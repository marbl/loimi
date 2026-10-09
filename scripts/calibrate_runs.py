#!/usr/bin/env python3
# SPDX-License-Identifier: MIT-0
"""Train technology-specific likelihoods only on the frozen training windows.

Truth labels and accepted votes come from the real HPC path. Each unique
placement locus contributes unit mass regardless of coverage. Estimate
P(observed length|true length,base,tier), not a truth-length prior: the label
sampling intentionally enriches long runs. Validation never updates models.
"""
import collections
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path

ROOT = None  # Set once by --data; all I/O stays in this dataset directory.
CAP = 100
MAX_TRUE = 50
MIN_LOCI = 30
MIN_DEPTH = 5
PRIOR_LOCI = 10
OUTLIER = .001
TIERS = ('LA', 'LAO', 'UL')


def checksum(path):
    """Hash large evidence files without retaining their contents."""
    digest = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def mode(h):
    """Use the same absolute-deviation/shorter-length tie break as C++."""
    peak = max(h.values())
    return min((sum(abs(k - v) * n for v, n in h.items()), k)
               for k, count in h.items() if count == peak)[1]


def labels(folder):
    """Map zero-based fixture HPC indexes to independent truth labels."""
    result = {}
    with (folder / 'labels.tsv').open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            run, length = int(row['run']), int(row['truth_length'])
            if run < 0 or run in result or length < 1 or row['base'] not in 'ACGT':
                raise ValueError(f'Invalid or duplicate truth label in {folder}: {row}')
            result[run] = (row['base'], length)
    return result


def evidence(folder):
    """Stream one histogram per run/tier/placement class in native export order."""
    truth = labels(folder)
    with (folder / 'consensus.run-evidence.tsv').open() as f:
        rows = csv.DictReader(f, delimiter='\t')
        key = lambda r: (int(r['run']), r['tier'], r['placement_class'])
        for (run, tier, kind), group in itertools.groupby(rows, key):
            # Truth masks may exclude additional ambiguous associations after
            # evidence collection. Unlabeled exports cannot become training data.
            if run not in truth:
                continue
            base, length = truth[run]
            h = {}
            for r in group:
                observed, count = int(r['length']), int(r['count'])
                if r['base'] != base or observed < 1 or count < 1 or observed in h:
                    raise ValueError(f'Invalid evidence row in {folder}: {r}')
                if tier not in TIERS or kind not in ('UNIQUE', 'AMBIGUOUS'):
                    raise ValueError(f'Unknown tier or placement class: {r}')
                h[observed] = count
            yield run, base, length, tier, kind, h


def train(windows):
    """Fit conditional distributions and freeze their provenance before scoring."""
    if (ROOT / 'models/metadata.json').exists():
        raise SystemExit('Models already frozen; use a new dataset directory to refit')
    counts = collections.defaultdict(collections.Counter)
    loci = collections.Counter()
    reads = collections.Counter()
    for w in windows:
        if w['split'] != 'train':
            continue
        for run, base, truth, tier, kind, h in evidence(ROOT / w['name']):
            if kind != 'UNIQUE' or truth > MAX_TRUE:
                continue
            key = (tier, base, truth)
            n = sum(h.values())
            loci[key] += 1
            reads[key] += n
            for obs, count in h.items():
                counts[key][min(obs, CAP)] += count / n
    modeldir = ROOT / 'models'
    modeldir.mkdir(exist_ok=True)
    metadata = {
        'format': 'HPC_RUN_MODEL_V1',
        'training_windows': [w['name'] for w in windows if w['split'] == 'train'],
        'validation_windows': [w['name'] for w in windows if w['split'] == 'validation'],
        'observed_max': CAP,
        'minimum_prediction_depth': MIN_DEPTH,
        'minimum_training_loci': MIN_LOCI,
        'maximum_true_length': MAX_TRUE,
        'prior_equivalent_loci': PRIOR_LOCI,
        'outlier_fraction': OUTLIER,
        'training_weight': 'one unit per unique-placement locus, not per read',
        'sampling': 'all labeled truth runs >=10; 1/100 shorter runs; no inferred length prior',
        'models': {},
    }
    metadata['manifest_sha256'] = checksum(ROOT / 'manifest.json')
    metadata['dataset_manifest'] = json.loads((ROOT / 'manifest.json').read_text())
    metadata['trainer_sha256'] = checksum(Path(__file__))
    metadata['training_input_sha256'] = {}
    for w in windows:
        if w['split'] == 'train':
            for filename in ('labels.tsv', 'consensus.run-evidence.tsv'):
                relative = f'{w["name"]}/{filename}'
                metadata['training_input_sha256'][relative] = checksum(ROOT / relative)
    for tier in TIERS:
        ranges = {}
        details = {}
        path = modeldir / f'{tier}.run-model.tsv'
        with path.open('w') as f:
            f.write('# Conditional HPC run-length model; experimental, chemistry/dataset-specific.\n')
            f.write(f'HPC_RUN_MODEL_V1 {CAP} {MIN_DEPTH}\n')
            for base in 'ACGT':
                upper = 0
                for truth in range(1, MAX_TRUE + 1):
                    key = (tier, base, truth)
                    n = loci[key]
                    if n < MIN_LOCI:
                        break
                    family = 'AT' if base in 'AT' else 'CG'
                    pool = collections.Counter()
                    for b in family:
                        pool.update(counts[(tier, b, truth)])
                    total = sum(pool.values())
                    probabilities = []
                    for obs in range(CAP + 1):
                        empirical = (counts[key][obs] + PRIOR_LOCI * pool[obs] / total) / (n + PRIOR_LOCI)
                        probabilities.append((1 - OUTLIER) * empirical + OUTLIER / (CAP + 1))
                    if abs(sum(probabilities) - 1) >= 1e-8:
                        raise ValueError('Fitted probabilities do not sum to one')
                    f.write(f'{base} {truth} {n} ' +
                            ' '.join(format(p, '.17g') for p in probabilities) + '\n')
                    upper = truth
                ranges[base] = upper
                details[base] = {
                    str(t): {'loci': loci[(tier, base, t)], 'observations': reads[(tier, base, t)]}
                    for t in range(1, MAX_TRUE + 1) if loci[(tier, base, t)]
                }
        metadata['models'][tier] = {
            'path': str(path), 'sha256': checksum(path),
            'supported_lengths_1_to': ranges, 'training_counts': details,
        }
        if not any(ranges.values()):
            raise ValueError(f'No sufficiently supported model rows for {tier}')
    (modeldir / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')


class Predictor:
    """Match RunModel::estimate, including depth, range, overflow and tie rules."""

    def __init__(self, path):
        self.rows = collections.defaultdict(list)
        self.cache = {}
        with path.open() as f:
            header = next(line for line in f if not line.startswith('#')).split()
            if header[0] != 'HPC_RUN_MODEL_V1':
                raise ValueError('Unsupported model version')
            self.cap, self.depth = map(int, header[1:])
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                fields = line.split()
                base, length = fields[0], int(fields[1])
                if length != len(self.rows[base]) + 1:
                    raise ValueError('Non-contiguous model lengths')
                self.rows[base].append([math.log(float(p)) for p in fields[3:]])

    def call(self, base, h, fallback):
        rows = self.rows[base]
        if sum(h.values()) < self.depth or not 1 <= fallback <= len(rows):
            return fallback
        key = (base, tuple(sorted(h.items())))
        if key not in self.cache:
            scores = [sum(n * p[min(o, self.cap)] for o, n in h.items()) for p in rows]
            self.cache[key] = max(range(len(scores)), key=scores.__getitem__) + 1
        return self.cache[key]


def record(s, truth, old, new, weight):
    """Keep unweighted counts distinct from inverse-sampling-weighted estimates."""
    s['loci'] += 1
    s['mode_wrong'] += old != truth
    s['model_wrong'] += new != truth
    s['mode_absolute_error'] += abs(old - truth)
    s['model_absolute_error'] += abs(new - truth)
    s['fixes_mode'] += old != truth and new == truth
    s['breaks_mode'] += old == truth and new != truth
    s['changed_calls'] += old != new
    s['sampling_weighted_loci'] += weight
    s['sampling_weighted_mode_wrong'] += weight * (old != truth)
    s['sampling_weighted_model_wrong'] += weight * (new != truth)


def validate(windows):
    """Score frozen models; neither distributions nor decision rules are fitted."""
    metadata = json.loads((ROOT / 'models/metadata.json').read_text())
    if checksum(ROOT / 'manifest.json') != metadata['manifest_sha256']:
        raise ValueError('Dataset manifest changed after fitting')
    for key, split in (('training_windows', 'train'), ('validation_windows', 'validation')):
        if metadata[key] != [w['name'] for w in windows if w['split'] == split]:
            raise ValueError('Dataset split changed after fitting')
    for tier in TIERS:
        path = ROOT / 'models' / f'{tier}.run-model.tsv'
        if checksum(path) != metadata['models'][tier]['sha256']:
            raise ValueError(f'Model changed after fitting: {tier}')
    predictors = {t: Predictor(ROOT / 'models' / f'{t}.run-model.tsv') for t in TIERS}
    report = {'windows': {}, 'by_tier': {}, 'priority_pipelines': {}}
    total = collections.defaultdict(collections.Counter)
    pipelines = collections.defaultdict(collections.Counter)
    for w in windows:
        if w['split'] != 'validation':
            continue
        folder = ROOT / w['name']
        per = collections.defaultdict(collections.Counter)
        pipe = collections.defaultdict(collections.Counter)
        observed_runs = set()
        for run, rows in itertools.groupby(evidence(folder), lambda x: x[0]):
            observed_runs.add(run)
            calls = {}
            truth = None
            for _, base, truth, tier, kind, h in rows:
                old = mode(h)
                new = predictors[tier].call(base, h, old)
                depth = sum(h.values())
                weight = 1 if truth >= 10 else 100
                depth_bin = 'depth5plus' if depth >= 5 else 'depth1to4'
                length_bin = ('length1to9' if truth < 10 else
                              'length10to19' if truth < 20 else 'length20plus')
                keys = (f'{tier}/{kind}/all', f'{tier}/{kind}/{depth_bin}', f'{tier}/{kind}/{length_bin}')
                for key in keys:
                    record(per[key], truth, old, new, weight)
                    record(total[key], truth, old, new, weight)
                    applicable = (depth >= predictors[tier].depth and
                                  1 <= old <= len(predictors[tier].rows[base]))
                    per[key]['model_eligible_loci'] += applicable
                    total[key]['model_eligible_loci'] += applicable
                calls[(tier, kind)] = (old, new)
            for name, order in [('LA-first', TIERS), ('LAO-first', ('LAO', 'LA', 'UL'))]:
                choices = [calls[(tier, kind)]
                           for kind in ('UNIQUE', 'AMBIGUOUS')
                           for tier in order if (tier, kind) in calls]
                if choices:
                    old, new = choices[0]
                    record(pipe[name], truth, old, new, weight)
                    record(pipelines[name], truth, old, new, weight)
        # Include labelable sites with no accepted evidence in the complete
        # priority simulation. Production's default fallback is length one;
        # omitting these sites would overstate accuracy in difficult regions.
        missing = 0
        for run, (base, truth) in labels(folder).items():
            if run in observed_runs:
                continue
            missing += 1
            weight = 1 if truth >= 10 else 100
            for name in ('LA-first', 'LAO-first'):
                record(pipe[name], truth, 1, 1, weight)
                record(pipelines[name], truth, 1, 1, weight)
        report['windows'][w['name']] = {
            'by_tier': per, 'priority_pipelines': pipe,
            'sampled_no_evidence_loci': missing,
            'label_summary': json.loads((folder / 'label-summary.json').read_text()),
        }
    report['by_tier'] = total
    report['priority_pipelines'] = pipelines
    report['caveats'] = [
        'HPC-labeled sampled runs, not full-sequence QV.',
        'Training uses only unique placements; ambiguous validation is separately reported.',
        'Weighted counts correct 1/100 short-run sampling only; excluded mappings remain excluded.',
        'No independence assumption for statistical significance; reads/loci are correlated.',
        'One individual and chemistry mixture; these models are not universal.']
    (ROOT / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True,
                        help='Dataset with manifest.json and window folders')
    parser.add_argument('action', choices=['train', 'validate'])
    args = parser.parse_args()
    ROOT = args.data.resolve()
    windows = json.loads((ROOT / 'manifest.json').read_text())['windows']
    names = [w['name'] for w in windows]
    if (len(names) != len(set(names)) or any(Path(n).name != n or n in ('.', '..') for n in names)
            or any(w['split'] not in ('train', 'validation') for w in windows)):
        raise SystemExit('Invalid window names or splits')
    if {w['split'] for w in windows} != {'train', 'validation'}:
        raise SystemExit('Both training and validation windows are required')
    split = 'train' if args.action == 'train' else 'validation'
    required = [w for w in windows if w['split'] == split]
    for w in required:
        if not all((ROOT / w['name'] / suffix).exists()
                   for suffix in ('consensus.summary.json', 'consensus.run-evidence.tsv')):
            raise SystemExit('Missing window ' + w['name'])
    {'train': train, 'validate': validate}[args.action](windows)
