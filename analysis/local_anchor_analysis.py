#!/usr/bin/env python3
# SPDX-License-Identifier: MIT-0
"""Collect bounded-site native diagnostics and replay local-anchor eligibility.

Pilot: training w00/w02, validation w01/w11. Uniform deterministic subsample
of existing truth-labelable loci plus a separately reported, biased diagnostic
set of native 1-kb-changed calls. This is run-length evidence replay, not native
polished assembly quality. No truth is used to select read observations.
"""
import argparse
import collections
import concurrent.futures
import csv
import json
import subprocess
from pathlib import Path
from corroboration import c, select, record, REPO

POLICIES = ('baseline', 'span1000', 'anchor2', 'anchor5', 'anchor10', 'anchor20')
WINDOWS = {'w00': 'train', 'w02': 'train', 'w01': 'validation', 'w11': 'validation'}


def selected_sites(root, name, all_labels=False):
    labels = c.labels(root/'calibration20'/name)
    # Multiplicative integer mixing avoids periodic artifacts in sampled HPC indexes.
    random = set(labels) if all_labels else {r for r in labels if ((r*2654435761) & 0xffffffff) % 10 == 0}
    changed = set()
    native = root/'flank-retention10-resumed'/'LA-LAO-UL'/name
    if (native/'flank1000.run-calls.tsv').exists():
        with (native/'baseline.run-calls.tsv').open() as x, (native/'flank1000.run-calls.tsv').open() as y:
            next(x); next(y)
            for a, b in zip(x, y, strict=True):
                aa, bb = a.split('\t'), b.split('\t')
                assert aa[1] == bb[1]
                run = int(aa[1])
                if run in labels and aa[5] != bb[5]:
                    changed.add(run)
    return labels, random, changed


def collect(root, out, exe, name, all_labels=False):
    folder = out/name; folder.mkdir(exist_ok=True)
    labels, random, changed = selected_sites(root, name, all_labels)
    sites = random | changed
    primary = 'labelable' if all_labels else 'random'
    (folder/'selection.json').write_text(json.dumps({primary: sorted(random), 'span_changed': sorted(changed)})+'\n')
    (folder/'sites.tsv').write_text(''.join(f'tig00000109\t{r}\n' for r in sorted(sites)))
    cmd = [str(exe), '--import', str(root/'calibration20'/name/'input.cnspack'),
           '--template-hpc', str(root/'calibration20'/name/'template.fasta'),
           '--tier', 'LA=hifi_m*', '--tier', 'LAO=hifi_[0-9a-ln-z]*', '--tier', 'UL=ont_*',
           '--tier-model', f'LA={REPO}/models/revio.run-model.tsv', '--threads', '4',
           '--no-extend-ends', '--experiment-retain-higher-unique',
           '--evidence-sites', str(folder/'sites.tsv'), '--evidence-read-votes',
           '--output', str(folder/'trace')]
    (folder/'command.json').write_text(json.dumps(cmd, indent=2)+'\n')
    if not (folder/'trace.summary.json').exists():
        with (folder/'trace.log').open('w') as log:
            subprocess.run(cmd, stdout=log, stderr=log, check=True)
    print(name, 'diagnostics complete', flush=True)


def analyze(root, out, name):
    labels = c.labels(root/'calibration20'/name)
    selection = json.loads((out/name/'selection.json').read_text())
    groups = {k: set(v) for k, v in selection.items()}
    hist = {p: collections.defaultdict(lambda: collections.defaultdict(collections.Counter)) for p in POLICIES}
    strands = collections.defaultdict(lambda: [collections.Counter(), collections.Counter()])
    with (out/name/'trace.run-read-votes.tsv').open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['status'] != 'ELIGIBLE':
                continue
            run = int(row['run']); length = int(row['length'])
            kind = 'UNIQUE' if int(row['placements']) == 1 else 'AMBIGUOUS'
            key = row['tier'], kind
            left, right = int(row['left_anchor_runs']), int(row['right_anchor_runs'])
            if key == ('LA', 'UNIQUE'):
                strands[run][int(row['reverse'])][length] += 1
            for policy in POLICIES:
                eligible = policy == 'baseline'
                if policy == 'span1000':
                    eligible = min(left, right) >= 1 and min(int(row['left_aligned_bases']), int(row['right_aligned_bases'])) >= 1000
                if policy.startswith('anchor'):
                    eligible = min(left, right) >= int(policy[6:])
                if eligible:
                    hist[policy][run][key][length] += 1
    model = c.Predictor(REPO/'models/revio.run-model.tsv')
    def call(hs, order, base):
        tier, kind, h = select(hs, order)
        modal = c.mode(h) if h else 1
        call = model.call(base, h, modal) if tier == 'LA' else modal
        if kind == 'UNIQUE':
            combined = collections.Counter(h)
            for t in order[:order.index(tier)]:
                combined.update(hs.get((t, 'UNIQUE'), {}))
            if combined != h:
                call = c.mode(combined)
        return call, tier, kind, sum(h.values())
    stats = collections.defaultdict(collections.Counter)
    changes = []
    selected = set().union(*groups.values())
    native = {}
    with (out/name/'trace.run-calls.tsv').open() as f:
        next(f)
        for line in f:
            fields = line.split('\t')
            run = int(fields[1])
            if run in selected:
                native[run] = int(fields[5])
    for run in selected:
        base, truth = labels[run]
        f, r = strands[run]
        discordant = sum(f.values()) >= 5 and sum(r.values()) >= 5 and c.mode(f) != c.mode(r)
        for config, order in [('LA-LAO-UL', ('LA', 'LAO', 'UL')), ('LA-UL', ('LA', 'UL'))]:
            old, oldtier, oldkind, olddepth = call(hist['baseline'][run], order, base)
            if config == 'LA-LAO-UL' and old != native[run]:
                raise AssertionError(f'Native/replay mismatch {name}:{run}: {native[run]} vs {old}')
            for policy in POLICIES:
                new, tier, kind, depth = call(hist[policy][run], order, base)
                for group, runs in groups.items():
                    if run not in runs:
                        continue
                    s = stats[f'{config}/{group}/{policy}']
                    record(s, truth, old, new)
                    s['tier_switches'] += tier != oldtier
                    s['lost_LA_depth5'] += oldtier == 'LA' and oldkind == 'UNIQUE' and olddepth >= 5 and (tier != 'LA' or kind != 'UNIQUE' or depth < 5)
                    s['strand_discordant_loci'] += discordant
                if old != new:
                    changes.append(dict(run=run, config=config, policy=policy, old=old, new=new,
                                        truth=truth, tier_before=oldtier, tier_after=tier,
                                        depth_before=olddepth, depth_after=depth, strand_discordant=discordant))
    (out/name/'analysis.json').write_text(json.dumps(dict(stats=stats, changes=changes), indent=2)+'\n')
    return stats


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--executable', type=Path, required=True)
    p.add_argument('--validation-all', action='store_true',
                   help='Evaluate fixed alternatives on all ten historical validation windows; do not select a rule from validation')
    p.add_argument('--all-labels', action='store_true',
                   help='Use all existing truth labels, without the pilot subsample')
    a = p.parse_args()
    root, out, exe = a.input_root.resolve(), a.output.resolve(), a.executable.resolve()
    out.mkdir(parents=True, exist_ok=True)
    windows = WINDOWS
    if a.validation_all:
        source = json.loads((root/'calibration20'/'manifest.json').read_text())
        windows = {w['name']: w['split'] for w in source['windows'] if w['split'] == 'validation'}
    manifest = dict(windows=windows, policies=POLICIES, binary_sha256=c.checksum(exe),
                    model_sha256=c.checksum(REPO/'models/revio.run-model.tsv'),
                    script_sha256=c.checksum(Path(__file__)),
                    source_manifest_sha256=c.checksum(root/'calibration20'/'manifest.json'),
                    scope=__doc__, all_labels=a.all_labels,
                    validation_only=a.validation_all,
                    selection='Fixed alternatives; no validation-based selection' if a.validation_all else 'Training selection')
    path = out/'manifest.json'
    if path.exists():
        assert json.loads(path.read_text()) == json.loads(json.dumps(manifest)), 'Changed provenance'
    else:
        path.write_text(json.dumps(manifest, indent=2)+'\n')
    # Complete training and freeze selection before validation analysis.
    totals = {}; chosen = {}
    for split in ('train', 'validation'):
        names = [n for n, s in windows.items() if s == split]
        if not names:
            continue
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            list(pool.map(lambda n: collect(root, out, exe, n, a.all_labels), names))
        total = collections.defaultdict(collections.Counter)
        partitions = {'pilot_windows': collections.defaultdict(collections.Counter),
                      'additional_windows': collections.defaultdict(collections.Counter)}
        for name in names:
            for key, value in analyze(root, out, name).items():
                total[key].update(value)
                partition = 'pilot_windows' if name in ('w01', 'w11') else 'additional_windows'
                partitions[partition][key].update(value)
            print(name, 'replay and native parity complete', flush=True)
        totals[split] = total
        if split == 'train':
            for config in ('LA-LAO-UL', 'LA-UL'):
                primary = 'labelable' if a.all_labels else 'random'
                chosen[config] = min(POLICIES, key=lambda p: (total[f'{config}/{primary}/{p}']['weighted_wrong'], total[f'{config}/{primary}/{p}']['changes']))
            (out/'frozen.json').write_text(json.dumps(dict(chosen=chosen, training=total), indent=2)+'\n')
        (out/'results.json').write_text(json.dumps(dict(chosen=chosen, totals=totals, validation_partitions=partitions if split == 'validation' else {}), indent=2)+'\n')
    print('LOCAL ANCHOR ANALYSIS COMPLETE', chosen, flush=True)


if __name__ == '__main__':
    main()
