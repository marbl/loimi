#!/usr/bin/env python3
# SPDX-License-Identifier: MIT-0
"""Streaming, offline tier-corroboration and confidence experiment.

Never changes assemblies/models. Train selects one rule before validation is
opened. Historical validation is exploratory, not a pristine holdout. One lower
tier is consulted, so corrected/raw ONT are never pooled as independent votes.
Short truth runs were sampled 1/100; weighted counts are estimates, not QV.
"""
import argparse
import collections
import functools
import importlib.util
import itertools
import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('calibration', REPO/'scripts/calibrate_runs.py')
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


def select(hist, order):
    for pooled in (False, True):
        for tier in order:
            h = collections.Counter(hist.get((tier, 'UNIQUE'), {}))
            if pooled:
                h.update(hist.get((tier, 'AMBIGUOUS'), {}))
            if sum(h.values()) >= 5:
                return tier, 'POOLED' if pooled else 'UNIQUE', h
    for tier in order:
        h = collections.Counter(hist.get((tier, 'UNIQUE'), {}))
        h.update(hist.get((tier, 'AMBIGUOUS'), {}))
        if h:
            return tier, 'SPARSE', h
    return None, 'NONE', {}


def record(s, truth, old, new):
    weight = 100 if truth < 10 else 1
    s['loci'] += 1
    s['wrong'] += new != truth
    s['absolute_error'] += abs(new-truth)
    s['weighted_wrong'] += weight*(new != truth)
    s['weighted_loci'] += weight
    s['changes'] += new != old
    s['fixes'] += old != truth and new == truth
    s['breaks'] += old == truth and new != truth


def bucket(value, cuts):
    return str(next((x for x in cuts if value < x), 'inf'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data', type=Path, required=True)
    p.add_argument('--lao-model', type=Path, required=True)
    p.add_argument('--ul-model', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    paths = {'LA': REPO/'models/revio.run-model.tsv', 'LAO': a.lao_model, 'UL': a.ul_model}
    models = {t: c.Predictor(path) for t, path in paths.items()}

    @functools.lru_cache(maxsize=8192)
    def scores(tier, base, items):
        model = models[tier]
        return tuple(sum(n*row[min(obs, model.cap)] for obs, n in items)
                     for row in model.rows[base])

    def call(tier, base, h):
        modal = c.mode(h) if h else 1
        model = models.get(tier)
        if model and sum(h.values()) >= model.depth and 1 <= modal <= len(model.rows[base]):
            s = scores(tier, base, tuple(sorted(h.items())))
            return max(range(len(s)), key=s.__getitem__)+1
        return modal

    rules = {'baseline': None}
    for margin, strength, fraction in itertools.product((1., 3., 10.), (.1, .5), (.6, .8)):
        rules[f'gap{margin}-lower{strength}-fraction{fraction}'] = (margin, strength, fraction)
    manifest = json.loads((a.data/'manifest.json').read_text())
    provenance = dict(models={t: c.checksum(path) for t, path in paths.items()},
                      manifest_sha256=c.checksum(a.data/'manifest.json'), rules=rules,
                      selection='training minimum weighted wrong; ties fewer changed calls',
                      caveat='Existing model training windows reused for confidence fitting; historical validation previously inspected. No full sequence/phase/QV claim.')
    (a.output/'design.json').write_text(json.dumps(provenance, indent=2)+'\n')
    totals = {}; calibration = {}; chosen = {}; window_results = {}
    for split in ('train', 'validation'):
        stats = collections.defaultdict(collections.Counter)
        bins = collections.defaultdict(collections.Counter)
        for w in manifest['windows']:
            if w['split'] != split:
                continue
            ws = collections.defaultdict(collections.Counter)
            for run, rows in itertools.groupby(c.evidence(a.data/w['name']), lambda r: r[0]):
                hist = {}
                for _, base, truth, tier, kind, h in rows:
                    hist[tier, kind] = h
                for config, order in [('LA-LAO-UL', ('LA', 'LAO', 'UL')), ('LA-UL', ('LA', 'UL'))]:
                    tier, kind, h = select(hist, order)
                    old = call('LA' if tier == 'LA' else None, base, h)
                    # Match additive retention: a mixed-technology histogram
                    # uses mode rather than pretending all reads are one tech.
                    if kind == 'UNIQUE':
                        pooled = collections.Counter(h)
                        for higher in order[:order.index(tier)]:
                            pooled.update(hist.get((higher, 'UNIQUE'), {}))
                        if pooled != h:
                            old = c.mode(pooled)
                    updates = {}; feature = None
                    if tier == 'LA' and kind == 'UNIQUE' and 1 <= c.mode(h) <= len(models['LA'].rows[base]):
                        s = scores('LA', base, tuple(sorted(h.items())))
                        ranked = sorted(range(len(s)), key=lambda i: (-s[i], i))
                        first, second = ranked[:2]; gap = s[first]-s[second]
                        assert old == first+1
                        modal = c.mode(h)
                        margin_bin = bucket(gap, (1, 3, 10, 30, 100))
                        feature = (f'{config}/gap<{margin_bin}',
                                   f'{config}/{base}/depth<{bucket(sum(h.values()), (10,20,40))}/length<{bucket(old,(10,20,40))}/gap<{margin_bin}/disagree{old!=modal}')
                        confidence = 1/sum(math.exp(x-s[first]) for x in s)
                        for key in feature:
                            b = bins[key]; b['loci'] += 1; b['wrong'] += old != truth
                            b['predicted_error_sum'] += 1-confidence
                            b['naive_brier_sum'] += ((1-confidence)-(old != truth))**2
                            if split == 'validation':
                                fitted = calibration['train'].get(feature[1], {})
                                if fitted.get('loci', 0) < 100:
                                    fitted = calibration['train'].get(feature[0], {})
                                risk = (fitted.get('wrong', 0)+1)/(fitted.get('loci', 0)+2)
                                b['fitted_error_sum'] += risk
                                b['fitted_brier_sum'] += (risk-(old != truth))**2
                        for lower in order[1:]:
                            lh = hist.get((lower, 'UNIQUE'), {})
                            if sum(lh.values()) < 5 or max(first, second) >= len(models[lower].rows[base]):
                                continue
                            ls = scores(lower, base, tuple(sorted(lh.items())))
                            strength = (ls[second]-ls[first])/sum(lh.values())
                            fraction = lh.get(second+1, 0)/sum(lh.values())
                            for label, rule in rules.items():
                                if rule and gap <= rule[0] and strength >= rule[1] and fraction >= rule[2]:
                                    updates[label] = second+1
                            # Never add UL to LAO: they may be the same molecule.
                            break
                    for label in rules:
                        new = updates.get(label, old)
                        key = config+'/'+label
                        record(stats[key], truth, old, new); record(ws[key], truth, old, new)
            window_results[w['name']] = ws
            print(split, w['name'], 'complete', flush=True)
        totals[split] = stats; calibration[split] = bins
        if split == 'train':
            for config in ('LA-LAO-UL', 'LA-UL'):
                chosen[config] = min(rules, key=lambda rule: (stats[config+'/'+rule]['weighted_wrong'], stats[config+'/'+rule]['changes']))
            # Freeze before reading validation evidence/labels.
            (a.output/'frozen.json').write_text(json.dumps(dict(chosen=chosen, calibration=bins), indent=2)+'\n')
        (a.output/'results.json').write_text(json.dumps(dict(chosen=chosen, totals=totals, windows=window_results, confidence=calibration), indent=2)+'\n')
    print('Selected on training:', chosen, flush=True)
    for config, label in chosen.items():
        print(config, 'baseline', dict(totals['validation'][config+'/baseline']), 'selected', dict(totals['validation'][config+'/'+label]), flush=True)


if __name__ == '__main__':
    main()
