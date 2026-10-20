#!/usr/bin/env python3
"""Frozen-label replay of conservative LA model overrides; no model fitting.

This is separate from borrowing: compare unchanged model, raw mode, a 100:1
likelihood requirement to override mode, and protection of >=60% modal support.
These hypothesis-driven fixed rules are exploratory, not independent validation.
"""
import argparse
import collections
import importlib.util
import itertools
import json
import math
from pathlib import Path

SPEC = importlib.util.spec_from_file_location('calibration', Path(__file__).resolve().parents[1]/'scripts/calibrate_runs.py')
calibration = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(calibration)


def select(histograms, order):
    """First qualifying unique tier, then qualifying pooled, then sparse."""
    for pooled in (False, True):
        for tier in order:
            h = collections.Counter(histograms.get((tier, 'UNIQUE'), {}))
            if pooled:
                h.update(histograms.get((tier, 'AMBIGUOUS'), {}))
            if sum(h.values()) >= 5:
                return tier, h
    for tier in order:
        h = collections.Counter(histograms.get((tier, 'UNIQUE'), {}))
        h.update(histograms.get((tier, 'AMBIGUOUS'), {}))
        if h:
            return tier, h
    return None, {}


def calls(model, base, h, tier):
    mode = calibration.mode(h) if h else 1
    current = model.call(base, h, mode) if tier == 'LA' else mode
    confident = protected = current
    if current != mode:
        rows = model.rows[base]
        delta = sum(n*(rows[current-1][min(k, model.cap)]-rows[mode-1][min(k, model.cap)]) for k, n in h.items())
        if delta < math.log(100):
            confident = mode
        if h[mode]/sum(h.values()) >= .6:
            protected = mode
    return dict(baseline=current, mode=mode, odds100=confident, protect_mode60=protected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--model', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    model = calibration.Predictor(args.model)
    total = collections.defaultdict(collections.Counter)
    windows = {}
    for w in json.loads((args.data/'manifest.json').read_text())['windows']:
        if w['split'] != 'validation':
            continue
        stats = collections.defaultdict(collections.Counter)
        seen = set()

        def score(base, truth, histograms):
            for name, order in [('LA_LAO_UL', ('LA', 'LAO', 'UL')), ('LA_UL', ('LA', 'UL'))]:
                tier, h = select(histograms, order)
                predictions = calls(model, base, h, tier)
                for method, value in predictions.items():
                    key = name+'/'+method
                    for target in (stats[key], total[key]):
                        target['loci'] += 1
                        target['wrong'] += value != truth
                        target['absolute_error'] += abs(value-truth)
                        target['fixes'] += predictions['baseline'] != truth and value == truth
                        target['breaks'] += predictions['baseline'] == truth and value != truth
                        target['weighted_wrong'] += (100 if truth < 10 else 1)*(value != truth)
                        target['changes'] += value != predictions['baseline']

        for run, rows in itertools.groupby(calibration.evidence(args.data/w['name']), lambda row: row[0]):
            histograms = {}
            seen.add(run)
            for _, base, truth, tier, kind, h in rows:
                histograms[(tier, kind)] = h
            score(base, truth, histograms)
        for run, (base, truth) in calibration.labels(args.data/w['name']).items():
            if run not in seen:
                score(base, truth, {})
        windows[w['name']] = stats
        model.cache.clear()
        print(w['name'], 'complete', flush=True)
    result = dict(total=total, windows=windows,
                  model_sha256=calibration.checksum(args.model),
                  caveat='Sampled labelable runs, not whole-sequence QV. No borrowing. Short runs sampled 1/100; weighted estimates are separate. Rules motivated by inspected diagnostic sites; require independent validation.')
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(total, indent=2))


if __name__ == '__main__':
    main()
