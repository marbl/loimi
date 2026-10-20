#!/usr/bin/env python3
"""Classify paired native Q100 differences without changing production calling.

Coordinates are zero-based, half-open unless a field says otherwise. HP and
short-repeat association describe sequence context, not a proven error cause.
Overlapping reference alignment blocks in either arm are excluded. Input files
are the ten-window runner's manifest, FASTAs, coords and show-snps outputs.
"""
import argparse
import bisect
import collections
import json
from pathlib import Path


def fasta(path):
    with path.open() as source:
        return ''.join(line.strip() for line in source if not line.startswith('>')).upper()


def unique_spans(path):
    events = collections.Counter()
    for line in path.read_text().splitlines():
        a, b = sorted(map(int, line.split()[:2]))
        events[a-1] += 1
        events[b] -= 1
    depth = 0
    last = 0
    result = []
    for p, delta in sorted(events.items()):
        if depth == 1 and last < p:
            result.append((last, p))
        last = p
        depth += delta
    return result


def contains(spans, starts, p):
    i = bisect.bisect_right(starts, p)-1
    return i >= 0 and p < spans[i][1]


def run_bounds(sequence, p):
    if not 0 <= p < len(sequence):
        raise ValueError('Position outside sequence')
    a, b = p, p+1
    while a and sequence[a-1] == sequence[p]:
        a -= 1
    while b < len(sequence) and sequence[b] == sequence[p]:
        b += 1
    return a, b


def tandem_period(sequence, p):
    """Smallest period 2..6 with >=3 exact copies overlapping the base.

    Restrict the search to an 18-base radius. Pure homopolymer motifs are
    excluded here because run-length context is reported separately.
    """
    for period in range(2, 7):
        for a in range(max(0, p-3*period+1), p+1):
            motif = sequence[a:a+period]
            if len(set(motif)) > 1 and sequence[a:a+3*period] == motif*3:
                return period
    return None


def classify(root):
    rows = []
    for window in json.loads((root/'manifest.json').read_text())['windows']:
        folder = root/window['name']
        masks = [unique_spans(folder/(arm+'.coords')) for arm in ('baseline', 'borrow')]
        starts = [[a for a, b in mask] for mask in masks]
        truth = fasta(folder/'truth.fasta')
        query = fasta(folder/'borrow.core.fasta')
        for line in (folder/'borrow.snps').read_text().splitlines():
            fields = line.split()
            rp, qp = int(fields[0])-1, int(fields[3])-1
            if not all(contains(mask, start, rp) for mask, start in zip(masks, starts)):
                continue
            rb, qb = fields[1:3]
            kind = 'insertion' if rb == '.' else 'deletion' if qb == '.' else 'substitution'
            seq, pos = (truth, rp) if kind == 'deletion' else (query, qp)
            a, b = run_bounds(seq, pos)
            hp = kind != 'substitution' and b-a >= 2
            period = tandem_period(seq, pos)
            category = 'hp_associated' if hp else 'short_repeat' if period else 'other'
            rows.append(dict(window=window['name'],truth_pos=window['truth_raw_begin']+rp,
                             truth_local=rp,query_pos=qp,kind=kind,category=category,
                             run_length=b-a,repeat_period=period))
    counts = collections.Counter((r['kind'], r['category']) for r in rows)
    return dict(counts={f'{k}/{c}':n for (k,c),n in sorted(counts.items())},
                total=len(rows),records=rows,
                caveat='Truth-selected diagnostics; context association is not causal classification. No switch assessment.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native-root', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = classify(args.native_root)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result['counts'], indent=2))


if __name__ == '__main__':
    main()
