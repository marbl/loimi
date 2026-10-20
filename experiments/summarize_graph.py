#!/usr/bin/env python3
"""Summarize diagnostic graphs without treating selected cases as genome QV."""
import argparse
import collections
import json
from pathlib import Path


def summarize(rows):
    summary = {}
    for kind in ('error', 'control'):
        cases = [row for row in rows if row['case']['kind'] == kind]
        groups = {}
        for arm in ('LA_LAO_UL', 'LA_UL'):
            counts = collections.Counter()
            for row in cases:
                policy = row['policies'][arm]
                before = policy.get('calibrated_correct')
                after = policy.get('graph_correct')
                counts['cases'] += 1
                counts['calibrated_correct'] += before is True
                counts['graph_correct'] += after is True
                counts['graph_fixes'] += before is False and after is True
                counts['graph_regressions'] += before is True and after is False
                counts['unresolved_graph'] += after is None
                counts['rejected_read_chunks'] += len(policy['rejected_reads'])
                # Descriptive strand-mode discordance, not a proposed filter.
                modes = []
                for histogram in policy['strand_histograms'].values():
                    if sum(histogram.values()) >= 5:
                        peak = max(histogram.values())
                        modes.append({int(k) for k, n in histogram.items() if n == peak})
                if len(modes) == 2:
                    counts['both_strands_depth5'] += 1
                    counts['disjoint_strand_modes'] += not bool(modes[0] & modes[1])
            groups[arm] = dict(counts)
        summary[kind] = dict(cases=len(cases),
                             mature_correct=sum(row['case']['mature_correct'] is True for row in cases),
                             policies=groups)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    result = summarize(json.loads((args.directory/'results.json').read_text()))
    (args.directory/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
