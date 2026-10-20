#!/usr/bin/env python3
"""Truth-selected, matched-read diagnostic; never edits an assembly.

Select up to two anchored HP errors and two nearby correct controls per window.
Compare raw mode, calibrated LA HPC calls and the external Verkko graph on
identical eligible unique read chunks. LAO removal is an explicit second arm.
Read chunks come from the uploaded CRAM, not a fresh native consensus alignment;
thus this diagnoses algorithms/evidence, not full pipeline quality or switches.
"""
import argparse
import collections
import json
import re
import subprocess
from pathlib import Path
from classify_remaining import fasta, run_bounds

def anchored(sequence, left, right):
    """Return the allele between unique exact anchors, or None if unresolved."""
    a, b = (sequence.find(left), sequence.find(right))
    if a < 0 or b < a + len(left) or b - a > 1000:
        return None
    if sequence.find(left, a + 1) >= 0 or sequence.find(right, b + 1) >= 0:
        return None
    return (a + len(left), b)

def chunks(cram, begin, end):
    """Extract forward-reference-oriented chunks via SAM CIGAR coordinates.

    SAM SEQ is already reverse-complemented for reverse-strand alignments.
    Do not complement it again; retain FLAG strand separately for diagnostics.
    """
    sam = subprocess.check_output(['samtools', 'view', str(cram), f'tig00000109:{begin + 1}-{end}'], text=True)
    reads = {}
    strands = {}
    for line in sam.splitlines():
        f = line.split('\t')
        flag = int(f[1])
        if int(f[4]) != 60 or flag & (4 | 256 | 2048) or f[9] == '*':
            continue
        tier = 'LA' if f[0].startswith('hifi_m') else 'LAO' if f[0].startswith('hifi_') else 'UL' if f[0].startswith('ont_') else None
        if tier is None:
            continue
        rp, qp = (int(f[3]) - 1, 0)
        mapped = {}
        for n, op in re.findall('(\\d+)([MIDNSHP=X])', f[5]):
            n = int(n)
            if op in 'M=X':
                for boundary in (begin - 1, end):
                    if rp <= boundary < rp + n:
                        mapped[boundary] = qp + boundary - rp
                rp += n
                qp += n
            elif op in 'IS':
                qp += n
            elif op in 'DN':
                rp += n
        if len(mapped) == 2:
            reads[f[0]] = (f[0], tier, f[9][mapped[begin - 1] + 1:mapped[end]])
            strands[f[0]] = 'reverse' if flag & 16 else 'forward'
    return (list(reads.values()), strands)

def evaluate(args, case, folder, reads, reference, strands):
    """Run both tier configurations and compare estimators on selected IDs."""
    output = {}
    for arm in ('LA_LAO_UL', 'LA_UL'):
        chosen = [r for r in reads if arm != 'LA_UL' or r[1] != 'LAO']
        path = folder / (arm + '.reads.tsv')
        path.write_text(''.join(('\t'.join(r) + '\n' for r in chosen)))
        command = [str(args.probe), str(folder / 'reference.txt'), str(path), str(case['begin']), str(case['end'])]
        text = subprocess.check_output(command, text=True)
        (folder / (arm + '.probe.tsv')).write_text(text)
        observations = {}
        names = []
        graph = None
        rejected = []
        for line in text.splitlines():
            f = line.split('\t')
            if f[0] == 'READ':
                observations[f[1]] = dict(tier=f[2], length=int(f[3]), allele=f[4])
            elif f[0] == 'REJECT':
                rejected.append(f[1])
            elif f[0] == 'SELECT' and f[1] == 'LAO_exclusive':
                names.append(f[2])
            elif f[0] == 'GRAPH' and f[1] == 'LAO_exclusive':
                graph = f[2]
        hist = collections.Counter((observations[n]['length'] for n in names))
        result = dict(depth=len(names), tiers=dict(collections.Counter((observations[n]['tier'] for n in names))), histogram=dict(hist), rejected_reads=rejected)
        result['strand_histograms'] = {strand: dict(collections.Counter((observations[n]['length'] for n in names if strands[n] == strand))) for strand in ('forward', 'reverse')}
        if hist:
            mode = min(hist, key=lambda k: (-hist[k], sum((abs(k - v) * n for v, n in hist.items())), k))
            predicted = mode
            if observations[names[0]]['tier'] == 'LA':
                line = case['base'] + ' ' + ' '.join((f'{k}:{v}' for k, v in sorted(hist.items()))) + '\n'
                predicted = int(subprocess.check_output([str(args.probe.parent / 'apply-run-model'), str(args.model)], input=line, text=True).split()[1])
            result.update(mode=mode, calibrated=predicted, calibrated_correct=predicted == case['truth_length'])
        if graph:
            bounds = anchored(graph, reference[case['begin'] - 40:case['begin']], reference[case['end']:case['end'] + 40])
            allele = graph[bounds[0]:bounds[1]] if bounds else None
            result.update(graph_allele=allele, graph_correct=allele == case['base'] * case['truth_length'] if allele is not None else None)
        output[arm] = result
    return output

def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('native-root', 'input-root', 'probe', 'model', 'output'):
        p.add_argument('--' + name, required=True, type=Path)
    p.add_argument('--require-mature-agreement', action='store_true', help='Diagnostic subset where mature Verkko and Q100 agree at the focal run')
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    records = json.loads((args.native_root / 'remaining-classification.json').read_text())['records']
    uploaded = fasta(args.input_root / 'part002_fast.fasta')
    mature = fasta(args.input_root / 'part002.fasta')
    report = []
    for w in json.loads((args.native_root / 'manifest.json').read_text())['windows']:
        folder = args.native_root / w['name']
        seq = fasta(folder / 'borrow.core.fasta')
        truth = fasta(folder / 'truth.fasta')
        candidates = []
        seen = set()
        for r in records:
            if r['window'] != w['name'] or r['category'] != 'hp_associated':
                continue
            a, b = run_bounds(seq, r['query_pos'])
            if b - a < 2 or a in seen:
                continue
            seen.add(a)
            candidates.append((a, b, 'error'))
        controls = []
        for a, b, _ in candidates[:30]:
            for match in re.finditer('A{10,}|C{10,}|G{10,}|T{10,}', seq[b + 100:b + 5000]):
                controls.append((b + 100 + match.start(), b + 100 + match.end(), 'control'))
                break
        counts = collections.Counter()
        used = []
        for a, b, kind in candidates + controls:
            if counts[kind] >= 2 or any((abs(a - v) < 100 for v in used)):
                continue
            left, right = (seq[a - 40:a], seq[b:b + 40])
            tb = anchored(truth, left, right)
            if not tb:
                continue
            allele = truth[tb[0]:tb[1]]
            if not allele or set(allele) != {seq[a]} or (len(allele) == b - a) != (kind == 'control'):
                continue
            ub = anchored(uploaded, left, right)
            if not ub or not uploaded[ub[0]:ub[1]] or set(uploaded[ub[0]:ub[1]]) != {seq[a]}:
                continue
            mb = anchored(mature, left, right)
            if args.require_mature_agreement and (not mb or mature[mb[0]:mb[1]] != allele):
                continue
            begin, end = (ub[0] - 500, ub[1] + 500)
            if begin < 1 or end >= len(uploaded):
                continue
            reads, strands = chunks(args.input_root / 'part002_fast.cram', begin, end)
            if len(reads) < 5:
                continue
            name = f"{w['name']}-{kind}-{counts[kind] + 1}"
            dest = args.output / name
            dest.mkdir(exist_ok=True)
            reference = uploaded[begin:end]
            (dest / 'reference.txt').write_text(reference + '\n')
            case = dict(name=name, kind=kind, window=w['name'], base=seq[a], truth_length=len(allele), native_length=b - a, uploaded_length=ub[1] - ub[0], uploaded_begin=ub[0], truth_begin=w['truth_raw_begin'] + tb[0], begin=500, end=500 + ub[1] - ub[0], mature_correct=mature[mb[0]:mb[1]] == allele if mb else None)
            (dest / 'strands.json').write_text(json.dumps(strands, indent=2) + '\n')
            result = dict(case=case, policies=evaluate(args, case, dest, reads, reference, strands))
            report.append(result)
            counts[kind] += 1
            used.append(a)
            print(name, {k: (v.get('calibrated_correct'), v.get('graph_correct'), v['tiers']) for k, v in result['policies'].items()}, flush=True)
            (args.output / 'results.json').write_text(json.dumps(report, indent=2) + '\n')
    print('Completed', len(report), 'diagnostic cases', flush=True)
if __name__ == '__main__':
    main()
