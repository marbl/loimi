"""Multi-tig regression: overhangs, small samples, and zero unique reads."""
import pathlib
import random
import subprocess
import sys
import tempfile

exe = pathlib.Path(sys.argv[1])
tools = pathlib.Path(sys.argv[2]) / 'build/bin'


def run(args):
    # Legacy tiny fixtures cannot supply the production 500-bp flank.
    if str(args[0]) == str(exe):
        args = [*args, "--no-hp-flank-filter", "--no-hp-realign"]
    p = subprocess.run(list(map(str, args)), capture_output=True, text=True, cwd=root)
    assert p.returncode == 0, p.stdout + p.stderr
    return p.stdout + p.stderr


with tempfile.TemporaryDirectory(prefix='hpc-offset-') as temporary:
    root = pathlib.Path(temporary)
    rng = random.Random(913)
    bases = ['A']
    for _ in range(11999):
        bases.append(rng.choice([b for b in 'ACGT' if b != bases[-1]]))
    sequence = ''.join(bases)
    reference = sequence[3000:9000]
    # Sparse cross-tier reads have no long exact anchors. Their only useful
    # evidence is a noisy enclosing overlap, so the joint fallback is required.
    mutated = list(sequence)
    for i in range(125, len(mutated) - 1, 250):
        mutated[i] = next(b for b in 'ACGT'
                          if b not in (mutated[i-1], mutated[i], mutated[i+1]))
    mutated = ''.join(mutated)
    # Distinct read IDs across tigs preserve package-wide uniqueness. The
    # second tig has five agreeing reads, one offset outlier, and one read
    # with no exact internal anchors. It must not require six agreeing votes.
    noise = ''.join(rng.choice('ACGT') for _ in range(12000))
    # One unique enclosing LA read must not be vetoed by two conflicting
    # placements of a second LA molecule. The template has no unique 128-base
    # markers, and both reads have too many distributed errors for long chunks.
    unit = list(sequence[:2662])
    unit[-1] = next(b for b in 'ACGT' if b not in (unit[-2], unit[0]))
    unit = ''.join(unit)
    reference_e = (unit * 3)[:6278]
    prefix = list(sequence[8000:8208])
    prefix[-1] = next(b for b in 'ACGT' if b not in (prefix[-2], unit[0]))
    sparse_read = list(''.join(prefix) + (unit * 8)[:20378])
    for i in range(209, len(sparse_read)-1):
        noisy = i < 208 + len(reference_e) and (i - 208) % 600 == 0
        noisy |= i >= 208 + len(reference_e) and (i - 208) % 31 == 0
        if noisy:
            sparse_read[i] = next(b for b in 'ACGT' if b not in
                                  (sparse_read[i-1], sparse_read[i], sparse_read[i+1]))
    sparse_read = ''.join(sparse_read)
    ambiguous_read = sparse_read[:208] + reference_e + sparse_read[208+len(reference_e):]
    # Six unique short UL reads and two multiply placed LA molecules on a
    # perfectly periodic template: all ordinary optimal locations are tied.
    unit_f = list(sequence[:732])
    unit_f[-1] = next(b for b in 'ACGT' if b not in (unit_f[-2], unit_f[0]))
    unit_f = ''.join(unit_f)
    reference_f = unit_f * 4
    repeat_read = unit_f * 2
    (root/'reads.fa').write_text(''.join(f'>LA{group}{i}\n{sequence}\n'
                                        for group in ('a', 'b', 'c') for i in range(6)) +
                                 '>LAnoise\n' + noise + '\n' +
                                 f'>LAd0\n{mutated}\n>ULd1\n{mutated}\n' +
                                 f'>LAe0\n{sparse_read}\n>LAe1\n{ambiguous_read}\n' +
                                 ''.join(f'>ULf{i}\n{repeat_read}\n' for i in range(6)) +
                                 ''.join(f'>LAf{i}\n{repeat_read}\n' for i in range(2)))
    (root/'template.fa').write_text(''.join(f'>{i}\n{reference}\n' for i in (1, 2, 3, 4)) +
                                    f'>5\n{reference_e}\n>6\n{reference_f}\n')
    layouts = []
    for group in ('a', 'b', 'c'):
        rows = []
        for i in range(6):
            begin = 160 if group == 'b' and i == 5 else 100
            placements = 2 if group == 'c' else 1
            rows.append(f'reads_LA{group}{i}\t{begin}\t{begin+12000}\t{placements}\n')
        if group == 'b':
            rows.append('reads_LAnoise\t100\t12100\t1\n')
        layouts.append(f'tig\tfixture{group}\nlen\t12160\nrds\t{len(rows)}\n' +
                       ''.join(rows) + 'end\n')
    layouts.append('tig\tfixtured\nlen\t12100\nrds\t2\n'
                   'reads_LAd0\t100\t12100\t1\n'
                   'reads_ULd1\t100\t12100\t1\nend\n')
    layouts.append('tig\tfixturee\nlen\t28686\nrds\t3\n'
                   'reads_LAe0\t100\t20686\t1\n'
                   'reads_LAe1\t4100\t24686\t2\n'
                   'reads_LAe1\t8100\t28686\t2\nend\n')
    repeat_rows = []
    for i in range(6):
        begin = 100 + (i // 2) * 732
        repeat_rows.append(f'reads_ULf{i}\t{begin}\t{begin+1464}\t1\n')
    for i in range(2):
        for begin in (2296, 3760):
            repeat_rows.append(f'reads_LAf{i}\t{begin}\t{begin+1464}\t2\n')
    layouts.append('tig\tfixturef\nlen\t5224\nrds\t10\n' + ''.join(repeat_rows) + 'end\n')
    (root/'layout').write_text(''.join(layouts))
    run([tools/'sqStoreCreate', '-o', root/'store', '-minlength', '1',
         '-raw', '-pacbio', 'reads', root/'reads.fa'])
    package_log = run([tools/'layoutToPackage', '-layout', root/'layout', '-store', root/'store',
         '-partition', '10', '1', '2', '-output', root/'part#.cnspack'])
    # Native partition numbering may also create an empty partition zero.
    packages = [p for p in root.glob('*.cnspack') if p.stat().st_size]
    assert len(packages) == 1, package_log
    package = packages[0]
    for label, flags in [('full', []), ('region', ['--region', 'tig00000002:0-1000']),
                         ('unique_only', ['--unique-only']),
                         ('sparse_tiers', ['--region', 'tig00000004:0-1000']),
                         ('zero_unique', ['--region', 'tig00000003:0-1000']),
                         ('sparse_unique', ['--region', 'tig00000005:0-1000']),
                         ('periodic', ['--region', 'tig00000006:0-1000'])]:
        log = run([exe, '--import', package, '--template-hpc', root/'template.fa',
                   '--tier', 'LA=reads_LA*', '--tier', 'UL=reads_UL*', '--no-extend-ends',
                   '--output', root/label, *flags])
        # layoutToPackage normalizes away the shared 100-base layout prefix.
        assert ('offset rescue: internal anchors, offset 3000' in log or
                'offset joint rescue: offset=3000' in log or
                'offset rescue: unique template anchors, offset 3000' in log or
                'offset rescue: sparse unique overlaps, offset 208' in log or
                'offset rescue: periodic sequence proposals, offset 732' in log), log
        observed = ''.join(s for s in (root/(label+'.fasta')).read_text().splitlines()
                           if not s.startswith('>'))
        expected = (reference * 4 + reference_e + reference_f if label in ('full', 'unique_only') else
                    reference_e[:1000] if label == 'sparse_unique' else
                    reference_f[:1000] if label == 'periodic' else reference[:1000])
        assert observed == expected
        if label == 'sparse_unique':
            assert 'offset rescue: sparse unique overlaps, offset 208' in log, log
        if label == 'periodic':
            assert 'offset rescue: periodic sequence proposals, offset 732' in log, log
            assert 'REPEAT_LAYOUT' in (root/(label+'.template-offsets.tsv')).read_text()
        if label in ('full', 'unique_only', 'zero_unique', 'sparse_tiers'):
            methods = (root/(label+'.template-offsets.tsv')).read_text()
            assert 'LAYOUT_OVERLAP' in methods or 'TEMPLATE_ANCHOR' in methods, methods
