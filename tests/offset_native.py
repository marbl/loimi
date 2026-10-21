"""Multi-tig regression: overhangs, small samples, and zero unique reads."""
import pathlib
import random
import subprocess
import sys
import tempfile

exe = pathlib.Path(sys.argv[1])
tools = pathlib.Path(sys.argv[2]) / 'build/bin'


def run(args):
    p = subprocess.run(list(map(str, args)), capture_output=True, text=True)
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
    (root/'reads.fa').write_text(''.join(f'>LA{group}{i}\n{sequence}\n'
                                        for group in ('a', 'b', 'c') for i in range(6)) +
                                 '>LAnoise\n' + noise + '\n' +
                                 f'>LAd0\n{mutated}\n>ULd1\n{mutated}\n')
    (root/'template.fa').write_text(''.join(f'>{i}\n{reference}\n' for i in (1, 2, 3, 4)))
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
                         ('zero_unique', ['--region', 'tig00000003:0-1000'])]:
        log = run([exe, '--import', package, '--template-hpc', root/'template.fa',
                   '--tier', 'LA=reads_LA*', '--tier', 'UL=reads_UL*', '--no-extend-ends',
                   '--output', root/label, *flags])
        # layoutToPackage normalizes away the shared 100-base layout prefix.
        assert ('offset rescue: internal anchors, offset 3000' in log or
                'offset joint rescue: offset=3000' in log), log
        observed = ''.join(s for s in (root/(label+'.fasta')).read_text().splitlines()
                           if not s.startswith('>'))
        assert observed == (reference * 4 if label in ('full', 'unique_only') else reference[:1000])
        if label in ('full', 'unique_only', 'zero_unique', 'sparse_tiers'):
            assert 'LAYOUT_OVERLAP' in (root/(label+'.template-offsets.tsv')).read_text()
