"""Native regression: all reads overhang both ends of a shorter HPC template."""
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
    return p.stderr


with tempfile.TemporaryDirectory(prefix='hpc-offset-') as temporary:
    root = pathlib.Path(temporary)
    rng = random.Random(913)
    bases = ['A']
    for _ in range(11999):
        bases.append(rng.choice([b for b in 'ACGT' if b != bases[-1]]))
    sequence = ''.join(bases)
    reference = sequence[3000:9000]
    (root/'reads.fa').write_text(''.join(f'>LA{i}\n{sequence}\n' for i in range(6)))
    (root/'template.fa').write_text('>1\n' + reference + '\n')
    (root/'layout').write_text('tig\tfixture\nlen\t12100\nrds\t6\n' +
                              ''.join(f'reads_LA{i}\t100\t12100\t1\n'
                                      for i in range(6)) + 'end\n')
    run([tools/'sqStoreCreate', '-o', root/'store', '-minlength', '1',
         '-raw', '-pacbio', 'reads', root/'reads.fa'])
    run([tools/'layoutToPackage', '-layout', root/'layout', '-store', root/'store',
         '-output', root/'part#.cnspack'])
    package = next(root.glob('*.cnspack'))
    for label, flags in [('full', []), ('region', ['--region', 'tig00000001:0-1000'])]:
        log = run([exe, '--import', package, '--template-hpc', root/'template.fa',
                   '--tier', 'LA=reads_LA*', '--no-extend-ends',
                   '--output', root/label, *flags])
        # layoutToPackage normalizes away the shared 100-base layout prefix.
        assert 'offset rescue: internal anchors, offset 3000' in log, log
        observed = ''.join(s for s in (root/(label+'.fasta')).read_text().splitlines()
                           if not s.startswith('>'))
        assert observed == (reference if label == 'full' else reference[:1000])
