"""Native integration: corroboration changes only uncertain first-tier calls.

Also exercises HP realignment with unchanged controls and optional final CRAM.
Synthetic files live only in a temporary directory.
"""
import pathlib
import random
import shutil
import subprocess
import sys
import tempfile

exe = pathlib.Path(sys.argv[1])
tools = pathlib.Path(sys.argv[2])/'build/bin'


def run(args):
    result = subprocess.run(list(map(str, args)), capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


with tempfile.TemporaryDirectory(prefix='hpc-refinement-') as temp:
    root = pathlib.Path(temp)
    rng = random.Random(17)
    def flank(n):
        s = 'C'
        while len(s) < n:
            s += rng.choice([b for b in 'ACGT' if b != s[-1]])
        return s
    left, right = flank(600), 'T' + flank(600)
    if left[-1] == 'G':
        left += 'C'
    template = left + 'G' + right
    reads = [(f'LA{i}', 2 if i < 3 else 3) for i in range(5)]
    reads += [(f'LOW{i}', 3) for i in range(5)]
    (root/'reads.fa').write_text(''.join(f'>{n}\n{left}{"G"*length}{right}\n' for n,length in reads))
    (root/'template.fa').write_text(f'>1\n{template}\n')
    (root/'layout').write_text(f'tig\tfixture\nlen\t{len(template)}\nrds\t10\n'+''.join(
        f'reads_{n}\t0\t{len(template)}\t1\n' for n,_ in reads)+'end\n')
    for name, p, q in [('la', '.55', '.42'), ('low', '.80', '.17')]:
        (root/(name+'.model')).write_text('HPC_RUN_MODEL_V1 4 5\nG 1 1 .01 .96 .01 .01 .01\n'
            f'G 2 1 .01 .01 {p} {q} .01\nG 3 1 .01 .01 {q} {p} .01\n')
    run([tools/'sqStoreCreate','-o',root/'store','-minlength','1','-raw','-pacbio','reads',root/'reads.fa'])
    run([tools/'layoutToPackage','-layout',root/'layout','-store',root/'store','-output',root/'part#.cnspack'])
    for arm, flags in [('baseline', []), ('corroboration', ['--experiment-corroboration-model', f'LOW={root}/low.model']),
                       ('realign', ['--experiment-hp-realign']),
                       ('combined', ['--experiment-hp-realign', '--experiment-corroboration-model', f'LOW={root}/low.model']),
                       ('flank500', ['--experiment-hp-flank-bases', '500']),
                       ('all-three', ['--experiment-hp-realign', '--experiment-corroboration-model', f'LOW={root}/low.model',
                                      '--experiment-hp-flank-bases', '500']),
                       ('all-three-linked', ['--experiment-hp-realign', '--experiment-corroboration-model', f'LOW={root}/low.model',
                                      '--experiment-hp-flank-bases', '500'])]:
        cram = ['--cram', root/(arm+'.cram')] if shutil.which('samtools') else []
        grouping = [] if arm == 'all-three-linked' else ['--experiment-master-events']
        run([exe,'--import',next(root.glob('*.cnspack')),'--template-hpc',root/'template.fa',
             '--tier','LA=reads_LA*','--tier','LOW=reads_LOW*','--tier-model', f'LA={root}/la.model',
             '--polish-non-hpc','--no-extend-ends',*grouping,
             '--output',root/arm,*flags,*cram])
        sequence = ''.join(s for s in (root/(arm+'.fasta')).read_text().splitlines() if not s.startswith('>'))
        assert sequence == left + ('GGG' if arm in ('corroboration','combined','all-three','all-three-linked') else 'GG') + right, (arm, sequence)
        if cram:
            run(['samtools','quickcheck',root/(arm+'.cram')])
    assert '\t2\t3\n' in (root/'corroboration.corroboration.tsv').read_text()
    assert len((root/'realign.hp-realignment.tsv').read_text().splitlines()) == 1

    # Actual post-draft HP edit: uncertain high-tier votes lack the 64-base
    # realignment context, so complete lower-tier UNIQUE observations qualify.
    # This exercises changed-reference CRAM realignment, not only no-op controls.
    short = root/'short'
    read_start = len(left) - 20
    short.mkdir()
    (short/'reads.fa').write_text(''.join(
        f'>{n}\n{left[read_start:] if n.startswith("LA") else left}{"G"*length}{right}\n'
        for n,length in reads))
    (short/'layout').write_text(f'tig\tfixture\nlen\t{len(template)}\nrds\t10\n'+''.join(
        f'reads_{n}\t{read_start if n.startswith("LA") else 0}\t{len(template)}\t1\n'
        for n,_ in reads)+'end\n')
    run([tools/'sqStoreCreate','-o',short/'store','-minlength','1','-raw','-pacbio','reads',short/'reads.fa'])
    run([tools/'layoutToPackage','-layout',short/'layout','-store',short/'store','-output',short/'part#.cnspack'])
    cram = ['--cram', short/'changed.cram'] if shutil.which('samtools') else []
    run([exe,'--import',next(short.glob('*.cnspack')),'--template-hpc',root/'template.fa',
         '--tier','LA=reads_LA*','--tier','LOW=reads_LOW*','--polish-non-hpc','--no-extend-ends',
         '--experiment-hp-realign','--experiment-master-events','--output',short/'changed',*cram])
    seq = lambda suffix: ''.join(s for s in (short/('changed'+suffix)).read_text().splitlines()
                                if not s.startswith('>'))
    assert seq('.draft.fasta') == left+'GG'+right
    assert seq('.fasta') == left+'GGG'+right
    assert '\t2\t3\tLOW\t5\n' in (short/'changed.hp-realignment.tsv').read_text()
    if cram:
        run(['samtools','quickcheck',short/'changed.cram'])
        records = subprocess.check_output(['samtools','view',str(short/'changed.cram')], text=True)
        assert len(records.splitlines()) == 10
