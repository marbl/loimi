"""Native regression: censored observations cannot qualify a higher tier."""
import csv
import pathlib
import subprocess
import sys
import tempfile

exe = pathlib.Path(sys.argv[1])
tools = pathlib.Path(sys.argv[2])/'build/bin'


def run(args):
    result = subprocess.run(list(map(str, args)), capture_output=True, text=True)
    assert result.returncode == 0, result.stdout+result.stderr


with tempfile.TemporaryDirectory(prefix='hpc-anchor-') as temp:
    for side in ('left', 'right'):
        root = pathlib.Path(temp)/side
        root.mkdir()
        reads = [(f'LA{i}', 'ACGTGGGACT', 0, 8) for i in range(4)]
        reads += [('LAterminal', 'GACT', 4, 8)] if side == 'left' else [('LAterminal', 'ACGTG', 0, 5)]
        reads += [(f'LOW{i}', 'ACGTGGACT', 0, 8) for i in range(5)]
        (root/'reads.fa').write_text(''.join(f'>{n}\n{s}\n' for n,s,_,_ in reads))
        (root/'template.fa').write_text('>1\nACGTGACT\n')
        (root/'sites.tsv').write_text('tig00000001\t4\n')
        (root/'layout').write_text('tig\tfixture\nlen\t8\nrds\t10\n'+''.join(
            f'reads_{n}\t{b}\t{e}\t1\n' for n,_,b,e in reads)+'end\n')
        run([tools/'sqStoreCreate','-o',root/'store','-minlength','1','-raw','-pacbio','reads',root/'reads.fa'])
        run([tools/'layoutToPackage','-layout',root/'layout','-store',root/'store','-output',root/'part#.cnspack'])
        for label, flags in [('baseline', []), ('anchored', ['--experiment-terminal-hp-filter']),
                             ('span', ['--experiment-hp-flank-bases', '1'])]:
            run([exe,'--import',next(root.glob('*.cnspack')),'--template-hpc',root/'template.fa',
                 '--tier','LA=reads_LA*','--tier','LOW=reads_LOW*','--output',root/label,
                 '--evidence-sites',root/'sites.tsv','--evidence-read-votes',*flags])
        def sequence(label):
            return ''.join(s for s in (root/(label+'.fasta')).read_text().splitlines() if not s.startswith('>'))
        assert sequence('baseline') == 'ACGTGGGACT'
        assert sequence('anchored') == 'ACGTGGACT'
        assert sequence('span') == 'ACGTGGACT'
        with (root/'anchored.run-read-votes.tsv').open() as source:
            rows = list(csv.DictReader(source, delimiter='\t'))
        assert next(r for r in rows if r['name']=='reads_LAterminal')['status']=='UNFLANKED_RUN'
        with (root/'anchored.run-calls.tsv').open() as source:
            call = next(r for r in csv.DictReader(source, delimiter='\t') if r['run']=='4')
        assert call['tier']=='LOW' and call['depth']=='5'
