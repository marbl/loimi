"""Native package regression: ambiguous evidence changes an uncertain run only."""
import pathlib
import csv
import collections
import subprocess
import sys
import tempfile

exe=pathlib.Path(sys.argv[1]);tools=pathlib.Path(sys.argv[2])/'build/bin'
def run(args):
    p=subprocess.run(list(map(str,args)),capture_output=True,text=True)
    assert p.returncode==0,p.stdout+p.stderr
with tempfile.TemporaryDirectory(prefix='hpc-borrow-') as tmp:
    root=pathlib.Path(tmp)
    reads=[(f'u{i}',2 if i<6 else 3,1) for i in range(10)]
    reads += [(f'a{i}',3,2) for i in range(20)]
    (root/'reads.fa').write_text(''.join(f'>{n}\nACGT'+ 'G'*length+'ACT\n' for n,length,_ in reads))
    (root/'template.fa').write_text('>1\nACGTGACT\n')
    (root/'layout').write_text('tig\tfixture\nlen\t8\nrds\t30\n'+''.join(f'reads_{n}\t0\t8\t{copies}\n' for n,_,copies in reads)+'end\n')
    run([tools/'sqStoreCreate','-o',root/'store','-minlength','1','-raw','-pacbio','reads',root/'reads.fa'])
    run([tools/'layoutToPackage','-layout',root/'layout','-store',root/'store','-output',root/'part#.cnspack'])
    for label,extra in [('baseline',[]),('borrow',['--experiment-consistent-ambiguous'])]:
        run([exe,'--import',next(root.glob('*.cnspack')),'--template-hpc',root/'template.fa',
             '--tier','LA=reads_*','--output',root/label,*extra])
    seq=lambda name: ''.join(x for x in (root/(name+'.fasta')).read_text().splitlines() if not x.startswith('>'))
    assert seq('baseline')=='ACGTGGACT'
    assert seq('borrow')=='ACGTGGGACT'
    assert 'UNIQUE_BORROWED' in (root/'borrow.run-calls.tsv').read_text()
    assert '\t2\t3\n' in (root/'borrow.borrowing.tsv').read_text()
    (root/'sites.tsv').write_text('tig00000001\t4\n')
    run([exe,'--import',next(root.glob('*.cnspack')),'--template-hpc',root/'template.fa',
         '--tier','LA=reads_*','--output',root/'trace','--evidence-sites',root/'sites.tsv',
         '--evidence-read-votes'])
    assert seq('trace') == seq('baseline')
    with (root/'trace.run-read-votes.tsv').open() as source:
        rows=list(csv.DictReader(source,delimiter='\t'))
    assert len(rows)==30 and all(r['status']=='ELIGIBLE' for r in rows)
    assert collections.Counter(int(r['length']) for r in rows if r['placements']=='1')=={2:6,3:4}
    assert all(int(r['query_raw_end'])-int(r['query_raw_begin'])==int(r['length']) for r in rows)
