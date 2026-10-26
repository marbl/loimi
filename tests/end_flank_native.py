"""500-bp interior protection must coexist with template-boundary extension.

Real native packages, both orientations and ends; short interior read support
must stay excluded even though template-clipped overhangs supply end context.
"""
import csv
import json
import pathlib
import random
import shutil
import subprocess
import sys
import tempfile

exe = pathlib.Path(sys.argv[1])
tools = pathlib.Path(sys.argv[2])/'build/bin'


def run(args):
    p = subprocess.run(list(map(str,args)),capture_output=True,text=True)
    assert p.returncode == 0, p.stdout+p.stderr
    return p


def sequence(path):
    return ''.join(x for x in path.read_text().splitlines() if not x.startswith('>'))


with tempfile.TemporaryDirectory(prefix='hpc-end-flank-') as tmp:
    root = pathlib.Path(tmp)
    rng = random.Random(89)
    bases = 'A'
    for _ in range(1799):
        bases += rng.choice([b for b in 'ACGT' if b != bases[-1]])
    left, core, right = bases[:300], bases[300:1500], bases[1500:]
    double = lambda s: ''.join(b*2 for b in s)
    reads = [(f'anchor{i}',double(core),1000,2200) for i in range(6)]
    reads += [(f'left{i}',double(left+core),700,2200) for i in range(6)]
    reads += [(f'right{i}',double(core+right),1000,2500) for i in range(6)]
    # Interior short reads would incorrectly lengthen a run if the fix waived
    # flank filtering globally instead of accounting only for boundary clips.
    short = double(core[500:600])+core[600]*3+double(core[601:700])
    reads += [(f'short{i}',short,1500,1700) for i in range(6)]
    reads = [(n,s.translate(str.maketrans('ACGT','TGCA'))[::-1],e,b)
             if n in ('left1','left3','right1','right3') else (n,s,b,e) for n,s,b,e in reads]
    (root/'reads.fa').write_text(''.join(f'>{n}\n{s}\n' for n,s,_,_ in reads))
    (root/'template.fa').write_text('>1\n'+core+'\n')
    (root/'layout').write_text('tig\tends\nlen\t2500\nrds\t24\n'+''.join(
        f'ends_{n}\t{b}\t{e}\t1\n' for n,s,b,e in reads)+'end\n')
    run([tools/'sqStoreCreate','-o',root/'store','-minlength','1','-raw','-pacbio','ends',root/'reads.fa'])
    run([tools/'layoutToPackage','-layout',root/'layout','-store',root/'store','-output',root/'part#.cnspack'])
    args = [exe,'--import',next(root.glob('*.cnspack')),'--template-hpc',root/'template.fa',
            '--tier','short=ends_short*','--tier','anchors=ends_anchor*','--tier','tails=*',
            '--polish-non-hpc','--experiment-hp-realign','--experiment-hp-flank-bases','500']
    for arm, extra, expected in [
        ('disabled',['--no-extend-ends'],double(core)),
        ('both',[],double(left[1:]+core+right[:-1])),
        ('left',['--end-extension-side','left'],double(left[1:]+core)),
        ('right',['--end-extension-side','right'],double(core+right[:-1]))]:
        cram = ['--cram',root/(arm+'.cram')] if arm=='both' and shutil.which('samtools') else []
        run([*args,*extra,*cram,'--output',root/arm])
        assert sequence(root/(arm+'.fasta')) == expected, arm
        with (root/(arm+'.run-calls.tsv')).open() as stream:
            calls = list(csv.DictReader(stream,delimiter='\t'))
        assert all(int(c['length'])==2 and c['tier']!='NO_EVIDENCE' for c in calls), arm
        assert calls[600]['tier']=='anchors' and calls[600]['depth']=='6', arm
        if arm=='both':
            summary=json.loads((root/(arm+'.summary.json')).read_text())
            assert summary['left_extension_bases']==598 and summary['right_extension_bases']==598
        if cram:
            run(['samtools','quickcheck',root/(arm+'.cram')])
            result=run(['samtools','calmd',root/(arm+'.cram'),root/(arm+'.fasta')])
            assert 'different' not in result.stderr
