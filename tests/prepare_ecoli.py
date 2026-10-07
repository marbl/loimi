"""Recover layout origins by observing the original deterministic layout generator."""
import pathlib
import runpy
import sys

asm, verkko, output = map(pathlib.Path, sys.argv[1:])
output.mkdir(parents=True, exist_ok=True)
script = verkko / 'src/scripts/get_layout_from_mbg.py'
layout = asm / '6-layoutContigs'
origins = {}
def trace(frame, event, arg):
    if frame.f_code.co_filename == str(script) and event == 'line' and frame.f_lineno == 620:
        v = frame.f_locals
        origins[v['contig']] = (v['start_pos'], v['end_pos'])
    return trace
sys.path.insert(0, str(script.parent))
sys.argv = [str(script)] + [str(layout / f) for f in ['combined-nodemap.txt', 'combined-edges.gfa', 'combined-alignments.gaf', 'consensus_paths.txt', 'nodelens.txt']] + [str(output / 'regenerated.layout'), str(output / 'regenerated.scfmap')]
sys.settrace(trace)
runpy.run_path(str(script), run_name='__main__')
sys.settrace(None)
with (output / 'origins.tsv').open('w') as out:
    for name, (start, end) in origins.items():
        out.write(f'{name}\t{start}\t{end}\n')
with (asm / 'assembly.homopolymer-compressed.gfa').open() as inp, (output / 'template.fa').open('w') as out:
    for line in inp:
        fields = line.rstrip().split('\t')
        if fields[0] == 'S' and fields[1] == 'utig4-0':
            out.write('>piece000001\n' + fields[2] + '\n')
print(origins)
