#!/usr/bin/env python3
"""Small data-free checks for matched evidence and additive unique retention."""
import pathlib
import subprocess
import sys
import tempfile

with tempfile.TemporaryDirectory(prefix='hpc-graph-test-') as tmp:
    root = pathlib.Path(tmp)
    left, right = 'CGTCGTCGTCGTCGT', 'CGTCGTCGTCGTCGT'
    (root/'reference').write_text(left+'AAAAAA'+right+'\n')
    (root/'reads').write_text(''.join(
        f'{tier}{i}\t{tier}\t{left}{"A"*length}{right}\n'
        for tier, depth, length in [('LA', 3, 6), ('UL', 5, 4), ('LAO', 6, 6)]
        for i in range(depth)))
    result = subprocess.check_output([sys.argv[1], str(root/'reference'), str(root/'reads'),
                                      str(len(left)), str(len(left)+6)], text=True)
    selected = {}; graphs = {}
    for line in result.splitlines():
        f = line.split('\t')
        if f[0] == 'READ':
            assert int(f[3]) == (4 if f[2] == 'UL' else 6)
            assert f[4] == 'A'*int(f[3])
        elif f[0] == 'SELECT': selected.setdefault(f[1], set()).add(f[2])
        elif f[0] == 'GRAPH': graphs[f[1]] = f[2]
    assert len(selected['UL_exclusive']) == 5
    assert len(selected['UL_additive']) == 8
    assert {'LA0', 'LA1', 'LA2'} <= selected['UL_additive']
    assert len(selected['LAO_additive']) == 9
    assert len(selected['all_unique']) == 14
    assert graphs['UL_exclusive'].startswith(left+'AAAA'+right[:-1])
    assert graphs['LAO_additive'].startswith(left+'AAAAAA'+right[:-1])
    # With LAO absent, LA/UL still uses the first qualifying tier. If neither
    # qualifies, preserve sparse LA rather than silently choosing the last tier.
    for la_depth, ul_depth, expected in [(3, 5, 'UL'), (3, 2, 'LA')]:
        (root/'reads').write_text(''.join(
            f'{tier}{i}\t{tier}\t{left}{"A"*length}{right}\n'
            for tier, depth, length in [('LA', la_depth, 6), ('UL', ul_depth, 4)]
            for i in range(depth)))
        sparse = subprocess.check_output([sys.argv[1], str(root/'reference'), str(root/'reads'),
                                          str(len(left)), str(len(left)+6)], text=True)
        names = [line.split('\t')[2] for line in sparse.splitlines()
                 if line.startswith('SELECT\tLAO_exclusive\t')]
        assert names and all(name.startswith(expected) for name in names)
    with (root/'reads').open('a') as out:
        out.write('bad\tLA\t'+'N'*len(left+'AAAAAA'+right)+'\n')
    rejected = subprocess.check_output([sys.argv[1], str(root/'reference'), str(root/'reads'),
                                       str(len(left)), str(len(left)+6)], text=True)
    assert 'REJECT\tbad\tAlignment exceeds maximum error rate' in rejected
    assert 'SELECT\tLAO_exclusive\tLA0' in rejected
print('Matched-read graph and additive-selection checks passed')
