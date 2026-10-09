# SPDX-License-Identifier: MIT-0
"""Synthetic train/held-out separation and frozen model provenance checks."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

script = Path(__file__).resolve().parents[1] / 'scripts/calibrate_runs.py'

with tempfile.TemporaryDirectory(prefix='hpc-model-test-') as td:
    root = Path(td)
    (root / 'manifest.json').write_text(json.dumps({'windows': [
        {'name': 'train', 'split': 'train'}, {'name': 'heldout', 'split': 'validation'}]}))
    for name in ('train', 'heldout'):
        folder = root / name
        folder.mkdir()
        (folder / 'consensus.summary.json').write_text('{}')
        (folder / 'label-summary.json').write_text('{}')
        labels = ['run\tbase\ttruth_length']
        evidence = ['tig\trun\tbase\ttier\tplacement_class\tlength\tcount']
        for run in range(120):
            length = run // 40 + 1
            labels.append(f'{run}\tA\t{length}')
            for tier in ('LA', 'LAO', 'UL'):
                evidence.append(f'tig00000001\t{run}\tA\t{tier}\tUNIQUE\t{length}\t10')
        labels.append('999\tA\t3')  # Missing evidence must remain in pipeline denominators.
        evidence.append('tig00000001\t998\tA\tLA\tUNIQUE\t50\t100')  # No truth label: exclude.
        (folder / 'labels.tsv').write_text('\n'.join(labels) + '\n')
        (folder / 'consensus.run-evidence.tsv').write_text('\n'.join(evidence) + '\n')
    heldout = root / 'heldout/consensus.run-evidence.tsv'
    valid = heldout.read_bytes()
    heldout.write_text('This validation evidence must not be read during fitting.\n')
    command = [sys.executable, str(script), '--data', str(root)]
    subprocess.run(command + ['train'], check=True, capture_output=True)
    metadata = json.loads((root / 'models/metadata.json').read_text())
    assert metadata['training_windows'] == ['train']
    assert metadata['models']['LA']['supported_lengths_1_to']['A'] == 3
    heldout.write_bytes(valid)
    subprocess.run(command + ['validate'], check=True, capture_output=True)
    result = json.loads((root / 'validation.json').read_text())
    assert result['windows']['heldout']['sampled_no_evidence_loci'] == 1
    assert result['priority_pipelines']['LA-first']['loci'] == 121
    assert result['priority_pipelines']['LA-first']['model_wrong'] == 1
    for tier in ('LA', 'LAO', 'UL'):
        row = result['by_tier'][tier + '/UNIQUE/all']
        assert row['loci'] == 120 and row['model_wrong'] == 0
    assert subprocess.run(command + ['train'], capture_output=True).returncode != 0
    with (root / 'models/LA.run-model.tsv').open('a') as out:
        out.write('# changed after freezing\n')
    assert subprocess.run(command + ['validate'], capture_output=True).returncode != 0
