# SPDX-License-Identifier: MIT-0
"""Data-free contracts for offline evidence selection and scoring."""
import collections
import unittest
import tempfile
from pathlib import Path
from unittest.mock import patch
from corroboration import select, record, bucket
from local_anchor_analysis import selected_sites


class ReplayContracts(unittest.TestCase):
    def test_all_labels_selection(self):
        labels = {1: ('A', 3), 2: ('G', 12), 3: ('T', 4)}
        with tempfile.TemporaryDirectory() as td, patch('local_anchor_analysis.c.labels', return_value=labels):
            returned, sites, changed = selected_sites(Path(td), 'w01', all_labels=True)
        self.assertEqual(returned, labels)
        self.assertEqual(sites, set(labels))
        self.assertEqual(changed, set())
    def test_unique_tier_precedes_pooled(self):
        h = {('LA', 'UNIQUE'): {10: 2}, ('LA', 'AMBIGUOUS'): {11: 100},
             ('LAO', 'UNIQUE'): {10: 5}}
        tier, kind, hist = select(h, ('LA', 'LAO', 'UL'))
        self.assertEqual((tier, kind, dict(hist)), ('LAO', 'UNIQUE', {10: 5}))

    def test_pooled_fallback(self):
        tier, kind, hist = select({('LA', 'UNIQUE'): {10: 2}, ('LA', 'AMBIGUOUS'): {10: 3}}, ('LA', 'UL'))
        self.assertEqual((tier, kind, dict(hist)), ('LA', 'POOLED', {10: 5}))

    def test_sparse_and_absent(self):
        self.assertEqual(select({}, ('LA', 'UL')), (None, 'NONE', {}))
        self.assertEqual(select({('UL', 'UNIQUE'): {4: 1}}, ('LA', 'UL'))[:2], ('UL', 'SPARSE'))

    def test_fixes_breaks_and_sampling_weights(self):
        s = collections.Counter()
        record(s, 12, 11, 12)
        record(s, 2, 2, 3)
        self.assertEqual((s['fixes'], s['breaks'], s['wrong'], s['weighted_wrong']), (1, 1, 1, 100))
        self.assertEqual(bucket(3, (1, 3, 10)), '10')


if __name__ == '__main__':
    unittest.main()
