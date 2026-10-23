"""Data-free coordinate and repeat-context checks."""
import unittest
from classify_remaining import contains, run_bounds, tandem_period


class ContextTests(unittest.TestCase):
    def test_runs_and_edges(self):
        self.assertEqual(run_bounds('AACTTT', 0), (0, 2))
        self.assertEqual(run_bounds('AACTTT', 5), (3, 6))
        with self.assertRaises(ValueError):
            run_bounds('A', 1)

    def test_period(self):
        self.assertEqual(tandem_period('GATATATC', 3), 2)
        self.assertIsNone(tandem_period('GAAAAAAC', 3))
        self.assertIsNone(tandem_period('ACGT', 2))

    def test_half_open(self):
        self.assertTrue(contains([(2, 5)], [2], 2))
        self.assertFalse(contains([(2, 5)], [2], 5))


if __name__ == '__main__':
    unittest.main()
