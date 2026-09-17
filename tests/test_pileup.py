import gzip
from io import StringIO
from pathlib import Path
import tempfile
import unittest

from cnv_pileup import COLUMNS, align_catalog, native_rows, training_summaries
from test_allelic import policy


def catalog():
    return [dict(chromosome='chr1', position=i, ref=ref, alt=alt, kind=kind)
            for i, ref, alt, kind in [(1, 'A', 'C', 'SNP'), (2, '.', '.', 'DEPTH'),
                                      (3, 'G', 'T', 'SNP')]]


class PileupTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.paths = {'a': self.root/'a.csv', 'b': self.root/'b.csv.gz'}
        self.paths['a'].write_text(','.join(COLUMNS)+'\nchr1,1,A,C,0,0,5,0\nchr1,2,.,.,4,0,0,0\n')
        with gzip.open(self.paths['b'], 'wt') as stream:
            stream.write(','.join(COLUMNS)+'\nchr1,1,A,C,10,10,0,1\n')

    def test_missing_zero_depth_and_pseudo_are_separate(self):
        rows = list(training_summaries(iter(catalog()), self.paths, {'chr1': 10},
                                       ['a', 'b'], policy()))
        self.assertEqual(rows[0]['zero_depth_normals'], 1)
        self.assertEqual(rows[0]['missing_normals'], 0)
        self.assertEqual(rows[0]['native_error_count'], 5)
        self.assertEqual(rows[0]['native_deletion_count'], 1)
        self.assertEqual(rows[1]['status'], 'DEPTH_ONLY')
        self.assertEqual(rows[1]['missing_normals'], 1)
        self.assertIsNone(rows[1]['correction_bias'])
        self.assertEqual(rows[2]['missing_normals'], 2)
        self.assertIsNone(rows[2]['median_observed_depth'])
        self.assertIsNone(rows[2]['native_error_count'])

    def test_alleles_cannot_be_silently_swapped(self):
        loci = catalog()
        loci[0].update(ref='C', alt='A')
        with self.assertRaisesRegex(ValueError, 'allele mismatch'):
            list(align_catalog(loci, self.paths, {'chr1': 10}))

    def test_extra_positions_do_not_disappear(self):
        with self.assertRaisesRegex(ValueError, 'trailing'):
            list(align_catalog(catalog()[:1], self.paths, {'chr1': 10}))

    def test_rejects_duplicates_invalid_counts_and_multisample_header(self):
        header = ','.join(COLUMNS)+'\n'
        for body in ('chr1,1,A,C,1,2,0,0\nchr1,1,A,C,1,2,0,0\n',
                     'chr1,1,A,C,-1,2,0,0\n', 'chr1,2,.,.,4,1,0,0\n'):
            with self.subTest(body=body), self.assertRaises(ValueError):
                list(native_rows(StringIO(header+body), {'chr1': 10}))
        with self.assertRaisesRegex(ValueError, 'schema'):
            list(native_rows(StringIO(header.strip()+',File2R\n'), {'chr1': 10}))

    def test_rejects_nonpon_and_reused_paths(self):
        with self.assertRaisesRegex(ValueError, 'membership'):
            list(training_summaries(catalog(), self.paths, {'chr1': 10}, ['a'], policy()))
        with self.assertRaisesRegex(ValueError, 'reused'):
            list(align_catalog(catalog(), {'a': self.paths['a'], 'b': self.paths['a']}, {'chr1': 10}))


if __name__ == '__main__':
    unittest.main()
