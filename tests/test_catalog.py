import unittest
from cnv_catalog import merge_depth_grid


def snp(chrom, pos):
    return dict(chromosome=chrom, position=pos, ref='A', alt='C', kind='SNP')


class CatalogTests(unittest.TestCase):
    def test_native_grid_boundaries_and_snp_precedence(self):
        rows = list(merge_depth_grid([snp('chr1', 1), snp('chr1', 10), snp('chr1', 21)],
                                    {'chr1': 21}, 10))
        self.assertEqual([row['position'] for row in rows], [1, 10, 20, 21])
        self.assertEqual([row['kind'] for row in rows], ['SNP', 'SNP', 'DEPTH', 'SNP'])

    def test_reference_order_not_lexical_and_short_contigs(self):
        rows = list(merge_depth_grid([snp('chr2', 4), snp('chr10', 1)],
                                    {'chr2': 5, 'chr10': 20}, 10))
        self.assertEqual([(row['chromosome'], row['position']) for row in rows],
                         [('chr2', 4), ('chr10', 1), ('chr10', 10), ('chr10', 20)])

    def test_no_observed_snps_does_not_remove_grid(self):
        rows = list(merge_depth_grid(iter(()), {'chr1': 20}, 10))
        self.assertEqual([row['position'] for row in rows], [10, 20])

    def test_duplicate_unsorted_out_of_bounds_rejected(self):
        for rows in ([snp('chr1', 1), snp('chr1', 1)],
                     [snp('chr1', 3), snp('chr1', 2)], [snp('chr1', 11)]):
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                list(merge_depth_grid(rows, {'chr1': 10}, 10))

    def test_stream_does_not_materialize_entire_reference(self):
        iterator = merge_depth_grid([], {'chr1': 300000000}, 10)
        self.assertEqual(next(iterator)['position'], 10)
        self.assertEqual(next(iterator)['position'], 20)
        iterator.close()


if __name__ == '__main__':
    unittest.main()
