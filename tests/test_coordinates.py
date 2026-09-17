from pathlib import Path
import tempfile
import unittest
from cnv_coordinates import bed_intervals, fai_lengths


class CoordinateTests(unittest.TestCase):
    def test_zero_based_territory_and_contig_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed = Path(tmp)/'targets.bed'
            bed.write_text('chr1\t0\t10\nchr1\t10\t20\n')
            self.assertEqual(sum(e-s for _,s,e in bed_intervals(bed, {'chr1': 20})), 20)
            for content in ('1\t0\t10\n', 'chr1\t0\t21\n', 'chr1\t10\t10\n',
                            'chr1\t10\t20\nchr1\t0\t10\n', 'chr1\t0\t15\nchr1\t10\t20\n'):
                bed.write_text(content)
                with self.subTest(content=content), self.assertRaises(ValueError):
                    bed_intervals(bed, {'chr1': 20})

    def test_fai_rejects_duplicate_reference_contig(self):
        with tempfile.TemporaryDirectory() as tmp:
            fai = Path(tmp)/'reference.fa.fai'
            fai.write_text('chr1\t20\t6\t20\t21\n')
            self.assertEqual(fai_lengths(fai), {'chr1': 20})
            fai.write_text(fai.read_text()*2)
            with self.assertRaises(ValueError):
                fai_lengths(fai)


if __name__ == '__main__':
    unittest.main()
