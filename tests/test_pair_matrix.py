from pathlib import Path
import tempfile
import unittest
from cnv_pair_matrix import pair_matrix,native_rcmat_rows
from cnv_pileup import COLUMNS
from cnv_territory import Territory
from test_contract import fixture


class PairMatrixTests(unittest.TestCase):
    def setUp(self):
        temp=tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        root=Path(temp.name)
        self.paths={'study':root/'study.csv','pon':root/'pon.csv'}
        header=','.join(COLUMNS)+'\n'
        self.paths['study'].write_text(header+'chr1,1,A,C,0,0,5,0\nchr1,2,.,.,7,0,0,0\n')
        self.paths['pon'].write_text(header+'chr1,1,A,C,10,10,0,0\nchr1,3,G,T,10,5,0,0\n')
        self.manifest=fixture()
        for row in self.manifest['samples']:
            row['reference_sex']='Female'
        self.loci=[dict(chromosome='chr1',position=p,ref=r,alt=a,kind=k)
                   for p,r,a,k in [(1,'A','C','SNP'),(2,'.','.','DEPTH'),(3,'G','T','SNP')]]
        self.mask=Territory([('chr1',0,1)],{'chr1':10})

    def rows(self,reference='pon'):
        return list(pair_matrix(self.manifest,'study',reference,self.loci,self.paths,
                                {'chr1':10},mask=self.mask))

    def test_native_zero_and_absent_rows_remain_different(self):
        rows=self.rows()
        self.assertEqual(rows[0]['TUM.DP'],0)
        self.assertEqual(rows[0]['TUM.ERROR'],5)
        self.assertTrue(rows[0]['native_input_eligible'])
        self.assertTrue(rows[0]['masked'])
        self.assertIsNone(rows[1]['NOR.DP'])
        self.assertIn('REFERENCE_PILEUP_ROW_ABSENT',rows[1]['reasons'])
        self.assertIsNone(rows[2]['TUM.DP'])
        self.assertEqual(rows[2]['NOR.DP'],15)
        native=list(native_rcmat_rows(iter(rows)))
        self.assertEqual(len(native),1)
        self.assertEqual(native[0]['catalog_row'],1)
        self.assertEqual(native[0]['TUM.DP'],0)

    def test_no_self_or_calibration_reference(self):
        for reference in ('study','calibration'):
            with self.subTest(reference=reference),self.assertRaisesRegex(ValueError,'independent'):
                self.rows(reference)

    def test_reference_sex_must_match_and_be_resolved(self):
        for sex in ('Unknown','Male'):
            self.manifest['samples'][1]['reference_sex']=sex
            with self.subTest(sex=sex),self.assertRaisesRegex(ValueError,'sex'):
                self.rows()

    def test_foreign_inputs_rejected(self):
        self.paths['foreign']=self.paths['study']
        with self.assertRaisesRegex(ValueError,'explicit'):
            self.rows()
