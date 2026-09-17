from pathlib import Path
import tempfile
import unittest
from cnv_manta import manta_evidence,support_pair


class MantaEvidenceTests(unittest.TestCase):
    def setUp(self):
        temp=tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.path=Path(temp.name)/'tumorSV.vcf'

    def parse(self,records,sample='bam-sm'):
        self.path.write_text('##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample+'\n'+records)
        return list(manta_evidence(self.path,'study-id',expected_vcf_sample='bam-sm',lengths={'chr1':100,'chr2':100}))

    def test_imprecise_deletion_retains_missing_split_support(self):
        rows=self.parse('chr1\t10\tdel\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=40;IMPRECISE;CIPOS=-2,2\tPR\t20,3\n')
        row=rows[0]
        self.assertEqual(row['paired_read_support'],[20,3])
        self.assertIsNone(row['split_read_support'])
        self.assertTrue(row['source_imprecise'])
        self.assertFalse(row['somatic_origin_classified'])
        self.assertIsNone(row['event_count'])

    def test_breakend_pair_remains_two_records_not_two_events(self):
        text='chr1\t10\tb1\tA\tA]chr2:20]\t.\tPASS\tSVTYPE=BND;MATEID=b2;EVENT=e\tPR:SR\t20,3:10,2\nchr2\t20\tb2\tA\tA]chr1:10]\t.\tPASS\tSVTYPE=BND;MATEID=b1;EVENT=e\tPR:SR\t20,3:10,2\n'
        rows=self.parse(text)
        self.assertEqual(len(rows),2)
        self.assertEqual(rows[0]['source_record_raw'],text.splitlines()[0])
        self.assertEqual(rows[0]['source_mateid'],'b2')
        self.assertTrue(all(row['event_count'] is None for row in rows))

    def test_partial_missing_counts_and_omitted_trailing_field(self):
        row=self.parse('chr1\t10\td\tA\t<DEL>\t.\tLowQual\tSVTYPE=DEL\tPR:SR\t.,0\n')[0]
        self.assertEqual(row['paired_read_support'],[None,0])
        self.assertIsNone(row['split_read_support'])
        self.assertEqual(row['filter_raw'],'LowQual')

    def test_wrong_sample_and_bad_counts_rejected(self):
        with self.assertRaisesRegex(ValueError,'header'):
            self.parse('',sample='wrong-sm')
        for value in ('1,-2','3','1,2,3'):
            with self.assertRaises(ValueError):
                support_pair(value)

    def test_empty_but_valid_vcf_stays_empty(self):
        self.assertEqual(self.parse(''),[])
