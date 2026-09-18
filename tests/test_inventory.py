import unittest
from cnv_inventory import reconcile


class InventoryTests(unittest.TestCase):
    def setUp(self):
        self.link = dict(individual_id='i', sample_id='s', person_id='p', analysis_id='a',
                         fap='f', assay_id='kit-full', assay_source='exact',
                         bam_path='/server/a.bam', bai_path='/server/a.bam.bai')
        self.row = dict(individual_id='i', analysis_id='a', fap='f', assay_id='kit',
                        BAM_path='/server/a.bam', batch_id='b', sex='Female')
        self.assay = dict(assay_id='kit-full', assay_short_id='kit',
                          assay_name_master='Test', proposed_pon_stratum='Test',
                          pon_group_equivalence_basis='declared',
                          capture_bed_identity_status='NOT_VERIFIED')

    def run_case(self):
        return reconcile([self.link], [self.row], [self.assay], 1)

    def test_unknown_is_not_negative(self):
        row = self.run_case()[0]
        self.assertIsNone(row['sex_inferred'])
        self.assertIsNone(row['mean_target_depth'])
        self.assertEqual(row['umi_status'], 'UNKNOWN')
        self.assertEqual(row['status'], 'HOLD')

    def test_identity_conflict(self):
        self.row['fap'] = 'other'
        with self.assertRaisesRegex(ValueError, 'identity'):
            self.run_case()

    def test_assay_conflict(self):
        self.row['assay_id'] = 'other'
        with self.assertRaisesRegex(ValueError, 'assay conflict'):
            self.run_case()

    def test_duplicate_people(self):
        other = dict(self.link, individual_id='i2', sample_id='s2', analysis_id='a2', bam_path='/b')
        with self.assertRaisesRegex(ValueError, 'person_id'):
            reconcile([self.link, other], [self.row], [self.assay], 2)

    def test_path_difference_retained(self):
        self.row['BAM_path'] = '/old/a.bam'
        row = self.run_case()[0]
        self.assertIn('REVIEW_PATH_DIFFERENCE', row['hold_reasons'])
        self.assertEqual(row['workbook_bam_path'], '/old/a.bam')

    def test_wrong_membership(self):
        self.row['individual_id'] = 'other'
        with self.assertRaisesRegex(ValueError, 'membership'):
            self.run_case()

    def test_missing_batch(self):
        self.row['batch_id'] = ''
        self.assertIn('BATCH_UNKNOWN', self.run_case()[0]['hold_reasons'])
