import copy
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from cnv_training import verified_training_summaries
from cnv_pileup import COLUMNS
from cnv_territory import Territory
from test_contract import fixture
from test_allelic import policy


class TrainingTests(unittest.TestCase):
    def setUp(self):
        self.manifest = fixture()
        second = copy.deepcopy(self.manifest['samples'][1])
        second.update(sample_id='pon2',person_id='pon2')
        second['bam']['path']='/example/pon2.bam'
        self.manifest['samples'].append(second)
        for sample in self.manifest['samples']:
            sample['reference_sex']='Female'
        temp=tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        root=Path(temp.name)
        self.registry={}
        for sid in ('pon','pon2'):
            path=root/(sid+'.csv')
            path.write_text(','.join(COLUMNS)+'\nchr1,1,A,C,20,20,0,0\nchr1,2,A,C,0,0,5,0\n')
            self.registry[sid]=dict(output=str(path),task_directory=str(root/sid))
        self.catalog=[dict(chromosome='chr1',position=p,ref='A',alt='C',kind='SNP') for p in (1,2,3)]
        self.mask=Territory([],{'chr1':10})
        self.policy=dict(mapq=20,baseq=20,pseudo_spacing=10,max_depth=10000,
                         count_orphans=False,min_read_count=0)

    def run_rows(self):
        return list(verified_training_summaries(self.manifest,'S','Female',self.catalog,
            self.registry,{'chr1':10},mask=self.mask,executable='/example/native',
            pileup_policy=self.policy,marker_policy=policy(),rehash_inputs=True))

    @patch('cnv_training.verify_completed_task',return_value={'fingerprint':'test-only'})
    def test_training_members_and_missingness(self,verify):
        rows=self.run_rows()
        self.assertEqual(verify.call_count,2)
        self.assertEqual(rows[0]['correction_bias'],0)
        self.assertEqual(rows[0]['pon_sample_ids'],['pon','pon2'])
        self.assertEqual(rows[1]['zero_depth_normals'],2)
        self.assertEqual(rows[2]['missing_normals'],2)
        self.assertIsNone(rows[2]['median_observed_depth'])
        self.assertFalse(rows[0]['sensitivity_calibrated'])

    @patch('cnv_training.verify_completed_task',return_value={'fingerprint':'test-only'})
    def test_mask_preserves_summary_but_blocks_correction(self,verify):
        self.mask=Territory([('chr1',0,1)],{'chr1':10})
        row=self.run_rows()[0]
        self.assertEqual(row['raw_median_bias'],0)
        self.assertIsNone(row['correction_bias'])
        self.assertEqual(row['status'],'MASKED_TRAINING_LOCUS')

    def test_foreign_or_omitted_member_rejected_before_native_access(self):
        for name in ('study','calibration','omitted'):
            with self.subTest(name=name):
                original=self.registry.copy()
                if name=='omitted':
                    self.registry.pop('pon2')
                else:
                    self.registry[name]=self.registry['pon']
                with self.assertRaisesRegex(ValueError,'exactly all eligible'):
                    self.run_rows()
                self.registry=original

    def test_unknown_sex_not_silently_excluded(self):
        self.manifest['samples'][-1]['reference_sex']='Unknown'
        with self.assertRaisesRegex(ValueError,'resolve all'):
            self.run_rows()

    @patch('cnv_training.verify_completed_task',side_effect=ValueError('bad receipt'))
    def test_receipt_failure_stops_training(self,verify):
        with self.assertRaisesRegex(ValueError,'bad receipt'):
            self.run_rows()
