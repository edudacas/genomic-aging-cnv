from dataclasses import asdict
import json
import unittest
from unittest.mock import patch
from cnv_train_panel import train
from cnv_joint_inputs import verify_joint_inputs
from cnv_allelic import summarize_marker
from test_allelic import policy
import test_prepare_pair


class JointInputTests(unittest.TestCase):
    def setUp(self):
        self.f=test_prepare_pair.PreparePairTests()
        self.f.setUp()
        self.addCleanup(self.f.doCleanups)
        self.f.run_preparation()
        self.panel=self.f.root/'panel'
        self.panel_config=self.f.root/'panel.json'
        config=json.loads(self.f.config.read_text())
        del config['sample_id']
        del config['reference_id']
        config.update(stratum='S',reference_sex='Female',marker_policy=asdict(policy()),
                      native_chromosome_map={'chr1':1})
        self.panel_config.write_text(json.dumps(config))
        self.row=dict(chromosome='chr1',position=1,ref='A',alt='C',kind='SNP',
            **summarize_marker({'a':(20,20),'b':(20,20)},['a','b'],policy()))
        self.build_panel()

    def build_panel(self):
        # Artifact-link checks only. Native role/receipt validation has separate tests.
        with patch('cnv_train_panel.verified_training_summaries',return_value=iter([self.row])):
            train(self.f.manifest,self.panel_config,self.panel)

    def verify(self):
        return verify_joint_inputs(self.f.manifest,self.f.output,self.f.config,self.panel,self.panel_config)

    def test_matching_artifacts_without_fit_claim(self):
        receipt=self.verify()
        self.assertEqual(receipt['status'],'JOINT_INPUT_ARTIFACTS_VERIFIED')
        self.assertFalse(receipt['native_joint_fit_validated'])
        self.assertFalse(receipt['upstream_bams_rehashed_this_check'])

    def test_changed_output_rejected(self):
        (self.panel/'allelic_bias.tsv').write_text('changed')
        with self.assertRaisesRegex(ValueError,'output bytes'):
            self.verify()

    def test_started_receipt_edit_rejected(self):
        path=self.panel/'STARTED.json'
        path.write_text(path.read_text()+' ')
        with self.assertRaisesRegex(ValueError,'stage receipt'):
            self.verify()

    def test_failure_marker_rejected(self):
        (self.panel/'FAILURE.json').write_text('{}')
        with self.assertRaisesRegex(ValueError,'failed attempt'):
            self.verify()

    def test_valid_but_wrong_sex_panel_rejected(self):
        self.panel=self.f.root/'male-panel'
        config=json.loads(self.panel_config.read_text())
        config['reference_sex']='Male'
        self.panel_config=self.f.root/'male-panel.json'
        self.panel_config.write_text(json.dumps(config))
        self.build_panel()
        with self.assertRaisesRegex(ValueError,'reference sex'):
            self.verify()
