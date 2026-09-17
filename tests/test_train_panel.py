from dataclasses import asdict
import json
import unittest
from unittest.mock import patch

from cnv_train_panel import train
from cnv_allelic import summarize_marker
from test_allelic import policy
import test_prepare_pair


class TrainPanelTests(unittest.TestCase):
    def setUp(self):
        self.fixture=test_prepare_pair.PreparePairTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.config=json.loads(self.fixture.config.read_text())
        del self.config['sample_id']
        del self.config['reference_id']
        self.config.update(stratum='S',reference_sex='Female',marker_policy=asdict(policy()),
                           native_chromosome_map={'chr1':1})
        self.fixture.config.write_text(json.dumps(self.config))
        self.row=dict(chromosome='chr1',position=1,ref='A',alt='C',kind='SNP',
                      **summarize_marker({'a':(20,20),'b':(20,20)},['a','b'],policy()))

    def run_stage(self,rows):
        # Stage-only tests: native receipt and role checks are tested elsewhere.
        with patch('cnv_train_panel.verified_training_summaries',return_value=rows):
            return train(self.fixture.manifest,self.fixture.config,self.fixture.output)

    def test_failed_stream_cannot_publish_complete(self):
        def broken():
            yield self.row
            raise ValueError('trailing invalid locus')
        with self.assertRaisesRegex(ValueError,'trailing invalid'):
            self.run_stage(broken())
        self.assertTrue((self.fixture.output/'FAILURE.json').exists())
        self.assertFalse((self.fixture.output/'COMPLETE.json').exists())

    def test_input_mutation_cannot_publish_complete(self):
        def changed():
            yield self.row
            self.fixture.config.write_text('{}')
        with self.assertRaisesRegex(ValueError,'changed during'):
            self.run_stage(changed())
        self.assertFalse((self.fixture.output/'COMPLETE.json').exists())

    def test_no_stable_markers_is_explicit_not_fabricated(self):
        row=dict(self.row,status='EXCLUDED_TRAINING_MARKER',correction_bias=None)
        receipt=self.run_stage(iter([row]))
        self.assertEqual(receipt['native_bias_rows'],0)
        self.assertFalse(receipt['has_native_bias_rows'])
        self.assertFalse(receipt['sensitivity_calibrated'])
        self.assertEqual(len((self.fixture.output/'allelic_bias.tsv').read_text().splitlines()),1)

    def test_invalid_native_map_rejected_before_output(self):
        self.config['native_chromosome_map']={'chr1':True}
        self.fixture.config.write_text(json.dumps(self.config))
        with self.assertRaisesRegex(ValueError,'chromosome map'):
            self.run_stage(iter([self.row]))
        self.assertFalse(self.fixture.output.exists())
