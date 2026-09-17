import copy
from pathlib import Path
import tempfile
import unittest
from cnv_calibration import summarize_calibration
from cnv_execution import sha
from test_contract import fixture


class CalibrationTests(unittest.TestCase):
    def setUp(self):
        temp=tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        evidence=Path(temp.name)/'evidence.json'
        evidence.write_text('{"scope":"synthetic test evidence only"}')
        self.artifact=dict(path=str(evidence),sha256=sha(evidence))
        self.case=dict(case_id='a',source_sample_id='calibration',truth_source='SIMULATED',
            event_class='LOSS',clone_fraction=.2,size_band='predeclared-size',depth_band='predeclared-depth',
            informative_loci_band='predeclared-density',truth_artifact=self.artifact)
        self.plan=dict(protocol_sha256='a'*64,cases=[self.case])

    def outcome(self,cid,status):
        return dict(case_id=cid,status=status,protocol_sha256='a'*64,result_artifact=self.artifact)

    def summarize(self,outcomes):
        return summarize_calibration(fixture(),self.plan,outcomes,check_artifact_files=True)

    def test_missing_run_is_not_missed_event(self):
        result=self.summarize([])
        cell=result['cells'][0]
        self.assertEqual(cell['counts'],{'NOT_RUN':1})
        self.assertIsNone(cell['conditional_detection_fraction'])
        self.assertFalse(result['calibration_gate_closed'])

    def test_non_evaluable_excluded_but_visible_in_denominator_accounting(self):
        self.plan['cases'] += [dict(self.case,case_id='b'),dict(self.case,case_id='c')]
        cell=self.summarize([self.outcome('a','DETECTED'),self.outcome('b','MISSED'),
                             self.outcome('c','NOT_EVALUABLE')])['cells'][0]
        self.assertEqual(cell['conditional_detection_fraction'],.5)
        self.assertEqual(cell['evaluable_cases'],2)
        self.assertEqual(cell['planned_cases'],3)
        self.assertEqual(cell['unique_source_samples'],1)
        self.assertFalse(cell['all_planned_cases_evaluable'])
        self.assertFalse(cell['independent_replicates_assumed'])

    def test_simulated_and_orthogonal_results_not_pooled(self):
        self.plan['cases'].append(dict(self.case,case_id='b',truth_source='ORTHOGONAL'))
        cells=self.summarize([self.outcome('a','DETECTED'),self.outcome('b','MISSED')])['cells']
        self.assertEqual(len(cells),2)
        self.assertEqual({c['conditional_detection_fraction'] for c in cells},{0,1})

    def test_unknown_orthogonal_fraction_remains_missing(self):
        self.case.update(truth_source='ORTHOGONAL',clone_fraction=None)
        cell=self.summarize([self.outcome('a','DETECTED')])['cells'][0]
        self.assertIsNone(cell['clone_fraction'])

    def test_null_control_has_separate_false_positive_fraction(self):
        self.case.update(event_class='NULL_CONTROL',clone_fraction=None)
        cell=self.summarize([self.outcome('a','FALSE_POSITIVE')])['cells'][0]
        self.assertEqual(cell['conditional_false_positive_fraction'],1)
        self.assertIsNone(cell['conditional_detection_fraction'])

    def test_study_or_pon_donor_is_rejected(self):
        for sample in ('study','pon'):
            self.case['source_sample_id']=sample
            with self.assertRaisesRegex(ValueError,'study or PoN'):
                self.summarize([])

    def test_protocol_mismatch_and_duplicate_rejected(self):
        outcome=self.outcome('a','DETECTED')
        outcome['protocol_sha256']='b'*64
        with self.assertRaisesRegex(ValueError,'different calibration protocol'):
            self.summarize([outcome])
        with self.assertRaisesRegex(ValueError,'duplicate'):
            self.summarize([self.outcome('a','DETECTED')]*2)

    def test_changed_truth_artifact_rejected(self):
        Path(self.artifact['path']).write_text('changed')
        with self.assertRaisesRegex(ValueError,'checksum mismatch'):
            self.summarize([])
