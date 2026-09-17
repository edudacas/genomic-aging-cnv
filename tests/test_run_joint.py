import json
from pathlib import Path
import unittest
from unittest.mock import patch
from cnv_execution import sha
from cnv_r_library import snapshot_library
from cnv_run_joint import run_joint,verify_preparation
import test_prepare_joint


class RunJointTests(unittest.TestCase):
    def setUp(self):
        self.fixture=test_prepare_joint.PrepareJointTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.fixture.prepare()
        self.root=self.fixture.output.parent
        library=self.root/'R-library'
        for name in ('facetsCH','pctGCdata'):
            (library/name/'Meta').mkdir(parents=True)
            (library/name/'DESCRIPTION').write_text('synthetic metadata only')
            (library/name/'Meta/package.rds').write_bytes(b'synthetic')
        self.snapshot=self.root/'runtime.json'
        self.snapshot.write_text(json.dumps(snapshot_library(library,required_packages=['facetsCH','pctGCdata'])))
        self.binary=self.root/'fake-Rscript'
        self.binary.write_text('not executable; native execution is mocked')
        self.artifact=dict(path=str(self.binary),sha256=sha(self.binary))

    def invoke(self):
        return run_joint(self.fixture.output,self.snapshot,self.artifact,self.root/'run',timeout_seconds=30)

    def test_changed_preparation_rejected_before_execution(self):
        (self.fixture.output/'fit_config.tsv').write_text('changed')
        with patch('cnv_run_joint.execute') as native,self.assertRaisesRegex(ValueError,'configuration changed'):
            self.invoke()
        native.assert_not_called()

    def test_native_failure_is_not_completion(self):
        with patch('cnv_run_joint.execute',side_effect=RuntimeError('native failure')):
            with self.assertRaisesRegex(RuntimeError,'native failure'):
                self.invoke()
        self.assertTrue((self.root/'run/FAILURE.json').exists())
        self.assertFalse((self.root/'run/COMPLETE.json').exists())

    def test_complete_contract_and_declared_scope(self):
        def fake(command,*args,**kwargs):
            destination=Path(command.outputs[-1]).parent
            destination.mkdir(parents=True)
            for path in command.outputs:
                Path(path).write_text('synthetic native artifact')
            (destination/'segments.tsv').write_text('chromosome\tstart_0based\tend_0based_exclusive\tn_markers\tn_informative_allelic\tmean_log2_ratio\tevent_class\tevaluable_bp\tsomatic_origin_classified\nchr1\t0\t1\t1\t1\t0\tUNCLASSIFIED\tNA\tFALSE\n')
            (destination/'locus_audit.tsv').write_text('catalog_row\tchromosome\tposition_1based\tsignal_status\tnormal_depth\tsample_depth\n1\tchr1\t1\tNATIVE_SIGNAL_RETAINED\t20\t0\n')
            (destination/'signal.tsv').write_text('chrom\tmaploc\tcnlr\tvalor\thet\tsource_catalog_row\n1\t1\t0\t0\t1\t1\n')
            (destination/'exclusions.tsv').write_text('input_row\tchromosome\tposition_1based\treason\n')
            Path(command.outputs[-1]).write_text('JOINT_SEGMENTATION_CANDIDATES\nsensitivity_calibrated=FALSE\nsomatic_origin_classified=FALSE\ncopy_neutral_loh_confirmed=FALSE\n')
            return dict(status='NATIVE_TASK_COMPLETE',outputs={p:dict(sha256=sha(p)) for p in command.outputs})
        with patch('cnv_run_joint.execute',side_effect=fake):
            result=self.invoke()
        self.assertEqual(result['status'],'JOINT_NATIVE_EXECUTION_COMPLETE')
        self.assertFalse(result['sensitivity_calibrated'])
        self.assertFalse(result['system_runtime_dependency_closure_verified'])
        self.assertTrue((self.root/'run/COMPLETE.json').exists())
