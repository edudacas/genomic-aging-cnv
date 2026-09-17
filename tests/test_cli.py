import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from cnv_cli import main
from test_contract import fixture


class CliTests(unittest.TestCase):
    def test_comparison_export_does_not_launch_callers(self):
        with tempfile.TemporaryDirectory() as temp:
            artifact=Path(temp)/'artifact.json'
            artifact.write_text(json.dumps({'path':'/receipt','sha256':'a'*64}))
            with patch('cnv_compare_export.export_comparison',return_value={'status':'TEST'}) as export, \
                    patch('cnv_cbs.run_cbs') as cbs,patch('cnv_run_joint.run_joint') as joint, \
                    contextlib.redirect_stdout(io.StringIO()):
                main(['export-comparison','--receipt-artifact',str(artifact),
                      '--sample-id','s','--output','out','--batch-size','2'])
            export.assert_called_once_with({'path':'/receipt','sha256':'a'*64},'s','out',batch_size=2)
            cbs.assert_not_called()
            joint.assert_not_called()

    def test_cnvkit_stages_require_execution_before_reading_inputs(self):
        for stage in ('run-reference','run-normalize','run-cbs'):
            args=[stage,'--config','/missing','--output','/missing',
                  '--implementation-file','/missing','--timeout-seconds','30']
            if stage=='run-cbs':
                args+=['--library-snapshot','/missing','--rscript-artifact','/missing',
                       '--processes','1','--threshold','.001']
            with self.subTest(stage=stage),contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit) as err:
                main(args)
            self.assertEqual(err.exception.code,2)

    def test_normalization_command_does_not_chain_cbs(self):
        with tempfile.TemporaryDirectory() as temp:
            config=Path(temp)/'config.json'
            data=dict(manifest_path='manifest',sample_id='study',sample_registry_path='sample-registry',
                pon_registry_path='pon-registry',reference_directory='reference',
                executable='cnvkit.py',processes=1,min_mapq=20,minimum_pon_samples=2)
            config.write_text(json.dumps(data))
            with patch('cnv_normalize.normalize_sample',return_value={'status':'TEST'}) as run, \
                    patch('cnv_cbs.run_cbs') as cbs,contextlib.redirect_stdout(io.StringIO()):
                main(['run-normalize','--config',str(config),'--output','output',
                      '--implementation-file','caller','--timeout-seconds','30','--execute'])
            run.assert_called_once_with(**data,output_directory='output',
                implementation_files=['caller'],timeout_seconds=30.0)
            cbs.assert_not_called()

    def test_cnvkit_config_cannot_override_output(self):
        with tempfile.TemporaryDirectory() as temp:
            config=Path(temp)/'config.json'
            config.write_text(json.dumps({'output_directory':'different'}))
            with patch('cnv_reference.build_reference') as run,contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):
                main(['run-reference','--config',str(config),'--output','output',
                      '--implementation-file','caller','--timeout-seconds','30','--execute'])
            run.assert_not_called()

    def test_help_does_not_import_or_launch_runtime(self):
        output=io.StringIO()
        with contextlib.redirect_stdout(output),self.assertRaises(SystemExit) as error:
            main(['--help'])
        self.assertEqual(error.exception.code,0)
        self.assertIn('prepare-joint',output.getvalue())

    def test_run_requires_explicit_execution_before_input_read(self):
        args=['run-joint','--prepared','/missing','--library-snapshot','/missing',
              '--rscript-artifact','/missing','--output','/missing','--timeout-seconds','30']
        with patch('cnv_run_joint.run_joint') as run,contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit) as error:
            main(args)
        self.assertEqual(error.exception.code,2)
        run.assert_not_called()

    def test_manifest_validation_reports_scope(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'manifest.json'
            path.write_text(json.dumps(fixture()))
            output=io.StringIO()
            with contextlib.redirect_stdout(output):
                self.assertEqual(main(['validate-manifest','--manifest',str(path),
                                      '--validation-stage','reference_build']),0)
            result=json.loads(output.getvalue())
            self.assertFalse(result['files_verified'])
            self.assertFalse(result['cohort_executed'])

    def test_pair_preparation_does_not_chain_native_fit(self):
        with (patch('cnv_prepare_pair.prepare',return_value={'status':'TEST_PREPARATION'}) as prepare,
              patch('cnv_run_joint.run_joint') as run):
            with contextlib.redirect_stdout(io.StringIO()):
                main(['prepare-pair','--manifest','m','--config','c','--output','o'])
            prepare.assert_called_once_with('m','c','o')
            run.assert_not_called()
