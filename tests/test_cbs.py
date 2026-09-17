import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from cnv_cbs import run_cbs
from cnv_execution import sha
from cnv_r_library import snapshot_library


class CBSTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.fai = self.root/'ref.fai'
        self.fai.write_text('chr1\t100\t6\t100\t101\n')
        self.manifest = self.root/'manifest.json'
        self.manifest.write_text(json.dumps({'strata':{'S':{'reference_fai':{'path':str(self.fai)}}}}))
        self.norm = self.root/'norm'
        self.norm.mkdir()
        for name in ('STARTED.json','COMPLETE.json','sample.cnr'):
            (self.norm/name).write_text('synthetic normalization fixture')
        library = self.root/'library'
        (library/'DNAcopy/Meta').mkdir(parents=True)
        (library/'DNAcopy/DESCRIPTION').write_text('synthetic metadata')
        (library/'DNAcopy/Meta/package.rds').write_text('synthetic metadata')
        self.snapshot = self.root/'runtime.json'
        self.snapshot.write_text(json.dumps(snapshot_library(library,required_packages=['DNAcopy'])))
        self.rscript = self.root/'Rscript'
        self.rscript.write_text('fake executable; execution mocked')
        self.config = self.root/'config.json'
        self.config.write_text(json.dumps(dict(manifest_path=str(self.manifest),sample_id='study',
            sample_registry_path='unused',pon_registry_path='unused',reference_directory='unused',
            directory=str(self.norm),executable='cnvkit.py',processes=1,min_mapq=20,minimum_pon_samples=2)))
        self.verified = dict(sample_id='study',stratum='S',normalized_table=str(self.norm/'sample.cnr'),
                             manifest_sha256=sha(self.manifest),sample_bam_sha256='0'*64,reference_build='GRCh38',
                             inputs=[self.manifest,self.fai])

    def invoke(self):
        return run_cbs(self.config,self.snapshot,dict(path=str(self.rscript),sha256=sha(self.rscript)),
            self.root/'run',processes=1,threshold=.001,implementation_files=[__file__],timeout_seconds=10)

    def fake_execute(self,command,*args,**kwargs):
        path = Path(command.outputs[0])
        path.parent.mkdir(parents=True)
        path.write_text('chromosome\tstart\tend\tlog2\tprobes\nchr1\t0\t100\t-1\t10\n')
        self.assertIn('--rscript-path',command.argv)
        self.assertEqual(kwargs['runtime_environment']['R_PROFILE_USER'],'/dev/null')
        return dict(outputs={str(path):dict(sha256=sha(path))})

    def test_complete_and_provenance_rechecked(self):
        with patch('cnv_cbs.verify_normalization',return_value=self.verified) as check, \
                patch('cnv_cbs.execute',side_effect=self.fake_execute):
            result = self.invoke()
        self.assertEqual(check.call_count,2)
        self.assertEqual(result['segment_rows'],1)
        self.assertFalse(result['somatic_origin_classified'])
        self.assertTrue((self.root/'run/COMPLETE.json').exists())

    def test_native_failure_cannot_complete(self):
        with patch('cnv_cbs.verify_normalization',return_value=self.verified), \
                patch('cnv_cbs.execute',side_effect=RuntimeError('CBS failed')):
            with self.assertRaisesRegex(RuntimeError,'CBS failed'):
                self.invoke()
        self.assertTrue((self.root/'run/FAILURE.json').exists())
        self.assertFalse((self.root/'run/COMPLETE.json').exists())

    def test_bad_normalization_stops_before_native_launch(self):
        with patch('cnv_cbs.verify_normalization',side_effect=ValueError('invalid normalization')), \
                patch('cnv_cbs.execute') as native:
            with self.assertRaisesRegex(ValueError,'invalid normalization'):
                self.invoke()
        native.assert_not_called()
        self.assertFalse((self.root/'run').exists())

    def test_out_of_bounds_segment_preserved_as_failure(self):
        def outside(command,*args,**kwargs):
            result = self.fake_execute(command,*args,**kwargs)
            Path(command.outputs[0]).write_text('chromosome\tstart\tend\tlog2\tprobes\nchr1\t0\t101\t-1\t10\n')
            return result
        with patch('cnv_cbs.verify_normalization',return_value=self.verified), \
                patch('cnv_cbs.execute',side_effect=outside):
            with self.assertRaises(ValueError):
                self.invoke()
        self.assertTrue((self.root/'run/FAILURE.json').exists())
        self.assertFalse((self.root/'run/COMPLETE.json').exists())
