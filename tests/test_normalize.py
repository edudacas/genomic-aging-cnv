from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import unittest
from unittest.mock import patch

from cnv_commands import coverage
from cnv_execution import sha
from cnv_normalize import normalize_sample, validate_normalized_table
import test_reference


class NormalizeTests(unittest.TestCase):
    def setUp(self):
        self.f = test_reference.ReferenceTests()
        self.f.setUp()
        self.addCleanup(self.f.doCleanups)
        self.root = self.f.root
        self.sample = self.f.manifest['samples'][0]
        for key in ('bam','bam_index'):
            self.sample[key] = test_reference.artifact(self.root/('study-'+key), 'synthetic '+key)
        self.f.make_reference_receipt()
        self.task = self.root/'sample-coverage'
        self.task.mkdir()
        table = self.task/'target.cnn'
        table.write_text(self.f.table.read_text())
        bed = self.f.manifest['strata']['S']['target_bed']
        command = coverage(self.sample['bam']['path'],bed['path'],table,
                           executable='cnvkit.py',processes=1,min_mapq=20)
        frozen = dict(command=asdict(command),
            inputs={x['path']:x['sha256'] for x in (self.sample['bam'],self.sample['bam_index'],bed)},
            implementation={str(self.f.implementation):sha(self.f.implementation)})
        fingerprint = hashlib.sha256(json.dumps(frozen,sort_keys=True).encode()).hexdigest()
        (self.task/'STARTED.json').write_text(json.dumps(dict(frozen,fingerprint=fingerprint)))
        (self.task/'COMPLETE.json').write_text(json.dumps(dict(status='NATIVE_TASK_COMPLETE',
            fingerprint=fingerprint,outputs={str(table):dict(sha256=sha(table),bytes=table.stat().st_size)})))
        self.registry = self.root/'sample-registry.json'
        self.registry.write_text(json.dumps({self.sample['sample_id']:{'target':str(self.task)}}))

    def invoke(self, sample_id=None):
        return normalize_sample(self.root/'manifest.json',sample_id or self.sample['sample_id'],
            self.registry,self.root/'registry.json',self.f.reference_dir,self.root/'fix',
            **self.f.options, implementation_files=[self.f.implementation],timeout_seconds=30)

    def test_verified_reference_and_sample_feed_fix(self):
        def native(command,*args,**kwargs):
            Path(command.outputs[0]).write_text('chromosome\tstart\tend\tlog2\tweight\nchr1\t0\t4\t0\t1\n')
            return {'status':'NATIVE_TASK_COMPLETE'}
        with patch('cnv_normalize.execute',side_effect=native) as call:
            result = self.invoke()
        argv = call.call_args.args[0].argv
        self.assertEqual(argv[1],'fix')
        self.assertEqual(argv[4],str(self.f.reference_dir/'reference.cnn'))
        self.assertEqual((self.root/'fix/off_target_disabled.cnn').read_text(),
                         'chromosome\tstart\tend\tgene\tlog2\tdepth\n')
        self.assertFalse(result['sex_chromosome_reference_policy_validated'])
        self.assertFalse(result['somatic_origin_classified'])

    def test_empty_normalization_is_not_evaluable(self):
        path=self.root/'empty.cnr'
        path.write_text('chromosome\tstart\tend\tlog2\tweight\n')
        with self.assertRaisesRegex(ValueError,'retained no bins'):
            validate_normalized_table(path)

    def test_pon_not_normalized_as_study(self):
        with patch('cnv_normalize.execute') as call,self.assertRaisesRegex(ValueError,'study/calibration'):
            self.invoke('pon')
        call.assert_not_called()

    def test_changed_reference_blocks_fix(self):
        (self.f.reference_dir/'reference.cnn').write_text('changed')
        with patch('cnv_normalize.execute') as call,self.assertRaises(ValueError):
            self.invoke()
        call.assert_not_called()

    def test_extra_sample_registry_blocks_fix(self):
        self.registry.write_text(json.dumps({'other':{'target':str(self.task)}}))
        with patch('cnv_normalize.execute') as call,self.assertRaisesRegex(ValueError,'registry'):
            self.invoke()
        call.assert_not_called()
