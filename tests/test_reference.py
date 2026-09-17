import copy
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from cnv_commands import coverage
from cnv_execution import sha
from cnv_reference import build_reference, reference_inputs, verify_reference
from test_contract import fixture


def artifact(path, text):
    path.write_text(text)
    return {'path': str(path), 'sha256': sha(path)}


class ReferenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.manifest = fixture()
        spec = self.manifest['strata']['S']
        spec['reference_fasta'] = artifact(self.root/'reference.fa', '>chr1\nACGT\n')
        spec['reference_fai'] = artifact(self.root/'reference.fa.fai', 'chr1\t4\t6\t4\t5\n')
        spec['target_bed'] = artifact(self.root/'target.bed', 'chr1\t0\t4\tGENE\n')
        self.implementation = self.root/'implementation'
        self.implementation.write_text('synthetic receipt fixture, not a caller')
        sample = self.manifest['samples'][1]
        sample['bam'] = artifact(self.root/'pon.bam', 'synthetic BAM binding')
        sample['bam_index'] = artifact(self.root/'pon.bai', 'synthetic index binding')
        self.task = self.root/'coverage'
        self.task.mkdir()
        self.table = self.task/'pon.targetcoverage.cnn'
        self.table.write_text('chromosome\tstart\tend\tgene\tlog2\tdepth\nchr1\t0\t4\tGENE\t0\t0\n')
        command = coverage(sample['bam']['path'], spec['target_bed']['path'], self.table,
                           executable='cnvkit.py', processes=1, min_mapq=20)
        self.frozen = dict(command=asdict(command),
                           inputs={x['path']: x['sha256'] for x in
                                   (sample['bam'], sample['bam_index'], spec['target_bed'])},
                           implementation={str(self.implementation): sha(self.implementation)})
        self.seal()
        self.registry = {'pon': {'target': str(self.task)}}

    def seal(self):
        fingerprint = hashlib.sha256(json.dumps(self.frozen, sort_keys=True).encode()).hexdigest()
        (self.task/'STARTED.json').write_text(json.dumps(dict(self.frozen, fingerprint=fingerprint)))
        (self.task/'COMPLETE.json').write_text(json.dumps(dict(
            status='NATIVE_TASK_COMPLETE', fingerprint=fingerprint,
            outputs={str(self.table): {'sha256': sha(self.table), 'bytes': self.table.stat().st_size}})))

    def run_build(self, **overrides):
        manifest, registry = self.root/'manifest.json', self.root/'registry.json'
        manifest.write_text(json.dumps(self.manifest))
        registry.write_text(json.dumps(self.registry))
        options = dict(executable='cnvkit.py', processes=1, min_mapq=20,
                       minimum_pon_samples=1, implementation_files=[self.implementation],
                       timeout_seconds=120)
        options.update(overrides)
        return build_reference(manifest, registry, 'S', self.root/'reference', **options)

    def test_only_pon_enters_reference_zero_depth_bins_retained(self):
        with patch('cnv_reference.execute', return_value={'status': 'NATIVE_TASK_COMPLETE'}) as run:
            receipt = self.run_build()
        self.assertEqual(receipt['reference_members'], ['pon'])
        self.assertEqual(receipt['calibration_samples_used'], 0)
        argv = run.call_args.args[0].argv
        self.assertIn(str(self.table), argv)
        self.assertNotIn('/example/study.bam', argv)

    def test_rejects_study_and_calibration_registry_members(self):
        for role in ('study', 'calibration'):
            self.registry[role] = {'target': str(self.task)}
            with self.subTest(role=role), self.assertRaisesRegex(ValueError, 'exactly'):
                self.run_build()
            del self.registry[role]

    def test_rejects_changed_bam_binding_and_policy(self):
        self.manifest['samples'][1]['bam']['sha256'] = 'f'*64
        with self.assertRaisesRegex(ValueError, 'binding'):
            self.run_build()
        with self.assertRaisesRegex(ValueError, 'policy'):
            self.run_build(min_mapq=30)

    def test_rejects_rehashed_wrong_bins(self):
        self.table.write_text(self.table.read_text().replace('chr1\t0\t4', 'chr1\t1\t4'))
        self.seal()
        with self.assertRaisesRegex(ValueError, 'bins differ'):
            self.run_build()

    def test_rejects_mutated_table_and_insufficient_pon(self):
        with self.assertRaisesRegex(ValueError, 'insufficient'):
            self.run_build(minimum_pon_samples=2)
        self.table.write_text(self.table.read_text() + '\n')
        with self.assertRaisesRegex(ValueError, 'output changed'):
            self.run_build()

    def test_rejects_reference_sample_with_other_assay(self):
        self.manifest['samples'][1]['capture_id'] = 'OTHER'
        with self.assertRaisesRegex(ValueError, 'mismatch'):
            self.run_build()

    def make_reference_receipt(self):
        # Receipt contract fixture, not a native CNVkit execution claim.
        with patch('cnv_reference.execute', return_value={'status':'NATIVE_TASK_COMPLETE'}):
            self.run_build()
        self.reference_dir = self.root/'reference'
        self.reference_dir.mkdir()
        self.options = dict(executable='cnvkit.py', processes=1, min_mapq=20,
                            minimum_pon_samples=1)
        command, inputs, _ = reference_inputs(self.root/'manifest.json',
            self.root/'registry.json', 'S', self.reference_dir, **self.options)
        output = Path(command.outputs[0])
        output.write_text('chromosome\tstart\tend\tgene\tlog2\nchr1\t0\t4\tGENE\t0\n')
        frozen = dict(command=asdict(command), inputs={str(p):sha(p) for p in inputs},
                      implementation={str(self.implementation):sha(self.implementation)})
        fingerprint = hashlib.sha256(json.dumps(frozen,sort_keys=True).encode()).hexdigest()
        (self.reference_dir/'STARTED.json').write_text(json.dumps(dict(frozen,fingerprint=fingerprint)))
        (self.reference_dir/'COMPLETE.json').write_text(json.dumps(dict(
            status='NATIVE_TASK_COMPLETE',fingerprint=fingerprint,
            outputs={str(output):dict(sha256=sha(output),bytes=output.stat().st_size)})))

    def verify_reference_fixture(self):
        return verify_reference(self.root/'manifest.json', self.root/'registry.json',
                                'S', self.reference_dir, **self.options)

    def test_read_only_reference_verification(self):
        self.make_reference_receipt()
        with patch('cnv_reference.execute') as caller:
            result = self.verify_reference_fixture()
        caller.assert_not_called()
        self.assertEqual(result['reference_members'],['pon'])
        self.assertFalse(result['clinical_reference_suitability_verified'])
        self.assertFalse(result['sex_chromosome_reference_policy_validated'])

    def test_reference_changed_after_build_rejected(self):
        self.make_reference_receipt()
        output = self.reference_dir/'reference.cnn'
        output.write_text(output.read_text().replace('GENE\t0','GENE\t1'))
        with self.assertRaisesRegex(ValueError,'output bytes changed'):
            self.verify_reference_fixture()

    def test_reference_cannot_rebind_changed_registry(self):
        self.make_reference_receipt()
        path = self.root/'registry.json'
        path.write_text(path.read_text()+'\n')
        with self.assertRaisesRegex(ValueError,'input binding'):
            self.verify_reference_fixture()


if __name__ == '__main__':
    unittest.main()
