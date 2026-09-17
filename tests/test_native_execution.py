from pathlib import Path
import sys
import subprocess
import tempfile
import unittest
from cnv_commands import Command, snp_pileup, facets_ch, segment, fix
from cnv_execution import execute, verify_completed_task


class NativeTests(unittest.TestCase):
    def test_runtime_environment_is_frozen_and_does_not_mutate_parent(self):
        import os
        before = os.environ.get('R_LIBS_USER')
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            output = root/'environment.txt'
            command = Command((sys.executable,'-c',
                'import os; from pathlib import Path; Path("environment.txt").write_text(os.environ["R_LIBS_USER"])'),
                (str(output),))
            kwargs = dict(implementation_files=[__file__],timeout_seconds=10,
                          runtime_environment={'R_LIBS_USER':'/isolated/test-library'})
            execute(command,root,[],**kwargs)
            self.assertEqual(output.read_text(),'/isolated/test-library')
            self.assertEqual(os.environ.get('R_LIBS_USER'),before)
            self.assertTrue(verify_completed_task(command,root,{},rehash_inputs=True)['outputs_verified'])
            kwargs['runtime_environment']={'R_LIBS_USER':'/different/library'}
            with self.assertRaisesRegex(ValueError,'inputs/implementation changed'):
                execute(command,root,[],**kwargs)

    def test_runtime_environment_rejects_undeclared_or_secret_variables(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            command = Command((sys.executable,'-c','pass'),(str(root/'out'),))
            with self.assertRaisesRegex(ValueError,'runtime environment'):
                execute(command,root,[],implementation_files=[__file__],timeout_seconds=10,
                        runtime_environment={'API_TOKEN':'do-not-record'})

    def test_cnvkit_fix_reference_is_positional(self):
        command = fix('target.cnn', 'antitarget.cnn', 'reference.cnn', 'sample.cnr',
                      executable='cnvkit.py')
        self.assertEqual(command.argv, ('cnvkit.py', 'fix', 'target.cnn',
                         'antitarget.cnn', 'reference.cnn', '-o', 'sample.cnr'))

    def test_timeout_does_not_issue_completion(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            command = Command((sys.executable, '-c', 'import time; time.sleep(20)'),
                              (str(root/'output'),))
            with self.assertRaises(subprocess.TimeoutExpired):
                execute(command, root, [], implementation_files=[__file__], timeout_seconds=0.1)
            self.assertFalse((root/'COMPLETE.json').exists())
            self.assertTrue((root/'FAILURE.json').exists())

    def test_pileup_preserves_overlap_detection_and_native_order(self):
        command = snp_pileup('sample.bam', 'sites.vcf.gz', 'counts.gz', executable='snp-pileup',
                            mapq=20, baseq=20, pseudo_spacing=100, max_depth=10000,
                            count_orphans=False, min_read_count=0)
        self.assertEqual(command.argv[-3:], ('sites.vcf.gz', 'counts.gz', 'sample.bam'))
        self.assertNotIn('-x', command.argv)
        self.assertNotIn('-A', command.argv)

    def test_unknown_sex_does_not_default_to_female(self):
        with self.assertRaisesRegex(ValueError, 'resolved reference sex'):
            facets_ch('pu', 'db', 'out', 'sample', 'WES', 'Unknown', rscript='Rscript',
                      wrapper='wrapper.R', bin_size=100, cval=25, gamma=0.25)

    def test_invalid_parameters_fail(self):
        for threshold in (float('nan'), -1, 1, True):
            with self.subTest(threshold=threshold), self.assertRaises(ValueError):
                segment('in', 'out', executable='cnvkit.py', processes=1, threshold=threshold)

    def test_execution_reuses_only_identical_validated_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root/'input'; source.write_text('data')
            output = root/'task'/'table.tsv'
            command = Command((sys.executable, '-c',
                'from pathlib import Path; Path("table.tsv").write_text("start\\tend\\n0\\t10\\n")'),
                (str(output),), ('start', 'end'))
            kwargs = dict(implementation_files=[__file__], timeout_seconds=10)
            first = execute(command, root/'task', [source], **kwargs)
            self.assertFalse(first['reused'])
            self.assertTrue(execute(command, root/'task', [source], **kwargs)['reused'])
            output.write_text('start\tend\n1\t10\n')
            with self.assertRaisesRegex(ValueError, 'outputs changed'):
                execute(command, root/'task', [source], **kwargs)

    def test_success_exit_without_output_fails_and_preserves_attempt(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            command = Command((sys.executable, '-c', 'pass'), (str(root/'missing.tsv'),))
            kwargs = dict(implementation_files=[__file__], timeout_seconds=10)
            with self.assertRaisesRegex(ValueError, 'missing/empty'):
                execute(command, root, [], **kwargs)
            self.assertTrue((root/'FAILURE.json').exists())
            self.assertFalse((root/'COMPLETE.json').exists())
            with self.assertRaisesRegex(ValueError, 'incomplete task'):
                execute(command, root, [], **kwargs)


if __name__ == '__main__':
    unittest.main()
