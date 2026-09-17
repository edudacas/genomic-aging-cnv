from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from cnv_commands import joint_fit


class JointCommandTests(unittest.TestCase):
    def test_explicit_argv_and_all_outputs(self):
        cmd=joint_fit('/a/matrix','/a/bias','/a/config','/a/map','/a/task/result',
            rscript='/usr/bin/Rscript',wrapper='/a/facets_run.R',library='/a/R library')
        self.assertEqual(cmd.argv[1],'--vanilla')
        self.assertIn('/a/R library',cmd.argv)
        self.assertNotIn('--bestnormal',cmd.argv)
        self.assertEqual(len(cmd.outputs),7)
        self.assertTrue(all(Path(p).parent==Path('/a/task/result') for p in cmd.outputs))

    def test_relative_paths_rejected(self):
        with self.assertRaisesRegex(ValueError,'absolute'):
            joint_fit('relative','/b','/c','/d','/e',rscript='/r',wrapper='/w',library='/l')

    @unittest.skipUnless(shutil.which('Rscript'),'Rscript unavailable')
    def test_cli_help_is_independent_of_cwd_and_native_install(self):
        script=Path(__file__).resolve().parents[1]/'facets_run.R'
        with tempfile.TemporaryDirectory() as directory:
            result=subprocess.run(['Rscript','--vanilla',str(script),'--help'],cwd=directory,
                capture_output=True,text=True,timeout=30)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertIn('--chromosome-map',result.stdout)

    @unittest.skipUnless(shutil.which('Rscript'),'Rscript unavailable')
    def test_cli_rejects_old_upstream_flags(self):
        script=Path(__file__).resolve().parents[1]/'facets_run.R'
        result=subprocess.run(['Rscript','--vanilla',str(script),'--bestnormal','TRUE'],
            capture_output=True,text=True,timeout=30)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('exact named path arguments',result.stderr)
