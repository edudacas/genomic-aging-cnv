from pathlib import Path
import shutil
import subprocess
import unittest


class FacetsAdapterTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which('Rscript'), 'Rscript unavailable')
    def test_file_boundary_with_injected_engine(self):
        root = Path(__file__).resolve().parents[1]
        completed = subprocess.run(['Rscript', 'tests/test_facets_io.R'], cwd=root,
                                   capture_output=True, text=True, timeout=60)
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)

    @unittest.skipUnless(shutil.which('Rscript'), 'Rscript unavailable')
    def test_segment_export_coordinates_and_missingness(self):
        root = Path(__file__).resolve().parents[1]
        completed = subprocess.run(['Rscript', 'tests/test_facets_segments.R'], cwd=root,
                                   capture_output=True, text=True, timeout=60)
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)

    @unittest.skipUnless(shutil.which('Rscript'), 'Rscript unavailable')
    def test_adapter_contract_with_injected_engine(self):
        root = Path(__file__).resolve().parents[1]
        completed = subprocess.run(['Rscript', 'tests/test_facets_adapter.R'], cwd=root,
                                   capture_output=True, text=True, timeout=60)
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        self.assertIn('native_engine_tested=FALSE', completed.stdout)
