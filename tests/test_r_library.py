from pathlib import Path
import tempfile
import unittest
from cnv_r_library import snapshot_library,verify_library


class RLibraryTests(unittest.TestCase):
    def setUp(self):
        temp=tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root=Path(temp.name)/'library'
        (self.root/'example/Meta').mkdir(parents=True)
        (self.root/'example/DESCRIPTION').write_text('Package: example\nVersion: 1.0\n')
        (self.root/'example/Meta/package.rds').write_bytes(b'not a real R package; fingerprint fixture only')

    def snapshot(self):
        return snapshot_library(self.root,required_packages=['example'])

    def test_snapshot_does_not_claim_native_loading(self):
        result=verify_library(self.snapshot())
        self.assertEqual(result['file_count'],2)
        self.assertFalse(result['native_load_verified'])

    def test_changed_or_added_files_invalidate_snapshot(self):
        snapshot=self.snapshot()
        path=self.root/'example/extra'
        path.write_bytes(b'new')
        with self.assertRaisesRegex(ValueError,'differs'):
            verify_library(snapshot)
        snapshot=self.snapshot()
        path.write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError,'differs'):
            verify_library(snapshot)

    def test_lock_and_missing_package_rejected(self):
        with self.assertRaisesRegex(ValueError,'metadata missing'):
            snapshot_library(self.root,required_packages=['facetsCH'])
        (self.root/'00LOCK-example').mkdir()
        with self.assertRaisesRegex(ValueError,'lock present'):
            self.snapshot()

    def test_external_symlink_cannot_hide_dependency(self):
        external=self.root.parent/'external'
        external.write_text('external bytes')
        (self.root/'example/external').symlink_to(external)
        with self.assertRaisesRegex(ValueError,'external'):
            self.snapshot()

    def test_internal_file_symlink_tracked(self):
        (self.root/'example/link').symlink_to('DESCRIPTION')
        snapshot=self.snapshot()
        self.assertEqual(snapshot['files']['example/link']['link_target'],'example/DESCRIPTION')
        verify_library(snapshot)
