import copy
import hashlib
from pathlib import Path
import tempfile
import unittest
from cnv_contract import RESOURCE_KEYS, DERIVED_RESOURCE_KEYS, validate_manifest


def fixture():
    artifact = {'path': '/example/resource', 'sha256': '0' * 64}
    spec = dict(reference_build='GRCh38', assay_id='A', capture_id='C',
                snp_ld_pruned=False, off_target_mode='disabled')
    spec.update({key: copy.deepcopy(artifact) for key in RESOURCE_KEYS})
    samples = [dict(sample_id=role, person_id=role, stratum='S', role=role,
                    assay_id='A', capture_id='C',
                    bam=dict(artifact, path='/example/' + role + '.bam'),
                    bam_index=dict(artifact)) for role in ('study', 'pon', 'calibration')]
    return dict(schema_version=1, reference_build='GRCh38', strata={'S': spec}, samples=samples)


class ContractTests(unittest.TestCase):
    def test_reference_construction_does_not_require_its_own_outputs(self):
        manifest = fixture()
        for key in DERIVED_RESOURCE_KEYS:
            del manifest['strata']['S'][key]
        receipt = validate_manifest(manifest, stage='reference_build')
        self.assertEqual(receipt['validation_stage'], 'reference_build')
        self.assertFalse(receipt['derived_resources_required'])
        with self.assertRaises(ValueError):
            validate_manifest(manifest, stage='analysis')
        manifest['strata']['S']['stable_snp_list'] = {'path': '/unverified'}
        with self.assertRaises(ValueError):
            validate_manifest(manifest, stage='reference_build')

    def test_reference_build_still_enforces_role_separation(self):
        manifest = fixture()
        manifest['samples'][1]['person_id'] = manifest['samples'][0]['person_id']
        with self.assertRaisesRegex(ValueError, 'person reused'):
            validate_manifest(manifest, stage='reference_build')

    def test_valid_is_not_execution_or_coordinate_validation(self):
        result = validate_manifest(fixture())
        self.assertFalse(result['cohort_executed'])
        self.assertFalse(result['coordinate_compatibility_verified'])

    def test_rejects_leakage_and_mismatched_assays(self):
        for key, value in [('person_id', 'study'), ('sample_id', 'study'),
                           ('assay_id', 'B'), ('capture_id', 'D')]:
            manifest = fixture()
            manifest['samples'][1][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                validate_manifest(manifest)

    def test_rejects_reused_bam(self):
        manifest = fixture()
        manifest['samples'][1]['bam'] = manifest['samples'][0]['bam']
        with self.assertRaises(ValueError):
            validate_manifest(manifest)

    def test_rejects_missing_heldout_or_pon(self):
        for role in ('pon', 'calibration'):
            manifest = fixture()
            manifest['samples'] = [s for s in manifest['samples'] if s['role'] != role]
            with self.subTest(role=role), self.assertRaises(ValueError):
                validate_manifest(manifest)

    def test_resource_bytes_and_build(self):
        with tempfile.TemporaryDirectory() as temp:
            manifest = fixture()
            objects = list(manifest['strata']['S'][key] for key in RESOURCE_KEYS)
            objects += [s[key] for s in manifest['samples'] for key in ('bam', 'bam_index')]
            for idx, item in enumerate(objects):
                path = Path(temp) / str(idx)
                path.write_bytes(str(idx).encode())
                item.update(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            self.assertTrue(validate_manifest(manifest, True)['files_verified'])
            path.write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError, 'checksum'):
                validate_manifest(manifest, True)
        manifest = fixture()
        manifest['strata']['S']['reference_build'] = 'GRCh37'
        with self.assertRaises(ValueError):
            validate_manifest(manifest)


if __name__ == '__main__':
    unittest.main()
