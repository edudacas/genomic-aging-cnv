"""Read-only validation of the mosaic CNV input manifest; no caller execution."""
import argparse
import hashlib
import json
from pathlib import Path
import re


ROLES = {'study', 'pon', 'calibration'}
BASE_RESOURCE_KEYS = ('reference_fasta', 'reference_fai', 'target_bed', 'dense_snp_vcf')
DERIVED_RESOURCE_KEYS = ('stable_snp_list', 'allelic_bias_table')
RESOURCE_KEYS = BASE_RESOURCE_KEYS + DERIVED_RESOURCE_KEYS


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate_artifact(item, check_files):
    require(isinstance(item, dict), 'artifact must declare path and sha256')
    require(isinstance(item.get('path'), str) and bool(item['path']), 'artifact path missing')
    digest = item.get('sha256')
    require(isinstance(digest, str) and re.fullmatch('[0-9a-f]{64}', digest),
            'artifact sha256 missing/invalid')
    if check_files:
        path = Path(item['path'])
        require(path.is_file(), 'artifact file missing: ' + str(path))
        h = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
                h.update(block)
        require(h.hexdigest() == digest, 'artifact checksum mismatch: ' + str(path))


def validate_manifest(manifest, check_files=False, *, stage='analysis'):
    require(stage in ('reference_build', 'analysis'), 'unsupported validation stage')
    require(manifest.get('schema_version') == 1, 'unsupported schema_version')
    require(manifest.get('reference_build') == 'GRCh38', 'GRCh38 declaration required')
    strata = manifest.get('strata')
    require(isinstance(strata, dict) and bool(strata), 'strata missing')
    samples = manifest.get('samples')
    require(isinstance(samples, list) and bool(samples), 'samples missing')
    people, sample_ids, bam_paths = set(), set(), set()
    counts = {key: dict.fromkeys(ROLES, 0) for key in strata}
    for sid, spec in strata.items():
        require(bool(sid) and isinstance(spec, dict), 'invalid stratum')
        require(bool(spec.get('capture_id')) and bool(spec.get('assay_id')), 'capture/assay missing')
        require(spec.get('reference_build') == 'GRCh38', 'stratum build mismatch')
        require(spec.get('snp_ld_pruned') is False, 'dense SNP evidence must not be LD-pruned')
        for key in BASE_RESOURCE_KEYS:
            validate_artifact(spec.get(key), check_files)
        for key in DERIVED_RESOURCE_KEYS:
            if stage == 'analysis' or key in spec:
                validate_artifact(spec.get(key), check_files)
        require(spec.get('off_target_mode') in ('supported', 'disabled'), 'off-target policy missing')
        if spec['off_target_mode'] == 'supported':
            validate_artifact(spec.get('antitarget_bed'), check_files)
    for sample in samples:
        require(isinstance(sample, dict), 'invalid sample row')
        for key in ('sample_id', 'person_id', 'stratum', 'role', 'capture_id', 'assay_id'):
            require(isinstance(sample.get(key), str) and bool(sample[key].strip()), key + ' missing')
        require(sample['sample_id'] not in sample_ids, 'duplicate sample_id')
        require(sample['person_id'] not in people, 'person reused across samples/roles')
        require(sample['role'] in ROLES, 'invalid participant role')
        require(sample['stratum'] in strata, 'unknown stratum')
        stratum = strata[sample['stratum']]
        for key in ('capture_id', 'assay_id'):
            require(sample[key] == stratum[key], 'sample/stratum ' + key + ' mismatch')
        for key in ('bam', 'bam_index'):
            validate_artifact(sample.get(key), check_files)
        resolved = str(Path(sample['bam']['path']).resolve())
        require(resolved not in bam_paths, 'BAM reused across samples/roles')
        sample_ids.add(sample['sample_id'])
        people.add(sample['person_id'])
        bam_paths.add(resolved)
        counts[sample['stratum']][sample['role']] += 1
    for key, count in counts.items():
        require(count['pon'] > 0, 'stratum lacks independent PoN: ' + key)
        require(count['calibration'] > 0, 'stratum lacks held-out calibration samples: ' + key)
        require(count['study'] > 0, 'stratum lacks study samples: ' + key)
    return {'status': 'INPUT_CONTRACT_PASS', 'files_verified': check_files,
            'validation_stage': stage,
            'derived_resources_required': stage == 'analysis',
            'counts': counts, 'caller_validation_performed': False,
            'coordinate_compatibility_verified': False, 'cohort_executed': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--check-files', action='store_true')
    parser.add_argument('--stage', choices=('reference_build', 'analysis'), default='analysis')
    args = parser.parse_args()
    print(json.dumps(validate_manifest(json.loads(Path(args.manifest).read_text()),
                                       args.check_files, stage=args.stage), sort_keys=True))


if __name__ == '__main__':
    main()
