"""Build CNVkit references from receipt-verified PoN coverage, never study data.

Coverage registry maps sample_id to target/antitarget task directories. Native
STARTED/COMPLETE receipts bind each table to the manifest BAM/index and BED.
This layer verifies those bindings; it does not establish PoN clinical fitness.
"""
from dataclasses import asdict
import csv
import hashlib
import json
import math
from pathlib import Path

from cnv_commands import coverage, pooled_reference, identifier, integer
from cnv_contract import validate_manifest, validate_artifact
from cnv_coordinates import fai_lengths, bed_intervals
from cnv_execution import execute, sha, verify_completed_task


def verified_coverage(directory, sample, bed, lengths, *, executable, processes, min_mapq):
    directory = Path(directory).resolve()
    started_path, complete_path = directory/'STARTED.json', directory/'COMPLETE.json'
    started, complete = (json.loads(path.read_text()) for path in (started_path, complete_path))
    frozen = {key: started[key] for key in ('command', 'inputs', 'implementation')}
    fingerprint = hashlib.sha256(json.dumps(frozen, sort_keys=True).encode()).hexdigest()
    if (started.get('fingerprint') != fingerprint or
            complete.get('fingerprint') != fingerprint or
            complete.get('status') != 'NATIVE_TASK_COMPLETE'):
        raise ValueError('coverage receipt fingerprint/status mismatch')
    outputs = frozen['command']['outputs']
    if len(outputs) != 1:
        raise ValueError('coverage must have exactly one declared table')
    table = Path(outputs[0]).resolve()
    if not table.is_relative_to(directory):
        raise ValueError('coverage table outside task directory')
    expected = coverage(sample['bam']['path'], bed['path'], table,
                        executable=executable, processes=processes, min_mapq=min_mapq)
    if json.loads(json.dumps(asdict(expected))) != frozen['command']:
        raise ValueError('coverage command/policy mismatch')
    for artifact in (sample['bam'], sample['bam_index'], bed):
        key = str(Path(artifact['path']).resolve())
        if frozen['inputs'].get(key) != artifact['sha256']:
            raise ValueError('coverage BAM/index/BED binding mismatch')
    if not frozen['implementation']:
        raise ValueError('coverage implementation provenance missing')
    for path, digest in frozen['implementation'].items():
        if sha(path) != digest:
            raise ValueError('coverage implementation changed')
    observed = {str(table): {'sha256': sha(table), 'bytes': table.stat().st_size}}
    if complete.get('outputs') != observed:
        raise ValueError('coverage output changed')
    expected_bins = bed_intervals(bed['path'], lengths)
    observed_bins = []
    with table.open() as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not set(expected.table_columns) <= set(reader.fieldnames or []):
            raise ValueError('coverage table schema mismatch')
        for row in reader:
            observed_bins.append((row['chromosome'], int(row['start']), int(row['end'])))
            depth, log2 = float(row['depth']), float(row['log2'])
            if not math.isfinite(depth) or depth < 0 or not math.isfinite(log2):
                raise ValueError('invalid coverage depth/log2')
    if observed_bins != expected_bins:
        raise ValueError('coverage bins differ from frozen BED')
    return table, [started_path, complete_path, table]


def reference_inputs(manifest_path, registry_path, stratum_id, output_directory, *,
                     executable, processes, min_mapq, minimum_pon_samples):
    """Reconstruct the exact allowed reference command from independent inputs."""
    integer(minimum_pon_samples, 'minimum_pon_samples')
    identifier(stratum_id)
    manifest_path, registry_path = Path(manifest_path).resolve(), Path(registry_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    validate_manifest(manifest, stage='reference_build')
    if stratum_id not in manifest['strata']:
        raise ValueError('unknown stratum')
    spec = manifest['strata'][stratum_id]
    registry = json.loads(registry_path.read_text())
    normals = sorted((sample for sample in manifest['samples']
                      if sample['stratum'] == stratum_id and sample['role'] == 'pon'),
                     key=lambda sample: sample['sample_id'])
    if len(normals) < minimum_pon_samples:
        raise ValueError('insufficient independent PoN samples')
    if set(registry) != {sample['sample_id'] for sample in normals}:
        raise ValueError('registry must contain exactly the stratum PoN; no study/calibration inputs')
    resources = ['reference_fasta', 'reference_fai', 'target_bed']
    if spec['off_target_mode'] == 'supported':
        resources.append('antitarget_bed')
    for key in resources:
        validate_artifact(spec[key], True)
    lengths = fai_lengths(spec['reference_fai']['path'])
    kinds = {'target': 'target_bed'}
    if spec['off_target_mode'] == 'supported':
        kinds['antitarget'] = 'antitarget_bed'
        # A base cannot contribute twice to the pooled reference.
        targets = bed_intervals(spec['target_bed']['path'], lengths)
        anti = bed_intervals(spec['antitarget_bed']['path'], lengths)
        rank = {chrom: i for i, chrom in enumerate(lengths)}
        merged = sorted(targets + anti, key=lambda x: (rank[x[0]], x[1]))
        if any(a[0] == b[0] and b[1] < a[2] for a, b in zip(merged, merged[1:])):
            raise ValueError('target/antitarget overlap')
    tables, inputs = [], [manifest_path, registry_path]
    inputs.extend(spec[key]['path'] for key in resources)
    for sample in normals:
        if set(registry[sample['sample_id']]) != set(kinds):
            raise ValueError('coverage registry does not match off-target mode')
        for kind, bed_key in kinds.items():
            table, receipt_inputs = verified_coverage(
                registry[sample['sample_id']][kind], sample, spec[bed_key], lengths,
                executable=executable, processes=processes, min_mapq=min_mapq)
            tables.append(table)
            inputs.extend(receipt_inputs)
    if len(set(tables)) != len(tables):
        raise ValueError('coverage table reused across PoN samples')
    output = Path(output_directory).resolve()/'reference.cnn'
    command = pooled_reference(tables, spec['reference_fasta']['path'], output,
                               executable=executable)
    return command, inputs, [sample['sample_id'] for sample in normals]


def build_reference(manifest_path, registry_path, stratum_id, output_directory, *,
                    executable, processes, min_mapq, minimum_pon_samples,
                    implementation_files, timeout_seconds):
    """Explicit invocation only. Caller chooses minimum PoN size before analysis."""
    command, inputs, members = reference_inputs(manifest_path, registry_path,
        stratum_id, output_directory, executable=executable, processes=processes,
        min_mapq=min_mapq, minimum_pon_samples=minimum_pon_samples)
    output = Path(command.outputs[0])
    receipt = execute(command, output.parent, inputs,
                      implementation_files=[__file__, *implementation_files],
                      timeout_seconds=timeout_seconds)
    return dict(receipt, stratum=stratum_id,
                reference_members=members,
                study_samples_used=0, calibration_samples_used=0,
                clinical_reference_suitability_verified=False)


def verify_reference(manifest_path, registry_path, stratum_id, directory, *,
                     executable, processes, min_mapq, minimum_pon_samples):
    """Read-only reference provenance check for downstream fix/segmentation.

Does not rerun a caller or infer membership from reference filename. Rechecks
the complete expected PoN coverage registry and reference task bindings.
"""
    command, inputs, members = reference_inputs(manifest_path, registry_path,
        stratum_id, directory, executable=executable, processes=processes,
        min_mapq=min_mapq, minimum_pon_samples=minimum_pon_samples)
    verified = verify_completed_task(command, directory,
        {str(Path(path).resolve()):sha(path) for path in inputs}, rehash_inputs=True)
    return dict(status='CNVKIT_REFERENCE_PROVENANCE_VERIFIED',
        reference=command.outputs[0], reference_sha256=sha(command.outputs[0]),
        stratum=stratum_id, reference_members=members, native_task=verified,
        clinical_reference_suitability_verified=False,
        sex_chromosome_reference_policy_validated=False)
