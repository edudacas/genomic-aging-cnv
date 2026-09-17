"""Normalize one sample with a independently verified PoN reference.

Explicit low-level stage API, not a cohort launcher or sex-chromosome policy.
Native log ratios are candidate evidence, never copy-number or somatic calls.
"""
import json
import csv
import math
from pathlib import Path

from cnv_commands import fix
from cnv_contract import validate_manifest, validate_artifact
from cnv_coordinates import fai_lengths
from cnv_execution import execute, sha, verify_completed_task
from cnv_reference import verified_coverage, verify_reference


def validate_normalized_table(path):
    count = 0
    with Path(path).open() as stream:
        reader = csv.DictReader(stream,delimiter='\t')
        if not {'chromosome','start','end','log2','weight'} <= set(reader.fieldnames or []):
            raise ValueError('normalized table schema mismatch')
        for row in reader:
            if (None in row or any(v is None for v in row.values()) or
                    int(row['start']) < 0 or int(row['end']) <= int(row['start']) or
                    not math.isfinite(float(row['log2'])) or
                    not math.isfinite(float(row['weight'])) or float(row['weight']) < 0):
                raise ValueError('invalid normalized bin')
            count += 1
    if not count:
        raise ValueError('normalization retained no bins; not evaluable for CBS')
    return count


def normalization_inputs(manifest_path, sample_id, sample_registry_path,
                     pon_registry_path, reference_directory, output_directory, *,
                     executable, processes, min_mapq, minimum_pon_samples,
                     create_placeholder=False):
    manifest_path = Path(manifest_path).resolve()
    sample_registry_path = Path(sample_registry_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    validate_manifest(manifest, stage='reference_build')
    matches = [s for s in manifest['samples'] if s['sample_id'] == sample_id]
    if len(matches) != 1 or matches[0]['role'] not in ('study','calibration'):
        raise ValueError('normalization requires one declared study/calibration sample')
    sample = matches[0]
    spec = manifest['strata'][sample['stratum']]
    kinds = {'target':'target_bed'}
    if spec['off_target_mode'] == 'supported':
        kinds['antitarget'] = 'antitarget_bed'
    registry = json.loads(sample_registry_path.read_text())
    if set(registry) != {sample_id} or set(registry[sample_id]) != set(kinds):
        raise ValueError('sample coverage registry differs from sample/off-target mode')
    for key in ['reference_fai', *kinds.values()]:
        validate_artifact(spec[key], True)
    lengths = fai_lengths(spec['reference_fai']['path'])
    tables = {}
    inputs = [manifest_path, sample_registry_path, Path(pon_registry_path).resolve(),
              spec['reference_fai']['path']]
    for kind, bed_key in kinds.items():
        tables[kind], receipts = verified_coverage(registry[sample_id][kind],
            sample, spec[bed_key], lengths, executable=executable,
            processes=processes, min_mapq=min_mapq)
        inputs.extend([*receipts, spec[bed_key]['path']])
    reference = verify_reference(manifest_path, pon_registry_path, sample['stratum'],
        reference_directory, executable=executable, processes=processes,
        min_mapq=min_mapq, minimum_pon_samples=minimum_pon_samples)
    reference_directory = Path(reference_directory).resolve()
    inputs.extend([reference['reference'], reference_directory/'STARTED.json',
                   reference_directory/'COMPLETE.json'])
    output_directory = Path(output_directory).resolve()
    # CNVkit requires a second coverage argument even for target-only analysis.
    # A header-only table explicitly means off-target disabled, not zero depth.
    if 'antitarget' not in tables:
        if create_placeholder:
            output_directory.mkdir(parents=True, exist_ok=True)
        empty = output_directory/'off_target_disabled.cnn'
        expected = 'chromosome\tstart\tend\tgene\tlog2\tdepth\n'
        if empty.exists():
            if empty.read_text() != expected:
                raise ValueError('target-only placeholder changed')
        else:
            if not create_placeholder:
                raise ValueError('missing target-only placeholder in completed normalization')
            with empty.open('x') as stream:
                stream.write(expected)
        tables['antitarget'] = empty
        inputs.append(empty)
    command = fix(tables['target'], tables['antitarget'], reference['reference'],
        output_directory/'sample.cnr', executable=executable)
    metadata = dict(sample_id=sample_id,role=sample['role'],stratum=sample['stratum'],
                    manifest_sha256=sha(manifest_path),sample_bam_sha256=sample['bam']['sha256'],
                    reference_build=manifest['reference_build'],
                    reference_provenance=reference,off_target_mode=spec['off_target_mode'])
    return command, inputs, metadata


def normalize_sample(manifest_path, sample_id, sample_registry_path,
                     pon_registry_path, reference_directory, output_directory, *,
                     executable, processes, min_mapq, minimum_pon_samples,
                     implementation_files, timeout_seconds):
    command, inputs, metadata = normalization_inputs(manifest_path,sample_id,
        sample_registry_path,pon_registry_path,reference_directory,output_directory,
        executable=executable,processes=processes,min_mapq=min_mapq,
        minimum_pon_samples=minimum_pon_samples,create_placeholder=True)
    native = execute(command, output_directory, inputs,
        implementation_files=[*Path(__file__).resolve().parent.glob('cnv_*.py'),
                              *implementation_files], timeout_seconds=timeout_seconds)
    retained_bins = validate_normalized_table(command.outputs[0])
    return dict(metadata,status='CNVKIT_NORMALIZATION_COMPLETE',native_task=native,
        normalized_table=command.outputs[0],retained_bins=retained_bins,sensitivity_calibrated=False,
        upstream_bam_bytes_rehashed_this_stage=False,
        sex_chromosome_reference_policy_validated=False, somatic_origin_classified=False)


def verify_normalization(manifest_path, sample_id, sample_registry_path,
                         pon_registry_path, reference_directory, directory, *,
                         executable, processes, min_mapq, minimum_pon_samples):
    command, inputs, metadata = normalization_inputs(manifest_path,sample_id,
        sample_registry_path,pon_registry_path,reference_directory,directory,
        executable=executable,processes=processes,min_mapq=min_mapq,
        minimum_pon_samples=minimum_pon_samples)
    native = verify_completed_task(command,directory,
        {str(Path(p).resolve()):sha(p) for p in inputs},rehash_inputs=True)
    retained_bins = validate_normalized_table(command.outputs[0])
    return dict(metadata,status='CNVKIT_NORMALIZATION_VERIFIED',native_task=native,
        normalized_table=command.outputs[0],retained_bins=retained_bins,inputs=inputs,
        upstream_bam_bytes_rehashed_this_stage=False,
        sex_chromosome_reference_policy_validated=False)
