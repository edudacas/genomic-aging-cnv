"""Explicit CBS stage over a receipt-verified normalization; no event calls."""
import json
from pathlib import Path

from cnv_commands import segment
from cnv_contract import validate_artifact
from cnv_coordinates import fai_lengths
from cnv_execution import execute, atomic_json, sha
from cnv_normalize import verify_normalization
from cnv_r_library import verify_library
from cnv_reconcile import read_segments
from cnv_territory import Territory


def run_cbs(normalization_config, library_snapshot, rscript, output, *,
            processes, threshold, implementation_files, timeout_seconds):
    config_path, snapshot_path = (Path(p).resolve() for p in
                                 (normalization_config,library_snapshot))
    config = json.loads(config_path.read_text())
    required = {'manifest_path','sample_id','sample_registry_path','pon_registry_path',
        'reference_directory','directory','executable','processes','min_mapq','minimum_pon_samples'}
    if set(config) != required:
        raise ValueError('exact normalization verification config required')
    verified = verify_normalization(**config)
    validate_artifact(rscript,True)
    snapshot = json.loads(snapshot_path.read_text())
    if 'DNAcopy' not in snapshot['required_packages']:
        raise ValueError('CBS library snapshot must require DNAcopy')
    verify_library(snapshot)
    manifest = json.loads(Path(config['manifest_path']).read_text())
    fai = manifest['strata'][verified['stratum']]['reference_fai']['path']
    territory = Territory([],fai_lengths(fai))
    directory = Path(config['directory']).resolve()
    inputs = [config_path,snapshot_path,*verified['inputs'],verified['normalized_table'],
              directory/'STARTED.json',directory/'COMPLETE.json']
    frozen = {str(Path(p).resolve()):sha(p) for p in inputs}
    runtime = {k:snapshot['directory'] for k in ('R_LIBS','R_LIBS_USER','R_LIBS_SITE')}
    runtime.update({k:'/dev/null' for k in
                    ('R_PROFILE','R_PROFILE_USER','R_ENVIRON','R_ENVIRON_USER')})
    output = Path(output).resolve()
    command = segment(verified['normalized_table'],output/'native/sample.cns',
        executable=config['executable'],processes=processes,threshold=threshold,
        rscript=rscript['path'])
    output.mkdir(parents=True,exist_ok=False)
    atomic_json(output/'STARTED.json',dict(inputs=frozen,rscript=rscript,
        processes=processes,threshold=threshold,runtime_environment=runtime))
    try:
        native = execute(command,output/'native',inputs,
            implementation_files=[*Path(__file__).resolve().parent.glob('cnv_*.py'),
                rscript['path'],*implementation_files,
                *(Path(snapshot['directory'])/name for name in snapshot['files'])],
            timeout_seconds=timeout_seconds,runtime_environment=runtime)
        count = 0
        ends = {}
        for row in read_segments(command.outputs[0],'cnvkit'):
            chrom, start, end = row['chromosome'],row['start_0based'],row['end_0based_exclusive']
            territory.validate(chrom,start,end)
            if start < ends.get(chrom,0) or row['mean_log2_ratio'] is None:
                raise ValueError('overlapping/unsorted or missing-signal CBS segments')
            probes = float(row['source_native_fields']['probes'])
            if not probes.is_integer() or probes < 1:
                raise ValueError('CBS segment requires positive integer probe count')
            ends[chrom] = end
            count += 1
        if not count:
            raise ValueError('CBS produced no segments')
        verify_library(snapshot)
        verify_normalization(**config)
        if any(sha(p)!=d for p,d in frozen.items()):
            raise ValueError('CBS inputs changed during execution')
        receipt = dict(status='CNVKIT_CBS_COMPLETE',sample_id=verified['sample_id'],
            sample_provenance={key:verified[key] for key in
                ('sample_id','stratum','manifest_sha256','sample_bam_sha256','reference_build')},
            stratum=verified['stratum'],segment_rows=count,native_task=native,
            normalized_provenance_verified=True,outputs=native['outputs'],
            started_sha256=sha(output/'STARTED.json'),
            sensitivity_calibrated=False,somatic_origin_classified=False,
            sex_chromosome_reference_policy_validated=False,
            full_runtime_dependency_closure_verified=False)
        atomic_json(output/'COMPLETE.json',receipt)
        return receipt
    except Exception as error:
        atomic_json(output/'FAILURE.json',dict(error=type(error).__name__,message=str(error)))
        raise
