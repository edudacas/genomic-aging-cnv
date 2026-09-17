"""Execute one explicitly prepared native joint fit; never launches a cohort."""
import json
from pathlib import Path
from cnv_commands import joint_fit
from cnv_contract import validate_artifact
from cnv_execution import atomic_json,execute,sha
from cnv_r_library import verify_library
from cnv_joint_output import validate_joint_output


def verify_preparation(directory):
    directory=Path(directory).resolve()
    receipt=json.loads((directory/'COMPLETE.json').read_text())
    if ((directory/'FAILURE.json').exists() or receipt.get('status')!='JOINT_FIT_PREPARED_NOT_LAUNCHED'
            or receipt.get('started_sha256')!=sha(directory/'STARTED.json')):
        raise ValueError('invalid joint preparation receipt')
    started=json.loads((directory/'STARTED.json').read_text())
    if any(started[key]!=receipt[key] for key in ('inputs','implementation')):
        raise ValueError('joint preparation bindings differ')
    for path,digest in (receipt['inputs']|receipt['implementation']).items():
        if sha(path)!=digest:
            raise ValueError('prepared joint inputs or implementation changed')
    if set(receipt['outputs'])!={'fit_config.tsv','chromosome_map.tsv'}:
        raise ValueError('unexpected joint preparation output set')
    for name,digest in receipt['outputs'].items():
        if sha(directory/name)!=digest:
            raise ValueError('prepared fit configuration changed')
    return receipt


def run_joint(prepared,library_snapshot_path,rscript,output,*,timeout_seconds):
    prepared=Path(prepared).resolve()
    snapshot_path=Path(library_snapshot_path).resolve()
    receipt=verify_preparation(prepared)
    validate_artifact(rscript,True)
    snapshot=json.loads(snapshot_path.read_text())
    if not {'facetsCH','pctGCdata'}<=set(snapshot['required_packages']):
        raise ValueError('runtime snapshot must require facetsCH and pctGCdata')
    verify_library(snapshot)
    frozen_snapshot=sha(snapshot_path)
    source=Path(__file__).resolve().parent
    output=Path(output).resolve()
    output.mkdir(parents=True,exist_ok=False)
    atomic_json(output/'STARTED.json',dict(preparation_sha256=sha(prepared/'COMPLETE.json'),
        library_snapshot_sha256=frozen_snapshot,rscript=rscript))
    try:
        result_dir=output/'native'/'result'
        command=joint_fit(receipt['verified_inputs']['pair_matrix'],
            receipt['verified_inputs']['allelic_bias'],prepared/'fit_config.tsv',
            prepared/'chromosome_map.tsv',result_dir,rscript=rscript['path'],
            wrapper=source/'facets_run.R',library=snapshot['directory'])
        implementation=[*receipt['implementation'],rscript['path'],
            *(Path(snapshot['directory'])/name for name in snapshot['files'])]
        inputs=[prepared/'STARTED.json',prepared/'COMPLETE.json',snapshot_path,
                *(prepared/name for name in receipt['outputs']),*receipt['inputs']]
        native=execute(command,output/'native',inputs,implementation_files=implementation,
                       timeout_seconds=timeout_seconds)
        expected=['JOINT_SEGMENTATION_CANDIDATES','sensitivity_calibrated=FALSE',
                  'somatic_origin_classified=FALSE','copy_neutral_loh_confirmed=FALSE']
        if (result_dir/'FAILURE.txt').exists() or (result_dir/'COMPLETE.txt').read_text().splitlines()!=expected:
            raise ValueError('native completion contract failed')
        output_validation=validate_joint_output(result_dir,receipt['verified_inputs']['native_rows'],
                                               input_matrix=receipt['verified_inputs']['pair_matrix'])
        verify_preparation(prepared)
        verify_library(snapshot)
        if sha(snapshot_path)!=frozen_snapshot:
            raise ValueError('runtime snapshot changed during execution')
        result=dict(status='JOINT_NATIVE_EXECUTION_COMPLETE',native_receipt=native,output_validation=output_validation,
            sample_provenance={key:receipt['verified_inputs'][key] for key in
                ('sample_id','stratum','manifest_sha256','sample_bam_sha256','reference_build')},
            prepared_receipt_sha256=sha(prepared/'COMPLETE.json'),library_snapshot_sha256=frozen_snapshot,
            started_sha256=sha(output/'STARTED.json'),outputs=native['outputs'],
            sensitivity_calibrated=False,somatic_origin_classified=False,
            copy_neutral_loh_confirmed=False,system_runtime_dependency_closure_verified=False)
        atomic_json(output/'COMPLETE.json',result)
        return result
    except Exception as error:
        atomic_json(output/'FAILURE.json',dict(error=type(error).__name__,message=str(error)))
        raise
