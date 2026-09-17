"""Reconcile hash-pinned completed caller stages for one exact manifest sample."""
import json
from pathlib import Path
from cnv_commands import Command
from cnv_contract import validate_artifact, validate_manifest
from cnv_coordinates import fai_lengths
from cnv_execution import atomic_json, sha, verify_completed_task
from cnv_reconcile import read_segments, reconcile_segments


def verified_source(artifact, status, provenance, relative_segment):
    validate_artifact(artifact,True)
    receipt_path=Path(artifact['path']).resolve()
    root=receipt_path.parent
    receipt=json.loads(receipt_path.read_text())
    if (receipt_path.name!='COMPLETE.json' or (root/'FAILURE.json').exists() or
            receipt.get('status')!=status or receipt.get('started_sha256')!=sha(root/'STARTED.json')):
        raise ValueError('invalid completed caller stage')
    if receipt.get('sample_provenance')!=provenance:
        raise ValueError('caller sample/BAM/manifest provenance mismatch')
    native=root/'native'
    started=json.loads((native/'STARTED.json').read_text())
    command=Command(**started['command'])
    verify_completed_task(command,native,started['inputs'],rehash_inputs=True)
    native_complete=json.loads((native/'COMPLETE.json').read_text())
    if native_complete['outputs']!=receipt['outputs']:
        raise ValueError('stage and native output registries differ')
    segment=(root/relative_segment).resolve()
    if str(segment) not in receipt['outputs']:
        raise ValueError('expected segment table not declared by caller')
    return segment


def compare_callers(manifest_artifact,sample_id,cbs_artifact,joint_artifact,output,*,
                    logr_direction_threshold):
    validate_artifact(manifest_artifact,True)
    manifest=json.loads(Path(manifest_artifact['path']).read_text())
    validate_manifest(manifest,stage='reference_build')
    samples=[s for s in manifest['samples'] if s['sample_id']==sample_id]
    if len(samples)!=1 or samples[0]['role'] not in ('study','calibration'):
        raise ValueError('one declared study/calibration sample required')
    sample=samples[0]
    provenance=dict(sample_id=sample_id,stratum=sample['stratum'],
        manifest_sha256=manifest_artifact['sha256'],sample_bam_sha256=sample['bam']['sha256'],
        reference_build=manifest['reference_build'])
    cbs=verified_source(cbs_artifact,'CNVKIT_CBS_COMPLETE',provenance,'native/sample.cns')
    facets=verified_source(joint_artifact,'JOINT_NATIVE_EXECUTION_COMPLETE',provenance,
                           'native/result/segments.tsv')
    fai=manifest['strata'][sample['stratum']]['reference_fai']
    validate_artifact(fai,True)
    lengths=fai_lengths(fai['path'])
    artifacts=[manifest_artifact,cbs_artifact,joint_artifact,fai,
               dict(path=str(cbs),sha256=sha(cbs)),dict(path=str(facets),sha256=sha(facets))]
    output=Path(output).resolve()
    output.mkdir(parents=True,exist_ok=False)
    atomic_json(output/'STARTED.json',dict(sample_provenance=provenance,inputs=artifacts,
        logr_direction_threshold=logr_direction_threshold))
    try:
        counts={}
        def sourced(path,caller):
            for row in read_segments(path,caller):
                row['source_receipt_verified']=True
                yield row
        result=output/'comparison.jsonl'
        with result.open('x') as stream:
            for row in reconcile_segments(sourced(cbs,'cnvkit'),sourced(facets,'facets'),lengths,
                                           logr_direction_threshold=logr_direction_threshold):
                row['sample_id']=sample_id
                stream.write(json.dumps(row,sort_keys=True,allow_nan=False)+'\n')
                key=row['comparison_status']
                counts[key]=counts.get(key,0)+1
        for artifact in artifacts:
            validate_artifact(artifact,True)
        receipt=dict(status='CALLER_COMPARISON_COMPLETE',sample_provenance=provenance,
            comparison_intervals=sum(counts.values()),comparison_status_counts=counts,
            outputs={str(result):dict(sha256=sha(result),bytes=result.stat().st_size)},
            started_sha256=sha(output/'STARTED.json'),event_count_estimated=False,
            sensitivity_calibrated=False,somatic_origin_classified=False)
        atomic_json(output/'COMPLETE.json',receipt)
        return receipt
    except Exception as error:
        atomic_json(output/'FAILURE.json',dict(error=type(error).__name__,message=str(error)))
        raise
