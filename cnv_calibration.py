"""Descriptive calibration accounting, not automatic LoD or biological validation.

Outcomes must be adjudicated upstream under a frozen matching/review protocol.
This module verifies declared sample roles and optional artifact bytes, not the
truth of an orthogonal assay or the adequacy of a simulation model.
"""
import math
import re
from cnv_contract import validate_manifest,validate_artifact


def summarize_calibration(manifest,plan,outcomes,*,check_artifact_files):
    if type(check_artifact_files) is not bool:
        raise ValueError('explicit artifact verification policy required')
    validate_manifest(manifest,stage='reference_build')
    if (set(plan)!={'protocol_sha256','cases'} or not isinstance(plan['protocol_sha256'],str) or
            not re.fullmatch('[0-9a-f]{64}',plan['protocol_sha256']) or not plan['cases']):
        raise ValueError('nonempty frozen calibration plan and protocol digest required')
    samples={s['sample_id']:s for s in manifest['samples']}
    cases={}
    fields={'case_id','source_sample_id','truth_source','event_class','clone_fraction',
            'size_band','depth_band','informative_loci_band','truth_artifact'}
    for case in plan['cases']:
        if set(case)!=fields:
            raise ValueError('exact calibration case schema required')
        cid=case['case_id']
        if not isinstance(cid,str) or not cid or cid in cases:
            raise ValueError('unique nonempty calibration case ID required')
        sample=samples.get(case['source_sample_id'])
        if sample is None or sample['role']!='calibration':
            raise ValueError('calibration cannot use study or PoN donors')
        if case['truth_source'] not in ('SIMULATED','ORTHOGONAL'):
            raise ValueError('simulation and orthogonal evidence must be explicit')
        if case['event_class'] not in ('LOSS','GAIN','CNLOH','NULL_CONTROL'):
            raise ValueError('unsupported calibration event class')
        f=case['clone_fraction']
        if case['event_class']=='NULL_CONTROL':
            if f is not None:
                raise ValueError('null-control clone fraction must be unavailable, not fabricated')
        elif f is None and case['truth_source']=='ORTHOGONAL':
            pass  # Confirmed event need not have a measured clonal fraction.
        elif type(f) not in (int,float) or not math.isfinite(f) or not 0<f<=1:
            raise ValueError('positive calibration needs an explicit clone fraction')
        for key in ('size_band','depth_band','informative_loci_band'):
            if not isinstance(case[key],str) or not case[key]:
                raise ValueError('explicit prespecified calibration bands required')
        validate_artifact(case['truth_artifact'],check_artifact_files)
        cases[cid]=case
    observed={}
    for outcome in outcomes:
        cid=outcome['case_id']
        if cid not in cases or cid in observed:
            raise ValueError('foreign or duplicate calibration outcome')
        if outcome.get('protocol_sha256')!=plan['protocol_sha256']:
            raise ValueError('outcome uses a different calibration protocol')
        allowed=({'FALSE_POSITIVE','CLEAN'} if cases[cid]['event_class']=='NULL_CONTROL'
                 else {'DETECTED','MISSED'}) | {'NOT_EVALUABLE','NOT_RUN'}
        if outcome.get('status') not in allowed:
            raise ValueError('outcome incompatible with calibration truth class')
        if outcome['status']!='NOT_RUN':
            validate_artifact(outcome.get('result_artifact'),check_artifact_files)
        observed[cid]=outcome
    cells={}
    for cid,case in cases.items():
        sample=samples[case['source_sample_id']]
        key=(sample['stratum'],case['truth_source'],case['event_class'],case['clone_fraction'],
             case['size_band'],case['depth_band'],case['informative_loci_band'])
        cell=cells.setdefault(key,dict(counts={},sources=set(),case_ids=[]))
        status=observed.get(cid,{'status':'NOT_RUN'})['status']
        cell['counts'][status]=cell['counts'].get(status,0)+1
        cell['sources'].add(case['source_sample_id'])
        cell['case_ids'].append(cid)
    result=[]
    for key,cell in cells.items():
        counts=cell['counts']
        negative=key[2]=='NULL_CONTROL'
        success='FALSE_POSITIVE' if negative else 'DETECTED'
        other='CLEAN' if negative else 'MISSED'
        n=counts.get(success,0)+counts.get(other,0)
        fraction=counts.get(success,0)/n if n else None
        result.append(dict(zip(('stratum','truth_source','event_class','clone_fraction',
            'size_band','depth_band','informative_loci_band'),key),
            planned_cases=len(cell['case_ids']),unique_source_samples=len(cell['sources']),
            case_ids=sorted(cell['case_ids']),counts=counts,evaluable_cases=n,
            conditional_detection_fraction=None if negative else fraction,
            conditional_false_positive_fraction=fraction if negative else None,
            all_planned_cases_evaluable=n==len(cell['case_ids']),
            independent_replicates_assumed=False))
    return dict(status='CALIBRATION_DESCRIPTIVE_SUMMARY',protocol_sha256=plan['protocol_sha256'],
        cells=result,planned_cases=len(cases),
        unexecuted_cases=sum(c['counts'].get('NOT_RUN',0) for c in result),
        non_evaluable_cases=sum(c['counts'].get('NOT_EVALUABLE',0) for c in result),
        artifact_bytes_verified=check_artifact_files,truth_validity_verified=False,
        matching_or_review_protocol_validated=False,detection_limit_established=False,
        calibration_gate_closed=False)
