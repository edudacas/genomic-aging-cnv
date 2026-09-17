"""Prepare adapted FACETS inputs. No R invocation, package install or cohort run."""
import csv
import json
import math
from pathlib import Path
from cnv_execution import atomic_json, sha
from cnv_joint_inputs import verify_joint_inputs


NUMERIC_FIELDS = ('ndepth','ndepthmax','het_threshold','snp_neighborhood','initial_cval',
                  'final_cval','gamma','min_nhet','seed','sample_het_min_depth','delta_cn')


def validate_fit_policy(policy):
    if set(policy)!=set(NUMERIC_FIELDS)|{'het_scale'}:
        raise ValueError('exact explicit fit policy required; estimator/map come from PoN')
    for field in NUMERIC_FIELDS:
        value=policy[field]
        if type(value) not in (int,float) or not math.isfinite(value):
            raise ValueError('invalid numeric fit policy: '+field)
    if type(policy['het_scale']) is not bool:
        raise ValueError('explicit het_scale boolean required')
    if not (0<policy['ndepth']<policy['ndepthmax'] and 0<policy['het_threshold']<.5
            and policy['snp_neighborhood']>=0 and 0<policy['initial_cval']<=policy['final_cval']
            and policy['gamma']>0 and policy['min_nhet']>=1 and policy['delta_cn']>=0
            and 0<=policy['seed']<=2147483647 and policy['seed']==int(policy['seed'])
            and policy['sample_het_min_depth']==50):
        raise ValueError('fit policy differs from supported native parameter ranges')


def prepare_joint(manifest_path,pair_directory,pair_config,panel_directory,panel_config,
                  fit_policy_path,output):
    paths=[Path(p).resolve() for p in (manifest_path,pair_config,panel_config,fit_policy_path)]
    for directory in (pair_directory,panel_directory):
        paths.extend(Path(directory).resolve()/name for name in ('STARTED.json','COMPLETE.json'))
    frozen={str(p):sha(p) for p in paths}
    verified=verify_joint_inputs(manifest_path,pair_directory,pair_config,panel_directory,panel_config)
    policy=json.loads(Path(fit_policy_path).read_text())
    validate_fit_policy(policy)
    if abs(verified['allelic_pseudocount']-1/6)>8*2.220446049250313e-16:
        raise ValueError('training estimator differs from audited native primitive')
    if verified['native_rows']==0:
        raise ValueError('no numeric pair rows; joint fit is not evaluable')
    if policy['het_scale'] and verified['native_bias_rows']==0:
        raise ValueError('no trained allelic markers; separate CNV-only review required')
    for path in (verified['pair_matrix'],verified['allelic_bias']):
        frozen[path]=sha(path)
    source=Path(__file__).resolve().parent
    implementation={str(p):sha(p) for p in [*source.glob('cnv_*.py'),*source.glob('facets_*.R')]}
    output=Path(output).resolve()
    output.mkdir(parents=True,exist_ok=False)
    atomic_json(output/'STARTED.json',dict(inputs=frozen,implementation=implementation))
    try:
        config=dict(policy,reference_build='GRCh38',allelic_pseudocount=1/6)
        with (output/'fit_config.tsv').open('x',newline='') as stream:
            writer=csv.writer(stream,delimiter='\t')
            writer.writerow(['key','value'])
            for key,value in sorted(config.items()):
                writer.writerow([key,str(value).upper() if type(value) is bool else value])
        with (output/'chromosome_map.tsv').open('x',newline='') as stream:
            writer=csv.writer(stream,delimiter='\t')
            writer.writerow(['source','native'])
            for chrom,code in sorted(verified['native_chromosome_map'].items(),key=lambda item:item[1]):
                writer.writerow([chrom,{23:'X',24:'Y'}.get(code,str(code))])
        if any(sha(p)!=d for p,d in (frozen|implementation).items()):
            raise ValueError('inputs changed during joint preparation')
        receipt=dict(status='JOINT_FIT_PREPARED_NOT_LAUNCHED',verified_inputs=verified,
            inputs=frozen,implementation=implementation,started_sha256=sha(output/'STARTED.json'),
            outputs={name:sha(output/name) for name in ('fit_config.tsv','chromosome_map.tsv')},
            facets_ch_executed=False,runtime_validated=False,sensitivity_calibrated=False)
        atomic_json(output/'COMPLETE.json',receipt)
        return receipt
    except Exception as error:
        atomic_json(output/'FAILURE.json',dict(error=type(error).__name__,message=str(error)))
        raise
