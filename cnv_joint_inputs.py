"""Verify prepared training and pair artifacts before a joint fit. Read-only."""
import json
from pathlib import Path
from cnv_contract import validate_manifest
from cnv_execution import sha


STAGES = {
    'PAIR_PREPARATION_COMPLETE': ('PAIR_PREPARATION', {'pair.evidence.jsonl.gz','pair.rcmat.tsv'}),
    'PON_MARKER_TRAINING_COMPLETE': ('PON_MARKER_TRAINING',
        {'training.audit.jsonl.gz','allelic_bias.tsv','stable_markers.tsv'}),
}


def verify_prepared(directory, status, manifest_path, config_path):
    directory = Path(directory).resolve()
    stage, names = STAGES[status]
    if (directory/'FAILURE.json').exists():
        raise ValueError('failed attempt cannot supply joint inputs')
    started = json.loads((directory/'STARTED.json').read_text())
    complete = json.loads((directory/'COMPLETE.json').read_text())
    if (complete.get('status')!=status or started.get('stage')!=stage or
            complete.get('started_sha256')!=sha(directory/'STARTED.json')):
        raise ValueError('invalid prepared stage receipt')
    for path in (manifest_path,config_path):
        if started['inputs'].get(str(Path(path).resolve()))!=sha(path):
            raise ValueError('prepared manifest/config binding changed')
    config = json.loads(Path(config_path).read_text())
    if config!=started['config']:
        raise ValueError('prepared inline configuration differs')
    if not started.get('implementation') or any(sha(p)!=d for p,d in started['implementation'].items()):
        raise ValueError('prepared implementation changed')
    if set(complete.get('outputs',{}))!=names:
        raise ValueError('prepared output set differs')
    for name,metadata in complete['outputs'].items():
        path=directory/name
        if (not path.is_file() or path.resolve().parent!=directory or
                path.stat().st_size!=metadata['bytes'] or sha(path)!=metadata['sha256']):
            raise ValueError('prepared output bytes changed')
    return config, complete


def verify_joint_inputs(manifest_path,pair_directory,pair_config,panel_directory,panel_config):
    manifest=json.loads(Path(manifest_path).read_text())
    validate_manifest(manifest,stage='reference_build')
    pair, pair_receipt=verify_prepared(pair_directory,'PAIR_PREPARATION_COMPLETE',manifest_path,pair_config)
    panel, panel_receipt=verify_prepared(panel_directory,'PON_MARKER_TRAINING_COMPLETE',manifest_path,panel_config)
    samples={s['sample_id']:s for s in manifest['samples']}
    sample=samples[pair['sample_id']]
    if sample['stratum']!=panel['stratum'] or sample.get('reference_sex')!=panel['reference_sex']:
        raise ValueError('pair/panel stratum or reference sex differs')
    for field in ('catalog','mask_bed','pileup_policy','executable'):
        if pair[field]!=panel[field]:
            raise ValueError('pair/panel resource or filter mismatch: '+field)
    if pair['reference_id'] not in panel['registry']:
        raise ValueError('pair reference not in trained PoN')
    if pair['registry'][pair['reference_id']]!=panel['registry'][pair['reference_id']]:
        raise ValueError('reference task differs from training task')
    return dict(status='JOINT_INPUT_ARTIFACTS_VERIFIED',sample_id=pair['sample_id'],
        stratum=sample['stratum'],manifest_sha256=sha(manifest_path),
        sample_bam_sha256=sample['bam']['sha256'],reference_build=manifest['reference_build'],
        pair_matrix=str(Path(pair_directory).resolve()/'pair.rcmat.tsv'),
        allelic_bias=str(Path(panel_directory).resolve()/'allelic_bias.tsv'),
        native_chromosome_map=panel['native_chromosome_map'],
        allelic_pseudocount=panel['marker_policy']['allelic_pseudocount'],
        native_rows=pair_receipt['native_rows'],native_bias_rows=panel_receipt['native_bias_rows'],
        upstream_bams_rehashed_this_check=False,native_joint_fit_validated=False,
        sensitivity_calibrated=False)
