"""Write PoN training artifacts only; no study calling or calibration claim."""
import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path

from cnv_allelic import MarkerPolicy
from cnv_contract import validate_manifest, validate_artifact
from cnv_coordinates import fai_lengths, bed_intervals
from cnv_execution import sha, atomic_json
from cnv_prepare_pair import catalog_rows
from cnv_territory import Territory
from cnv_training import verified_training_summaries


def train(manifest_path, config_path, output):
    manifest_path, config_path = Path(manifest_path).resolve(), Path(config_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    config = json.loads(config_path.read_text())
    validate_manifest(manifest, stage='reference_build')
    required = {'stratum','reference_sex','catalog','mask_bed','registry','executable',
                'pileup_policy','marker_policy','rehash_inputs','native_chromosome_map'}
    if set(config) != required:
        raise ValueError('explicit exact training config required')
    policy = MarkerPolicy(**config['marker_policy'])
    policy.validate()
    spec = manifest['strata'][config['stratum']]
    for artifact in (config['catalog'],config['mask_bed'],spec['reference_fai']):
        validate_artifact(artifact,True)
    lengths = fai_lengths(spec['reference_fai']['path'])
    mapping = config['native_chromosome_map']
    if (not isinstance(mapping,dict) or not mapping or not set(mapping)<=set(lengths) or
            any(type(v) is not int or not 1<=v<=24 for v in mapping.values()) or
            len(set(mapping.values()))!=len(mapping)):
        raise ValueError('explicit one-to-one native chromosome map required (X=23,Y=24)')
    mask = Territory(bed_intervals(config['mask_bed']['path'],lengths,allow_empty=True),lengths)
    paths = [manifest_path,config_path,config['catalog']['path'],config['mask_bed']['path'],
             spec['reference_fai']['path']]
    for entry in config['registry'].values():
        paths.extend([entry['output'],Path(entry['task_directory'])/'STARTED.json',
                      Path(entry['task_directory'])/'COMPLETE.json'])
    frozen = {str(Path(p).resolve()):sha(p) for p in paths}
    implementation = {str(p):sha(p) for p in Path(__file__).resolve().parent.glob('cnv_*.py')}
    output = Path(output).resolve()
    output.mkdir(parents=True,exist_ok=False)
    atomic_json(output/'STARTED.json',dict(stage='PON_MARKER_TRAINING',inputs=frozen,
                                         implementation=implementation,config=config))
    try:
        ledger, bias, stable = (output/name for name in
            ('training.audit.jsonl.gz','allelic_bias.tsv','stable_markers.tsv'))
        statuses = Counter()
        bias_rows = 0
        with gzip.open(ledger,'wt') as audit, bias.open('x',newline='') as bf, stable.open('x',newline='') as sf:
            bias_writer = csv.writer(bf,delimiter='\t')
            bias_writer.writerow(['snp_index','correction_bias','bias_unit','allelic_pseudocount'])
            stable_writer = csv.writer(sf,delimiter='\t')
            stable_writer.writerow(['chromosome','position','ref','alt'])
            for row in verified_training_summaries(manifest,config['stratum'],config['reference_sex'],
                    catalog_rows(config['catalog']['path']),config['registry'],lengths,
                    mask=mask,executable=config['executable'],pileup_policy=config['pileup_policy'],
                    marker_policy=policy,rehash_inputs=config['rehash_inputs']):
                statuses[row['status']] += 1
                native_chrom = mapping.get(row['chromosome'])
                export_bias = row['status']=='STABLE_TRAINING_MARKER' and native_chrom is not None and native_chrom<=23
                audit.write(json.dumps(dict(row,native_bias_exported=export_bias),sort_keys=True,allow_nan=False)+'\n')
                if row['status']=='STABLE_TRAINING_MARKER':
                    stable_writer.writerow([row[k] for k in ('chromosome','position','ref','alt')])
                if export_bias:
                    bias_writer.writerow([f"{native_chrom}_{row['position']}",row['correction_bias'],
                                          row['bias_unit'],row['allelic_pseudocount']])
                    bias_rows += 1
        if not statuses:
            raise ValueError('empty training catalog')
        if any(sha(p)!=digest for p,digest in (frozen|implementation).items()):
            raise ValueError('input/implementation changed during training')
        receipt = dict(status='PON_MARKER_TRAINING_COMPLETE',stratum=config['stratum'],
            started_sha256=sha(output/'STARTED.json'),
            reference_sex=config['reference_sex'],status_counts=dict(statuses),native_bias_rows=bias_rows,
            outputs={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in (ledger,bias,stable)},
            has_native_bias_rows=bias_rows>0,sensitivity_calibrated=False,facets_ch_executed=False,
            upstream_input_bytes_rehashed=config['rehash_inputs'])
        atomic_json(output/'COMPLETE.json',receipt)
        return receipt
    except Exception as error:
        atomic_json(output/'FAILURE.json',dict(error=type(error).__name__,message=str(error)))
        raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',required=True)
    parser.add_argument('--config',required=True)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    print(json.dumps(train(args.manifest,args.config,args.output),sort_keys=True))
