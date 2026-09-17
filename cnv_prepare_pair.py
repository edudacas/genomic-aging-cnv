"""Prepare verified pair count artifacts only; never launches FACETS-CH."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import time

from cnv_contract import validate_manifest,validate_artifact
from cnv_coordinates import fai_lengths,bed_intervals
from cnv_execution import sha,atomic_json
from cnv_pair_matrix import verified_pair_matrix,native_rcmat_rows
from cnv_territory import Territory

CATALOG_COLUMNS=('chromosome','position','ref','alt','kind')
NATIVE_COLUMNS=('Chromosome','Position','NOR.DP','NOR.RD','TUM.DP','TUM.RD',
                'catalog_kind','masked','catalog_row')


def catalog_rows(path):
    with Path(path).open(newline='') as stream:
        reader=csv.DictReader(stream,delimiter='\t')
        if tuple(reader.fieldnames or ())!=CATALOG_COLUMNS:
            raise ValueError('frozen catalog TSV schema mismatch')
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError('ragged frozen catalog')
            row['position']=int(row['position'])
            yield row


def prepare(manifest_path,config_path,output):
    manifest_path,config_path=Path(manifest_path).resolve(),Path(config_path).resolve()
    manifest=json.loads(manifest_path.read_text())
    config=json.loads(config_path.read_text())
    validate_manifest(manifest,stage='reference_build')
    required={'sample_id','reference_id','catalog','mask_bed','registry',
              'executable','pileup_policy','rehash_inputs'}
    if set(config)!=required:
        raise ValueError('pair preparation requires an exact explicit config')
    samples={s['sample_id']:s for s in manifest['samples']}
    if config['sample_id'] not in samples:
        raise ValueError('unknown sample')
    stratum=manifest['strata'][samples[config['sample_id']]['stratum']]
    for artifact in (config['catalog'],config['mask_bed'],stratum['reference_fai']):
        validate_artifact(artifact,True)
    lengths=fai_lengths(stratum['reference_fai']['path'])
    mask=Territory(bed_intervals(config['mask_bed']['path'],lengths,allow_empty=True),lengths)
    inputs=[manifest_path,config_path,config['catalog']['path'],config['mask_bed']['path'],
            stratum['reference_fai']['path']]
    for entry in config['registry'].values():
        inputs.extend([entry['output'],Path(entry['task_directory'])/'STARTED.json',
                       Path(entry['task_directory'])/'COMPLETE.json'])
    frozen={str(Path(path).resolve()):sha(path) for path in inputs}
    implementation={str(path):sha(path) for path in Path(__file__).resolve().parent.glob('cnv_*.py')}
    output=Path(output).resolve()
    # A new attempt only. Existing results are never overwritten or repaired
    # implicitly; native task resume remains independently hash-verified.
    output.mkdir(parents=True,exist_ok=False)
    started=time.monotonic()
    atomic_json(output/'STARTED.json',dict(stage='PAIR_PREPARATION',inputs=frozen,
                                         implementation=implementation,config=config))
    try:
        total=complete=missing=masked=0
        ledger=output/'pair.evidence.jsonl.gz'
        native=output/'pair.rcmat.tsv'
        with gzip.open(ledger,'wt') as full,native.open('x',newline='') as tab:
            writer=csv.DictWriter(tab,fieldnames=NATIVE_COLUMNS,delimiter='\t')
            writer.writeheader()
            for row in verified_pair_matrix(manifest,config['sample_id'],config['reference_id'],
                    catalog_rows(config['catalog']['path']),config['registry'],lengths,
                    mask=mask,executable=config['executable'],pileup_policy=config['pileup_policy'],
                    rehash_inputs=config['rehash_inputs']):
                full.write(json.dumps(row,sort_keys=True,allow_nan=False)+'\n')
                total+=1
                masked+=int(row['masked'])
                if row['native_input_eligible']:
                    native_row=next(native_rcmat_rows([row]))
                    native_row['masked']='TRUE' if native_row['masked'] else 'FALSE'
                    writer.writerow(native_row)
                    complete+=1
                else:
                    missing+=1
        if not total:
            raise ValueError('empty prepared evidence')
        if any(sha(path)!=digest for path,digest in (frozen|implementation).items()):
            raise ValueError('input/implementation changed during preparation')
        receipt=dict(status='PAIR_PREPARATION_COMPLETE',sample_id=config['sample_id'],
            started_sha256=sha(output/'STARTED.json'),
            reference_id=config['reference_id'],catalog_loci=total,native_rows=complete,
            incomplete_loci=missing,masked_loci=masked,
            has_complete_numeric_rows=complete>0,
            native_joint_fit_validated=False,
            outputs={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in (ledger,native)},
            upstream_input_bytes_rehashed=config['rehash_inputs'],
            facets_ch_executed=False,elapsed_seconds=time.monotonic()-started)
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
    print(json.dumps(prepare(args.manifest,args.config,args.output),sort_keys=True))
