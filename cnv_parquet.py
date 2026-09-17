"""Typed candidate Parquet with exact canonical JSONL round-trip verification.

PyArrow is optional until invocation; no automatic installation. These private
evidence exports are not adjudicated calls and are not de-identification tools.
API reference: https://arrow.apache.org/docs/python/parquet.html
"""
import hashlib
import json
from pathlib import Path
from cnv_execution import atomic_json,sha
from cnv_export_contract import FIELDS,validated_candidates,validate_candidate


def canonical(row):
    return json.dumps(row,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n'


def candidate_schema(pa):
    integers={'schema_version','start_0based','end_0based_exclusive','evaluable_bp'}
    booleans={'somatic_origin_classified','copy_neutral_loh_confirmed','sensitivity_calibrated'}
    nullable={'chromosome','start_0based','end_0based_exclusive','evaluable_bp'}
    return pa.schema([pa.field(name,pa.int64() if name in integers else
        pa.bool_() if name in booleans else pa.string(),nullable=name in nullable) for name in FIELDS],
        metadata={b'cnv_candidate_schema_version':b'1',b'scientific_scope':b'UNCLASSIFIED_EVIDENCE'})


def export_candidates(records,output,*,source_artifacts,batch_size):
    if type(batch_size) is not int or batch_size<1:
        raise ValueError('explicit positive batch size required')
    # Fail before creating output if the optional runtime is absent.
    import pyarrow as pa
    import pyarrow.parquet as pq
    if not isinstance(source_artifacts,dict) or not source_artifacts:
        raise ValueError('source digest-to-path registry required')
    for digest,path in source_artifacts.items():
        if sha(path)!=digest:
            raise ValueError('source artifact checksum mismatch')
    schema=candidate_schema(pa)
    output=Path(output).resolve()
    output.mkdir(parents=True,exist_ok=False)
    atomic_json(output/'STARTED.json',dict(stage='CANDIDATE_PARQUET_EXPORT',
        source_artifacts=source_artifacts,pyarrow_version=pa.__version__,batch_size=batch_size))
    parquet=output/'candidates.parquet'
    jsonl=output/'candidates.jsonl'
    try:
        expected=hashlib.sha256()
        count=0
        with pq.ParquetWriter(parquet,schema,compression='snappy') as writer,jsonl.open('x') as stream:
            batch=[]
            for row in validated_candidates(records):
                if row['source_sha256'] not in source_artifacts:
                    raise ValueError('candidate source missing from artifact registry')
                text=canonical(row)
                expected.update(text.encode('utf-8'))
                stream.write(text)
                batch.append(row)
                count+=1
                if len(batch)==batch_size:
                    writer.write_table(pa.Table.from_pylist(batch,schema=schema))
                    batch=[]
            if batch:
                writer.write_table(pa.Table.from_pylist(batch,schema=schema))
        observed=hashlib.sha256()
        observed_count=0
        with pq.ParquetFile(parquet) as reader:
            if not reader.schema_arrow.equals(schema,check_metadata=True):
                raise ValueError('Parquet schema changed during serialization')
            for batch in reader.iter_batches(batch_size=batch_size):
                for row in batch.to_pylist():
                    validate_candidate(row)
                    observed.update(canonical(row).encode('utf-8'))
                    observed_count+=1
        if observed_count!=count or observed.digest()!=expected.digest() or sha(jsonl)!=expected.hexdigest():
            raise ValueError('Parquet/JSONL round-trip mismatch')
        for digest,path in source_artifacts.items():
            if sha(path)!=digest:
                raise ValueError('source artifact changed during export')
        receipt=dict(status='CANDIDATE_PARQUET_EXPORT_PASS',rows=count,
            canonical_records_sha256=expected.hexdigest(),roundtrip_fields_match=True,
            source_artifact_bytes_verified=True,payload_derivation_verified=False,
            pyarrow_version=pa.__version__,outputs={p.name:sha(p) for p in (parquet,jsonl)},
            scientific_validation=False,deidentification_performed=False)
        atomic_json(output/'COMPLETE.json',receipt)
        return receipt
    except Exception as error:
        atomic_json(output/'FAILURE.json',dict(error=type(error).__name__,message=str(error)))
        raise
