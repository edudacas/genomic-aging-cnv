"""Export a hash-pinned comparison as unclassified candidate evidence."""
import json
from pathlib import Path
from cnv_contract import validate_artifact
from cnv_execution import sha, atomic_json
from cnv_export_contract import candidate_record
from cnv_parquet import export_candidates


def export_comparison(receipt_artifact, sample_id, output, *, batch_size):
    validate_artifact(receipt_artifact,True)
    path=Path(receipt_artifact['path']).resolve()
    root=path.parent
    receipt=json.loads(path.read_text())
    if (path.name!='COMPLETE.json' or (root/'FAILURE.json').exists() or
            receipt.get('status')!='CALLER_COMPARISON_COMPLETE' or
            receipt.get('started_sha256')!=sha(root/'STARTED.json')):
        raise ValueError('invalid completed comparison')
    started=json.loads((root/'STARTED.json').read_text())
    if (receipt['sample_provenance']!=started['sample_provenance'] or
            receipt['sample_provenance']['sample_id']!=sample_id):
        raise ValueError('comparison sample identity mismatch')
    table=root/'comparison.jsonl'
    if set(receipt['outputs'])!={str(table)}:
        raise ValueError('unexpected comparison output registry')
    metadata=receipt['outputs'][str(table)]
    if sha(table)!=metadata['sha256'] or table.stat().st_size!=metadata['bytes']:
        raise ValueError('comparison output changed')
    artifacts=[receipt_artifact,dict(path=str(root/'STARTED.json'),sha256=sha(root/'STARTED.json')),
               dict(path=str(table),sha256=sha(table)),*started['inputs']]
    for item in artifacts:
        validate_artifact(item,True)
    sources={item['sha256']:item['path'] for item in artifacts}
    def records():
        counts={}
        with table.open() as stream:
            for index,line in enumerate(stream,1):
                row=json.loads(line)
                if row.get('sample_id')!=sample_id:
                    raise ValueError('comparison row belongs to another sample')
                if (row.get('event_class')!='UNCLASSIFIED' or row.get('evaluable_bp','missing') is not None or
                    any(row.get(k) is not False for k in
                        ('somatic_origin_classified','copy_neutral_loh_confirmed','sensitivity_calibrated'))):
                    raise ValueError('comparison contains unsupported biological claim')
                key=row['comparison_status']
                counts[key]=counts.get(key,0)+1
                yield candidate_record(sample_id,f'comparison:{index}','SEGMENT',row,metadata['sha256'],
                    chromosome=row['chromosome'],start_0based=row['start_0based'],
                    end_0based_exclusive=row['end_0based_exclusive'])
        if counts!=receipt['comparison_status_counts'] or sum(counts.values())!=receipt['comparison_intervals']:
            raise ValueError('comparison rows differ from completion counts')
    result=export_candidates(records(),output,source_artifacts=sources,batch_size=batch_size)
    lineage=dict(status='COMPARISON_CANDIDATE_EXPORT_PASS',sample_id=sample_id,
        comparison_receipt_sha256=receipt_artifact['sha256'],comparison_payload_preserved=True,
        rows=result['rows'],roundtrip_fields_match=result['roundtrip_fields_match'],
        candidate_export_receipt_sha256=sha(Path(output)/'COMPLETE.json'),
        scientific_validation=False,event_count_estimated=False,deidentification_performed=False)
    atomic_json(Path(output)/'COMPARISON_EXPORT.json',lineage)
    return lineage
