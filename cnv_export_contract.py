"""Candidate evidence export contract, separate from adjudicated event exports.

This envelope preserves full source evidence as canonical JSON and exposes typed
identity/coordinate fields. It does not verify source hashes against files,
de-identify payloads, assign a negative result or approve clinical release.
"""
import json
import re


EVIDENCE_TYPES={'SEGMENT','REGION','FRAGMENT_PHASE','STRUCTURAL_VARIANT'}
FIELDS=('schema_version','sample_id','record_id','evidence_type','chromosome',
        'start_0based','end_0based_exclusive','source_sha256','payload_json',
        'event_class','evaluable_bp','somatic_origin_classified',
        'copy_neutral_loh_confirmed','sensitivity_calibrated')


def validate_candidate(record):
    if set(record)!=set(FIELDS) or type(record['schema_version']) is not int or record['schema_version']!=1:
        raise ValueError('candidate export schema mismatch')
    if any(not isinstance(record[key],str) or not record[key] for key in ('sample_id','record_id')):
        raise ValueError('explicit study and record identifiers required')
    if record['evidence_type'] not in EVIDENCE_TYPES:
        raise ValueError('unknown evidence type')
    if not isinstance(record['source_sha256'],str) or not re.fullmatch('[0-9a-f]{64}',record['source_sha256']):
        raise ValueError('source artifact digest required')
    chrom,start,end=(record[key] for key in ('chromosome','start_0based','end_0based_exclusive'))
    if not (chrom is None and start is None and end is None):
        if not isinstance(chrom,str) or not chrom or type(start) is not int or type(end) is not int or not 0<=start<end:
            raise ValueError('coordinates must be a complete half-open interval or all null')
    if record['event_class']!='UNCLASSIFIED' or record['evaluable_bp'] is not None:
        raise ValueError('candidate evidence cannot supply adjudicated events/denominators')
    if any(record[key] is not False for key in ('somatic_origin_classified','copy_neutral_loh_confirmed','sensitivity_calibrated')):
        raise ValueError('candidate evidence cannot claim biological validation')
    if not isinstance(record['payload_json'],str):
        raise ValueError('canonical JSON evidence required')
    payload=json.loads(record['payload_json'])
    if not isinstance(payload,dict) or json.dumps(payload,sort_keys=True,separators=(',',':'),allow_nan=False)!=record['payload_json']:
        raise ValueError('source evidence must be canonical finite JSON object')
    return record


def candidate_record(sample_id,record_id,evidence_type,payload,source_sha256,*,
                      chromosome=None,start_0based=None,end_0based_exclusive=None):
    record=dict(schema_version=1,sample_id=sample_id,record_id=record_id,
        evidence_type=evidence_type,chromosome=chromosome,start_0based=start_0based,
        end_0based_exclusive=end_0based_exclusive,source_sha256=source_sha256,
        payload_json=json.dumps(payload,sort_keys=True,separators=(',',':'),allow_nan=False),
        event_class='UNCLASSIFIED',evaluable_bp=None,somatic_origin_classified=False,
        copy_neutral_loh_confirmed=False,sensitivity_calibrated=False)
    return validate_candidate(record)


def validated_candidates(records):
    seen=set()
    for record in records:
        validate_candidate(record)
        key=(record['sample_id'],record['evidence_type'],record['record_id'])
        if key in seen:
            raise ValueError('duplicate candidate record identity')
        seen.add(key)
        yield record
