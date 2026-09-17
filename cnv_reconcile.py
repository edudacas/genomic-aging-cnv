"""Reconcile caller segment geometry and directions, not biological event calls.

Intervals are split at every boundary from either method. Original IDs and
segment-level signals remain attached; no marker counts are redistributed to
subintervals. A near-zero logR band is NOT proof of copy neutrality or cnLOH.
"""
import math
import csv
from cnv_territory import Territory


def read_segments(path, caller):
    """Read native CNVkit CNS or adapted FACETS segment TSV, retaining fields.

    Caller output hashes/receipts must be verified by orchestration before use.
    This reader validates schema/types, not source provenance or event validity.
    """
    schemas={
        'cnvkit':('chromosome','start','end','log2'),
        'facets':('chromosome','start_0based','end_0based_exclusive','mean_log2_ratio'),
    }
    if caller not in schemas:
        raise ValueError('unsupported segment source')
    chrom,start,end,ratio=schemas[caller]
    with open(path,newline='') as stream:
        reader=csv.DictReader(stream,delimiter='\t')
        fields=reader.fieldnames or []
        if len(fields)!=len(set(fields)) or not set(schemas[caller])<=set(fields):
            raise ValueError('segment source schema mismatch')
        for index,row in enumerate(reader,1):
            if None in row or any(v is None for v in row.values()):
                raise ValueError('ragged segment source')
            value=None if row[ratio]=='NA' else float(row[ratio])
            if value is not None and not math.isfinite(value):
                raise ValueError('nonfinite segment signal; missing must be explicit NA')
            yield dict(segment_id=f'{caller}:{index}',chromosome=row[chrom],
                start_0based=int(row[start]),end_0based_exclusive=int(row[end]),
                mean_log2_ratio=value,source_native_fields=row,
                source_receipt_verified=False)


def signal_direction(value, threshold):
    if value is None:
        return 'NOT_EVALUATED'
    if type(value) not in (int,float) or not math.isfinite(value):
        raise ValueError('segment log2 ratio must be finite or explicitly missing')
    if value < -threshold:
        return 'LOSS_DIRECTION'
    if value > threshold:
        return 'GAIN_DIRECTION'
    return 'WITHIN_DECLARED_LOGR_BAND'


def reconcile_segments(cnvkit, facets, lengths, *, logr_direction_threshold):
    if (type(logr_direction_threshold) not in (int,float) or
            not math.isfinite(logr_direction_threshold) or logr_direction_threshold<=0):
        raise ValueError('explicit positive finite direction threshold required')
    reference = Territory([],lengths)
    events = {chrom:{} for chrom in lengths}
    for caller, rows in (('cnvkit',cnvkit),('facets',facets)):
        seen = set()
        # Validate within-caller overlap separately, without assuming input order.
        previous = {}
        for row in rows:
            sid = row['segment_id']
            if not isinstance(sid,str) or not sid or sid in seen:
                raise ValueError('unique nonempty source segment IDs required per caller')
            seen.add(sid)
            chrom,start,end = row['chromosome'],row['start_0based'],row['end_0based_exclusive']
            reference.validate(chrom,start,end)
            if chrom in previous and start<previous[chrom]:
                raise ValueError('source segments overlap or are unsorted within chromosome')
            previous[chrom]=end
            direction=signal_direction(row['mean_log2_ratio'],logr_direction_threshold)
            source=dict(row,direction=direction)
            events[chrom].setdefault(start,[]).append((caller,source))
            events[chrom].setdefault(end,[]).append((caller,None))
    for chrom,chrom_events in events.items():
        boundaries=sorted(chrom_events)
        active={'cnvkit':None,'facets':None}
        for index,start in enumerate(boundaries[:-1]):
            # Remove ended segments before adding adjacent segments at the same base.
            for caller,row in sorted(chrom_events[start],key=lambda event:event[1] is not None):
                active[caller]=row
            end=boundaries[index+1]
            a,b=active['cnvkit'],active['facets']
            if a is None and b is None:
                continue
            da=a['direction'] if a is not None else 'NOT_EVALUATED'
            db=b['direction'] if b is not None else 'NOT_EVALUATED'
            if a is None:
                agreement='FACETS_ONLY'
            elif b is None:
                agreement='CNVKIT_ONLY'
            elif 'NOT_EVALUATED' in (da,db):
                agreement='DIRECTION_NOT_COMPARABLE'
            elif da==db:
                agreement='CONCORDANT_DIRECTION' if da!='WITHIN_DECLARED_LOGR_BAND' else 'BOTH_WITHIN_DECLARED_LOGR_BAND'
            else:
                agreement='DISCORDANT_DIRECTION'
            yield dict(chromosome=chrom,start_0based=start,end_0based_exclusive=end,
                cnvkit_source=a,facets_source=b,comparison_status=agreement,
                logr_direction_threshold=logr_direction_threshold,
                signal_scope='ORIGINAL_SOURCE_SEGMENT_NOT_REESTIMATED_IN_SUBINTERVAL',
                event_class='UNCLASSIFIED',evaluable_bp=None,
                copy_neutral_loh_confirmed=False,somatic_origin_classified=False,
                sensitivity_calibrated=False)
