"""Descriptive aggregation of an externally reviewed event ledger.

This function consumes review decisions, never classifies caller segments.
Provenance/calibration of the review and territory must be verified separately
before using these descriptive quantities as inferential-model endpoints.
"""
from cnv_regions import merged_intervals
from cnv_territory import Territory


def reviewed_burden(sample_id,events,*,review_status,lengths,evaluable_territory):
    if not isinstance(sample_id,str) or not sample_id or review_status not in ('COMPLETE','INCOMPLETE'):
        raise ValueError('explicit sample and review status required')
    bounds=Territory([],lengths)
    if evaluable_territory is not None and evaluable_territory.lengths!=lengths:
        raise ValueError('evaluable territory reference mismatch')
    seen=set(); held=[]; rejected=[]; accepted=[]; outside=[]
    for event in events:
        if event['sample_id']!=sample_id:
            raise ValueError('foreign sample event')
        key=event['event_id']
        if not isinstance(key,str) or not key or key in seen:
            raise ValueError('duplicate or missing reviewed event ID')
        seen.add(key)
        state=event['review_decision']
        if state not in ('ACCEPTED','HOLD','REJECTED'):
            raise ValueError('explicit event review decision required')
        intervals=event['intervals']
        if not intervals:
            raise ValueError('event intervals required')
        for chrom,start,end in intervals:
            bounds.validate(chrom,start,end)
        union=merged_intervals(intervals,lengths)
        if state=='HOLD':
            held.append(key); continue
        if state=='REJECTED':
            rejected.append(key); continue
        if event['event_class'] not in ('LOSS','GAIN','CNLOH'):
            raise ValueError('accepted event requires reviewed class')
        if event.get('somatic_origin_reviewed') is not True:
            raise ValueError('accepted somatic burden event requires origin review')
        if event['event_class']=='CNLOH' and event.get('copy_neutral_evidence_reviewed') is not True:
            raise ValueError('CNLOH requires reviewed copy-neutral evidence')
        overlap=(sum(evaluable_territory.overlap(c,a,b) for c,a,b in union)
                 if evaluable_territory is not None else None)
        if overlap==0:
            outside.append(key); continue
        accepted.append(dict(event_id=key,event_class=event['event_class'],intervals=union,
                             evaluable_overlap_bp=overlap))
    denominator=(sum(b-a for _,a,b in evaluable_territory.intervals)
                 if evaluable_territory is not None else None)
    complete=review_status=='COMPLETE' and not held and denominator is not None and denominator>0
    intervals=merged_intervals([interval for event in accepted for interval in event['intervals']],lengths)
    span=(sum(evaluable_territory.overlap(c,a,b) for c,a,b in intervals)
          if evaluable_territory is not None else None)
    return dict(sample_id=sample_id,status='DESCRIPTIVE_REVIEW_COMPLETE' if complete else 'NOT_EVALUABLE_OR_REVIEW_INCOMPLETE',
        event_count=len(accepted) if complete else None,
        altered_evaluable_union_bp=span if complete else None,
        declared_evaluable_bp=denominator,
        reviewed_accepted_events=accepted,held_event_ids=held,rejected_event_ids=rejected,
        accepted_outside_scope_event_ids=outside,
        class_event_counts={kind:sum(e['event_class']==kind for e in accepted) if complete else None
                            for kind in ('LOSS','GAIN','CNLOH')},
        review_provenance_verified=False,territory_calibration_verified=False,
        inferential_analysis_authorized=False)
