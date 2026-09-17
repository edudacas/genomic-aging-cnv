"""Region-level segment evidence (e.g. supplied TP53/17p BED definitions).

No built-in gene coordinates, transcript choice or centromere boundaries.
Region definitions may overlap one another; their totals must not be summed
as independent genomic territory. Absence of a segment is not a negative call.
"""
from cnv_territory import Territory


def merged_intervals(intervals,lengths):
    rank={chrom:i for i,chrom in enumerate(lengths)}
    output=[]
    for chrom,start,end in sorted(intervals,key=lambda row:(rank[row[0]],row[1],row[2])):
        if output and output[-1][0]==chrom and start<=output[-1][2]:
            output[-1]=(chrom,output[-1][1],max(end,output[-1][2]))
        else:
            output.append((chrom,start,end))
    return output


def regional_evidence(segments,regions,assay,*,coverage_eligible=None):
    """regions maps exact resource region ID to disjoint sorted BED intervals.

    Source segments must already be normalized and provenance-verified upstream.
    Outputs retain source IDs; segment overlaps are geometrically unioned, not
    counted as multiple events. Missing per-sample coverage mask remains None.
    """
    if not isinstance(regions,dict) or not regions:
        raise ValueError('explicit named region definitions required')
    if coverage_eligible is not None:
        coverage_eligible.require_subset(assay)
    segments=list(segments)
    seen=set()
    for segment in segments:
        sid=segment['segment_id']
        if not isinstance(sid,str) or not sid or sid in seen:
            raise ValueError('unique segment IDs required within regional input')
        seen.add(sid)
        assay.validate(segment['chromosome'],segment['start_0based'],segment['end_0based_exclusive'])
    for name,intervals in regions.items():
        if not isinstance(name,str) or not name:
            raise ValueError('region identifier required')
        territory=Territory(intervals,assay.lengths)
        if not territory.intervals:
            raise ValueError('empty region definition')
        intersections=[]
        hits=[]
        for segment in segments:
            clipped=[]
            for chrom,start,end in territory.intervals:
                if chrom==segment['chromosome']:
                    a=max(start,segment['start_0based'])
                    b=min(end,segment['end_0based_exclusive'])
                    if a<b:
                        clipped.append((chrom,a,b))
            if clipped:
                hits.append(dict(segment_id=segment['segment_id'],
                    region_intersection_bp=sum(end-start for _,start,end in clipped),
                    source_segment=segment))
                intersections.extend(clipped)
        union=merged_intervals(intersections,assay.lengths)
        assay_bp=sum(assay.overlap(c,a,b) for c,a,b in territory.intervals)
        covered_bp=(sum(coverage_eligible.overlap(c,a,b) for c,a,b in territory.intervals)
                    if coverage_eligible is not None else None)
        if assay_bp==0:
            status='REGION_OUTSIDE_ASSAY'
        elif covered_bp==0:
            status='NO_COVERAGE_ELIGIBLE_REGION_BASES'
        elif not hits:
            status='NO_SEGMENT_EVIDENCE'
        else:
            status='REGIONAL_SEGMENT_EVIDENCE'
        yield dict(region_id=name,status=status,source_segments=hits,
            region_definition_bp=sum(b-a for _,a,b in territory.intervals),
            assay_overlap_bp=assay_bp,coverage_eligible_overlap_bp=covered_bp,
            segment_span_union_overlap_bp=sum(b-a for _,a,b in union),
            assay_overlap_with_segment_spans_bp=sum(assay.overlap(c,a,b) for c,a,b in union),
            event_count=None,evaluable_bp=None,event_class='UNCLASSIFIED',
            region_resource_provenance_verified=False,negative_call_assigned=False,
            somatic_origin_classified=False,copy_neutral_loh_confirmed=False)
