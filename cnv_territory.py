"""Exact half-open territory accounting; no calibration or event inference."""
from bisect import bisect_right


class Territory:
    """Disjoint frozen intervals with O(log n) overlap queries."""
    def __init__(self, intervals, lengths):
        if not lengths or any(type(n) is not int or n <= 0 for n in lengths.values()):
            raise ValueError('positive reference lengths required')
        self.lengths = dict(lengths)
        rank = {chrom: i for i, chrom in enumerate(lengths)}
        self.intervals = []
        by_chrom = {chrom: [] for chrom in lengths}
        previous = None
        for chrom, start, end in intervals:
            self.validate(chrom, start, end)
            if previous and (rank[chrom] < rank[previous[0]] or
                             (chrom == previous[0] and start < previous[2])):
                raise ValueError('territory unsorted or overlapping')
            previous = chrom, start, end
            self.intervals.append(previous)
            by_chrom[chrom].append((start, end))
        self.index = {}
        for chrom, rows in by_chrom.items():
            starts, ends, prefix = [], [], [0]
            for start, end in rows:
                starts.append(start)
                ends.append(end)
                prefix.append(prefix[-1]+end-start)
            self.index[chrom] = starts, ends, prefix

    def validate(self, chrom, start, end):
        if (chrom not in self.lengths or type(start) is not int or type(end) is not int or
                not 0 <= start < end <= self.lengths[chrom]):
            raise ValueError('invalid half-open reference interval')

    def overlap(self, chrom, start, end):
        self.validate(chrom, start, end)
        starts, ends, prefix = self.index[chrom]
        def cumulative(x):
            i = bisect_right(starts, x)-1
            if i < 0:
                return 0
            return prefix[i] + max(0, min(x, ends[i])-starts[i])
        return cumulative(end)-cumulative(start)

    def require_subset(self, parent):
        if self.lengths != parent.lengths:
            raise ValueError('territory reference mismatch')
        if any(parent.overlap(chrom,start,end) != end-start for chrom,start,end in self.intervals):
            raise ValueError('child territory extends outside parent')


def annotate_segment_territory(segments, assay, *, coverage_eligible=None,
                               calibrated=None):
    """Each optional territory is externally supplied, not inferred from signal.

    A provided empty mask yields measured overlap zero. An unavailable mask
    yields None. A calibrated mask must be derived by a separate validated
    calibration stage; this function checks geometry only, not its provenance.
    """
    if coverage_eligible is not None:
        coverage_eligible.require_subset(assay)
    if calibrated is not None:
        calibrated.require_subset(coverage_eligible if coverage_eligible is not None else assay)
    for segment in segments:
        chrom = segment['chromosome']
        start, end = segment['start_0based'], segment['end_0based_exclusive']
        assay.validate(chrom,start,end)
        result = dict(segment)
        result.update(marker_bounding_span_bp=end-start,
            assay_overlap_bp=assay.overlap(chrom,start,end),
            coverage_eligible_overlap_bp=(coverage_eligible.overlap(chrom,start,end)
                                         if coverage_eligible is not None else None),
            calibrated_overlap_bp=(calibrated.overlap(chrom,start,end)
                                   if calibrated is not None else None),
            territory_calibration_provenance_verified=False)
        # Do not promote geometry to a validated evaluable-burden denominator.
        result['evaluable_bp'] = None
        yield result
