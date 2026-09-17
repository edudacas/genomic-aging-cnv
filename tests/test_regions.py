import unittest
from cnv_regions import regional_evidence
from cnv_territory import Territory


def seg(name,start,end):
    return dict(segment_id=name,chromosome='chr17',start_0based=start,end_0based_exclusive=end)


class RegionalTests(unittest.TestCase):
    def setUp(self):
        self.assay=Territory([('chr17',10,20),('chr17',40,50)],{'chr17':100})
        # Synthetic coordinates, never real TP53/cytoband definitions.
        self.regions={'gene_fixture':[('chr17',5,25)],'arm_fixture':[('chr17',0,60)]}

    def test_gene_and_arm_overlap_are_separate_not_additive(self):
        rows=list(regional_evidence([seg('a',0,60)],self.regions,self.assay))
        self.assertEqual(rows[0]['assay_overlap_bp'],10)
        self.assertEqual(rows[1]['assay_overlap_bp'],20)
        self.assertEqual(rows[1]['segment_span_union_overlap_bp'],60)
        self.assertIsNone(rows[1]['evaluable_bp'])
        self.assertIsNone(rows[1]['event_count'])

    def test_overlap_union_not_multiple_event_bases(self):
        row=list(regional_evidence([seg('a',0,30),seg('b',20,60)],self.regions,self.assay))[1]
        self.assertEqual(row['segment_span_union_overlap_bp'],60)
        self.assertEqual(len(row['source_segments']),2)
        self.assertEqual(row['assay_overlap_with_segment_spans_bp'],20)

    def test_no_segment_is_not_negative(self):
        row=next(regional_evidence([],self.regions,self.assay))
        self.assertEqual(row['status'],'NO_SEGMENT_EVIDENCE')
        self.assertFalse(row['negative_call_assigned'])
        self.assertIsNone(row['coverage_eligible_overlap_bp'])

    def test_empty_coverage_is_different_from_unknown(self):
        empty=Territory([],self.assay.lengths)
        row=next(regional_evidence([seg('a',0,60)],self.regions,self.assay,coverage_eligible=empty))
        self.assertEqual(row['coverage_eligible_overlap_bp'],0)
        self.assertEqual(row['status'],'NO_COVERAGE_ELIGIBLE_REGION_BASES')

    def test_outside_assay_retains_source_intersection(self):
        row=next(regional_evidence([seg('a',70,90)],{'outside':[('chr17',70,90)]},self.assay))
        self.assertEqual(row['status'],'REGION_OUTSIDE_ASSAY')
        self.assertEqual(row['segment_span_union_overlap_bp'],20)
        self.assertEqual(row['assay_overlap_with_segment_spans_bp'],0)

    def test_empty_regions_or_duplicate_sources_rejected(self):
        for segments,regions in (([],{}),([],{'x':[]}),([seg('a',0,10),seg('a',10,20)],self.regions)):
            with self.assertRaises(ValueError):
                list(regional_evidence(segments,regions,self.assay))
