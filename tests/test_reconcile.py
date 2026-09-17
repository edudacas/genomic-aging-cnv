import unittest
from pathlib import Path
import tempfile
from cnv_reconcile import reconcile_segments,read_segments


def segment(sid,start,end,ratio):
    return dict(segment_id=sid,chromosome='chr1',start_0based=start,
                end_0based_exclusive=end,mean_log2_ratio=ratio)


class ReconcileTests(unittest.TestCase):
    def test_native_table_conversion_keeps_coordinates_and_missing(self):
        with tempfile.TemporaryDirectory() as temp:
            cnvkit=Path(temp)/'sample.cns'
            facets=Path(temp)/'segments.tsv'
            cnvkit.write_text('chromosome\tstart\tend\tlog2\tprobes\nchr1\t9\t20\t-0.3\t10\n')
            facets.write_text('chromosome\tstart_0based\tend_0based_exclusive\tmean_log2_ratio\tn_informative_allelic\nchr1\t9\t20\tNA\t0\n')
            rows=self.rows(read_segments(cnvkit,'cnvkit'),read_segments(facets,'facets'))
            self.assertEqual(rows[0]['start_0based'],9)
            self.assertEqual(rows[0]['comparison_status'],'DIRECTION_NOT_COMPARABLE')
            self.assertEqual(rows[0]['cnvkit_source']['source_native_fields']['probes'],'10')
            self.assertFalse(rows[0]['facets_source']['source_receipt_verified'])

    def rows(self,a,b):
        return list(reconcile_segments(a,b,{'chr1':100},logr_direction_threshold=.1))

    def test_boundaries_and_opposite_direction(self):
        rows=self.rows([segment('a',0,20,-.4)],[segment('b',10,30,.3)])
        self.assertEqual([(r['start_0based'],r['end_0based_exclusive']) for r in rows],[(0,10),(10,20),(20,30)])
        self.assertEqual([r['comparison_status'] for r in rows],['CNVKIT_ONLY','DISCORDANT_DIRECTION','FACETS_ONLY'])
        self.assertEqual(rows[1]['cnvkit_source']['start_0based'],0)
        self.assertFalse(rows[1]['somatic_origin_classified'])

    def test_adjacent_segments_end_before_start(self):
        rows=self.rows([segment('a',0,10,-.4),segment('b',10,20,.3)],[segment('f',0,20,-.3)])
        self.assertEqual(rows[0]['comparison_status'],'CONCORDANT_DIRECTION')
        self.assertEqual(rows[1]['cnvkit_source']['segment_id'],'b')
        self.assertEqual(rows[1]['comparison_status'],'DISCORDANT_DIRECTION')

    def test_gaps_not_filled_and_missing_not_zero(self):
        rows=self.rows([segment('a',0,10,None)],[segment('f',20,30,0)])
        self.assertEqual(len(rows),2)
        self.assertIsNone(rows[0]['cnvkit_source']['mean_log2_ratio'])
        self.assertIsNone(rows[0]['facets_source'])

    def test_near_zero_never_confirms_cnloh(self):
        rows=self.rows([segment('a',0,10,0)],[segment('f',0,10,.01)])
        self.assertEqual(rows[0]['comparison_status'],'BOTH_WITHIN_DECLARED_LOGR_BAND')
        self.assertEqual(rows[0]['event_class'],'UNCLASSIFIED')
        self.assertFalse(rows[0]['copy_neutral_loh_confirmed'])
        self.assertIsNone(rows[0]['evaluable_bp'])

    def test_overlaps_duplicates_and_nonfinite_rejected(self):
        for bad in ([segment('a',0,20,0),segment('b',10,30,0)],
                    [segment('a',0,10,0),segment('a',10,20,0)],
                    [segment('a',0,10,float('nan'))]):
            with self.assertRaises(ValueError):
                self.rows(bad,[])

    def test_basewise_source_membership_oracle(self):
        a=[segment('a',0,10,-.2),segment('b',15,60,.2)]
        b=[segment('c',3,22,-.2),segment('d',40,70,0)]
        rows=self.rows(a,b)
        for pos in range(100):
            expected=[next((s['segment_id'] for s in source if s['start_0based']<=pos<s['end_0based_exclusive']),None) for source in (a,b)]
            matches=[r for r in rows if r['start_0based']<=pos<r['end_0based_exclusive']]
            self.assertEqual(len(matches),int(any(s is not None for s in expected)))
            if matches:
                actual=[matches[0][key]['segment_id'] if matches[0][key] else None for key in ('cnvkit_source','facets_source')]
                self.assertEqual(actual,expected)
