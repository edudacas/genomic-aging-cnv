import unittest
from cnv_burden import reviewed_burden
from cnv_territory import Territory


def event(name,intervals,**kwargs):
    return dict(dict(sample_id='s',event_id=name,intervals=intervals,event_class='LOSS',
                     review_decision='ACCEPTED',somatic_origin_reviewed=True),**kwargs)


class BurdenTests(unittest.TestCase):
    lengths={'chr1':100}
    def run_summary(self,events,status='COMPLETE',territory=None):
        return reviewed_burden('s',events,review_status=status,lengths=self.lengths,
                               evaluable_territory=territory)

    def test_no_evaluation_is_not_zero(self):
        self.assertIsNone(self.run_summary([])['event_count'])
        empty=Territory([],self.lengths)
        self.assertIsNone(self.run_summary([],territory=empty)['event_count'])

    def test_review_complete_no_events_has_descriptive_zero(self):
        result=self.run_summary([],territory=Territory([('chr1',0,20)],self.lengths))
        self.assertEqual(result['event_count'],0)
        self.assertEqual(result['altered_evaluable_union_bp'],0)
        self.assertFalse(result['inferential_analysis_authorized'])

    def test_event_count_is_not_segment_count_or_base_sum(self):
        events=[event('a',[('chr1',0,15),('chr1',10,25)]),event('b',[('chr1',20,40)])]
        territory=Territory([('chr1',0,10),('chr1',20,30)],self.lengths)
        result=self.run_summary(events,territory=territory)
        self.assertEqual(result['event_count'],2)
        self.assertEqual(result['altered_evaluable_union_bp'],20)

    def test_hold_prevents_complete_endpoint(self):
        result=self.run_summary([event('a',[('chr1',0,10)],review_decision='HOLD')],
            territory=Territory([('chr1',0,20)],self.lengths))
        self.assertIsNone(result['event_count'])
        self.assertEqual(result['held_event_ids'],['a'])

    def test_cnloh_requires_copy_neutral_review(self):
        with self.assertRaisesRegex(ValueError,'copy-neutral'):
            self.run_summary([event('a',[('chr1',0,10)],event_class='CNLOH')])

    def test_duplicate_event_or_foreign_sample_rejected(self):
        row=event('a',[('chr1',0,10)])
        for events in ([row,row],[dict(row,sample_id='other')]):
            with self.assertRaises(ValueError):
                self.run_summary(events)
