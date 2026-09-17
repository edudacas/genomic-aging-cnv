import unittest
from cnv_phasing import PhaseReadPolicy,fragment_evidence,summarize_fragment_evidence


def obs(name,locus,allele,**kwargs):
    return dict(dict(sample_id='s',read_group='rg',query_name=name,locus_id=locus,
                     allele=allele,mapq=60,baseq=30,flag=0),**kwargs)


class PhasingTests(unittest.TestCase):
    policy=PhaseReadPolicy(20,20,4|256|2048|1024)

    def summarize(self,rows):
        return summarize_fragment_evidence(rows,'s',['a','b'],policy=self.policy)

    def test_overlapping_mates_count_once(self):
        rows=[obs('q','a','ALT'),obs('q','a','ALT'),obs('q','b','ALT'),obs('q','b','ALT')]
        result=self.summarize(rows)
        self.assertEqual(result['linked_fragments'],1)
        self.assertEqual(result['haplotype_fragment_counts']['ALT_ALT'],1)
        self.assertIsNone(result['phase_assignment'])
        self.assertFalse(result['somatic_second_hit_confirmed'])

    def test_disagreeing_mates_not_majority_voted(self):
        rows=[obs('q','a','ALT'),obs('q','a','ALT'),obs('q','a','REF'),obs('q','b','ALT')]
        result=self.summarize(rows)
        self.assertEqual(result['linked_fragments'],0)
        self.assertEqual(result['fragment_status_counts']['AMBIGUOUS_FRAGMENT'],1)

    def test_no_link_is_not_trans(self):
        result=self.summarize([obs('q1','a','ALT'),obs('q2','b','ALT')])
        self.assertEqual(result['status'],'PHASE_NOT_EVALUATED')
        self.assertIsNone(result['phase_assignment'])

    def test_filters_and_unknown_mapq(self):
        for bad in (dict(flag=1024),dict(flag=256),dict(flag=2048),dict(mapq=255),dict(mapq=10),dict(baseq=10)):
            rows=[obs('q','a','REF'),obs('q','b','ALT',**bad)]
            with self.subTest(bad=bad):
                self.assertEqual(self.summarize(rows)['linked_fragments'],0)

    def test_read_group_prevents_template_collision(self):
        rows=[obs('q','a','ALT',read_group='r1'),obs('q','b','ALT',read_group='r2')]
        self.assertEqual(self.summarize(rows)['linked_fragments'],0)

    def test_reappearing_template_or_foreign_identity_rejected(self):
        for rows in ([obs('q2','a','REF'),obs('q1','b','ALT')],
                     [obs('q','a','REF',sample_id='other')],
                     [obs('q','c','REF')]):
            with self.assertRaises(ValueError):
                self.summarize(rows)

    def test_other_allele_is_ambiguous_not_reference(self):
        result=self.summarize([obs('q','a','OTHER'),obs('q','b','ALT')])
        self.assertEqual(result['fragment_status_counts']['AMBIGUOUS_FRAGMENT'],1)
        self.assertEqual(result['linked_fragments'],0)

    def test_empty_is_not_evaluated(self):
        self.assertEqual(self.summarize([])['status'],'PHASE_NOT_EVALUATED')
