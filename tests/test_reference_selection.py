import copy
import unittest
from cnv_reference_selection import select_reference
from test_contract import fixture


def manifest():
    value = fixture()
    second = copy.deepcopy(value['samples'][1])
    second.update(sample_id='pon2',person_id='pon2')
    second['bam']['path']='/example/pon2.bam'
    value['samples'].append(second)
    for sample in value['samples']:
        sample['reference_sex']='Female'
    return value


def select(value, signals, loci=('a','b','c'), **kwargs):
    options=dict(minimum_common_loci=2,minimum_common_fraction=0.5)
    options.update(kwargs)
    return select_reference(value,'study',signals,loci,**options)


class ReferenceSelectionTests(unittest.TestCase):
    def test_all_candidates_use_same_marker_set(self):
        result=select(manifest(),{'pon':{'a':0.1,'b':0.1,'c':4},
                                  'pon2':{'a':0.2,'b':0.2,'c':None}})
        self.assertEqual(result['selected_reference'],'pon')
        self.assertEqual(result['common_loci'],2)
        self.assertAlmostEqual(result['candidate_scores']['pon'],0.02)
        self.assertFalse(result['biological_validation_performed'])

    def test_missing_common_evidence_is_not_zero(self):
        with self.assertRaisesRegex(ValueError,'insufficient common'):
            select(manifest(),{'pon':{'a':0.1},'pon2':{'b':0.1}})

    def test_unknown_sex_never_defaults(self):
        value=manifest();value['samples'][0]['reference_sex']='Unknown'
        with self.assertRaisesRegex(ValueError,'unresolved'):
            select(value,{'pon':{},'pon2':{}})

    def test_cannot_omit_normal_or_insert_study(self):
        for signals in ({'pon':{}},{'pon':{},'pon2':{},'study':{}}):
            with self.subTest(signals=signals),self.assertRaisesRegex(ValueError,'exactly'):
                select(manifest(),signals)

    def test_stable_tie_and_same_sex_only(self):
        signals={'pon':{'a':0,'b':0},'pon2':{'a':0,'b':0}}
        self.assertEqual(select(manifest(),signals)['selected_reference'],'pon')
        value=manifest();value['samples'][-1]['reference_sex']='Male'
        self.assertEqual(select(value,{'pon':signals['pon']})['selected_reference'],'pon')

    def test_rejects_nonfinite_and_undeclared_loci(self):
        for signal in ({'a':float('nan')},{'undeclared':0}):
            with self.subTest(signal=signal),self.assertRaises(ValueError):
                select(manifest(),{'pon':signal,'pon2':{}})
