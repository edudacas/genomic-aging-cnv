import json
import unittest
from cnv_export_contract import candidate_record,validate_candidate,validated_candidates


class ExportContractTests(unittest.TestCase):
    def record(self,payload=None):
        return candidate_record('s','r','FRAGMENT_PHASE',payload or {'missing':None,'zero':0},'a'*64)

    def test_lossless_missing_zero_and_nested_data(self):
        payload={'missing':None,'zero':0,'nested':{'alleles':['ALT',None],'supported':False}}
        record=self.record(payload)
        self.assertEqual(json.loads(record['payload_json']),payload)
        self.assertIsNone(record['start_0based'])
        self.assertIsNone(record['evaluable_bp'])

    def test_candidate_cannot_be_promoted_to_call(self):
        for key,value in [('event_class','LOSS'),('evaluable_bp',0),('somatic_origin_classified',True),('sensitivity_calibrated',True)]:
            with self.subTest(key=key),self.assertRaises(ValueError):
                validate_candidate(dict(self.record(),**{key:value}))

    def test_no_partial_or_boolean_coordinates(self):
        for interval in (dict(chromosome='chr1'),dict(chromosome='chr1',start_0based=True,end_0based_exclusive=2)):
            with self.assertRaises(ValueError):
                candidate_record('s','r','SEGMENT',{},'a'*64,**interval)

    def test_nonfinite_evidence_rejected(self):
        with self.assertRaises(ValueError):
            self.record({'log2':float('nan')})

    def test_duplicates_not_silently_dropped(self):
        with self.assertRaisesRegex(ValueError,'duplicate'):
            list(validated_candidates([self.record(),self.record()]))

    def test_valid_interval_preserved(self):
        record=candidate_record('s','r','SEGMENT',{},'a'*64,chromosome='chr1',
                                start_0based=9,end_0based_exclusive=20)
        self.assertEqual(record['start_0based'],9)
        validate_candidate(record)
