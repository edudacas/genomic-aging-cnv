import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from cnv_compare_export import export_comparison
from cnv_execution import sha


@unittest.skipUnless(importlib.util.find_spec('pyarrow'),'PyArrow required')
class CompareExportTests(unittest.TestCase):
    def setUp(self):
        temp=tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root=Path(temp.name)
        self.table=self.root/'comparison.jsonl'
        self.row=dict(sample_id='s',chromosome='chr1',start_0based=0,end_0based_exclusive=10,
            comparison_status='CNVKIT_ONLY',cnvkit_source={'mean_log2_ratio':0},facets_source=None,
            event_class='UNCLASSIFIED',evaluable_bp=None,somatic_origin_classified=False,
            copy_neutral_loh_confirmed=False,sensitivity_calibrated=False)
        self.started=self.root/'STARTED.json'
        self.started.write_text(json.dumps(dict(sample_provenance={'sample_id':'s'},inputs=[])))
        self.receipt=self.root/'COMPLETE.json'
        self.seal()

    def seal(self,count=1):
        self.table.write_text(json.dumps(self.row)+'\n')
        self.receipt.write_text(json.dumps(dict(status='CALLER_COMPARISON_COMPLETE',
            sample_provenance={'sample_id':'s'},started_sha256=sha(self.started),
            outputs={str(self.table):dict(sha256=sha(self.table),bytes=self.table.stat().st_size)},
            comparison_status_counts={'CNVKIT_ONLY':count},comparison_intervals=count)))
        self.artifact=dict(path=str(self.receipt),sha256=sha(self.receipt))

    def invoke(self,sample='s'):
        return export_comparison(self.artifact,sample,self.root/'export',batch_size=1)

    def test_payload_null_and_observed_zero_survive(self):
        result=self.invoke()
        self.assertTrue(result['comparison_payload_preserved'])
        row=json.loads((self.root/'export/candidates.jsonl').read_text())
        self.assertEqual(json.loads(row['payload_json']),self.row)
        self.assertIsNone(row['evaluable_bp'])
        self.assertFalse(result['event_count_estimated'])

    def test_wrong_sample_rejected_before_export(self):
        with self.assertRaisesRegex(ValueError,'identity mismatch'):
            self.invoke('other')
        self.assertFalse((self.root/'export').exists())

    def test_changed_comparison_rejected(self):
        self.table.write_text(self.table.read_text()+'\n')
        with self.assertRaisesRegex(ValueError,'output changed'):
            self.invoke()

    def test_count_mismatch_cannot_complete(self):
        self.seal(count=2)
        with self.assertRaisesRegex(ValueError,'completion counts'):
            self.invoke()
        self.assertTrue((self.root/'export/FAILURE.json').exists())
        self.assertFalse((self.root/'export/COMPLETE.json').exists())

    def test_biological_claim_cannot_enter_candidate_export(self):
        self.row['somatic_origin_classified']=True
        self.seal()
        with self.assertRaisesRegex(ValueError,'biological claim'):
            self.invoke()
        self.assertFalse((self.root/'export/COMPARISON_EXPORT.json').exists())
