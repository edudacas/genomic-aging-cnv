import gzip
import json
from pathlib import Path
import tempfile
import unittest
from collect_snv_qc import collect, sha


class QCCollectionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.sample=self.root/'qc/S'
        self.sample.mkdir(parents=True)
        outputs={}
        for name in ('qc.json','coverage.summary.json','coverage.compartments.tsv'):
            p=self.sample/name
            p.write_text('{}\n')
            outputs[name]=sha(p)
        (self.sample/'COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',fingerprint='f',outputs=outputs)))
        (self.root/'QC_COMPLETE.json').write_text(json.dumps(dict(status='QC_COMPLETE',fingerprint='f',
            samples={'S':dict(receipt_sha256=sha(self.sample/'COMPLETE.json'))})))

    def test_valid_collection(self):
        out=self.root/'out.gz'
        collect(self.root,out,1)
        with gzip.open(out,'rt') as f:
            data=json.load(f)
        self.assertEqual(data['selected_samples'],1)
        self.assertFalse(data['sex_inferred'])
        self.assertFalse(data['bams_scanned'])

    def test_artifact_tamper(self):
        (self.sample/'qc.json').write_text('changed')
        with self.assertRaisesRegex(ValueError,'artifact changed'):
            collect(self.root,self.root/'out.gz',1)

    def test_receipt_tamper(self):
        (self.sample/'COMPLETE.json').write_text('{}')
        with self.assertRaisesRegex(ValueError,'receipt mismatch'):
            collect(self.root,self.root/'out.gz',1)

    def test_expected_count(self):
        with self.assertRaisesRegex(ValueError,'membership'):
            collect(self.root,self.root/'out.gz',2)
