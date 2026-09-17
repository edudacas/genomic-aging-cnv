import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from cnv_execution import sha
from cnv_export_contract import candidate_record
from cnv_parquet import export_candidates


@unittest.skipUnless(importlib.util.find_spec('pyarrow'),'PyArrow unavailable: native export NOT validated')
class ParquetTests(unittest.TestCase):
    def setUp(self):
        temp=tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root=Path(temp.name)
        self.source=self.root/'source.json'
        self.source.write_text('{"missing": null, "zero": 0}')
        self.digest=sha(self.source)
        self.records=[candidate_record('s',str(i),'SEGMENT',json.loads(self.source.read_text()),
            self.digest,chromosome='chr1',start_0based=i,end_0based_exclusive=i+1) for i in range(3)]

    def export(self,records):
        return export_candidates(records,self.root/'out',source_artifacts={self.digest:str(self.source)},batch_size=2)

    def test_roundtrip_nullable_fields_and_multiple_batches(self):
        receipt=self.export(iter(self.records))
        self.assertEqual(receipt['rows'],3)
        self.assertTrue(receipt['roundtrip_fields_match'])
        self.assertFalse(receipt['scientific_validation'])
        with self.assertRaises(FileExistsError):
            self.export(self.records)

    def test_empty_typed_export(self):
        self.assertEqual(self.export([])['rows'],0)

    def test_duplicate_leaves_failure_not_completion(self):
        with self.assertRaisesRegex(ValueError,'duplicate'):
            self.export([self.records[0],self.records[0]])
        self.assertTrue((self.root/'out/FAILURE.json').exists())
        self.assertFalse((self.root/'out/COMPLETE.json').exists())
