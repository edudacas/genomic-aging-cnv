import csv
import gzip
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from cnv_coordinates import bed_intervals
from cnv_execution import sha
from cnv_prepare_pair import prepare
import test_pair_matrix


class PreparePairTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        fixture = test_pair_matrix.PairMatrixTests()
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        self.rows = fixture.rows()
        manifest = fixture.manifest
        def artifact(name, text):
            path = self.root / name
            path.write_text(text)
            return dict(path=str(path), sha256=sha(path))
        manifest['strata']['S']['reference_fai'] = artifact('ref.fai', 'chr1\t10\t0\t10\t11\n')
        catalog = artifact('catalog.tsv', 'chromosome\tposition\tref\talt\tkind\nchr1\t1\tA\tC\tSNP\n')
        mask = artifact('mask.bed', '')
        self.manifest = self.root / 'manifest.json'
        self.manifest.write_text(json.dumps(manifest))
        registry = {}
        for sample in ('study', 'pon'):
            directory = self.root / sample
            directory.mkdir()
            for name in ('STARTED.json', 'COMPLETE.json', 'pileup.csv'):
                (directory / name).write_text('{}')
            registry[sample] = dict(task_directory=str(directory), output=str(directory / 'pileup.csv'))
        self.config = self.root / 'config.json'
        self.config.write_text(json.dumps(dict(sample_id='study', reference_id='pon',
            catalog=catalog, mask_bed=mask, registry=registry, executable='/mock/native',
            pileup_policy={}, rehash_inputs=True)))
        self.output = self.root / 'prepared'

    def run_preparation(self, rows=None):
        # This tests stage serialization, not native receipt or scientific validity.
        with patch('cnv_prepare_pair.verified_pair_matrix', return_value=iter(self.rows) if rows is None else rows):
            return prepare(self.manifest, self.config, self.output)

    def test_preserves_missing_and_observed_zero(self):
        receipt = self.run_preparation()
        self.assertEqual(receipt['catalog_loci'], 3)
        self.assertEqual(receipt['native_rows'], 1)
        self.assertEqual(receipt['incomplete_loci'], 2)
        self.assertFalse(receipt['facets_ch_executed'])
        self.assertFalse(receipt['native_joint_fit_validated'])
        with gzip.open(self.output / 'pair.evidence.jsonl.gz', 'rt') as stream:
            rows = [json.loads(line) for line in stream]
        self.assertEqual(rows[0]['TUM.DP'], 0)
        self.assertIsNone(rows[1]['NOR.DP'])
        with (self.output / 'pair.rcmat.tsv').open() as stream:
            native = list(csv.DictReader(stream, delimiter='\t'))
        self.assertEqual(native[0]['masked'], 'TRUE')
        self.assertEqual(native[0]['TUM.DP'], '0')
        for name, metadata in receipt['outputs'].items():
            self.assertEqual(metadata['sha256'], sha(self.output / name))
        with self.assertRaises(FileExistsError):
            self.run_preparation()

    def test_midstream_failure_cannot_complete(self):
        def broken():
            yield self.rows[0]
            raise ValueError('stream validation failed')
        with self.assertRaisesRegex(ValueError, 'stream validation'):
            self.run_preparation(broken())
        self.assertTrue((self.output / 'FAILURE.json').exists())
        self.assertFalse((self.output / 'COMPLETE.json').exists())

    def test_input_mutation_cannot_complete(self):
        def changing():
            yield from self.rows
            self.manifest.write_text('{}')
        with self.assertRaisesRegex(ValueError, 'changed during'):
            self.run_preparation(changing())
        self.assertFalse((self.output / 'COMPLETE.json').exists())

    def test_empty_mask_requires_explicit_allowance(self):
        with self.assertRaises(ValueError):
            bed_intervals(self.root / 'mask.bed', {'chr1': 10})
        self.assertEqual(bed_intervals(self.root / 'mask.bed', {'chr1': 10}, allow_empty=True), [])
