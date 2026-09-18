import csv
import gzip
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from collect_mloy_intervals import collect


FIELDS = [
    "chrom", "start", "end", "parent_start", "parent_end", "compartment",
    "territory_bases", "depth_sum", "covered_bases", "bases_ge_10",
    "bases_ge_20", "bases_ge_30", "bases_ge_50", "bases_ge_80",
    "bases_ge_100", "bases_ge_200",
]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def interval(chrom, start, end, parent_start, parent_end, depth, compartment="OLD"):
    territory = end - start
    return {
        "chrom": chrom, "start": start, "end": end,
        "parent_start": parent_start, "parent_end": parent_end,
        "compartment": compartment, "territory_bases": territory,
        "depth_sum": depth * territory, "covered_bases": territory,
        "bases_ge_10": territory, "bases_ge_20": territory,
        "bases_ge_30": territory, "bases_ge_50": territory,
        "bases_ge_80": 0, "bases_ge_100": 0, "bases_ge_200": 0,
    }


class MloyCollectorTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.stages = self.root / "stages"
        self.mask = self.root / "mask.json"
        self.mask.write_text(json.dumps({
            "assembly": "GRCh38", "coordinate_system": "0-based-half-open",
            "par": {"chrY": [[0, 10]]},
            "xtr": {"chrY": [[30, 40], [50, 60]],
                    "status": "FROZEN_FOR_FUTURE_RECOMPUTATION_NOT_USED_IN_CURRENT_SEX_GATE"},
        }))
        self.fingerprint = "fp"

    def tearDown(self):
        self.tmp.cleanup()

    def make_stage(self, sample="S1", corrupt=False):
        directory = self.stages / "qc" / sample
        directory.mkdir(parents=True)
        rows = [
            interval("chr1", 0, 100, 0, 100, 100),
            interval("chrY", 10, 20, 10, 20, 30),
            interval("chrY", 25, 30, 25, 35, 30),
            interval("chrY", 30, 35, 25, 35, 30),
            interval("chrY", 40, 50, 40, 50, 20),
        ]
        path = directory / "coverage.intervals.tsv.gz"
        with gzip.open(path, "wt", newline="") as stream:
            writer = csv.DictWriter(stream, FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)
        digest = "0" * 64 if corrupt else sha(path)
        complete = directory / "COMPLETE.json"
        complete.write_text(json.dumps({
            "status": "COMPLETE", "fingerprint": self.fingerprint,
            "outputs": {"coverage.intervals.tsv.gz": digest},
        }))
        self.stages.mkdir(exist_ok=True)
        marker = self.stages / "QC_COMPLETE.json"
        marker.write_text(json.dumps({
            "status": "QC_COMPLETE", "fingerprint": self.fingerprint,
            "samples": {sample: {"receipt_sha256": sha(complete)}},
        }))

    def test_collects_eligible_parents_and_excludes_boundary_parent(self):
        self.make_stage()
        output = self.root / "out.json.gz"
        result = collect(self.stages, self.mask, output, 1, 1)
        self.assertEqual(result["eligible_y_intervals"], 2)
        with gzip.open(output, "rt") as stream:
            payload = json.load(stream)
        sample = payload["samples"][0]
        self.assertEqual(sample["autosomal"]["mean_depth"], 100)
        self.assertEqual(set(sample["eligible_y_intervals"]), {"chrY:10-20", "chrY:40-50"})
        excluded = {r["interval_id"]: r for r in payload["y_interval_catalog"]}
        self.assertEqual(excluded["chrY:25-35"]["status"], "EXCLUDED_XTR")
        self.assertEqual(excluded["chrY:25-35"]["conservative_overexcluded_bases"], 5)
        self.assertFalse(payload["bams_scanned"])
        self.assertFalse(payload["mloy_called"])

    def test_changed_interval_artifact_fails(self):
        self.make_stage(corrupt=True)
        with self.assertRaisesRegex(ValueError, "artifact changed"):
            collect(self.stages, self.mask, self.root / "out.json.gz", 1, 1)

    def test_membership_mismatch_fails(self):
        self.make_stage()
        with self.assertRaisesRegex(ValueError, "membership"):
            collect(self.stages, self.mask, self.root / "out.json.gz", 2, 1)


if __name__ == "__main__":
    unittest.main()
