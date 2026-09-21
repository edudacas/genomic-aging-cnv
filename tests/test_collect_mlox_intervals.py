import csv
import gzip
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from collect_mlox_intervals import collect


FIELDS = [
    "chrom", "start", "end", "parent_start", "parent_end", "compartment",
    "territory_bases", "depth_sum", "covered_bases", "bases_ge_10",
    "bases_ge_20", "bases_ge_30", "bases_ge_50", "bases_ge_80",
    "bases_ge_100", "bases_ge_200",
]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def interval(chrom, start, end, parent_start, parent_end, depth):
    territory = end - start
    return {
        "chrom": chrom, "start": start, "end": end,
        "parent_start": parent_start, "parent_end": parent_end,
        "compartment": "OLD", "territory_bases": territory,
        "depth_sum": depth * territory, "covered_bases": territory,
        "bases_ge_10": territory, "bases_ge_20": territory,
        "bases_ge_30": territory, "bases_ge_50": territory,
        "bases_ge_80": 0, "bases_ge_100": 0, "bases_ge_200": 0,
    }


class MloxCollectorTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.stages = self.root / "stages"
        self.mask = self.root / "mask.json"
        self.mask.write_text(json.dumps({
            "assembly": "GRCh38", "coordinate_system": "0-based-half-open",
            "par": {"chrX": [[0, 10]]},
            "xtr": {
                "chrX": [[30, 40]],
                "status": "FROZEN_FOR_FUTURE_RECOMPUTATION_NOT_USED_IN_CURRENT_SEX_GATE",
            },
        }))
        self.fingerprint = "fp"

    def tearDown(self):
        self.tmp.cleanup()

    def make_stage(self, corrupt=False):
        directory = self.stages / "qc" / "S1"
        directory.mkdir(parents=True)
        rows = [
            interval("chr1", 0, 100, 0, 100, 100),
            interval("chrX", 10, 20, 10, 20, 100),
            interval("chrX", 25, 30, 25, 35, 100),
            interval("chrX", 30, 35, 25, 35, 100),
            interval("chrX", 40, 50, 40, 50, 80),
        ]
        path = directory / "coverage.intervals.tsv.gz"
        with gzip.open(path, "wt", newline="") as stream:
            writer = csv.DictWriter(stream, FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        digest = "0" * 64 if corrupt else sha(path)
        complete = directory / "COMPLETE.json"
        complete.write_text(json.dumps({
            "status": "COMPLETE", "fingerprint": self.fingerprint,
            "outputs": {"coverage.intervals.tsv.gz": digest},
        }))
        marker = self.stages / "QC_COMPLETE.json"
        marker.write_text(json.dumps({
            "status": "QC_COMPLETE", "fingerprint": self.fingerprint,
            "samples": {"S1": {"receipt_sha256": sha(complete)}},
        }))

    def test_collects_x_nonpar_nonxtr_without_bam_scan(self):
        self.make_stage()
        output = self.root / "x.json.gz"
        result = collect(self.stages, self.mask, output, 1, 1)
        self.assertEqual(result["eligible_x_intervals"], 2)
        with gzip.open(output, "rt") as stream:
            payload = json.load(stream)
        self.assertFalse(payload["bams_scanned"])
        self.assertFalse(payload["mlox_called"])
        sample = payload["samples"][0]
        self.assertEqual(set(sample["eligible_x_intervals"]), {"chrX:10-20", "chrX:40-50"})
        excluded = {r["interval_id"]: r for r in payload["x_interval_catalog"]}
        self.assertEqual(excluded["chrX:25-35"]["status"], "EXCLUDED_XTR")
        self.assertEqual(excluded["chrX:25-35"]["conservative_overexcluded_bases"], 5)

    def test_changed_artifact_fails_closed(self):
        self.make_stage(corrupt=True)
        with self.assertRaisesRegex(ValueError, "artifact changed"):
            collect(self.stages, self.mask, self.root / "x.json.gz", 1, 1)


if __name__ == "__main__":
    unittest.main()
