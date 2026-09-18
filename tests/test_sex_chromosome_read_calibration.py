import importlib.util
from pathlib import Path
import tempfile
import unittest

from cnv_execution import sha
from sex_chromosome_read_calibration import scan_sample


@unittest.skipUnless(importlib.util.find_spec("pysam"), "pysam required for native BAM tests")
class SexChromosomeReadCalibrationTests(unittest.TestCase):
    def setUp(self):
        import pysam
        self.pysam = pysam
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.bam = self.root / "sample.bam"
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [
                {"SN": "chr1", "LN": 1000},
                {"SN": "chrX", "LN": 1000},
                {"SN": "chrY", "LN": 1000},
            ],
            "RG": [{"ID": "rg", "SM": "sample"}],
        }
        with pysam.AlignmentFile(str(self.bam), "wb", header=header) as stream:
            for reference_id in (0, 1, 2):
                for index in range(40):
                    read = pysam.AlignedSegment()
                    read.query_name = f"r{reference_id}_{index}"
                    read.query_sequence = "A" * 50
                    read.query_qualities = pysam.qualitystring_to_array("I" * 50)
                    read.flag = 0
                    read.reference_id = reference_id
                    read.reference_start = 100 + index * 5
                    read.mapping_quality = 60
                    read.cigarstring = "50M"
                    read.set_tag("RG", "rg")
                    stream.write(read)
        pysam.index(str(self.bam))

    def binding(self):
        index = Path(str(self.bam) + ".bai")
        return {
            "sample_id": "S",
            "assay": "A",
            "bam_sample_name": "sample",
            "bam": {"path": str(self.bam), "sha256": sha(self.bam)},
            "bam_index": {"path": str(index), "sha256": sha(index)},
        }

    def frozen_depth(self, contig):
        total = 0
        with self.pysam.AlignmentFile(str(self.bam), "rb") as stream:
            for read in stream.fetch(contig, 100, 400):
                total += sum(max(0, min(right, 400) - max(left, 100)) for left, right in read.get_blocks())
        return total

    def test_mloy_target_only_scan_exactly_reconstructs_fraction_zero(self):
        depth_sum = self.frozen_depth("chrY")
        sample = {
            "sample_id": "S",
            "autosomal": {"mean_depth": 10.0},
            "eligible_y_intervals": {"chrY:100-400": {"territory_bases": 300, "depth_sum": depth_sum}},
        }
        metric = {
            "sample_id": "S", "analysis_id": "1", "assay": "A",
            "normalization_median": -1.0,
            "raw_log2_y_auto_weighted": -1.0, "mLRR_Y_zero_aware": 0.0,
            "weak_signal_threshold_2sigma": -0.04, "extreme_signal_threshold_3sigma": -0.08,
        }
        panel = {"intervals": [{"interval_id": "chrY:100-400", "selected": True}]}
        rows, receipt = scan_sample(
            "mloy", self.binding(),
            {"sample_id": "S", "assay": "A", "depth_stratum": "MID"},
            sample, metric, panel, fractions=(0.0, 0.5), seeds=(17,),
        )
        self.assertEqual(receipt["fraction_zero_filter_contract"], "EXACT_MATCH")
        self.assertFalse(receipt["pseudo_bams_written"])
        self.assertEqual(rows[0]["negative_call"], "NOT_EMITTED")
        self.assertLess(rows[1]["read_level_raw_log2_ratio_median"], rows[0]["read_level_raw_log2_ratio_median"])

    def test_mlox_mismatch_fails_closed(self):
        sample = {
            "sample_id": "S", "autosomal": {"mean_depth": 10.0},
            "eligible_x_intervals": {"chrX:100-400": {"territory_bases": 300, "depth_sum": self.frozen_depth("chrX") + 1}},
        }
        metric = {
            "sample_id": "S", "assay": "A", "female_normalization_median": 0.0,
            "male_dosage_median": -1.0, "weak_signal_threshold_2sigma": -0.04,
            "extreme_signal_threshold_3sigma": -0.08,
        }
        panel = {"intervals": [{"interval_id": "chrX:100-400", "selected": True}]}
        with self.assertRaisesRegex(ValueError, "filter contract mismatch"):
            scan_sample(
                "mlox", self.binding(),
                {"sample_id": "S", "assay": "A", "depth_stratum": "MID"},
                sample, metric, panel, fractions=(0.0, 0.5), seeds=(17,),
            )


if __name__ == "__main__":
    unittest.main()
