import importlib.util
from pathlib import Path
import tempfile
import unittest

from cnv_execution import sha
from mloy_calibration import (
    run_coverage_calibration,
    thin_y_bam,
    y_keep_probability,
    y_template_retained,
)


class MloyCoverageCalibrationTests(unittest.TestCase):
    def fixture(self):
        collection_samples = []
        metrics = []
        depths = [30, 40, 50, 60, 70, 80, 90, 100, 110]
        for index, depth in enumerate(depths):
            sample_id = f"M{index}"
            intervals = {
                "chrY:10-110": {
                    "territory_bases": 100,
                    "depth_sum": depth * 50,
                    "depth_over_autosome": 0.5,
                },
                "chrY:200-300": {
                    "territory_bases": 100,
                    "depth_sum": depth * 50,
                    "depth_over_autosome": 0.5,
                },
            }
            collection_samples.append({
                "sample_id": sample_id,
                "autosomal": {"mean_depth": depth},
                "eligible_y_intervals": intervals,
            })
            metrics.append({
                "sample_id": sample_id,
                "analysis_id": str(index),
                "assay": "A",
                "study_group": "TP53_negative",
                "mloy_reference_eligible": True,
                "sex_chromosome_anomaly_review": False,
                "metric_status": "ELIGIBLE",
                "screen_tier_LOY": "WITHIN_REFERENCE_UNCALIBRATED",
                "normalization_median": -1.0,
                "raw_log2_y_auto_weighted": -1.0,
                "mLRR_Y_zero_aware": 0.0,
                "weak_signal_threshold_2sigma": -0.04,
                "extreme_signal_threshold_3sigma": -0.08,
            })
        collection = {
            "status": "MLOY_INTERVAL_COLLECTION_COMPLETE",
            "samples": collection_samples,
        }
        panel = {
            "assays": {
                "A": {
                    "intervals": [
                        {"interval_id": "chrY:10-110", "selected": True},
                        {"interval_id": "chrY:200-300", "selected": True},
                    ]
                }
            }
        }
        return collection, metrics, panel

    def test_cell_fraction_maps_to_expected_y_retention(self):
        self.assertEqual(y_keep_probability(0), 1)
        self.assertEqual(y_keep_probability(1), 0)
        self.assertAlmostEqual(y_keep_probability(0.4), 0.6)
        for invalid in (True, -0.1, 1.1, float("nan")):
            with self.assertRaises(ValueError):
                y_keep_probability(invalid)

    def test_template_decisions_are_repeatable_and_nested(self):
        mild = {i for i in range(1000) if y_template_retained("rg", str(i), 17, 0.2)}
        strong = {i for i in range(1000) if y_template_retained("rg", str(i), 17, 0.8)}
        self.assertTrue(strong < mild)
        self.assertEqual(mild, {i for i in range(1000) if y_template_retained("rg", str(i), 17, 0.2)})

    def test_coverage_calibration_stratifies_depth_and_never_calls_events(self):
        replicates, summary, plan = run_coverage_calibration(
            *self.fixture(),
            fractions=(0.0, 0.4),
            seeds=(17, 29),
            per_depth_stratum=1,
        )
        self.assertEqual(len(plan["donors"]), 3)
        self.assertEqual({row["depth_stratum"] for row in plan["donors"]}, {"LOW", "MID", "HIGH"})
        zero = [row for row in replicates if row["cell_fraction"] == 0]
        altered = [row for row in replicates if row["cell_fraction"] == 0.4]
        self.assertTrue(all(row["screen_tier_LOY"] == "WITHIN_REFERENCE_UNCALIBRATED" for row in zero))
        self.assertTrue(all(row["borderline_or_extreme"] for row in altered))
        self.assertTrue(all(row["event_call_LOY"] == "NOT_CALIBRATED" for row in replicates))
        self.assertFalse(summary["event_calls_emitted"])
        self.assertFalse(summary["negative_calls_emitted"])
        self.assertTrue(summary["bam_validation_required"])
        self.assertEqual(summary["provisional_borderline_lod_by_assay"]["A"], 0.4)

    def test_donor_selection_rejects_candidate_signal_as_background(self):
        collection, metrics, panel = self.fixture()
        for row in metrics:
            row["screen_tier_LOY"] = "BORDERLINE_CANDIDATE_UNCALIBRATED"
        with self.assertRaisesRegex(ValueError, "no eligible|insufficient"):
            run_coverage_calibration(
                collection,
                metrics,
                panel,
                fractions=(0.0, 0.4),
                seeds=(17,),
                per_depth_stratum=1,
            )


@unittest.skipUnless(importlib.util.find_spec("pysam"), "pysam required for native BAM tests")
class MloyBamThinningTests(unittest.TestCase):
    def setUp(self):
        import pysam
        self.pysam = pysam
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.bam = self.root / "source.bam"

    def make_bam(self):
        p = self.pysam
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 10000}, {"SN": "chrY", "LN": 10000}],
            "RG": [{"ID": "rg", "SM": "sample"}],
        }
        with p.AlignmentFile(str(self.bam), "wb", header=header) as stream:
            for reference_id in (0, 1):
                for index in range(200):
                    read = p.AlignedSegment()
                    read.query_name = f"r{reference_id}_{index}"
                    read.query_sequence = "A" * 50
                    read.query_qualities = p.qualitystring_to_array("I" * 50)
                    read.flag = 0
                    read.reference_id = reference_id
                    read.reference_start = index * 20
                    read.mapping_quality = 60
                    read.cigarstring = "50M"
                    read.set_tag("RG", "rg")
                    stream.write(read)
        p.index(str(self.bam))

    def test_bam_mode_preserves_all_non_y_and_thins_y(self):
        self.make_bam()
        receipt = thin_y_bam(
            {"path": str(self.bam), "sha256": sha(self.bam)},
            "sample",
            self.root / "result",
            seed=17,
            cell_fraction=0.5,
        )
        self.assertEqual(receipt["non_y_alignment_records"], 200)
        self.assertEqual(receipt["retained_non_y_alignment_records"], 200)
        self.assertEqual(receipt["non_y_records_removed"], 0)
        self.assertGreater(receipt["retained_y_alignment_records"], 70)
        self.assertLess(receipt["retained_y_alignment_records"], 130)
        self.assertFalse(receipt["event_calls_authorized"])
        self.assertTrue((self.root / "result" / "mloy.depth.bam.bai").exists())


if __name__ == "__main__":
    unittest.main()
