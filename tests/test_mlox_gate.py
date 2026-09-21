import unittest

from mlox_gate import build_gate


def interval(depth, territory=100):
    return {
        "territory_bases": territory,
        "depth_over_autosome": depth,
        "bases_ge_20": territory,
    }


class MloxGateTests(unittest.TestCase):
    def fixture(self):
        samples = []
        sex = []
        metadata = {}
        female_depths = [0.98, 1.00, 1.02, 1.01, 0.99, 1.03, 0.97]
        male_depths = [0.48, 0.50, 0.52, 0.49, 0.51, 0.47]
        index = 0
        for label, depths in (("Female", female_depths), ("Male", male_depths)):
            for depth in depths:
                sample_id = f"S{index}"
                samples.append({
                    "sample_id": sample_id,
                    "eligible_x_intervals": {
                        "chrX:10-110": interval(depth),
                        "chrX:200-300": interval(depth * 1.01),
                    },
                })
                sex.append({
                    "sample_id": sample_id, "analysis_id": str(index),
                    "sex_inferred": label, "sex_gate_status": "PASS",
                    "sex_reported": label, "proposed_capture_family": "A",
                    "batch_id": "B", "batch_name": "batch",
                })
                metadata[str(index)] = {
                    "age_at_collection": 10 + index,
                    "study_group": "TP53_negative",
                }
                index += 1
        samples.append({
            "sample_id": "CASE",
            "eligible_x_intervals": {
                "chrX:10-110": interval(0.80),
                "chrX:200-300": interval(0.808),
            },
        })
        sex.append({
            "sample_id": "CASE", "analysis_id": "CASE",
            "sex_inferred": "Female", "sex_gate_status": "PASS",
            "sex_reported": "Female", "proposed_capture_family": "A",
            "batch_id": "B", "batch_name": "batch",
        })
        metadata["CASE"] = {"age_at_collection": 60, "study_group": "LFS_R337H"}
        return ({
            "status": "MLOX_INTERVAL_COLLECTION_COMPLETE",
            "selected_samples": len(samples),
            "mask_sha256": "mask",
            "samples": samples,
        }, sex, metadata)

    def test_retains_diffuse_x_loss_without_event_call(self):
        rows, summary, panel = build_gate(*self.fixture())
        case = next(row for row in rows if row["sample_id"] == "CASE")
        self.assertEqual(panel["assays"]["A"]["selected_interval_n"], 2)
        self.assertGreater(case["f_LOX_interpretive"], 0.30)
        self.assertEqual(case["screen_tier_LOX"], "EXTREME_CANDIDATE_UNCALIBRATED")
        self.assertEqual(case["event_call_LOX"], "NOT_CALIBRATED")
        self.assertFalse(summary["event_calls_emitted"])
        self.assertFalse(summary["negative_calls_emitted"])

    def test_sex_anomaly_is_retained(self):
        collection, sex, metadata = self.fixture()
        collection["samples"].append({
            "sample_id": "AMB", "eligible_x_intervals": {
                "chrX:10-110": interval(0.65), "chrX:200-300": interval(0.66),
            },
        })
        collection["selected_samples"] += 1
        sex.append({
            "sample_id": "AMB", "analysis_id": "AMB", "sex_inferred": None,
            "sex_gate_status": "HOLD", "sex_reported": "Unknown",
            "proposed_capture_family": "A", "batch_id": "B", "batch_name": "batch",
        })
        metadata["AMB"] = {"age_at_collection": 50, "study_group": "TP53_negative"}
        rows, summary, _ = build_gate(collection, sex, metadata)
        case = next(row for row in rows if row["sample_id"] == "AMB")
        self.assertEqual(case["screen_tier_LOX"], "SEX_CHROMOSOME_ANOMALY_REVIEW")
        self.assertEqual(summary["sex_chromosome_anomaly_review_samples"], 1)

    def test_collection_membership_mismatch_fails(self):
        collection, sex, metadata = self.fixture()
        collection["samples"].pop()
        with self.assertRaisesRegex(ValueError, "sample mismatch"):
            build_gate(collection, sex, metadata)

    def test_missing_age_is_retained_and_excluded_only_from_young_null(self):
        collection, sex, metadata = self.fixture()
        metadata["CASE"]["age_at_collection"] = None
        rows, summary, _ = build_gate(collection, sex, metadata)
        case = next(row for row in rows if row["sample_id"] == "CASE")
        self.assertIsNone(case["age_at_collection"])
        self.assertEqual(case["age_status"], "MISSING_RETAINED")
        self.assertIsNotNone(case["mLRR_X"])
        self.assertFalse(summary["negative_calls_emitted"])

    def test_held_out_control_is_scored_but_excluded_from_training(self):
        rows, summary, panel = build_gate(
            *self.fixture(), training_excluded_sample_ids={"S0"}
        )
        held_out = next(row for row in rows if row["sample_id"] == "S0")
        self.assertTrue(held_out["held_out_specificity"])
        self.assertEqual(held_out["female_normalization_control_n"], 6)
        self.assertEqual(panel["assays"]["A"]["reference_female_n"], 7)
        self.assertEqual(summary["training_excluded_sample_n"], 1)


if __name__ == "__main__":
    unittest.main()
