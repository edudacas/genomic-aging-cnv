import unittest

from mloy_gate import build_gate


def interval(depth, territory=100, covered_fraction=1.0):
    covered = int(territory * covered_fraction)
    return {
        "territory_bases": territory,
        "depth_over_autosome": depth,
        "bases_ge_20": covered,
    }


class MloyGateTests(unittest.TestCase):
    def fixture(self):
        samples = []
        sex = []
        metadata = {}
        control_depths = [0.46, 0.48, 0.50, 0.52, 0.54, 0.49, 0.51, 0.53]
        for index in range(9):
            sample_id = f"S{index}"
            depth = control_depths[index] if index < 8 else 0.25
            samples.append({
                "sample_id": sample_id,
                "eligible_y_intervals": {
                    "chrY:10-110": interval(depth),
                    "chrY:200-300": interval(depth, covered_fraction=0.0),
                },
            })
            sex.append({
                "sample_id": sample_id, "analysis_id": str(index),
                "sex_inferred": "Male", "sex_gate_status": "PASS",
                "sex_reported": "Male", "proposed_capture_family": "A",
                "batch_id": "B", "batch_name": "batch",
            })
            metadata[str(index)] = {
                "age_at_collection": 10 + index,
                "study_group": "TP53_negative" if index < 8 else "LFS_R337H",
            }
        collection = {
            "status": "MLOY_INTERVAL_COLLECTION_COMPLETE",
            "selected_samples": len(samples),
            "mask_sha256": "mask",
            "samples": samples,
        }
        return collection, sex, metadata

    def test_freezes_panel_normalizes_and_blocks_event_calls(self):
        rows, summary, panel = build_gate(*self.fixture())
        self.assertEqual(panel["assays"]["A"]["selected_interval_n"], 1)
        case = next(row for row in rows if row["study_group"] == "LFS_R337H")
        self.assertLess(case["mLRR_Y"], -0.9)
        self.assertGreater(case["f_LOY_interpretive"], 0.45)
        self.assertEqual(case["screen_tier_LOY"], "EXTREME_CANDIDATE_UNCALIBRATED")
        self.assertGreater(case["fraction_targets_below_control"], 0.0)
        self.assertEqual(case["event_call_LOY"], "NOT_CALIBRATED")
        self.assertFalse(summary["event_calls_emitted"])
        self.assertFalse(summary["negative_calls_emitted"])

    def test_missing_collection_sample_fails(self):
        collection, sex, metadata = self.fixture()
        collection["samples"].pop()
        collection["selected_samples"] -= 1
        with self.assertRaisesRegex(ValueError, "sample mismatch"):
            build_gate(collection, sex, metadata)

    def test_absent_y_is_manual_review_not_event(self):
        collection, sex, metadata = self.fixture()
        collection["samples"][-1]["eligible_y_intervals"]["chrY:10-110"] = interval(0.0)
        rows, _, _ = build_gate(collection, sex, metadata)
        case = next(row for row in rows if row["study_group"] == "LFS_R337H")
        self.assertEqual(case["metric_status"], "Y_SIGNAL_ABSENT_MANUAL_REVIEW")
        self.assertEqual(case["screen_flag_LOY"], "Y_SIGNAL_ABSENT_MANUAL_REVIEW")
        self.assertEqual(case["event_call_LOY"], "NOT_CALIBRATED")

    def test_borderline_signal_is_retained_not_negative(self):
        collection, sex, metadata = self.fixture()
        collection["samples"][-1]["eligible_y_intervals"]["chrY:10-110"] = interval(0.47)
        rows, _, _ = build_gate(collection, sex, metadata)
        case = next(row for row in rows if row["study_group"] == "LFS_R337H")
        self.assertEqual(case["screen_tier_LOY"], "BORDERLINE_CANDIDATE_UNCALIBRATED")
        self.assertNotIn("NEGATIVE", case["screen_tier_LOY"])

    def test_sex_gate_hold_is_retained_for_anomaly_review(self):
        collection, sex, metadata = self.fixture()
        collection["samples"].append({
            "sample_id": "SX",
            "eligible_y_intervals": {
                "chrY:10-110": interval(0.40),
                "chrY:200-300": interval(0.40, covered_fraction=0.0),
            },
        })
        collection["selected_samples"] += 1
        sex.append({
            "sample_id": "SX", "analysis_id": "X",
            "sex_inferred": None, "sex_gate_status": "HOLD",
            "sex_reported": "Unknown", "proposed_capture_family": "A",
            "batch_id": "B", "batch_name": "batch",
        })
        metadata["X"] = {"age_at_collection": 50, "study_group": "TP53_negative"}
        rows, summary, _ = build_gate(collection, sex, metadata)
        case = next(row for row in rows if row["sample_id"] == "SX")
        self.assertEqual(case["screen_tier_LOY"], "SEX_CHROMOSOME_ANOMALY_REVIEW")
        self.assertFalse(case["mloy_reference_eligible"])
        self.assertEqual(summary["sex_chromosome_anomaly_review_samples"], 1)

    def test_held_out_control_is_scored_but_excluded_from_training(self):
        rows, summary, panel = build_gate(
            *self.fixture(), training_excluded_sample_ids={"S0"}
        )
        held_out = next(row for row in rows if row["sample_id"] == "S0")
        self.assertTrue(held_out["held_out_specificity"])
        self.assertEqual(held_out["normalization_control_n"], 7)
        self.assertEqual(panel["assays"]["A"]["reference_male_n"], 8)
        self.assertEqual(summary["training_excluded_sample_n"], 1)


if __name__ == "__main__":
    unittest.main()
