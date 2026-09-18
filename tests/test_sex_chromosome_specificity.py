import unittest

from sex_chromosome_specificity import cross_validate_specificity


def interval(depth, territory=100):
    return {
        "territory_bases": territory,
        "depth_over_autosome": depth,
        "bases_ge_20": territory,
    }


class SexChromosomeSpecificityTests(unittest.TestCase):
    def fixture(self, axis):
        samples, sex, metadata = [], [], {}
        target_sex = "Male" if axis == "mloy" else "Female"
        target_contig = "y" if axis == "mloy" else "x"
        base = 0.5 if axis == "mloy" else 1.0
        count = 18
        for index in range(count):
            sample_id = f"T{index:02d}"
            depth = 40 + index
            ratio = base + ((index % 7) - 3) * 0.005
            samples.append({
                "sample_id": sample_id,
                "autosomal": {"mean_depth": depth},
                f"eligible_{target_contig}_intervals": {
                    f"chr{target_contig.upper()}:10-110": interval(ratio),
                    f"chr{target_contig.upper()}:200-300": interval(ratio * 1.002),
                },
            })
            sex.append({
                "sample_id": sample_id,
                "analysis_id": str(index),
                "sex_inferred": target_sex,
                "sex_reported": target_sex,
                "sex_gate_status": "PASS",
                "proposed_capture_family": "A",
                "batch_id": f"B{index % 2}",
                "batch_name": "batch",
            })
            metadata[str(index)] = {
                "age_at_collection": 15 + index,
                "study_group": "TP53_negative",
            }
        if axis == "mlox":
            for offset in range(8):
                index = count + offset
                sample_id = f"M{offset:02d}"
                depth = 45 + offset
                ratio = 0.5 + ((offset % 5) - 2) * 0.005
                samples.append({
                    "sample_id": sample_id,
                    "autosomal": {"mean_depth": depth},
                    "eligible_x_intervals": {
                        "chrX:10-110": interval(ratio),
                        "chrX:200-300": interval(ratio * 1.002),
                    },
                })
                sex.append({
                    "sample_id": sample_id,
                    "analysis_id": str(index),
                    "sex_inferred": "Male",
                    "sex_reported": "Male",
                    "sex_gate_status": "PASS",
                    "proposed_capture_family": "A",
                    "batch_id": f"B{offset % 2}",
                    "batch_name": "batch",
                })
                metadata[str(index)] = {
                    "age_at_collection": 15 + index,
                    "study_group": "TP53_negative",
                }
        collection = {
            "status": "MLOY_INTERVAL_COLLECTION_COMPLETE" if axis == "mloy" else "MLOX_INTERVAL_COLLECTION_COMPLETE",
            "selected_samples": len(samples),
            "mask_sha256": "mask",
            "samples": samples,
        }
        return collection, sex, metadata

    def test_mloy_controls_are_scored_once_without_fold_leakage(self):
        rows, summary, plan = cross_validate_specificity(
            "mloy", *self.fixture("mloy"), folds=3
        )
        self.assertEqual(len(rows), 18)
        self.assertEqual(len({row["sample_id"] for row in rows}), 18)
        self.assertTrue(all(row["negative_call"] == "NOT_EMITTED" for row in rows))
        self.assertTrue(summary["callable_panel_refit_without_held_out_fold"])
        self.assertFalse(summary["event_calls_emitted"])
        self.assertEqual({row["fold"] for row in plan["controls"]}, {0, 1, 2})

    def test_mlox_controls_are_scored_once_without_truth_leakage(self):
        rows, summary, _ = cross_validate_specificity(
            "mlox", *self.fixture("mlox"), folds=3
        )
        self.assertEqual(len(rows), 18)
        self.assertFalse(summary["orthogonal_truth_used_for_fitting"])
        assay = summary["held_out_reference_behavior_by_assay"]["A"]
        self.assertIn("retained_signal_rate_wilson_95", assay)
        self.assertFalse(assay["false_positive_rate_interpretable"])


if __name__ == "__main__":
    unittest.main()
