import unittest

from cnv_sex_gate import calibrate_and_apply


def row(sample, assay, reported, x, y=0.0):
    return {
        "sample_id": sample,
        "person_id": "P:" + sample,
        "proposed_capture_family": assay,
        "sex_reported": reported,
        "sex_inferred": None,
        "sex_inference_status": "NOT_ASSESSED_GATE_REQUIRED",
        "x_auto_ratio": x,
        "y_auto_ratio": y,
        "hold_reasons": ["PROVISIONAL_SEX_COMPARTMENT_MASK"],
    }


class SexGateTests(unittest.TestCase):
    def cohort(self):
        return [
            row("m1", "A", "Male", 0.48, 0.3),
            row("m2", "A", "Male", 0.52, 0.0),
            row("f1", "A", "Female", 0.96, 0.0),
            row("f2", "A", "Female", 1.04, 0.01),
            row("u1", "A", "Unknown", 0.50, 0.0),
            row("u2", "A", "Unknown", 1.00, 0.0),
            row("u3", "A", "Unknown", 0.70, 0.30),
        ]

    def test_known_concordance_and_unknown_assignment(self):
        rows, summary, thresholds = calibrate_and_apply(self.cohort(), "abc")
        by_id = {r["sample_id"]: r for r in rows}
        self.assertEqual(by_id["u1"]["sex_inferred"], "Male")
        self.assertEqual(by_id["u2"]["sex_inferred"], "Female")
        self.assertIsNone(by_id["u3"]["sex_inferred"])
        self.assertEqual(by_id["u3"]["sex_gate_status"], "HOLD")
        self.assertFalse(by_id["u3"]["sex_inference_uses_y"])
        self.assertEqual(summary["unknown_inference_counts"], {"Ambiguous": 1, "Female": 1, "Male": 1})
        self.assertEqual(thresholds["assays"]["A"]["leave_one_out_accuracy"], 1.0)

    def test_reported_sex_is_never_overwritten(self):
        rows, _, _ = calibrate_and_apply(self.cohort(), "abc")
        self.assertEqual([r["sex_reported"] for r in rows], [r["sex_reported"] for r in self.cohort()])

    def test_input_is_not_mutated(self):
        cohort = self.cohort()
        calibrate_and_apply(cohort, "abc")
        self.assertEqual(cohort[0]["sex_inference_status"], "NOT_ASSESSED_GATE_REQUIRED")

    def test_overlap_fails_closed(self):
        cohort = self.cohort()
        cohort[1]["x_auto_ratio"] = 1.00
        with self.assertRaisesRegex(ValueError, "overlaps"):
            calibrate_and_apply(cohort, "abc")

    def test_duplicate_sample_fails_closed(self):
        cohort = self.cohort()
        cohort[-1]["sample_id"] = cohort[0]["sample_id"]
        with self.assertRaisesRegex(ValueError, "unique"):
            calibrate_and_apply(cohort, "abc")


if __name__ == "__main__":
    unittest.main()
