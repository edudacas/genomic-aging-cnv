"""GRCm38 masks, configurable control group and TSV metadata for the mLOY/mLOX path."""

import csv
import gzip
import hashlib
import json
import math
import tempfile
import unittest
from pathlib import Path

import collect_mlox_intervals
import collect_mloy_intervals
from collect_mloy_intervals import gc_log2_correction, validate_mask
from interval_gc import gc_fraction
from mloy_gate import (
    build_control_sets,
    build_gate as build_mloy_gate,
    choose_control_set,
    load_metadata,
    select_depth_metric,
)
from mlox_gate import build_gate as build_mlox_gate
try:  # `python -m unittest discover -s tests` puts tests/ on sys.path
    import test_collect_mloy_intervals as collector_fixtures
    import test_mloy_gate as mloy_fixtures
    import test_mlox_gate as mlox_fixtures
except ImportError:  # run as tests.test_mouse_grcm38 from the repository root
    from tests import test_collect_mloy_intervals as collector_fixtures
    from tests import test_mloy_gate as mloy_fixtures
    from tests import test_mlox_gate as mlox_fixtures

# Module references (not class imports) so unittest does not collect those suites twice.
FIELDS, interval = collector_fixtures.FIELDS, collector_fixtures.interval

RESOURCE = Path(__file__).resolve().parents[1] / "resources" / "sex_compartments.GRCm38.v1.json"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def mouse_mask(**extra):
    mask = {
        "assembly": "GRCm38", "coordinate_system": "0-based-half-open",
        "par": {"chrX": [[900, 1000]], "chrY": [[900, 1000]]},
        "xtr": {"chrX": [], "chrY": [],
                "status": "FROZEN_FOR_FUTURE_RECOMPUTATION_NOT_USED_IN_CURRENT_SEX_GATE"},
        "additional_exclusions": {"chrY": [{"start": 100, "end": 900, "label": "AMPLICONIC_YQ"}]},
    }
    mask.update(extra)
    return mask


class MouseResourceTests(unittest.TestCase):
    def test_shipped_resource_is_valid(self):
        mask = json.loads(RESOURCE.read_text())
        autosomes = validate_mask(mask)
        self.assertEqual(len(autosomes), 19)
        self.assertNotIn("chr20", autosomes)
        y = collect_mloy_intervals.mask_exclusions(mask, "chrY")
        self.assertEqual([label for *_, label in y], ["AMPLICONIC_YQ", "PAR"])
        x = collect_mloy_intervals.mask_exclusions(mask, "chrX")
        self.assertEqual([label for *_, label in x], ["PAR"])

    def test_unknown_assembly_and_wrong_autosomes_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "assembly"):
            validate_mask(mouse_mask(assembly="GRCm39"))
        with self.assertRaisesRegex(ValueError, "autosomes"):
            validate_mask(mouse_mask(autosomes=[f"chr{i}" for i in range(1, 23)]))

    def test_additional_exclusion_label_must_be_safe(self):
        mask = mouse_mask(additional_exclusions={"chrY": [{"start": 1, "end": 2, "label": "PAR"}]})
        with self.assertRaisesRegex(ValueError, "label"):
            collect_mloy_intervals.mask_exclusions(mask, "chrY")


class MouseCollectorTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.stages = self.root / "stages"
        self.mask = self.root / "mask.json"
        self.mask.write_text(json.dumps(mouse_mask()))

    def tearDown(self):
        self.tmp.cleanup()

    def make_stage(self):
        directory = self.stages / "qc" / "M1"
        directory.mkdir(parents=True)
        rows = [
            interval("chr19", 0, 100, 0, 100, 100),
            interval("chr20", 0, 100, 0, 100, 900),   # not a mouse autosome: ignored
            interval("chrX", 10, 60, 10, 60, 50),
            interval("chrX", 950, 960, 950, 960, 50),  # PAR
            interval("chrY", 10, 60, 10, 60, 50),
            interval("chrY", 200, 260, 200, 260, 5),  # ampliconic Yq
        ]
        path = directory / "coverage.intervals.tsv.gz"
        with gzip.open(path, "wt", newline="") as stream:
            writer = csv.DictWriter(stream, FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)
        complete = directory / "COMPLETE.json"
        complete.write_text(json.dumps({"status": "COMPLETE", "fingerprint": "fp",
                                        "outputs": {"coverage.intervals.tsv.gz": sha(path)}}))
        (self.stages / "QC_COMPLETE.json").write_text(json.dumps({
            "status": "QC_COMPLETE", "fingerprint": "fp",
            "samples": {"M1": {"receipt_sha256": sha(complete)}}}))

    def test_y_collection_uses_mouse_autosomes_and_excludes_yq(self):
        self.make_stage()
        out = self.root / "y.json.gz"
        collect_mloy_intervals.collect(self.stages, self.mask, out, 1, 1)
        payload = json.loads(gzip.decompress(out.read_bytes()))
        sample = payload["samples"][0]
        self.assertEqual(sample["autosomal"]["mean_depth"], 100)
        self.assertEqual(set(sample["eligible_y_intervals"]), {"chrY:10-60"})
        status = {r["interval_id"]: r["status"] for r in payload["y_interval_catalog"]}
        self.assertEqual(status["chrY:200-260"], "EXCLUDED_AMPLICONIC_YQ")

    def test_x_collection_excludes_mouse_par(self):
        self.make_stage()
        out = self.root / "x.json.gz"
        collect_mlox_intervals.collect(self.stages, self.mask, out, 1, 1)
        payload = json.loads(gzip.decompress(out.read_bytes()))
        self.assertEqual(set(payload["samples"][0]["eligible_x_intervals"]), {"chrX:10-60"})
        self.assertEqual(payload["samples"][0]["autosomal"]["mean_depth"], 100)


class GcCorrectionTests(unittest.TestCase):
    """Autosomes at 40% GC read 2x deeper than at 60% GC; sex-chromosome depth is corrected."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.stages = self.root / "stages"
        self.mask = self.root / "mask.json"
        self.mask.write_text(json.dumps(mouse_mask()))
        self.gc = self.root / "gc.tsv"

    def tearDown(self):
        self.tmp.cleanup()

    def make_stage(self):
        rows, gc_rows = [], []
        for i in range(600):
            start = i * 100
            gc, depth = (0.40, 100) if i % 2 == 0 else (0.60, 50)
            rows.append(interval("chr1", start, start + 100, start, start + 100, depth))
            gc_rows.append(("chr1", start, start + 100, gc))
        for chrom in ("chrX", "chrY"):
            rows.append(interval(chrom, 10, 50, 10, 50, 100))   # 40% GC, inflated
            rows.append(interval(chrom, 50, 90, 50, 90, 50))    # 60% GC, depressed
            gc_rows += [(chrom, 10, 50, 0.40), (chrom, 50, 90, 0.60)]
        directory = self.stages / "qc" / "M1"
        directory.mkdir(parents=True)
        path = directory / "coverage.intervals.tsv.gz"
        with gzip.open(path, "wt", newline="") as stream:
            writer = csv.DictWriter(stream, FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)
        complete = directory / "COMPLETE.json"
        complete.write_text(json.dumps({"status": "COMPLETE", "fingerprint": "fp",
                                        "outputs": {"coverage.intervals.tsv.gz": sha(path)}}))
        (self.stages / "QC_COMPLETE.json").write_text(json.dumps({
            "status": "QC_COMPLETE", "fingerprint": "fp",
            "samples": {"M1": {"receipt_sha256": sha(complete)}}}))
        self.gc.write_text("chrom\tstart\tend\tgc\n" + "".join(f"{c}\t{a}\t{b}\t{g}\n" for c, a, b, g in gc_rows))

    def test_collectors_equalise_gc_biased_intervals(self):
        self.make_stage()
        for module, key in ((collect_mloy_intervals, "eligible_y_intervals"),
                            (collect_mlox_intervals, "eligible_x_intervals")):
            out = self.root / f"{key}.json.gz"
            module.collect(self.stages, self.mask, out, 1, 1, self.gc)
            payload = json.loads(gzip.decompress(out.read_bytes()))
            self.assertIn("gc_correction", payload)
            values = payload["samples"][0][key]
            low, high = values[f"chr{key[9].upper()}:10-50"], values[f"chr{key[9].upper()}:50-90"]
            self.assertAlmostEqual(low["depth_over_autosome"] / high["depth_over_autosome"], 2.0)
            self.assertAlmostEqual(low["depth_over_autosome_gc"], high["depth_over_autosome_gc"])
            self.assertAlmostEqual(low["gc_fraction"], 0.40)
            self.assertAlmostEqual(payload["samples"][0]["gc_curve"]["range_log2"], 1.0)

    def test_without_gc_table_output_is_unchanged(self):
        self.make_stage()
        out = self.root / "y.json.gz"
        collect_mloy_intervals.collect(self.stages, self.mask, out, 1, 1)
        payload = json.loads(gzip.decompress(out.read_bytes()))
        self.assertNotIn("gc_correction", payload)
        self.assertNotIn("gc_curve", payload["samples"][0])
        self.assertNotIn("depth_over_autosome_gc", payload["samples"][0]["eligible_y_intervals"]["chrY:10-50"])

    def test_sex_interval_missing_from_gc_table_fails(self):
        self.make_stage()
        self.gc.write_text("\n".join(l for l in self.gc.read_text().splitlines() if not l.startswith("chrY\t50")) + "\n")
        with self.assertRaisesRegex(ValueError, "GC table lacks"):
            collect_mloy_intervals.collect(self.stages, self.mask, self.root / "y.json.gz", 1, 1, self.gc)

    def test_curve_interpolates_and_is_flat_outside(self):
        curve = ([0.4, 0.6], [0.5, -0.5])
        self.assertEqual(gc_log2_correction(curve, 0.3), 0.5)
        self.assertEqual(gc_log2_correction(curve, 0.9), -0.5)
        self.assertAlmostEqual(gc_log2_correction(curve, 0.5), 0.0)

    def test_gc_fraction_ignores_n(self):
        self.assertAlmostEqual(gc_fraction("ggccNNatAT"), 0.5)
        self.assertIsNone(gc_fraction("NNNN"))

    def test_gate_gc_metric_requires_corrected_collection(self):
        collection, sex, metadata = mloy_fixtures.MloyGateTests().fixture()
        with self.assertRaisesRegex(ValueError, "no GC correction"):
            build_mloy_gate(collection, sex, metadata, depth_metric="gc")
        collection["gc_correction"] = {"method": "test"}
        for sample in collection["samples"]:
            for values in sample["eligible_y_intervals"].values():
                values["depth_over_autosome_gc"] = values["depth_over_autosome"] * 2
        swapped = select_depth_metric(collection, "eligible_y_intervals", "gc")
        first = swapped["samples"][0]["eligible_y_intervals"]["chrY:10-110"]
        self.assertEqual(first["depth_over_autosome"], 2 * first["depth_over_autosome_raw"])
        self.assertNotIn("depth_over_autosome_raw", collection["samples"][0]["eligible_y_intervals"]["chrY:10-110"])
        rows, _, panel = build_mloy_gate(collection, sex, metadata, depth_metric="gc")
        self.assertEqual(panel["depth_metric"], "gc")


class ReferenceStratumTests(unittest.TestCase):
    def test_control_set_order(self):
        controls = [{"batch_id": "B1", "reference_stratum": "embryo"} for _ in range(5)]
        controls += [{"batch_id": "B1", "reference_stratum": "postnatal"} for _ in range(3)]
        controls += [{"batch_id": "B2", "reference_stratum": "postnatal"} for _ in range(3)]
        sets = build_control_sets(controls, "A")
        self.assertEqual(choose_control_set(sets, {"batch_id": "B1", "reference_stratum": "embryo"}, "A"),
                         ("BATCH_STRATUM", "B1|embryo"))
        self.assertEqual(choose_control_set(sets, {"batch_id": "B1", "reference_stratum": "postnatal"}, "A"),
                         ("BATCH", "B1"))
        self.assertEqual(choose_control_set(sets, {"batch_id": "B2", "reference_stratum": "postnatal"}, "A"),
                         ("STRATUM", "postnatal"))
        self.assertEqual(choose_control_set(sets, {"batch_id": "B2", "reference_stratum": None}, "A"),
                         ("ASSAY", "A"))

    def test_stratum_offset_is_removed(self):
        collection, sex, metadata = mloy_fixtures.MloyGateTests().fixture()
        # Six more controls in an "embryo" stratum sitting 10% lower than postnatal controls.
        for index, depth in zip(range(9, 15), (0.45, 0.44, 0.46, 0.43, 0.47, 0.45)):
            sid = f"S{index}"
            collection["samples"].append({"sample_id": sid, "eligible_y_intervals": {
                "chrY:10-110": mloy_fixtures.interval(depth),
                "chrY:200-300": mloy_fixtures.interval(depth, covered_fraction=0.0)}})
            sex.append({**sex[0], "sample_id": sid, "analysis_id": str(index)})
            metadata[str(index)] = {"age_at_collection": 1, "study_group": "TP53_negative",
                                    "reference_stratum": "embryo"}
        collection["selected_samples"] = len(collection["samples"])
        for key in map(str, range(9)):
            metadata[key]["reference_stratum"] = "postnatal"
        rows, summary, _ = build_mloy_gate(collection, sex, metadata)
        embryo = next(r for r in rows if r["sample_id"] == "S9")
        self.assertEqual(embryo["normalization_scope"], "BATCH_STRATUM_TP53_NEGATIVE_MALES")
        self.assertAlmostEqual(embryo["mLRR_Y"], 0.0)
        self.assertIn("BATCH_STRATUM_TP53_NEGATIVE_MALES", summary["normalization_scope_counts"])


class IntervalResidualScoreTests(unittest.TestCase):
    """Intervals with different capture efficiency: the residual score removes it."""

    def fixture(self):
        capture = {"chrY:0-100": 0.3, "chrY:100-200": 0.6, "chrY:200-300": 1.2}
        sample_scale = [0.96, 0.98, 1.00, 1.02, 1.04, 0.99, 1.01, 0.97]
        collection, sex, metadata = mloy_fixtures.MloyGateTests().fixture()
        for index, sample in enumerate(collection["samples"]):
            scale = sample_scale[index] if index < 8 else 0.8 * 1.0   # S8: 20% of cells lost Y
            noise = (1.0, 1.0, 1.0) if index >= 8 else (1.0 + 0.01 * (index % 3), 1.0, 1.0 - 0.01 * (index % 2))
            sample["eligible_y_intervals"] = {
                interval_id: mloy_fixtures.interval(0.5 * factor * scale * n)
                for (interval_id, factor), n in zip(capture.items(), noise)
            }
        return collection, sex, metadata

    def test_mloy_residual_score_with_leave_one_out(self):
        rows, summary, panel = build_mloy_gate(*self.fixture(), score_metric="interval_residual")
        self.assertEqual(panel["score_metric"], "interval_residual")
        self.assertEqual(summary["score_metric"], "interval_residual")
        case = next(r for r in rows if r["sample_id"] == "S8")
        self.assertAlmostEqual(case["mLRR_Y"], math.log2(0.8), delta=0.02)
        self.assertAlmostEqual(case["f_LOY_interpretive"], 0.2, delta=0.015)
        self.assertFalse(case["normalization_leave_one_out"])
        control = next(r for r in rows if r["sample_id"] == "S0")
        self.assertTrue(control["normalization_leave_one_out"])
        self.assertEqual(case["screen_tier_LOY"], "EXTREME_CANDIDATE_UNCALIBRATED")

    def test_default_score_is_unchanged(self):
        rows_default, _, panel = build_mloy_gate(*self.fixture())
        rows_explicit, _, _ = build_mloy_gate(*self.fixture(), score_metric="median_ratio")
        self.assertEqual(panel["score_metric"], "median_ratio")
        self.assertEqual([r["mLRR_Y"] for r in rows_default], [r["mLRR_Y"] for r in rows_explicit])

    def test_unknown_score_metric_fails(self):
        with self.assertRaisesRegex(ValueError, "score metric"):
            build_mloy_gate(*self.fixture(), score_metric="mean")

    def test_mlox_residual_score_matches_read_substitution_fraction(self):
        collection, sex, metadata = mlox_fixtures.MloxGateTests().fixture()
        ratio_rows, _, _ = build_mlox_gate(collection, sex, metadata)
        rows, _, panel = build_mlox_gate(collection, sex, metadata, score_metric="interval_residual")
        self.assertEqual(panel["score_metric"], "interval_residual")
        case = next(r for r in rows if r["sample_id"] == "CASE")
        ratio_case = next(r for r in ratio_rows if r["sample_id"] == "CASE")
        self.assertAlmostEqual(case["f_LOX_interpretive"], ratio_case["f_LOX_interpretive"], delta=0.02)
        self.assertFalse(case["normalization_leave_one_out"])
        self.assertTrue(next(r for r in rows if r["sample_id"] == "S0")["normalization_leave_one_out"])


class ControlGroupTests(unittest.TestCase):
    @staticmethod
    def relabel(metadata, old, new):
        for info in metadata.values():
            if info["study_group"] == old:
                info["study_group"] = new
        return metadata

    def test_mloy_gate_uses_named_control_group(self):
        collection, sex, metadata = mloy_fixtures.MloyGateTests().fixture()
        self.relabel(metadata, "TP53_negative", "WT")
        with self.assertRaisesRegex(ValueError, "TP53_negative"):
            build_mloy_gate(collection, sex, metadata)
        rows, summary, panel = build_mloy_gate(collection, sex, metadata, control_group="WT")
        case = next(r for r in rows if r["study_group"] == "LFS_R337H")
        self.assertEqual(case["screen_tier_LOY"], "EXTREME_CANDIDATE_UNCALIBRATED")
        self.assertEqual(case["normalization_scope"], "BATCH_WT_MALES")
        self.assertEqual(panel["control_group"], "WT")

    def test_several_control_groups(self):
        collection, sex, metadata = mloy_fixtures.MloyGateTests().fixture()
        for index, info in enumerate(metadata.values()):
            if info["study_group"] == "TP53_negative":
                info["study_group"] = "WT" if index % 2 else "Het"
        with self.assertRaisesRegex(ValueError, "insufficient"):
            build_mloy_gate(collection, sex, metadata, control_group="WT")
        rows, _, panel = build_mloy_gate(collection, sex, metadata, control_group="WT, Het")
        self.assertEqual(panel["control_group"], "Het+WT")
        case = next(r for r in rows if r["study_group"] == "LFS_R337H")
        self.assertEqual(case["normalization_control_n"], 8)

    def test_mloy_gate_retains_missing_age_outside_young_null(self):
        collection, sex, metadata = mloy_fixtures.MloyGateTests().fixture()
        metadata["8"]["age_at_collection"] = None
        rows, summary, _ = build_mloy_gate(collection, sex, metadata)
        case = next(r for r in rows if r["sample_id"] == "S8")
        self.assertEqual(case["age_status"], "MISSING_RETAINED")
        self.assertEqual(case["screen_tier_LOY"], "EXTREME_CANDIDATE_UNCALIBRATED")
        self.assertEqual(summary["missing_age_sample_n"], 1)

    def test_mlox_gate_uses_named_control_group(self):
        collection, sex, metadata = mlox_fixtures.MloxGateTests().fixture()
        self.relabel(metadata, "TP53_negative", "WT")
        rows, _, panel = build_mlox_gate(collection, sex, metadata, control_group="WT")
        case = next(r for r in rows if r["sample_id"] == "CASE")
        self.assertTrue(case["normalization_scope"].endswith("_WT"))
        self.assertEqual(panel["control_group"], "WT")


class MetadataTsvTests(unittest.TestCase):
    def test_tsv_metadata_with_blank_age(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "meta.tsv"
            path.write_text("analysis_id\tage_at_collection\tstudy_group\treference_stratum\n"
                            "A\t79.5\tWT\tpostnatal\nB\t\tHet\t\n")
            metadata, key, source = load_metadata(None, path)
        self.assertEqual(key, "metadata_tsv_sha256")
        self.assertEqual(metadata["A"]["age_at_collection"], 79.5)
        self.assertEqual(metadata["A"]["reference_stratum"], "postnatal")
        self.assertIsNone(metadata["B"]["age_at_collection"])
        self.assertIsNone(metadata["B"]["reference_stratum"])

    def test_exactly_one_metadata_source(self):
        with self.assertRaisesRegex(ValueError, "exactly one"):
            load_metadata(None, None)


if __name__ == "__main__":
    unittest.main()
