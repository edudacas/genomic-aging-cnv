#!/usr/bin/env python3
"""Build a sensitive, non-diagnostic mLOX screening layer from X intervals.

Sex inference is context, never an exclusion criterion for reported females or
sex-chromosome anomalies. No negative or event call is emitted before BAM-level
downsampling calibrates sensitivity and technical target dropout.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


GATE_VERSION = "mlox-sensitive-retention-v1"
CALLABLE_BASE_FRACTION = 0.80
CALLABLE_SAMPLE_FRACTION = 0.95
CALLABLE_DEPTH_THRESHOLD = 20
MIN_RUN_CONTROLS = 5
MIN_SAMPLE_PANEL_COMPLETENESS = 0.95
EXTREME_MIN_INTERPRETIVE_FRACTION = 0.10
BORDERLINE_MIN_INTERPRETIVE_FRACTION = 0.05
WEAK_SIGNAL_SIGMA = 2.0
EXTREME_SIGNAL_SIGMA = 3.0


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _median(values: list[float], label: str) -> float:
    if not values:
        raise ValueError(f"empty median: {label}")
    return statistics.median(values)


def _load_workbook_metadata(path: Path) -> dict[str, dict[str, Any]]:
    import openpyxl

    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    sheet = workbook["Review"]
    iterator = sheet.iter_rows(values_only=True)
    header = next(iterator)
    required = {"analysis_id", "age_at_collection", "study_group"}
    if not required.issubset(header):
        raise ValueError("workbook metadata columns are missing")
    result: dict[str, dict[str, Any]] = {}
    for values in iterator:
        row = dict(zip(header, values))
        key = str(row["analysis_id"])
        if key in result:
            raise ValueError(f"duplicate analysis_id in workbook: {key}")
        result[key] = {
            "analysis_id": key,
            "age_at_collection": row["age_at_collection"],
            "study_group": row["study_group"],
        }
    workbook.close()
    return result


def _raw_metrics(intervals: dict[str, dict[str, Any]], selected: list[str]) -> dict[str, Any]:
    positive = [
        math.log2(float(intervals[interval_id]["depth_over_autosome"]))
        for interval_id in selected
        if float(intervals[interval_id]["depth_over_autosome"]) > 0
    ]
    numerator = sum(
        float(intervals[interval_id]["depth_over_autosome"])
        * int(intervals[interval_id]["territory_bases"])
        for interval_id in selected
    )
    denominator = sum(int(intervals[interval_id]["territory_bases"]) for interval_id in selected)
    weighted_ratio = numerator / denominator
    completeness = len(positive) / len(selected)
    if not positive:
        status = "X_SIGNAL_ABSENT_MANUAL_REVIEW"
        raw_median = None
    elif completeness < MIN_SAMPLE_PANEL_COMPLETENESS:
        status = "INCOMPLETE_X_PANEL_HOLD"
        raw_median = statistics.median(positive)
    else:
        status = "ELIGIBLE"
        raw_median = statistics.median(positive)
    return {
        "callable_interval_n": len(selected),
        "positive_interval_n": len(positive),
        "zero_interval_n": len(selected) - len(positive),
        "panel_completeness": completeness,
        "raw_log2_x_auto_median": raw_median,
        "raw_x_auto_weighted_ratio": weighted_ratio,
        "raw_log2_x_auto_weighted": math.log2(weighted_ratio) if weighted_ratio > 0 else None,
        "metric_status": status,
    }


def build_gate(
    collection: dict[str, Any],
    sex_rows: list[dict[str, Any]],
    metadata: dict[str, dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    if collection.get("status") != "MLOX_INTERVAL_COLLECTION_COMPLETE":
        raise ValueError("mLOX interval collection is incomplete")
    collected = {row["sample_id"]: row for row in collection["samples"]}
    sex = {row["sample_id"]: row for row in sex_rows}
    if len(collected) != collection.get("selected_samples") or set(collected) != set(sex):
        raise ValueError("sex-gate/collection sample mismatch")

    rows: list[dict[str, Any]] = []
    all_reference_males: list[dict[str, Any]] = []
    for sample_id, sex_row in sorted(sex.items()):
        analysis_id = str(sex_row.get("analysis_id"))
        if analysis_id not in metadata:
            raise ValueError(f"{sample_id}: workbook linkage missing")
        info = metadata[analysis_id]
        age = info.get("age_at_collection")
        if not isinstance(age, (int, float)) or not math.isfinite(float(age)) or age < 0:
            raise ValueError(f"{sample_id}: invalid age_at_collection")
        common = {
            "sample_id": sample_id,
            "analysis_id": analysis_id,
            "assay": sex_row["proposed_capture_family"],
            "batch_id": str(sex_row["batch_id"]),
            "batch_name": sex_row.get("batch_name"),
            "sex_reported": sex_row["sex_reported"],
            "sex_inferred": sex_row.get("sex_inferred"),
            "study_group": info.get("study_group"),
            "age_at_collection": float(age),
            "collection": collected[sample_id],
        }
        reported_female = sex_row.get("sex_reported") == "Female"
        inferred_female = sex_row.get("sex_inferred") == "Female"
        reported_male = sex_row.get("sex_reported") == "Male"
        inferred_male = sex_row.get("sex_inferred") == "Male"
        anomaly = sex_row.get("sex_gate_status") != "PASS" or sex_row.get("sex_inferred") is None
        if reported_female or inferred_female or anomaly:
            rows.append({
                **common,
                "mlox_reference_eligible": bool(
                    sex_row.get("sex_gate_status") == "PASS" and (reported_female or inferred_female)
                ),
                "sex_chromosome_anomaly_review": bool(anomaly),
            })
        if sex_row.get("sex_gate_status") == "PASS" and (reported_male or inferred_male):
            all_reference_males.append(common)
    if not rows:
        raise ValueError("no mLOX candidates or sex-chromosome anomaly reviews")

    by_assay: dict[str, list[dict[str, Any]]] = defaultdict(list)
    males_by_assay: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_assay[row["assay"]].append(row)
    for row in all_reference_males:
        males_by_assay[row["assay"]].append(row)

    panel_receipt: dict[str, Any] = {
        "schema_version": 1,
        "gate_version": GATE_VERSION,
        "callability_rule": {
            "per_sample_interval_rule": f"bases_ge_{CALLABLE_DEPTH_THRESHOLD}/territory_bases >= {CALLABLE_BASE_FRACTION}",
            "cohort_rule": f"per-sample rule true in >= {CALLABLE_SAMPLE_FRACTION} of reference-eligible females within assay",
            "sample_metric": "median(log2(interval_mean_depth/autosomal_mean_depth))",
            "minimum_sample_panel_completeness": MIN_SAMPLE_PANEL_COMPLETENESS,
        },
        "cell_fraction_model": "female*(1-f)+male*f on linear X/autosome depth",
        "retention_policy": "reported females and sex-chromosome anomalies are never excluded by sex inference",
        "event_call_policy": "BLOCKED_PENDING_BAM_READ_SUBSTITUTION_CALIBRATION",
        "assays": {},
    }
    output: list[dict[str, Any]] = []
    for assay, assay_rows in sorted(by_assay.items()):
        reference_females = [r for r in assay_rows if r["mlox_reference_eligible"]]
        reference_males = males_by_assay.get(assay, [])
        if not reference_females or not reference_males:
            raise ValueError(f"{assay}: missing reference females or males")
        interval_sets = [set(r["collection"]["eligible_x_intervals"]) for r in reference_females]
        if any(value != interval_sets[0] for value in interval_sets[1:]):
            raise ValueError(f"{assay}: interval membership differs across samples")
        interval_ids = sorted(interval_sets[0])
        selected: list[str] = []
        interval_stats: list[dict[str, Any]] = []
        for interval_id in interval_ids:
            values = [r["collection"]["eligible_x_intervals"][interval_id] for r in reference_females]
            territory = int(values[0]["territory_bases"])
            callable_n = sum(
                value[f"bases_ge_{CALLABLE_DEPTH_THRESHOLD}"] / territory >= CALLABLE_BASE_FRACTION
                for value in values
            )
            fraction = callable_n / len(values)
            is_selected = fraction >= CALLABLE_SAMPLE_FRACTION
            if is_selected:
                selected.append(interval_id)
            interval_stats.append({
                "interval_id": interval_id,
                "territory_bases": territory,
                "callable_reference_female_n": callable_n,
                "reference_female_n": len(values),
                "callable_reference_female_fraction": fraction,
                "selected": is_selected,
            })
        if not selected:
            raise ValueError(f"{assay}: callable X panel is empty")

        for row in assay_rows + reference_males:
            row.update(_raw_metrics(row["collection"]["eligible_x_intervals"], selected))
        eligible_females = [r for r in assay_rows if r["metric_status"] == "ELIGIBLE"]
        female_controls = [
            r for r in eligible_females
            if r["study_group"] == "TP53_negative" and r["mlox_reference_eligible"]
        ]
        male_controls = [
            r for r in reference_males
            if r["metric_status"] == "ELIGIBLE" and r["study_group"] == "TP53_negative"
        ]
        if len(female_controls) < 6 or len(male_controls) < 6:
            raise ValueError(f"{assay}: insufficient TP53-negative sex-specific controls")
        female_by_batch: dict[str, list[dict[str, Any]]] = defaultdict(list)
        male_by_batch: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in female_controls:
            female_by_batch[row["batch_id"]].append(row)
        for row in male_controls:
            male_by_batch[row["batch_id"]].append(row)

        for row in assay_rows:
            row.update({
                "normalization_scope": None,
                "female_normalization_control_n": None,
                "male_dosage_control_n": None,
                "mLRR_X": None,
                "f_LOX_interpretive": None,
                "target_residual_median": None,
                "fraction_targets_below_control": None,
                "fraction_targets_below_minus_0_05": None,
                "fraction_targets_below_minus_0_10": None,
            })
            if row["metric_status"] != "ELIGIBLE":
                continue
            female_reference = female_by_batch[row["batch_id"]]
            male_reference = male_by_batch[row["batch_id"]]
            if len(female_reference) < MIN_RUN_CONTROLS:
                female_reference = female_controls
            if len(male_reference) < MIN_RUN_CONTROLS:
                male_reference = male_controls
            female_baseline = _median(
                [r["raw_log2_x_auto_median"] for r in female_reference],
                f"{assay} female baseline",
            )
            male_baseline = _median(
                [r["raw_log2_x_auto_median"] for r in male_reference],
                f"{assay} male baseline",
            )
            female_linear = 2.0 ** female_baseline
            male_linear = 2.0 ** male_baseline
            observed_linear = 2.0 ** row["raw_log2_x_auto_median"]
            if not female_linear > male_linear > 0:
                raise ValueError(f"{assay}: invalid female/male X dosage baselines")
            row["normalization_scope"] = (
                "BATCH_SEX_SPECIFIC_TP53_NEGATIVE"
                if female_reference is female_by_batch[row["batch_id"]]
                and male_reference is male_by_batch[row["batch_id"]]
                else "ASSAY_SEX_SPECIFIC_TP53_NEGATIVE"
            )
            row["female_normalization_control_n"] = len(female_reference)
            row["male_dosage_control_n"] = len(male_reference)
            row["female_normalization_median"] = female_baseline
            row["male_dosage_median"] = male_baseline
            row["mLRR_X"] = row["raw_log2_x_auto_median"] - female_baseline
            row["f_LOX_interpretive"] = max(
                0.0, min(1.0, (female_linear - observed_linear) / (female_linear - male_linear))
            )

            residuals: list[float] = []
            below_zero = below_005 = below_010 = 0
            intervals = row["collection"]["eligible_x_intervals"]
            for interval_id in selected:
                observed = float(intervals[interval_id]["depth_over_autosome"])
                reference = _median([
                    float(value["collection"]["eligible_x_intervals"][interval_id]["depth_over_autosome"])
                    for value in female_reference
                ], f"{assay}/{interval_id} female controls")
                if reference <= 0:
                    raise ValueError(f"{assay}/{interval_id}: non-positive female control median")
                if observed > 0:
                    residual = math.log2(observed / reference)
                    residuals.append(residual)
                    below_zero += residual < 0
                    below_005 += residual < -0.05
                    below_010 += residual < -0.10
                else:
                    below_zero += 1
                    below_005 += 1
                    below_010 += 1
            row["target_residual_median"] = statistics.median(residuals) if residuals else None
            row["fraction_targets_below_control"] = below_zero / len(selected)
            row["fraction_targets_below_minus_0_05"] = below_005 / len(selected)
            row["fraction_targets_below_minus_0_10"] = below_010 / len(selected)

        young_n = math.ceil(len(female_controls) / 3)
        young = sorted(female_controls, key=lambda r: (r["age_at_collection"], r["sample_id"]))[:young_n]
        null_values = [r["mLRR_X"] for r in young]
        null_median = _median(null_values, f"{assay} young female null")
        null_mad = _median([abs(v - null_median) for v in null_values], f"{assay} young female MAD")
        robust_sigma = 1.4826 * null_mad
        if robust_sigma <= 0:
            raise ValueError(f"{assay}: degenerate young female null")
        weak_threshold = null_median - WEAK_SIGNAL_SIGMA * robust_sigma
        extreme_threshold = null_median - EXTREME_SIGNAL_SIGMA * robust_sigma

        for row in assay_rows:
            row["young_null_median"] = null_median
            row["young_null_robust_sigma"] = robust_sigma
            row["weak_signal_threshold_2sigma"] = weak_threshold
            row["extreme_signal_threshold_3sigma"] = extreme_threshold
            row["event_call_LOX"] = "NOT_CALIBRATED"
            if row["metric_status"] == "X_SIGNAL_ABSENT_MANUAL_REVIEW":
                tier = "X_SIGNAL_ABSENT_MANUAL_REVIEW"
            elif row["metric_status"] == "INCOMPLETE_X_PANEL_HOLD":
                tier = "TECHNICAL_NOT_EVALUABLE_RETAINED"
            elif row["sex_chromosome_anomaly_review"]:
                tier = "SEX_CHROMOSOME_ANOMALY_REVIEW"
            elif row["metric_status"] == "ELIGIBLE":
                if (
                    row["mLRR_X"] < extreme_threshold
                    and row["f_LOX_interpretive"] >= EXTREME_MIN_INTERPRETIVE_FRACTION
                ):
                    tier = "EXTREME_CANDIDATE_UNCALIBRATED"
                elif row["f_LOX_interpretive"] >= BORDERLINE_MIN_INTERPRETIVE_FRACTION:
                    tier = "BORDERLINE_CANDIDATE_UNCALIBRATED"
                elif row["mLRR_X"] < weak_threshold:
                    tier = "WEAK_SIGNAL_RETAINED"
                else:
                    tier = "WITHIN_REFERENCE_UNCALIBRATED"
            else:
                tier = "NOT_EVALUABLE_RETAINED"
            row["screen_tier_LOX"] = tier

        panel_receipt["assays"][assay] = {
            "analysis_candidate_n": len(assay_rows),
            "reference_female_n": len(reference_females),
            "reference_male_n": len(reference_males),
            "eligible_female_control_n": len(female_controls),
            "eligible_male_control_n": len(male_controls),
            "selected_interval_n": len(selected),
            "selected_territory_bases": sum(
                next(item["territory_bases"] for item in interval_stats if item["interval_id"] == interval_id)
                for interval_id in selected
            ),
            "young_null_n": len(young),
            "young_null_max_age": max(r["age_at_collection"] for r in young),
            "young_null_median": null_median,
            "young_null_robust_sigma": robust_sigma,
            "weak_signal_threshold_2sigma": weak_threshold,
            "extreme_signal_threshold_3sigma": extreme_threshold,
            "intervals": interval_stats,
        }
        output.extend(assay_rows)

    clean_rows = []
    for row in sorted(output, key=lambda value: value["sample_id"]):
        clean = {key: value for key, value in row.items() if key != "collection"}
        clean["mlox_gate_version"] = GATE_VERSION
        clean_rows.append(clean)
    tiers = Counter(r["screen_tier_LOX"] for r in clean_rows)
    retained = [
        r for r in clean_rows
        if r["screen_tier_LOX"] in {
            "EXTREME_CANDIDATE_UNCALIBRATED",
            "BORDERLINE_CANDIDATE_UNCALIBRATED",
            "WEAK_SIGNAL_RETAINED",
        }
    ]
    summary = {
        "status": "MLOX_GATE_METRICS_COMPLETE_EVENT_CALLS_BLOCKED",
        "gate_version": GATE_VERSION,
        "mlox_analysis_samples": len(clean_rows),
        "reference_female_samples": sum(r["mlox_reference_eligible"] for r in clean_rows),
        "sex_chromosome_anomaly_review_samples": sum(r["sex_chromosome_anomaly_review"] for r in clean_rows),
        "metric_status_counts": dict(sorted(Counter(r["metric_status"] for r in clean_rows).items())),
        "screen_tier_counts": dict(sorted(tiers.items())),
        "retained_candidate_by_assay_and_group": {
            f"{assay}|{group}": count
            for (assay, group), count in sorted(Counter((r["assay"], r["study_group"]) for r in retained).items())
        },
        "event_calls_emitted": False,
        "negative_calls_emitted": False,
        "borderline_signals_retained": True,
        "screening_is_not_an_event_call": True,
        "calibration_required": "BAM-level X-read substitution by assay, depth and fraction",
        "collection_mask_sha256": collection["mask_sha256"],
    }
    return clean_rows, summary, panel_receipt


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, ensure_ascii=False)
        stream.write("\n")
    os.chmod(path, 0o600)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--collection", type=Path, required=True)
    parser.add_argument("--sex-cohort", type=Path, required=True)
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-sha256")
    args = parser.parse_args()
    observed = sha256_file(args.collection)
    if args.expected_sha256 and observed != args.expected_sha256.lower():
        raise SystemExit("collection SHA-256 mismatch")
    with gzip.open(args.collection, "rt", encoding="utf-8") as stream:
        collection = json.load(stream)
    sex_rows = json.loads(args.sex_cohort.read_text(encoding="utf-8"))
    metadata = _load_workbook_metadata(args.workbook)
    rows, summary, panel = build_gate(collection, sex_rows, metadata)
    summary.update({
        "collection_sha256": observed,
        "sex_cohort_sha256": sha256_file(args.sex_cohort),
        "workbook_sha256": sha256_file(args.workbook),
    })
    panel.update({
        "collection_sha256": observed,
        "sex_cohort_sha256": summary["sex_cohort_sha256"],
        "workbook_sha256": summary["workbook_sha256"],
    })
    _write_json(args.output_dir / "mlox.metrics.PRIVATE.json", rows)
    _write_json(args.output_dir / "mlox.summary.json", summary)
    _write_json(args.output_dir / "mlox.callable_panel.json", panel)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
