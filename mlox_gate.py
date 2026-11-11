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

from mloy_gate import (
    SCORE_METRICS,
    IntervalReferences,
    build_control_sets,
    interval_residual_score,
    choose_control_set,
    load_metadata,
    parse_control_groups,
    select_depth_metric,
)


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
    *,
    training_excluded_sample_ids: set[str] | None = None,
    control_group: str = "TP53_negative",
    depth_metric: str = "raw",
    score_metric: str = "median_ratio",
) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    training_excluded_sample_ids = set(training_excluded_sample_ids or ())
    if score_metric not in SCORE_METRICS:
        raise ValueError(f"unknown score metric: {score_metric}")
    if collection.get("status") != "MLOX_INTERVAL_COLLECTION_COMPLETE":
        raise ValueError("mLOX interval collection is incomplete")
    control_groups, control_label = parse_control_groups(control_group)
    collection = select_depth_metric(collection, "eligible_x_intervals", depth_metric)
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
        age_value = (
            float(age)
            if isinstance(age, (int, float)) and math.isfinite(float(age)) and age >= 0
            else None
        )
        common = {
            "sample_id": sample_id,
            "analysis_id": analysis_id,
            "assay": sex_row["proposed_capture_family"],
            "batch_id": str(sex_row["batch_id"]),
            "batch_name": sex_row.get("batch_name"),
            "sex_reported": sex_row["sex_reported"],
            "sex_inferred": sex_row.get("sex_inferred"),
            "study_group": info.get("study_group"),
            "age_at_collection": age_value,
            "age_status": "AVAILABLE" if age_value is not None else "MISSING_RETAINED",
            "reference_stratum": info.get("reference_stratum"),
            "collection": collected[sample_id],
            "held_out_specificity": sample_id in training_excluded_sample_ids,
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
        "control_group": control_label,
        "normalization_rule": (
            "female and male baselines each from the most specific control set with >=5 "
            "eligible controls: batch+reference stratum, batch, reference stratum, else assay"
        ),
        "depth_metric": depth_metric,
        "score_metric": score_metric,
        "score_rule": (
            "mLRR_X = median(log2(interval_mean_depth/autosomal_mean_depth)) - female control median"
            if score_metric == "median_ratio" else
            "mLRR_X = median over panel intervals of log2(depth_over_autosome / female control "
            "interval median), sample left out of its own control set; f from the male/female "
            "baseline ratio under the same read-substitution model"
        ),
        "gc_correction": collection.get("gc_correction") if depth_metric == "gc" else None,
        "assays": {},
    }
    output: list[dict[str, Any]] = []
    for assay, assay_rows in sorted(by_assay.items()):
        reference_females = [
            r for r in assay_rows
            if r["mlox_reference_eligible"] and not r["held_out_specificity"]
        ]
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
            if (
                r["study_group"] in control_groups
                and r["mlox_reference_eligible"]
                and not r["held_out_specificity"]
            )
        ]
        male_controls = [
            r for r in reference_males
            if r["metric_status"] == "ELIGIBLE" and r["study_group"] in control_groups
        ]
        if len(female_controls) < 6 or len(male_controls) < 6:
            raise ValueError(f"{assay}: insufficient {control_label} sex-specific controls")
        female_control_sets = build_control_sets(female_controls, assay)
        references = IntervalReferences("eligible_x_intervals", selected)
        male_control_sets = build_control_sets(male_controls, assay)

        female_log_baselines = {
            key: _median(
                [row["raw_log2_x_auto_median"] for row in control_rows],
                f"{assay}/{key} female baseline",
            )
            for key, control_rows in female_control_sets.items()
        }
        male_log_baselines = {
            key: _median(
                [row["raw_log2_x_auto_median"] for row in control_rows],
                f"{assay}/{key} male baseline",
            )
            for key, control_rows in male_control_sets.items()
        }
        female_interval_baselines: dict[tuple[str, str], dict[str, float]] = {}
        for key, control_rows in female_control_sets.items():
            female_interval_baselines[key] = {
                interval_id: _median(
                    [
                        float(control["collection"]["eligible_x_intervals"][interval_id]["depth_over_autosome"])
                        for control in control_rows
                    ],
                    f"{assay}/{key}/{interval_id} female controls",
                )
                for interval_id in selected
            }

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
            female_key = choose_control_set(female_control_sets, row, assay)
            male_key = choose_control_set(male_control_sets, row, assay)
            female_reference = female_control_sets[female_key]
            male_reference = male_control_sets[male_key]
            female_baseline = female_log_baselines[female_key]
            male_baseline = male_log_baselines[male_key]
            female_linear = 2.0 ** female_baseline
            male_linear = 2.0 ** male_baseline
            observed_linear = 2.0 ** row["raw_log2_x_auto_median"]
            if not female_linear > male_linear > 0:
                raise ValueError(f"{assay}: invalid female/male X dosage baselines")
            row["normalization_scope"] = (
                f"FEMALE_{female_key[0]}_MALE_{male_key[0]}_{control_label.upper()}"
            )
            row["female_normalization_set"] = female_key[1]
            row["male_dosage_set"] = male_key[1]
            row["female_normalization_control_n"] = len(female_reference)
            row["male_dosage_control_n"] = len(male_reference)
            row["female_normalization_median"] = female_baseline
            row["male_dosage_median"] = male_baseline
            row["normalization_leave_one_out"] = False
            if score_metric == "interval_residual":
                interval_refs = references.medians(female_key, female_reference, row["sample_id"])
                score = interval_residual_score(row["collection"]["eligible_x_intervals"], interval_refs)
                row["mLRR_X"] = score
                row["normalization_leave_one_out"] = any(
                    r["sample_id"] == row["sample_id"] for r in female_reference
                )
                male_relative = male_linear / female_linear
                row["f_LOX_interpretive"] = max(0.0, min(1.0, (1.0 - 2.0 ** score) / (1.0 - male_relative)))
            else:
                row["mLRR_X"] = row["raw_log2_x_auto_median"] - female_baseline
                row["f_LOX_interpretive"] = max(
                    0.0, min(1.0, (female_linear - observed_linear) / (female_linear - male_linear))
                )

            residuals: list[float] = []
            below_zero = below_005 = below_010 = 0
            intervals = row["collection"]["eligible_x_intervals"]
            interval_references = female_interval_baselines[female_key]
            for interval_id in selected:
                observed = float(intervals[interval_id]["depth_over_autosome"])
                reference = interval_references[interval_id]
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

        age_eligible_controls = [
            row for row in female_controls if row["age_at_collection"] is not None
        ]
        if len(age_eligible_controls) < 6:
            raise ValueError(f"{assay}: insufficient age-eligible female controls")
        young_n = math.ceil(len(age_eligible_controls) / 3)
        young = sorted(
            age_eligible_controls,
            key=lambda r: (r["age_at_collection"], r["sample_id"]),
        )[:young_n]
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
            "female_controls_missing_age_n": len(female_controls) - len(age_eligible_controls),
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
        "age_status_counts": dict(sorted(Counter(r["age_status"] for r in clean_rows).items())),
        "control_group": control_label,
        "depth_metric": depth_metric,
        "score_metric": score_metric,
        "normalization_scope_counts": dict(sorted(
            Counter(r["normalization_scope"] for r in clean_rows if r.get("normalization_scope")).items()
        )),
        "metric_status_counts": dict(sorted(Counter(r["metric_status"] for r in clean_rows).items())),
        "screen_tier_counts": dict(sorted(tiers.items())),
        "retained_candidate_by_assay_and_group": {
            f"{assay}|{group}": count
            for (assay, group), count in sorted(Counter((r["assay"], r["study_group"]) for r in retained).items())
        },
        "event_calls_emitted": False,
        "negative_calls_emitted": False,
        "training_excluded_sample_n": len(training_excluded_sample_ids),
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
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--workbook", type=Path)
    source.add_argument("--metadata-tsv", type=Path,
                        help="TSV with analysis_id, age_at_collection, study_group")
    parser.add_argument("--control-group", default="TP53_negative",
                        help="study_group value used as the control reference")
    parser.add_argument("--depth-metric", choices=("raw", "gc"), default="raw",
                        help="gc uses depth_over_autosome_gc from a collection built with --gc-table")
    parser.add_argument("--score-metric", choices=SCORE_METRICS, default="median_ratio",
                        help="interval_residual scores each interval against its control median (leave-one-out)")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-sha256")
    args = parser.parse_args()
    observed = sha256_file(args.collection)
    if args.expected_sha256 and observed != args.expected_sha256.lower():
        raise SystemExit("collection SHA-256 mismatch")
    with gzip.open(args.collection, "rt", encoding="utf-8") as stream:
        collection = json.load(stream)
    sex_rows = json.loads(args.sex_cohort.read_text(encoding="utf-8"))
    metadata, metadata_key, metadata_path = load_metadata(args.workbook, args.metadata_tsv)
    rows, summary, panel = build_gate(
        collection, sex_rows, metadata, control_group=args.control_group, depth_metric=args.depth_metric,
        score_metric=args.score_metric,
    )
    summary.update({
        "collection_sha256": observed,
        "sex_cohort_sha256": sha256_file(args.sex_cohort),
        metadata_key: sha256_file(metadata_path),
    })
    panel.update({
        "collection_sha256": observed,
        "sex_cohort_sha256": summary["sex_cohort_sha256"],
        metadata_key: summary[metadata_key],
    })
    _write_json(args.output_dir / "mlox.metrics.PRIVATE.json", rows)
    _write_json(args.output_dir / "mlox.summary.json", summary)
    _write_json(args.output_dir / "mlox.callable_panel.json", panel)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
