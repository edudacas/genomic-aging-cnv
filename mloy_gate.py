#!/usr/bin/env python3
"""Freeze assay-specific callable-Y panels and calculate exploratory mLRR-Y.

Event calls remain blocked until BAM-level downsampling calibrates sensitivity.
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


GATE_VERSION = "mloy-sensitive-retention-v2"
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
    return result


def _load_metadata_tsv(path: Path) -> dict[str, dict[str, Any]]:
    """Tab-separated alternative to the workbook: analysis_id, age_at_collection, study_group.

    A blank age is kept as missing (the sample is retained but never enters the young null).
    An optional reference_stratum column (for example embryo / postnatal) enables
    stratum-matched control baselines; blank means no stratum.
    """
    import csv

    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"analysis_id", "age_at_collection", "study_group"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("metadata TSV columns are missing")
        result: dict[str, dict[str, Any]] = {}
        for row in reader:
            key = str(row["analysis_id"])
            if key in result:
                raise ValueError(f"duplicate analysis_id in metadata TSV: {key}")
            age_text = (row["age_at_collection"] or "").strip()
            result[key] = {
                "analysis_id": key,
                "age_at_collection": float(age_text) if age_text else None,
                "study_group": row["study_group"],
                "reference_stratum": (row.get("reference_stratum") or "").strip() or None,
            }
    return result


def parse_control_groups(control_group: str) -> tuple[frozenset[str], str]:
    """Comma-separated study_group values that form the control reference."""
    groups = frozenset(part.strip() for part in control_group.split(",") if part.strip())
    if not groups:
        raise ValueError("empty control group")
    return groups, "+".join(sorted(groups))


def select_depth_metric(collection: dict[str, Any], intervals_key: str, depth_metric: str) -> dict[str, Any]:
    """Return a collection whose depth_over_autosome is the raw or GC-corrected value."""
    if depth_metric == "raw":
        return collection
    if depth_metric != "gc":
        raise ValueError(f"unknown depth metric: {depth_metric}")
    if not collection.get("gc_correction"):
        raise ValueError("collection has no GC correction; rerun the collector with --gc-table")
    samples = []
    for sample in collection["samples"]:
        intervals = {}
        for interval_id, values in sample[intervals_key].items():
            if "depth_over_autosome_gc" not in values:
                raise ValueError(f"{sample['sample_id']}/{interval_id}: GC-corrected depth missing")
            intervals[interval_id] = {
                **values,
                "depth_over_autosome_raw": values["depth_over_autosome"],
                "depth_over_autosome": values["depth_over_autosome_gc"],
            }
        samples.append({**sample, intervals_key: intervals})
    return {**collection, "samples": samples}


def build_control_sets(controls: list[dict[str, Any]], assay: str) -> dict[tuple[str, str], list[dict[str, Any]]]:
    """Control sets with at least MIN_RUN_CONTROLS members; the assay-wide set always exists.

    Keys: ("BATCH_STRATUM", "batch|stratum"), ("BATCH", batch), ("STRATUM", stratum), ("ASSAY", assay).
    """
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in controls:
        grouped[("BATCH", row["batch_id"])].append(row)
        stratum = row.get("reference_stratum")
        if stratum:
            grouped[("BATCH_STRATUM", f"{row['batch_id']}|{stratum}")].append(row)
            grouped[("STRATUM", stratum)].append(row)
    sets = {key: rows for key, rows in grouped.items() if len(rows) >= MIN_RUN_CONTROLS}
    sets[("ASSAY", assay)] = controls
    return sets


def choose_control_set(
    sets: dict[tuple[str, str], list[dict[str, Any]]], row: dict[str, Any], assay: str
) -> tuple[str, str]:
    """Most specific available set: batch+stratum, batch, stratum, then assay."""
    stratum = row.get("reference_stratum")
    order = []
    if stratum:
        order.append(("BATCH_STRATUM", f"{row['batch_id']}|{stratum}"))
    order.append(("BATCH", row["batch_id"]))
    if stratum:
        order.append(("STRATUM", stratum))
    order.append(("ASSAY", assay))
    return next(key for key in order if key in sets)


SCORE_METRICS = ("median_ratio", "interval_residual")


class IntervalReferences:
    """Per-interval control medians of depth_over_autosome, cached per control set.

    For the interval_residual score the sample itself is left out of its own reference.
    """

    def __init__(self, intervals_key: str, selected: list[str]):
        self.intervals_key = intervals_key
        self.selected = selected
        self.cache: dict[tuple[tuple[str, str], str | None], dict[str, float]] = {}

    def medians(self, key: tuple[str, str], controls: list[dict[str, Any]], exclude: str | None) -> dict[str, float]:
        members = [row for row in controls if row["sample_id"] != exclude]
        cache_key = (key, exclude if len(members) < len(controls) else None)
        if cache_key not in self.cache:
            if not members:
                raise ValueError(f"{key}: empty interval reference")
            result = {}
            for interval_id in self.selected:
                value = statistics.median(
                    float(row["collection"][self.intervals_key][interval_id]["depth_over_autosome"])
                    for row in members
                )
                if value <= 0:
                    raise ValueError(f"{key}/{interval_id}: non-positive control median")
                result[interval_id] = value
            self.cache[cache_key] = result
        return self.cache[cache_key]


def interval_residual_score(values: dict[str, dict[str, Any]], references: dict[str, float]) -> float | None:
    """Median over panel intervals of log2(observed / control median); zero-depth intervals skipped."""
    residuals = [
        math.log2(float(values[interval_id]["depth_over_autosome"]) / reference)
        for interval_id, reference in references.items()
        if float(values[interval_id]["depth_over_autosome"]) > 0
    ]
    return statistics.median(residuals) if residuals else None


def load_metadata(workbook: Path | None, metadata_tsv: Path | None) -> tuple[dict[str, dict[str, Any]], str, Path]:
    """Load exactly one metadata source; returns (metadata, receipt key, path)."""
    if (workbook is None) == (metadata_tsv is None):
        raise ValueError("provide exactly one of --workbook or --metadata-tsv")
    if workbook is not None:
        return _load_workbook_metadata(workbook), "workbook_sha256", workbook
    return _load_metadata_tsv(metadata_tsv), "metadata_tsv_sha256", metadata_tsv


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
    if collection.get("status") != "MLOY_INTERVAL_COLLECTION_COMPLETE":
        raise ValueError("mLOY interval collection is incomplete")
    control_groups, control_label = parse_control_groups(control_group)
    collection = select_depth_metric(collection, "eligible_y_intervals", depth_metric)
    collected = {row["sample_id"]: row for row in collection["samples"]}
    if len(collected) != collection.get("selected_samples"):
        raise ValueError("collection membership mismatch")
    sex = {row["sample_id"]: row for row in sex_rows}
    if len(sex) != len(sex_rows) or set(sex) != set(collected):
        raise ValueError("sex-gate/collection sample mismatch")

    candidates: list[dict[str, Any]] = []
    for sample_id, row in sorted(sex.items()):
        reported_male = row.get("sex_reported") == "Male"
        inferred_male = row.get("sex_inferred") == "Male"
        anomaly_review = row.get("sex_gate_status") != "PASS" or row.get("sex_inferred") is None
        if not (reported_male or inferred_male or anomaly_review):
            continue
        analysis_id = str(row.get("analysis_id"))
        if analysis_id not in metadata:
            raise ValueError(f"{sample_id}: workbook linkage missing")
        age = metadata[analysis_id].get("age_at_collection")
        if age is None:
            age_value = None
        elif not isinstance(age, (int, float)) or not math.isfinite(float(age)) or age < 0:
            raise ValueError(f"{sample_id}: invalid age_at_collection")
        else:
            age_value = float(age)
        group = metadata[analysis_id].get("study_group")
        if not group:
            raise ValueError(f"{sample_id}: study_group missing")
        candidates.append(
            {
                "sample_id": sample_id,
                "analysis_id": analysis_id,
                "assay": row["proposed_capture_family"],
                "batch_id": str(row["batch_id"]),
                "batch_name": row.get("batch_name"),
                "sex_reported": row["sex_reported"],
                "sex_inferred": row["sex_inferred"],
                "study_group": group,
                "age_at_collection": age_value,
                "age_status": "AVAILABLE" if age_value is not None else "MISSING_RETAINED",
                "reference_stratum": metadata[analysis_id].get("reference_stratum"),
                "mloy_reference_eligible": bool(
                    row.get("sex_gate_status") == "PASS" and (reported_male or inferred_male)
                ),
                "held_out_specificity": sample_id in training_excluded_sample_ids,
                "sex_chromosome_anomaly_review": bool(anomaly_review),
                "collection": collected[sample_id],
            }
        )
    if not candidates:
        raise ValueError("no mLOY candidates or sex-chromosome anomaly reviews")

    by_assay: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in candidates:
        by_assay[row["assay"]].append(row)

    panel_receipt: dict[str, Any] = {
        "schema_version": 1,
        "gate_version": GATE_VERSION,
        "callability_rule": {
            "per_sample_interval_rule": f"bases_ge_{CALLABLE_DEPTH_THRESHOLD}/territory_bases >= {CALLABLE_BASE_FRACTION}",
            "cohort_rule": f"per-sample rule true in >= {CALLABLE_SAMPLE_FRACTION} of reference-eligible males within assay",
            "sample_metric": "median(log2(interval_mean_depth/autosomal_mean_depth))",
            "zero_aware_metric": "log2(territory-weighted mean interval depth/autosomal depth), including zero-depth intervals",
            "minimum_sample_panel_completeness": MIN_SAMPLE_PANEL_COMPLETENESS,
        },
        "normalization_rule": (
            f"{control_label} male median of the most specific control set with >=5 eligible "
            "controls: batch+reference stratum, batch, reference stratum, else assay "
            "(the stratum levels apply only when metadata supplies reference_stratum)"
        ),
        "control_group": control_label,
        "depth_metric": depth_metric,
        "score_metric": score_metric,
        "score_rule": (
            "mLRR_Y = median(log2(interval_mean_depth/autosomal_mean_depth)) - control median"
            if score_metric == "median_ratio" else
            "mLRR_Y = median over panel intervals of log2(depth_over_autosome / control interval "
            "median), with the sample left out of its own control set"
        ),
        "gc_correction": collection.get("gc_correction") if depth_metric == "gc" else None,
        "retention_policy": (
            "sex inference never removes reported males or sex-chromosome anomalies; "
            "no sample is labeled negative before calibration"
        ),
        "event_call_policy": "BLOCKED_PENDING_BAM_DOWNSAMPLING_CALIBRATION",
        "assays": {},
    }

    output_rows: list[dict[str, Any]] = []
    for assay, assay_rows in sorted(by_assay.items()):
        reference_rows = [
            row for row in assay_rows
            if row["mloy_reference_eligible"] and not row["held_out_specificity"]
        ]
        if not reference_rows:
            raise ValueError(f"{assay}: no reference-eligible males")
        interval_sets = [set(row["collection"]["eligible_y_intervals"]) for row in reference_rows]
        if any(values != interval_sets[0] for values in interval_sets[1:]):
            raise ValueError(f"{assay}: interval membership differs across samples")
        interval_ids = sorted(interval_sets[0])
        selected: list[str] = []
        interval_stats: list[dict[str, Any]] = []
        for interval_id in interval_ids:
            values = [row["collection"]["eligible_y_intervals"][interval_id] for row in reference_rows]
            territory = int(values[0]["territory_bases"])
            if any(int(value["territory_bases"]) != territory for value in values):
                raise ValueError(f"{assay}/{interval_id}: territory mismatch")
            callable_n = sum(
                value[f"bases_ge_{CALLABLE_DEPTH_THRESHOLD}"] / territory >= CALLABLE_BASE_FRACTION
                for value in values
            )
            fraction = callable_n / len(values)
            is_selected = fraction >= CALLABLE_SAMPLE_FRACTION
            if is_selected:
                selected.append(interval_id)
            interval_stats.append(
                {
                    "interval_id": interval_id,
                    "territory_bases": territory,
                    "callable_male_n": callable_n,
                    "reference_male_n": len(values),
                    "callable_male_fraction": fraction,
                    "selected": is_selected,
                }
            )
        if not selected:
            raise ValueError(f"{assay}: callable panel is empty")

        for row in assay_rows:
            values = row["collection"]["eligible_y_intervals"]
            positive = [
                math.log2(float(values[interval_id]["depth_over_autosome"]))
                for interval_id in selected
                if float(values[interval_id]["depth_over_autosome"]) > 0
            ]
            completeness = len(positive) / len(selected)
            weighted_numerator = sum(
                float(values[interval_id]["depth_over_autosome"])
                * int(values[interval_id]["territory_bases"])
                for interval_id in selected
            )
            weighted_denominator = sum(
                int(values[interval_id]["territory_bases"]) for interval_id in selected
            )
            weighted_ratio = weighted_numerator / weighted_denominator
            row["callable_interval_n"] = len(selected)
            row["positive_interval_n"] = len(positive)
            row["zero_interval_n"] = len(selected) - len(positive)
            row["panel_completeness"] = completeness
            row["raw_y_auto_weighted_ratio"] = weighted_ratio
            row["raw_log2_y_auto_weighted"] = math.log2(weighted_ratio) if weighted_ratio > 0 else None
            if not positive:
                row["raw_log2_y_auto_median"] = None
                row["metric_status"] = "Y_SIGNAL_ABSENT_MANUAL_REVIEW"
            elif completeness < MIN_SAMPLE_PANEL_COMPLETENESS:
                row["raw_log2_y_auto_median"] = statistics.median(positive)
                row["metric_status"] = "INCOMPLETE_Y_PANEL_HOLD"
            else:
                row["raw_log2_y_auto_median"] = statistics.median(positive)
                row["metric_status"] = "ELIGIBLE"

        eligible = [row for row in assay_rows if row["metric_status"] == "ELIGIBLE"]
        controls = [
            row for row in eligible
            if (
                row["study_group"] in control_groups
                and row["mloy_reference_eligible"]
                and not row["held_out_specificity"]
            )
        ]
        if len(controls) < 6:
            raise ValueError(f"{assay}: insufficient eligible {control_label} male controls")
        control_sets = build_control_sets(controls, assay)
        references = IntervalReferences("eligible_y_intervals", selected)

        for row in assay_rows:
            row["normalization_scope"] = None
            row["normalization_control_n"] = None
            row["normalization_median"] = None
            row["mLRR_Y"] = None
            row["mLRR_Y_zero_aware"] = None
            row["f_LOY_interpretive"] = None
            row["f_LOY_zero_aware_interpretive"] = None
            row["target_residual_median"] = None
            row["target_residual_q10"] = None
            row["target_residual_q90"] = None
            row["fraction_targets_below_control"] = None
            row["fraction_targets_below_minus_0_05"] = None
            row["fraction_targets_below_minus_0_10"] = None
            if row["metric_status"] != "ELIGIBLE":
                continue
            key = choose_control_set(control_sets, row, assay)
            interval_controls = control_sets[key]
            baseline = _median(
                [value["raw_log2_y_auto_median"] for value in interval_controls],
                f"{assay}/{key} controls",
            )
            scope = f"{key[0]}_{control_label.upper()}_MALES"
            control_n = len(interval_controls)
            row["normalization_leave_one_out"] = False
            if score_metric == "interval_residual":
                interval_refs = references.medians(key, interval_controls, row["sample_id"])
                mlrr = interval_residual_score(row["collection"]["eligible_y_intervals"], interval_refs)
                row["normalization_leave_one_out"] = any(
                    value["sample_id"] == row["sample_id"] for value in interval_controls
                )
            else:
                mlrr = row["raw_log2_y_auto_median"] - baseline
            weighted_baseline = _median(
                [value["raw_log2_y_auto_weighted"] for value in interval_controls],
                f"{assay}/{key} zero-aware controls",
            )
            zero_aware_mlrr = row["raw_log2_y_auto_weighted"] - weighted_baseline
            row["normalization_scope"] = scope
            row["normalization_set"] = key[1]
            row["normalization_control_n"] = control_n
            row["normalization_median"] = baseline
            row["mLRR_Y"] = mlrr
            row["mLRR_Y_zero_aware"] = zero_aware_mlrr
            row["f_LOY_interpretive"] = max(0.0, min(1.0, 1.0 - 2.0**mlrr))
            row["f_LOY_zero_aware_interpretive"] = max(
                0.0, min(1.0, 1.0 - 2.0**zero_aware_mlrr)
            )

            residuals: list[float] = []
            below_zero = below_005 = below_010 = 0
            values = row["collection"]["eligible_y_intervals"]
            for interval_id in selected:
                observed = float(values[interval_id]["depth_over_autosome"])
                reference = _median(
                    [
                        float(value["collection"]["eligible_y_intervals"][interval_id]["depth_over_autosome"])
                        for value in interval_controls
                    ],
                    f"{assay}/{interval_id} controls",
                )
                if reference <= 0:
                    raise ValueError(f"{assay}/{interval_id}: non-positive control median")
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
            ordered = sorted(residuals)
            row["target_residual_median"] = statistics.median(residuals) if residuals else None
            row["target_residual_q10"] = ordered[int(0.10 * (len(ordered) - 1))] if ordered else None
            row["target_residual_q90"] = ordered[int(0.90 * (len(ordered) - 1))] if ordered else None
            row["fraction_targets_below_control"] = below_zero / len(selected)
            row["fraction_targets_below_minus_0_05"] = below_005 / len(selected)
            row["fraction_targets_below_minus_0_10"] = below_010 / len(selected)

        age_controls = [row for row in controls if row["age_at_collection"] is not None]
        if len(age_controls) < 6:
            raise ValueError(f"{assay}: insufficient age-annotated {control_label} male controls")
        young_n = math.ceil(len(age_controls) / 3)
        young = sorted(age_controls, key=lambda row: (row["age_at_collection"], row["sample_id"]))[:young_n]
        null_values = [row["mLRR_Y"] for row in young]
        null_median = _median(null_values, f"{assay} young null")
        null_mad = _median([abs(value - null_median) for value in null_values], f"{assay} young MAD")
        robust_sigma = 1.4826 * null_mad
        if robust_sigma <= 0:
            raise ValueError(f"{assay}: degenerate young null")
        weak_threshold = null_median - WEAK_SIGNAL_SIGMA * robust_sigma
        extreme_threshold = null_median - EXTREME_SIGNAL_SIGMA * robust_sigma

        for row in assay_rows:
            row["young_null_median"] = null_median
            row["young_null_robust_sigma"] = robust_sigma
            row["weak_signal_threshold_2sigma"] = weak_threshold
            row["extreme_signal_threshold_3sigma"] = extreme_threshold
            row["screen_tier_LOY"] = "NOT_EVALUABLE_RETAINED"
            row["event_call_LOY"] = "NOT_CALIBRATED"
            if row["metric_status"] == "Y_SIGNAL_ABSENT_MANUAL_REVIEW":
                tier = "Y_SIGNAL_ABSENT_MANUAL_REVIEW"
            elif row["metric_status"] == "INCOMPLETE_Y_PANEL_HOLD":
                tier = "TECHNICAL_NOT_EVALUABLE_RETAINED"
            elif row["sex_chromosome_anomaly_review"]:
                tier = "SEX_CHROMOSOME_ANOMALY_REVIEW"
            elif row["metric_status"] == "ELIGIBLE":
                if (
                    row["mLRR_Y"] < extreme_threshold
                    and row["f_LOY_interpretive"] >= EXTREME_MIN_INTERPRETIVE_FRACTION
                ):
                    tier = "EXTREME_CANDIDATE_UNCALIBRATED"
                elif row["f_LOY_interpretive"] >= BORDERLINE_MIN_INTERPRETIVE_FRACTION:
                    tier = "BORDERLINE_CANDIDATE_UNCALIBRATED"
                elif row["mLRR_Y"] < weak_threshold:
                    tier = "WEAK_SIGNAL_RETAINED"
                else:
                    tier = "WITHIN_REFERENCE_UNCALIBRATED"
            else:
                tier = "NOT_EVALUABLE_RETAINED"
            row["screen_tier_LOY"] = tier
            row["screen_flag_LOY"] = tier

        panel_receipt["assays"][assay] = {
            "analysis_candidate_n": len(assay_rows),
            "reference_male_n": len(reference_rows),
            "eligible_control_male_n": len(controls),
            "selected_interval_n": len(selected),
            "selected_territory_bases": sum(
                next(item["territory_bases"] for item in interval_stats if item["interval_id"] == interval_id)
                for interval_id in selected
            ),
            "young_null_n": len(young),
            "young_null_max_age": max(row["age_at_collection"] for row in young),
            "young_null_median": null_median,
            "young_null_mad": null_mad,
            "young_null_robust_sigma": robust_sigma,
            "weak_signal_threshold_2sigma": weak_threshold,
            "extreme_signal_threshold_3sigma": extreme_threshold,
            "intervals": interval_stats,
        }
        output_rows.extend(assay_rows)

    clean_rows: list[dict[str, Any]] = []
    for row in sorted(output_rows, key=lambda value: value["sample_id"]):
        cleaned = {key: value for key, value in row.items() if key != "collection"}
        cleaned["mloy_gate_version"] = GATE_VERSION
        clean_rows.append(cleaned)

    statuses = Counter(row["metric_status"] for row in clean_rows)
    screens = Counter(row["screen_tier_LOY"] for row in clean_rows)
    retained_candidates = [
        row for row in clean_rows
        if row["screen_tier_LOY"] in {
            "EXTREME_CANDIDATE_UNCALIBRATED",
            "BORDERLINE_CANDIDATE_UNCALIBRATED",
            "WEAK_SIGNAL_RETAINED",
        }
    ]
    summary = {
        "status": "MLOY_GATE_METRICS_COMPLETE_EVENT_CALLS_BLOCKED",
        "gate_version": GATE_VERSION,
        "mloy_analysis_samples": len(clean_rows),
        "reference_male_samples": sum(row["mloy_reference_eligible"] for row in clean_rows),
        "sex_chromosome_anomaly_review_samples": sum(
            row["sex_chromosome_anomaly_review"] for row in clean_rows
        ),
        "metric_status_counts": dict(sorted(statuses.items())),
        "screen_tier_counts": dict(sorted(screens.items())),
        "retained_candidate_by_assay_and_group": {
            f"{assay}|{group}": count
            for (assay, group), count in sorted(
                Counter((row["assay"], row["study_group"]) for row in retained_candidates).items()
            )
        },
        "retained_candidate_age_range": (
            [
                min(row["age_at_collection"] for row in aged_candidates),
                max(row["age_at_collection"] for row in aged_candidates),
            ]
            if (aged_candidates := [r for r in retained_candidates if r["age_at_collection"] is not None])
            else None
        ),
        "missing_age_sample_n": sum(row["age_at_collection"] is None for row in clean_rows),
        "control_group": control_label,
        "depth_metric": depth_metric,
        "score_metric": score_metric,
        "normalization_scope_counts": dict(sorted(
            Counter(row["normalization_scope"] for row in clean_rows if row["normalization_scope"]).items()
        )),
        "manual_review_sample_ids": sorted(
            row["sample_id"]
            for row in clean_rows
            if row["metric_status"] == "Y_SIGNAL_ABSENT_MANUAL_REVIEW"
        ),
        "incomplete_panel_sample_ids": sorted(
            row["sample_id"]
            for row in clean_rows
            if row["metric_status"] == "INCOMPLETE_Y_PANEL_HOLD"
        ),
        "study_group_counts": dict(sorted(Counter(row["study_group"] for row in clean_rows).items())),
        "assay_counts": dict(sorted(Counter(row["assay"] for row in clean_rows).items())),
        "event_calls_emitted": False,
        "screening_is_not_an_event_call": True,
        "negative_calls_emitted": False,
        "training_excluded_sample_n": len(training_excluded_sample_ids),
        "borderline_signals_retained": True,
        "calibration_required": "BAM-level Y non-PAR downsampling by assay and depth",
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
    summary.update(
        {
            "collection_sha256": observed,
            "sex_cohort_sha256": sha256_file(args.sex_cohort),
            metadata_key: sha256_file(metadata_path),
        }
    )
    panel.update(
        {
            "collection_sha256": observed,
            "sex_cohort_sha256": summary["sex_cohort_sha256"],
            metadata_key: summary[metadata_key],
        }
    )
    _write_json(args.output_dir / "mloy.metrics.PRIVATE.json", rows)
    _write_json(args.output_dir / "mloy.summary.json", summary)
    _write_json(args.output_dir / "mloy.callable_panel.json", panel)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
