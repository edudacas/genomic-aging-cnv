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


GATE_VERSION = "mloy-callable-y-v1"
CALLABLE_BASE_FRACTION = 0.80
CALLABLE_SAMPLE_FRACTION = 0.95
CALLABLE_DEPTH_THRESHOLD = 20
MIN_RUN_CONTROLS = 5
MIN_SAMPLE_PANEL_COMPLETENESS = 0.95


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


def build_gate(
    collection: dict[str, Any],
    sex_rows: list[dict[str, Any]],
    metadata: dict[str, dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    if collection.get("status") != "MLOY_INTERVAL_COLLECTION_COMPLETE":
        raise ValueError("mLOY interval collection is incomplete")
    collected = {row["sample_id"]: row for row in collection["samples"]}
    if len(collected) != collection.get("selected_samples"):
        raise ValueError("collection membership mismatch")
    sex = {row["sample_id"]: row for row in sex_rows}
    if len(sex) != len(sex_rows) or set(sex) != set(collected):
        raise ValueError("sex-gate/collection sample mismatch")

    males: list[dict[str, Any]] = []
    for sample_id, row in sorted(sex.items()):
        if row.get("sex_inferred") != "Male" or row.get("sex_gate_status") != "PASS":
            continue
        analysis_id = str(row.get("analysis_id"))
        if analysis_id not in metadata:
            raise ValueError(f"{sample_id}: workbook linkage missing")
        age = metadata[analysis_id].get("age_at_collection")
        if not isinstance(age, (int, float)) or not math.isfinite(float(age)) or age < 0:
            raise ValueError(f"{sample_id}: invalid age_at_collection")
        group = metadata[analysis_id].get("study_group")
        if not group:
            raise ValueError(f"{sample_id}: study_group missing")
        males.append(
            {
                "sample_id": sample_id,
                "analysis_id": analysis_id,
                "assay": row["proposed_capture_family"],
                "batch_id": str(row["batch_id"]),
                "batch_name": row.get("batch_name"),
                "sex_reported": row["sex_reported"],
                "sex_inferred": row["sex_inferred"],
                "study_group": group,
                "age_at_collection": float(age),
                "collection": collected[sample_id],
            }
        )
    if not males:
        raise ValueError("no eligible male samples")

    by_assay: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in males:
        by_assay[row["assay"]].append(row)

    panel_receipt: dict[str, Any] = {
        "schema_version": 1,
        "gate_version": GATE_VERSION,
        "callability_rule": {
            "per_sample_interval_rule": f"bases_ge_{CALLABLE_DEPTH_THRESHOLD}/territory_bases >= {CALLABLE_BASE_FRACTION}",
            "cohort_rule": f"per-sample rule true in >= {CALLABLE_SAMPLE_FRACTION} of eligible males within assay",
            "sample_metric": "median(log2(interval_mean_depth/autosomal_mean_depth))",
            "minimum_sample_panel_completeness": MIN_SAMPLE_PANEL_COMPLETENESS,
        },
        "normalization_rule": (
            "TP53-negative male batch median when >=5 eligible controls; "
            "otherwise TP53-negative male assay median"
        ),
        "event_call_policy": "BLOCKED_PENDING_BAM_DOWNSAMPLING_CALIBRATION",
        "assays": {},
    }

    output_rows: list[dict[str, Any]] = []
    for assay, assay_rows in sorted(by_assay.items()):
        interval_sets = [set(row["collection"]["eligible_y_intervals"]) for row in assay_rows]
        if any(values != interval_sets[0] for values in interval_sets[1:]):
            raise ValueError(f"{assay}: interval membership differs across samples")
        interval_ids = sorted(interval_sets[0])
        selected: list[str] = []
        interval_stats: list[dict[str, Any]] = []
        for interval_id in interval_ids:
            values = [row["collection"]["eligible_y_intervals"][interval_id] for row in assay_rows]
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
                    "eligible_male_n": len(values),
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
            row["callable_interval_n"] = len(selected)
            row["positive_interval_n"] = len(positive)
            row["panel_completeness"] = completeness
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
        controls = [row for row in eligible if row["study_group"] == "TP53_negative"]
        if len(controls) < 6:
            raise ValueError(f"{assay}: insufficient eligible TP53-negative male controls")
        assay_control_median = _median(
            [row["raw_log2_y_auto_median"] for row in controls], f"{assay} controls"
        )
        controls_by_batch: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in controls:
            controls_by_batch[row["batch_id"]].append(row)

        for row in assay_rows:
            row["normalization_scope"] = None
            row["normalization_control_n"] = None
            row["normalization_median"] = None
            row["mLRR_Y"] = None
            row["f_LOY_interpretive"] = None
            if row["metric_status"] != "ELIGIBLE":
                continue
            batch_controls = controls_by_batch.get(row["batch_id"], [])
            if len(batch_controls) >= MIN_RUN_CONTROLS:
                baseline = _median(
                    [value["raw_log2_y_auto_median"] for value in batch_controls],
                    f"{assay}/{row['batch_id']} controls",
                )
                scope = "BATCH_TP53_NEGATIVE_MALES"
                control_n = len(batch_controls)
            else:
                baseline = assay_control_median
                scope = "ASSAY_TP53_NEGATIVE_MALES"
                control_n = len(controls)
            mlrr = row["raw_log2_y_auto_median"] - baseline
            row["normalization_scope"] = scope
            row["normalization_control_n"] = control_n
            row["normalization_median"] = baseline
            row["mLRR_Y"] = mlrr
            row["f_LOY_interpretive"] = max(0.0, min(1.0, 1.0 - 2.0**mlrr))

        young_n = math.ceil(len(controls) / 3)
        young = sorted(controls, key=lambda row: (row["age_at_collection"], row["sample_id"]))[:young_n]
        null_values = [row["mLRR_Y"] for row in young]
        null_median = _median(null_values, f"{assay} young null")
        null_mad = _median([abs(value - null_median) for value in null_values], f"{assay} young MAD")
        robust_sigma = 1.4826 * null_mad
        if robust_sigma <= 0:
            raise ValueError(f"{assay}: degenerate young null")
        screen_threshold = null_median - 3.0 * robust_sigma

        for row in assay_rows:
            row["young_null_median"] = null_median
            row["young_null_robust_sigma"] = robust_sigma
            row["screen_threshold_3sigma"] = screen_threshold
            row["screen_flag_LOY"] = "NOT_EVALUABLE"
            row["event_call_LOY"] = "NOT_CALIBRATED"
            if row["metric_status"] == "Y_SIGNAL_ABSENT_MANUAL_REVIEW":
                row["screen_flag_LOY"] = "Y_SIGNAL_ABSENT_MANUAL_REVIEW"
            elif row["metric_status"] == "ELIGIBLE":
                is_screen = (
                    row["mLRR_Y"] < screen_threshold
                    and row["f_LOY_interpretive"] >= 0.10
                )
                row["screen_flag_LOY"] = (
                    "SCREEN_POSITIVE_UNCALIBRATED" if is_screen else "SCREEN_NEGATIVE_UNCALIBRATED"
                )

        panel_receipt["assays"][assay] = {
            "eligible_male_n": len(assay_rows),
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
            "screen_threshold_3sigma": screen_threshold,
            "intervals": interval_stats,
        }
        output_rows.extend(assay_rows)

    clean_rows: list[dict[str, Any]] = []
    for row in sorted(output_rows, key=lambda value: value["sample_id"]):
        cleaned = {key: value for key, value in row.items() if key != "collection"}
        cleaned["mloy_gate_version"] = GATE_VERSION
        clean_rows.append(cleaned)

    statuses = Counter(row["metric_status"] for row in clean_rows)
    screens = Counter(row["screen_flag_LOY"] for row in clean_rows)
    screen_positive = [
        row for row in clean_rows if row["screen_flag_LOY"] == "SCREEN_POSITIVE_UNCALIBRATED"
    ]
    summary = {
        "status": "MLOY_GATE_METRICS_COMPLETE_EVENT_CALLS_BLOCKED",
        "gate_version": GATE_VERSION,
        "male_samples": len(clean_rows),
        "metric_status_counts": dict(sorted(statuses.items())),
        "screen_flag_counts": dict(sorted(screens.items())),
        "screen_positive_by_assay_and_group": {
            f"{assay}|{group}": count
            for (assay, group), count in sorted(
                Counter((row["assay"], row["study_group"]) for row in screen_positive).items()
            )
        },
        "screen_positive_age_range": (
            [
                min(row["age_at_collection"] for row in screen_positive),
                max(row["age_at_collection"] for row in screen_positive),
            ]
            if screen_positive else None
        ),
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
    summary.update(
        {
            "collection_sha256": observed,
            "sex_cohort_sha256": sha256_file(args.sex_cohort),
            "workbook_sha256": sha256_file(args.workbook),
        }
    )
    panel.update(
        {
            "collection_sha256": observed,
            "sex_cohort_sha256": summary["sex_cohort_sha256"],
            "workbook_sha256": summary["workbook_sha256"],
        }
    )
    _write_json(args.output_dir / "mloy.metrics.PRIVATE.json", rows)
    _write_json(args.output_dir / "mloy.summary.json", summary)
    _write_json(args.output_dir / "mloy.callable_panel.json", panel)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
