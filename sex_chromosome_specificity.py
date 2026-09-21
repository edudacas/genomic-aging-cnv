#!/usr/bin/env python3
"""Cross-validated reference-control behavior for mLOX/mLOY screens.

Every TP53-negative control is scored exactly once after excluding its fold
from callable-panel construction, normalization, and young-null thresholds.
The output is a screening-specificity receipt; it never emits event or negative
calls and it does not use orthogonal truth cases for fitting.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any

import mlox_gate
import mloy_gate


AXES = {
    "mlox": {
        "collection_status": "MLOX_INTERVAL_COLLECTION_COMPLETE",
        "target_sex": "Female",
        "tier_field": "screen_tier_LOX",
        "build_gate": mlox_gate.build_gate,
        "load_metadata": mlox_gate._load_workbook_metadata,
    },
    "mloy": {
        "collection_status": "MLOY_INTERVAL_COLLECTION_COMPLETE",
        "target_sex": "Male",
        "tier_field": "screen_tier_LOY",
        "build_gate": mloy_gate.build_gate,
        "load_metadata": mloy_gate._load_workbook_metadata,
    },
}
SIGNAL_TIERS = {
    "EXTREME_CANDIDATE_UNCALIBRATED",
    "BORDERLINE_CANDIDATE_UNCALIBRATED",
    "WEAK_SIGNAL_RETAINED",
}
BORDERLINE_TIERS = {
    "EXTREME_CANDIDATE_UNCALIBRATED",
    "BORDERLINE_CANDIDATE_UNCALIBRATED",
}
VERSION = "sex-chromosome-null-specificity-v1"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _wilson(successes: int, total: int, z: float = 1.96) -> list[float]:
    if total <= 0:
        raise ValueError("Wilson interval requires observations")
    proportion = successes / total
    denominator = 1 + z * z / total
    centre = (proportion + z * z / (2 * total)) / denominator
    margin = z * math.sqrt(
        proportion * (1 - proportion) / total + z * z / (4 * total * total)
    ) / denominator
    return [max(0.0, centre - margin), min(1.0, centre + margin)]


def _eligible_controls(
    axis: str,
    collection: dict[str, Any],
    sex_rows: list[dict[str, Any]],
    metadata: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    config = AXES[axis]
    collected = {row["sample_id"]: row for row in collection["samples"]}
    controls: list[dict[str, Any]] = []
    for row in sex_rows:
        sample_id = row["sample_id"]
        analysis_id = str(row.get("analysis_id"))
        info = metadata.get(analysis_id)
        if info is None or info.get("study_group") != "TP53_negative":
            continue
        target = config["target_sex"]
        if row.get("sex_gate_status") != "PASS" or not (
            row.get("sex_reported") == target or row.get("sex_inferred") == target
        ):
            continue
        sample = collected.get(sample_id)
        if sample is None:
            raise ValueError(f"collection missing null control {sample_id}")
        depth = float(sample["autosomal"]["mean_depth"])
        if not math.isfinite(depth) or depth <= 0:
            raise ValueError(f"invalid autosomal depth for {sample_id}")
        controls.append({
            "sample_id": sample_id,
            "analysis_id": analysis_id,
            "assay": row["proposed_capture_family"],
            "batch_id": str(row["batch_id"]),
            "autosomal_mean_depth": depth,
        })
    if not controls:
        raise ValueError("no eligible TP53-negative null controls")
    return controls


def _assign_folds(controls: list[dict[str, Any]], folds: int) -> list[dict[str, Any]]:
    if type(folds) is not int or folds < 2:
        raise ValueError("at least two folds are required")
    by_assay: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in controls:
        by_assay[row["assay"]].append(row)
    prepared: list[dict[str, Any]] = []
    for assay, rows in sorted(by_assay.items()):
        if len(rows) - math.ceil(len(rows) / folds) < 6:
            raise ValueError(f"{assay}: too few training controls after fold exclusion")
        ordered = sorted(rows, key=lambda row: (row["autosomal_mean_depth"], row["sample_id"]))
        for index, row in enumerate(ordered):
            stratum = ("LOW", "MID", "HIGH")[min(2, math.floor(3 * index / len(ordered)))]
            prepared.append({**row, "depth_stratum": stratum})
    strata: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in prepared:
        strata[(row["assay"], row["batch_id"], row["depth_stratum"])].append(row)
    assigned: list[dict[str, Any]] = []
    for key, rows in sorted(strata.items()):
        rows = sorted(
            rows,
            key=lambda row: hashlib.sha256(row["sample_id"].encode()).hexdigest(),
        )
        offset = int.from_bytes(hashlib.sha256(json.dumps(key).encode()).digest()[:4], "big") % folds
        for index, row in enumerate(rows):
            assigned.append({**row, "fold": (offset + index) % folds})
    return sorted(assigned, key=lambda row: row["sample_id"])


def cross_validate_specificity(
    axis: str,
    collection: dict[str, Any],
    sex_rows: list[dict[str, Any]],
    metadata: dict[str, dict[str, Any]],
    *,
    folds: int = 5,
) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    if axis not in AXES:
        raise ValueError("axis must be mlox or mloy")
    config = AXES[axis]
    if collection.get("status") != config["collection_status"]:
        raise ValueError(f"{axis} interval collection is incomplete")
    controls = _assign_folds(
        _eligible_controls(axis, collection, sex_rows, metadata), folds
    )
    expected_ids = {row["sample_id"] for row in controls}
    scored: list[dict[str, Any]] = []
    fold_receipts: list[dict[str, Any]] = []
    for fold in range(folds):
        excluded = {row["sample_id"] for row in controls if row["fold"] == fold}
        if not excluded:
            fold_receipts.append({"fold": fold, "held_out_n": 0, "status": "EMPTY_FOLD"})
            continue
        metrics, gate_summary, panel = config["build_gate"](
            collection,
            sex_rows,
            metadata,
            training_excluded_sample_ids=excluded,
        )
        by_id = {row["sample_id"]: row for row in metrics}
        missing = excluded - set(by_id)
        if missing:
            raise ValueError(f"fold {fold}: held-out samples not scored: {sorted(missing)}")
        training_by_assay = Counter(
            row["assay"] for row in controls if row["sample_id"] not in excluded
        )
        for design in (row for row in controls if row["fold"] == fold):
            result = by_id[design["sample_id"]]
            tier = result[config["tier_field"]]
            scored.append({
                **design,
                "screen_tier": tier,
                "retained_signal": tier in SIGNAL_TIERS,
                "borderline_or_extreme": tier in BORDERLINE_TIERS,
                "metric_status": result["metric_status"],
                "normalization_scope": result.get("normalization_scope"),
                "training_control_n_assay": training_by_assay[design["assay"]],
                "event_call": "NOT_CALIBRATED",
                "negative_call": "NOT_EMITTED",
            })
        fold_receipts.append({
            "fold": fold,
            "held_out_n": len(excluded),
            "training_n": len(controls) - len(excluded),
            "gate_status": gate_summary["status"],
            "panel_gate_version": panel["gate_version"],
            "status": "COMPLETE",
        })
    if len(scored) != len(expected_ids) or {row["sample_id"] for row in scored} != expected_ids:
        raise ValueError("each null control must be scored exactly once")
    by_assay: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in scored:
        by_assay[row["assay"]].append(row)
    specificity: dict[str, Any] = {}
    for assay, rows in sorted(by_assay.items()):
        retained = sum(row["retained_signal"] for row in rows)
        borderline = sum(row["borderline_or_extreme"] for row in rows)
        specificity[assay] = {
            "held_out_null_n": len(rows),
            "retained_signal_n": retained,
            "retained_signal_rate": retained / len(rows),
            "retained_signal_rate_wilson_95": _wilson(retained, len(rows)),
            "borderline_or_extreme_n": borderline,
            "borderline_or_extreme_rate": borderline / len(rows),
            "borderline_or_extreme_rate_wilson_95": _wilson(borderline, len(rows)),
            "screen_tier_counts": dict(sorted(Counter(row["screen_tier"] for row in rows).items())),
            "false_positive_rate_interpretable": False,
        }
    plan = {
        "status": "NULL_SPECIFICITY_FOLDS_FROZEN",
        "version": VERSION,
        "axis": axis,
        "folds": folds,
        "stratification": "assay, batch, assay-specific autosomal-depth tertile",
        "assignment": "deterministic SHA-256 order with stratum-specific fold offset",
        "controls": controls,
    }
    summary = {
        "status": "NULL_SPECIFICITY_CROSS_VALIDATION_COMPLETE_EVENT_CALLS_BLOCKED",
        "version": VERSION,
        "axis": axis,
        "held_out_null_n": len(scored),
        "fold_receipts": fold_receipts,
        "held_out_reference_behavior_by_assay": specificity,
        "callable_panel_refit_without_held_out_fold": True,
        "normalization_refit_without_held_out_fold": True,
        "young_null_thresholds_refit_without_held_out_fold": True,
        "orthogonal_truth_used_for_fitting": False,
        "event_calls_emitted": False,
        "negative_calls_emitted": False,
        "interpretation": "cross-validated reference-control retained-signal rate; TP53-negative status is not orthogonal absence of age-related mosaicism, so this is not a false-positive rate",
    }
    return scored, summary, plan


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    os.chmod(path, 0o600)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--axis", choices=sorted(AXES), required=True)
    parser.add_argument("--collection", type=Path, required=True)
    parser.add_argument("--sex-cohort", type=Path, required=True)
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    with gzip.open(args.collection, "rt", encoding="utf-8") as stream:
        collection = json.load(stream)
    sex_rows = json.loads(args.sex_cohort.read_text(encoding="utf-8"))
    metadata = AXES[args.axis]["load_metadata"](args.workbook)
    rows, summary, plan = cross_validate_specificity(
        args.axis, collection, sex_rows, metadata, folds=args.folds
    )
    provenance = {
        "collection_sha256": sha256_file(args.collection),
        "sex_cohort_sha256": sha256_file(args.sex_cohort),
        "workbook_sha256": sha256_file(args.workbook),
        "implementation_sha256": sha256_file(Path(__file__)),
    }
    summary.update(provenance)
    plan.update(provenance)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    _write_json(args.output_dir / f"{args.axis}.specificity.rows.PRIVATE.json", rows)
    _write_json(args.output_dir / f"{args.axis}.specificity.summary.json", summary)
    _write_json(args.output_dir / f"{args.axis}.specificity.plan.PRIVATE.json", plan)
    outputs = {
        path.name: {"sha256": sha256_file(path), "bytes": path.stat().st_size}
        for path in sorted(args.output_dir.iterdir())
    }
    _write_json(args.output_dir / "COMPLETE.json", {
        "status": summary["status"],
        "outputs": outputs,
        **provenance,
    })
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
