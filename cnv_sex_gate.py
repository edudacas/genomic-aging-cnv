#!/usr/bin/env python3
"""Calibrate and apply an assay-stratified coverage sex gate.

The gate is deliberately based only on X_NONPAR/autosome coverage.  Y coverage
is retained for audit but is not used because the imported cohort summaries
were produced with a provisional Y-XTR mask.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable


KNOWN_SEXES = {"Female", "Male"}
GATE_VERSION = "x-auto-assay-support-v1"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _finite_positive(value: Any, field: str, sample_id: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{sample_id}: invalid {field}: {value!r}") from exc
    if not math.isfinite(number) or number <= 0:
        raise ValueError(f"{sample_id}: {field} must be finite and > 0")
    return number


def _thresholds(rows: list[dict[str, Any]]) -> dict[str, Any]:
    males = [r for r in rows if r["sex_reported"] == "Male"]
    females = [r for r in rows if r["sex_reported"] == "Female"]
    if len(males) < 2 or len(females) < 2:
        raise ValueError("each assay requires at least two known males and females")

    male_upper = max(r["_x"] for r in males)
    female_lower = min(r["_x"] for r in females)
    if male_upper >= female_lower:
        raise ValueError(
            "known-sex X/autosome support overlaps; assay-specific gate cannot be frozen"
        )

    midpoint = math.sqrt(male_upper * female_lower)
    loo = Counter()
    for held_out in males + females:
        training = [r for r in males + females if r is not held_out]
        train_male_upper = max(
            r["_x"] for r in training if r["sex_reported"] == "Male"
        )
        train_female_lower = min(
            r["_x"] for r in training if r["sex_reported"] == "Female"
        )
        if train_male_upper >= train_female_lower:
            raise ValueError("leave-one-out training support overlaps")
        train_midpoint = math.sqrt(train_male_upper * train_female_lower)
        predicted = "Male" if held_out["_x"] < train_midpoint else "Female"
        loo[f"{held_out['sex_reported']}->{predicted}"] += 1

    correct = loo["Male->Male"] + loo["Female->Female"]
    total = len(males) + len(females)
    return {
        "known_female_n": len(females),
        "known_male_n": len(males),
        "male_support_upper_inclusive": male_upper,
        "female_support_lower_inclusive": female_lower,
        "empirical_gap_width": female_lower - male_upper,
        "binary_midpoint_for_leave_one_out_only": midpoint,
        "leave_one_out_confusion": dict(sorted(loo.items())),
        "leave_one_out_correct": correct,
        "leave_one_out_total": total,
        "leave_one_out_accuracy": correct / total,
    }


def calibrate_and_apply(
    input_rows: Iterable[dict[str, Any]], mask_sha256: str
) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    rows = [copy.deepcopy(row) for row in input_rows]
    if not rows:
        raise ValueError("empty cohort")

    sample_ids = [str(r.get("sample_id", "")) for r in rows]
    if any(not value for value in sample_ids) or len(sample_ids) != len(set(sample_ids)):
        raise ValueError("sample_id must be present and unique")

    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        sample_id = row["sample_id"]
        assay = str(row.get("proposed_capture_family") or "")
        if not assay:
            raise ValueError(f"{sample_id}: proposed_capture_family is missing")
        if row.get("sex_reported") not in KNOWN_SEXES | {"Unknown"}:
            raise ValueError(f"{sample_id}: unsupported sex_reported")
        row["_x"] = _finite_positive(row.get("x_auto_ratio"), "x_auto_ratio", sample_id)
        grouped[assay].append(row)

    thresholds = {assay: _thresholds(group) for assay, group in sorted(grouped.items())}
    if any(t["leave_one_out_accuracy"] != 1.0 for t in thresholds.values()):
        raise ValueError("leave-one-out validation did not achieve 100% concordance")

    inference_counts = Counter()
    reported_concordance = Counter()
    assay_counts: dict[str, Counter[str]] = defaultdict(Counter)
    ambiguous_ids: list[str] = []

    for row in rows:
        assay = row["proposed_capture_family"]
        gate = thresholds[assay]
        x_ratio = row.pop("_x")
        if x_ratio <= gate["male_support_upper_inclusive"]:
            inferred = "Male"
        elif x_ratio >= gate["female_support_lower_inclusive"]:
            inferred = "Female"
        else:
            inferred = None

        row["sex_inferred"] = inferred
        row["sex_inference_gate_version"] = GATE_VERSION
        row["sex_inference_mask_sha256"] = mask_sha256
        row["sex_inference_feature"] = "mean_depth(X_NONPAR)/mean_depth(AUTOSOMES)"
        row["sex_inference_uses_y"] = False
        row["y_auto_ratio_interpretation"] = (
            "AUDIT_ONLY_NOT_USED_PROVISIONAL_Y_XTR_MASK"
        )
        row["sex_gate_hold_reasons"] = []

        if inferred is None:
            row["sex_inference_status"] = "AMBIGUOUS_EMPIRICAL_X_GAP"
            row["sex_gate_status"] = "HOLD"
            row["sex_gate_hold_reasons"].append("SEX_INFERENCE_AMBIGUOUS_X_AUTO")
            ambiguous_ids.append(row["sample_id"])
            inference_counts["Ambiguous"] += 1
            assay_counts[assay]["Ambiguous"] += 1
        else:
            row["sex_inference_status"] = "INFERRED_ASSAY_EMPIRICAL_SUPPORT"
            row["sex_gate_status"] = "PASS"
            inference_counts[inferred] += 1
            assay_counts[assay][inferred] += 1

        reported = row["sex_reported"]
        if reported in KNOWN_SEXES:
            if inferred == reported:
                row["sex_reported_inferred_concordance"] = "CONCORDANT"
                reported_concordance["CONCORDANT"] += 1
            else:
                row["sex_reported_inferred_concordance"] = "DISCORDANT"
                row["sex_gate_status"] = "HOLD"
                row["sex_gate_hold_reasons"].append("SEX_REPORTED_INFERRED_DISCORDANT")
                reported_concordance["DISCORDANT"] += 1
        else:
            row["sex_reported_inferred_concordance"] = "NOT_APPLICABLE_UNREPORTED"
            reported_concordance["NOT_APPLICABLE_UNREPORTED"] += 1

    unknown_rows = [r for r in rows if r["sex_reported"] == "Unknown"]
    unknown_counts = Counter(
        r["sex_inferred"] if r["sex_inferred"] is not None else "Ambiguous"
        for r in unknown_rows
    )
    known_rows = [r for r in rows if r["sex_reported"] in KNOWN_SEXES]
    known_discordant = [
        r["sample_id"]
        for r in known_rows
        if r["sex_reported_inferred_concordance"] != "CONCORDANT"
    ]

    summary = {
        "status": "SEX_GATE_PASS_WITH_HOLDS" if ambiguous_ids else "SEX_GATE_PASS",
        "gate_version": GATE_VERSION,
        "samples": len(rows),
        "known_sex_samples": len(known_rows),
        "unknown_reported_samples": len(unknown_rows),
        "known_reported_inferred_discordant": known_discordant,
        "unknown_inference_counts": dict(sorted(unknown_counts.items())),
        "all_inference_counts": dict(sorted(inference_counts.items())),
        "reported_concordance_counts": dict(sorted(reported_concordance.items())),
        "ambiguous_sample_ids": sorted(ambiguous_ids),
        "assay_inference_counts": {
            assay: dict(sorted(counts.items())) for assay, counts in sorted(assay_counts.items())
        },
        "classification_basis": (
            "X_NONPAR/autosome coverage calibrated independently within capture family; "
            "values inside the empty empirical male-female support gap are held"
        ),
        "y_policy": (
            "Y/autosome retained for audit only and never used for classification because "
            "the imported summaries used a provisional Y-XTR mask"
        ),
        "mask_sha256": mask_sha256,
    }
    threshold_receipt = {
        "schema_version": 1,
        "gate_version": GATE_VERSION,
        "feature": "mean_depth(X_NONPAR)/mean_depth(AUTOSOMES)",
        "stratification": "proposed_capture_family",
        "decision_rule": (
            "Male at or below known-male maximum; Female at or above known-female "
            "minimum; otherwise Ambiguous/HOLD"
        ),
        "y_used": False,
        "mask_sha256": mask_sha256,
        "assays": thresholds,
    }
    return rows, summary, threshold_receipt


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, ensure_ascii=False)
        stream.write("\n")
    os.chmod(path, 0o600)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--mask", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected", type=int, default=1329)
    args = parser.parse_args()

    cohort = json.loads(args.cohort.read_text(encoding="utf-8"))
    if not isinstance(cohort, list) or len(cohort) != args.expected:
        raise SystemExit(f"expected {args.expected} cohort rows, observed {len(cohort)}")
    mask = json.loads(args.mask.read_text(encoding="utf-8"))
    if mask.get("sex_inference_status") != "FROZEN_X_NONPAR_GATE":
        raise SystemExit("mask is not frozen for X_NONPAR sex inference")

    mask_hash = sha256_file(args.mask)
    rows, summary, thresholds = calibrate_and_apply(cohort, mask_hash)
    summary["input_cohort_sha256"] = sha256_file(args.cohort)
    summary["mask_file"] = args.mask.name
    thresholds["input_cohort_sha256"] = summary["input_cohort_sha256"]

    _write_json(args.output_dir / "cohort.sex_gate.PRIVATE.json", rows)
    _write_json(args.output_dir / "sex_gate.summary.json", summary)
    _write_json(args.output_dir / "sex_gate.thresholds.json", thresholds)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
