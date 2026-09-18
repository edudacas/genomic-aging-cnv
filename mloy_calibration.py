#!/usr/bin/env python3
"""Calibrate coverage-based mLOY screening without emitting event calls.

Coverage mode performs Monte Carlo binomial thinning of observed Y read-base
counts while holding each donor's autosomal denominator fixed. BAM mode thins
mapped Y alignments and preserves every non-Y alignment. Neither mode is a
substitute for biological or orthogonal validation.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
from typing import Any

from cnv_contract import validate_artifact
from cnv_execution import atomic_json, sha


CALIBRATION_VERSION = "mloy-coverage-thinning-v1"
DEFAULT_FRACTIONS = (0.0, 0.01, 0.025, 0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.40)
DEFAULT_SEEDS = (17, 29, 43, 71, 101)
DEPTH_STRATA = ("LOW", "MID", "HIGH")
SIGNAL_TIERS = {
    "EXTREME_CANDIDATE_UNCALIBRATED",
    "BORDERLINE_CANDIDATE_UNCALIBRATED",
    "WEAK_SIGNAL_RETAINED",
}
BORDERLINE_TIERS = {
    "EXTREME_CANDIDATE_UNCALIBRATED",
    "BORDERLINE_CANDIDATE_UNCALIBRATED",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _finite_fraction(value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("cell fraction must be numeric")
    value = float(value)
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError("cell fraction must be finite in [0,1]")
    return value


def y_keep_probability(cell_fraction: float) -> float:
    """Expected Y-depth retention when fraction f of male cells lose Y."""
    return 1.0 - _finite_fraction(cell_fraction)


def y_template_retained(read_group: str, query_name: str, seed: int, cell_fraction: float) -> bool:
    probability = y_keep_probability(cell_fraction)
    if type(seed) is not int or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    if not isinstance(read_group, str) or not read_group:
        raise ValueError("read group is required")
    if not isinstance(query_name, str) or not query_name:
        raise ValueError("query name is required")
    key = json.dumps([seed, read_group, query_name], separators=(",", ":")).encode()
    draw = int.from_bytes(hashlib.sha256(key).digest()[:8], "big")
    cutoff = Fraction(str(probability)) * (1 << 64)
    return draw * cutoff.denominator < cutoff.numerator


def thin_y_bam(
    bam_artifact: dict[str, str],
    expected_bam_sample_name: str,
    output: Path,
    *,
    seed: int,
    cell_fraction: float,
    y_contig: str = "chrY",
) -> dict[str, Any]:
    """Thin mapped Y records only and preserve all non-Y records."""
    import pysam

    fraction = _finite_fraction(cell_fraction)
    keep_probability = y_keep_probability(fraction)
    validate_artifact(bam_artifact, True)
    if not isinstance(expected_bam_sample_name, str) or not expected_bam_sample_name:
        raise ValueError("explicit BAM sample name required")
    if not isinstance(y_contig, str) or not y_contig:
        raise ValueError("Y contig is required")
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    atomic_json(output / "STARTED.json", {
        "stage": "MLOY_Y_DEPTH_THINNING",
        "bam": bam_artifact,
        "expected_bam_sample_name": expected_bam_sample_name,
        "seed": seed,
        "cell_fraction": fraction,
        "y_keep_probability": keep_probability,
        "y_contig": y_contig,
        "pysam_version": pysam.__version__,
        "implementation_sha256": sha(__file__),
    })
    destination = output / "mloy.depth.bam"
    total = y_total = y_kept = non_y_total = non_y_kept = 0
    try:
        with pysam.AlignmentFile(bam_artifact["path"], "rb") as source:
            header = source.header.to_dict()
            if header.get("HD", {}).get("SO") != "coordinate":
                raise ValueError("coordinate-sorted BAM declaration required")
            if y_contig not in {entry.get("SN") for entry in header.get("SQ", [])}:
                raise ValueError("declared Y contig absent from BAM header")
            groups = header.get("RG", [])
            if (
                not groups
                or any(not group.get("ID") or group.get("SM") != expected_bam_sample_name for group in groups)
                or len({group["ID"] for group in groups}) != len(groups)
            ):
                raise ValueError("BAM read-group/sample identity mismatch")
            known = {group["ID"] for group in groups}
            with pysam.AlignmentFile(str(destination), "wb", header=source.header) as target:
                for read in source.fetch(until_eof=True):
                    total += 1
                    if any(read.has_tag(tag) for tag in ("RX", "MI", "UR", "UB")):
                        raise ValueError("UMI-tagged input requires molecule-aware Y thinning")
                    read_group = read.get_tag("RG") if read.has_tag("RG") else None
                    if read_group not in known:
                        raise ValueError("read lacks a declared read group")
                    on_y = not read.is_unmapped and source.get_reference_name(read.reference_id) == y_contig
                    if on_y:
                        y_total += 1
                        if y_template_retained(read_group, read.query_name, seed, fraction):
                            target.write(read)
                            y_kept += 1
                    else:
                        target.write(read)
                        non_y_total += 1
                        non_y_kept += 1
        pysam.index(str(destination))
        validate_artifact(bam_artifact, True)
        outputs = (destination, Path(str(destination) + ".bai"))
        receipt = {
            "status": "MLOY_Y_DEPTH_THINNING_COMPLETE",
            "calibration_version": CALIBRATION_VERSION,
            "seed": seed,
            "cell_fraction": fraction,
            "expected_y_record_retention": keep_probability,
            "observed_y_record_retention": y_kept / y_total if y_total else None,
            "input_alignment_records": total,
            "y_alignment_records": y_total,
            "retained_y_alignment_records": y_kept,
            "non_y_alignment_records": non_y_total,
            "retained_non_y_alignment_records": non_y_kept,
            "non_y_records_removed": 0,
            "depth_component_simulated": True,
            "event_calls_authorized": False,
            "intended_use": "Y_DEPTH_QC_CALIBRATION_ONLY_NOT_VARIANT_CALLING",
            "cross_chromosome_pair_consistency_guaranteed": False,
            "outputs": {
                item.name: {"sha256": sha(item), "bytes": item.stat().st_size}
                for item in outputs
            },
            "started_sha256": sha(output / "STARTED.json"),
        }
        atomic_json(output / "COMPLETE.json", receipt)
        return receipt
    except Exception as error:
        atomic_json(output / "FAILURE.json", {"error": type(error).__name__, "message": str(error)})
        raise


def _evenly_spaced(rows: list[dict[str, Any]], count: int) -> list[dict[str, Any]]:
    if len(rows) < count:
        raise ValueError("insufficient donors in depth stratum")
    indexes = [min(len(rows) - 1, math.floor((index + 0.5) * len(rows) / count)) for index in range(count)]
    if len(set(indexes)) != count:
        raise ValueError("donor selection produced duplicate indexes")
    return [rows[index] for index in indexes]


def select_donors(
    collection_by_id: dict[str, dict[str, Any]],
    metrics: list[dict[str, Any]],
    *,
    per_depth_stratum: int,
) -> list[dict[str, Any]]:
    if type(per_depth_stratum) is not int or per_depth_stratum < 1:
        raise ValueError("per-depth-stratum donor count must be positive")
    by_assay: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in metrics:
        if (
            row.get("mloy_reference_eligible") is True
            and row.get("sex_chromosome_anomaly_review") is False
            and row.get("study_group") == "TP53_negative"
            and row.get("metric_status") == "ELIGIBLE"
            and row.get("screen_tier_LOY") == "WITHIN_REFERENCE_UNCALIBRATED"
        ):
            sample_id = row["sample_id"]
            sample = collection_by_id.get(sample_id)
            if sample is None:
                raise ValueError(f"collection missing calibration donor {sample_id}")
            depth = float(sample["autosomal"]["mean_depth"])
            if not math.isfinite(depth) or depth <= 0:
                raise ValueError(f"invalid autosomal depth for {sample_id}")
            by_assay[row["assay"]].append({**row, "autosomal_mean_depth": depth})
    selected: list[dict[str, Any]] = []
    for assay, rows in sorted(by_assay.items()):
        rows.sort(key=lambda row: (row["autosomal_mean_depth"], row["sample_id"]))
        if len(rows) < 3 * per_depth_stratum:
            raise ValueError(f"{assay}: insufficient within-reference male controls for depth strata")
        strata: dict[str, list[dict[str, Any]]] = {name: [] for name in DEPTH_STRATA}
        for index, row in enumerate(rows):
            stratum = DEPTH_STRATA[min(2, math.floor(3 * index / len(rows)))]
            strata[stratum].append(row)
        for stratum in DEPTH_STRATA:
            for row in _evenly_spaced(strata[stratum], per_depth_stratum):
                selected.append({**row, "depth_stratum": stratum})
    if not selected:
        raise ValueError("no eligible calibration donors")
    return selected


def _classify(row: dict[str, Any], mlrr: float, estimated_fraction: float) -> str:
    if mlrr < float(row["extreme_signal_threshold_3sigma"]) and estimated_fraction >= 0.10:
        return "EXTREME_CANDIDATE_UNCALIBRATED"
    if estimated_fraction >= 0.05:
        return "BORDERLINE_CANDIDATE_UNCALIBRATED"
    if mlrr < float(row["weak_signal_threshold_2sigma"]):
        return "WEAK_SIGNAL_RETAINED"
    return "WITHIN_REFERENCE_UNCALIBRATED"


def _wilson(successes: int, total: int, z: float = 1.96) -> list[float]:
    if total <= 0:
        raise ValueError("Wilson interval requires observations")
    proportion = successes / total
    denominator = 1 + z * z / total
    centre = (proportion + z * z / (2 * total)) / denominator
    margin = z * math.sqrt(proportion * (1 - proportion) / total + z * z / (4 * total * total)) / denominator
    return [max(0.0, centre - margin), min(1.0, centre + margin)]


def run_coverage_calibration(
    collection: dict[str, Any],
    metrics: list[dict[str, Any]],
    panel: dict[str, Any],
    *,
    fractions: tuple[float, ...] = DEFAULT_FRACTIONS,
    seeds: tuple[int, ...] = DEFAULT_SEEDS,
    per_depth_stratum: int = 3,
) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    import numpy as np

    if collection.get("status") != "MLOY_INTERVAL_COLLECTION_COMPLETE":
        raise ValueError("mLOY interval collection is incomplete")
    if not metrics or not panel.get("assays"):
        raise ValueError("frozen mLOY metrics and callable panel are required")
    fractions = tuple(_finite_fraction(value) for value in fractions)
    if not fractions or fractions[0] != 0 or tuple(sorted(set(fractions))) != fractions:
        raise ValueError("fractions must be unique, sorted and start at zero")
    if not seeds or any(type(seed) is not int or seed < 0 for seed in seeds) or len(set(seeds)) != len(seeds):
        raise ValueError("seeds must be unique nonnegative integers")
    collection_by_id = {row["sample_id"]: row for row in collection["samples"]}
    donors = select_donors(collection_by_id, metrics, per_depth_stratum=per_depth_stratum)
    replicates: list[dict[str, Any]] = []
    for donor in donors:
        sample = collection_by_id[donor["sample_id"]]
        assay_panel = panel["assays"].get(donor["assay"])
        if assay_panel is None:
            raise ValueError(f"callable panel missing assay {donor['assay']}")
        interval_ids = [item["interval_id"] for item in assay_panel["intervals"] if item["selected"]]
        intervals = sample["eligible_y_intervals"]
        if any(interval_id not in intervals for interval_id in interval_ids):
            raise ValueError(f"{donor['sample_id']}: callable interval missing")
        depth_sums = np.asarray([int(intervals[key]["depth_sum"]) for key in interval_ids], dtype=np.int64)
        territories = np.asarray([int(intervals[key]["territory_bases"]) for key in interval_ids], dtype=np.float64)
        autosomal_depth = float(donor["autosomal_mean_depth"])
        baseline = float(donor["normalization_median"])
        weighted_baseline = float(donor["raw_log2_y_auto_weighted"]) - float(donor["mLRR_Y_zero_aware"])
        for fraction in fractions:
            probability = y_keep_probability(fraction)
            for seed in seeds:
                key = json.dumps([donor["sample_id"], fraction, seed], separators=(",", ":")).encode()
                derived_seed = int.from_bytes(hashlib.sha256(key).digest()[:8], "big")
                rng = np.random.default_rng(derived_seed)
                simulated_sums = rng.binomial(depth_sums, probability)
                ratios = simulated_sums / territories / autosomal_depth
                positive = ratios[ratios > 0]
                completeness = len(positive) / len(ratios)
                weighted_ratio = float(simulated_sums.sum() / territories.sum() / autosomal_depth)
                if completeness < 0.95 or weighted_ratio <= 0:
                    tier = "TECHNICAL_NOT_EVALUABLE_RETAINED"
                    raw_median = mlrr = estimated_fraction = None
                    zero_mlrr = zero_fraction = None
                else:
                    raw_median = float(np.median(np.log2(positive)))
                    mlrr = raw_median - baseline
                    estimated_fraction = max(0.0, min(1.0, 1.0 - 2.0**mlrr))
                    zero_mlrr = math.log2(weighted_ratio) - weighted_baseline
                    zero_fraction = max(0.0, min(1.0, 1.0 - 2.0**zero_mlrr))
                    tier = _classify(donor, mlrr, estimated_fraction)
                replicates.append({
                    "sample_id": donor["sample_id"],
                    "analysis_id": donor.get("analysis_id"),
                    "assay": donor["assay"],
                    "depth_stratum": donor["depth_stratum"],
                    "autosomal_mean_depth": autosomal_depth,
                    "cell_fraction": fraction,
                    "seed": seed,
                    "derived_seed": derived_seed,
                    "y_keep_probability": probability,
                    "expected_log2_shift": math.log2(probability) if probability > 0 else None,
                    "simulated_raw_log2_y_auto_median": raw_median,
                    "simulated_mLRR_Y": mlrr,
                    "estimated_f_LOY": estimated_fraction,
                    "simulated_mLRR_Y_zero_aware": zero_mlrr,
                    "estimated_f_LOY_zero_aware": zero_fraction,
                    "screen_tier_LOY": tier,
                    "retained_signal": tier in SIGNAL_TIERS,
                    "borderline_or_extreme": tier in BORDERLINE_TIERS,
                    "event_call_LOY": "NOT_CALIBRATED",
                })
    grouped: dict[tuple[str, str, float], list[dict[str, Any]]] = defaultdict(list)
    for row in replicates:
        grouped[(row["assay"], row["depth_stratum"], row["cell_fraction"])].append(row)
    strata_results: list[dict[str, Any]] = []
    for (assay, depth_stratum, fraction), rows in sorted(grouped.items()):
        retained = sum(row["retained_signal"] for row in rows)
        borderline = sum(row["borderline_or_extreme"] for row in rows)
        estimates = [row["estimated_f_LOY"] for row in rows if row["estimated_f_LOY"] is not None]
        strata_results.append({
            "assay": assay,
            "depth_stratum": depth_stratum,
            "cell_fraction": fraction,
            "replicates": len(rows),
            "retained_signal_n": retained,
            "retained_signal_rate": retained / len(rows),
            "borderline_or_extreme_n": borderline,
            "borderline_or_extreme_rate": borderline / len(rows),
            "borderline_or_extreme_wilson_95": _wilson(borderline, len(rows)),
            "median_estimated_f_LOY": statistics.median(estimates) if estimates else None,
        })
    assays = sorted({row["assay"] for row in replicates})
    lod_by_assay: dict[str, float | None] = {}
    false_positive_by_assay: dict[str, dict[str, Any]] = {}
    for assay in assays:
        zero_by_donor: dict[str, bool] = {}
        for row in replicates:
            if row["assay"] == assay and row["cell_fraction"] == 0:
                zero_by_donor[row["sample_id"]] = zero_by_donor.get(row["sample_id"], False) or row["retained_signal"]
        false_n = sum(zero_by_donor.values())
        false_positive_by_assay[assay] = {
            "donors": len(zero_by_donor),
            "retained_signal_donors": false_n,
            "retained_signal_rate": false_n / len(zero_by_donor),
        }
        lod = None
        for fraction in fractions[1:]:
            strata = [row for row in strata_results if row["assay"] == assay and row["cell_fraction"] == fraction]
            pooled = [row for row in replicates if row["assay"] == assay and row["cell_fraction"] == fraction]
            pooled_rate = sum(row["borderline_or_extreme"] for row in pooled) / len(pooled)
            if pooled_rate >= 0.95 and all(row["borderline_or_extreme_rate"] >= 0.90 for row in strata):
                lod = fraction
                break
        lod_by_assay[assay] = lod
    donor_plan = {
        "status": "MLOY_CALIBRATION_DONORS_FROZEN",
        "calibration_version": CALIBRATION_VERSION,
        "selection": "TP53-negative reference males within uncalibrated reference; deterministic depth tertiles",
        "per_assay_per_depth_stratum": per_depth_stratum,
        "donors": [
            {
                "sample_id": row["sample_id"],
                "analysis_id": row.get("analysis_id"),
                "assay": row["assay"],
                "depth_stratum": row["depth_stratum"],
                "autosomal_mean_depth": row["autosomal_mean_depth"],
            }
            for row in donors
        ],
    }
    summary = {
        "status": "MLOY_COVERAGE_THINNING_CALIBRATION_COMPLETE_BAM_VALIDATION_REQUIRED",
        "calibration_version": CALIBRATION_VERSION,
        "simulation_unit": "interval depth_sum binomial thinning",
        "autosomal_denominator": "held fixed per donor",
        "cell_fraction_model": "Y depth retention = 1 - f",
        "fractions": list(fractions),
        "seeds": list(seeds),
        "donor_n": len(donors),
        "replicate_n": len(replicates),
        "strata": strata_results,
        "false_positive_by_assay": false_positive_by_assay,
        "provisional_borderline_lod_by_assay": lod_by_assay,
        "lod_rule": "lowest fraction with pooled borderline-or-extreme sensitivity >=0.95 and every depth stratum >=0.90",
        "event_calls_emitted": False,
        "negative_calls_emitted": False,
        "clinical_validation_complete": False,
        "bam_validation_required": True,
        "interpretation": "analytical coverage calibration only; true or orthogonally confirmed mLOY remains required for biological validation",
    }
    return replicates, summary, donor_plan


def _write_private_gzip(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", compresslevel=6) as stream:
        json.dump(payload, stream, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    os.chmod(path, 0o600)


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    os.chmod(path, 0o600)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="mode", required=True)
    coverage = subparsers.add_parser("coverage", help="Monte Carlo thinning of collected interval depth")
    coverage.add_argument("--collection", type=Path, required=True)
    coverage.add_argument("--metrics", type=Path, required=True)
    coverage.add_argument("--panel", type=Path, required=True)
    coverage.add_argument("--output-dir", type=Path, required=True)
    coverage.add_argument("--collection-sha256")
    coverage.add_argument("--fractions", default=",".join(map(str, DEFAULT_FRACTIONS)))
    coverage.add_argument("--seeds", default=",".join(map(str, DEFAULT_SEEDS)))
    coverage.add_argument("--per-depth-stratum", type=int, default=3)
    bam = subparsers.add_parser("thin-bam", help="Create one depth-only Y-thinned BAM")
    bam.add_argument("--bam", type=Path, required=True)
    bam.add_argument("--bam-sha256", required=True)
    bam.add_argument("--sample-name", required=True)
    bam.add_argument("--output-dir", type=Path, required=True)
    bam.add_argument("--cell-fraction", type=float, required=True)
    bam.add_argument("--seed", type=int, required=True)
    bam.add_argument("--y-contig", default="chrY")
    args = parser.parse_args()
    if args.mode == "thin-bam":
        receipt = thin_y_bam(
            {"path": str(args.bam.resolve()), "sha256": args.bam_sha256.lower()},
            args.sample_name,
            args.output_dir,
            seed=args.seed,
            cell_fraction=args.cell_fraction,
            y_contig=args.y_contig,
        )
        print(json.dumps(receipt, indent=2, sort_keys=True))
        return 0
    collection_hash = sha256_file(args.collection)
    if args.collection_sha256 and collection_hash != args.collection_sha256.lower():
        raise SystemExit("collection SHA-256 mismatch")
    with gzip.open(args.collection, "rt", encoding="utf-8") as stream:
        collection = json.load(stream)
    metrics = json.loads(args.metrics.read_text(encoding="utf-8"))
    panel = json.loads(args.panel.read_text(encoding="utf-8"))
    fractions = tuple(float(value) for value in args.fractions.split(","))
    seeds = tuple(int(value) for value in args.seeds.split(","))
    replicates, summary, donor_plan = run_coverage_calibration(
        collection,
        metrics,
        panel,
        fractions=fractions,
        seeds=seeds,
        per_depth_stratum=args.per_depth_stratum,
    )
    provenance = {
        "collection_sha256": collection_hash,
        "metrics_sha256": sha256_file(args.metrics),
        "panel_sha256": sha256_file(args.panel),
        "implementation_sha256": sha256_file(Path(__file__)),
    }
    summary.update(provenance)
    donor_plan.update(provenance)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    _write_private_gzip(args.output_dir / "mloy.calibration.replicates.PRIVATE.json.gz", replicates)
    _write_json(args.output_dir / "mloy.calibration.summary.json", summary)
    _write_json(args.output_dir / "mloy.calibration.plan.PRIVATE.json", donor_plan)
    receipt = {
        "status": summary["status"],
        "outputs": {
            path.name: {"sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in sorted(args.output_dir.iterdir())
        },
        **provenance,
    }
    _write_json(args.output_dir / "COMPLETE.json", receipt)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
