#!/usr/bin/env python3
"""Target-only read-level mLOX/mLOY depth calibration.

The scanner reads only frozen callable X or Y intervals, keeps the original
autosomal denominator, and applies deterministic template thinning in memory.
It writes no pseudo-BAMs. Fraction zero must exactly reconstruct every frozen
interval depth_sum or the run fails closed as a filter-contract mismatch.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
from typing import Any

from cnv_contract import validate_artifact
from mlox_calibration import (
    BORDERLINE_TIERS as MLOX_BORDERLINE,
    SIGNAL_TIERS as MLOX_SIGNAL,
    _classify as classify_mlox,
    x_keep_probability,
    x_template_retained,
)
from mloy_calibration import (
    BORDERLINE_TIERS as MLOY_BORDERLINE,
    SIGNAL_TIERS as MLOY_SIGNAL,
    _classify as classify_mloy,
    y_keep_probability,
    y_template_retained,
)


VERSION = "sex-chromosome-target-read-calibration-v1"
DEFAULT_SEEDS = (17, 29, 43)
DEFAULT_FRACTIONS = {
    "mlox": {
        "MGI_IDT_EXOME": (0.0, 0.10, 0.15, 0.20),
        "TWIST_WES_OS_T7": (0.0, 0.025, 0.05, 0.075),
    },
    "mloy": {
        "MGI_IDT_EXOME": (0.0, 0.05, 0.075, 0.10),
        "TWIST_WES_OS_T7": (0.0, 0.075, 0.10, 0.15),
    },
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _parse_interval(interval_id: str, expected_contig: str) -> tuple[int, int]:
    chrom, coordinates = interval_id.split(":", 1)
    start, end = (int(value) for value in coordinates.split("-", 1))
    if chrom != expected_contig or not 0 <= start < end:
        raise ValueError(f"invalid frozen interval: {interval_id}")
    return start, end


def _overlap_bases(read: Any, start: int, end: int) -> int:
    return sum(max(0, min(right, end) - max(left, start)) for left, right in read.get_blocks())


def _validate_fractions(fractions: tuple[float, ...]) -> tuple[float, ...]:
    values = tuple(float(value) for value in fractions)
    if (
        not values
        or values[0] != 0
        or tuple(sorted(set(values))) != values
        or any(not math.isfinite(value) or not 0 <= value <= 1 for value in values)
    ):
        raise ValueError("fractions must be unique, sorted, finite and start at zero")
    return values


def scan_sample(
    axis: str,
    binding: dict[str, Any],
    donor: dict[str, Any],
    sample: dict[str, Any],
    metric: dict[str, Any],
    assay_panel: dict[str, Any],
    *,
    fractions: tuple[float, ...],
    seeds: tuple[int, ...] = DEFAULT_SEEDS,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    import pysam

    if axis not in ("mlox", "mloy"):
        raise ValueError("axis must be mlox or mloy")
    fractions = _validate_fractions(fractions)
    if not seeds or len(set(seeds)) != len(seeds) or any(type(seed) is not int or seed < 0 for seed in seeds):
        raise ValueError("seeds must be unique nonnegative integers")
    if donor["sample_id"] != sample["sample_id"] or metric["sample_id"] != sample["sample_id"]:
        raise ValueError("donor, collection and metric sample mismatch")
    if donor["assay"] != metric["assay"] or binding["assay"] != donor["assay"]:
        raise ValueError("assay mismatch")
    validate_artifact(binding["bam"], True)
    validate_artifact(binding["bam_index"], True)
    bam_path = Path(binding["bam"]["path"])
    expected_index = Path(str(bam_path) + ".bai")
    if Path(binding["bam_index"]["path"]).resolve() != expected_index.resolve():
        raise ValueError("BAM index must be the adjacent .bam.bai artifact")

    contig = "chrX" if axis == "mlox" else "chrY"
    interval_field = "eligible_x_intervals" if axis == "mlox" else "eligible_y_intervals"
    intervals = sample[interval_field]
    interval_ids = [row["interval_id"] for row in assay_panel["intervals"] if row["selected"]]
    if not interval_ids or any(interval_id not in intervals for interval_id in interval_ids):
        raise ValueError("frozen callable interval missing from sample")
    combinations = [(fraction, seed) for fraction in fractions for seed in seeds]
    simulated = {key: {interval_id: 0 for interval_id in interval_ids} for key in combinations}
    observed = {interval_id: 0 for interval_id in interval_ids}
    records_seen = filtered_records = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as stream:
        header = stream.header.to_dict()
        if header.get("HD", {}).get("SO") != "coordinate":
            raise ValueError("coordinate-sorted BAM declaration required")
        if contig not in {row.get("SN") for row in header.get("SQ", [])}:
            raise ValueError(f"{contig} absent from BAM header")
        groups = header.get("RG", [])
        expected_name = binding.get("bam_sample_name")
        if (
            not expected_name
            or not groups
            or any(not row.get("ID") or row.get("SM") != expected_name for row in groups)
        ):
            raise ValueError("BAM read-group/sample identity mismatch")
        known_groups = {row["ID"] for row in groups}
        for interval_id in interval_ids:
            start, end = _parse_interval(interval_id, contig)
            for read in stream.fetch(contig, start, end):
                records_seen += 1
                if read.is_unmapped or read.is_secondary or read.is_qcfail or read.is_duplicate:
                    filtered_records += 1
                    continue
                if any(read.has_tag(tag) for tag in ("RX", "MI", "UR", "UB")):
                    raise ValueError("UMI-tagged input requires molecule-aware thinning")
                read_group = read.get_tag("RG") if read.has_tag("RG") else None
                if read_group not in known_groups:
                    raise ValueError("read lacks a declared read group")
                overlap = _overlap_bases(read, start, end)
                if not overlap:
                    continue
                observed[interval_id] += overlap
                for fraction, seed in combinations:
                    retained = (
                        x_template_retained(read_group, read.query_name, seed, fraction)
                        if axis == "mlox"
                        else y_template_retained(read_group, read.query_name, seed, fraction)
                    )
                    if retained:
                        simulated[(fraction, seed)][interval_id] += overlap

    mismatches = [
        {
            "interval_id": interval_id,
            "frozen_depth_sum": int(intervals[interval_id]["depth_sum"]),
            "read_level_depth_sum": observed[interval_id],
        }
        for interval_id in interval_ids
        if observed[interval_id] != int(intervals[interval_id]["depth_sum"])
    ]
    if mismatches:
        preview = mismatches[:10]
        raise ValueError(
            f"fraction-zero filter contract mismatch in {len(mismatches)}/{len(interval_ids)} intervals: {preview}"
        )
    autosomal_depth = float(sample["autosomal"]["mean_depth"])
    results: list[dict[str, Any]] = []
    for fraction, seed in combinations:
        sums = simulated[(fraction, seed)]
        ratios = [
            sums[interval_id] / int(intervals[interval_id]["territory_bases"]) / autosomal_depth
            for interval_id in interval_ids
        ]
        positive = [value for value in ratios if value > 0]
        completeness = len(positive) / len(ratios)
        tier = "TECHNICAL_NOT_EVALUABLE_RETAINED"
        raw_median = mlrr = estimated = None
        zero_mlrr = zero_estimated = None
        if completeness >= 0.95:
            raw_median = statistics.median(math.log2(value) for value in positive)
            if axis == "mlox":
                female_baseline = float(metric["female_normalization_median"])
                male_baseline = float(metric["male_dosage_median"])
                female_linear, male_linear = 2.0**female_baseline, 2.0**male_baseline
                mlrr = raw_median - female_baseline
                estimated = max(0.0, min(1.0, (female_linear - 2.0**raw_median) / (female_linear - male_linear)))
                tier = classify_mlox(metric, mlrr, estimated)
            else:
                baseline = float(metric["normalization_median"])
                mlrr = raw_median - baseline
                estimated = max(0.0, min(1.0, 1.0 - 2.0**mlrr))
                territory_sum = sum(int(intervals[key]["territory_bases"]) for key in interval_ids)
                weighted_ratio = sum(sums.values()) / territory_sum / autosomal_depth
                frozen_weighted_baseline = float(metric["raw_log2_y_auto_weighted"]) - float(metric["mLRR_Y_zero_aware"])
                if weighted_ratio > 0:
                    zero_mlrr = math.log2(weighted_ratio) - frozen_weighted_baseline
                    zero_estimated = max(0.0, min(1.0, 1.0 - 2.0**zero_mlrr))
                tier = classify_mloy(metric, mlrr, estimated)
        signal = MLOX_SIGNAL if axis == "mlox" else MLOY_SIGNAL
        borderline = MLOX_BORDERLINE if axis == "mlox" else MLOY_BORDERLINE
        results.append({
            "axis": axis,
            "sample_id": sample["sample_id"],
            "analysis_id": metric.get("analysis_id"),
            "assay": donor["assay"],
            "depth_stratum": donor["depth_stratum"],
            "cell_fraction": fraction,
            "seed": seed,
            "keep_probability": x_keep_probability(fraction) if axis == "mlox" else y_keep_probability(fraction),
            "panel_completeness": completeness,
            "read_level_raw_log2_ratio_median": raw_median,
            "read_level_mLRR": mlrr,
            "estimated_cell_fraction": estimated,
            "read_level_zero_aware_mLRR": zero_mlrr,
            "zero_aware_estimated_cell_fraction": zero_estimated,
            "screen_tier": tier,
            "retained_signal": tier in signal,
            "borderline_or_extreme": tier in borderline,
            "event_call": "NOT_CALIBRATED",
            "negative_call": "NOT_EMITTED",
        })
    receipt = {
        "status": "TARGET_READ_CALIBRATION_SAMPLE_COMPLETE",
        "version": VERSION,
        "axis": axis,
        "sample_id": sample["sample_id"],
        "assay": donor["assay"],
        "depth_stratum": donor["depth_stratum"],
        "selected_interval_n": len(interval_ids),
        "fraction_zero_exact_interval_matches": len(interval_ids),
        "fraction_zero_filter_contract": "EXACT_MATCH",
        "alignment_records_seen_across_target_fetches": records_seen,
        "alignment_records_removed_by_frozen_flag_filter": filtered_records,
        "flag_filter": "UNMAP|SECONDARY|QCFAIL|DUP",
        "minimum_mapping_quality": 0,
        "minimum_base_quality": 0,
        "overlapping_mates_counted_independently": True,
        "autosomal_denominator_source": "frozen interval collection; BAM autosomes not rescanned",
        "pseudo_bams_written": False,
        "event_calls_emitted": False,
        "negative_calls_emitted": False,
        "bam_sha256": binding["bam"]["sha256"],
        "bam_index_sha256": binding["bam_index"]["sha256"],
    }
    return results, receipt


def run_grid(
    axis: str,
    collection: dict[str, Any],
    metrics: list[dict[str, Any]],
    panel: dict[str, Any],
    donor_plan: dict[str, Any],
    bam_manifest: list[dict[str, Any]],
    *,
    fractions_by_assay: dict[str, tuple[float, ...]],
    seeds: tuple[int, ...] = DEFAULT_SEEDS,
) -> tuple[list[dict[str, Any]], dict[str, Any], list[dict[str, Any]]]:
    donors = {row["sample_id"]: row for row in donor_plan["donors"]}
    samples = {row["sample_id"]: row for row in collection["samples"]}
    metric_by_id = {row["sample_id"]: row for row in metrics}
    bindings = {row["sample_id"]: row for row in bam_manifest}
    if len(bindings) != len(bam_manifest) or not bindings:
        raise ValueError("BAM manifest must contain unique samples")
    if not set(bindings).issubset(donors):
        raise ValueError("BAM manifest contains samples outside frozen donor plan")
    results: list[dict[str, Any]] = []
    receipts: list[dict[str, Any]] = []
    for sample_id, binding in sorted(bindings.items()):
        donor = donors[sample_id]
        assay = donor["assay"]
        if assay not in fractions_by_assay:
            raise ValueError(f"no read-level fraction grid for assay {assay}")
        rows, receipt = scan_sample(
            axis,
            binding,
            donor,
            samples[sample_id],
            metric_by_id[sample_id],
            panel["assays"][assay],
            fractions=fractions_by_assay[assay],
            seeds=seeds,
        )
        results.extend(rows)
        receipts.append(receipt)
        print(json.dumps({
            "stage": "TARGET_READ_CALIBRATION",
            "axis": axis,
            "completed": len(receipts),
            "total": len(bindings),
            "sample_id": sample_id,
        }), flush=True)
    grouped: dict[tuple[str, float], list[dict[str, Any]]] = defaultdict(list)
    for row in results:
        grouped[(row["assay"], row["cell_fraction"])].append(row)
    strata = []
    for (assay, fraction), rows in sorted(grouped.items()):
        borderline = sum(row["borderline_or_extreme"] for row in rows)
        strata.append({
            "assay": assay,
            "cell_fraction": fraction,
            "replicate_n": len(rows),
            "borderline_or_extreme_n": borderline,
            "borderline_or_extreme_rate": borderline / len(rows),
            "retained_signal_n": sum(row["retained_signal"] for row in rows),
        })
    summary = {
        "status": "TARGET_READ_CALIBRATION_COMPLETE_EVENT_CALLS_BLOCKED",
        "version": VERSION,
        "axis": axis,
        "sample_n": len(bindings),
        "replicate_n": len(results),
        "fractions_by_assay": {key: list(value) for key, value in sorted(fractions_by_assay.items())},
        "seeds": list(seeds),
        "strata": strata,
        "all_fraction_zero_filter_contracts_exact": all(
            row["fraction_zero_filter_contract"] == "EXACT_MATCH" for row in receipts
        ),
        "pseudo_bams_written": False,
        "event_calls_emitted": False,
        "negative_calls_emitted": False,
        "biological_specificity_validated": False,
        "interpretation": "target-only read-level analytical depth validation; orthogonal truth remains separate",
    }
    return results, summary, receipts


def _write_json(path: Path, payload: Any) -> None:
    with path.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    os.chmod(path, 0o600)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--axis", choices=("mlox", "mloy"), required=True)
    parser.add_argument("--collection", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--donor-plan", type=Path, required=True)
    parser.add_argument("--bam-manifest", type=Path, required=True)
    parser.add_argument("--fractions-json", type=Path)
    parser.add_argument("--seeds", default=",".join(map(str, DEFAULT_SEEDS)))
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    with gzip.open(args.collection, "rt", encoding="utf-8") as stream:
        collection = json.load(stream)
    metrics = json.loads(args.metrics.read_text(encoding="utf-8"))
    panel = json.loads(args.panel.read_text(encoding="utf-8"))
    donor_plan = json.loads(args.donor_plan.read_text(encoding="utf-8"))
    bam_manifest = json.loads(args.bam_manifest.read_text(encoding="utf-8"))
    raw_fractions = (
        json.loads(args.fractions_json.read_text(encoding="utf-8"))
        if args.fractions_json
        else DEFAULT_FRACTIONS[args.axis]
    )
    fractions = {key: _validate_fractions(tuple(value)) for key, value in raw_fractions.items()}
    seeds = tuple(int(value) for value in args.seeds.split(","))
    results, summary, receipts = run_grid(
        args.axis,
        collection,
        metrics,
        panel,
        donor_plan,
        bam_manifest,
        fractions_by_assay=fractions,
        seeds=seeds,
    )
    provenance = {
        "collection_sha256": sha256_file(args.collection),
        "metrics_sha256": sha256_file(args.metrics),
        "panel_sha256": sha256_file(args.panel),
        "donor_plan_sha256": sha256_file(args.donor_plan),
        "bam_manifest_sha256": sha256_file(args.bam_manifest),
        "implementation_sha256": sha256_file(Path(__file__)),
    }
    summary.update(provenance)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    _write_json(args.output_dir / f"{args.axis}.read-calibration.rows.PRIVATE.json", results)
    _write_json(args.output_dir / f"{args.axis}.read-calibration.sample-receipts.PRIVATE.json", receipts)
    _write_json(args.output_dir / f"{args.axis}.read-calibration.summary.json", summary)
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
