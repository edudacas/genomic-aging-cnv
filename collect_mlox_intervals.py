#!/usr/bin/env python3
"""Collect verified autosomal and X non-PAR/non-XTR interval summaries.

This read-only collector reuses completed SNV-QC interval artifacts. It does
not scan BAMs, recompute depth, infer sex, or call mLOX. Parent intervals that
touch frozen PAR or XTR boundaries are excluded whole so aggregated depth is
never interpolated across a boundary.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from collect_mloy_intervals import (
    AUTOSOMES,
    MEASURES,
    SAFE_SAMPLE,
    _add,
    _find_output,
    _overlap_length,
    _zero_totals,
    sha256_bytes,
    sha256_file,
)


def _load_x_mask(path: Path) -> tuple[dict[str, Any], list[tuple[int, int, str]]]:
    mask = json.loads(path.read_text(encoding="utf-8"))
    if mask.get("assembly") != "GRCh38" or mask.get("coordinate_system") != "0-based-half-open":
        raise ValueError("mask assembly/coordinate mismatch")
    if mask.get("xtr", {}).get("status") != "FROZEN_FOR_FUTURE_RECOMPUTATION_NOT_USED_IN_CURRENT_SEX_GATE":
        raise ValueError("XTR resource is not frozen")
    regions: list[tuple[int, int, str]] = []
    for start, end in mask["par"]["chrX"]:
        regions.append((int(start), int(end), "PAR"))
    for start, end in mask["xtr"]["chrX"]:
        regions.append((int(start), int(end), "XTR"))
    regions.sort()
    if any(not 0 <= start < end for start, end, _ in regions):
        raise ValueError("invalid X exclusion interval")
    if any(a[1] > b[0] for a, b in zip(regions, regions[1:])):
        raise ValueError("overlapping X exclusions")
    return mask, regions


def _process_sample(
    stages: Path,
    sample_id: str,
    binding: dict[str, Any],
    fingerprint: str,
    x_exclusions: list[tuple[int, int, str]],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if not SAFE_SAMPLE.fullmatch(sample_id):
        raise ValueError("unsafe sample identifier")
    directory = stages / "qc" / sample_id
    complete = directory / "COMPLETE.json"
    if directory.is_symlink() or complete.is_symlink():
        raise ValueError(f"{sample_id}: symlinked QC path")
    if sha256_file(complete) != binding["receipt_sha256"]:
        raise ValueError(f"{sample_id}: QC receipt mismatch")
    task = json.loads(complete.read_text(encoding="utf-8"))
    if task.get("status") != "COMPLETE" or task.get("fingerprint") != fingerprint:
        raise ValueError(f"{sample_id}: QC status/fingerprint mismatch")

    name, expected_digest = _find_output(task, "coverage.intervals.tsv.gz")
    relative = Path(name)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"{sample_id}: unsafe interval output path")
    interval_path = directory / relative
    if interval_path.is_symlink() or directory.resolve() not in interval_path.resolve().parents:
        raise ValueError(f"{sample_id}: unsafe interval artifact")
    compressed = interval_path.read_bytes()
    observed_digest = sha256_bytes(compressed)
    if observed_digest != expected_digest:
        raise ValueError(f"{sample_id}: interval artifact changed")
    try:
        text = gzip.decompress(compressed).decode("utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise ValueError(f"{sample_id}: invalid interval gzip") from exc

    autosomal = _zero_totals()
    parents: dict[tuple[int, int], dict[str, int]] = {}
    parent_segments: dict[tuple[int, int], list[tuple[int, int]]] = {}
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    required = {"chrom", "start", "end", "parent_start", "parent_end", *MEASURES}
    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        raise ValueError(f"{sample_id}: interval schema mismatch")
    for row in reader:
        chrom = row["chrom"]
        if chrom in AUTOSOMES:
            _add(autosomal, row)
            continue
        if chrom != "chrX":
            continue
        start, end = int(row["start"]), int(row["end"])
        parent_start, parent_end = int(row["parent_start"]), int(row["parent_end"])
        if not 0 <= parent_start <= start < end <= parent_end:
            raise ValueError(f"{sample_id}: invalid X parent/segment coordinates")
        key = (parent_start, parent_end)
        parents.setdefault(key, _zero_totals())
        parent_segments.setdefault(key, []).append((start, end))
        _add(parents[key], row)

    if autosomal["territory_bases"] <= 0 or autosomal["depth_sum"] <= 0:
        raise ValueError(f"{sample_id}: empty autosomal denominator")
    auto_mean = autosomal["depth_sum"] / autosomal["territory_bases"]
    catalog: list[dict[str, Any]] = []
    intervals: dict[str, dict[str, Any]] = {}
    for (start, end), totals in sorted(parents.items()):
        segments = sorted(parent_segments[(start, end)])
        if any(a[1] != b[0] for a, b in zip(segments, segments[1:])):
            raise ValueError(f"{sample_id}: noncontiguous split X parent")
        if segments[0][0] != start or segments[-1][1] != end:
            raise ValueError(f"{sample_id}: incomplete split X parent")
        if totals["territory_bases"] != end - start:
            raise ValueError(f"{sample_id}: X territory mismatch")
        overlap, labels = _overlap_length(start, end, x_exclusions)
        interval_id = f"chrX:{start}-{end}"
        status = "EXCLUDED_" + "_AND_".join(labels) if labels else "ELIGIBLE_NONPAR_NONXTR"
        catalog.append({
            "interval_id": interval_id,
            "chrom": "chrX",
            "start": start,
            "end": end,
            "territory_bases": end - start,
            "status": status,
            "masked_overlap_bases": overlap,
            "conservative_overexcluded_bases": (end - start - overlap) if overlap else 0,
        })
        if not labels:
            mean_depth = totals["depth_sum"] / totals["territory_bases"]
            intervals[interval_id] = {
                **totals,
                "mean_depth": mean_depth,
                "depth_over_autosome": mean_depth / auto_mean,
            }

    return ({
        "sample_id": sample_id,
        "qc_receipt_sha256": binding["receipt_sha256"],
        "coverage_intervals_sha256": observed_digest,
        "autosomal": {**autosomal, "mean_depth": auto_mean},
        "eligible_x_intervals": intervals,
    }, catalog)


def collect(stages: Path, mask_path: Path, destination: Path, expected: int, workers: int) -> dict[str, Any]:
    stages = stages.resolve()
    marker = stages / "QC_COMPLETE.json"
    initial_marker_hash = sha256_file(marker)
    receipt = json.loads(marker.read_text(encoding="utf-8"))
    samples = receipt.get("samples", {})
    if receipt.get("status") != "QC_COMPLETE" or len(samples) != expected:
        raise ValueError("QC completion/membership mismatch")
    mask, x_exclusions = _load_x_mask(mask_path)
    mask_hash = sha256_file(mask_path)

    results: dict[str, dict[str, Any]] = {}
    reference_catalog: list[dict[str, Any]] | None = None
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                _process_sample, stages, sample_id, binding,
                receipt.get("fingerprint"), x_exclusions,
            ): sample_id
            for sample_id, binding in sorted(samples.items())
        }
        completed = 0
        for future in as_completed(futures):
            sample_id = futures[future]
            sample, catalog = future.result()
            if reference_catalog is None:
                reference_catalog = catalog
            elif catalog != reference_catalog:
                raise ValueError(f"{sample_id}: X interval catalog differs across samples")
            results[sample_id] = sample
            completed += 1
            print(json.dumps({
                "stage": "MLOX_INTERVAL_COLLECTION", "completed": completed,
                "total": expected, "sample_id": sample_id,
            }), flush=True)

    if sha256_file(marker) != initial_marker_hash:
        raise ValueError("QC completion receipt changed during collection")
    assert reference_catalog is not None
    eligible = [r for r in reference_catalog if r["status"] == "ELIGIBLE_NONPAR_NONXTR"]
    excluded = [r for r in reference_catalog if r["status"] != "ELIGIBLE_NONPAR_NONXTR"]
    payload = {
        "status": "MLOX_INTERVAL_COLLECTION_COMPLETE",
        "schema_version": 1,
        "selected_samples": len(results),
        "qc_complete_sha256": initial_marker_hash,
        "qc_fingerprint": receipt.get("fingerprint"),
        "mask": mask,
        "mask_sha256": mask_hash,
        "aggregation_policy": "WHOLE_PARENT_INTERVAL_EXCLUDED_IF_ANY_PAR_OR_XTR_OVERLAP",
        "bams_scanned": False,
        "depth_recomputed": False,
        "mlox_called": False,
        "x_interval_catalog": reference_catalog,
        "eligible_x_interval_count": len(eligible),
        "eligible_x_territory_bases": sum(r["territory_bases"] for r in eligible),
        "excluded_x_parent_count": len(excluded),
        "excluded_x_parent_territory_bases": sum(r["territory_bases"] for r in excluded),
        "conservative_overexcluded_bases": sum(r["conservative_overexcluded_bases"] for r in excluded),
        "samples": [results[sid] for sid in sorted(results)],
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("xb") as raw:
        os.chmod(destination, 0o600)
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as stream:
            stream.write(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8"))
    result = {
        "status": payload["status"],
        "samples": len(results),
        "eligible_x_intervals": len(eligible),
        "eligible_x_territory_bases": payload["eligible_x_territory_bases"],
        "conservative_overexcluded_bases": payload["conservative_overexcluded_bases"],
        "output": str(destination),
        "sha256": sha256_file(destination),
    }
    print(json.dumps(result, sort_keys=True), flush=True)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stages", type=Path, required=True)
    parser.add_argument("--mask", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected", type=int, default=1329)
    parser.add_argument("--workers", type=int, default=4, choices=range(1, 5))
    args = parser.parse_args()
    collect(args.stages, args.mask, args.output, args.expected, args.workers)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
