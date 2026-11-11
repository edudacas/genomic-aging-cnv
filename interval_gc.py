#!/usr/bin/env python3
"""Write the GC fraction of every BED interval for the collectors' --gc-table option.

GC = (G + C) / (A + C + G + T), case-insensitive; N and other symbols are ignored,
and an interval with no A/C/G/T is written as NA. Coordinates are copied from the
BED (0-based half-open), so they match the coverage.intervals.tsv.gz rows built from
the same BED.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


def gc_fraction(sequence: str) -> float | None:
    sequence = sequence.upper()
    gc = sequence.count("G") + sequence.count("C")
    acgt = gc + sequence.count("A") + sequence.count("T")
    return gc / acgt if acgt else None


def write_gc_table(bed: Path, fasta: Path, output: Path) -> dict[str, object]:
    import pysam

    reference = pysam.FastaFile(str(fasta))
    rows = 0
    missing = 0
    with bed.open(encoding="utf-8") as source, output.open("x", encoding="utf-8") as sink:
        sink.write("chrom\tstart\tend\tgc\n")
        for line in source:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end = line.split("\t")[:3]
            value = gc_fraction(reference.fetch(chrom, int(start), int(end)))
            missing += value is None
            sink.write(f"{chrom}\t{int(start)}\t{int(end)}\t{'NA' if value is None else f'{value:.6f}'}\n")
            rows += 1
    return {
        "rows": rows,
        "no_acgt_rows": missing,
        "output": str(output),
        "sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bed", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True, help="indexed reference FASTA")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(write_gc_table(args.bed, args.fasta, args.output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
