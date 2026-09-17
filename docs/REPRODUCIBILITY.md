# Reproducibility and production boundary

## Portable verification

The default test command exercises pure Python contracts, evidence semantics,
receipts, resume behavior, reconciliation and exports:

```bash
python -m pip install -r requirements-test.txt
python -m unittest discover -s tests -q
```

The fixtures are generated in temporary directories and contain no study data.

## Native verification

Native smoke scripts require explicit paths to external software. The historical
receipts in `validation/` describe the versions used in the development session,
but the temporary runtime itself is not a distributed or reproducible release.

Reproduction requires independently installing and recording:

- CNVkit 0.9.12 and pysam 0.23.3;
- samtools and the pinned `snp-pileup` source revision;
- R, DNAcopy, FACETS-CH and pctGCdata;
- PyArrow for exact JSONL/Parquet round trips;
- system and shared-library versions.

Do not interpret a synthetic smoke pass as assay validation. The artificial
reference, BAMs and events test interfaces and evidence flow rather than the
study's sensitivity or specificity.

## Required study-side inputs

Private study execution requires resources that are intentionally absent:

- authorized participant/sample linkage and BAM locations;
- exact GRCh38 reference and capture territories;
- assay-stratified PoN membership;
- SNP, depth-grid and mask resources;
- calibration truth and outcome manifests;
- clinical covariates and review decisions.

These inputs should be mounted at execution time and represented by hashes in
receipts. They must never be copied into this repository.

## Release criterion

A production release requires a frozen container digest, resource manifest,
successful native tests, a representative real-BAM canary, calibrated event
rules, measured CPU/RAM/I/O/disk requirements, failure/restart testing and an
explicit authorization receipt. None is implied by this repository snapshot.
