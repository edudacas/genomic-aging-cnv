# Genomic Aging: mosaic CNV/LOH pipeline

Research software accompanying the Genomic Aging in Li-Fraumeni Syndrome
study. This repository contains the developing WES mosaic/somatic CNV and LOH
workflow, separate from the constitutional Rainbow CNV analysis and from the
staged SNV/indel pipeline in `genomic-aging-snv`.

## Status and validation boundary

This is a **development snapshot, not a production release**. No CNV/LOH
cohort run has been launched for the 1,329 study samples from this code.

With the pinned Python dependencies installed, the current portable regression
suite contains 195 passing tests. Historical
synthetic native tests also exercised CNVkit 0.9.12/DNAcopy, FACETS-CH,
`snp-pileup`, paired-caller comparison and Parquet export; their receipts are
retained under `validation/`. Those tests demonstrate implementation and file
boundary behavior. They do not establish clinical validity, study sensitivity,
specificity, a limit of detection or fitness of a study panel of normals.

The isolated native runtime used for those tests was temporary and is not
distributed here. Production readiness still requires a reproducible container,
assay-matched resources and PoNs, real-BAM canaries, event-aware calibration and
frozen classification thresholds. See `DEVELOPMENT_STATE.md`.

## Intended workflow

1. Validate GRCh38 BAM/index identity, sample roles, assay/capture strata and
   frozen resources.
2. Generate CNVkit target/antitarget coverage and FACETS SNP/pseudo-depth counts.
3. Build independent, assay-matched PoN references and allelic-bias summaries
   without study or held-out calibration leakage.
4. Run CNVkit normalization/CBS and the audited FACETS-CH GRCh38 adapter.
5. Reconcile caller segments while retaining caller-only, discordant and
   not-evaluable evidence.
6. Add TP53/17p regional and fragment evidence; retain Manta as exploratory
   structural-variant evidence.
7. Calibrate loss, gain and cnLOH sensitivity by event size, cell fraction,
   depth and informative-marker density before biological classification.
8. Export native evidence plus typed JSONL/Parquet for review and analysis.

Missing or unevaluated evidence is never converted to zero. Caller concordance
is not orthogonal truth, and candidate segments are not automatically assigned
somatic origin or a TP53 second hit.

## Repository layout

- `cnv_contract.py`, `cnv_coordinates.py`: input, role, coordinate and resource
  contracts.
- `cnv_pileup.py`, `cnv_catalog.py`: SNP/pseudo-depth evidence handling.
- `cnv_training.py`, `cnv_train_panel.py`: PoN marker and allelic-bias training.
- `facets_*.R`, `cnv_run_joint.py`: audited FACETS-CH preparation and execution.
- `cnv_reference.py`, `cnv_normalize.py`, `cnv_cbs.py`: CNVkit/CBS stages.
- `cnv_compare.py`, `cnv_compare_export.py`: provenance-bound reconciliation
  and evidence export.
- `cnv_downsample.py`, `cnv_simulation.py`, `cnv_calibration.py`: depth stress,
  ideal diploid mixture expectations and calibration accounting.
- `cnv_regions.py`, `cnv_phase_*.py`, `cnv_manta.py`: regional, fragment and
  exploratory SV evidence.
- `tests/`: synthetic unit, contract and optional native smoke tests.
- `validation/`: historical synthetic native-test receipts; no patient data.

## Run the portable tests

Python 3.11 or 3.12 is recommended:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-test.txt
python -m unittest discover -s tests -q
```

Native smoke tests additionally require separately installed and licensed
CNVkit, samtools, `snp-pileup`, R, DNAcopy, FACETS-CH and `pctGCdata`. Their
paths are explicit command-line inputs; the test suite does not install them.
See `docs/REPRODUCIBILITY.md` and `UPSTREAM_AUDIT.md`.

## Production boundary

Before a 1,329-sample run, the study team must freeze and validate:

- the exact sample/assay manifest and exclusion or HOLD reasons;
- capture-specific target/antitarget territory and X/Y/PAR policy;
- independent assay-matched PoNs and stable-marker resources;
- masks, reference FASTA, SNP catalog and pseudo-depth grid;
- QC, contamination, duplication and UMI policy;
- event classification rules and evaluable denominators;
- simulated and, where available, orthogonally confirmed calibration cases;
- a representative real-BAM canary, resource budget and resumable cohort
  orchestration.

The explicit CLI runs individual development stages only. It intentionally does
not auto-chain or authorize a cohort launch.

## Data, privacy and citation

No workbook, patient identifier, FAP, name, date of birth, BAM/VCF, clinical
linkage or cohort result is included. Synthetic fixtures are not patient data.
Do not commit study data or server configuration to this repository.

For review, cite the exact Git commit together with the relevant validation
receipt. Publication DOI, archival release and software licensing remain author
decisions. No third-party binaries, packages or reference datasets are
redistributed here.
