# Sex-chromosome validation before cohort interpretation

This stage validates analytical screening behavior for mLOX and mLOY. It does
not authorize event, negative, somatic-origin or clinical calls.

## 1. Held-out reference-control analysis

Run `sex_chromosome_specificity.py` independently for mLOX and mLOY using the
frozen interval collection, sex-gate cohort and private study workbook. The
program creates deterministic folds stratified by assay, batch and
assay-specific autosomal-depth tertile. For every fold it reconstructs the
callable panel, normalization and young-null thresholds without the held-out
samples and scores each held-out sample once.

The output is a retained-signal rate among held-out TP53-negative reference
controls, with Wilson intervals. It is not a false-positive rate: TP53-negative
people may carry genuine age-related sex-chromosome mosaicism.

```bash
python3 sex_chromosome_specificity.py \
  --axis mlox \
  --collection mlox-x-intervals.GRCh38.v1.PRIVATE.json.gz \
  --sex-cohort sex-gate.cohort.PRIVATE.json \
  --workbook cohort.PRIVATE.xlsx \
  --folds 5 \
  --output-dir mlox_specificity_cv.PRIVATE

python3 sex_chromosome_specificity.py \
  --axis mloy \
  --collection mloy-y-intervals.GRCh38.v1.PRIVATE.json.gz \
  --sex-cohort sex-gate.cohort.PRIVATE.json \
  --workbook cohort.PRIVATE.xlsx \
  --folds 5 \
  --output-dir mloy_specificity_cv.PRIVATE
```

## 2. Target-only read-level boundary calibration

The BAM manifest is private JSON. Each row contains only pseudonymous linkage
and immutable artifacts:

```json
{
  "sample_id": "GA_...",
  "assay": "MGI_IDT_EXOME",
  "bam_sample_name": "declared BAM SM",
  "bam": {"path": "/absolute/sample.bam", "sha256": "..."},
  "bam_index": {"path": "/absolute/sample.bam.bai", "sha256": "..."}
}
```

Samples must be a subset of the frozen interval-summary donor plan. Use one
representative donor from each low/middle/high depth stratum per assay. Default
grids are:

- mLOX Twist: 0, 2.5, 5 and 7.5%; mLOX IDT: 0, 10, 15 and 20%;
- mLOY IDT: 0, 5, 7.5 and 10%; mLOY Twist: 0, 7.5, 10 and 15%;
- seeds 17, 29 and 43.

The scanner fetches only callable X or Y targets and performs template-level
thinning in memory. It does not create pseudo-BAMs. The original autosomal
denominator is retained. Fraction zero must exactly reproduce every frozen
interval `depth_sum`; disagreement blocks the sample and indicates a filter or
input-contract mismatch.

```bash
python3 sex_chromosome_read_calibration.py \
  --axis mloy \
  --collection mloy-y-intervals.GRCh38.v1.PRIVATE.json.gz \
  --metrics mloy.metrics.PRIVATE.json \
  --panel mloy.callable_panel.json \
  --donor-plan mloy.calibration.plan.PRIVATE.json \
  --bam-manifest mloy.read-calibration.bams.PRIVATE.json \
  --output-dir mloy_read_calibration.PRIVATE
```

Repeat with the corresponding mLOX inputs and `--axis mlox`.

## 3. Orthogonal WES-array cases

Internal WES-array cases should remain external to threshold fitting. Include
both positive and negative cases when available. Transfer a private,
de-identified manifest with:

- pseudonymous WES sample ID and axis (`mlox` or `mloy`);
- array platform, version and genome build;
- positive, negative or indeterminate truth label and interpretation date;
- event chromosome and 0-based half-open coordinates for positives;
- array logR/BAF or estimated mosaic fraction when available;
- tissue and days between WES and array collections;
- source artifact SHA-256 and whether review was blinded to WES;
- assay, batch and WES depth stratum.

Do not include name, date of birth, FAP or other direct identifiers. Array
negatives are only strong negative truth within the array's validated mosaic
sensitivity; an array-negative low-fraction WES signal remains indeterminate,
not automatically false positive.

## Interpretation boundary

After cross-validation, target-read calibration and orthogonal review, freeze
assay-specific research thresholds. Below-limit signals remain
`SCREEN_RETAINED_BELOW_LOD`; above-limit signals may become
`CANDIDATE_ABOVE_ASSAY_LOD`. Neither label is a confirmed somatic event.
