# mLOX analytical coverage calibration

`mlox_calibration.py` separates two validation layers that must not be conflated.

## Interval-summary Monte Carlo calibration

The `coverage` mode uses observed X interval `depth_sum` values from eligible,
within-reference TP53-negative females. Donors are selected deterministically
within low, middle and high autosomal-depth tertiles for each assay. For each
prespecified cell fraction and seed, the X read-base count is binomially thinned
with probability:

`P(retain X) = 1 - cell_fraction / 2`.

The autosomal denominator is held fixed. This reproduces the expected depth
component of losing one X copy in a fraction of otherwise diploid female cells.
The frozen mLOX baselines and thresholds are then applied without refitting.
The output reports sensitivity for retained signal and for the
borderline-or-extreme endpoint, false-positive donors at fraction zero, and a
provisional assay-specific limit of detection.

This mode preserves real interval-to-interval and donor-to-donor coverage
variation, while adding stochastic thinning. It does not reproduce read-level
mapping effects, duplicate handling, base-level threshold transitions or
haplotype-specific allelic imbalance. Its status therefore ends in
`BAM_VALIDATION_REQUIRED`.

## BAM depth calibration

The `thin-bam` mode preserves every non-X alignment and deterministically thins
mapped X alignments by read-group/query-name. Run the existing compact coverage
summarizer and the frozen mLOX gate on each resulting BAM. The BAM output is for
depth QC only: cross-chromosome pairs can become orphaned, and the BAM must not
be reused for generic variant calling.

The BAM mode rejects UMI-tagged input because molecule-aware thinning requires a
separate design. It records the input checksum, implementation checksum,
retained X/non-X counts, output checksums and an explicit statement that
haplotype-specific imbalance was not simulated.

Example interval calibration:

```bash
python3 mlox_calibration.py coverage \
  --collection mlox-x-intervals.GRCh38.v1.PRIVATE.json.gz \
  --metrics mlox.metrics.PRIVATE.json \
  --panel mlox.callable_panel.json \
  --output-dir mlox_coverage_calibration.PRIVATE \
  --collection-sha256 98be0cfb9c2ef7a9988be71785b48308fa9f36477da995e6d522d0b1ba597c1a
```

Example one-BAM depth simulation:

```bash
python3 mlox_calibration.py thin-bam \
  --bam sample.bam \
  --bam-sha256 SHA256 \
  --sample-name BAM_SM \
  --cell-fraction 0.05 \
  --seed 17 \
  --output-dir sample.f005.seed17.PRIVATE
```

## FULL1329 interval-summary result

The frozen FULL1329 X collection was calibrated with 18 deterministic donors:
three donors from each autosomal-depth tertile in each of the two assays. Ten
cell fractions and five seeds produced 900 replicates. At fraction zero, no
donor entered a retained-signal tier in either assay.

Using the prespecified rule (pooled borderline-or-extreme sensitivity at least
95% and every depth stratum at least 90%), the provisional coverage limits were:

- `TWIST_WES_OS_T7`: 5% cell fraction (44/45 borderline-or-extreme replicates);
- `MGI_IDT_EXOME`: 15% cell fraction (45/45 replicates).

At 5%, IDT reached 34/45 (75.6%); at 7.5% and 10%, it reached 40/45 (88.9%).
This assay difference must be retained in downstream interpretation. These are
interval-summary limits, not final BAM-validated limits. Aggregate provenance is
recorded in `validation/mlox_coverage_calibration_20260918.json`.

## Interpretation boundary

No mode emits negative or mLOX event calls. Coverage calibration estimates
analytical detection behavior only. Haplotype-aware simulation or true and
orthogonally confirmed mLOX samples remain necessary to validate allelic
imbalance, biological specificity and clinical sensitivity.

## Leakage-controlled reference behavior

`sex_chromosome_specificity.py --axis mlox` assigns TP53-negative reference
females to deterministic folds stratified by assay, batch and autosomal-depth
tertile. For each fold, callable intervals, normalization baselines and the
young-null thresholds are rebuilt without that fold; held-out samples are then
scored exactly once. A TP53-negative donor can still carry real age-related
mLOX, so the retained-signal rate is not labeled a false-positive rate.

## Target-only read-level validation

`sex_chromosome_read_calibration.py --axis mlox` scans only frozen callable X
targets, applies deterministic template thinning in memory and retains the
frozen autosomal denominator. It writes no pseudo-BAMs. Before any altered
fraction is accepted, fraction zero must reproduce every frozen target
`depth_sum` exactly. The boundary grid is Twist 2.5/5/7.5% and IDT 10/15/20%,
plus fraction zero and three prespecified seeds.
