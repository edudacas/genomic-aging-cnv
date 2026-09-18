# mLOY analytical coverage calibration

`mloy_calibration.py` separates interval-summary calibration from BAM-level
depth validation. Neither mode emits event or negative calls.

## Interval-summary Monte Carlo calibration

The `coverage` mode selects eligible, within-reference TP53-negative males
deterministically within low, middle and high autosomal-depth tertiles for each
assay. For every prespecified fraction and seed, observed Y interval
`depth_sum` is binomially thinned with:

`P(retain Y) = 1 - cell_fraction`.

The autosomal denominator remains fixed. Frozen mLOY baselines and thresholds
are applied without refitting. Outputs include retained-signal and
borderline-or-extreme sensitivity, fraction-zero false-positive donors, and a
provisional assay-specific limit of detection.

The summary mode preserves real donor and interval coverage heterogeneity but
does not reproduce read-level mapping, duplicate handling or base-level
threshold transitions. Its result remains provisional until BAM validation.

## BAM depth calibration

The `thin-bam` mode preserves every non-Y alignment and deterministically thins
mapped Y records by read-group/query-name. Its BAMs are restricted to depth-QC
calibration. Cross-chromosome pairs can become orphaned, so these BAMs must not
be reused for generic variant calling.

UMI-tagged inputs are rejected pending a molecule-aware design. Each attempt
records input and implementation hashes, retained Y/non-Y counts and output
hashes.

Example interval calibration:

```bash
python3 mloy_calibration.py coverage \
  --collection mloy-y-intervals.GRCh38.v1.PRIVATE.json.gz \
  --metrics mloy.metrics.PRIVATE.json \
  --panel mloy.callable_panel.json \
  --output-dir mloy_coverage_calibration.PRIVATE \
  --collection-sha256 74c019f2a95f24bceeeeaec824eb2a0614e7173a26787eb94212c7bf04b2d8de
```

Example one-BAM depth simulation:

```bash
python3 mloy_calibration.py thin-bam \
  --bam sample.bam \
  --bam-sha256 SHA256 \
  --sample-name BAM_SM \
  --cell-fraction 0.05 \
  --seed 17 \
  --output-dir sample.f005.seed17.PRIVATE
```

## FULL1329 interval-summary result

The frozen FULL1329 Y collection was calibrated with 18 deterministic donors:
three donors from each autosomal-depth tertile in each of the two assays. Ten
cell fractions and five seeds produced 900 replicates. At fraction zero, no
donor entered a retained-signal tier in either assay.

Using the prespecified rule (pooled borderline-or-extreme sensitivity at least
95% and every depth stratum at least 90%), the provisional coverage limits were:

- `MGI_IDT_EXOME`: 7.5% cell fraction (45/45 borderline-or-extreme replicates);
- `TWIST_WES_OS_T7`: 10% cell fraction (45/45 replicates).

At 5%, IDT reached 22/45 borderline-or-extreme replicates while all 45 entered
at least a retained-signal tier. Twist reached 17/45 at 5% and 32/45 at 7.5%.
Thus 5% remains a sensitive retention zone, not a reliable event threshold.
These are interval-summary limits, not final BAM-validated limits. Aggregate
provenance is recorded in
`validation/mloy_coverage_calibration_20260918.json`.

## Interpretation boundary

The interval-summary limit is analytical and coverage-specific. BAM validation
is required before it can unlock a research event threshold. True or
orthogonally confirmed mLOY samples remain necessary for biological specificity
and clinical sensitivity. The substantially smaller callable Y territory means
uncertainty must be reported by assay and depth stratum.
