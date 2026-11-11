# Mouse (GRCm38) mLOY / mLOX screening

The interval collectors and screening gates run on mouse WES with three generic
extensions; human behaviour and defaults are unchanged.

## What changed

- **Assembly from the mask.** `collect_mloy_intervals.py` and
  `collect_mlox_intervals.py` accept masks for `GRCh38` or `GRCm38`; the assembly
  fixes the autosomal denominator (chr1–22 or chr1–19). A mask may list
  `autosomes` explicitly, and it must then match the assembly.
- **Labelled additional exclusions.** A mask may carry
  `additional_exclusions.{chrX,chrY}` as `{start, end, label}` entries. They are
  applied like PAR/XTR (whole parent interval excluded, status
  `EXCLUDED_<LABEL>`).
- **Configurable controls and metadata.** Both gates take `--control-group`
  (default `TP53_negative`) and either `--workbook` or `--metadata-tsv`
  (`analysis_id`, `age_at_collection`, `study_group`). The receipt records
  `metadata_tsv_sha256` or `workbook_sha256` accordingly.
- **Missing age in the mLOY gate.** As already done in the mLOX gate, a sample
  without `age_at_collection` is retained (`age_status=MISSING_RETAINED`) and only
  excluded from the young null.

## GRCm38 resource

`resources/sex_compartments.GRCm38.v1.json`:

- **PAR:** chrX ≥ 169,969,758 and chrY 90,745,844–91,644,698 (Ensembl GRCm38,
  0-based half-open).
- **XTR:** none; the mouse genome has no X-transposed region.
- **`AMPLICONIC_YQ`:** chrY 4,000,000–90,745,844 is excluded. Yq is multicopy
  (Sly, Ssty, Srsy) and its copy number varies by strain, so it is not a
  single-copy readout. Multicopy Yp clusters (for example Rbmy) fall out through
  the gates' per-interval callability rule.

## GC correction (optional)

`interval_gc.py` writes the GC fraction of every BED interval. Passing that table
to either collector with `--gc-table` adds `gc_fraction` and
`depth_over_autosome_gc` to every eligible interval. It also adds a per-sample
`gc_curve`, and records the table hash and method under `gc_correction`.

The method:

1. Take the sample's autosomal intervals with relative depth ≥ 0.1.
2. Compute log2(interval depth / autosomal mean) and take its median within
   each 1% GC bin. Only bins with ≥ 200 intervals are kept.
3. Centre the curve on the overall median.
4. Divide each sex-chromosome interval by `2**curve(GC)`. The curve is
   interpolated linearly and stays flat beyond the outermost bins.

Without `--gc-table` the output is unchanged. The gates use the corrected value
only with `--depth-metric gc`, and they refuse an uncorrected collection.

## Reference strata (optional)

`--metadata-tsv` may carry a `reference_stratum` column, for example `embryo`
or `postnatal`. Baselines then come from the most specific control set that
has at least 5 eligible controls, checked in this order:

1. batch + stratum
2. batch
3. stratum
4. assay

Without strata the order collapses to the original batch-then-assay rule. Each
row records `normalization_scope` and the set used. For mLOX, the female and
male baselines are chosen independently.

## Per-interval residual score (optional)

With `--score-metric interval_residual`, mLRR is the median over panel
intervals of log2(depth_over_autosome / control interval median). The control
set is the one chosen above, and a sample that belongs to that set is left out
of its own reference (`normalization_leave_one_out`).

This removes interval-to-interval capture differences that the default
`median_ratio` score keeps. The interpretive fractions follow the same models:

- LOY: `f = 1 - 2**mLRR`.
- LOX: `f = (1 - 2**mLRR) / (1 - male/female baseline ratio)`.

The default (`median_ratio`) is unchanged.

## Several control groups

`--control-group` accepts a comma-separated list, for example
`WT,Het,Homo` for all non-tumour muscle.

## Example

```bash
python3 interval_gc.py --bed capture.bed --fasta genome.fa --output interval_gc.tsv
python3 collect_mloy_intervals.py --stages STAGES --mask resources/sex_compartments.GRCm38.v1.json \
  --gc-table interval_gc.tsv --output mloy_collection.json.gz --expected N
python3 mloy_gate.py --collection mloy_collection.json.gz --sex-cohort sex_cohort.json \
  --metadata-tsv metadata.tsv --control-group WT,Het,Homo --depth-metric gc \
  --score-metric interval_residual --output-dir mloy/
```

The same pattern applies to the mLOX collector and gate.

## Mouse cohort check (108 WES samples, batches of 95 and 13)

The cohort was screened against an interval-level calibration done outside the
repository: GC correction, leave-one-out stage-matched references, and a
thinning-based LoD95 of 5% of cells for LOY and 3% for LOX in batch 1.

| Configuration | LOY: robust sigma of non-loss males, batch 1 | LOX: same, non-loss females, batch 1 |
|---|---|---|
| raw depth, WT controls, no strata, median_ratio | 0.064 | 0.028 |
| GC + strata, WT + Het + Homo muscle controls, median_ratio | 0.030 | 0.008 |
| GC + strata, WT + Het + Homo muscle controls, interval_residual | 0.018 | 0.005 |

With GC + strata + interval_residual:

- **mLOY:** every calibrated loss is flagged.
- **mLOX:** every calibrated batch-1 loss is flagged. The exception is a
  sex-discordant sample, which is reported as sex-chromosome anomaly review
  instead.
- **Extra signals:** the only flags beyond the calibrated set are weak signals,
  a low-depth sample already excluded from controls, and batch-2 samples whose
  batch has fewer than 5 controls per sex.
- **Fractions:** above ~5% of cells, the gate's fraction estimates agree with
  the calibrated ones to within about 1 percentage point.

## Known limitations

- **Autosomal denominator.** The autosomes enter as one mean depth. A per-interval
  autosomal ratio is not available, so tumour aneuploidy or strong residual GC
  shape can shift small estimates (under ~3% of cells) by a point or two.
- **Small batches.** A batch with fewer than 5 controls per sex still falls
  back to a wider set, so its own capture bias remains.
- **Young null.** The young null is the youngest third of controls. With
  strata and few young controls its sigma can be very small, which inflates
  `WEAK_SIGNAL_RETAINED`.
- **Calibration modules.** `mloy_calibration.py`, `mlox_calibration.py` and
  `sex_chromosome_specificity.py` still assume `TP53_negative` controls, the
  workbook and raw depth.
