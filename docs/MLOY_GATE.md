# Sensitive mLOY screening gate

`mloy_gate.py` freezes an assay-specific callable-Y panel and emits continuous,
uncalibrated mLRR-Y screening metrics. It never emits negative or event calls.

## Preservation of borderline biology

Sex inference is context, not an exclusion filter. Reported males, inferred
males and samples held for a sex-chromosome anomaly remain in the output. Only
sex-gate-passing male-like samples contribute to the callable panel and control
reference.

The pre-calibration tiers are:

- `EXTREME_CANDIDATE_UNCALIBRATED`: below the young-control three-sigma
  threshold with interpretive loss fraction at least 10%;
- `BORDERLINE_CANDIDATE_UNCALIBRATED`: interpretive loss fraction at least 5%;
- `WEAK_SIGNAL_RETAINED`: below the young-control two-sigma threshold;
- `WITHIN_REFERENCE_UNCALIBRATED`: not flagged by the exploratory thresholds;
- explicit technical, absent-signal and sex-chromosome-anomaly review states.

`WITHIN_REFERENCE_UNCALIBRATED` is not a negative result. Every row retains the
continuous metrics and `event_call_LOY=NOT_CALIBRATED`.

## Target breadth and zero-aware metrics

The gate retains the legacy median interval metric for continuity and adds:

- a territory-weighted zero-aware Y/autosome metric;
- the number of zero-depth targets;
- interval-specific residuals against batch controls when available;
- the fraction of targets below the control median, -0.05 log2 and -0.10 log2;
- residual median and 10th/90th percentiles.

This distinguishes chromosome-wide depression from a few failing targets.
Samples with incomplete panels are retained as technically non-evaluable, and
complete absence of Y remains a mandatory manual-review state.

## FULL1329 sensitive rerun

The verified existing interval collection was re-evaluated without reading
BAMs:

- 273 reference male samples plus one sex-chromosome anomaly review;
- 272 eligible metrics, one incomplete panel and one absent-Y manual review;
- nine extreme candidates;
- 26 borderline candidates between the current exploratory fraction tiers;
- eight additional weak signals retained;
- no negative or mLOY event calls.

The candidate distribution includes very young controls. Therefore the tiers
remain technical screening categories and cannot be interpreted as prevalence.
BAM-level Y non-PAR downsampling by assay, depth and low fraction is required.
