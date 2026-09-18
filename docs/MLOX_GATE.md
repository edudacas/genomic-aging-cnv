# Sensitive mLOX screening gate

`collect_mlox_intervals.py` reads the immutable, completed SNV-QC interval
summaries and collects X non-PAR/non-XTR depth without scanning BAMs or
recomputing coverage. Parent intervals touching frozen PAR or XTR boundaries
are excluded whole; aggregated depth is never interpolated.

`mlox_gate.py` then builds an assay-specific callable-X panel and preserves:

- every reported or inferred female;
- every sample held for a sex-chromosome anomaly;
- continuous X/autosome, interval-breadth and technical-completeness metrics.

The interpretive mosaic fraction follows the prespecified read-substitution
model on linear depth:

`observed X = female_reference * (1 - f) + male_reference * f`.

Female and male baselines are batch-specific when each has at least five
TP53-negative controls and otherwise assay-specific. The same sensitive,
non-negative tiers used for mLOY are emitted for mLOX. `event_call_LOX` always
remains `NOT_CALIBRATED`.

## FULL1329 sensitive run

The receipt-verified X collection completed for all 1,329 samples with 6,656
eligible parent intervals covering 1,261,299 bases. Assay-specific callable
panels retained 6,304 IDT intervals and 6,312 Twist intervals.

The screening gate retained 1,056 female-like or sex-chromosome-anomaly rows:

- four extreme candidates;
- 18 borderline candidates;
- 41 additional weak signals;
- one sex-chromosome anomaly review;
- 992 within the current uncalibrated reference distribution;
- no negative or mLOX event calls.

All 1,056 rows had evaluable interval metrics. One missing collection age was
retained in the metrics and excluded only from the young-null definition. The
four extreme signals were distributed across approximately 63–80% of callable
X targets rather than being confined to a few intervals. Borderline signals
include very young controls, so the tiers cannot be interpreted as prevalence.

`mlox_calibration.py` now implements prespecified interval-summary Monte Carlo
thinning and depth-only BAM X thinning. Interval-summary results remain
provisional until representative BAMs are thinned and rerun through compact
coverage collection. Neither mode models haplotype-specific allelic imbalance;
true or orthogonally confirmed mLOX samples remain necessary for biological and
clinical validation. See `docs/MLOX_CALIBRATION.md`.
