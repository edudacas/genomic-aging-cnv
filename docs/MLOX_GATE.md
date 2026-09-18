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

The collector must be run on the server because the downloaded Y collection
does not contain interval-level X values. A real-cohort mLOX result therefore
remains pending. BAM-level X read substitution at low fractions is required to
calibrate sensitivity and distinguish diffuse mLOX from target dropout.
