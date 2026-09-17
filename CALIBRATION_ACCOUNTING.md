# Calibration accounting contract

`cnv_calibration.summarize_calibration` summarizes a frozen set of planned
cases and externally adjudicated outcomes. It does not generate truth, match
caller intervals to truth, or close the study gate.

Every source sample must have role `calibration` in the validated manifest.
Study and PoN donors are excluded. A plan supplies a protocol SHA-256 and cases
with unique case ID, source sample ID, truth source (SIMULATED or ORTHOGONAL),
event class (LOSS, GAIN, CNLOH or NULL_CONTROL), clone fraction, prespecified
size/depth/informative-locus bands and a truth-artifact path/hash.

Outcomes carry case ID, the same protocol digest, status and result artifact.
With file checks enabled, truth/result bytes must match their declared hashes.
This does not prove the truth assay or adjudication protocol is valid.

## Outcomes and denominators

- Positive cases: DETECTED or MISSED.
- Null controls: FALSE_POSITIVE or CLEAN.
- Either: NOT_EVALUABLE or NOT_RUN.
- A missing outcome remains NOT_RUN, never MISSED or CLEAN.

For each exact stratum/truth-source/event/fraction/size/depth/density cell, the
summary reports planned cases, evaluable cases, status counts, case IDs and
unique source samples. The detection proportion uses DETECTED/(DETECTED+MISSED);
the false-positive proportion uses FALSE_POSITIVE/(FALSE_POSITIVE+CLEAN).
Both are conditional on evaluability, which is reported separately. No
evaluable denominator means null, not zero.

Simulated and orthogonal results are not pooled. Orthogonally confirmed events
may have an unknown measured clone fraction; this remains null. Simulated
positive cases require a declared fraction. Null controls do not receive a
fabricated clone fraction.

Multiple fractions/downsampling replicates from one donor are not assumed
independent. No confidence intervals or detection-limit selection are inferred
from these raw proportions. Matching rules, event-aware simulation validity,
held-out adequacy, uncertainty analysis and gate decisions remain separate work.
