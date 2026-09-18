# Native runtime and cohort resource freeze

Isolated local installation and synthetic tests are authorized. No cohort launch.
Required software: CNVkit 0.9.12, pysam 0.23.3, PyArrow 25.0.1,
samtools/HTSlib, R/DNAcopy, FACETS-CH and pctGCdata.

Require a dependency lock, source archive checksums, native LOSS/GAIN/CNLOH
tests, a container build and its image digest before marking the runtime ready.
Synthetic tests do not establish assay sensitivity or study readiness.

Study resources must be frozen separately by assay. Missing, unresolved or
discordant metadata remain HOLD; no automatic pooling of assay IDs.
