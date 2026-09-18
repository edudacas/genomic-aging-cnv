# Runtime and resource preparation, 2026-09-18

This change adds a persistent local native runtime bootstrap, an explicit
container build candidate, synthetic event validation, and a private inventory
builder. It does not authorize a CNV cohort run.

## Runtime

`runtime/bootstrap_local.py ROOT` installs CNVkit 0.9.12, pysam 0.23.3,
PyArrow 25.0.1 and numpy 1.26.4 into a new project-local Python environment.
It copies available R dependencies into an isolated library and installs missing
ones there, plus DNAcopy 1.86.0. Existing system packages are not upgraded.
FACETS-CH and pctGCdata archives are checked against the audited SHA-256 values;
the FACETS package-name metadata repair is recorded. snp-pileup is compiled from
the same pinned FACETS-CH commit. Sources, installation logs and resolved Python
versions remain in ROOT.

`runtime/validate_native.py --root ROOT --output NEW_DIRECTORY` runs portable
tests and the native synthetic LOSS, GAIN and CNLOH fixtures. Each uses the same
synthetic BAMs in both callers. `runtime/check_event_signals.py NEW_DIRECTORY`
additionally checks log-depth and allelic shifts against the ideal mixture.
These tests are not an empirical mosaic detection limit.

`runtime/freeze_runtime.py` records installed files and resolved ELF dependencies.
The local installation remains host-bound: system R, samtools, Python and shared
libraries are required. It is not a portable image.

The Docker candidate requires an explicit Ubuntu 24.04 digest through
`runtime/build_image.py --base ubuntu@sha256:DIGEST --output NEW_DIRECTORY`.
Its build runs the native tests and exports receipts and image ID. The current
recipe resolves apt/CRAN/Python dependencies at build time. Full source rebuild
locking remains open; an unbuilt recipe is not a versioned production image.

## Private inventory

`cnv_inventory.py --project STUDY_DIRECTORY --output NEW_PRIVATE_DIRECTORY`
checks exact workbook/linkage membership, person/sample/BAM uniqueness, FAP and
analysis equality, assay aliases and target BED digests. It preserves reported
sex separately from inferred sex and leaves unevaluated fields null/UNKNOWN.
No names, dates of birth or private study files belong in this repository.

Current local inventory: 1,329 samples, 1,208 Twist and 121 IDT, 66 assay/batch
combinations, no missing batch IDs, and 234 unknown reported sexes. The previously
confirmed Twist naming equivalence is retained, but capture provenance and CNV
normalization compatibility still require verification. Assay/batch combinations
are inventory strata, not an instruction to construct 66 independent PoNs.

The existing 80 external donors are retained as established study resources.
Reusing their donor eligibility does not create a CNV coverage reference, stable
marker panel or allelic-bias table. Those derived resources must be built and
validated for this workflow. Exact analysis-ID overlap is zero; the inventory
does not claim new person-identity validation from that check alone.

`collect_snv_qc.py --stages STAGES --output NEW_FILE.json.gz` reads completed QC
receipts, checks sample receipt bindings and selected artifact hashes, and exports
the compact summaries. It does not read BAMs or infer sex. Importing this bundle,
establishing UMI/library metadata, freezing masks/SNP catalogs/pseudo-depth grid,
verifying BAM/FASTA/BED identity, and defining X/Y/PAR handling remain gates.
