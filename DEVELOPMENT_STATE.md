# Development state

The complete manuscript/figure objective remains open. This folder now has
input contracts, coordinate checks, native command constructors and a locked
task runner with hash-verified resume and failure receipts.

## Repository snapshot, 2026-09-17

The sanitized repository passes 195 Python unit/contract tests with PyArrow
25.0.1 and pysam 0.23.3. The ideal diploid mixture helper now covers loss,
gain and cnLOH expectations jointly for depth and haplotype balance. This is a
mathematical fixture model, not an empirical WES error or sensitivity model.

The previously authorized native runtime lived under `/tmp` and is no longer
available. Historical synthetic receipts remain under `validation/`, but native
gain/cnLOH reruns were interrupted and are not claimed as completed evidence.
A permanent container, assay-matched resources, real-BAM canary, calibration,
classification rules and cohort orchestration remain required before production.
No CNV cohort run was launched.

## Verified

### Calibration denominator accounting

Implemented `summarize_calibration` over frozen planned cases and externally
adjudicated outcomes. It enforces calibration-only donors, protocol matching,
unique cases/outcomes and declared artifact hashes. Simulation/orthogonal
evidence and null/positive trials are separated; missing runs are NOT_RUN,
non-evaluable results remain explicit, and proportions have explicit evaluable
denominators. Repeated cases from one source are not independent by assumption.
189 tests pass. This is descriptive bookkeeping, not a completed calibration
experiment, truth validation, event-aware simulator or automatic gate.

### Whole-BAM template depth subsampling

Re-read manuscript calibration and mosaic-CNV/LOH sections: minimum event size,
density, amplitude and quality must be prespecified, and clonal sensitivity
must be measured rather than inferred from implementation tests.
Added deterministic whole-BAM template subsampling with hash-pinned input,
single-sample RG checks, shared mate decisions, nested fractions, output index
and receipts. UMI-tagged BAMs fail closed pending molecule-aware design.
Native pysam tests preserve both mate records, reject wrong sample/UMI and
exercise output indexing. 181 tests pass. This is depth stress testing only,
not event-aware simulation or calibration completion. Independent positive
calibration resources were requested; none are assumed to exist.

### Comparison-to-Parquet integration

Added `export_comparison` with pinned receipt/output/input checks, explicit
sample identity, row/status-count checks and rejection of unsupported biological
claims. Full nested comparison payloads and null evaluability survive typed
Parquet/canonical JSONL export. All six intervals from the real paired synthetic
run passed round-trip field verification. Separate compare-callers and
export-comparison CLI commands do not launch native callers. 176 tests pass.
Evidence: `validation/native_comparison_export.json`. This is not clinical
classification, calibration, de-identification or study execution.

### Provenance-bound paired-caller comparison

Both native stage receipts now bind sample, stratum, BAM digest, manifest digest
and build. `cnv_compare.compare_callers` verifies externally hash-pinned stage
receipts, native fingerprints/inputs/output bytes and implementation files
before comparison. The same-BAM native fixture passed this API: four CNVkit-only
boundary intervals, one concordant loss-direction interval and one within-band
interval; no event counts inferred. 170 tests pass, including same-ID/different-
BAM rejection and changed input/output rejection. Old receipts lack this new
binding and are not accepted by the comparison stage. See
`validation/native_verified_comparison.json`.

### Same-BAM native paired-caller integration

The synthetic joint-chain test now optionally runs CNVkit reference/fix/CBS on
the same study and two PoN BAMs and reconciles its CNS with native FACETS output.
An initial all-A reference exposed complete GC-based bin rejection. That failed
attempt is preserved; balanced-GC reference and matching read sequences retain
the planned count signal and pass. Normalization now explicitly rejects empty
or malformed/nonfinite normalized tables, including at later verification.
Both callers produced two segments; six comparison intervals include concordant
loss and within-band signal, without event calls or evaluable denominators.
165 tests pass. Evidence: `validation/native_paired_callers.json`.

### CNVkit command-line stage interface

Added execution-gated run-reference/run-normalize/run-cbs commands, exact config
keys, and repeatable implementation-file inputs. No stage automatically launches
the next. The real native fixture now invokes CBS through cnv_cli.py and passes
both target-only and hybrid modes. 164 unit/contract tests pass, including
pre-input-read execution gates and rejection of configuration output overrides.
`CNVKIT_STAGES.md` documents schemas and explicitly synthetic example values.

### Verified normalization to CBS stage

Added read-only `verify_normalization` by sharing reconstruction logic with
`normalize_sample`. Added `cnv_cbs.run_cbs`: checks normalization/reference,
uses explicit hashed R runtime, validates CNS bounds/overlap/probe counts,
rechecks inputs and writes completion only after those checks. Native CNVkit
smoke now uses the stage API and passes target-only and hybrid modes.
161 tests pass; new contracts cover failure preservation and prelaunch gates.
Evidence: `validation/native_cbs_stage.json`. This completes the single-sample
CBS execution component, not paired-caller validation or study calibration.

### Explicit CBS runtime fingerprint

`execute` accepts a restricted optional R runtime environment, includes it in
the task fingerprint, and passes it only to the child. Legacy task fingerprints
remain unchanged when no override is supplied. Contract tests reject unsafe
variable keys and changed-runtime resume; parent environment remains unchanged.
The native CNVkit/CBS fixture passes both coverage modes with explicit Rscript,
startup profiles disabled and a before/after hash snapshot of the R library.
157 tests pass. This is not full shared-library closure or cohort validation.
Evidence: `validation/native_cbs_runtime.json`.

### Native CNVkit normalization validation, 2026-09-16

With user authorization, CNVkit 0.9.12 and pysam 0.23.3 were restored to the
isolated Python environment. The native smoke now invokes `normalize_sample`
instead of calling fix directly. Both target-only (100 bins) and hybrid
(200 bins) passed actual coverage, receipt-bound PoN construction, verified
normalization and DNAcopy CBS. Known synthetic loss/neutral medians are -1/0.
155 tests pass with all current dependencies. This supersedes the earlier
new-normalization-native-test-pending note. It does not validate X/Y, a
study-suitable PoN, low-clonal-fraction sensitivity, or paired-caller concordance.
See `validation/native_cnvkit_normalization.json` and the dated Python freeze.

### CNVkit reference-to-normalization integration

Added read-only `verify_reference`, reconstructing expected PoN coverage
membership and verifying reference command/input/output receipts. Added
`normalize_sample` to bind one study/calibration coverage registry to that
verified reference before native fix. Four normalization and three reference
tests reject changed references, rebound registries and PoN-as-study misuse.
155 unit/contract tests pass with PyArrow enabled. The new normalization tests
mock the native fix process; they do not establish native integration or a
sex-chromosome reference policy. No cohort launched.

### Full synthetic preparation chain through native FACETS

`tests/native_joint_chain_smoke.py` now runs three synthetic BAMs through
actual snp-pileup, a two-normal synthetic training panel, pair preparation,
joint preparation and `run_joint`, with no mocked execution or verification.
The enriched fixture passes: 601 catalog loci, 600 numeric input rows,
599 native signal/segmented loci, two segments. One absent locus stays null;
one masked locus retains its audit evidence but is not used in segmentation.
The held-out manifest row is structural only and is not measured. Reference
sequence and parameter choices are synthetic, not a study resource release.
See `validation/native_joint_chain.json`. This supersedes the earlier
preparation-chain-not-tested limitation for this synthetic FACETS branch only.

### Real file-boundary FACETS test

The operational R entrypoint was tested with site/user libraries excluded,
after copying 29 existing dependency packages into the isolated library
(originals unchanged). `tests/native_joint_files_smoke.py` passed actual TSV
input, hg38 GC correction, native segmentation/export, exact input-locus checks,
task receipt rehashing and reuse of an identical completed task. The synthetic
600-locus fixture yielded two segments. Evidence is recorded in
`validation/native_joint_files.json`. This exercises native task execution,
not the full prepare-pair/train-panel/prepare-joint chain or study calibration.

### 2026-09-16 authorized isolated dependency refresh

Installed PyArrow 25.0.1 and the pinned FACETS-CH/pctGCdata with R dependencies
under `/tmp/cnv-native-deps.An66yG`; no system upgrade or server action.
All 148 Python tests pass, including the three previously skipped Parquet tests.
Actual snp-pileup synthetic paired/unpaired tests, DNAcopy two-segment smoke,
and `tests/native_facets_smoke.R` pass. Native FACETS exposed the pctGCdata
search-path requirement; the adapter now attaches that package temporarily.
The temporary FACETS source has one packaging repair (DESCRIPTION package name).
See `validation/native_dependencies_authorized.json` for versions and hashes.
This supersedes older missing-runtime notes below, not the pending study
calibration, full dependency freeze, or end-to-end pipeline validation.

- Synthetic unit tests cover identity/PoN separation, assay matching, coordinate
  boundaries, native command construction, failed outputs and resume integrity.
- pysam was exercised on synthetic BAM headers and reference/SNP allele checks.
- FACETS-CH upstream interfaces and several GRCh38 adaptation issues were read
  directly at the commit recorded in UPSTREAM_AUDIT.md.

## Native verification completed

- Historical installation session 31999 exited successfully with CNVkit 0.9.12
  and pysam 0.23.3 in `/tmp/cnv-native-testenv`. That temporary environment was
  subsequently lost; these results do not establish current runtime availability.
- Native coverage/reference/fix smoke passed on 20 generated target bins.
  Study depths were 2 for five bins and 4 for fifteen bins; normalized median
  log2 was -1.0 and 0.0, respectively. This is a gross synthetic loss, not
  a measurement of low-cell-fraction mosaic sensitivity.
- Native execution exposed and fixed an incorrect `fix -r reference`
  argument. CNVkit expects its reference as the third positional argument.
  A regression test now protects this interface; 61 unit/contract tests pass.
- The observed Python dependency set is recorded in
  `requirements-native-smoke.lock.txt`; full environment recreation and
  FACETS-CH execution are not yet verified.
- Native off-target coverage/reference/fix passed with 40 bins (20 targets,
  20 antitargets). Both modes preserved loss log2=-1 and neutral log2=0.
- User explicitly approved isolated DNAcopy installation. R 4.6.1 with
  DNAcopy 1.86.0, BiocVersion 3.23.1 and BiocManager 1.30.27 is available
  through `R_LIBS_USER=/tmp/cnv-r-library.0Hs1Qb`. No system R changes.
- Native CBS completed but did not split the initial 20-target fixture at
  alpha=0.0001 (single log2=-0.25 segment). A denser 100-target fixture
  passed the loss/neutral separation checks at the SAME alpha, both target-only
  (100 bins) and with antitargets (200 bins). This does not establish clinical
  marker-density thresholds or low-fraction sensitivity.
- Reproduce with `PATH=/tmp/cnv-native-testenv/bin:$PATH
  R_LIBS_USER=/tmp/cnv-r-library.0Hs1Qb
  /tmp/cnv-native-testenv/bin/python tests/native_cnvkit_smoke.py --cbs`
  as one shell command.
- Reference-build validation now permits missing derived marker/bias resources;
  analysis validation still requires them. Role separation remains enforced.
- cnv_reference.py now builds the native reference using only exact stratum
  PoN membership from the manifest and hash-verified coverage task receipts.
  It validates BAM/index/BED bindings, caller filter policy, unchanged output,
  exact bin coordinates, target/antitarget non-overlap and explicit minimum
  PoN sample count. Study/calibration registry entries are rejected.
- The actual native smoke now exercises this manifest/registry integration
  (not the pooled-reference command alone). Both 100-target and 200 combined
  bin coverage/reference/fix/CBS runs passed with loss=-1 and neutral=0.
  One synthetic normal is used ONLY to test the interface.
- Completed reading FACETS-CH utils, procreads and PoN notebook at the pinned
  commit. Bias uses log2 without pseudocount in the notebook, but sample signal
  uses natural log with 1/6 pseudocount. Recorded both differences and other
  marker/window/coordinate issues in UPSTREAM_AUDIT.md.
- cnv_allelic.py now provides per-marker training summaries, explicit log units,
  declared pseudocount and bias subtraction; missing/zero/insufficient hets are
  separate. New tests verify signed conversion, estimator consistency, PoN-only
  IDs and unavailable correction. No WES thresholds are declared calibrated.
- cnv_pileup.py now streams native single-sample CSV/BGZF against a frozen
  SNP/pseudo-depth catalog, retains all-missing loci, distinguishes native
  zero ref+alt rows with errors, and rejects allele swaps/duplicates/foreign
  IDs/off-catalog loci. Training summaries preserve missingness and depth-only
  semantics.
- cnv_catalog.py streams a frozen whole-reference depth grid merged with
  biallelic SNPs, keeping SNP precedence at grid overlaps. The VCF adapter
  verifies REF against FASTA with pysam; the pure grid merger has unit tests.
  The catalog denotes expected coordinates, not evaluated territory.
- The earlier temporary Python environment disappeared after session change;
  its prior test results remain historical, not proof of current availability.
  Native snp-pileup tests now use installed samtools, without pysam.
- User authorized isolated compilation of pinned bundled snp-pileup.
  Build at /tmp/cnv-snp-pileup.ODBXqD/snp-pileup succeeded against HTSlib 1.19.
  The native smoke passed unpaired and fully-overlapping proper-pair fixtures,
  with exact ref/alt/error counts and no overlapping-fragment double count.
  Provenance is in validation/native_pileup_build.json. Revalidate /tmp paths
  on future sessions. Full FACETS-CH R execution remains untested.
- facets_grch38.R implements the native primitive composition with hg38,
  explicit chromosome mapping/masks, bias correction before segsnps, and
  exclusion rows. Native fixed 1/6 pseudocount and sample-het depth 50 are
  asserted. R contract tests pass with an injected fake engine ONLY.
  The current environment has R 4.6.1 but no facetsCH installed.
- Verified segsnps signature uses delta, not deltaCN. Its scaling divides by
  mad(diff(logR)); zero-dispersion and no-allele/hetscale inputs now fail
  explicitly pending separate CNV-only review, never generate negative calls.
- The adapter now emits a full locus_audit retaining original coordinates,
  input row linkage, sample/normal depth and distinct corrected/depth-only/
  bias-unavailable/not-evaluated states. R contract tests verify the ledger.
- Prepared tests/native_facets_smoke.R for actual package execution, still
  UNRUN. Installation plan is validation/facets_install_plan.json. User was
  asked to authorize isolated FACETS-CH/pctGCdata/dependency installation;
  no response is recorded in this file. Recheck conversation before installing.
- facets_segments.R summarizes final native jointseg into explicitly labelled
  1-based native coordinates and 0-based half-open marker-bounding intervals.
  Singleton span is one base; missing allelic evidence remains NA/CNV_ONLY.
  Evaluable bp stays NA until territory/calibration integration; event classes
  remain UNCLASSIFIED. Tests include gaps, singletons and invalid segment IDs.
  The unrun native FACETS smoke now includes this export conversion.
- cnv_territory.py provides indexed exact half-open overlap accounting for
  assay, coverage-eligible and externally supplied calibrated territories.
  Missing masks remain None, known empty masks yield zero overlap. Geometric
  overlap alone never marks calibration verified or fills evaluable_bp.
  Unit tests compare all small intervals with a basewise oracle and test
  actual R-summary TSV to Python-territory integration across a capture gap.
- cnv_reference_selection.py ranks every eligible same-stratum/reference-sex
  PoN on the exact common finite locus set, with explicit minimum common count/
  fraction, deterministic tie handling and common-locus hash. Unknown sex,
  candidate omission, undeclared loci and inadequate intersection are rejected.
  Study/calibration samples cannot become references. Tests verify that
  differing missingness does not unfairly favor a reference. Signal receipt
  verification and orchestration remain to be connected.
- cnv_pair_matrix.py converts exact sample/PoN matched pileups to R count
  rows while retaining the complete absence/error/mask ledger. It rejects
  self-reference, calibration references, sex mismatches and foreign inputs.
  Native complete numeric rows retain original catalog IDs; facets_grch38.R
  now propagates those IDs into signal and locus audits. Upstream native
  receipt verification is still an orchestration prerequisite, not claimed
  by this conversion function.
- verify_completed_task now checks STARTED/COMPLETE fingerprints, expected
  command/filter policy, manifest input bindings, current implementation
  hashes and output schemas/hashes. Input-byte rehashing is an explicit option,
  separately reported; metadata checks do not imply rereading BAMs.
- verified_pair_matrix integrates those checks for both sample and reference
  before conversion. Tests reject changed commands, manifest bindings, status
  and output bytes. Native snp-pileup smoke was rerun successfully with receipt
  verification and input-byte rehashing for both unpaired/overlapping-pair cases.

## Next implementation work

Latest local integration evidence:
- `cnv_burden.py` adds descriptive aggregation of externally reviewed event IDs,
  separately from candidate segments. Multiple intervals per event count once;
  altered bases use geometric union within the supplied territory. Missing/empty
  territory, incomplete review and HOLD events do not produce complete zero
  endpoints. Accepted CNLOH requires an explicit external copy-neutral review
  declaration; the function itself does not validate that declaration or authorize
  inference. Current suite: 148 discovered, 145 passed, 3 PyArrow tests skipped.
  Review receipt validation, event-specific calibrated territory and model-ready
  endpoints remain outstanding.
- `cnv_parquet.py` implements candidate Parquet plus canonical JSONL, explicit
  nullable schema, batching, source-byte checks and full field round-trip digest
  verification. Its three native integration tests cover multiple batches,
  empty typed output and duplicate-record failure. They are currently SKIPPED
  because PyArrow is unavailable: 142 tests discovered, 139 passed, 3 skipped.
  Do not report Parquet as validated until those tests run with real PyArrow.
  This export does not verify evidence derivation, de-identify payloads or replace
  adjudicated-event exports. Isolated dependency installation remains pending.
- `cnv_export_contract.py` defines typed candidate identity/coordinate fields and
  canonical full-evidence JSON, preserving nulls, zeros and nested source data.
  Duplicate identities, nonfinite values, partial coordinates and premature
  biological/denominator claims are rejected. Current suite: 139 tests passing.
  This is the candidate-record contract ONLY, not implemented/validated Parquet
  serialization or the separate adjudicated-event export. PyArrow is absent from
  current Linux Python; bundled Windows Python could not start through WSL.
  Isolated PyArrow installation approval has been requested; none was attempted.
- Joint execution now compares each audited catalog ID, chromosome, position and
  sample/reference depth against the submitted numeric matrix, not only its row
  count. Tests reject changed loci and depths. Standalone structural validation
  without a matrix explicitly reports submitted_matrix_rows_verified=false;
  the execution path always supplies the matrix. Current suite: 133 tests passing.
- `cnv_cli.py` exposes validation, pair/panel preparation, joint preparation,
  isolated-library snapshotting and single-fit execution as separate commands.
  Native execution requires explicit --execute and timeout. Preparation does not
  chain to R or launch a cohort. Current suite: 132 tests passing, including CLI
  execution gating and preparation routing. This is a development interface,
  not a production-ready workflow or native FACETS validation.
- `cnv_manta.py` retains single-sample exploratory SV records, original fields,
  separate PR/SR counts, missing evidence, IMPRECISE, MATEID and EVENT. BND pairs
  remain records, not independent event counts; PASS never assigns somatic origin
  or CNV class. Exact VCF sample names are required separately from study IDs.
  Current suite: 128 tests passing. This is synthetic VCF parsing only, not native
  Manta execution or source-provenance validation.
- `cnv_joint_output.py` now gates the final joint receipt on tabular structure,
  complete numeric-input audit, unique source-catalog linkage, finite signals,
  segment/count bounds and absence of premature biological claims. File presence
  and COMPLETE.txt alone no longer suffice. Current suite: 123 tests passing;
  malformed audit/linkage and unsupported CNLOH labels are rejected. Validation
  is explicitly structural only. Real FACETS output remains to be exercised after
  the separately requested isolated installation is authorized.
- `cnv_run_joint.py` now links joint preparation, explicit R executable/library
  fingerprints and the native task runner. It verifies preparation/runtime before
  and after execution, validates the native completion text and emits a separate
  final receipt only on success. Failed/incomplete attempts are preserved. Current
  suite: 119 tests; new orchestration tests MOCK native execution and do not prove
  actual FACETS operation. System shared-library dependency closure, native FACETS
  execution and scientific calibration remain unverified. Nothing was launched.
- `cnv_regions.py` maps normalized segments to explicitly supplied named regions
  (gene/exon/arm definitions), retaining each original segment and unioning
  overlapping geometry. Assay overlap, missing versus zero coverage eligibility,
  absent segment evidence and outside-assay regions remain distinct. Overlapping
  gene/arm summaries are not independent territory and do not create event counts.
  No built-in TP53/17p coordinates or negative calls are inferred. Current suite:
  116 tests passing, using synthetic region coordinates. Actual frozen regional
  resources, source-receipt wiring and calibrated regional interpretation remain
  outstanding.
- Runtime recheck confirms facetsCH, pctGCdata, data.table, argparse, R.utils
  and DNAcopy remain absent from the current vanilla R 4.6.1 environment.
  No installation was attempted; additional FACETS installation approval is
  pending. `cnv_r_library.py` fingerprints all isolated-library files, detects
  additions/changes, rejects installation locks and external/directory symlinks.
  This is file stability only, not native loading or dependency closure. Current
  suite: 110 tests passing. Native FACETS execution remains unverified.
- `cnv_phase_verified.py` connects phase evidence to the role/assay manifest,
  explicit BAM SM, reference FAI and both native extraction receipts. It verifies
  BAM/index/FASTA bindings, exact requested loci, REF alleles and the SAM/FAI
  sequence dictionary. The real synthetic BAM test now covers this entire path
  and rejects a different locus, changed BAM binding, incorrect SM and dictionary
  mismatch. 105 tests pass with expanded coverage. Provenance is explicitly
  limited to declared manifest/native-task bindings, not clinical identity or
  calibrated phase confidence. No cohort execution was performed.
- Phase extraction now requires separate study ID and expected BAM SM. Native
  tests use different identifiers and reject a mismatched SM, without relabeling
  BAM data or guessing identity from filenames. The SNP reference-check path
  now runs real samtools faidx on a synthetic indexed FASTA and rejects REF
  mismatches/missing loci. This validates alleles against the supplied reference,
  not its GRCh38 provenance. The existing 105-test suite passes with these
  expanded native checks. Manifest-level linkage remains to be integrated.
- `cnv_phase_sam.py` adds indexed samtools multi-region SNP extraction, CIGAR
  coordinate mapping and exact read-group/SM checks, followed by template-grouped
  observations. A real synthetic BAM test confirms that two overlapping mates
  spanning both loci yield two SAM records, four base observations and one linked
  ALT_ALT fragment. CIGAR clipping/insertion/deletion offsets are tested. Current
  suite: 105 tests. Existing samtools 1.19.2 was reused; no package installation.
  Only biallelic SNP pairs are supported; complex alleles are rejected. FASTA REF
  validation, study manifest receipt orchestration, UMI handling and calibrated
  phase confidence remain incomplete. Synthetic identifiers are not patient data.
- `cnv_phasing.py` aggregates two-locus observations by explicit sample,
  read-group and template name. Overlapping mates count once; conflicting or
  non-REF/ALT alleles remain ambiguous, never majority-voted. Read/base quality
  and flag policy is explicit; unavailable MAPQ=255 is excluded. No linking
  template means PHASE_NOT_EVALUATED, not trans. Four haplotype counts are raw
  fragment evidence, not calibrated phase or a somatic second hit. This does not
  yet extract observations from BAM or collapse UMIs; provenance/phase-confidence
  calibration remain separate requirements. Current suite: 102 tests passing.
- `cnv_reconcile.py` reads native CNVkit CNS and adapted FACETS TSV and splits
  their interval union at every caller boundary. Original source fields and
  segment signals remain attached; marker counts are not redistributed to the
  new subintervals. Explicit direction thresholds distinguish discordance,
  concordance, single-method evidence and incomparable missing signal. The
  within-logR-band state never establishes copy neutrality or cnLOH. Gaps,
  absence and evaluability are not filled with zero. Current suite: 94 tests,
  including a basewise geometry oracle. Threshold calibration, source-receipt
  orchestration and reviewed biological classification remain outstanding.
- `cnv_prepare_joint.py` binds the verified pair/panel to an explicit fit policy
  and writes R-readable config/map TSVs. The estimator and chromosome mapping
  cannot be overridden by fit policy. No-row and no-bias/het-scaling cases stop
  explicitly; the latter requires CNV-only review rather than manufactured bias.
  Preparation hashes Python/R implementations and inputs, preserves existing
  attempts and returns JOINT_FIT_PREPARED_NOT_LAUNCHED. Current suite: 87 passing
  tests. Native execution/runtime verification and calibrated review are still
  separate, unfinished requirements.
- `facets_run.R` and `cnv_commands.joint_fit` define the current native adapted
  invocation, with absolute input/output paths, R --vanilla, explicit isolated
  library, script-relative source loading and no legacy implicit best-normal
  selection. All seven result artifacts are declared for outer hashing. The
  output is a fresh child directory, separate from the task runner's receipts.
  Current suite: 84 tests passing, including CLI help outside the project cwd
  and rejection of legacy flags. Full native invocation remains untested until
  FACETS-CH installation is authorized and completed. The command constructor
  does not itself verify the prepared artifacts or hash installed dependencies;
  those orchestration requirements remain open.
- `facets_io.R` provides typed TSV input, data-only configuration parsing and
  native-result/segment/audit export. Empty bias tables preserve numeric column
  types, failures do not write COMPLETE, and existing attempts are preserved.
  The R file-boundary test uses an injected fake engine; it is NOT native FACETS
  validation. Current suite: 80 tests passing. TSV round-trip testing exposed
  decimal rounding of 1/6; bias pseudocount comparison now permits only eight
  double-machine epsilons of rounding, rather than exact floating equality.
  Provenance verification still belongs to the Python orchestrator; the R text
  completion marker alone is not a verified pipeline receipt.
- `cnv_joint_inputs.py` verifies prepared pair/panel receipt bindings, output
  hashes, current implementation, matching manifest/config, stratum/reference sex,
  catalog/mask/filter policy and the exact reference pileup task. Prepared stages
  now bind COMPLETE to the STARTED receipt hash; older unbound prepared artifacts
  are intentionally rejected. This verifies artifact compatibility, not joint-fit
  success, clinical identity or calibrated sensitivity; it does not rehash BAMs.
  The current suite passes 79 tests. Native pileup/pair/panel smoke remains passing;
  the new cross-stage compatibility tests use synthetic mocked stage computations.
- `cnv_training.py` now verifies native pileup receipts and manifest input
  bindings for exactly all same-stratum/resolved-sex PoN members before streaming
  marker summaries. Study/held-out inclusion, omitted eligible normals, unknown
  PoN sex and failed receipts are rejected. Masks retain raw summaries but disable
  correction; pseudo-depth rows cannot acquire allelic correction.
- The current suite passes 70 tests. Native snp-pileup smoke also exercises this
  PoN training path with two generated BAMs in a separate synthetic-role manifest,
  for both unpaired and overlapping paired reads. This is an interface/count test,
  not a clinically suitable PoN or calibration of the synthetic fixture thresholds.
  Full FACETS integration remains needed.
- `cnv_train_panel.py` now writes `training.audit.jsonl.gz`, `stable_markers.tsv`
  and `allelic_bias.tsv`, with new-attempt-only output and STARTED/FAILURE/COMPLETE
  receipts. Native bias indices use an explicit numeric chromosome map (X=23);
  Y and unmapped loci remain in the audit rather than entering the native bias
  table. Empty stable-marker panels are explicitly reported, not fabricated.
  Completion means successful training artifact construction, not sensitivity
  calibration. The integrated native smoke verifies the output hashes and
  rejection of overwrite attempts. The current suite passes 74 tests, including
  interrupted streams, changed inputs, invalid maps and zero stable markers.
- `cnv_prepare_pair.py` now writes a full gzip JSONL evidence ledger and a
  complete-numeric-row R matrix, with STARTED/COMPLETE/FAILURE receipts and
  pre/post input/implementation hashes. Existing attempts cannot be overwritten.
  An empty mask is permitted explicitly; an empty target remains invalid.
- The pure local suite passes 65 tests. New tests cover observed zero versus
  missing data, serialization, midstream failure, input mutation and preservation
  of existing attempts. These stage unit tests mock native verification.
- Separately, `tests/native_pileup_smoke.py` runs the real compiled snp-pileup
  for distinct synthetic sample/reference BAMs and calls pair preparation with
  real native receipts and upstream-byte rehashing. Both unpaired and overlapping
  paired cases pass: 20 catalog loci, 4 complete numeric rows, 16 incomplete rows.
  Masked error-only zero depth is retained and absent loci are excluded only from
  the numeric matrix, never converted to zero in the ledger. The held-out manifest
  row is structural only and does not represent executed calibration.
- Pair preparation does not execute FACETS-CH or establish that a joint fit is
  possible. Native joint fitting, calibration and biological classification remain
  unverified; no clinical/sample identity claim follows from these synthetic BAMs.

1. Extend the manifest/receipt integration from reference building to coverage
   preflight and study fix/CBS orchestration, including explicit sex/PAR policy.
2. Implement the actual audited GRCh38 FACETS-CH adaptation and native tests;
   the command constructor alone does not provide that adaptation.
3. Wire reference construction, marker/allelic-bias derivation and preflight to
   role/assay manifests, then stage orchestration and calibrated review outputs.
4. Complete TP53/17p, Manta exploratory evidence, Parquet/export, calibration and
   models against all seven requirements in README.md before claiming completion.

There has been no CNV cohort launch or change to the running SNV pipeline.
