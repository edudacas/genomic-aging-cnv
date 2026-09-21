# Representative real-BAM canary commands

These commands are prepared for a deliberately small, stratified pilot. They do
not authorize a 1,329-sample launch and must be run only after the private
resource manifest and independent PoN references pass hash validation.

Select de-identified canaries covering Twist and IDT, low/middle/high depth,
both inferred reference-sex groups where evaluable, and at least one known
technical edge case. Do not select cases by CNV outcome.

## 1. Bind paths and fail closed

```bash
set -euo pipefail
CNV_REPO=/home/eduardo/06_Projeto_de_pesquisa/Phase0_DataFoundation/06_VARIANT_ANNOTATION_LAYER/output/TP53_parquet/projects/Genomic_aging_LFS/github_genomic_aging_cnv
CNV_RUNTIME=/home/eduardo/06_Projeto_de_pesquisa/Phase0_DataFoundation/06_VARIANT_ANNOTATION_LAYER/output/TP53_parquet/projects/Genomic_aging_LFS/30_CNV_NATIVE_RUNTIME_v1
CNV_PRIVATE_ROOT=/absolute/private/cnv_canary_20260921
CNV_MANIFEST=/absolute/private/cnv_resources/frozen_manifest.json
CNV_PY="${CNV_RUNTIME}/python/bin/python"
CNVKIT="${CNV_RUNTIME}/python/bin/cnvkit.py"
RSCRIPT=/usr/bin/Rscript

test -d "${CNV_REPO}"
test -x "${CNV_PY}"
test -x "${CNVKIT}"
test -x "${RSCRIPT}"
test -f "${CNV_MANIFEST}"
test ! -e "${CNV_PRIVATE_ROOT}"
mkdir -m 700 "${CNV_PRIVATE_ROOT}"
cd "${CNV_REPO}"

"${CNV_PY}" cnv_cli.py validate-manifest \
  --manifest "${CNV_MANIFEST}" \
  --validation-stage analysis \
  --check-files
sha256sum --quiet -c SHA256SUMS
```

Replace every `/absolute/private/...` value only with a reviewed private path.
The guard intentionally refuses an existing output root.

## 2. Verify the frozen container and resource receipt

```bash
CNV_IMAGE_REF=/absolute/private/container/IMAGE_REFERENCE.txt
CNV_RESOURCE_RECEIPT=/absolute/private/cnv_resources/COMPLETE.json

test -s "${CNV_IMAGE_REF}"
test -s "${CNV_RESOURCE_RECEIPT}"
docker image inspect "$(cat "${CNV_IMAGE_REF}")" >/dev/null
"${CNV_PY}" -m json.tool "${CNV_RESOURCE_RECEIPT}" >/dev/null
```

The resource receipt must say that GRCh38 FASTA/FAI, dense SNP catalog, masks,
capture-specific target/antitarget BEDs, stable-marker tables, allelic-bias
tables and independent PoN membership are frozen for both assay families.

## 3. Build one CNVkit reference per frozen assay stratum

Create reviewed JSON configurations matching `CNVKIT_STAGES.md` and run:

```bash
for ASSAY in TWIST_WES_OS_T7 MGI_IDT_EXOME; do
  CONFIG="/absolute/private/cnv_resources/${ASSAY}.reference.json"
  OUTPUT="${CNV_PRIVATE_ROOT}/reference_${ASSAY}"
  test -f "${CONFIG}"
  test ! -e "${OUTPUT}"
  "${CNV_PY}" cnv_cli.py run-reference \
    --config "${CONFIG}" \
    --output "${OUTPUT}" \
    --implementation-file "${CNVKIT}" \
    --timeout-seconds 14400 \
    --execute
done
```

PoN registries must contain exactly the independent external donors for that
stratum. Study and held-out calibration samples must not appear in them.

## 4. Run CNVkit normalization and CBS for each canary

For every canary, prepare exact normalization and CBS verification JSON files.
Use a new output directory for every attempt.

```bash
CANARY_ID=REPLACE_WITH_PSEUDONYMOUS_SAMPLE
NORMALIZE_CONFIG="/absolute/private/cnv_canaries/${CANARY_ID}.normalize.json"
CBS_CONFIG="/absolute/private/cnv_canaries/${CANARY_ID}.cbs.json"
R_LIBRARY_SNAPSHOT=/absolute/private/cnv_resources/R_LIBRARY_SNAPSHOT.json
RSCRIPT_ARTIFACT=/absolute/private/cnv_resources/RSCRIPT_ARTIFACT.json

case "${CANARY_ID}" in
  REPLACE_*|"") echo "Set a reviewed pseudonymous CANARY_ID" >&2; exit 2 ;;
esac

"${CNV_PY}" cnv_cli.py run-normalize \
  --config "${NORMALIZE_CONFIG}" \
  --output "${CNV_PRIVATE_ROOT}/${CANARY_ID}/cnvkit_normalize" \
  --implementation-file "${CNVKIT}" \
  --timeout-seconds 7200 \
  --execute

"${CNV_PY}" cnv_cli.py run-cbs \
  --config "${CBS_CONFIG}" \
  --library-snapshot "${R_LIBRARY_SNAPSHOT}" \
  --rscript-artifact "${RSCRIPT_ARTIFACT}" \
  --output "${CNV_PRIVATE_ROOT}/${CANARY_ID}/cnvkit_cbs" \
  --implementation-file "${CNVKIT}" \
  --processes 1 \
  --threshold 0.0001 \
  --timeout-seconds 7200 \
  --execute
```

The CBS threshold above reproduces the validated smoke-test interface; it is not
a production event threshold. Freeze the pilot policy before interpreting
segments.

## 5. Prepare and run the allelic/FACETS side

```bash
PAIR_CONFIG="/absolute/private/cnv_canaries/${CANARY_ID}.pair.json"
PANEL_CONFIG=/absolute/private/cnv_resources/PON_PANEL_CONFIG.json
PAIR_DIR="${CNV_PRIVATE_ROOT}/${CANARY_ID}/pair"
PANEL_DIR=/absolute/private/cnv_resources/PON_PANEL
FIT_POLICY=/absolute/private/cnv_resources/FIT_POLICY.json
PREPARED_DIR="${CNV_PRIVATE_ROOT}/${CANARY_ID}/joint_prepared"
JOINT_DIR="${CNV_PRIVATE_ROOT}/${CANARY_ID}/joint_native"

"${CNV_PY}" cnv_cli.py prepare-pair \
  --manifest "${CNV_MANIFEST}" \
  --config "${PAIR_CONFIG}" \
  --output "${PAIR_DIR}"

"${CNV_PY}" cnv_cli.py prepare-joint \
  --manifest "${CNV_MANIFEST}" \
  --pair-directory "${PAIR_DIR}" \
  --pair-config "${PAIR_CONFIG}" \
  --panel-directory "${PANEL_DIR}" \
  --panel-config "${PANEL_CONFIG}" \
  --fit-policy "${FIT_POLICY}" \
  --output "${PREPARED_DIR}"

"${CNV_PY}" cnv_cli.py run-joint \
  --prepared "${PREPARED_DIR}" \
  --library-snapshot /absolute/private/cnv_resources/FACETS_R_LIBRARY_SNAPSHOT.json \
  --rscript-artifact "${RSCRIPT_ARTIFACT}" \
  --output "${JOINT_DIR}" \
  --timeout-seconds 14400 \
  --execute
```

## 6. Reconcile callers and export evidence

```bash
"${CNV_PY}" cnv_cli.py compare-callers \
  --manifest-artifact "${CNV_MANIFEST}" \
  --sample-id "${CANARY_ID}" \
  --cbs-artifact "${CNV_PRIVATE_ROOT}/${CANARY_ID}/cnvkit_cbs/COMPLETE.json" \
  --joint-artifact "${JOINT_DIR}/COMPLETE.json" \
  --output "${CNV_PRIVATE_ROOT}/${CANARY_ID}/comparison" \
  --logr-direction-threshold 0.05

"${CNV_PY}" cnv_cli.py export-comparison \
  --receipt-artifact "${CNV_PRIVATE_ROOT}/${CANARY_ID}/comparison/COMPLETE.json" \
  --sample-id "${CANARY_ID}" \
  --output "${CNV_PRIVATE_ROOT}/${CANARY_ID}/export" \
  --batch-size 1024
```

The numerical thresholds are pilot configuration inputs, not biological truth.
Any caller-only, discordant or unevaluable evidence must remain explicit.

## 7. Sex-chromosome BAM-level validation

Run the already frozen target-read calibration independently for mLOY and mLOX:

```bash
"${CNV_PY}" sex_chromosome_read_calibration.py \
  --axis mloy \
  --collection /absolute/private/mloy-y-intervals.GRCh38.v1.PRIVATE.json.gz \
  --metrics /absolute/private/mloy.metrics.PRIVATE.json \
  --panel /absolute/private/mloy.callable_panel.json \
  --donor-plan /absolute/private/mloy.calibration.plan.PRIVATE.json \
  --bam-manifest /absolute/private/mloy.read-calibration.bams.PRIVATE.json \
  --output-dir "${CNV_PRIVATE_ROOT}/mloy_read_calibration"

"${CNV_PY}" sex_chromosome_read_calibration.py \
  --axis mlox \
  --collection /absolute/private/mlox-x-intervals.GRCh38.v1.PRIVATE.json.gz \
  --metrics /absolute/private/mlox.metrics.PRIVATE.json \
  --panel /absolute/private/mlox.callable_panel.json \
  --donor-plan /absolute/private/mlox.calibration.plan.PRIVATE.json \
  --bam-manifest /absolute/private/mlox.read-calibration.bams.PRIVATE.json \
  --output-dir "${CNV_PRIVATE_ROOT}/mlox_read_calibration"
```

Fraction-zero reconstruction must be exact. A mismatch blocks the sample.

## Stop/go criteria

Do not scale beyond the canary set unless all attempts have hash-valid
`COMPLETE.json` receipts, no silent sample drops, expected GRCh38 contigs,
stable resource hashes, no PoN/study overlap, bounded resource use and reviewed
CNVkit/FACETS concordance. Candidate segments remain research evidence, not
confirmed somatic events or TP53 second hits.
