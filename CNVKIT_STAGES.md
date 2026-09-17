# Explicit CNVkit stage interface

Development interface; not authorization to run the study. All paths below are
placeholders. Use absolute paths in real configurations. There is no automatic
cohort launch, segmentation, classification or downstream analysis.

## Preconditions

The manifest must pass `reference_build` validation. Native coverage tasks must
already exist with hash-valid STARTED/COMPLETE receipts. The PoN registry has
exactly the stratum's PoN members; study/calibration samples never enter it.
Coverage registries map a sample ID to `target` and, only when the manifest
declares supported off-target coverage, `antitarget` task directories.

The coverage process count and mapping-quality cutoff must reproduce the
original native coverage command. Minimum PoN sample count is an explicit
study-design input; the software does not choose an adequate clinical PoN size.

## Reference configuration

Exact JSON keys for `run-reference --config reference.json`:

```json
{
  "manifest_path": "/absolute/manifest.json",
  "registry_path": "/absolute/pon-coverage-registry.json",
  "stratum_id": "ASSAY_STRATUM",
  "executable": "/absolute/isolated/bin/cnvkit.py",
  "processes": 1,
  "min_mapq": 20,
  "minimum_pon_samples": 2
}
```

Numbers illustrate the schema only, not recommended study parameters.

```bash
python cnv_cli.py run-reference --config /absolute/reference.json \
  --output /absolute/new-reference-task \
  --implementation-file /absolute/isolated/bin/cnvkit.py \
  --timeout-seconds 600 --execute
```

`--implementation-file` is repeatable. Include the executable and relevant
runtime implementation files. Hashing an entrypoint alone does not freeze the
Python environment, native dependencies or reference suitability.

## Normalization configuration

For `run-normalize`, replace `registry_path` and `stratum_id` above with:

```json
{
  "sample_id": "STUDY_SAMPLE",
  "sample_registry_path": "/absolute/sample-coverage-registry.json",
  "pon_registry_path": "/absolute/pon-coverage-registry.json",
  "reference_directory": "/absolute/completed-reference-task"
}
```

Keep `manifest_path`, `executable`, `processes`, `min_mapq` and
`minimum_pon_samples`. Pass the same command-level flags as reference, using
`run-normalize` and a separate output directory. Only a declared study or
calibration sample is accepted. This stage does not automatically invoke CBS.

## CBS configuration and runtime

The CBS `--config` has the same exact keys as normalization plus `directory`,
pointing to the completed normalization task (not the destination CBS task).

Create the isolated R library separately. Then freeze it:

```bash
python cnv_cli.py snapshot-r-library --library /absolute/isolated/R-library \
  --require-package DNAcopy --output /absolute/new-r-library-snapshot.json
```

The `--rscript-artifact` JSON has `path` and `sha256` of the actual Rscript
executable. Both are verified. CBS process count and threshold are separate
from the coverage policy in the normalization verification config.

```bash
python cnv_cli.py run-cbs --config /absolute/cbs-normalization.json \
  --library-snapshot /absolute/r-library-snapshot.json \
  --rscript-artifact /absolute/rscript-artifact.json \
  --output /absolute/new-cbs-attempt \
  --implementation-file /absolute/isolated/bin/cnvkit.py \
  --processes 1 --threshold 0.0001 --timeout-seconds 600 --execute
```

The threshold shown is a synthetic-test value, not a calibrated production
setting. CBS requires a fresh attempt directory; incomplete attempts are
preserved. It produces candidate segments, not copy-number or somatic-origin
calls. X/Y policy, study calibration, reviewed classification, evaluable burden
and paired-caller comparison are separate unfinished requirements.
