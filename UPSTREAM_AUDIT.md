# FACETS-CH adaptation audit

Inspected upstream commit `696113be3e31dedde130d46be0c9b2b89c41499b` at
https://github.com/papaemmelab/Gao_NC_CH.
Files read: `facets-ch/README.md`, `facets-ch/facets-ch.R`,
`facets-ch/R/facets-wrapper.R`, complete `facets-ch/R/facets-ch-utils.R`,
`facets-ch/R/facets-procreads.R` and all source cells of
`facets-ch/create_pon_files.ipynb`.

## Findings requiring implementation or calibration

1. `preProcSample` supports `hg38`, but its first/default genome is `hg19`.
   `get_logr` does not supply `gbuild`. GRCh38 must be passed explicitly and its
   GC resources checked against reference contigs/coordinates.
2. `get_logr` contains chromosome-specific coordinate filters on chromosomes
   2, 6 and 7. Do not reinterpret these coordinates as GRCh38; replace them with
   a documented, hash-bound mask derived for this assay/build.
3. `facets-ch.R` chooses Female PoN samples for every sex value outside its male
   list. Unknown/discordant sex must be held explicitly; the adapter now rejects
   unknown rather than invoking that branch.
4. `get_logr` fixes `ndepth=50`, `ndepthmax=1500`, `het.thresh=0.1`, initial
   `cval=25`, and random seed 0. Final `cval` is a different parameter.
   These settings require explicit configuration and WES calibration.
5. Allelic correction uses `bias_median/1.5` with a source comment about log-base
   conversion. Audit the PoN notebook's actual stored units before choosing a
   conversion; do not silently modify the formula and call it equivalent.
6. `call_aberrant` applies default `p_threshold=1e-8`, requires `p_y<0.05`, and
   recurses until the number of aberrant segments stabilizes. Evaluate absent
   heterozygotes, too few normal segments, NaNs, small background windows, and
   oscillation; add bounded convergence with explicit unresolved outputs.
7. The nearest-model classifier considers loss/gain/cnLOH on a cell-fraction
   grid 0–0.90 with 0.01 spacing and a fixed residual cutoff. Preserve native
   evidence, but do not promote a nearest-model LOH label to confirmed cnLOH
   without calibrated copy-neutral support and informative markers.
8. Reference selection minimizes squared logR. Candidate reference comparison
   must use comparable marker sets. Self-reference must be excluded by exact
   identity, not the upstream test that merely compares all alternate counts.

## PoN notebook and native allele-estimator findings

- The notebook builds bias from `log2(File1A/File1R)`. Native unmatched
  `counts2logROR` instead uses
  `ln(alt + 1/6) - ln(ref + 1/6)`, with positive sign for alt enrichment.
  Thus there are TWO estimator differences: log base and pseudocount.
  Multiplying a legacy log2 bias by ln(2) fixes units only; it does not make
  that bias identical to an estimate rebuilt with the native pseudocount.
- New `cnv_allelic.py` declares ln_alt_ref units and requires an explicit
  pseudocount. Rebuilding training bias with the same estimator as sample
  signal is an adaptation requiring calibration, not upstream equivalence.
- Notebook heterozygote filters are depth >50 and <2000 and BAF (0.1,0.9);
  median coverage cutoff 125, log2 median deviation 0.4 and SD 0.4 are
  panel-specific choices, not study defaults. Its whitelist accepts n=1;
  the new summary requires a configured minimum >=2 to estimate dispersion.
- Missing normal rows must not become zero depth. The summary reports observed,
  missing, zero-depth and informative-heterozygote counts separately; absent
  or rejected correction is None, never an invented zero.
- Native unmatched heterozygosity also hardcodes sample depth >=50 inside
  procSnps, independent of the outer normal-depth threshold.
- procSnps excludes Y by default. Any X/Y extension needs explicit ploidy/
  compartment handling, not only a GRCh38 switch.
- get_seg uses SNP positions as start/end, with size=end-start; these are not
  BED coordinates or continuous evaluable territory lengths. Exports must
  retain native coordinates and explicit converted coordinate semantics.
- calc_sem_x currently averages n_marker+1 values per window via i:(i+n_marker).
  Its 1:(length-n_marker) construction is unsafe when background is too short.
  Bounds and estimator intent must be addressed in the adapted caller.

## Native snp-pileup interface

Inspected `mskcc/facets/inst/extcode/snp-pileup.cpp` on upstream HEAD (not yet a
release-pinned dependency). Arguments are VCF, output, BAM(s). `-v` is verbose,
not a VCF-valued option. `-x` disables paired-read overlap detection and is
absent from the adapter. `-A` counts anomalous pairs and requires explicit
configuration. `-P` generates pseudo-SNP depth positions. The adapter exposes
MQ, BQ, depth cap and minimum output count explicitly. Cross-sample marker
completeness and zero-depth denominator handling still need end-to-end tests.

The bundled C++ source at the pinned Gao_NC_CH commit was additionally checked
for output semantics. Pseudo-SNP rows use Ref=Alt='.' and place their passing
read count in File1R, with File1A/E/D zero. SNP output requires at least one
ref/alt/error count across files; zero rows can therefore be omitted even with
minimum read count zero. No absent row is converted to measured zero by the
streaming adapter. An error-only SNP row is observed with ref+alt depth zero,
distinct from missing. Pseudo counting does not apply the same explicit
deletion branch as SNP counting; these remain separately labelled depth
evidence pending native tests and WES calibration.

## Current execution boundary

The pinned DESCRIPTION declares Package: FACETS-CH (0.0.1, GPL >=2), while
the pinned facets-ch.R CLI loads facetsCH. The installation plan records this
packaging discrepancy; any repair must be explicitly recorded. pctGCdata is
pinned in that plan, but has not yet been installed or validated here.

Command builders and the native-task runner are implemented. The runner hashes
declared inputs and implementations, validates declared output schemas, and
rejects changed outputs on resume. CNVkit 0.9.12 coverage/reference/fix now pass
on generated synthetic BAMs, including a preserved log2=-1 coverage reduction.
This exposed the positional reference argument required by the native fix CLI.
CBS with DNAcopy 1.86.0 has now been executed: loss/neutral checks passed on
a 100-target fixture, also with 100 antitargets, at alpha=0.0001. The initial
20-target fixture did not split at that alpha. These are integration results,
not sensitivity estimates. FACETS-CH has not been executed. The GRCh38 wrapper and complete
dependency environment remain to be implemented and verified before production
caller execution is offered.
# Exploratory Manta interface audit

Source: https://github.com/Illumina/manta/blob/v1.6.0/docs/userGuide/README.md

The v1.6.0 guide distinguishes tumor-only `tumorSV.vcf.gz` from scored paired
somatic output. Unpaired calling has no quality scoring model; paired/split
support and selected filters remain available, and SR may be absent. Do not
zero-fill absent support, add PR and SR as independent fragments, or treat PASS
as a somatic classification. The exome flag changes depth filtering and remains
necessary with callRegions. These interface facts do not establish suitability
or sensitivity for low-fraction blood mosaicism. Native Manta remains untested.
