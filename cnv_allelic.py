"""Explicit-unit PoN marker summaries for the GRCh38 FACETS-CH adaptation.

These are training summaries, not genotypes or somatic calls. Thresholds must
be frozen/calibrated outside this module. None of the upstream panel-specific
thresholds are silently reused for WES. Missing evidence is never zero depth.
"""
from dataclasses import dataclass
import math
import statistics


def to_natural_log(value, unit):
    if unit not in ('log2_alt_ref', 'ln_alt_ref'):
        raise ValueError('explicit alt/ref log unit required')
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError('nonfinite log ratio')
    return float(value) * math.log(2) if unit == 'log2_alt_ref' else float(value)


def subtract_bias(raw_ln_alt_ref, bias, *, bias_unit):
    """Preserve missing bias rather than treating it as no correction."""
    raw = to_natural_log(raw_ln_alt_ref, 'ln_alt_ref')
    bias_ln = to_natural_log(bias, bias_unit)
    if raw is None or bias_ln is None:
        return None
    return raw - bias_ln


def allelic_log_ratio(ref, alt, *, pseudocount):
    if any(type(value) is not int or value < 0 for value in (ref, alt)):
        raise ValueError('counts must be nonnegative integers')
    if (isinstance(pseudocount, bool) or not isinstance(pseudocount, (int, float))
            or not math.isfinite(pseudocount) or pseudocount < 0):
        raise ValueError('explicit nonnegative pseudocount required')
    if ref + alt == 0 or min(ref + pseudocount, alt + pseudocount) <= 0:
        return None
    return math.log(alt + pseudocount) - math.log(ref + pseudocount)


@dataclass(frozen=True)
class MarkerPolicy:
    min_depth: int
    max_depth: int
    het_baf_lower: float
    het_baf_upper: float
    min_heterozygous_normals: int
    min_observed_fraction: float
    min_median_depth: float
    max_abs_median_ln_bias: float
    max_sd_ln_bias: float
    allelic_pseudocount: float

    def validate(self):
        if (type(self.min_depth) is not int or type(self.max_depth) is not int or
                not 0 < self.min_depth <= self.max_depth):
            raise ValueError('invalid depth range')
        if type(self.min_heterozygous_normals) is not int or self.min_heterozygous_normals < 2:
            raise ValueError('at least two heterozygous normals required to estimate dispersion')
        fields = (self.het_baf_lower, self.het_baf_upper, self.min_observed_fraction,
                  self.min_median_depth, self.max_abs_median_ln_bias, self.max_sd_ln_bias,
                  self.allelic_pseudocount)
        if any(isinstance(x, bool) or not isinstance(x, (int, float)) or
               not math.isfinite(x) for x in fields):
            raise ValueError('invalid/nonfinite marker policy')
        if not 0 < self.het_baf_lower < 0.5 < self.het_baf_upper < 1:
            raise ValueError('heterozygosity bounds must bracket 0.5')
        if not 0 < self.min_observed_fraction <= 1:
            raise ValueError('invalid observed fraction')
        if min(self.min_median_depth, self.max_abs_median_ln_bias, self.max_sd_ln_bias,
               self.allelic_pseudocount) < 0:
            raise ValueError('negative marker threshold')


def summarize_marker(counts_by_sample, pon_sample_ids, policy):
    """counts are (reference, alternate), or None when not evaluated.

    pon_sample_ids must come from the role-validated assay manifest. Extra IDs
    are rejected; omitted IDs count as missing, not as zero. A zero/zero row
    contributes depth zero to the observed median but never to allelic bias.
    """
    policy.validate()
    ids = list(pon_sample_ids)
    if not ids or len(ids) != len(set(ids)) or any(not isinstance(x, str) or not x for x in ids):
        raise ValueError('unique nonempty PoN sample IDs required')
    if set(counts_by_sample) - set(ids):
        raise ValueError('non-PoN sample supplied')
    depths, ln_ratios = [], []
    eligible_depth_samples = 0
    zero_depth_samples = 0
    for sample in ids:
        counts = counts_by_sample.get(sample)
        if counts is None:
            continue
        if not isinstance(counts, (list, tuple)) or len(counts) != 2:
            raise ValueError('expected ref/alt count pair')
        ref, alt = counts
        if any(type(value) is not int or value < 0 for value in (ref, alt)):
            raise ValueError('counts must be nonnegative integers')
        depth = ref + alt
        depths.append(depth)
        zero_depth_samples += depth == 0
        if not policy.min_depth <= depth <= policy.max_depth:
            continue
        eligible_depth_samples += 1
        baf = alt / depth
        if policy.het_baf_lower < baf < policy.het_baf_upper:
            ln_ratios.append(allelic_log_ratio(ref, alt, pseudocount=policy.allelic_pseudocount))
    observed_fraction = len(depths)/len(ids)
    median_depth = statistics.median(depths) if depths else None
    median_bias = statistics.median(ln_ratios) if ln_ratios else None
    sd_bias = statistics.stdev(ln_ratios) if len(ln_ratios) >= 2 else None
    reasons = []
    if observed_fraction < policy.min_observed_fraction:
        reasons.append('INSUFFICIENT_OBSERVED_NORMALS')
    if median_depth is None or median_depth < policy.min_median_depth:
        reasons.append('INSUFFICIENT_MEDIAN_DEPTH')
    if len(ln_ratios) < policy.min_heterozygous_normals:
        reasons.append('INSUFFICIENT_HETEROZYGOUS_NORMALS')
    if median_bias is not None and abs(median_bias) > policy.max_abs_median_ln_bias:
        reasons.append('EXCESS_MEDIAN_ALLELIC_BIAS')
    if sd_bias is not None and sd_bias > policy.max_sd_ln_bias:
        reasons.append('EXCESS_ALLELIC_DISPERSION')
    stable = not reasons
    return dict(status='STABLE_TRAINING_MARKER' if stable else 'EXCLUDED_TRAINING_MARKER',
                reasons=reasons, total_normals=len(ids), observed_normals=len(depths),
                missing_normals=len(ids)-len(depths), zero_depth_normals=zero_depth_samples,
                observed_fraction=observed_fraction, median_observed_depth=median_depth,
                eligible_depth_normals=eligible_depth_samples, heterozygous_normals=len(ln_ratios),
                bias_unit='ln_alt_ref', allelic_pseudocount=policy.allelic_pseudocount,
                raw_median_bias=median_bias, raw_sd_bias=sd_bias,
                correction_bias=median_bias if stable else None,
                sensitivity_calibrated=False, germline_genotype_assigned=False)
