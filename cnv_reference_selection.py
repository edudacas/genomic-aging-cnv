"""PoN reference ranking on a common, predeclared set of corrected-logR loci."""
import hashlib
import json
import math

from cnv_contract import validate_manifest
from cnv_commands import integer


def select_reference(manifest, sample_id, candidate_signals, selection_loci, *,
                     minimum_common_loci, minimum_common_fraction):
    """Signals map PoN ID -> locus ID -> corrected log2 ratio or None.

    Requires all same-stratum, same-reference-sex PoN candidates. This prevents
    silent outcome-driven candidate omission. The signal-generation provenance
    must be verified by orchestration; this function verifies selection logic.
    """
    validate_manifest(manifest, stage='reference_build')
    integer(minimum_common_loci, 'minimum_common_loci')
    if (isinstance(minimum_common_fraction,bool) or
        not isinstance(minimum_common_fraction,(int,float)) or
        not math.isfinite(minimum_common_fraction) or not 0 < minimum_common_fraction <= 1):
        raise ValueError('invalid common-locus fraction')
    by_id = {sample['sample_id']:sample for sample in manifest['samples']}
    if sample_id not in by_id:
        raise ValueError('unknown sample')
    sample = by_id[sample_id]
    sex = sample.get('reference_sex')
    if sex not in ('Male','Female'):
        raise ValueError('reference sex unresolved; never defaults to Female')
    candidates = []
    for normal in manifest['samples']:
        if normal['role'] != 'pon' or normal['stratum'] != sample['stratum']:
            continue
        if normal['sample_id'] == sample_id:
            continue
        if normal.get('reference_sex') not in ('Male','Female'):
            raise ValueError('PoN reference sex unresolved')
        if normal['reference_sex'] == sex:
            candidates.append(normal['sample_id'])
    candidates.sort()
    if not candidates or set(candidate_signals) != set(candidates):
        raise ValueError('expected exactly all eligible independent PoN references')
    loci = list(selection_loci)
    if not loci or any(not isinstance(x,str) or not x for x in loci) or len(loci)!=len(set(loci)):
        raise ValueError('unique frozen selection loci required')
    locus_set = set(loci)
    valid_sets = []
    for candidate in candidates:
        signal = candidate_signals[candidate]
        if set(signal)-locus_set:
            raise ValueError('candidate signal contains undeclared selection loci')
        valid = set()
        for locus,value in signal.items():
            if value is None:
                continue
            if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
                raise ValueError('invalid corrected logR; use None for missing')
            valid.add(locus)
        valid_sets.append(valid)
    shared = set.intersection(*valid_sets)
    common = [locus for locus in loci if locus in shared]
    fraction = len(common)/len(loci)
    if len(common)<minimum_common_loci or fraction<minimum_common_fraction:
        raise ValueError('insufficient common loci for comparable reference ranking')
    scores = {candidate:math.fsum(candidate_signals[candidate][locus]**2 for locus in common)
              for candidate in candidates}
    if any(not math.isfinite(value) for value in scores.values()):
        raise ValueError('reference score overflow')
    chosen = min(candidates,key=lambda candidate:(scores[candidate],candidate))
    encoded = json.dumps(common,separators=(',',':')).encode()
    return dict(status='REFERENCE_SELECTED_ON_COMMON_LOCI',sample_id=sample_id,
        selected_reference=chosen,reference_sex=sex,candidate_scores=scores,
        common_loci=len(common),selection_loci=len(loci),common_fraction=fraction,
        common_loci_sha256=hashlib.sha256(encoded).hexdigest(),
        missing_by_candidate={c:len(loci)-len(v) for c,v in zip(candidates,valid_sets)},
        score='SUM_SQUARED_CORRECTED_LOG2_RATIO_ON_IDENTICAL_LOCI',
        tie_break='LEXICAL_SAMPLE_ID',signal_provenance_verified=False,
        biological_validation_performed=False)
