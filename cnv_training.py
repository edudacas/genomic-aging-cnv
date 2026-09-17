"""Manifest-bound PoN marker training, never study/held-out calibration data.

All same-stratum, same-resolved-sex PoN members are required. Sex is explicit
because pooling sex chromosomes across ploidies requires a separate policy.
This module does not infer sex or establish that a PoN is clinically suitable.
Exhaust the generator before publishing its outputs or completion receipt.
"""
from cnv_contract import validate_manifest
from cnv_commands import snp_pileup
from cnv_execution import verify_completed_task
from cnv_pileup import training_summaries


def verified_training_summaries(manifest, stratum_id, reference_sex, catalog,
                                registry, lengths, *, mask, executable,
                                pileup_policy, marker_policy, rehash_inputs):
    validate_manifest(manifest, stage='reference_build')
    marker_policy.validate()
    if reference_sex not in ('Male', 'Female'):
        raise ValueError('resolved explicit reference sex required')
    if stratum_id not in manifest['strata']:
        raise ValueError('unknown training stratum')
    if mask.lengths != lengths:
        raise ValueError('mask/reference dictionary mismatch')
    normals = [s for s in manifest['samples'] if s['role']=='pon' and s['stratum']==stratum_id]
    if any(s.get('reference_sex') not in ('Male','Female') for s in normals):
        raise ValueError('resolve all stratum PoN reference sexes before training')
    selected = {s['sample_id']:s for s in normals if s['reference_sex']==reference_sex}
    if len(selected) < marker_policy.min_heterozygous_normals:
        raise ValueError('insufficient eligible PoN members for marker policy')
    if set(registry) != set(selected):
        raise ValueError('registry must contain exactly all eligible PoN members')
    required = {'mapq','baseq','pseudo_spacing','max_depth','count_orphans','min_read_count'}
    if set(pileup_policy) != required:
        raise ValueError('complete explicit pileup policy required')
    vcf = manifest['strata'][stratum_id]['dense_snp_vcf']
    paths, fingerprints = {}, {}
    for sid in sorted(selected):
        sample, entry = selected[sid], registry[sid]
        command = snp_pileup(sample['bam']['path'],vcf['path'],entry['output'],
                            executable=executable,**pileup_policy)
        receipt = verify_completed_task(command,entry['task_directory'],
            {a['path']:a['sha256'] for a in (sample['bam'],sample['bam_index'],vcf)},
            rehash_inputs=rehash_inputs)
        fingerprints[sid] = receipt['fingerprint']
        paths[sid] = entry['output']
    for row in training_summaries(catalog,paths,lengths,sorted(selected),marker_policy):
        masked = mask.overlap(row['chromosome'],row['position']-1,row['position']) > 0
        if masked:
            row = dict(row,status='MASKED_TRAINING_LOCUS',correction_bias=None,
                       reasons=[*row.get('reasons',[]),'EXCLUDED_BY_MASK'])
        yield dict(row,masked=masked,training_stratum=stratum_id,
                   training_reference_sex=reference_sex,
                   pon_sample_ids=sorted(selected),pileup_task_fingerprints=fingerprints,
                   pileup_receipt_bindings_verified=True,
                   upstream_input_bytes_rehashed=rehash_inputs,
                   sensitivity_calibrated=False,clinical_identity_verified=False)
