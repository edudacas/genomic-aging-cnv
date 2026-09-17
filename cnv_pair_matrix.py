"""Role-checked streaming sample/reference count matrices for the R adapter.

The caller must verify native pileup receipts before calling this function.
It validates role/stratum/sex and data geometry, not upstream file provenance.
"""
from cnv_contract import validate_manifest
from cnv_pileup import align_catalog
from cnv_commands import snp_pileup
from cnv_execution import verify_completed_task


def verified_pair_matrix(manifest,sample_id,reference_id,catalog,registry,lengths,*,
                         mask,executable,pileup_policy,rehash_inputs):
    """Verify both native tasks before yielding any paired count row.

    registry maps each sample to {task_directory, output}. Receipt verification
    does not substitute for identity validation or reference-resource calibration.
    """
    validate_manifest(manifest,stage='reference_build')
    samples={row['sample_id']:row for row in manifest['samples']}
    if set(registry)!={sample_id,reference_id} or not set(registry)<=set(samples):
        raise ValueError('registry must match known sample/reference pair')
    required={'mapq','baseq','pseudo_spacing','max_depth','count_orphans','min_read_count'}
    if set(pileup_policy)!=required:
        raise ValueError('complete explicit pileup policy required')
    paths,receipts={},{}
    for sid,entry in registry.items():
        sample=samples[sid]
        vcf=manifest['strata'][sample['stratum']]['dense_snp_vcf']
        command=snp_pileup(sample['bam']['path'],vcf['path'],entry['output'],
                          executable=executable,**pileup_policy)
        artifacts=(sample['bam'],sample['bam_index'],vcf)
        receipts[sid]=verify_completed_task(command,entry['task_directory'],
            {item['path']:item['sha256'] for item in artifacts},rehash_inputs=rehash_inputs)
        paths[sid]=entry['output']
    for row in pair_matrix(manifest,sample_id,reference_id,catalog,paths,lengths,mask=mask):
        yield dict(row,pileup_receipt_bindings_verified=True,
                   upstream_input_bytes_rehashed=rehash_inputs,
                   pileup_task_fingerprints={sid:receipt['fingerprint']
                                             for sid,receipt in receipts.items()})


def pair_matrix(manifest, sample_id, reference_id, catalog, pileup_paths, lengths, *, mask):
    validate_manifest(manifest, stage='reference_build')
    samples = {row['sample_id']:row for row in manifest['samples']}
    if sample_id not in samples or reference_id not in samples:
        raise ValueError('unknown sample/reference')
    sample, normal = samples[sample_id], samples[reference_id]
    if normal['role']!='pon' or sample['person_id']==normal['person_id']:
        raise ValueError('independent PoN reference required')
    if sample['stratum']!=normal['stratum']:
        raise ValueError('sample/reference stratum mismatch')
    if (sample.get('reference_sex') not in ('Male','Female') or
            normal.get('reference_sex')!=sample['reference_sex']):
        raise ValueError('resolved matching reference sex required')
    if set(pileup_paths)!={sample_id,reference_id}:
        raise ValueError('pileup inputs must match the explicit sample/reference pair')
    if mask.lengths!=lengths:
        raise ValueError('mask/reference dictionary mismatch')
    for index,(locus,counts) in enumerate(align_catalog(catalog,pileup_paths,lengths),1):
        tumor, control = counts[sample_id],counts[reference_id]
        reasons = []
        if tumor is None:
            reasons.append('SAMPLE_PILEUP_ROW_ABSENT')
        if control is None:
            reasons.append('REFERENCE_PILEUP_ROW_ABSENT')
        masked = mask.overlap(locus['chromosome'],locus['position']-1,locus['position'])>0
        row = dict(catalog_row=index,sample_id=sample_id,reference_id=reference_id,
                   Chromosome=locus['chromosome'],Position=locus['position'],
                   Ref=locus['ref'],Alt=locus['alt'],catalog_kind=locus['kind'],masked=masked,
                   reasons=reasons,native_input_eligible=not reasons,
                   input_provenance_verified=False)
        for prefix,value in (('TUM',tumor),('NOR',control)):
            row[prefix+'.DP'] = value[0]+value[1] if value is not None else None
            row[prefix+'.RD'] = value[0] if value is not None else None
            row[prefix+'.AD'] = value[1] if value is not None else None
            row[prefix+'.ERROR'] = value[2] if value is not None else None
            row[prefix+'.DEL'] = value[3] if value is not None else None
        # Observed zero remains eligible for explicit native depth filtering;
        # missing rows go to the exclusion ledger and never into native R.
        yield row


def native_rcmat_rows(rows):
    """Select complete numeric inputs without deleting the separate full ledger."""
    columns = ('Chromosome','Position','NOR.DP','NOR.RD','TUM.DP','TUM.RD',
               'catalog_kind','masked','catalog_row')
    for row in rows:
        if row['native_input_eligible']:
            yield {key:row[key] for key in columns}
