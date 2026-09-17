"""Manifest and native-receipt binding for two-SNP fragment evidence."""
from cnv_contract import validate_manifest,validate_artifact
from cnv_coordinates import fai_lengths
from cnv_execution import verify_completed_task
from cnv_phase_sam import (phase_subset_command,phase_reference_command,
    verify_reference_bases,extract_snp_observations,validate_loci)
from cnv_phasing import fragment_evidence,summarize_fragment_evidence


def verified_phase_evidence(manifest,sample_id,loci,sam_task,reference_task,*,
                            samtools,policy,rehash_inputs):
    """Verify completed local tasks only. Does not launch samtools or a cohort.

    bam_sample_name is an explicit manifest field, distinct from sample_id.
    Matching declared SM is not independent proof of clinical identity.
    """
    validate_manifest(manifest,stage='reference_build')
    validate_loci(loci)
    policy.validate()
    samples={row['sample_id']:row for row in manifest['samples']}
    if sample_id not in samples:
        raise ValueError('unknown phase sample')
    sample=samples[sample_id]
    sm=sample.get('bam_sample_name')
    if not isinstance(sm,str) or not sm:
        raise ValueError('explicit manifest bam_sample_name required')
    spec=manifest['strata'][sample['stratum']]
    validate_artifact(spec['reference_fai'],True)
    lengths=fai_lengths(spec['reference_fai']['path'])
    for locus in loci:
        if locus['chromosome'] not in lengths or locus['position']>lengths[locus['chromosome']]:
            raise ValueError('phase locus outside declared reference')
    sam_command=phase_subset_command(sample['bam']['path'],sam_task['output'],loci,samtools=samtools)
    sam_verified=verify_completed_task(sam_command,sam_task['task_directory'],
        {a['path']:a['sha256'] for a in (sample['bam'],sample['bam_index'])},rehash_inputs=rehash_inputs)
    reference_command=phase_reference_command(spec['reference_fasta']['path'],reference_task['output'],
                                              loci,samtools=samtools)
    reference_verified=verify_completed_task(reference_command,reference_task['task_directory'],
        {a['path']:a['sha256'] for a in (spec['reference_fasta'],spec['reference_fai'])},
        rehash_inputs=rehash_inputs)
    verify_reference_bases(reference_task['output'],loci)
    observations,diagnostics=extract_snp_observations(sam_task['output'],sample_id,loci,
                                                   expected_bam_sample_name=sm,
                                                   expected_reference_lengths=lengths)
    ids=[x['locus_id'] for x in loci]
    summary=summarize_fragment_evidence(observations,sample_id,ids,policy=policy)
    fragments=list(fragment_evidence(observations,sample_id,ids,policy=policy))
    binding=dict(bam_provenance_verified=True,
        bam_provenance_scope='MANIFEST_NATIVE_TASK_INPUT_OUTPUT_BINDINGS',
        sam_task_fingerprint=sam_verified['fingerprint'],
        reference_task_fingerprint=reference_verified['fingerprint'],
        upstream_input_bytes_rehashed=rehash_inputs,clinical_identity_verified=False)
    return dict(summary=dict(summary,**binding),
                fragments=[dict(row,**binding) for row in fragments],
                extraction_diagnostics=diagnostics,loci=loci)
