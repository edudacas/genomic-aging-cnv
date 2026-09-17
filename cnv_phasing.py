"""Two-locus fragment evidence without population/statistical phase inference.

Input observations must be grouped by (read_group, query_name) in lexical order
and restricted to one explicit sample. This does not merge distinct templates
by UMI, assert unique molecular identity, or establish a somatic second hit.
Allele extraction/BAM provenance is a separate upstream prerequisite.
"""
from collections import Counter
from dataclasses import dataclass


@dataclass(frozen=True)
class PhaseReadPolicy:
    min_mapq: int
    min_baseq: int
    exclude_flags: int

    def validate(self):
        if any(type(x) is not int or x<0 for x in (self.min_mapq,self.min_baseq,self.exclude_flags)):
            raise ValueError('explicit nonnegative integer read policy required')
        if self.min_mapq>255 or self.min_baseq>93 or self.exclude_flags>65535:
            raise ValueError('read policy outside SAM ranges')
        # Unmapped, secondary and supplementary records cannot create links.
        if self.exclude_flags & (4|256|2048) != (4|256|2048):
            raise ValueError('phase evidence must exclude unmapped/secondary/supplementary')


def fragment_evidence(observations, sample_id, locus_ids, *, policy):
    policy.validate()
    if not isinstance(sample_id,str) or not sample_id:
        raise ValueError('explicit sample ID required')
    loci=tuple(locus_ids)
    if len(loci)!=2 or any(not isinstance(x,str) or not x for x in loci) or loci[0]==loci[1]:
        raise ValueError('two distinct ordered locus IDs required')
    current=None
    seen=Counter()
    accepted={locus:set() for locus in loci}
    filtered=Counter()
    def emit(key):
        calls={locus:(next(iter(values)) if len(values)==1 and 'OTHER' not in values else None)
               for locus,values in accepted.items()}
        ambiguous=any(len(v)>1 or 'OTHER' in v for v in accepted.values())
        if ambiguous:
            status='AMBIGUOUS_FRAGMENT'
        elif all(calls[locus] is not None for locus in loci):
            status='LINKED_FRAGMENT'
        elif any(calls[locus] is not None for locus in loci):
            status='SINGLE_LOCUS_FRAGMENT'
        else:
            status='NO_ELIGIBLE_BASES'
        return dict(sample_id=sample_id,read_group=key[0],query_name=key[1],
            locus_ids=list(loci),status=status,alleles=calls,
            haplotype='_'.join(calls[locus] for locus in loci) if status=='LINKED_FRAGMENT' else None,
            observed_records=dict(seen),accepted_allele_sets={k:sorted(v) for k,v in accepted.items()},
            filtered_records=dict(filtered),unique_molecule_verified=False,bam_provenance_verified=False)
    required={'sample_id','read_group','query_name','locus_id','allele','mapq','baseq','flag'}
    for observation in observations:
        if not required<=set(observation):
            raise ValueError('incomplete phase observation')
        if observation['sample_id']!=sample_id or observation['locus_id'] not in loci:
            raise ValueError('foreign sample or locus in phase evidence')
        key=(observation['read_group'],observation['query_name'])
        if any(not isinstance(x,str) or not x for x in key):
            raise ValueError('read-group and template identity required')
        if current is not None and key<current:
            raise ValueError('observations must be template-grouped and sorted')
        if current is not None and key!=current:
            yield emit(current)
            seen=Counter(); filtered=Counter(); accepted={locus:set() for locus in loci}
        current=key
        for field,maximum in (('mapq',255),('baseq',93),('flag',65535)):
            value=observation[field]
            if type(value) is not int or not 0<=value<=maximum:
                raise ValueError('invalid SAM observation '+field)
        if observation['allele'] not in ('REF','ALT','OTHER'):
            raise ValueError('explicit REF/ALT/OTHER observation required')
        locus=observation['locus_id']
        seen[locus]+=1
        if observation['flag'] & policy.exclude_flags:
            filtered['EXCLUDED_FLAG']+=1
        elif observation['mapq']==255:
            filtered['MAPQ_UNAVAILABLE']+=1
        elif observation['mapq']<policy.min_mapq:
            filtered['LOW_MAPQ']+=1
        elif observation['baseq']<policy.min_baseq:
            filtered['LOW_BASEQ']+=1
        else:
            accepted[locus].add(observation['allele'])
    if current is not None:
        yield emit(current)


def summarize_fragment_evidence(observations,sample_id,locus_ids,*,policy):
    locus_ids=tuple(locus_ids)
    statuses=Counter()
    haplotypes=Counter({key:0 for key in ('REF_REF','REF_ALT','ALT_REF','ALT_ALT')})
    for row in fragment_evidence(observations,sample_id,locus_ids,policy=policy):
        statuses[row['status']]+=1
        if row['haplotype'] is not None:
            haplotypes[row['haplotype']]+=1
    linked=sum(haplotypes.values())
    return dict(sample_id=sample_id,locus_ids=list(locus_ids),fragment_status_counts=dict(statuses),
        linked_fragments=linked,haplotype_fragment_counts=dict(haplotypes),
        status='FRAGMENT_LINKAGE_OBSERVED' if linked else 'PHASE_NOT_EVALUATED',
        phase_assignment=None,phase_confidence_calibrated=False,
        somatic_second_hit_confirmed=False,unique_molecules_counted=False,
        bam_provenance_verified=False)
