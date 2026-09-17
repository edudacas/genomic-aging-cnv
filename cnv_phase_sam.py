"""Extract biallelic SNP observations from an indexed BAM's SAM subset.

SAM SEQ is already alignment-oriented, including reverse-strand records.
The source BAM/reference alleles still require manifest and FASTA validation.
"""
from collections import Counter
import re
from cnv_commands import Command


def validate_loci(loci):
    if len(loci)!=2 or len({x['locus_id'] for x in loci})!=2:
        raise ValueError('two distinct SNP loci required')
    positions=set()
    for locus in loci:
        if (not isinstance(locus['locus_id'],str) or not locus['locus_id'] or
                not isinstance(locus['chromosome'],str) or
                not re.fullmatch(r'[A-Za-z0-9_.-]+',locus['chromosome']) or
                type(locus['position']) is not int or locus['position']<1 or
                locus['ref'] not in ('A','C','G','T') or locus['alt'] not in ('A','C','G','T') or
                locus['ref']==locus['alt']):
            raise ValueError('only explicit biallelic SNP loci supported')
        key=(locus['chromosome'],locus['position'])
        if key in positions:
            raise ValueError('phase loci must occupy distinct positions')
        positions.add(key)


def phase_subset_command(bam,output,loci,*,samtools):
    validate_loci(loci)
    # -M selects the region union without reporting a spanning alignment twice.
    # Reference: https://www.htslib.org/doc/samtools-view.html
    regions=[f"{x['chromosome']}:{x['position']}-{x['position']}" for x in loci]
    return Command((str(samtools),'view','-h','-M','-o',str(output),str(bam),*regions),(str(output),))


def query_positions(cigar,start,sequence_length,positions):
    operations=re.findall(r'([1-9][0-9]*)([MIDNSHP=X])',cigar)
    if ''.join(n+op for n,op in operations)!=cigar or not operations:
        raise ValueError('unsupported or invalid CIGAR')
    reference,query=start,0
    found={}
    for count,op in operations:
        count=int(count)
        if op in 'M=X':
            for position in positions:
                if reference<=position<reference+count:
                    found[position]=query+position-reference
            reference+=count; query+=count
        elif op in 'IS':
            query+=count
        elif op in 'DN':
            reference+=count
    if query!=sequence_length:
        raise ValueError('CIGAR/query length mismatch')
    return found


def extract_snp_observations(path,sample_id,loci,*,expected_bam_sample_name,
                             expected_reference_lengths=None):
    validate_loci(loci)
    if any(not isinstance(x,str) or not x for x in (sample_id,expected_bam_sample_name)):
        raise ValueError('explicit study ID and expected BAM SM required')
    groups={}
    sequences={}
    observations=[]
    counts=Counter()
    with open(path) as stream:
        for line in stream:
            fields=line.rstrip('\r\n').split('\t')
            if line.startswith('@'):
                if counts['alignment_records']:
                    raise ValueError('SAM header after alignments')
                if fields[0]=='@RG':
                    tags=dict(item.split(':',1) for item in fields[1:])
                    if not tags.get('ID') or tags.get('SM')!=expected_bam_sample_name or tags['ID'] in groups:
                        raise ValueError('SAM read-group/sample identity mismatch')
                    groups[tags['ID']]=tags['SM']
                elif fields[0]=='@SQ':
                    tags=dict(item.split(':',1) for item in fields[1:])
                    if not tags.get('SN') or tags['SN'] in sequences or int(tags.get('LN','0'))<=0:
                        raise ValueError('invalid SAM reference dictionary')
                    sequences[tags['SN']]=int(tags['LN'])
                continue
            if len(fields)<11:
                raise ValueError('truncated SAM record')
            counts['alignment_records']+=1
            flag=int(fields[1])
            if not 0<=flag<=65535 or fields[0] in ('','*'):
                raise ValueError('SAM flag or template identity invalid')
            if flag&4:
                counts['unmapped_records']+=1
                continue
            rg=[item[5:] for item in fields[11:] if item.startswith('RG:Z:')]
            if len(rg)!=1 or rg[0] not in groups:
                raise ValueError('mapped read lacks a unique declared read group')
            sequence,qualities=fields[9],fields[10]
            if sequence=='*' or qualities=='*':
                counts['missing_sequence_or_quality_records']+=1
                continue
            if len(sequence)!=len(qualities) or any(not 33<=ord(x)<=126 for x in qualities):
                raise ValueError('invalid SAM base qualities')
            selected=[x for x in loci if x['chromosome']==fields[2]]
            if int(fields[3])<1:
                raise ValueError('mapped SAM position must be positive')
            mapped=query_positions(fields[5],int(fields[3]),len(sequence),[x['position'] for x in selected])
            for locus in selected:
                if locus['position'] not in mapped:
                    counts['locus_without_aligned_base']+=1
                    continue
                q=mapped[locus['position']]
                base=sequence[q].upper()
                observations.append(dict(sample_id=sample_id,read_group=rg[0],query_name=fields[0],
                    locus_id=locus['locus_id'],allele='REF' if base==locus['ref'] else 'ALT' if base==locus['alt'] else 'OTHER',
                    mapq=int(fields[4]),baseq=ord(qualities[q])-33,flag=flag))
    if not groups:
        raise ValueError('SAM requires declared sample read groups')
    if expected_reference_lengths is not None and sequences!=expected_reference_lengths:
        raise ValueError('SAM/reference sequence dictionary mismatch')
    observations.sort(key=lambda row:(row['read_group'],row['query_name']))
    return observations,dict(counts)


def phase_reference_command(fasta,output,loci,*,samtools):
    validate_loci(loci)
    regions=[f"{x['chromosome']}:{x['position']}-{x['position']}" for x in loci]
    return Command((str(samtools),'faidx','-o',str(output),str(fasta),*regions),(str(output),))


def verify_reference_bases(path,loci):
    """Validate exact faidx single-base output; no contig aliasing or allele flip."""
    validate_loci(loci)
    expected={f"{x['chromosome']}:{x['position']}-{x['position']}":x['ref'] for x in loci}
    records={}
    current=None
    with open(path) as stream:
        for line in stream:
            value=line.strip()
            if value.startswith('>'):
                current=value[1:]
                if current not in expected or current in records:
                    raise ValueError('reference subset locus mismatch')
                records[current]=''
            elif value:
                if current is None:
                    raise ValueError('reference sequence without header')
                records[current]+=value.upper()
    if records!=expected:
        raise ValueError('SNP REF differs from reference FASTA or locus is missing')
    return dict(status='PHASE_SNP_REFERENCE_MATCH',loci=len(loci),reference_build_inferred=False)
