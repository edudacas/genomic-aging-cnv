"""Lossless exploratory single-sample Manta VCF evidence, not CNV event calls.

Pinned interface reference: Illumina/manta v1.6.0 docs/userGuide/README.md.
Tumor-only mode has no somatic scoring model. PR and SR are retained separately,
not summed as independent fragments. Breakend records are not counted as events.
"""
import gzip


def support_pair(value):
    if value in (None,'.'):
        return None
    fields=value.split(',')
    if len(fields)!=2:
        raise ValueError('Manta PR/SR requires ref,alt support pair')
    result=[]
    for field in fields:
        if field=='.':
            result.append(None)
        else:
            if not field.isdigit():
                raise ValueError('invalid Manta support count')
            result.append(int(field))
    return result


def manta_evidence(path,sample_id,*,expected_vcf_sample,lengths):
    if any(not isinstance(x,str) or not x for x in (sample_id,expected_vcf_sample)):
        raise ValueError('explicit study and VCF sample identities required')
    opener=gzip.open if str(path).endswith('.gz') else open
    header=None
    ordinal=0
    with opener(path,'rt') as stream:
        for line in stream:
            if line.startswith('##'):
                if header is not None:
                    raise ValueError('VCF metadata after column header')
                continue
            raw=line.rstrip('\r\n')
            fields=raw.split('\t')
            if line.startswith('#CHROM'):
                expected=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',expected_vcf_sample]
                if header is not None or fields!=expected:
                    raise ValueError('exact single-sample VCF header required')
                header=fields
                continue
            if header is None or len(fields)!=10:
                raise ValueError('VCF record/header mismatch')
            chrom,pos,identifier,ref,alt,qual,filters,info,fmt,sample=fields
            if chrom not in lengths or not pos.isdigit() or not 1<=int(pos)<=lengths[chrom]:
                raise ValueError('VCF coordinate outside reference')
            if not ref or not alt or ',' in alt:
                raise ValueError('expected single alternate Manta record')
            info_values={}
            if info!='.':
                for item in info.split(';'):
                    parts=item.split('=',1)
                    if not parts[0] or parts[0] in info_values:
                        raise ValueError('empty/duplicate INFO key')
                    info_values[parts[0]]=parts[1] if len(parts)==2 else True
            keys=[] if fmt=='.' else fmt.split(':')
            values=[] if sample=='.' and not keys else sample.split(':')
            if len(keys)!=len(set(keys)) or len(values)>len(keys):
                raise ValueError('invalid VCF FORMAT fields')
            # VCF trailing missing sample fields may be omitted; never zero-fill.
            sample_values=dict(zip(keys,values+['.']*(len(keys)-len(values))))
            ordinal+=1
            yield dict(sample_id=sample_id,vcf_sample_name=expected_vcf_sample,
                source_record_index=ordinal,source_record_id=identifier,
                chromosome=chrom,position_1based=int(pos),ref=ref,alt=alt,
                qual_raw=qual,filter_raw=filters,info_raw=info,format_raw=fmt,sample_raw=sample,
                source_record_raw=raw,source_info=info_values,source_sample_fields=sample_values,
                paired_read_support=support_pair(sample_values.get('PR')),
                split_read_support=support_pair(sample_values.get('SR')),
                source_svtype=info_values.get('SVTYPE'),source_mateid=info_values.get('MATEID'),
                source_event=info_values.get('EVENT'),source_imprecise='IMPRECISE' in info_values,
                evidence_scope='EXPLORATORY_STRUCTURAL_VARIANT_RECORD',
                event_count=None,cnv_event_class='UNCLASSIFIED',somatic_origin_classified=False,
                caller_source_provenance_verified=False,phase_assignment=None)
    if header is None:
        raise ValueError('missing VCF column header')
