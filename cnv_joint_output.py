"""Structural validation of adapted native joint outputs, not event validation."""
import csv
import math
from pathlib import Path


def rows(path,required):
    with Path(path).open(newline='') as stream:
        reader=csv.DictReader(stream,delimiter='\t')
        fields=reader.fieldnames or []
        if len(fields)!=len(set(fields)) or not set(required)<=set(fields):
            raise ValueError('joint output schema mismatch: '+Path(path).name)
        result=list(reader)
        if any(None in row or any(v is None for v in row.values()) for row in result):
            raise ValueError('ragged joint output')
        return result


def integer(value,minimum):
    number=float(value)
    if not math.isfinite(number) or number<minimum or number!=int(number):
        raise ValueError('invalid joint integer field')
    return int(number)


def validate_joint_output(directory,expected_input_rows,*,input_matrix=None):
    directory=Path(directory)
    audit=rows(directory/'locus_audit.tsv',('catalog_row','chromosome','position_1based','signal_status'))
    if len(audit)!=expected_input_rows:
        raise ValueError('joint locus audit does not cover prepared numeric rows')
    ids=[integer(row['catalog_row'],1) for row in audit]
    if len(set(ids))!=len(ids):
        raise ValueError('duplicate audited catalog row')
    for row in audit:
        integer(row['position_1based'],1)
    if input_matrix is not None:
        original=rows(input_matrix,('catalog_row','Chromosome','Position','NOR.DP','TUM.DP'))
        indexed={integer(row['catalog_row'],1):row for row in original}
        if len(indexed)!=len(original) or set(indexed)!=set(ids):
            raise ValueError('audit catalog rows differ from submitted matrix')
        for row in audit:
            source=indexed[integer(row['catalog_row'],1)]
            if (row['chromosome']!=source['Chromosome'] or
                    integer(row['position_1based'],1)!=integer(source['Position'],1) or
                    integer(row.get('normal_depth','nan'),0)!=integer(source['NOR.DP'],0) or
                    integer(row.get('sample_depth','nan'),0)!=integer(source['TUM.DP'],0)):
                raise ValueError('audit locus/counts differ from submitted matrix')
    signal=rows(directory/'signal.tsv',('chrom','maploc','cnlr','valor','het','source_catalog_row'))
    if not signal:
        raise ValueError('empty joint signal')
    signal_ids=[]
    for row in signal:
        signal_ids.append(integer(row['source_catalog_row'],1))
        integer(row['chrom'],1); integer(row['maploc'],1)
        if integer(row['het'],0) not in (0,1) or not math.isfinite(float(row['cnlr'])):
            raise ValueError('invalid joint signal')
        if row['valor']!='NA' and not math.isfinite(float(row['valor'])):
            raise ValueError('nonfinite joint allele signal')
    if len(signal_ids)!=len(set(signal_ids)) or not set(signal_ids)<=set(ids):
        raise ValueError('joint signal lacks unique input linkage')
    segments=rows(directory/'segments.tsv',('chromosome','start_0based','end_0based_exclusive',
        'n_markers','n_informative_allelic','mean_log2_ratio','event_class','evaluable_bp','somatic_origin_classified'))
    if not segments:
        raise ValueError('empty joint segments')
    markers=0
    for row in segments:
        start=integer(row['start_0based'],0); end=integer(row['end_0based_exclusive'],1)
        n=integer(row['n_markers'],1); allelic=integer(row['n_informative_allelic'],0)
        if start>=end or allelic>n or not math.isfinite(float(row['mean_log2_ratio'])):
            raise ValueError('invalid joint segment geometry or signal')
        if row['event_class']!='UNCLASSIFIED' or row['evaluable_bp']!='NA' or row['somatic_origin_classified']!='FALSE':
            raise ValueError('native candidate output contains unsupported biological claim')
        markers+=n
    if markers>len(signal):
        raise ValueError('segmented markers exceed native signal rows')
    rows(directory/'exclusions.tsv',('input_row','chromosome','position_1based','reason'))
    return dict(status='JOINT_OUTPUT_STRUCTURE_PASS',input_rows=len(audit),signal_rows=len(signal),
                segment_rows=len(segments),segmented_markers=markers,scientific_validation=False,
                submitted_matrix_rows_verified=input_matrix is not None)
