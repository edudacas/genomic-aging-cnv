"""Deterministic whole-BAM template subsampling, not mosaic-event simulation.

All records with the same read-group/query-name receive the same decision.
No region-specific thinning, allele selection or per-read independent sampling.
UMI-tagged inputs need a separate molecule-aware design and are rejected here.
"""
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from cnv_contract import validate_artifact
from cnv_execution import atomic_json,sha


def template_retained(read_group,query_name,seed,fraction):
    if (type(seed) is not int or seed<0 or type(fraction) not in (float,int) or
            not math.isfinite(fraction) or not 0<fraction<=1):
        raise ValueError('explicit nonnegative seed and fraction in (0,1] required')
    if not isinstance(read_group,str) or not read_group or not isinstance(query_name,str) or not query_name:
        raise ValueError('explicit read group and template name required')
    key=json.dumps([seed,read_group,query_name],separators=(',',':')).encode()
    draw=int.from_bytes(hashlib.sha256(key).digest()[:8],'big')
    cutoff=Fraction(str(fraction))*(1<<64)
    return draw*cutoff.denominator < cutoff.numerator


def downsample_bam(bam_artifact,expected_bam_sample_name,output,*,seed,fraction):
    import pysam
    template_retained('validation','validation',seed,fraction)
    if not isinstance(expected_bam_sample_name,str) or not expected_bam_sample_name:
        raise ValueError('explicit BAM sample name required')
    validate_artifact(bam_artifact,True)
    output=Path(output).resolve()
    output.mkdir(parents=True,exist_ok=False)
    atomic_json(output/'STARTED.json',dict(stage='WHOLE_BAM_TEMPLATE_SUBSAMPLING',
        bam=bam_artifact,expected_bam_sample_name=expected_bam_sample_name,seed=seed,
        fraction=fraction,pysam_version=pysam.__version__,implementation_sha256=sha(__file__)))
    destination=output/'subsampled.bam'
    total=kept=0
    try:
        with pysam.AlignmentFile(bam_artifact['path'],'rb') as source:
            header=source.header.to_dict()
            if header.get('HD',{}).get('SO')!='coordinate':
                raise ValueError('coordinate-sorted BAM declaration required')
            groups=header.get('RG',[])
            if (not groups or any(not g.get('ID') or g.get('SM')!=expected_bam_sample_name for g in groups) or
                    len({g['ID'] for g in groups})!=len(groups)):
                raise ValueError('BAM read-group/sample identity mismatch')
            known={g['ID'] for g in groups}
            with pysam.AlignmentFile(str(destination),'wb',header=source.header) as target:
                for read in source.fetch(until_eof=True):
                    total+=1
                    if any(read.has_tag(tag) for tag in ('RX','MI','UR','UB')):
                        raise ValueError('UMI-tagged input requires molecule-aware subsampling')
                    rg=read.get_tag('RG') if read.has_tag('RG') else None
                    if rg not in known:
                        raise ValueError('read lacks a declared read group')
                    if template_retained(rg,read.query_name,seed,fraction):
                        target.write(read)
                        kept+=1
        pysam.index(str(destination))
        validate_artifact(bam_artifact,True)
        receipt=dict(status='TEMPLATE_SUBSAMPLING_COMPLETE',seed=seed,fraction=fraction,
            bam_sample_name=expected_bam_sample_name,input_alignment_records=total,
            retained_alignment_records=kept,template_count_measured=False,
            outputs={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in
                     (destination,Path(str(destination)+'.bai'))},
            started_sha256=sha(output/'STARTED.json'),
            whole_bam_uniform_sampling=True,mosaic_event_simulated=False,sensitivity_calibrated=False)
        atomic_json(output/'COMPLETE.json',receipt)
        return receipt
    except Exception as error:
        atomic_json(output/'FAILURE.json',dict(error=type(error).__name__,message=str(error)))
        raise
