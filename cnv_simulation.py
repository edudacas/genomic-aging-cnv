"""Ideal diploid two-haplotype mixture expectations for synthetic tests.

Not an empirical WES error model. No sex-chromosome, aneuploid baseline,
subclonal hierarchy, mapping bias, molecule-count or assay validation is implied.
Haplotype A/B are phased copies, not universal REF/ALT assignments.
"""
import math


def diploid_mixture(event_class,cell_fraction,*,affected_haplotype):
    if event_class not in ('LOSS','GAIN','CNLOH') or affected_haplotype not in ('A','B'):
        raise ValueError('explicit supported event and affected haplotype required')
    if (type(cell_fraction) not in (int,float) or not math.isfinite(cell_fraction) or
            not 0<=cell_fraction<=1):
        raise ValueError('cell fraction must be finite in [0,1]')
    # LOSS removes affected haplotype; GAIN duplicates it; CNLOH loses it
    # and duplicates the other haplotype while retaining two total copies.
    altered={'LOSS':(0,1),'GAIN':(2,1),'CNLOH':(0,2)}[event_class]
    if affected_haplotype=='B':
        altered=altered[::-1]
    a=(1-cell_fraction)+cell_fraction*altered[0]
    b=(1-cell_fraction)+cell_fraction*altered[1]
    total=a+b
    return dict(event_class=event_class,cell_fraction=cell_fraction,
        affected_haplotype=affected_haplotype,clone_haplotype_copies=list(altered),
        expected_haplotype_A_copies=a,expected_haplotype_B_copies=b,
        expected_total_copies=total,expected_depth_ratio=total/2,
        expected_log2_ratio=math.log2(total/2),
        expected_haplotype_A_read_fraction=a/total,
        expected_haplotype_B_read_fraction=b/total,
        expected_clone_read_fraction=cell_fraction*sum(altered)/total,
        baseline_haplotype_copies=[1,1],assay_error_model_calibrated=False,
        sensitivity_calibrated=False)
