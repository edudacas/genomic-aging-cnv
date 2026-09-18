"""Prepare private cohort/resource inventory; never launches or approves callers."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import os
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def text(value):
    return '' if value is None else str(value).strip()


def keyed(rows, key):
    result = {}
    for row in rows:
        value = text(row.get(key))
        if not value or value in result:
            raise ValueError('missing or duplicate key: ' + key)
        result[value] = row
    return result


def read_tsv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def reconcile(linkage, workbook, assays, expected):
    links = keyed(linkage, 'individual_id')
    reviews = keyed(workbook, 'individual_id')
    keyed(linkage, 'sample_id')
    keyed(linkage, 'person_id')
    keyed(linkage, 'analysis_id')
    keyed(linkage, 'bam_path')
    if len(links) != expected or set(links) != set(reviews):
        raise ValueError('workbook/linkage membership differs from expected cohort')
    maps = keyed(assays, 'assay_id')
    output = []
    for identity, link in links.items():
        row = reviews[identity]
        for field in ('analysis_id', 'fap'):
            if not text(row.get(field)) or text(row[field]) != text(link[field]):
                raise ValueError('exact identity linkage mismatch: ' + field)
        assay = maps.get(link['assay_id'])
        if assay is None:
            raise ValueError('unmapped assay ID')
        reported = text(row.get('assay_id'))
        if reported and reported not in (link['assay_id'], assay['assay_short_id']):
            raise ValueError('workbook assay conflict')
        reasons = ['QC_RECEIPTS_NOT_IMPORTED', 'BAM_REFERENCE_BED_NOT_VERIFIED',
                   'UMI_STATUS_NOT_VERIFIED', 'CNV_RESOURCES_NOT_FROZEN']
        # Workbook paths may predate the server reconciliation: retain, never rewrite.
        workbook_bam = text(row.get('BAM_path'))
        path_status = 'MATCH' if workbook_bam == link['bam_path'] else 'REVIEW_PATH_DIFFERENCE'
        if path_status != 'MATCH':
            reasons.append(path_status)
        if not text(row.get('batch_id')):
            reasons.append('BATCH_UNKNOWN')
        output.append(dict(sample_id=link['sample_id'], person_id=link['person_id'],
            analysis_id=link['analysis_id'], assay_id=link['assay_id'],
            assay_source=link['assay_source'], assay_name=assay['assay_name_master'],
            proposed_capture_family=assay['proposed_pon_stratum'],
            capture_equivalence_basis=assay['pon_group_equivalence_basis'],
            capture_bed_identity_status=assay['capture_bed_identity_status'],
            batch_id=text(row.get('batch_id')), batch_name=text(row.get('batch_name')),
            sex_reported=text(row.get('sex')), sex_inferred=None,
            sex_inference_status='NOT_IMPORTED', laboratory=None,
            library_protocol=None, umi_status='UNKNOWN',
            mean_target_depth=None, coverage_qc_status='NOT_IMPORTED',
            bam_path=link['bam_path'], bai_path=link['bai_path'],
            workbook_bam_path=workbook_bam, bam_path_comparison=path_status,
            candidate_stratum=link['assay_id'] + ':' + (text(row.get('batch_id')) or 'UNKNOWN'),
            status='HOLD', hold_reasons=reasons))
    return output


def inventory(project, output, expected=1329):
    import openpyxl
    project = Path(project).resolve()
    paths = dict(
        linkage=project/'27_STAGED_CALLERS_V3_5_2_FINAL/study_linkage.PRIVATE.tsv',
        workbook=project.parent.parent/'TP53_LFS_WES_truthset.HUMAN_REVIEW.EN.FULL1329.private.simple.xlsx',
        assays=project/'06_PON_METADATA_20260909_v1/assay_mapping.v1.tsv',
        beds=project/'09_SERVER_INPUTS_v1/beds/manifest.json',
        pon=project/'13_PON_API_QC_v1/pon_index_v1/PON_QC_INDEX.PRIVATE.tsv')
    fingerprints = {k: dict(path=str(p), sha256=sha(p)) for k,p in paths.items()}
    wb = openpyxl.load_workbook(paths['workbook'], read_only=True, data_only=True)
    try:
        values = iter(wb['Review'].values)
        header = next(values)
        reviews = [dict(zip(header,r)) for r in values if any(v is not None for v in r)]
    finally:
        wb.close()
    links = read_tsv(paths['linkage'])
    rows = reconcile(links, reviews, read_tsv(paths['assays']), expected)
    pon = read_tsv(paths['pon'])
    keyed(pon, 'analysis_id')
    # Exact analysis overlap is a check, not proof of independent people.
    overlap = set(r['analysis_id'] for r in rows) & set(r['analysis_id'] for r in pon)
    if overlap:
        raise ValueError('study/PoN analysis overlap')
    bed_manifest = json.loads(paths['beds'].read_text())
    resources = {}
    for name, metadata in bed_manifest['outputs'].items():
        p = paths['beds'].parent/name
        digest = sha(p)
        if digest != metadata['sha256']:
            raise ValueError('BED checksum drift: ' + name)
        resources[name] = dict(path=str(p), sha256=digest, intervals=metadata['intervals'],
            bases=metadata['bases'], status='HASH_VERIFIED_CANDIDATE',
            capture_identity_verified=False, bam_reference_compatibility_verified=False)
    for k,p in paths.items():
        if sha(p) != fingerprints[k]['sha256']:
            raise ValueError('input changed during inventory')
    output = Path(output).resolve()
    if 'PRIVATE' not in output.name or output.is_relative_to(Path(__file__).resolve().parent):
        raise ValueError('use a PRIVATE directory outside the software checkout')
    output.mkdir(parents=True, exist_ok=False, mode=0o700)
    os.chmod(output, 0o700)
    def write(name, value):
        p = output/name
        with p.open('x', encoding='utf-8') as f:
            json.dump(value, f, indent=2, sort_keys=True)
            f.write('\n')
        os.chmod(p, 0o600)
    write('cohort.PRIVATE.json', rows)
    write('sources.PRIVATE.json', fingerprints)
    write('resources.PRIVATE.json', dict(status='HOLD', candidates=resources,
        required_unfrozen=['reference_fasta_and_dictionary', 'capture_identity',
            'antitarget_or_target_only_policy', 'problem_region_masks',
            'dense_snp_catalog', 'pseudo_depth_grid', 'CNV_PoN_reference',
            'stable_markers_and_allelic_bias', 'X_Y_PAR_policy', 'calibration_holdout'],
        pon_existing_donors=len(pon), pon_analysis_overlap=0,
        existing_snv_pon_is_cnv_reference=False,
        independent_people_verification='REQUIRES_EXACT_PERSON_LINKAGE'))
    summary = dict(status='INVENTORY_PREPARED_WITH_HOLDS', selected_samples=len(rows),
        assay_counts=dict(Counter(r['assay_name'] for r in rows)),
        capture_family_counts=dict(Counter(r['proposed_capture_family'] for r in rows)),
        sex_reported_counts=dict(Counter(r['sex_reported'] or 'UNKNOWN' for r in rows)),
        batch_missing=sum(not r['batch_id'] for r in rows),
        candidate_assay_batch_strata=len(set(r['candidate_stratum'] for r in rows)),
        workbook_path_differences=sum(r['bam_path_comparison'] != 'MATCH' for r in rows),
        hold_counts=dict(Counter(reason for r in rows for reason in r['hold_reasons'])),
        existing_pon_donors=len(pon), pon_analysis_overlap=0,
        resources_frozen=False, cohort_executed=False)
    write('summary.json', summary)
    write('manifest.json', {p.name: sha(p) for p in sorted(output.iterdir()) if p.is_file()})
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--expected', type=int, default=1329)
    args = parser.parse_args()
    print(json.dumps(inventory(args.project, args.output, args.expected), indent=2))
