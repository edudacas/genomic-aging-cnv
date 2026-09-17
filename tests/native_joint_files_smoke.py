"""Opt-in real R entrypoint/receipt test; synthetic counts, not clinical validation.

Run with --library pointing to a complete isolated FACETS-CH library and
--output pointing to a new directory. Leaves auditable native artifacts.
This does not test cohort linkage or upstream BAM-to-matrix derivation.
"""
import argparse
import csv
import json
from pathlib import Path
import shutil
import sys

SOURCE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SOURCE))
from cnv_commands import joint_fit
from cnv_execution import atomic_json, execute, sha, verify_completed_task
from cnv_joint_output import validate_joint_output
from cnv_r_library import snapshot_library, verify_library


def write_tsv(path, header, rows):
    with path.open('x', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--library', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--rscript', default=shutil.which('Rscript'))
    args = parser.parse_args()
    if not args.rscript:
        parser.error('Rscript required')
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    library = snapshot_library(args.library, required_packages=['facetsCH', 'pctGCdata'])
    atomic_json(root/'library_snapshot.json', library)
    matrix, bias, config, mapping = [root/name for name in
        ('matrix.tsv', 'bias.tsv', 'config.tsv', 'mapping.tsv')]
    rows = []
    for i in range(1, 601):
        normal = 120+i % 19
        loss = i <= 150
        depth = round(normal*(0.75 if loss else 1)*(1+(i % 7-3)/100))
        rows.append(('chr1', 1000000+i*1000, normal, round(normal/2), depth,
                     round(depth*(2/3 if loss else 0.5)), 'SNP', 'FALSE', i))
    write_tsv(matrix, ['Chromosome', 'Position', 'NOR.DP', 'NOR.RD', 'TUM.DP',
                      'TUM.RD', 'catalog_kind', 'masked', 'catalog_row'], rows)
    write_tsv(bias, ['snp_index', 'correction_bias', 'bias_unit', 'allelic_pseudocount'],
              [(f'1_{r[1]}', 0, 'ln_alt_ref', 1/6) for r in rows])
    # Explicit synthetic-test parameters; not study-calibrated defaults.
    policy = dict(reference_build='GRCh38', ndepth=20, ndepthmax=2000,
        het_threshold=0.1, snp_neighborhood=0, initial_cval=25, final_cval=50,
        gamma=0.25, min_nhet=15, seed=0, allelic_pseudocount=1/6,
        sample_het_min_depth=50, het_scale='TRUE', delta_cn=0)
    write_tsv(config, ['key', 'value'], sorted(policy.items()))
    write_tsv(mapping, ['source', 'native'], [('chr1', '1')])
    native = root/'native'
    result = native/'result'
    command = joint_fit(matrix, bias, config, mapping, result,
        rscript=str(Path(args.rscript).resolve()), wrapper=SOURCE/'facets_run.R',
        library=args.library.resolve())
    inputs = [matrix, bias, config, mapping, root/'library_snapshot.json']
    implementation = [Path(args.rscript).resolve(), Path(__file__).resolve(),
        *SOURCE.glob('facets_*.R'), *SOURCE.glob('cnv_*.py'),
        *(Path(library['directory'])/name for name in library['files'])]
    kwargs = dict(implementation_files=implementation, timeout_seconds=120)
    first = execute(command, native, inputs, **kwargs)
    assert not first['reused']
    validation = validate_joint_output(result, 600, input_matrix=matrix)
    with (result/'segments.tsv').open() as stream:
        segments = list(csv.DictReader(stream, delimiter='\t'))
    assert len(segments) >= 2
    assert sum(int(s['n_markers']) for s in segments) == 600
    assert (result/'COMPLETE.txt').read_text().splitlines() == [
        'JOINT_SEGMENTATION_CANDIDATES', 'sensitivity_calibrated=FALSE',
        'somatic_origin_classified=FALSE', 'copy_neutral_loh_confirmed=FALSE']
    verify_library(library)
    verified = verify_completed_task(command, native,
        {str(p): sha(p) for p in inputs}, rehash_inputs=True)
    resumed = execute(command, native, inputs, **kwargs)
    assert resumed['reused']
    receipt = dict(status='NATIVE_JOINT_FILE_SMOKE_PASS', native_validation=validation,
        native_task_verification=verified, resume_verified=True,
        segments=len(segments), loci=600, library_snapshot_sha256=sha(inputs[-1]),
        synthetic_counts_only=True, bam_derivation_tested=False,
        sensitivity_calibrated=False, somatic_origin_classified=False,
        full_pipeline_validated=False)
    atomic_json(root/'SMOKE_COMPLETE.json', receipt)
    print(json.dumps(receipt, sort_keys=True))


if __name__ == '__main__':
    main()
