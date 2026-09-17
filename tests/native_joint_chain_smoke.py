"""Synthetic BAM -> native pileups -> PoN/pair preparation -> native FACETS.

Opt-in integration fixture, not a study reference or sensitivity experiment.
Requires --pileup, --library and a fresh --output directory. Preserves artifacts.
"""
import argparse
import csv
import gzip
from dataclasses import asdict
import json
from pathlib import Path
import shutil
import subprocess
import sys

SOURCE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SOURCE))
from cnv_allelic import MarkerPolicy
from cnv_commands import snp_pileup
from cnv_execution import atomic_json, execute, sha
from cnv_prepare_pair import prepare, CATALOG_COLUMNS
from cnv_train_panel import train
from cnv_prepare_joint import prepare_joint
from cnv_run_joint import run_joint
from cnv_r_library import snapshot_library
from cnv_commands import coverage
from cnv_reference import build_reference
from cnv_normalize import normalize_sample
from cnv_cbs import run_cbs
from cnv_reconcile import read_segments, reconcile_segments
from cnv_compare import compare_callers
from cnv_simulation import diploid_mixture


def paired_cnvkit(root, manifest, samples, target, fasta, library, executable, rscript, length, facets,event_class):
    """Both callers see identical synthetic study/PoN BAM identities."""
    tasks = {}
    for sample in samples:
        if sample['role']=='calibration':
            continue
        sid = sample['sample_id']
        directory = root/'cnvkit'/('coverage-'+sid)
        output = directory/'target.cnn'
        command = coverage(sample['bam']['path'],target,output,executable=executable,
                           processes=1,min_mapq=20)
        execute(command,directory,[sample['bam']['path'],sample['bam_index']['path'],target],
                implementation_files=[executable,__file__],timeout_seconds=120)
        tasks[sid] = {'target':str(directory)}
    normal_registry, sample_registry = root/'cnvkit/pon.json',root/'cnvkit/study.json'
    atomic_json(normal_registry,{sid:tasks[sid] for sid in ('pon1','pon2')})
    atomic_json(sample_registry,{'study':tasks['study']})
    options = dict(executable=executable,processes=1,min_mapq=20,minimum_pon_samples=2,
                   implementation_files=[executable,__file__],timeout_seconds=120)
    reference = root/'cnvkit/reference'
    build_reference(manifest,normal_registry,'S',reference,**options)
    normalize_sample(manifest,'study',sample_registry,normal_registry,reference,root/'cnvkit/fix',**options)
    config, runtime = root/'cnvkit/cbs.json',root/'cnvkit/runtime.json'
    atomic_json(config,dict(manifest_path=str(manifest),sample_id='study',
        sample_registry_path=str(sample_registry),pon_registry_path=str(normal_registry),
        reference_directory=str(reference),directory=str(root/'cnvkit/fix'),
        executable=executable,processes=1,min_mapq=20,minimum_pon_samples=2))
    atomic_json(runtime,snapshot_library(library,required_packages=['DNAcopy']))
    cbs = run_cbs(config,runtime,artifact(Path(rscript).resolve()),root/'cnvkit/cbs',
        processes=1,threshold=.0001,implementation_files=[executable,__file__],timeout_seconds=120)
    for receipt in (cbs,facets):
        for path,metadata in receipt['outputs'].items():
            assert sha(path)==metadata['sha256']
    cnvkit_path = root/'cnvkit/cbs/native/sample.cns'
    facets_path = root/'fit/native/result/segments.tsv'
    comparison = compare_callers(artifact(manifest),'study',
        artifact(root/'cnvkit/cbs/COMPLETE.json'),artifact(root/'fit/COMPLETE.json'),
        root/'comparison',logr_direction_threshold=.1)
    with (root/'comparison/comparison.jsonl').open() as stream:
        rows = list(map(json.loads,stream))
    assert comparison['comparison_intervals']==len(rows)
    # Direction only; boundaries may differ and gaps are not evaluable territory.
    loss = [r for r in rows if r['comparison_status']=='CONCORDANT_DIRECTION'
            and r['cnvkit_source']['direction']=='LOSS_DIRECTION']
    neutral = [r for r in rows if r['comparison_status']=='BOTH_WITHIN_DECLARED_LOGR_BAND']
    gain = [r for r in rows if r['comparison_status']=='CONCORDANT_DIRECTION'
            and r['cnvkit_source']['direction']=='GAIN_DIRECTION']
    assert neutral
    if event_class=='LOSS':
        assert loss
    elif event_class=='GAIN':
        assert gain
    assert all(r['event_class']=='UNCLASSIFIED' and r['evaluable_bp'] is None for r in rows)
    with (root/'paired-reconciliation.jsonl').open('x') as stream:
        for row in rows:
            stream.write(json.dumps(row,sort_keys=True)+'\n')
    return dict(status='SAME_SYNTHETIC_BAMS_PAIRED_CALLERS_PASS',
        comparison_intervals=len(rows),concordant_loss_intervals=len(loss),
        concordant_gain_intervals=len(gain),simulated_event_class=event_class,
        within_band_intervals=len(neutral),cbs_segments=cbs['segment_rows'],
        reconciliation_sha256=sha(root/'paired-reconciliation.jsonl'),
        biological_events_validated=False,sensitivity_calibrated=False)


def artifact(path):
    return dict(path=str(path), sha256=sha(path))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pileup', type=Path, required=True)
    p.add_argument('--library', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--cnvkit',type=Path,help='Optional same-BAM CNVkit/CBS reconciliation test')
    p.add_argument('--event-class',choices=('LOSS','GAIN','CNLOH'),default='LOSS')
    a = p.parse_args()
    samtools, rscript = shutil.which('samtools'), shutil.which('Rscript')
    if not samtools or not rscript:
        p.error('samtools and Rscript required')
    executable = str(a.pileup.resolve())
    root = a.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    length = 1700000
    positions = [1000000+i*1000 for i in range(1, 601)]
    mixture=diploid_mixture(a.event_class,.5,affected_haplotype='B')
    atomic_json(root/'synthetic-mixture.json',mixture)
    absent_position = 1650000
    masked_position = positions[299]
    fasta, vcf = root/'synthetic.fa', root/'snps.vcf'
    # Artificial sequence; hg38 GC lookup exercises the interface only.
    sequence = list(('ACGT'*((length+3)//4))[:length])
    for pos in positions+[absent_position]:
        sequence[pos-1] = 'A'
    reference_sequence = ''.join(sequence)
    fasta.write_text('>chr1\n'+reference_sequence+'\n')
    subprocess.run([samtools, 'faidx', str(fasta)], check=True)
    with vcf.open('x') as f:
        f.write(f'##fileformat=VCFv4.2\n##contig=<ID=chr1,length={length}>\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for pos in positions+[absent_position]:
            f.write(f'chr1\t{pos}\t.\tA\tC\t.\t.\t.\n')
    catalog = root/'catalog.tsv'
    with catalog.open('x', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(CATALOG_COLUMNS)
        w.writerows(('chr1', pos, 'A', 'C', 'SNP') for pos in positions+[absent_position])
    target, mask = root/'target.bed', root/'mask.bed'
    target.write_text(''.join(f'chr1\t{pos-50}\t{pos+50}\tBIN{i}\n'
                              for i,pos in enumerate(positions,1)))
    mask.write_text(f'chr1\t{masked_position-1}\t{masked_position}\n')
    policy = dict(mapq=20, baseq=20, pseudo_spacing=2000000, max_depth=10000,
                  count_orphans=True, min_read_count=0)
    registry, samples = {}, []
    for sid, role in [('study', 'study'), ('pon1', 'pon'), ('pon2', 'pon')]:
        sam, bam = root/(sid+'.sam'), root/(sid+'.bam')
        with sam.open('x') as f:
            f.write(f'@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{length}\n')
            f.write(f'@RG\tID:{sid}\tSM:{sid}\n')
            for i, pos in enumerate(positions, 1):
                normal = 120+2*(i % 19)
                affected = role == 'study' and i <= 150
                depth = (round(normal*(mixture['expected_depth_ratio'] if affected else 1)*(1+(i % 7-3)/100))
                         if role == 'study' else normal)
                ref_count = round(depth*(mixture['expected_haplotype_A_read_fraction'] if affected else 0.5))
                for j in range(depth):
                    read = list(reference_sequence[pos-50:pos+50])
                    read[49] = 'A' if j < ref_count else 'C'
                    seq = ''.join(read)
                    f.write(f'{sid}_{i}_{j}\t0\tchr1\t{pos-49}\t60\t100M\t*\t0\t0\t'
                            f'{seq}\t'+ 'I'*100+f'\tRG:Z:{sid}\n')
        subprocess.run([samtools, 'view', '-b', '-o', str(bam), str(sam)], check=True)
        subprocess.run([samtools, 'index', str(bam)], check=True)
        output = root/('pileup-'+sid)/'counts.csv.gz'
        command = snp_pileup(bam, vcf, output, executable=executable, **policy)
        execute(command, output.parent, [bam, str(bam)+'.bai', vcf],
                implementation_files=[executable, __file__], timeout_seconds=120)
        registry[sid] = dict(task_directory=str(output.parent), output=str(output))
        samples.append(dict(sample_id=sid, person_id=sid, role=role, stratum='S',
            capture_id='synthetic', assay_id='synthetic', reference_sex='Female',
            bam_sample_name=sid, bam=artifact(bam), bam_index=artifact(str(bam)+'.bai')))
    samples.append(dict(samples[0], sample_id='heldout', person_id='heldout', role='calibration',
        bam=dict(path=str(root/'not-executed.bam'), sha256='0'*64),
        bam_index=dict(path=str(root/'not-executed.bam.bai'), sha256='0'*64)))
    spec = dict(reference_build='GRCh38', capture_id='synthetic', assay_id='synthetic',
        snp_ld_pruned=False, off_target_mode='disabled', reference_fasta=artifact(fasta),
        reference_fai=artifact(str(fasta)+'.fai'), target_bed=artifact(target), dense_snp_vcf=artifact(vcf))
    manifest = root/'manifest.json'
    atomic_json(manifest, dict(schema_version=1, reference_build='GRCh38', strata={'S':spec}, samples=samples))
    common = dict(catalog=artifact(catalog), mask_bed=artifact(mask), executable=executable,
                  pileup_policy=policy, rehash_inputs=True)
    pair_config, panel_config, fit_config = root/'pair.json', root/'panel.json', root/'fit.json'
    atomic_json(pair_config, dict(common, sample_id='study', reference_id='pon1',
                                 registry={k:registry[k] for k in ('study','pon1')}))
    atomic_json(panel_config, dict(common, stratum='S', reference_sex='Female',
        native_chromosome_map={'chr1':1}, registry={k:registry[k] for k in ('pon1','pon2')},
        marker_policy=asdict(MarkerPolicy(10,200,.1,.9,2,1,10,2,2,1/6))))
    pair = prepare(manifest, pair_config, root/'pair')
    panel = train(manifest, panel_config, root/'panel')
    assert pair['native_rows'] == 600 and panel['native_bias_rows'] == 599
    assert pair['catalog_loci'] == 601 and pair['incomplete_loci'] == 1
    with gzip.open(root/'pair/pair.evidence.jsonl.gz','rt') as f:
        ledger = {r['Position']:r for r in map(json.loads,f)}
    assert ledger[absent_position]['TUM.DP'] is None
    assert ledger[absent_position]['NOR.DP'] is None
    assert ledger[masked_position]['masked'] is True
    atomic_json(fit_config, dict(ndepth=20, ndepthmax=2000, het_threshold=.1,
        snp_neighborhood=0, initial_cval=25, final_cval=50, gamma=.25, min_nhet=15,
        seed=0, sample_het_min_depth=50, delta_cn=0, het_scale=True))
    prepare_joint(manifest, root/'pair', pair_config, root/'panel', panel_config,
                  fit_config, root/'prepared')
    runtime = root/'runtime.json'
    atomic_json(runtime, snapshot_library(a.library, required_packages=['facetsCH','pctGCdata']))
    result = run_joint(root/'prepared', runtime, artifact(Path(rscript).resolve()),
                       root/'fit', timeout_seconds=120)
    assert result['output_validation']['segmented_markers'] == 599
    with (root/'fit/native/result/signal.tsv').open() as f:
        signal_positions = {int(r['maploc']) for r in csv.DictReader(f,delimiter='\t')}
    assert masked_position not in signal_positions and absent_position not in signal_positions
    assert result['output_validation']['segment_rows'] >= 2
    receipt = dict(status='SYNTHETIC_BAM_TO_JOINT_SEGMENTATION_PASS',
        native_pileup_samples=3, independent_synthetic_pon_bams=2, catalog_loci=601,
        numeric_loci=600, segmented_loci=599, missing_locus_preserved_as_null=True,
        masked_locus_excluded_from_native_signal=True,
        native_output_validation=result['output_validation'],
        native_completion_sha256=sha(root/'fit/COMPLETE.json'),
        all_preparation_stages_executed=True, native_engine_mocked=False,
        synthetic_parameters_only=True, synthetic_reference_sequence=True,
        synthetic_mixture=mixture,
        heldout_calibration_executed=False, sensitivity_calibrated=False,
        cnvkit_branch_tested=False, full_pipeline_validated=False)
    if a.cnvkit:
        receipt['paired_callers'] = paired_cnvkit(root,manifest,samples,target,fasta,a.library,
            str(a.cnvkit.resolve()),rscript,length,result,a.event_class)
        receipt['cnvkit_branch_tested'] = True
    atomic_json(root/'SMOKE_COMPLETE.json', receipt)
    print(json.dumps(receipt, sort_keys=True))


if __name__ == '__main__':
    main()
