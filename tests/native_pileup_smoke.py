"""Opt-in actual snp-pileup test; no patient data, no sensitivity claim."""
import argparse
import csv
import gzip
from dataclasses import asdict
import json
from pathlib import Path
import sys
import shutil
import subprocess
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cnv_catalog import merge_depth_grid
from cnv_commands import snp_pileup
from cnv_execution import execute,verify_completed_task,sha
from cnv_pileup import align_catalog
from cnv_prepare_pair import prepare, CATALOG_COLUMNS
from cnv_training import verified_training_summaries
from cnv_allelic import MarkerPolicy
from cnv_territory import Territory
from cnv_train_panel import train


def prepare_native_pair(root, samtools, executable, bam, sam, fasta, vcf, output, snps, paired):
    """Distinct synthetic reference BAM and native task; no mocked verification."""
    reference_sam = root/'normal.sam'
    lines = sam.read_text().splitlines()
    with reference_sam.open('w') as stream:
        for line in lines:
            if line.startswith('@'):
                stream.write(line+'\n')
                continue
            fields = line.split('\t')
            index = int(fields[0][4:])
            fields[0] = 'normal'+str(index)
            sequence = list(fields[9])
            sequence[9] = 'C' if index >= 20 else 'A'
            fields[9] = ''.join(sequence)
            stream.write('\t'.join(fields)+'\n')
    reference_bam = root/'normal.bam'
    subprocess.run([samtools,'view','-b','-o',str(reference_bam),str(reference_sam)],check=True)
    subprocess.run([samtools,'index',str(reference_bam)],check=True)
    subprocess.run([samtools,'faidx',str(fasta)],check=True)
    normal_output = root/'normal-pileup'/'normal.csv.gz'
    policy = dict(mapq=20,baseq=20,pseudo_spacing=10,max_depth=10000,
                  count_orphans=not paired,min_read_count=0)
    command = snp_pileup(reference_bam,vcf,normal_output,executable=executable,**policy)
    execute(command,normal_output.parent,[reference_bam,str(reference_bam)+'.bai',vcf],
            implementation_files=[executable,__file__],timeout_seconds=120)
    def artifact(path):
        return dict(path=str(path),sha256=sha(path))
    target = root/'targets.bed'
    target.write_text('chr1\t0\t200\n')
    mask = root/'mask.bed'
    mask.write_text('chr1\t29\t30\n')
    catalog = root/'catalog.tsv'
    with catalog.open('w',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=CATALOG_COLUMNS,delimiter='\t')
        writer.writeheader()
        writer.writerows(merge_depth_grid(snps,{'chr1':200},10))
    spec = dict(reference_build='GRCh38',capture_id='synthetic',assay_id='synthetic',
                snp_ld_pruned=False,off_target_mode='disabled',reference_fasta=artifact(fasta),
                reference_fai=artifact(str(fasta)+'.fai'),target_bed=artifact(target),
                dense_snp_vcf=artifact(vcf))
    samples = []
    for role, source in [('study',bam),('pon',reference_bam)]:
        samples.append(dict(sample_id=role,person_id=role,role=role,stratum='S',
                            capture_id='synthetic',assay_id='synthetic',reference_sex='Female',
                            bam=artifact(source),bam_index=artifact(str(source)+'.bai')))
    # Held-out row is structural only: no calibration BAM is used or measured.
    samples.append(dict(samples[0],sample_id='calibration',person_id='calibration',role='calibration',
        bam=dict(path=str(root/'heldout-not-executed.bam'),sha256='0'*64),
        bam_index=dict(path=str(root/'heldout-not-executed.bam.bai'),sha256='0'*64)))
    manifest = root/'manifest.json'
    manifest.write_text(json.dumps(dict(schema_version=1,reference_build='GRCh38',strata={'S':spec},samples=samples)))
    config = root/'pair.json'
    config.write_text(json.dumps(dict(sample_id='study',reference_id='pon',catalog=artifact(catalog),
        mask_bed=artifact(mask),executable=executable,pileup_policy=policy,rehash_inputs=True,
        registry={sid:dict(task_directory=str(path.parent),output=str(path))
                  for sid,path in [('study',output),('pon',normal_output)]})))
    destination = root/'pair-prepared'
    receipt = prepare(manifest,config,destination)
    with gzip.open(destination/'pair.evidence.jsonl.gz','rt') as stream:
        rows = {row['Position']:row for row in map(json.loads,stream)}
    assert rows[20]['TUM.RD']==30 and rows[20]['NOR.RD']==20
    assert rows[30]['TUM.DP']==0 and rows[30]['masked']
    assert rows[190]['TUM.DP'] is None and rows[190]['NOR.DP'] is None
    assert all(row['pileup_receipt_bindings_verified'] for row in rows.values())
    with (destination/'pair.rcmat.tsv').open() as stream:
        native = {int(row['Position']):row for row in csv.DictReader(stream,delimiter='\t')}
    assert 190 not in native and native[30]['TUM.DP']=='0' and native[30]['masked']=='TRUE'
    assert receipt['catalog_loci']==20 and not receipt['facets_ch_executed']
    # Independent synthetic training exercise: both generated BAMs now serve
    # as PoN fixtures in a separate manifest, not study or held-out evidence.
    training_manifest=json.loads(manifest.read_text())
    training_manifest['samples'][0].update(sample_id='ponA',person_id='ponA',role='pon')
    training_manifest['samples'].append(dict(samples[0],
        bam=dict(path=str(root/'study-not-executed.bam'),sha256='0'*64),
        bam_index=dict(path=str(root/'study-not-executed.bam.bai'),sha256='0'*64)))
    training=list(verified_training_summaries(training_manifest,'S','Female',
        merge_depth_grid(snps,{'chr1':200},10),
        {sid:dict(task_directory=str(path.parent),output=str(path))
         for sid,path in [('ponA',output),('pon',normal_output)]},
        {'chr1':200},mask=Territory([('chr1',29,30)],{'chr1':200}),
        executable=executable,pileup_policy=policy,
        marker_policy=MarkerPolicy(10,200,.1,.9,2,1,10,2,2,1/6),
        rehash_inputs=True))
    trained={row['position']:row for row in training}
    assert trained[20]['status']=='STABLE_TRAINING_MARKER'
    assert trained[30]['status']=='MASKED_TRAINING_LOCUS'
    assert trained[190]['missing_normals']==2 and trained[190]['correction_bias'] is None
    assert trained[40]['status']=='DEPTH_ONLY' and trained[40]['correction_bias'] is None
    training_manifest_path=root/'training-manifest.json'
    training_manifest_path.write_text(json.dumps(training_manifest))
    training_config=root/'training-config.json'
    training_config.write_text(json.dumps(dict(stratum='S',reference_sex='Female',
        catalog=artifact(catalog),mask_bed=artifact(mask),executable=executable,
        pileup_policy=policy,marker_policy=asdict(MarkerPolicy(10,200,.1,.9,2,1,10,2,2,1/6)),
        rehash_inputs=True,native_chromosome_map={'chr1':1},
        registry={sid:dict(task_directory=str(path.parent),output=str(path))
                  for sid,path in [('ponA',output),('pon',normal_output)]})))
    training_destination=root/'panel'
    training_receipt=train(training_manifest_path,training_config,training_destination)
    assert training_receipt['native_bias_rows']==1
    with (training_destination/'allelic_bias.tsv').open() as stream:
        bias_rows=list(csv.DictReader(stream,delimiter='\t'))
    assert len(bias_rows)==1 and bias_rows[0]['snp_index']=='1_20'
    assert bias_rows[0]['bias_unit']=='ln_alt_ref'
    for name,metadata in training_receipt['outputs'].items():
        assert sha(training_destination/name)==metadata['sha256']
    try:
        train(training_manifest_path,training_config,training_destination)
    except FileExistsError:
        pass
    else:
        raise AssertionError('training output overwritten')
    print(json.dumps(dict(status='REAL_PON_TRAINING_SMOKE_PASS',paired=paired,
        catalog_loci=len(training),training_thresholds='SYNTHETIC_FIXTURE_ONLY',
        sensitivity_calibrated=False)))
    print(json.dumps(dict(status='REAL_PAIR_PREPARATION_SMOKE_PASS',paired=paired,
        catalog_loci=20,native_rows=receipt['native_rows'],
        incomplete_loci=receipt['incomplete_loci'],heldout_calibration_executed=False,
        facets_ch_executed=False,mosaic_sensitivity_measured=False)))


def main(executable, paired=False):
    executable = str(Path(executable).resolve())
    samtools = shutil.which('samtools')
    if not samtools:
        raise RuntimeError('samtools required to generate the synthetic BAM')
    with tempfile.TemporaryDirectory(prefix='cnv-native-pileup-') as temp:
        root = Path(temp)
        fasta = root/'reference.fa'
        fasta.write_text('>chr1\n' + 'A'*200 + '\n')
        vcf = root/'snps.vcf'
        vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=200>\n'
                       '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                       'chr1\t20\t.\tA\tC\t.\t.\t.\n'
                       'chr1\t30\t.\tA\tC\t.\t.\t.\n'
                       'chr1\t190\t.\tA\tC\t.\t.\t.\n')
        bam = root/'sample.bam'
        sam = root/'sample.sam'
        with sam.open('w') as stream:
            stream.write('@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:200\n')
            for index in range(40):
                sequence = list('A'*40)
                if index >= 30:
                    sequence[9] = 'C'  # 1-based position 20.
                sequence[19] = 'G'  # Error-only evidence at 1-based position 30.
                mates = [(99, '=', 11, 40), (147, '=', 11, -40)] if paired else [(0, '*', 0, 0)]
                for flag, mate_chrom, mate_pos, template_length in mates:
                    stream.write(f'read{index}\t{flag}\tchr1\t11\t60\t40M\t'
                                 f'{mate_chrom}\t{mate_pos}\t{template_length}\t'
                                 + ''.join(sequence) + '\t' + 'I'*40 + '\n')
        subprocess.run([samtools, 'view', '-b', '-o', str(bam), str(sam)], check=True)
        subprocess.run([samtools, 'index', str(bam)], check=True)
        output = root/'pileup'/'sample.csv.gz'
        command = snp_pileup(bam, vcf, output, executable=executable, mapq=20, baseq=20,
                             pseudo_spacing=10, max_depth=10000,
                             count_orphans=not paired, min_read_count=0)
        try:
            receipt = execute(command, output.parent, [bam, str(bam)+'.bai', vcf],
                              implementation_files=[executable, __file__], timeout_seconds=120)
        except Exception:
            log = output.parent/'commands.log'
            if log.exists():
                print(log.read_text(errors='replace')[-12000:], file=sys.stderr)
            raise
        assert receipt['status'] == 'NATIVE_TASK_COMPLETE'
        verified = verify_completed_task(command,output.parent,
            {str(p):sha(p) for p in (bam,str(bam)+'.bai',vcf)},rehash_inputs=True)
        assert verified['input_bytes_rehashed'] and verified['outputs_verified']
        snps = [dict(chromosome='chr1', position=p, ref='A', alt='C', kind='SNP')
                for p in (20, 30, 190)]
        rows = list(align_catalog(merge_depth_grid(snps, {'chr1': 200}, 10),
                                  {'sample': output}, {'chr1': 200}))
        by_position = {locus['position']: counts['sample'] for locus, counts in rows}
        assert len(rows) == 20
        assert by_position[20] == (30, 10, 0, 0), by_position
        assert by_position[30] == (0, 0, 40, 0), by_position
        assert by_position[40] == (40, 0, 0, 0), by_position
        assert by_position[190] is None, by_position
        prepare_native_pair(root,samtools,executable,bam,sam,fasta,vcf,output,snps,paired)
        print(json.dumps(dict(status='REAL_SNP_PILEUP_SMOKE_PASS',
                              catalog_loci=len(rows), counts_verified=True,
                              error_only_zero_depth_verified=True,
                              missing_locus_not_zero=True, pseudo_depth_verified=True,
                              native_receipt_verified=True,
                              pair_overlap_tested=paired, mosaic_sensitivity_measured=False)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', required=True)
    args = parser.parse_args()
    main(args.executable)
    main(args.executable, paired=True)
