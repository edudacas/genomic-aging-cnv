"""Opt-in real CNVkit coverage/reference/fix smoke test on generated data.

Run from module root with a Python environment containing CNVkit and pysam.
Optional --cbs exercises DNAcopy on a denser fixture. Neither mode calibrates
mosaic detection or exercises FACETS-CH.
"""
import argparse
import csv
import json
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
import tempfile

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cnv_commands import coverage, pooled_reference, fix, segment
from cnv_execution import execute as native_execute
from cnv_execution import sha
from cnv_reference import build_reference
from cnv_normalize import normalize_sample
from cnv_r_library import snapshot_library, verify_library
from cnv_cbs import run_cbs


def execute(command, directory, *args, **kwargs):
    try:
        return native_execute(command, directory, *args, **kwargs)
    except Exception:
        log = Path(directory)/'commands.log'
        if log.exists():
            print(log.read_text(errors='replace')[-12000:], file=sys.stderr)
        raise


def bam(path, target_copies, off_target=False, bins=20):
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
              'SQ': [{'SN': 'chr1', 'LN': bins*2000}]}
    with pysam.AlignmentFile(str(path), 'wb', header=header) as stream:
        for bin_number in range(bins):
            start = bin_number * 2000
            copies = target_copies if bin_number < bins//4 else 4
            for position in range(start, start + (2000 if off_target else 1000), 100):
                for copy in range(copies):
                    read = pysam.AlignedSegment()
                    read.query_name = f'read_{position}_{copy}'
                    read.reference_id = 0
                    read.reference_start = position
                    read.query_sequence = 'ACGT' * 25
                    read.query_qualities = pysam.qualitystring_to_array('I' * 100)
                    read.flag = 0
                    read.mapping_quality = 60
                    read.cigarstring = '100M'
                    stream.write(read)
    pysam.index(str(path))


def run_mode(off_target, cbs=False, r_library=None):
    bins = 100 if cbs else 20
    loss_end = bins//4 * 2000
    executable = shutil.which('cnvkit.py')
    if not executable:
        raise RuntimeError('CNVkit not available')
    with tempfile.TemporaryDirectory(prefix='cnvkit-real-smoke-') as tmp:
        root = Path(tmp)
        fasta = root/'reference.fa'
        fasta.write_text('>chr1\n' + 'ACGT'*(bins*500) + '\n')
        pysam.faidx(str(fasta))
        bed = root/'targets.bed'
        bed.write_text(''.join(f'chr1\t{i*2000}\t{i*2000+1000}\tG{i}\n' for i in range(bins)))
        anti_bed = root/'antitargets.bed'
        anti_bed.write_text(''.join(
            f'chr1\t{i*2000+1000}\t{i*2000+2000}\tAntitarget\n' for i in range(bins)))
        implementations = [executable, __file__]
        tables = {}
        for sample, copies in [('normal', 4), ('study', 2)]:
            native_bam = root/(sample+'.bam')
            bam(native_bam, copies, off_target, bins)
            task = root/sample
            output = task/(sample+'.targetcoverage.cnn')
            command = coverage(native_bam, bed, output, executable=executable,
                               processes=1, min_mapq=20)
            receipt = execute(command, task, [native_bam, str(native_bam)+'.bai', bed],
                              implementation_files=implementations, timeout_seconds=120)
            assert receipt['status'] == 'NATIVE_TASK_COMPLETE'
            assert execute(command, task, [native_bam, str(native_bam)+'.bai', bed],
                           implementation_files=implementations, timeout_seconds=120)['reused']
            tables[sample] = output
            if off_target:
                anti_task = root/(sample+'_anti')
                anti_output = anti_task/(sample+'.antitargetcoverage.cnn')
                anti_command = coverage(native_bam, anti_bed, anti_output,
                                        executable=executable, processes=1, min_mapq=20)
                execute(anti_command, anti_task, [native_bam, str(native_bam)+'.bai', anti_bed],
                        implementation_files=implementations, timeout_seconds=120)
                tables[sample+'_anti'] = anti_output
        reference = root/'reference'/'reference.cnn'
        def declared(path):
            return {'path': str(path), 'sha256': sha(path)}
        # This synthetic one-normal reference is an interface test, not a
        # sufficient clinical PoN. Calibration remains a disjoint participant.
        calibration_bam = root/'calibration.bam'
        bam(calibration_bam, 4, off_target, bins)
        vcf = root/'dense.vcf'
        vcf.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                       'chr1\t1\t.\tA\tC\t.\t.\t.\n')
        spec = dict(capture_id='synthetic', assay_id='synthetic', reference_build='GRCh38',
                    snp_ld_pruned=False, off_target_mode='supported' if off_target else 'disabled',
                    reference_fasta=declared(fasta), reference_fai=declared(str(fasta)+'.fai'),
                    target_bed=declared(bed), dense_snp_vcf=declared(vcf))
        if off_target:
            spec['antitarget_bed'] = declared(anti_bed)
        samples = [dict(sample_id=name, person_id=name, role=role, stratum='synthetic',
                        capture_id='synthetic', assay_id='synthetic',
                        bam=declared(root/(name+'.bam')),
                        bam_index=declared(root/(name+'.bam.bai')))
                   for name, role in [('normal', 'pon'), ('study', 'study'),
                                      ('calibration', 'calibration')]]
        manifest = root/'manifest.json'
        manifest.write_text(json.dumps(dict(schema_version=1, reference_build='GRCh38',
                                            strata={'synthetic': spec}, samples=samples)))
        registry = root/'coverage_registry.json'
        entries = {'target': str(tables['normal'].parent)}
        if off_target:
            entries['antitarget'] = str(tables['normal_anti'].parent)
        registry.write_text(json.dumps({'normal': entries}))
        try:
            reference_receipt = build_reference(
                manifest, registry, 'synthetic', reference.parent, executable=executable,
                processes=1, min_mapq=20, minimum_pon_samples=1,
                implementation_files=implementations, timeout_seconds=120)
        except Exception:
            log = reference.parent/'commands.log'
            if log.exists():
                print(log.read_text()[-12000:], file=sys.stderr)
            raise
        assert reference_receipt['reference_members'] == ['normal']
        sample_registry = root/'sample_registry.json'
        sample_entries = {'target':str(tables['study'].parent)}
        if off_target:
            sample_entries['antitarget'] = str(tables['study_anti'].parent)
        sample_registry.write_text(json.dumps({'study':sample_entries}))
        try:
            normalized = normalize_sample(manifest, 'study', sample_registry, registry,
                reference.parent, root/'fix', executable=executable, processes=1,
                min_mapq=20, minimum_pon_samples=1, implementation_files=implementations,
                timeout_seconds=120)
        except Exception:
            log = root/'fix/commands.log'
            if log.exists():
                print(log.read_text()[-12000:],file=sys.stderr)
            raise
        assert normalized['reference_provenance']['reference_members'] == ['normal']
        assert not normalized['somatic_origin_classified']
        cnr = Path(normalized['normalized_table'])
        with tables['study'].open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        assert len(rows) == bins
        assert all(float(row['depth']) == 2 for row in rows[:bins//4])
        assert all(float(row['depth']) == 4 for row in rows[bins//4:])
        with cnr.open() as stream:
            fixed = list(csv.DictReader(stream, delimiter='\t'))
        expected_bins = bins * (2 if off_target else 1)
        assert len(fixed) == expected_bins
        loss_log2 = statistics.median(float(row['log2']) for row in fixed if int(row['start']) < loss_end)
        neutral_log2 = statistics.median(float(row['log2']) for row in fixed if int(row['start']) >= loss_end)
        assert loss_log2 < neutral_log2 - 0.5, (loss_log2, neutral_log2)
        if cbs:
            if r_library is None:
                raise ValueError('--r-library is required for a fingerprinted CBS test')
            library = snapshot_library(r_library,required_packages=['DNAcopy'])
            snapshot = root/'cbs-library.json'
            snapshot.write_text(json.dumps(library))
            rscript = shutil.which('Rscript')
            if not rscript:
                raise RuntimeError('Rscript required')
            rscript = str(Path(rscript).resolve())
            cbs_config = root/'cbs-normalization.json'
            cbs_config.write_text(json.dumps(dict(manifest_path=str(manifest),sample_id='study',
                sample_registry_path=str(sample_registry),pon_registry_path=str(registry),
                reference_directory=str(reference.parent),directory=str(cnr.parent),
                executable=executable,processes=1,min_mapq=20,minimum_pon_samples=1)))
            rscript_artifact = root/'rscript.json'
            rscript_artifact.write_text(json.dumps(declared(rscript)))
            cli_command = [sys.executable,str(Path(__file__).resolve().parents[1]/'cnv_cli.py'),
                'run-cbs','--config',str(cbs_config),'--library-snapshot',str(snapshot),
                '--rscript-artifact',str(rscript_artifact),'--output',str(root/'cbs'),
                '--processes','1','--threshold','.0001','--timeout-seconds','120','--execute']
            for implementation in implementations:
                cli_command.extend(['--implementation-file',str(implementation)])
            launched = subprocess.run(cli_command,check=True,text=True,stdout=subprocess.PIPE)
            cbs_receipt = json.loads(launched.stdout)
            assert cbs_receipt['status']=='CNVKIT_CBS_COMPLETE'
            assert cbs_receipt['normalized_provenance_verified']
            cns = root/'cbs/native/sample.cns'
            verify_library(library)
            with cns.open() as stream:
                segments = list(csv.DictReader(stream, delimiter='\t'))
            assert len(segments) >= 2, segments
            loss = [row for row in segments if int(row['start']) < loss_end
                    and int(row['end']) <= loss_end]
            neutral = [row for row in segments if int(row['start']) >= loss_end]
            assert loss and neutral, segments
            assert all(float(row['log2']) < -0.5 for row in loss), segments
            assert all(abs(float(row['log2'])) < 0.25 for row in neutral), segments
        print(json.dumps(dict(status='REAL_CNVKIT_SMOKE_PASS', bins=expected_bins,
                              off_target_tested=off_target,
                              stages=['coverage', 'reference', 'fix'] + (['cbs'] if cbs else []),
                              loss_log2_median=loss_log2,
                              neutral_log2_median=neutral_log2,
                              segmentation_tested=cbs, facets_ch_tested=False,
                              normalization_provenance_verified=True,
                              explicit_R_runtime_fingerprinted=cbs,
                              cbs_stage_orchestrator_tested=cbs,
                              cbs_cli_tested=cbs,
                              mosaic_sensitivity_measured=False)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cbs', action='store_true', help='Require native R/DNAcopy CBS test')
    parser.add_argument('--r-library',type=Path,help='Isolated R library containing DNAcopy')
    args = parser.parse_args()
    if args.cbs and args.r_library is None:
        parser.error('--cbs requires --r-library')
    run_mode(False, args.cbs, args.r_library)
    run_mode(True, args.cbs, args.r_library)
