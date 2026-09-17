from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from cnv_phase_sam import (query_positions,extract_snp_observations,phase_subset_command,
                           validate_loci,phase_reference_command,verify_reference_bases)
from cnv_phasing import summarize_fragment_evidence,PhaseReadPolicy
from cnv_execution import execute,sha
from cnv_phase_verified import verified_phase_evidence
from test_contract import fixture


class PhaseSamTests(unittest.TestCase):
    def test_cigar_offsets_and_absent_deleted_bases(self):
        positions=query_positions('2S3M2I2M2D3M1H',10,12,[10,12,13,14,15,16,17,19])
        self.assertEqual(positions,{10:2,12:4,13:7,14:8,17:9,19:11})
        with self.assertRaises(ValueError):
            query_positions('2M1B',1,2,[1])

    def test_complex_allele_rejected(self):
        loci=[dict(locus_id='a',chromosome='chr1',position=20,ref='A',alt='AC'),
              dict(locus_id='b',chromosome='chr1',position=30,ref='A',alt='G')]
        with self.assertRaisesRegex(ValueError,'biallelic SNP'):
            validate_loci(loci)

    @unittest.skipUnless(shutil.which('samtools'),'samtools unavailable')
    def test_real_bam_subset_to_fragment_evidence(self):
        samtools=shutil.which('samtools')
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            sam=root/'reads.sam'; bam=root/'reads.bam'
            seq=list('A'*40); seq[9]='C'; seq[19]='G'
            with sam.open('w') as stream:
                stream.write('@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:200\n@RG\tID:r\tSM:s\n')
                for flag in (99,147):
                    stream.write('fragment\t'+str(flag)+'\tchr1\t11\t60\t40M\t=\t11\t0\t'+''.join(seq)+'\t'+'I'*40+'\tRG:Z:r\n')
            subprocess.run([samtools,'view','-b','-o',str(bam),str(sam)],check=True)
            subprocess.run([samtools,'index',str(bam)],check=True)
            loci=[dict(locus_id='a',chromosome='chr1',position=20,ref='A',alt='C'),
                  dict(locus_id='b',chromosome='chr1',position=30,ref='A',alt='G')]
            fasta=root/'reference.fa'
            fasta.write_text('>chr1\n'+'A'*200+'\n')
            subprocess.run([samtools,'faidx',str(fasta)],check=True)
            refout=root/'reference-subset'/'bases.fa'
            refcmd=phase_reference_command(fasta,refout,loci,samtools=samtools)
            execute(refcmd,refout.parent,[fasta,str(fasta)+'.fai'],
                    implementation_files=[samtools,__file__],timeout_seconds=30)
            self.assertEqual(verify_reference_bases(refout,loci)['loci'],2)
            with self.assertRaisesRegex(ValueError,'REF differs'):
                verify_reference_bases(refout,[dict(loci[0],ref='G'),loci[1]])
            output=root/'subset'/'reads.sam'
            cmd=phase_subset_command(bam,output,loci,samtools=samtools)
            receipt=execute(cmd,output.parent,[bam,str(bam)+'.bai'],
                implementation_files=[samtools,__file__],timeout_seconds=30)
            self.assertEqual(receipt['status'],'NATIVE_TASK_COMPLETE')
            observations,diagnostics=extract_snp_observations(output,'study-id',loci,expected_bam_sample_name='s')
            self.assertEqual(diagnostics['alignment_records'],2)
            self.assertEqual(len(observations),4)
            self.assertEqual({o['sample_id'] for o in observations},{'study-id'})
            result=summarize_fragment_evidence(observations,'study-id',['a','b'],
                policy=PhaseReadPolicy(20,20,4|256|2048|1024))
            self.assertEqual(result['haplotype_fragment_counts']['ALT_ALT'],1)
            self.assertEqual(result['linked_fragments'],1)
            with self.assertRaisesRegex(ValueError,'identity mismatch'):
                extract_snp_observations(output,'study-id',loci,expected_bam_sample_name='wrong-sm')
            with self.assertRaisesRegex(ValueError,'sequence dictionary'):
                extract_snp_observations(output,'study-id',loci,expected_bam_sample_name='s',
                                         expected_reference_lengths={'chr1':201})
            manifest=fixture()
            def artifact(path):
                return dict(path=str(path),sha256=sha(path))
            study=manifest['samples'][0]
            study.update(bam=artifact(bam),bam_index=artifact(str(bam)+'.bai'),bam_sample_name='s')
            manifest['strata']['S'].update(reference_fasta=artifact(fasta),
                                           reference_fai=artifact(str(fasta)+'.fai'))
            sam_task=dict(task_directory=str(output.parent),output=str(output))
            reference_task=dict(task_directory=str(refout.parent),output=str(refout))
            policy=PhaseReadPolicy(20,20,4|256|2048|1024)
            linked=verified_phase_evidence(manifest,'study',loci,sam_task,reference_task,
                samtools=samtools,policy=policy,rehash_inputs=True)
            self.assertTrue(linked['summary']['bam_provenance_verified'])
            self.assertFalse(linked['summary']['clinical_identity_verified'])
            self.assertEqual(linked['summary']['haplotype_fragment_counts']['ALT_ALT'],1)
            # A different requested locus cannot reuse the same extracted SAM.
            with self.assertRaisesRegex(ValueError,'command or filter'):
                verified_phase_evidence(manifest,'study',[dict(loci[0],position=21),loci[1]],
                    sam_task,reference_task,samtools=samtools,policy=policy,rehash_inputs=True)
            study['bam']['sha256']='0'*64
            with self.assertRaisesRegex(ValueError,'binding differs'):
                verified_phase_evidence(manifest,'study',loci,sam_task,reference_task,
                    samtools=samtools,policy=policy,rehash_inputs=True)
