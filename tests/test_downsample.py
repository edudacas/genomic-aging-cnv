import importlib.util
from pathlib import Path
import tempfile
import unittest
from cnv_downsample import downsample_bam,template_retained
from cnv_execution import sha


class TemplateDecisionTests(unittest.TestCase):
    def test_nested_fractions_and_repeatability(self):
        low={i for i in range(1000) if template_retained('rg',str(i),17,.25)}
        high={i for i in range(1000) if template_retained('rg',str(i),17,.75)}
        self.assertTrue(low<high)
        self.assertEqual(low,{i for i in range(1000) if template_retained('rg',str(i),17,.25)})
        self.assertTrue(all(template_retained('rg',str(i),17,1) for i in range(1000)))

    def test_invalid_sampling_policy(self):
        for fraction in (0,-1,1.1,True,float('nan')):
            with self.assertRaises(ValueError):
                template_retained('rg','read',17,fraction)


@unittest.skipUnless(importlib.util.find_spec('pysam'),'pysam required for native BAM tests')
class BamDownsampleTests(unittest.TestCase):
    def setUp(self):
        import pysam
        self.pysam=pysam
        temp=tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root=Path(temp.name)
        self.bam=self.root/'source.bam'

    def make_bam(self,umi=False):
        p=self.pysam
        header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':1000}],
                'RG':[{'ID':'rg','SM':'sample'}]}
        with p.AlignmentFile(str(self.bam),'wb',header=header) as stream:
            for i in range(100):
                for flag in (99,147):
                    read=p.AlignedSegment()
                    read.query_name=f'read{i}'
                    read.query_sequence='A'*50
                    read.query_qualities=p.qualitystring_to_array('I'*50)
                    read.flag=flag
                    read.reference_id=read.next_reference_id=0
                    read.reference_start=read.next_reference_start=i*5
                    read.template_length=50 if flag==99 else -50
                    read.mapping_quality=60
                    read.cigarstring='50M'
                    read.set_tag('RG','rg')
                    if umi:
                        read.set_tag('RX','ACGT')
                    stream.write(read)

    def invoke(self,sample='sample'):
        return downsample_bam({'path':str(self.bam),'sha256':sha(self.bam)},sample,
                              self.root/'result',seed=17,fraction=.5)

    def test_native_bam_retains_mates_together_and_preserves_records(self):
        self.make_bam()
        result=self.invoke()
        with self.pysam.AlignmentFile(str(self.root/'result/subsampled.bam'),'rb') as bam:
            records=list(bam.fetch(until_eof=True))
        names={r.query_name for r in records}
        self.assertTrue(0<len(names)<100)
        self.assertEqual(len(records),2*len(names))
        for name in names:
            self.assertEqual({r.flag for r in records if r.query_name==name},{99,147})
        self.assertFalse(result['mosaic_event_simulated'])
        self.assertTrue((self.root/'result/subsampled.bam.bai').exists())

    def test_umi_requires_separate_design(self):
        self.make_bam(umi=True)
        with self.assertRaisesRegex(ValueError,'molecule-aware'):
            self.invoke()
        self.assertFalse((self.root/'result/COMPLETE.json').exists())
        self.assertTrue((self.root/'result/FAILURE.json').exists())

    def test_wrong_sample_is_rejected(self):
        self.make_bam()
        with self.assertRaisesRegex(ValueError,'identity mismatch'):
            self.invoke('other')
