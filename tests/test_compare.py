import json
from pathlib import Path
import sys
import tempfile
import unittest
from cnv_commands import Command
from cnv_compare import verified_source
from cnv_execution import execute,sha


class CompareSourceTests(unittest.TestCase):
    def setUp(self):
        temp=tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root=Path(temp.name)
        self.input=self.root/'input'
        self.input.write_text('synthetic input')
        self.output=self.root/'native/sample.cns'
        command=Command((sys.executable,'-c',
            'from pathlib import Path; Path("sample.cns").write_text("chromosome\\tstart\\tend\\tlog2\\tprobes\\nchr1\\t0\\t10\\t0\\t1\\n")'),
            (str(self.output),))
        native=execute(command,self.root/'native',[self.input],implementation_files=[__file__],timeout_seconds=10)
        (self.root/'STARTED.json').write_text('{}')
        self.provenance=dict(sample_id='study',stratum='S',manifest_sha256='a'*64,
                             sample_bam_sha256='b'*64,reference_build='GRCh38')
        self.receipt=dict(status='CNVKIT_CBS_COMPLETE',sample_provenance=self.provenance,
            started_sha256=sha(self.root/'STARTED.json'),outputs=native['outputs'])
        self.path=self.root/'COMPLETE.json'
        self.seal()

    def seal(self):
        self.path.write_text(json.dumps(self.receipt))
        self.artifact=dict(path=str(self.path),sha256=sha(self.path))

    def verify(self):
        return verified_source(self.artifact,'CNVKIT_CBS_COMPLETE',self.provenance,'native/sample.cns')

    def test_matching_source_and_native_receipt(self):
        self.assertEqual(self.verify(),self.output)

    def test_same_sample_label_different_bam_rejected(self):
        self.receipt['sample_provenance']=dict(self.provenance,sample_bam_sha256='c'*64)
        self.seal()
        with self.assertRaisesRegex(ValueError,'provenance mismatch'):
            self.verify()

    def test_changed_segment_rejected(self):
        self.output.write_text('modified')
        with self.assertRaisesRegex(ValueError,'output bytes changed'):
            self.verify()

    def test_changed_native_input_rejected(self):
        self.input.write_text('changed')
        with self.assertRaisesRegex(ValueError,'input bytes changed'):
            self.verify()

    def test_failed_stage_rejected_even_with_completion(self):
        (self.root/'FAILURE.json').write_text('{}')
        with self.assertRaisesRegex(ValueError,'invalid completed'):
            self.verify()
