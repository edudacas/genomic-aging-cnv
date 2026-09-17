from dataclasses import replace
import json
from pathlib import Path
import sys
import tempfile
import unittest
from cnv_commands import Command
from cnv_execution import execute,sha,verify_completed_task


class ReceiptVerificationTests(unittest.TestCase):
    def setUp(self):
        temp=tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root=Path(temp.name)
        self.source=self.root/'input'
        self.source.write_text('input')
        self.output=self.root/'task'/'output.tsv'
        self.command=Command((sys.executable,'-c',
            'from pathlib import Path; Path("output.tsv").write_text("x\\n1\\n")'),
            (str(self.output),),('x',))
        self.inputs={str(self.source):sha(self.source)}
        execute(self.command,self.output.parent,[self.source],
                implementation_files=[__file__],timeout_seconds=10)

    def verify(self,rehash=True):
        return verify_completed_task(self.command,self.output.parent,self.inputs,
                                     rehash_inputs=rehash)

    def test_verified_completion_and_explicit_rehash_scope(self):
        self.assertTrue(self.verify()['input_bytes_rehashed'])
        self.source.write_text('changed')
        self.assertFalse(self.verify(False)['input_bytes_rehashed'])
        with self.assertRaisesRegex(ValueError,'input bytes'):
            self.verify(True)

    def test_manifest_binding_change_rejected_even_without_rehash(self):
        self.inputs[str(self.source)]='f'*64
        with self.assertRaisesRegex(ValueError,'binding'):
            self.verify(False)

    def test_changed_command_or_output_rejected(self):
        wrong=replace(self.command,argv=(*self.command.argv,'unexpected'))
        with self.assertRaisesRegex(ValueError,'command'):
            verify_completed_task(wrong,self.output.parent,self.inputs,rehash_inputs=False)
        self.output.write_text('x\n2\n')
        with self.assertRaisesRegex(ValueError,'output bytes'):
            self.verify()

    def test_forged_success_status_rejected(self):
        path=self.output.parent/'COMPLETE.json'
        receipt=json.loads(path.read_text());receipt['status']='FAILED'
        path.write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError,'completion'):
            self.verify()
