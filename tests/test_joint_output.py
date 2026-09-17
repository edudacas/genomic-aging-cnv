import unittest
from cnv_joint_output import validate_joint_output
import test_run_joint


class JointOutputTests(unittest.TestCase):
    def setUp(self):
        self.fixture=test_run_joint.RunJointTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.fixture.test_complete_contract_and_declared_scope()
        self.output=self.fixture.root/'run/native/result'

    def test_incomplete_audit_rejected(self):
        with self.assertRaisesRegex(ValueError,'does not cover'):
            validate_joint_output(self.output,2)

    def test_false_biological_claim_rejected(self):
        path=self.output/'segments.tsv'
        path.write_text(path.read_text().replace('UNCLASSIFIED','CNLOH'))
        with self.assertRaisesRegex(ValueError,'biological claim'):
            validate_joint_output(self.output,1)

    def test_broken_source_linkage_rejected(self):
        path=self.output/'signal.tsv'
        path.write_text(path.read_text().replace('\t1\t1\n','\t1\t2\n'))
        with self.assertRaisesRegex(ValueError,'input linkage'):
            validate_joint_output(self.output,1)

    def test_success_is_structural_only(self):
        result=validate_joint_output(self.output,1)
        self.assertEqual(result['status'],'JOINT_OUTPUT_STRUCTURE_PASS')
        self.assertFalse(result['scientific_validation'])

    def test_audit_must_match_original_count_matrix(self):
        matrix=self.fixture.fixture.f.f.output/'pair.rcmat.tsv'
        result=validate_joint_output(self.output,1,input_matrix=matrix)
        self.assertTrue(result['submitted_matrix_rows_verified'])
        path=self.output/'locus_audit.tsv'
        original=path.read_text()
        for changed in (original.replace('\tchr1\t1\t','\tchr1\t2\t'),
                        original.replace('\t20\t0\n','\t20\t1\n')):
            path.write_text(changed)
            with self.assertRaisesRegex(ValueError,'locus/counts differ'):
                validate_joint_output(self.output,1,input_matrix=matrix)
