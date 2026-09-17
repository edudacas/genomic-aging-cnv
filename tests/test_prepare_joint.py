import csv
import json
import unittest
from cnv_prepare_joint import prepare_joint,validate_fit_policy
import test_joint_inputs


class PrepareJointTests(unittest.TestCase):
    def setUp(self):
        self.f=test_joint_inputs.JointInputTests()
        self.f.setUp()
        self.addCleanup(self.f.doCleanups)
        self.policy=dict(ndepth=20,ndepthmax=2000,het_threshold=.1,snp_neighborhood=0,
            initial_cval=25,final_cval=150,gamma=.25,min_nhet=15,seed=0,
            sample_het_min_depth=50,delta_cn=0,het_scale=True)
        self.path=self.f.f.root/'fit.json'
        self.path.write_text(json.dumps(self.policy))
        self.output=self.f.f.root/'fit-prepared'

    def prepare(self):
        return prepare_joint(self.f.f.manifest,self.f.f.output,self.f.f.config,
                             self.f.panel,self.f.panel_config,self.path,self.output)

    def test_prepared_without_launch_and_roundtrip(self):
        receipt=self.prepare()
        self.assertFalse(receipt['facets_ch_executed'])
        self.assertFalse(receipt['runtime_validated'])
        with (self.output/'fit_config.tsv').open() as stream:
            rows={row['key']:row['value'] for row in csv.DictReader(stream,delimiter='\t')}
        self.assertEqual(rows['het_scale'],'TRUE')
        self.assertEqual(float(rows['allelic_pseudocount']),1/6)
        self.assertIn('chr1\t1',(self.output/'chromosome_map.tsv').read_text())
        with self.assertRaises(FileExistsError):
            self.prepare()

    def test_cannot_override_panel_estimator(self):
        self.policy['allelic_pseudocount']=0
        with self.assertRaisesRegex(ValueError,'estimator/map'):
            validate_fit_policy(self.policy)

    def test_invalid_parameters_rejected(self):
        for key,value in [('seed',.5),('sample_het_min_depth',20),('het_scale','TRUE'),('gamma',float('nan'))]:
            with self.subTest(key=key),self.assertRaises(ValueError):
                validate_fit_policy(dict(self.policy,**{key:value}))
