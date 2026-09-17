import unittest
from cnv_simulation import diploid_mixture


class MixtureTests(unittest.TestCase):
    def model(self,event,fraction,hap='A'):
        return diploid_mixture(event,fraction,affected_haplotype=hap)

    def test_loss_changes_depth_and_balance_jointly(self):
        m=self.model('LOSS',.4)
        self.assertAlmostEqual(m['expected_depth_ratio'],.8)
        self.assertAlmostEqual(m['expected_haplotype_A_read_fraction'],.375)
        self.assertAlmostEqual(m['expected_clone_read_fraction'],.25)

    def test_gain_read_fraction_is_not_cell_fraction(self):
        m=self.model('GAIN',.4)
        self.assertAlmostEqual(m['expected_depth_ratio'],1.2)
        self.assertAlmostEqual(m['expected_haplotype_A_read_fraction'],1.4/2.4)
        self.assertAlmostEqual(m['expected_clone_read_fraction'],.5)

    def test_cnloh_keeps_depth_but_changes_allele_fraction(self):
        m=self.model('CNLOH',.4)
        self.assertEqual(m['expected_depth_ratio'],1)
        self.assertAlmostEqual(m['expected_haplotype_A_read_fraction'],.3)
        self.assertAlmostEqual(m['expected_haplotype_B_read_fraction'],.7)

    def test_orientation_swaps_alleles_not_dosage(self):
        for event in ('LOSS','GAIN','CNLOH'):
            a,b=self.model(event,.2,'A'),self.model(event,.2,'B')
            self.assertEqual(a['expected_depth_ratio'],b['expected_depth_ratio'])
            self.assertEqual(a['expected_haplotype_A_read_fraction'],b['expected_haplotype_B_read_fraction'])

    def test_zero_and_full_clone_limits(self):
        for event in ('LOSS','GAIN','CNLOH'):
            m=self.model(event,0)
            self.assertEqual(m['expected_depth_ratio'],1)
            self.assertEqual(m['expected_haplotype_A_read_fraction'],.5)
        self.assertEqual(self.model('CNLOH',1)['expected_haplotype_A_read_fraction'],0)
        self.assertEqual(self.model('GAIN',1)['expected_depth_ratio'],1.5)

    def test_invalid_inputs(self):
        for f in (True,float('nan'),-1,1.1):
            with self.assertRaises(ValueError):
                self.model('LOSS',f)
