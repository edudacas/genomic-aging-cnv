import math
import unittest
from cnv_allelic import MarkerPolicy, summarize_marker, subtract_bias, to_natural_log, allelic_log_ratio


def policy():
    # Fixture settings only; not a proposed study configuration.
    return MarkerPolicy(10, 200, 0.1, 0.9, 2, 0.75, 10, 1, 0.5, 1/6)


class AllelicTests(unittest.TestCase):
    def test_native_pseudocount_and_training_estimator_agree(self):
        expected = math.log(20 + 1/6) - math.log(10 + 1/6)
        self.assertAlmostEqual(allelic_log_ratio(10, 20, pseudocount=1/6), expected)
        result = summarize_marker({'a': (10, 20), 'b': (10, 20)}, ['a', 'b'], policy())
        self.assertAlmostEqual(result['correction_bias'], expected)
        self.assertAlmostEqual(subtract_bias(expected, result['correction_bias'],
                                            bias_unit=result['bias_unit']), 0)
        self.assertIsNone(allelic_log_ratio(0, 0, pseudocount=1/6))

    def test_log_units_and_signed_correction(self):
        self.assertAlmostEqual(to_natural_log(1, 'log2_alt_ref'), math.log(2))
        self.assertAlmostEqual(subtract_bias(math.log(2), 1, bias_unit='log2_alt_ref'), 0)
        self.assertAlmostEqual(subtract_bias(-math.log(2), -1, bias_unit='log2_alt_ref'), 0)
        self.assertIsNone(subtract_bias(0.1, None, bias_unit='ln_alt_ref'))
        with self.assertRaises(ValueError):
            to_natural_log(1, 'unspecified')

    def test_missing_is_not_zero(self):
        result = summarize_marker({'a': (0, 0), 'b': None, 'c': (10, 10)},
                                  ['a', 'b', 'c', 'd'], policy())
        self.assertEqual(result['missing_normals'], 2)
        self.assertEqual(result['zero_depth_normals'], 1)
        self.assertEqual(result['median_observed_depth'], 10)
        self.assertEqual(result['heterozygous_normals'], 1)
        self.assertIsNone(result['raw_sd_bias'])
        self.assertIsNone(result['correction_bias'])

    def test_balanced_training_and_no_genotype_claim(self):
        result = summarize_marker({'a': (10, 10), 'b': (20, 20)}, ['a', 'b'], policy())
        self.assertEqual(result['status'], 'STABLE_TRAINING_MARKER')
        self.assertEqual(result['correction_bias'], 0)
        self.assertFalse(result['germline_genotype_assigned'])
        self.assertFalse(result['sensitivity_calibrated'])

    def test_rejects_foreign_samples_negative_or_fractional_counts(self):
        for counts in ({'study': (10, 10)}, {'a': (-1, 10)}, {'a': (2.1, 10)}, {'a': (True, 10)}):
            with self.subTest(counts=counts), self.assertRaises(ValueError):
                summarize_marker(counts, ['a', 'b'], policy())

    def test_single_het_cannot_prove_low_dispersion(self):
        result = summarize_marker({'a': (10, 10), 'b': (20, 0)}, ['a', 'b'], policy())
        self.assertIn('INSUFFICIENT_HETEROZYGOUS_NORMALS', result['reasons'])
        self.assertIsNone(result['correction_bias'])

    def test_empty_evidence_stays_missing(self):
        result = summarize_marker({}, ['a', 'b'], policy())
        self.assertIsNone(result['median_observed_depth'])
        self.assertIsNone(result['raw_median_bias'])
        self.assertIsNone(result['correction_bias'])


if __name__ == '__main__':
    unittest.main()
