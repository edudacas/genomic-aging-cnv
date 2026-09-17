import unittest
import csv
import io
from pathlib import Path
import shutil
import subprocess
from cnv_territory import Territory, annotate_segment_territory


class TerritoryTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which('Rscript'), 'Rscript unavailable')
    def test_r_segment_export_to_python_territory(self):
        code = """
        source('facets_segments.R')
        x <- data.frame(chrom=c(1,1),maploc=c(10,90),seg=c(1,1),
                        cnlr=c(0,0.1),valor=c(NA,NA),het=c(0,0))
        write.table(summarize_native_segments(x,c(chr1='1')),stdout(),
                    sep='\t',row.names=FALSE,quote=FALSE,na='NA')
        """
        run = subprocess.run(['Rscript','-e',code],
                             cwd=Path(__file__).resolve().parents[1],
                             capture_output=True,text=True,check=True,timeout=30)
        row = next(csv.DictReader(io.StringIO(run.stdout), delimiter='\t'))
        for key in ('start_0based','end_0based_exclusive'):
            row[key] = int(row[key])
        assay = Territory([('chr1',0,20),('chr1',80,100)], {'chr1':100})
        result = next(annotate_segment_territory([row],assay))
        self.assertEqual(result['marker_bounding_span_bp'],81)
        self.assertEqual(result['assay_overlap_bp'],21)
        self.assertEqual(result['evidence_class'],'CNV_ONLY_SIGNAL')
        self.assertIsNone(result['evaluable_bp'])

    def test_gaps_are_not_examined_bases(self):
        territory = Territory([('chr1',0,10),('chr1',90,100)], {'chr1':100})
        segment = dict(chromosome='chr1',start_0based=0,end_0based_exclusive=100)
        result = next(annotate_segment_territory([segment], territory))
        self.assertEqual(result['marker_bounding_span_bp'],100)
        self.assertEqual(result['assay_overlap_bp'],20)
        self.assertIsNone(result['coverage_eligible_overlap_bp'])
        self.assertIsNone(result['calibrated_overlap_bp'])
        self.assertIsNone(result['evaluable_bp'])

    def test_empty_mask_zero_is_different_from_unknown(self):
        assay = Territory([('chr1',0,100)], {'chr1':100})
        empty = Territory([], {'chr1':100})
        segment = dict(chromosome='chr1',start_0based=0,end_0based_exclusive=100)
        result = next(annotate_segment_territory([segment], assay, coverage_eligible=empty))
        self.assertEqual(result['coverage_eligible_overlap_bp'],0)
        self.assertIsNone(result['calibrated_overlap_bp'])

    def test_overlap_matches_basewise_oracle(self):
        rows = [('chr1',0,3),('chr1',3,7),('chr1',11,15),('chr1',18,20)]
        territory = Territory(rows, {'chr1':20})
        bases = {i for _,start,end in rows for i in range(start,end)}
        for start in range(20):
            for end in range(start+1,21):
                self.assertEqual(territory.overlap('chr1',start,end),
                                 len(bases.intersection(range(start,end))))

    def test_calibrated_geometry_does_not_prove_calibration(self):
        assay = Territory([('chr1',0,100)], {'chr1':100})
        eligible = Territory([('chr1',0,20)], {'chr1':100})
        calibrated = Territory([('chr1',5,10)], {'chr1':100})
        segment = dict(chromosome='chr1',start_0based=0,end_0based_exclusive=100)
        result = next(annotate_segment_territory([segment], assay,
                      coverage_eligible=eligible, calibrated=calibrated))
        self.assertEqual(result['calibrated_overlap_bp'],5)
        self.assertFalse(result['territory_calibration_provenance_verified'])
        self.assertIsNone(result['evaluable_bp'])
        outside = Territory([('chr1',90,100)], {'chr1':100})
        with self.assertRaisesRegex(ValueError,'outside parent'):
            list(annotate_segment_territory([segment],assay,
                 coverage_eligible=eligible,calibrated=outside))

    def test_invalid_or_overlapping_territory_rejected(self):
        for rows in ([('chr1',0,0)],[('chr1',0,11)],
                     [('chr1',0,5),('chr1',4,10)],[('1',0,5)]):
            with self.subTest(rows=rows),self.assertRaises(ValueError):
                Territory(rows, {'chr1':10})
