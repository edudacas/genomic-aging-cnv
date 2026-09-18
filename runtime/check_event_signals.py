"""Check direction and amplitude in ideal native synthetic fixtures, not study LOD."""
import argparse
import csv
import json
import math
from pathlib import Path
from statistics import median


def check(directory):
    directory=Path(directory)
    results={}
    for event in ('LOSS','GAIN','CNLOH'):
        root=directory/event
        mixture=json.loads((root/'synthetic-mixture.json').read_text())
        with (root/'fit/native/result/signal.tsv').open() as f:
            rows=list(csv.DictReader(f,delimiter='\t'))
        affected=[r for r in rows if int(r['maploc'])<=1150000]
        neutral=[r for r in rows if int(r['maploc'])>1150000]
        if not affected or not neutral:
            raise ValueError('missing planned synthetic territory')
        depth=median(float(r['cnlr']) for r in affected)-median(float(r['cnlr']) for r in neutral)
        allelic=median(float(r['valor']) for r in affected)-median(float(r['valor']) for r in neutral)
        expected_depth=math.log2(mixture['expected_depth_ratio'])
        a=mixture['expected_haplotype_A_read_fraction']
        expected_allelic=math.log((1-a)/a)
        if abs(depth-expected_depth)>.2 or abs(allelic-expected_allelic)>.15:
            raise ValueError('native synthetic signal differs from intended '+event)
        results[event]=dict(observed_log2_depth_delta=depth,expected_log2_depth_delta=expected_depth,
            observed_ln_alt_ref_delta=allelic,expected_ln_alt_ref_delta=expected_allelic,
            affected_loci=len(affected),neutral_loci=len(neutral))
    result=dict(status='IDEAL_NATIVE_EVENT_SIGNALS_PASS',events=results,
        clinical_event_classification_validated=False,sensitivity_calibrated=False)
    with (directory/'EVENT_SIGNALS.json').open('x') as f:
        json.dump(result,f,indent=2)
        f.write('\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory')
    check(p.parse_args().directory)
