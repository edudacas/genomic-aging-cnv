"""Run portable regression and LOSS/GAIN/CNLOH native integration on synthetic BAMs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import unittest
import sys


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(root, output):
    root=Path(root).resolve()
    output=Path(output).resolve()
    output.mkdir(parents=True,exist_ok=False)
    repo=Path(__file__).resolve().parents[1]
    env=dict(os.environ,R_LIBS_USER=str(root/'R-library'),R_LIBS_SITE=str(root/'R-library'),
             PATH=str(root/'python/bin')+os.pathsep+os.environ['PATH'])
    commands=[('portable',[root/'python/bin/python','-m','unittest','discover','-s','tests','-v'])]
    for event in ('LOSS','GAIN','CNLOH'):
        commands.append((event,[root/'python/bin/python',repo/'tests/native_joint_chain_smoke.py',
            '--pileup',root/'bin/snp-pileup','--library',root/'R-library',
            '--cnvkit',root/'python/bin/cnvkit.py','--event-class',event,
            '--output',output/event]))
    results={}
    for name, command in commands:
        print('RUN',name,flush=True)
        log=output/(name+'.log')
        with log.open('x') as stream:
            completed=subprocess.run(list(map(str,command)),cwd=repo,env=env,stdout=stream,stderr=stream)
        results[name]=dict(exit_code=completed.returncode,log_sha256=sha(log))
        if completed.returncode:
            (output/'FAILURE.json').write_text(json.dumps(results,indent=2)+'\n')
            raise SystemExit('Native validation failed: '+name+'; see '+str(log))
        if name!='portable':
            receipt=output/name/'SMOKE_COMPLETE.json'
            results[name]['receipt_sha256']=sha(receipt)
    result=dict(status='NATIVE_SYNTHETIC_SUITE_PASS',tests=results,
        container_validated=False,study_sensitivity_calibrated=False,cohort_executed=False)
    (output/'COMPLETE.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',required=True)
    p.add_argument('--output',required=True)
    a=p.parse_args()
    main(a.root,a.output)
