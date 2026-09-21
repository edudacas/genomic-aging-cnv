"""Read-only collection of completed SNV QC; no BAM scan or stage execution."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import re


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):
            h.update(b)
    return h.hexdigest()


def collect(stages, destination, expected):
    stages=Path(stages).resolve()
    destination=Path(destination)
    marker=stages/'QC_COMPLETE.json'
    initial=sha(marker)
    receipt=json.loads(marker.read_text())
    samples=receipt.get('samples',{})
    if receipt.get('status')!='QC_COMPLETE' or len(samples)!=expected:
        raise ValueError('QC completion/membership mismatch')
    result=[]
    for sid, binding in sorted(samples.items()):
        if not re.fullmatch(r'[A-Za-z0-9_-]+',sid):
            raise ValueError('unsafe sample identifier')
        directory=stages/'qc'/sid
        complete=directory/'COMPLETE.json'
        if directory.is_symlink() or sha(complete)!=binding['receipt_sha256']:
            raise ValueError('QC sample receipt mismatch')
        task=json.loads(complete.read_text())
        if task.get('status')!='COMPLETE' or task.get('fingerprint')!=receipt.get('fingerprint'):
            raise ValueError('QC task status/fingerprint mismatch')
        selected={}
        for name,digest in task['outputs'].items():
            relative=Path(name)
            if relative.is_absolute() or '..' in relative.parts:
                raise ValueError('unsafe QC artifact')
            path=directory/relative
            if path.is_symlink() or directory.resolve() not in path.resolve().parents:
                raise ValueError('unsafe QC artifact path')
            if path.name not in ('qc.json','coverage.summary.json','coverage.compartments.tsv'):
                continue
            if sha(path)!=digest:
                raise ValueError('QC artifact changed')
            selected[name]=dict(sha256=digest,content=path.read_text())
        found={Path(n).name for n in selected}
        if found!={'qc.json','coverage.summary.json','coverage.compartments.tsv'}:
            raise ValueError('missing compact QC outputs')
        result.append(dict(sample_id=sid,receipt_sha256=binding['receipt_sha256'],artifacts=selected))
    if sha(marker)!=initial:
        raise ValueError('QC receipt changed during collection')
    payload=dict(status='QC_SELECTED_ARTIFACTS_VERIFIED',samples=result,
        qc_complete_sha256=initial,selected_samples=len(result),
        interval_files_included=False,all_task_outputs_rehashed=False,
        sex_inferred=False,umi_assessed=False,bams_scanned=False)
    with destination.open('xb') as raw:
        os.chmod(destination,0o600)
        with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as f:
            f.write(json.dumps(payload,sort_keys=True).encode())
    print(json.dumps(dict(status=payload['status'],samples=len(result),
        output=str(destination),sha256=sha(destination))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stages',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--expected',type=int,default=1329)
    a=p.parse_args()
    collect(a.stages,a.output,a.expected)
