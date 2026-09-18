"""Fingerprint installed runtime and resolved ELF shared libraries (host-bound)."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):
            h.update(b)
    return h.hexdigest()


def freeze(root, output):
    root=Path(root).resolve()
    files={}
    elfs=set()
    for subtree in ('python','R-library','bin'):
        for p in sorted((root/subtree).rglob('*')):
            if not p.is_file() or '__pycache__' in p.parts or p.suffix=='.pyc':
                continue
            files[str(p.relative_to(root))]=dict(sha256=sha(p),bytes=p.stat().st_size,
                external_target=str(p.resolve()) if p.is_symlink() else None)
            with p.open('rb') as stream:
                if stream.read(4)==b'\x7fELF':
                    elfs.add(p.resolve())
    for p in (Path('/usr/bin/samtools'),Path('/usr/bin/Rscript'),Path('/usr/lib/R/bin/exec/R')):
        if p.exists():
            elfs.add(p.resolve())
    shared={}
    unresolved={}
    for p in sorted(elfs):
        proc=subprocess.run(['ldd',str(p)],capture_output=True,text=True)
        report=proc.stdout+proc.stderr
        if 'not found' in report:
            # A wheel's secondary ELF may rely on its parent's RUNPATH.
            # Retain this uncertainty rather than claiming closure from ldd.
            unresolved[str(p)]=[line.strip() for line in report.splitlines() if 'not found' in line]
        for match in re.findall(r'(/[^\s()]+)',report):
            lib=Path(match)
            if lib.is_file():
                shared[str(lib.resolve())]=sha(lib)
    result=dict(status='HOST_BOUND_RUNTIME_FINGERPRINT',root=str(root),files=files,
        external_elfs={str(p):sha(p) for p in sorted(elfs) if not p.is_relative_to(root)},
        resolved_shared_libraries=shared,unresolved_elf_probes=unresolved,container_built=False,
        portable_dependency_closure=False,cohort_executed=False)
    with Path(output).open('x') as stream:
        json.dump(result,stream,sort_keys=True,indent=2)
        stream.write('\n')
    print(json.dumps(dict(status=result['status'],runtime_files=len(files),
        shared_libraries=len(shared),manifest_sha256=sha(output))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',required=True)
    p.add_argument('--output',required=True)
    a=p.parse_args()
    freeze(a.root,a.output)
