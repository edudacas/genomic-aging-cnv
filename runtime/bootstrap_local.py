"""Install an authorized persistent native test environment, outside the software repo.

Uses the host R/samtools/compiler. This is NOT a portable container image.
Preserves downloaded sources, logs and resolved dependency versions.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import urllib.request


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def main(root):
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    if (root/'NATIVE_INSTALLED.json').exists():
        raise ValueError('runtime already installed; verify it instead of overwriting')
    repo = Path(__file__).resolve().parents[1]
    def run(command, **kw):
        print('RUN', ' '.join(map(str,command)), flush=True)
        with (root/'install.log').open('a') as log:
            subprocess.run(list(map(str,command)), check=True, stdout=log, stderr=log, **kw)
    if not (root/'sources').exists():
        run([sys.executable, Path(__file__).with_name('fetch_sources.py'), root/'sources'])
    sources = root/'sources'
    source = next((sources/'facets').iterdir())/'facets-ch'
    gcsource = next((sources/'pctGCdata').iterdir())
    python = root/'python/bin/python'
    if not python.exists():
        run([sys.executable, '-m','venv',root/'python'])
    run([python,'-m','pip','install','CNVkit==0.9.12','pysam==0.23.3',
         'pyarrow==25.0.1','numpy==1.26.4'])
    with (root/'python.freeze.txt').open('w') as f:
        subprocess.run([str(python),'-m','pip','freeze','--all'],stdout=f,check=True)
    library = root/'R-library'
    if not library.exists():
        library.mkdir()
        # Copy existing package bytes, never update the system library.
        system_libraries = subprocess.check_output(['Rscript','--vanilla','-e',
            'cat(.libPaths(),sep="\\n")'],text=True).splitlines()
        for system in system_libraries:
            for package in Path(system).iterdir():
                if package.is_dir() and (package/'DESCRIPTION').is_file() and not (library/package.name).exists():
                    shutil.copytree(package,library/package.name,symlinks=False)
    archives=root/'R-sources'
    archives.mkdir(exist_ok=True)
    env=dict(os.environ, R_LIBS_USER=str(library),R_LIBS_SITE=str(library),
             R_PROFILE_USER='/dev/null',R_ENVIRON_USER='/dev/null')
    run(['Rscript','--vanilla',Path(__file__).with_name('install_r.R'),library,archives],env=env)
    dna=archives/'DNAcopy_1.86.0.tar.gz'
    if not dna.exists():
        urllib.request.urlretrieve('https://bioconductor.org/packages/3.23/bioc/src/contrib/DNAcopy_1.86.0.tar.gz',dna)
    for package in [dna,gcsource,source]:
        run(['R','CMD','INSTALL','--library='+str(library),package],env=env)
    bindir=root/'bin'
    bindir.mkdir(exist_ok=True)
    run(['g++','-O2','-o',bindir/'snp-pileup',source/'inst/extcode/snp-pileup.cpp','-lhts'])
    run(['Rscript','--vanilla','-e',
        'stopifnot(all(vapply(c("DNAcopy","facetsCH","pctGCdata"),requireNamespace,logical(1),quietly=TRUE))); stopifnot(as.character(packageVersion("DNAcopy"))=="1.86.0")'],env=env)
    commands={name:shutil.which(name) for name in ('Rscript','R','samtools','g++')}
    manifest=dict(status='PERSISTENT_NATIVE_INSTALLED_NOT_CONTAINER',
        host_executables={k:dict(path=v,sha256=sha(v)) for k,v in commands.items()},
        source_archives={str(p.relative_to(root)):sha(p) for p in root.rglob('*.tar.gz')},
        python_lock_sha256=sha(root/'python.freeze.txt'),
        snp_pileup_sha256=sha(bindir/'snp-pileup'),
        native_tests_passed=False, container_built=False, cohort_executed=False,
        dependency_closure_portable=False)
    (root/'NATIVE_INSTALLED.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root')
    main(p.parse_args().root)
