"""Fingerprint an explicit isolated R library; never installs or loads packages.

This proves file stability, not successful native loading, dependency closure,
upstream authenticity, shared-library compatibility or scientific validation.
"""
from pathlib import Path
from cnv_execution import sha


def snapshot_library(directory, *, required_packages):
    root=Path(directory).resolve()
    if not root.is_dir():
        raise ValueError('isolated R library missing')
    required=list(required_packages)
    if (not required or len(required)!=len(set(required)) or
            any(not isinstance(p,str) or not p or '/' in p or '\\' in p or p in ('.','..') for p in required)):
        raise ValueError('explicit unique package names required')
    for package in required:
        if not (root/package/'DESCRIPTION').is_file() or not (root/package/'Meta/package.rds').is_file():
            raise ValueError('installed package metadata missing: '+package)
    files={}
    for path in sorted(root.rglob('*')):
        relative=path.relative_to(root).as_posix()
        if any(part.startswith('00LOCK') for part in path.relative_to(root).parts):
            raise ValueError('R installation lock present; library may be incomplete')
        if path.is_symlink():
            target=path.resolve(strict=True)
            if not target.is_relative_to(root) or target.is_dir():
                raise ValueError('external or directory symlink in isolated library')
        if path.is_file():
            files[relative]=dict(sha256=sha(path),bytes=path.stat().st_size,
                link_target=path.resolve().relative_to(root).as_posix() if path.is_symlink() else None)
        elif not path.is_dir():
            raise ValueError('unsupported library filesystem object')
    if not files:
        raise ValueError('empty isolated R library')
    return dict(status='R_LIBRARY_FILE_SNAPSHOT',directory=str(root),
        required_packages=sorted(required),files=files,
        native_load_verified=False,dependency_closure_verified=False)


def verify_library(snapshot):
    current=snapshot_library(snapshot['directory'],required_packages=snapshot['required_packages'])
    if current!=snapshot:
        raise ValueError('isolated R library differs from frozen snapshot')
    return dict(status='R_LIBRARY_FILES_VERIFIED',file_count=len(current['files']),
                native_load_verified=False,dependency_closure_verified=False)
