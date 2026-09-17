"""Locked, hash-bound native task execution. Never converts missing outputs to success."""
from dataclasses import asdict
from datetime import datetime, timezone
import csv
import fcntl
import gzip
import hashlib
import json
import math
import os
import signal
from pathlib import Path
import subprocess
import tempfile
import time


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    descriptor, name = tempfile.mkstemp(prefix='.'+path.name, dir=path.parent)
    try:
        with os.fdopen(descriptor, 'w') as stream:
            json.dump(value, stream, indent=2, allow_nan=False)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def output_metadata(command):
    result = {}
    for name in command.outputs:
        path = Path(name)
        if not path.is_file() or not path.stat().st_size:
            raise ValueError('missing/empty caller output: ' + name)
        if command.table_columns:
            opener = gzip.open if path.suffix == '.gz' else open
            with opener(path, 'rt', newline='') as stream:
                header = stream.readline().rstrip('\r\n')
                delimiter = '\t' if '\t' in header else ','
                columns = next(csv.reader([header], delimiter=delimiter))
                if not set(command.table_columns) <= set(columns):
                    raise ValueError('native output schema mismatch: ' + name)
                # Exhaust compressed stream to detect truncation/CRC errors.
                for row in csv.reader(stream, delimiter=delimiter):
                    if len(row) != len(columns):
                        raise ValueError('ragged native table: ' + name)
        result[str(path.resolve())] = {'sha256': sha(path), 'bytes': path.stat().st_size}
    return result


def verify_completed_task(command, directory, expected_inputs, *, rehash_inputs):
    """Read-only receipt verification against an independently expected command.

    expected_inputs maps resolved paths to hashes from the frozen manifest.
    rehash_inputs explicitly controls re-reading potentially very large BAMs.
    """
    if type(rehash_inputs) is not bool:
        raise ValueError('explicit input rehash policy required')
    directory = Path(directory).resolve()
    if any(not Path(p).is_absolute() or not Path(p).resolve().is_relative_to(directory)
           for p in command.outputs):
        raise ValueError('declared outputs must stay within the task directory')
    started = json.loads((directory/'STARTED.json').read_text())
    complete = json.loads((directory/'COMPLETE.json').read_text())
    frozen = {key:started[key] for key in ('command','inputs','implementation')}
    if 'runtime_environment' in started:
        frozen['runtime_environment'] = validate_runtime_environment(started['runtime_environment'])
    fingerprint = hashlib.sha256(json.dumps(frozen,sort_keys=True).encode()).hexdigest()
    if (started.get('fingerprint') != fingerprint or complete.get('fingerprint') != fingerprint or
            complete.get('status') != 'NATIVE_TASK_COMPLETE'):
        raise ValueError('invalid completion receipt or fingerprint')
    if frozen['command'] != json.loads(json.dumps(asdict(command))):
        raise ValueError('native command or filter policy differs from expectation')
    for path,digest in expected_inputs.items():
        if frozen['inputs'].get(str(Path(path).resolve())) != digest:
            raise ValueError('native input binding differs from manifest')
    if not frozen['implementation']:
        raise ValueError('missing implementation provenance')
    if any(sha(path)!=digest for path,digest in frozen['implementation'].items()):
        raise ValueError('native implementation changed')
    if rehash_inputs and any(sha(path)!=digest for path,digest in frozen['inputs'].items()):
        raise ValueError('native input bytes changed')
    if output_metadata(command) != complete.get('outputs'):
        raise ValueError('native output bytes changed')
    return dict(status='COMPLETED_TASK_VERIFIED',fingerprint=fingerprint,
                manifest_bindings_verified=True,input_bytes_rehashed=rehash_inputs,
                outputs_verified=True)


def validate_runtime_environment(value):
    allowed = {'R_LIBS','R_LIBS_USER','R_LIBS_SITE','R_PROFILE','R_PROFILE_USER',
               'R_ENVIRON','R_ENVIRON_USER'}
    if (not isinstance(value, dict) or not value or not set(value)<=allowed or
            any(not isinstance(v,str) or not v or '\x00' in v for v in value.values())):
        raise ValueError('invalid explicit R runtime environment')
    return dict(value)


def execute(command, directory, input_files, *, implementation_files, timeout_seconds,
            runtime_environment=None):
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    if isinstance(timeout_seconds, bool) or not isinstance(timeout_seconds, (int, float)) or not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
        raise ValueError('timeout must be positive')
    if any(not Path(p).is_absolute() for p in command.outputs):
        raise ValueError('output paths must be absolute')
    output_paths = [Path(p).resolve() for p in command.outputs]
    if len(set(output_paths)) != len(output_paths):
        raise ValueError('duplicate output paths')
    if any(not p.is_relative_to(directory) for p in output_paths):
        raise ValueError('outputs must stay inside task directory')
    reserved = {'LOCK', 'STARTED.json', 'COMPLETE.json', 'FAILURE.json', 'commands.log'}
    if any(p == directory / name for p in output_paths for name in reserved):
        raise ValueError('output collides with task bookkeeping')
    inputs = {str(Path(p).resolve()): sha(p) for p in input_files}
    implementation = {str(Path(p).resolve()): sha(p) for p in implementation_files}
    if not implementation:
        raise ValueError('caller implementation hashes required')
    if set(inputs) & {str(p) for p in output_paths}:
        raise ValueError('input/output collision')
    frozen = dict(command=asdict(command), inputs=inputs, implementation=implementation)
    environment = None
    if runtime_environment is not None:
        frozen['runtime_environment'] = validate_runtime_environment(runtime_environment)
        environment = dict(os.environ, **frozen['runtime_environment'])
    fingerprint = hashlib.sha256(json.dumps(frozen, sort_keys=True).encode()).hexdigest()
    with (directory/'LOCK').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        complete = directory/'COMPLETE.json'
        if complete.exists():
            receipt = json.loads(complete.read_text())
            if receipt.get('fingerprint') != fingerprint:
                raise ValueError('task inputs/implementation changed')
            if receipt.get('outputs') != output_metadata(command):
                raise ValueError('completed outputs changed')
            return dict(receipt, reused=True)
        # Explicitly preserve incomplete attempts; caller chooses a fresh directory.
        if (directory/'STARTED.json').exists() or any(p.exists() for p in output_paths):
            raise ValueError('incomplete task preserved; use a new attempt directory')
        started = time.monotonic()
        atomic_json(directory/'STARTED.json', dict(frozen, fingerprint=fingerprint,
                    started_utc=datetime.now(timezone.utc).isoformat()))
        try:
            with (directory/'commands.log').open('xb') as log:
                process = subprocess.Popen(command.argv, cwd=directory, stdout=log,
                                           stderr=subprocess.STDOUT, shell=False,
                                           env=environment,
                                           start_new_session=True)
                try:
                    returncode = process.wait(timeout=timeout_seconds)
                except subprocess.TimeoutExpired:
                    # Terminate the whole task group, including R/caller children.
                    try:
                        os.killpg(process.pid, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        pass
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    process.wait()
                    raise
                if returncode:
                    raise subprocess.CalledProcessError(returncode, command.argv)
            outputs = output_metadata(command)
            # Detect inputs modified while the external tool ran.
            if any(sha(path) != digest for path, digest in (inputs | implementation).items()):
                raise ValueError('input/implementation changed during execution')
            receipt = dict(status='NATIVE_TASK_COMPLETE', fingerprint=fingerprint,
                           outputs=outputs, elapsed_seconds=time.monotonic()-started)
            atomic_json(complete, receipt)
            return dict(receipt, reused=False)
        except Exception as error:
            atomic_json(directory/'FAILURE.json', dict(error=type(error).__name__,
                        message=str(error), elapsed_seconds=time.monotonic()-started))
            raise
