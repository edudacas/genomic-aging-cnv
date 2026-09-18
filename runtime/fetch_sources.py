"""Fetch audited archives and verify checksums before safe extraction."""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile
import urllib.request


def fetch(destination):
    root = Path(destination)
    root.mkdir(parents=True, exist_ok=False)
    lock = json.loads(Path(__file__).with_name('source-lock.json').read_text())
    for name, item in lock.items():
        archive = root/(name + '.tar.gz')
        urllib.request.urlretrieve(item['url'], archive)
        if hashlib.sha256(archive.read_bytes()).hexdigest() != item['sha256']:
            raise ValueError('archive checksum mismatch: ' + name)
        with tarfile.open(archive) as tf:
            tf.extractall(root/name, filter='data')
    source = next((root/'facets').iterdir())/'facets-ch'
    description = source/'DESCRIPTION'
    original = description.read_text()
    if original.count('Package: FACETS-CH') != 1:
        raise ValueError('unexpected upstream package declaration')
    description.write_text(original.replace('Package: FACETS-CH', 'Package: facetsCH', 1))
    (root/'packaging-repair.json').write_text(json.dumps(dict(
        change='Package: FACETS-CH -> facetsCH', reason='match upstream NAMESPACE and CLI',
        before_sha256=hashlib.sha256(original.encode()).hexdigest(),
        after_sha256=hashlib.sha256(description.read_bytes()).hexdigest(),
        algorithm_changed=False), indent=2)+'\n')
    print(source)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('destination')
    fetch(p.parse_args().destination)
