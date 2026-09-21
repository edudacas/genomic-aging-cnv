"""Build/test a local candidate image using an explicit Ubuntu 24.04 base digest.

Dependencies are resolved at build time. Image content is frozen by image ID;
this recipe alone is not a fully locked, byte-reproducible source rebuild.
"""
import argparse
import json
from pathlib import Path
import re
import subprocess


def main(base, output):
    if not re.fullmatch(r'ubuntu@sha256:[0-9a-f]{64}',base):
        raise ValueError('provide the resolved Ubuntu 24.04 digest')
    repo=Path(__file__).resolve().parents[1]
    out=Path(output).resolve()
    out.mkdir(parents=True,exist_ok=False)
    subprocess.run(['docker','build','--build-arg','BASE_IMAGE='+base,
        '-f',str(repo/'runtime/Dockerfile'),'--iidfile',str(out/'image.id'),str(repo)],check=True)
    iid=(out/'image.id').read_text().strip()
    inspected=json.loads(subprocess.check_output(['docker','image','inspect',iid],text=True))
    (out/'image.inspect.json').write_text(json.dumps(inspected,indent=2)+'\n')
    cid=subprocess.check_output(['docker','create',iid],text=True).strip()
    try:
        for source in ('/opt/cnv-validation/COMPLETE.json','/opt/cnv-validation/EVENT_SIGNALS.json',
                       '/opt/cnv-runtime/NATIVE_INSTALLED.json','/opt/cnv-runtime/FILES.json',
                       '/opt/cnv-runtime/python.freeze.txt','/opt/cnv-runtime/os-packages.tsv'):
            subprocess.run(['docker','cp',cid+':'+source,str(out/Path(source).name)],check=True)
    finally:
        # Remove only the never-started helper created above; no volumes used.
        subprocess.run(['docker','rm',cid],check=True)
    (out/'IMAGE_BUILD.json').write_text(json.dumps(dict(status='IMAGE_BUILT_NATIVE_TESTS_PASS',
        image_id=iid,base_image=base,cohort_executed=False,
        bit_reproducible_source_rebuild=False),indent=2)+'\n')
    print(iid)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base',required=True)
    p.add_argument('--output',required=True)
    a=p.parse_args()
    main(a.base,a.output)
