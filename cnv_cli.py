"""Explicit-stage local CNV development interface; no automatic cohort chaining."""
import argparse
import json
from pathlib import Path


def parser():
    root=argparse.ArgumentParser(description=__doc__)
    commands=root.add_subparsers(dest='stage',required=True)
    validate=commands.add_parser('validate-manifest',help='Read-only input contract validation')
    validate.add_argument('--manifest',required=True)
    validate.add_argument('--validation-stage',choices=('reference_build','analysis'),required=True)
    validate.add_argument('--check-files',action='store_true')
    for name in ('prepare-pair','train-panel'):
        command=commands.add_parser(name,help='Build evidence artifacts; never launch a caller')
        command.add_argument('--manifest',required=True)
        command.add_argument('--config',required=True)
        command.add_argument('--output',required=True)
    prepare=commands.add_parser('prepare-joint',help='Verify artifacts and write config; do not execute R')
    for name in ('manifest','pair-directory','pair-config','panel-directory','panel-config','fit-policy','output'):
        prepare.add_argument('--'+name,required=True)
    runtime=commands.add_parser('snapshot-r-library',help='Hash an existing library; never install packages')
    runtime.add_argument('--library',required=True)
    runtime.add_argument('--require-package',action='append',required=True)
    runtime.add_argument('--output',required=True)
    run=commands.add_parser('run-joint',help='Execute ONE prepared R fit; not a cohort runner')
    for name in ('prepared','library-snapshot','rscript-artifact','output'):
        run.add_argument('--'+name,required=True)
    run.add_argument('--timeout-seconds',required=True,type=float)
    run.add_argument('--execute',action='store_true',help='Required explicit acknowledgement of local native execution')
    for name in ('run-reference','run-normalize','run-cbs'):
        command=commands.add_parser(name,help='Execute one explicit CNVkit stage; no automatic chaining')
        command.add_argument('--config',required=True)
        command.add_argument('--output',required=True)
        command.add_argument('--implementation-file',action='append',required=True)
        command.add_argument('--timeout-seconds',type=float,required=True)
        command.add_argument('--execute',action='store_true')
        if name=='run-cbs':
            command.add_argument('--library-snapshot',required=True)
            command.add_argument('--rscript-artifact',required=True)
            command.add_argument('--processes',type=int,required=True)
            command.add_argument('--threshold',type=float,required=True)
    compare=commands.add_parser('compare-callers',help='Verify existing stage receipts and reconcile; no native launch')
    for name in ('manifest-artifact','sample-id','cbs-artifact','joint-artifact','output'):
        compare.add_argument('--'+name,required=True)
    compare.add_argument('--logr-direction-threshold',type=float,required=True)
    export=commands.add_parser('export-comparison',help='Export existing comparison evidence, not biological calls')
    for name in ('receipt-artifact','sample-id','output'):
        export.add_argument('--'+name,required=True)
    export.add_argument('--batch-size',type=int,required=True)
    return root


def main(argv=None):
    cli=parser()
    args=cli.parse_args(argv)
    if args.stage.startswith('run-') and not args.execute:
        cli.error(args.stage+' requires --execute; preparation never launches callers')
    if args.stage=='validate-manifest':
        from cnv_contract import validate_manifest
        result=validate_manifest(json.loads(Path(args.manifest).read_text()),args.check_files,
                                 stage=args.validation_stage)
    elif args.stage=='prepare-pair':
        from cnv_prepare_pair import prepare
        result=prepare(args.manifest,args.config,args.output)
    elif args.stage=='train-panel':
        from cnv_train_panel import train
        result=train(args.manifest,args.config,args.output)
    elif args.stage=='prepare-joint':
        from cnv_prepare_joint import prepare_joint
        result=prepare_joint(args.manifest,args.pair_directory,args.pair_config,args.panel_directory,
                             args.panel_config,args.fit_policy,args.output)
    elif args.stage=='snapshot-r-library':
        from cnv_r_library import snapshot_library
        result=snapshot_library(args.library,required_packages=args.require_package)
        # Exclusive creation: a frozen runtime snapshot cannot be silently replaced.
        with Path(args.output).open('x') as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write('\n')
        result=dict(status=result['status'],output=str(Path(args.output).resolve()),
                    native_load_verified=False,dependency_closure_verified=False)
    elif args.stage=='run-joint':
        if not args.execute:
            cli.error('run-joint requires --execute; preparation never launches R')
        from cnv_run_joint import run_joint
        result=run_joint(args.prepared,args.library_snapshot,
                         json.loads(Path(args.rscript_artifact).read_text()),args.output,
                         timeout_seconds=args.timeout_seconds)
    elif args.stage in ('run-reference','run-normalize'):
        config=json.loads(Path(args.config).read_text())
        required={'manifest_path','executable','processes','min_mapq','minimum_pon_samples'}
        required |= ({'registry_path','stratum_id'} if args.stage=='run-reference' else
                     {'sample_id','sample_registry_path','pon_registry_path','reference_directory'})
        if not isinstance(config,dict) or set(config)!=required:
            cli.error('exact configuration keys required for '+args.stage)
        if args.stage=='run-reference':
            from cnv_reference import build_reference as native_stage
        else:
            from cnv_normalize import normalize_sample as native_stage
        result=native_stage(**config,output_directory=args.output,
                            implementation_files=args.implementation_file,
                            timeout_seconds=args.timeout_seconds)
    elif args.stage=='run-cbs':
        from cnv_cbs import run_cbs
        result=run_cbs(args.config,args.library_snapshot,
            json.loads(Path(args.rscript_artifact).read_text()),args.output,
            processes=args.processes,threshold=args.threshold,
            implementation_files=args.implementation_file,timeout_seconds=args.timeout_seconds)
    elif args.stage=='compare-callers':
        from cnv_compare import compare_callers
        result=compare_callers(json.loads(Path(args.manifest_artifact).read_text()),args.sample_id,
            json.loads(Path(args.cbs_artifact).read_text()),json.loads(Path(args.joint_artifact).read_text()),
            args.output,logr_direction_threshold=args.logr_direction_threshold)
    elif args.stage=='export-comparison':
        from cnv_compare_export import export_comparison
        result=export_comparison(json.loads(Path(args.receipt_artifact).read_text()),args.sample_id,
                                 args.output,batch_size=args.batch_size)
    print(json.dumps(result,sort_keys=True,allow_nan=False))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
