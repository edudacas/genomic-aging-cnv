"""Explicit caller commands. Parameters are supplied by a frozen pilot config."""
from dataclasses import dataclass
from pathlib import Path
import math
import re


@dataclass(frozen=True)
class Command:
    argv: tuple
    outputs: tuple
    table_columns: tuple = ()


def positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or value <= 0:
        raise ValueError('invalid positive parameter: ' + name)
    return value


def integer(value, name, minimum=1):
    if type(value) is not int or value < minimum:
        raise ValueError('invalid integer parameter: ' + name)
    return value


def identifier(value):
    if not isinstance(value, str) or not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]*', value):
        raise ValueError('unsafe sample/panel identifier')
    return value


def coverage(bam, bed, output, *, executable, processes, min_mapq):
    integer(processes, 'processes')
    integer(min_mapq, 'min_mapq', 0)
    return Command((executable, 'coverage', str(bam), str(bed), '-o', str(output),
                    '-p', str(processes), '-q', str(min_mapq)), (str(output),),
                   ('chromosome', 'start', 'end', 'gene', 'log2', 'depth'))


def pooled_reference(normal_coverages, fasta, output, *, executable):
    files = [str(Path(p).resolve()) for p in normal_coverages]
    if not files or len(set(files)) != len(files):
        raise ValueError('normal coverage inputs must be nonempty and unique')
    return Command((executable, 'reference', *files, '--fasta', str(fasta), '-o', str(output)),
                   (str(output),), ('chromosome', 'start', 'end', 'gene', 'log2'))


def fix(target, antitarget, reference, output, *, executable):
    # CNVkit's native hybrid interface includes an antitarget table even if empty.
    return Command((executable, 'fix', str(target), str(antitarget), str(reference),
                    '-o', str(output)), (str(output),),
                   ('chromosome', 'start', 'end', 'gene', 'log2', 'weight'))


def segment(cnr, output, *, executable, processes, threshold, rscript=None):
    integer(processes, 'processes')
    positive(threshold, 'CBS threshold')
    if threshold >= 1:
        raise ValueError('CBS threshold must be below 1')
    args = (executable, 'segment', str(cnr), '--method', 'cbs',
            '--threshold', str(threshold), '-p', str(processes), '-o', str(output))
    if rscript is not None:
        if not isinstance(rscript,(str,Path)) or not Path(rscript).is_absolute():
            raise ValueError('explicit Rscript path must be absolute')
        args += ('--rscript-path',str(rscript))
    return Command(args,
                   (str(output),), ('chromosome', 'start', 'end', 'log2', 'probes'))


def snp_pileup(bam, vcf, output, *, executable, mapq, baseq, pseudo_spacing,
               max_depth, count_orphans, min_read_count):
    for name, value in [('mapq', mapq), ('baseq', baseq), ('min_read_count', min_read_count)]:
        integer(value, name, 0)
    integer(pseudo_spacing, 'pseudo_spacing')
    integer(max_depth, 'max_depth')
    if type(count_orphans) is not bool:
        raise ValueError('count_orphans must be an explicit boolean')
    # -x would disable overlap detection: deliberately absent.
    args = [executable, '-g', '-q', str(mapq), '-Q', str(baseq), '-P', str(pseudo_spacing),
            '-d', str(max_depth), '-r', str(min_read_count)]
    if count_orphans:
        args.append('-A')
    args.extend([str(vcf), str(output), str(bam)])
    return Command(tuple(args), (str(output),),
                   ('Chromosome', 'Position', 'Ref', 'Alt', 'File1R', 'File1A', 'File1E', 'File1D'))


def facets_ch(pileup, db, output_dir, sample_id, panel, sex, *, rscript, wrapper,
              bin_size, cval, gamma):
    identifier(sample_id)
    identifier(panel)
    if sex not in ('Male', 'Female'):
        raise ValueError('FACETS-CH requires resolved reference sex; unknown is not female')
    integer(bin_size, 'bin_size')
    positive(cval, 'cval')
    positive(gamma, 'gamma')
    # wrapper must be the audited GRCh38 adaptation, never the default upstream CLI.
    args = (rscript, str(wrapper), '--db', str(db), '--pu', str(pileup), '--sex', sex,
            '--outdir', str(output_dir), '--prefix', sample_id, '--panel', panel,
            '--bin', str(bin_size), '--cval', str(cval), '--gamma', str(gamma),
            '--blacklist', 'TRUE', '--bestnormal', 'TRUE', '--correctbaf', 'TRUE',
            '--snpdistance', 'TRUE')
    return Command(args, tuple(str(Path(output_dir)/(sample_id+suffix))
                               for suffix in ('_seg.tsv', '_sig.tsv', '.png')))


def joint_fit(matrix, bias, config, chromosome_map, output_dir, *, rscript, wrapper, library):
    """Current adapted R interface. Inputs must first pass verify_joint_inputs.

    output_dir is a new child of the native task directory, not the runner's
    receipt directory. R creates it and preserves incomplete attempts.
    """
    paths = (matrix,bias,config,chromosome_map,output_dir,rscript,wrapper,library)
    if any(not isinstance(p,(str,Path)) or not Path(p).is_absolute() for p in paths):
        raise ValueError('joint fit requires absolute paths')
    argv = (str(rscript),'--vanilla',str(wrapper),'--matrix',str(matrix),'--bias',str(bias),
            '--config',str(config),'--chromosome-map',str(chromosome_map),
            '--output',str(output_dir),'--library',str(library))
    return Command(argv,tuple(str(Path(output_dir)/name) for name in
        ('segments.tsv','locus_audit.tsv','exclusions.tsv','signal.tsv',
         'native_fit.rds','R_session.txt','COMPLETE.txt')))
