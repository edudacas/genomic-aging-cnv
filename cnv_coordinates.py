"""Coordinate checks against the actual supplied reference, without contig renaming."""
from pathlib import Path


def fai_lengths(path):
    lengths = {}
    with Path(path).open() as stream:
        for line in stream:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 5 or fields[0] in lengths or int(fields[1]) <= 0:
                raise ValueError('invalid/duplicate FASTA index entry')
            lengths[fields[0]] = int(fields[1])
    if not lengths:
        raise ValueError('empty FASTA index')
    return lengths


def bed_intervals(path, lengths, *, allow_empty=False):
    rank = {chrom: i for i, chrom in enumerate(lengths)}
    previous = None
    intervals = []
    with Path(path).open() as stream:
        for number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f'BED row {number}: missing fields')
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if chrom not in lengths or not 0 <= start < end <= lengths[chrom]:
                raise ValueError(f'BED row {number}: outside reference coordinates')
            if previous and (rank[chrom] < rank[previous[0]] or
                             (chrom == previous[0] and start < previous[2])):
                raise ValueError(f'BED row {number}: unsorted/overlapping intervals')
            previous = (chrom, start, end)
            intervals.append(previous)
    if not intervals and not allow_empty:
        raise ValueError('empty target territory')
    return intervals


def bam_header(bam_path, index_path, lengths, expected_sm):
    import pysam
    if not isinstance(expected_sm, str) or not expected_sm:
        raise ValueError('exact BAM SM required')
    with pysam.AlignmentFile(str(bam_path), 'rb', index_filename=str(index_path),
                             require_index=True) as bam:
        if not bam.check_index():
            raise ValueError('BAM index unavailable')
        header = bam.header.to_dict()
        if header.get('HD', {}).get('SO') != 'coordinate':
            raise ValueError('coordinate-sorted BAM declaration required')
        for name, length in zip(bam.references, bam.lengths):
            if lengths.get(name) != length:
                raise ValueError('BAM/reference dictionary mismatch: ' + name)
        groups = header.get('RG', [])
        if not groups or any(g.get('SM') != expected_sm for g in groups):
            raise ValueError('BAM read-group sample mismatch')
    return {'bam_header_checked': True, 'whole_bam_integrity_checked': False}


def dense_snp_reference(vcf_path, fasta_path, *, expected_contigs):
    import pysam
    rank = {chrom: i for i, chrom in enumerate(expected_contigs)}
    previous = None
    count = 0
    with pysam.FastaFile(str(fasta_path)) as fasta, pysam.VariantFile(str(vcf_path)) as vcf:
        for variant in vcf:
            if variant.contig not in rank:
                raise ValueError('SNP contig absent from reference')
            key = (rank[variant.contig], variant.pos)
            if previous is not None and key <= previous:
                raise ValueError('SNP positions unsorted or duplicated')
            previous = key
            alleles = (variant.ref,) + (variant.alts or ())
            if len(alleles) != 2 or any(a not in ('A', 'C', 'G', 'T') for a in alleles) or alleles[0] == alleles[1]:
                raise ValueError('dense SNP resource must be biallelic SNVs')
            if fasta.fetch(variant.contig, variant.start, variant.stop).upper() != variant.ref:
                raise ValueError('SNP REF does not match reference FASTA')
            count += 1
    if not count:
        raise ValueError('empty SNP resource')
    return {'biallelic_snps_checked': count, 'reference_alleles_verified': True,
            'marker_density_calibrated': False}
