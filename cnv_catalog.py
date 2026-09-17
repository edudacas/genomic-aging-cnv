"""Frozen whole-reference SNP/depth catalog, generated without sample data.

This defines expected coordinates, not measured/evaluable coverage. SNPs take
precedence when a SNP coincides with the depth grid. Grid positions are 1-based
multiples of the native snp-pileup -P spacing. No LD pruning or contig remapping.
"""
from cnv_commands import integer
from cnv_pileup import position_key


def merge_depth_grid(snps, lengths, spacing):
    integer(spacing, 'pseudo-SNP spacing')
    if not lengths or any(type(n) is not int or n <= 0 for n in lengths.values()):
        raise ValueError('nonempty positive reference lengths required')
    rank = {chrom: i for i, chrom in enumerate(lengths)}

    def checked_snps():
        previous = None
        for locus in snps:
            key = position_key(locus['chromosome'], locus['position'], lengths, rank)
            if previous is not None and key <= previous:
                raise ValueError('SNP catalog unsorted or duplicate')
            if (locus.get('kind') != 'SNP' or locus['ref'] not in ('A', 'C', 'G', 'T') or
                    locus['alt'] not in ('A', 'C', 'G', 'T') or locus['ref'] == locus['alt']):
                raise ValueError('biallelic SNP required')
            previous = key
            yield key, locus

    iterator = iter(checked_snps())
    snp = next(iterator, None)
    for chrom, length in lengths.items():
        for position in range(spacing, length + 1, spacing):
            grid_key = rank[chrom], position
            while snp is not None and snp[0] < grid_key:
                yield snp[1]
                snp = next(iterator, None)
            if snp is not None and snp[0] == grid_key:
                yield snp[1]
                snp = next(iterator, None)
            else:
                yield dict(chromosome=chrom, position=position, ref='.', alt='.', kind='DEPTH')
    while snp is not None:
        yield snp[1]
        snp = next(iterator, None)


def reference_catalog(vcf_path, fasta_path, spacing):
    """Validate SNP REF against the actual FASTA while streaming the catalog."""
    import pysam
    with pysam.FastaFile(str(fasta_path)) as fasta, pysam.VariantFile(str(vcf_path)) as vcf:
        lengths = dict(zip(fasta.references, fasta.lengths))

        def snps():
            for variant in vcf:
                alleles = variant.alleles or ()
                if len(alleles) != 2 or any(a not in ('A', 'C', 'G', 'T') for a in alleles):
                    raise ValueError('VCF must contain only biallelic SNPs')
                if variant.contig not in lengths or not 1 <= variant.pos <= lengths[variant.contig]:
                    raise ValueError('VCF coordinate outside reference')
                if fasta.fetch(variant.contig, variant.pos-1, variant.pos).upper() != variant.ref:
                    raise ValueError('VCF REF/reference mismatch')
                yield dict(chromosome=variant.contig, position=variant.pos,
                           ref=variant.ref, alt=variant.alts[0], kind='SNP')
        yield from merge_depth_grid(snps(), lengths, spacing)
