"""Streaming exact-catalog alignment of single-sample FACETS snp-pileup files.

The catalog includes SNP and depth-grid positions, in reference order. It is
not the union of observed rows: sites absent from every sample remain present.
Exhaust the iterator before declaring validation complete.
"""
from contextlib import ExitStack
import csv
import gzip
import statistics
from pathlib import Path

from cnv_allelic import summarize_marker


COLUMNS = ('Chromosome', 'Position', 'Ref', 'Alt', 'File1R', 'File1A', 'File1E', 'File1D')


def position_key(chrom, position, lengths, rank):
    if chrom not in lengths or type(position) is not int or not 1 <= position <= lengths[chrom]:
        raise ValueError('pileup/catalog coordinate outside reference')
    return rank[chrom], position


def native_rows(stream, lengths):
    """Keep native Ref=Alt='.' depth pseudo-SNPs distinct from allelic loci."""
    reader = csv.DictReader(stream)
    if tuple(reader.fieldnames or ()) != COLUMNS:
        raise ValueError('expected native single-sample snp-pileup CSV schema')
    rank = {chrom: i for i, chrom in enumerate(lengths)}
    previous = None
    for row in reader:
        if None in row or any(value is None for value in row.values()):
            raise ValueError('ragged snp-pileup row')
        key = position_key(row['Chromosome'], int(row['Position']), lengths, rank)
        if previous is not None and key <= previous:
            raise ValueError('unsorted or duplicate snp-pileup coordinate')
        previous = key
        ref, alt = row['Ref'], row['Alt']
        pseudo = ref == alt == '.'
        if not pseudo and (ref not in 'ACGT' or alt not in 'ACGT' or
                           len(ref) != 1 or len(alt) != 1 or ref == alt):
            raise ValueError('invalid SNP alleles')
        counts = tuple(int(row[column]) for column in COLUMNS[4:])
        if any(value < 0 for value in counts):
            raise ValueError('negative snp-pileup count')
        if pseudo and counts[1:] != (0, 0, 0):
            raise ValueError('invalid native depth pseudo-SNP counts')
        yield key, (ref, alt), counts


def align_catalog(catalog, pileup_paths, lengths):
    """Yield (locus, per-sample native counts/None); memory O(number of normals).

    Catalog rows: chromosome, position (1-based int), ref, alt, kind (SNP/DEPTH).
    Every observed row must belong to the full frozen catalog. No contig
    renaming, allele flipping, LD pruning, imputation or implicit zero-fill.
    """
    if not lengths or not pileup_paths:
        raise ValueError('reference and pileup inputs required')
    paths = [str(Path(path).resolve()) for path in pileup_paths.values()]
    if len(paths) != len(set(paths)):
        raise ValueError('pileup file reused across samples')
    rank = {chrom: i for i, chrom in enumerate(lengths)}
    with ExitStack() as stack:
        iterators, heads = {}, {}
        for sample, path in pileup_paths.items():
            opener = gzip.open if str(path).endswith('.gz') else open
            stream = stack.enter_context(opener(path, 'rt', newline=''))
            iterators[sample] = native_rows(stream, lengths)
            heads[sample] = next(iterators[sample], None)
        previous = None
        count = 0
        for locus in catalog:
            key = position_key(locus['chromosome'], locus['position'], lengths, rank)
            if previous is not None and key <= previous:
                raise ValueError('unsorted or duplicate catalog coordinate')
            previous = key
            alleles = locus['ref'], locus['alt']
            if locus['kind'] == 'DEPTH':
                if alleles != ('.', '.'):
                    raise ValueError('depth catalog requires native dot alleles')
            elif locus['kind'] == 'SNP':
                if any(a not in ('A', 'C', 'G', 'T') for a in alleles) or alleles[0] == alleles[1]:
                    raise ValueError('catalog must contain biallelic SNPs')
            else:
                raise ValueError('unknown catalog kind')
            counts = {}
            for sample, head in heads.items():
                if head is not None and head[0] < key:
                    raise ValueError('observed pileup locus outside frozen catalog')
                if head is not None and head[0] == key:
                    if head[1] != alleles:
                        raise ValueError('pileup/catalog allele mismatch')
                    counts[sample] = head[2]
                    heads[sample] = next(iterators[sample], None)
                else:
                    counts[sample] = None
            count += 1
            yield locus, counts
        if not count:
            raise ValueError('empty frozen catalog')
        if any(head is not None for head in heads.values()):
            raise ValueError('trailing pileup loci outside frozen catalog')


def training_summaries(catalog, pileup_paths, lengths, pon_sample_ids, policy):
    ids = list(pon_sample_ids)
    if len(ids) != len(set(ids)) or set(pileup_paths) != set(ids):
        raise ValueError('pileups must match exactly the validated PoN membership')
    policy.validate()
    for locus, counts in align_catalog(catalog, pileup_paths, lengths):
        observed = [row for row in counts.values() if row is not None]
        if locus['kind'] == 'DEPTH':
            summary = dict(status='DEPTH_ONLY', total_normals=len(ids),
                           observed_normals=len(observed), missing_normals=len(ids)-len(observed),
                           zero_depth_normals=sum(row[0] == 0 for row in observed),
                           median_observed_depth=statistics.median(row[0] for row in observed) if observed else None,
                           correction_bias=None)
        else:
            summary = summarize_marker(
                {sample: row[:2] if row is not None else None for sample, row in counts.items()},
                ids, policy)
        # Errors and deletions remain visible; they are not alternate evidence.
        yield dict(locus, **summary,
                   native_error_count=sum(row[2] for row in observed) if observed else None,
                   native_deletion_count=sum(row[3] for row in observed) if observed else None)
