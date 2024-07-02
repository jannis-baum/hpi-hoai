from pyensembl.gene import Gene


def chrom(gene: Gene):
    return f'chr{gene.contig}'
