from pyensembl.gene import Gene
from pyfaidx import Fasta
from pysam import VariantFile

from script.data_loading import find_path
from script.gene import chrom

class Annotator():
    _complement_base = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }

    def __init__(self, vcf_path: str, fasta_path: str = find_path('ch38.fa')):
        self._fasta = Fasta(fasta_path)
        self._vcf = VariantFile(vcf_path)

    # get sequence with variants and reference indices for alignment:
    # reference indices correspond to elements of sequence and denote which
    # index/indices they correspond to in the reference genome
    def get_seq(self, chromosome: str, start: int, end: int, rc: bool) -> tuple[str, list[list[int]]]:
        seq: str = self._fasta.get_seq(chromosome, start, end + 1).seq
        idx = [[i] for i in range(start, end + 1)]

        # sort records by position in descending order to not mess up indices
        records = sorted(self._vcf.fetch(chromosome, start, end), key=lambda x: x.pos, reverse=True)
        for record in records:
            ref: str = record.ref
            alt: str = record.alts[0] # assume first alternative allele

            r_start = record.pos - 1 - start  # -1 for 0-based list index
            r_end = r_start + len(ref)

            # insert variance into sequence
            seq = seq[:r_start] + alt + seq[r_end:]
            # keep track of corresponding indices
            # (we have to reshift r_start/end by the original start position)
            idx = idx[:r_start] + [[i for i in range(r_start + start, r_end + start)]] * len(alt) + idx[r_end:]

        if rc:
            seq = [Annotator._complement_base[base] for base in seq[::-1]]
            idx = idx[::-1]
        return (''.join(seq), idx)

    def get_gene_seq(self, gene: Gene, padding: int = 5000) -> tuple[str, list[list[int]]]:
        return self.get_seq(chrom(gene), gene.start - padding, gene.end + padding, (gene.strand == '-'))
