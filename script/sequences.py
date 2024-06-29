from pyfaidx import Fasta
from pysam import VariantFile

from script.data_loading import get_path
from script.gene import Gene

class Annotator():
    def __init__(self, vcf_path: str, fasta_path: str = get_path('ch38.fa')):
        self._fasta = Fasta(fasta_path)
        self._vcf = VariantFile(vcf_path)

    # get sequence with variants and reference indices for alignment:
    # reference indices correspond to elements of sequence and denote which
    # index/indices they correspond to in the reference genome
    def get_seq(self, chromosome: str, start: int, end: int) -> tuple[str, list[list[int]]]:
        seq = list(self._fasta[chromosome][start:end + 1].seq)
        idx = [[i] for i in range(start, end + 1)]

        # sort records by position in descending order to not mess up indices
        records = sorted(self._vcf.fetch(chromosome, start, end), key=lambda x: x.pos, reverse=True)
        for record in records:
            ref = record.ref
            alt = record.alts[0] # assume first alternative allele

            r_start = record.pos - 1 - start  # -1 for 0-based list index
            r_end = r_start + len(ref)

            # insert variance into sequence
            seq = seq[:r_start] + list(alt) + seq[r_end:]
            # keep track of corresponding indices
            # (we have to reshift r_start/end by the original start position)
            idx = idx[:r_start] + [[i for i in range(r_start + start, r_end + start)]] * len(alt) + idx[r_end:]

        return (''.join(seq), idx)

    def get_gene_seq(self, gene: Gene, padding: int = 5000) -> tuple[str, list[list[int]]]:
        return self.get_seq(gene.chromosome, gene.bp_start - padding, gene.bp_end + padding)

