from pyfaidx import Fasta
from pysam import VariantFile

class Annotator():
    def __init__(self, fasta_path, vcf_path):
        self._fasta = Fasta(fasta_path)
        self._vcf = VariantFile(vcf_path)

    def get_seq(self, chromosome: str, start: int, end: int) -> str:
        ref_seq = list(self._fasta[chromosome][start:end].seq)

        for record in self._vcf.fetch(chromosome, start, end):
            pos = record.pos - 1 - start  # -1 for 0-based list index
            ref = record.ref
            alt = record.alts[0] # assume first alternative allele
            # variant
            ref_seq = ref_seq[:pos] + list(alt) + ref_seq[pos + len(ref):]

        return ''.join(ref_seq)
