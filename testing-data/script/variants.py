from pyfaidx import Fasta
from pysam import VariantFile

def get_first_chromosome(vcf_path: str):
    vcf = VariantFile(vcf_path, 'r')
    first_record = next(vcf)
    return first_record.chrom

def filter_positions(vcf_in_path: str, vcf_out_path: str, whitelist: set[int]):
    vcf_in = VariantFile(vcf_in_path, 'r')
    vcf_out = VariantFile(vcf_out_path, 'w', header=vcf_in.header)
    for record in vcf_in:
        if record.pos in whitelist:
            vcf_out.write(record)

class Annotator():
    def __init__(self, fasta_path, vcf_path):
        self._fasta = Fasta(fasta_path)
        self._vcf = VariantFile(vcf_path)

    def get_seq(self, chromosome, start, end):
        ref_seq = list(self._fasta[chromosome][start:end].seq)

        for record in self._vcf.fetch(chromosome, start, end):
            pos = record.pos - 1 - start  # -1 for 0-based list index
            ref = record.ref
            alt = record.alts[0] # assume first alternative allele
            # variant
            ref_seq = ref_seq[:pos] + list(alt) + ref_seq[pos + len(ref):]

        return ''.join(ref_seq)
