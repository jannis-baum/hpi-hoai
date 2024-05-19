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
