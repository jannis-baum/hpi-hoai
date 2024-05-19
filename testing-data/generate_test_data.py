import argparse

import numpy as np

from script.splice_distributions import get_splice_distributions
from script.variants import get_first_chromosome, filter_positions

if __name__ == '__main__':
    parser = argparse.ArgumentParser('generate_test_data')
    parser.add_argument('gtf_in', help='Input GTF file from Stringtie')
    parser.add_argument('vcf_in', help='Single chromosome VCF file to filter for relevant variants')
    parser.add_argument('-D', '--distance', help='Maximum considered distance between variant and splice site', type=int, default=50)
    args = parser.parse_args()

    chrom = get_first_chromosome(args.vcf_in)

    p_acceptor, p_donor = get_splice_distributions(args.gtf_in, chrom)
    p_acceptor.to_csv('p_acceptor.csv')
    p_donor.to_csv('p_donor.csv')

    def _considered_range(pos) -> list[int]:
        return list(range(pos - args.distance, pos + args.distance))

    pos_whitelist = {
        considered_pos
        for pos in np.concatenate([p_acceptor.index, p_donor.index])
        for considered_pos in _considered_range(pos)
    }
    filter_positions(args.vcf_in, 'filtered.vcf', pos_whitelist)
