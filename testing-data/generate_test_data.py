import argparse
import os

import numpy as np

from script.splice_distributions import get_splice_distributions
from script.variants import get_first_chromosome, filter_positions

if __name__ == '__main__':
    parser = argparse.ArgumentParser('generate_test_data')
    parser.add_argument('gtf_in', help='Input GTF file from Stringtie')
    parser.add_argument('vcf_in', help='Single chromosome VCF file to filter for relevant variants')
    parser.add_argument('-D', '--distance', help='Maximum considered distance between variant and splice site, default 50', type=int, default=50)
    parser.add_argument('-o', '--output', help='Output directory, default cwd', default='.')
    args = parser.parse_args()

    out_dir = args.output
    os.makedirs(out_dir, exist_ok=True)

    chrom = get_first_chromosome(args.vcf_in)

    p_acceptor, p_donor = get_splice_distributions(args.gtf_in, chrom)
    p_acceptor.to_csv(os.path.join(out_dir, 'p_acceptor.csv'), index=False)
    p_donor.to_csv(os.path.join(out_dir, 'p_donor.csv'), index=False)

    def _considered_range(pos) -> list[int]:
        return list(range(pos - args.distance, pos + args.distance))

    pos_whitelist = {
        considered_pos
        for pos in np.concatenate([p_acceptor['5p'], p_donor['3p']])
        for considered_pos in _considered_range(pos)
    }
    filter_positions(args.vcf_in, os.path.join(out_dir, 'filtered-variants.vcf'), pos_whitelist)
