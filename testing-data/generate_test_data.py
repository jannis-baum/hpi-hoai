import argparse

from script.splice_distributions import get_splice_distributions

if __name__ == '__main__':
    parser = argparse.ArgumentParser('generate_test_data')
    parser.add_argument('gtf_in', help='Input GTF file from Stringtie')
    parser.add_argument('-chr', '--chromosome', type=str, help='output only data for the given chromosome (or, more accurately, GTF seqname)')
    args = parser.parse_args()

    p_acceptor, p_donor = get_splice_distributions(args.gtf_in, args.chromosome)
    p_acceptor.to_csv('p_acceptor.csv')
    p_donor.to_csv('p_donor.csv')
