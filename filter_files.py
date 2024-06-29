from argparse import ArgumentParser
import json
import os
import sys

import pandas as pd

def filetype(file) -> str | None:
    if 'file_name' not in file or 'data_format' not in file or 'experimental_strategy' not in file:
        return None
    # Splice Junctions & Gene Expression conditions
    if all([
        file['data_format'] == 'TSV',
        file['experimental_strategy'] == 'RNA-Seq'
    ]):
        if file['file_name'].endswith('.rna_seq.star_splice_junctions.tsv.gz'):
            return 'splice'
        if file['file_name'].endswith('.rna_seq.augmented_star_gene_counts.tsv'):
            return 'expression'
    # WGS conditions
    if not all([
        file['data_format'] == 'VFC',
        file['experimental_strategy'] == 'WGS', 
    ]):
        if file['file_name'].endswith('.wgs.BRASS.raw_structural_variation.vcf.gz'): return 'wgs_brass'
        if file['file_name'].endswith('.wgs.sanger_raw_pindel.raw_somatic_mutation.vcf.gz'): return 'wgs_pindel'
        if file['file_name'].endswith('.wgs.CaVEMan.raw_somatic_mutation.vcf.gz'): return 'wgs_caveman'

    return None

def case_id(file) -> str | None:
    cases = file['cases']
    if len(cases) != 1: return None
    return cases[0]['case_id']

if __name__ == '__main__':
    parser = ArgumentParser('filter_files', description='Filter GDC cart based on JSON')
    parser.add_argument('json', type=str)
    parser.add_argument('--out', type=str, help='Info file output')

    args = parser.parse_args()

    json_dir = os.path.dirname(os.path.abspath(args.json))
    out_info_path = args.out if args.out else os.path.join(json_dir, 'filtered_info.csv')

    # inspect JSON to find relevant files
    with open(args.json, 'rb') as fp:
        files = json.load(fp)
    
    info = pd.DataFrame(columns=['case_id', 'splice', 'expression', 'wgs_brass', 'wgs_pindel', 'wgs_caveman', 'total_size'])
    info.set_index('case_id', inplace=True)
    for file in files:
        ft = filetype(file)
        cid = case_id(file)
        if not ft or not cid: continue
        if cid not in info.index:
            info.loc[cid] = [pd.NA, pd.NA, pd.NA, pd.NA, pd.NA, 0]
        info.at[cid, ft] = file['file_name']
        info.at[cid, 'total_size'] = info.at[cid, 'total_size'] + file['file_size']

    info.dropna(how='any', inplace=True)

    csv = info.to_csv()
    if args.out:
        with open(args.out, 'w') as fp:
            fp.write(csv)
    else:
        print(csv, end='')

    info['total_size'] = info['total_size'] / 10**9
    # print to stderr so we can still redirect csv to file
    print('Total size (GB): ', info['total_size'].sum(), file=sys.stderr)
