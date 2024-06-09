from argparse import ArgumentParser
import json
import os

import pandas as pd

def filetype(file) -> str | None:
    if 'file_name' not in file or 'data_format' not in file or 'experimental_strategy' not in file:
        return None
    # RNA-Seq conditions
    if all([
        file['file_name'].endswith('.rna_seq.transcriptome.gdc_realn.bam'),
        file['data_format'] == 'BAM',
        file['experimental_strategy'] == 'RNA-Seq'
    ]):
        return 'rna_seq'
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
    parser = ArgumentParser('filter_files', description='Filter GDC manifest file based on cat JSON')
    parser.add_argument('manifest', type=str)
    parser.add_argument('json', type=str)

    args = parser.parse_args()

    json_dir = os.path.dirname(os.path.abspath(args.json))
    out_info_path = os.path.join(json_dir, 'filtered_info.csv')
    manifest_dir = os.path.dirname(os.path.abspath(args.manifest))
    out_manifest_path = os.path.join(manifest_dir, 'filtered_manifest.txt')

    with open(args.json, 'rb') as fp:
        files = json.load(fp)
    
    info = pd.DataFrame(columns=['case_id', 'rna_seq', 'wgs_brass', 'wgs_pindel', 'wgs_caveman'])
    info.set_index('case_id', inplace=True)
    for file in files:
        ft = filetype(file)
        cid = case_id(file)
        if not ft or not cid: continue
        if cid not in info.index:
            info.loc[cid] = [None, None, None, None]
        info.at[cid, ft] = file['file_name']

    info.to_csv(out_info_path)
