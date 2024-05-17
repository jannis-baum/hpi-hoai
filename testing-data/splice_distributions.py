import argparse

import numpy as np
import pandas as pd

# helper to parse attributes and add them as a new column
def _extract_attribute(df, attr):
    def _parse_attribute(col, attr):
        # (some attributes are optional, that's why we need NaN)
        try: value = col.split(f'{attr} "')[1].split('";')[0]
        except: return np.NaN

        try: return float(value)
        except: return value
    df[attr] = df['attribute'].apply(_parse_attribute, args=(attr, ))

if __name__ == '__main__':
    parser = argparse.ArgumentParser('splice_distributions')
    parser.add_argument('gtf_in', help='Input GTF file from Stringtie')
    args = parser.parse_args()

    gtf_df = pd.read_csv(
        args.gtf_in, sep="\t", comment='#', low_memory=False,
        names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    )

    # get transcripts
    transcripts = gtf_df[gtf_df['feature'] == 'transcript'].loc[:]

    # parse attributes
    for attr in ['gene_id', 'transcript_id', 'FPKM']:
        _extract_attribute(transcripts, attr)
    del transcripts['attribute']

    # get total FPKM score of gene (we would get the same result with TPM)
    transcripts['FPKM_gene_total'] = transcripts.groupby('gene_id')['FPKM'].transform('sum')
    # get relative expression for transcript, i.e. what is the probability of this transcript given the gene is transcribed?
    transcripts['transcript_expression'] = transcripts['FPKM'] / transcripts['FPKM_gene_total']

    # get exons
    exons = gtf_df[gtf_df['feature'] == 'exon'].loc[:]

    # parse attributes
    _extract_attribute(exons, 'transcript_id')
    del exons['attribute']

    # merge relative expression of transcript into exons
    exons = pd.merge(exons, transcripts[['transcript_id', 'transcript_expression']], on='transcript_id')
    # sum up relative expressions of transcripts to receive donor/acceptor probabilities & save as csvs
    exons.groupby('start')['transcript_expression'].sum().rename('p_acceptor').to_csv('p_acceptor.csv')
    exons.groupby('end')['transcript_expression'].sum().rename('p_donor').to_csv('p_donor.csv')
