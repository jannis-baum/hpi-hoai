import numpy as np
import pandas as pd

# helper to parse attributes and add them as a new column
def _extract_attribute(df, attr, cast = None):
    def _parse_attribute(col, attr):
        # (some attributes are optional, that's why we need NaN)
        try: value = col.split(f'{attr} "')[1].split('";')[0]
        except: return np.NaN

        if not cast: return value
        try: return cast(value)
        except: return np.NaN
    df[attr] = df['attribute'].apply(_parse_attribute, args=(attr, ))

def get_splice_distributions(gtf_path, seqname_filter = None):
    gtf_df = pd.read_csv(
        gtf_path, sep='\t', comment='#', low_memory=False,
        names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    )
    if seqname_filter:
        gtf_df = gtf_df[gtf_df.seqname == seqname_filter]

    # get transcripts
    transcripts = gtf_df[gtf_df['feature'] == 'transcript'].loc[:]

    # parse attributes
    _extract_attribute(transcripts, 'gene_id')
    _extract_attribute(transcripts, 'transcript_id')
    _extract_attribute(transcripts, 'FPKM', float)
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
    return (
        exons.groupby('start')['transcript_expression'].sum().rename('p_acceptor'),
        exons.groupby('end')['transcript_expression'].sum().rename('p_donor')
    )
