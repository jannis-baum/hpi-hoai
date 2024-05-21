from typing import Literal
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

# seqname_filter matches gtf seqname, e.g. chromosome
# read_length is average read_length from RNA-Seq
def get_splice_distributions(gtf_path, seqname_filter = None, read_length=75):
    # --------------------------------------------------------------------------
    # READING DATA -------------------------------------------------------------

    gtf_df = pd.read_csv(
        gtf_path, sep='\t', comment='#', low_memory=False,
        names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    )
    if seqname_filter:
        gtf_df = gtf_df[gtf_df.seqname == seqname_filter]

    # get transcripts and exons
    transcripts = gtf_df[gtf_df['feature'] == 'transcript'].loc[:]
    exons = gtf_df[gtf_df['feature'] == 'exon'].loc[:]

    # parse attributes
    _extract_attribute(transcripts, 'gene_id')
    _extract_attribute(transcripts, 'transcript_id')
    _extract_attribute(transcripts, 'FPKM', float)
    _extract_attribute(transcripts, 'cov', float)
    del transcripts['attribute']
    _extract_attribute(exons, 'transcript_id')
    _extract_attribute(exons, 'gene_id')
    del exons['attribute']

    # --------------------------------------------------------------------------
    # EXPRESSION PROBABILITIES -------------------------------------------------

    # get total FPKM score of gene (we would get the same result with TPM)
    transcripts['FPKM_gene_total'] = transcripts.groupby('gene_id')['FPKM'].transform('sum')
    # get relative expression for transcript, i.e. what is the probability of
    # this transcript given the gene is transcribed?
    transcripts['transcript_expression'] = transcripts['FPKM'] / transcripts['FPKM_gene_total']

    # --------------------------------------------------------------------------
    # GENE READ COUNTS ---------------------------------------------------------

    # compute length of exons
    exons['length'] = exons.end - exons.start
    # sum up exons to get post-splicing length of transcript
    transcripts = pd.merge(transcripts, exons.groupby('transcript_id').length.sum(), on='transcript_id')
    # compute estimated read count for transcripts according to stringtie's
    # prepDE formula
    # (https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual#emode)
    transcripts['t_read_count'] = np.round(transcripts['cov'] * transcripts['length'] / read_length).astype(int)
    # sum up transcript read counts to get gene read counts
    transcripts = pd.merge(transcripts, transcripts.groupby('gene_id')['t_read_count'].sum().rename('g_read_count'), on='gene_id')

    # --------------------------------------------------------------------------
    # FINAL DATA ACCUMULATION --------------------------------------------------

    # merge relative expression of transcript into exons
    exons = pd.merge(exons, transcripts[['transcript_id', 'transcript_expression']], on='transcript_id')
    # directionally independent sides
    exons['5p'] = np.where(exons['strand'] == '+', exons['start'], exons['end'])
    exons['3p'] = np.where(exons['strand'] == '+', exons['end'], exons['start'])

    # helper to avoid repeating the same code for both splicing sides
    def _get_splice_probs(site: Literal['acceptor', 'donor']):
        side = '5p' if site == 'acceptor' else '3p'
        # sum up transcript expression to get splice probability
        sites = exons.groupby([side, 'gene_id'])[['transcript_expression']].sum()
        # (keep position as column instead of just index to preserve it in merge)
        sites.insert(0, side, sites.index.get_level_values(side))
        # merge estimated gene read counts into sites & rename columns
        sites = pd.merge(
            sites,
            transcripts.drop_duplicates(subset=['gene_id'])[['gene_id', 'g_read_count']],
            on='gene_id'
        ).rename(columns={
            'transcript_expression': 'p',
            'g_read_count': 'n'
        })
        # check that we don't have duplicate sites
        dup_sites = sites[side].duplicated(keep=False)
        if sum(dup_sites) > 0:
            print(f'Warning: dropping duplicate {site} sites:')
            print(sites[dup_sites])
            sites = sites[~dup_sites]
        # delete unneeded gene_id column
        del sites['gene_id']
        return sites

    return (
        _get_splice_probs('acceptor'),
        _get_splice_probs('donor')
    )
