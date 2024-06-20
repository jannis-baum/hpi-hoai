import numpy as np
import pandas as pd

def get_brca1(df):
    brca1_start = 43044295
    brca1_end = 43170245
    condition = (df['#chromosome'] == 'chr17') \
        & (df['intron_start'] > brca1_start) & (df['intron_start'] < brca1_end) \
        & (df['intron_end'] > brca1_start) & (df['intron_end'] < brca1_end) \
        & (df['strand'] != 0)
    return df[condition]

def get_junctions(df):
    strand_conditions = [df['strand'] == 1, df['strand'] == 2]
    junctions = pd.DataFrame({
        'intron_donor': np.select(strand_conditions, [df['intron_start'], df['intron_end']], default=np.nan),
        'intron_acceptor': np.select(strand_conditions, [df['intron_end'], df['intron_start']], default=np.nan),
        'n': df['n_unique_map']
    }).set_index(['intron_donor', 'intron_acceptor'])
    return junctions
