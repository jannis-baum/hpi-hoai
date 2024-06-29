import pandas as pd

from script.data_loading import find_path, study_info
from script.gene import Gene

def get_junctions(case_id: str, gene: Gene) -> pd.DataFrame:
    case = study_info.loc[case_id]

    df = pd.read_csv(find_path(case['splice']), compression='gzip', sep='\t')
    # filter splice junction quantification for given gene
    df = df[(df['#chromosome'] == gene.chromosome)
        & (df['intron_start'] > gene.bp_start) & (df['intron_start'] < gene.bp_end) \
        & (df['intron_end'] > gene.bp_start) & (df['intron_end'] < gene.bp_end) \
        & (df['strand'] != 0)
    ]
    return pd.DataFrame({
        'intron_donor': df['intron_start'],
        'intron_acceptor': df['intron_end'],
        'n': df['n_unique_map']
    }).set_index(['intron_donor', 'intron_acceptor'])
