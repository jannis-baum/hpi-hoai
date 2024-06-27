import numpy as np
import pandas as pd

from script.data_loading import get_path, study_info
from script.gene import Gene

def junctions_and_count(case_id: str, gene: Gene) -> tuple[pd.DataFrame, int]:
    case = study_info.loc[case_id]

    def _get_junctions() -> pd.DataFrame:
        df = pd.read_csv(get_path(case['splice']), compression='gzip', sep='\t')
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

    def _get_count() -> int:
        df = pd.read_csv(get_path(case['expression']), sep='\t', skiprows=1, index_col='gene_id')
        return df.loc[gene.id]['unstranded']

    return (_get_junctions(), _get_count())
