import pandas as pd
from pyensembl.gene import Gene

from script.data_loading import find_path, study_info
from script.gene import chrom

def get_junctions(case_id: str, gene: Gene) -> pd.DataFrame:
    case = study_info.loc[case_id]

    df = pd.read_csv(find_path(case['splice']), compression='gzip', sep='\t')
    # filter splice junction quantification for given gene
    strand = 1 if gene.strand == '+' else 2
    df = df[(df['#chromosome'] == chrom(gene))
        & (df['intron_start'] > gene.start) & (df['intron_start'] < gene.end) \
        & (df['intron_end'] > gene.start) & (df['intron_end'] < gene.end) \
        & (df['strand'] == strand)
    ]
    return pd.DataFrame({
        'intron_donor': df['intron_start'],
        'intron_acceptor': df['intron_end'],
        'n': df['n_unique_map']
    }).set_index(['intron_donor', 'intron_acceptor'])
