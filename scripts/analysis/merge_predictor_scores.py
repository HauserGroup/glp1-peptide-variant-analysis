"""

"""
import pathlib
import pandas as pd 

""" Create a seaborn clustered heatmap for a dataframe.
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pathlib as pl
from matplotlib.lines import Line2D
# made svg edited
plt.rcParams['svg.fonttype'] = 'none'



def _merge_individual_scores(data_dir):
    dfs = [data_dir / (name + ".csv") for name in ['boltz2_iptm', 'dG_cross', 'alphamissense_pathogenicity', 'chai_iptm', 'af2_plddt', 'af2_rosettarelax', 'per_residue_energy_int', 'predicted_potency', 'FEP', 'dG_separated', 'old_FEP', 'ddG', 'proteocast', 'finetuned', 'mean_binding', 'af2_iptm', 'af3_iptm', 'mean_cAMP', 'mean_arrestin', 'sc_value']]
    missing = [p for p in dfs if not p.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing benchmark metric inputs: {missing}")
    cols_to_drop = ['mean_binding_factor', 'source', 'sequence', 'pdb_path']
    merged_df = None
    for p in dfs:
        df = pd.read_csv(p)
        print(f"Read {p}, columns: {df.columns}")
        if merged_df is None:
            merged_df = df
        else:
            df = df.drop(columns=[c for c in cols_to_drop if c in df.columns])
            df = df.drop(columns=[c for c in df.columns if 'unnamed' in c.lower()])
            print(f"Merging {p} with merged_df, on 'peptide', how='outer'")
            # print(f'old columns: {merged_df.columns}')
            # print(f'new columns: {df.columns}')
            # drop columns from df that are already in merged_df, except for 'peptide'
            df = df.drop(columns=[c for c in df.columns if c in merged_df.columns and c != 'peptide'])

            merged_df = merged_df.merge(df, on='peptide', how='outer')
            
    return merged_df


def run():
    ""
    repo_root = pathlib.Path(__file__).resolve().parents[2]
    data_dir = repo_root / 'data' / 'benchmark' / 'metric-inputs'
    out_p = repo_root / 'output' / 'analysis' / 'predictor-benchmark' / 'g1_metrics.csv'
    out_p.parent.mkdir(parents=True, exist_ok=True)
    sources_to_use = [
        'human_variants', 
        'alanine_scan', 
        ]
    max_peptide_length = 32

    # merge .csv files on 'peptide'
    merged_df = _merge_individual_scores(data_dir)
    print(merged_df['source'].value_counts())
    # drop on length
    print(f'Dropping peptides with length > {max_peptide_length}, n before={len(merged_df)}')
    merged_df = merged_df[merged_df['sequence'].str.len() <= max_peptide_length]
    # only keep some sources
    _sources_dropped = set(merged_df['source'].unique()) - set(sources_to_use)
    print(merged_df['source'].value_counts())
    print(f"Dropping sources: {_sources_dropped}")
    for _to_drop in _sources_dropped:
        print(f"\tDropping source: {_to_drop}, n={len(merged_df[merged_df['source'] == _to_drop])}")
        merged_df = merged_df[merged_df['source'] != _to_drop]
    for _source in sources_to_use:
        print(f"Keeping source: {_source}, n={len(merged_df[merged_df['source'] == _source])}")
        
    # drop reference
    ### assert 'HV_GLP1_7-37' in merged_df['peptide'].values, "Reference peptide 'HV_GLP1_7-37' not found in data."
    if 'HV_GLP1_7-37' in merged_df['peptide'].values:
        merged_df = merged_df[merged_df['peptide'] != 'HV_GLP1_7-37']
    if 'HV_GLP1 7-37' in merged_df['peptide'].values:
        merged_df = merged_df[merged_df['peptide'] != 'HV_GLP1 7-37']


    print(f'n after={len(merged_df)}')
    merged_df.to_csv(out_p, index=False)
    print(f"Merged data saved to {out_p}")


if __name__ == '__main__':
    run()
