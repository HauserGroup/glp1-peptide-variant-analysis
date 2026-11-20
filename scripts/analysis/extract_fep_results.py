#!/usr/bin/env python3
"""Extract FEP summary values from three PMX replicas per mutation and state."""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
MUTATIONS = [
    'E100Q', 'T102N', 'T102I', 'S105T', 'D106A', 'S108R', 'G113V',
    'Q114P', 'A115G', 'K117N', 'F119L', 'I120N', 'I120S', 'A121S',
    'A121V', 'W122C', 'V124L', 'G126S', 'G126C', 'R127L', 'R127Q',
]
CONV_THRESHOLD = 1.0


def parse_pmx_file(file_path):
    """Parse BAR, Jarzynski, and convergence values from a PMX results.txt."""
    content = file_path.read_text()
    patterns = {
        'bar_dg': r'BAR: dG =\s+([\d\.-]+) kJ/mol',
        'jarz_dg': r'JARZ: dG Mean\s+=\s+([\d\.-]+) kJ/mol',
        'conv': r'BAR: Conv =\s+([\d\.-]+)',
    }
    values = {}
    for name, pattern in patterns.items():
        match = re.search(pattern, content)
        if match is None:
            raise ValueError(f'Missing {name} in {file_path}')
        values[name] = float(match.group(1))
    return values


def get_detailed_state_data(base_path, mutation, state):
    """Gather three replica results and calculate statistics for one state."""
    results = [
        parse_pmx_file(base_path / state / mutation / f'analyse{i}' / 'results.txt')
        for i in range(1, 4)
    ]
    bar_dgs = [result['bar_dg'] for result in results]
    return {
        'raw_bar': bar_dgs,
        'bar_avg': np.mean(bar_dgs),
        'bar_sem': np.std(bar_dgs, ddof=1) / np.sqrt(len(bar_dgs)),
        'jarz_avg': np.mean([result['jarz_dg'] for result in results]),
        'conv_max': max(abs(result['conv']) for result in results),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--data-root', type=Path,
        default=REPO_ROOT / 'zenodo-data/in-silico/03-FEP',
        help='Folder containing peptide/ and complex/ PMX replica directories',
    )
    parser.add_argument(
        '--output', type=Path,
        default=REPO_ROOT / 'output/analysis/mutation_final_report.csv',
    )
    parser.add_argument(
        '--plot', type=Path,
        default=REPO_ROOT / 'output/figures/fep_ddg_diagnostic.png',
    )
    args = parser.parse_args()
    results = []

    for mut in MUTATIONS:
        try:
            pep = get_detailed_state_data(args.data_root, mut, 'peptide')
            com = get_detailed_state_data(args.data_root, mut, 'complex')
        except (FileNotFoundError, ValueError) as exc:
            parser.error(str(exc))

        ddg_bar = com['bar_avg'] - pep['bar_avg']
        ddg_jarz = com['jarz_avg'] - pep['jarz_avg']
        total_err = np.sqrt(pep['bar_sem']**2 + com['bar_sem']**2)
        max_conv = max(pep['conv_max'], com['conv_max'])
        warning = 'OK' if max_conv <= CONV_THRESHOLD else f'HIGH CONV ({max_conv:.2f})'

        results.append({
            'ResNum': int(mut[1:-1]), 'WT': mut[0], 'Mut': mut[-1],
            'Mutation': mut,
            'Pep_Rep1': pep['raw_bar'][0], 'Pep_Rep2': pep['raw_bar'][1],
            'Pep_Rep3': pep['raw_bar'][2],
            'SEM_Peptide': round(pep['bar_sem'], 3),
            'Com_Rep1': com['raw_bar'][0], 'Com_Rep2': com['raw_bar'][1],
            'Com_Rep3': com['raw_bar'][2],
            'SEM_Complex': round(com['bar_sem'], 3),
            'BAR_ddG': round(ddg_bar, 3),
            'JARZ_ddG': round(ddg_jarz, 3),
            'Total_Propagated_Error': round(total_err, 3),
            'Max_Abs_Conv': round(max_conv, 3),
            'Status': warning,
        })

    df = pd.DataFrame(results)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)

    fig, ax = plt.subplots(figsize=(12, 6))
    colors = ['red' if status != 'OK' else 'skyblue' for status in df['Status']]
    ax.bar(df['Mutation'], df['BAR_ddG'], yerr=df['Total_Propagated_Error'],
           capsize=5, color=colors, label='BAR ddG')
    ax.set_ylabel(r'$\Delta\Delta G$ (kJ/mol)')
    ax.set_title('Binding Free Energy: Individual Errors Propagated to Total Error')
    ax.axhline(0, color='black', linewidth=0.8)
    args.plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.plot, bbox_inches='tight')
    plt.close(fig)
    print(f'Wrote {len(df)} mutations to {args.output} and {args.plot}')


if __name__ == '__main__':
    main()
