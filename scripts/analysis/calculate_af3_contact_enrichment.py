from pathlib import Path
import re

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D

import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import MaxNLocator
from scipy.stats import norm, rankdata
from statsmodels.stats.multitest import multipletests

plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8
plt.rcParams['svg.fonttype'] = 'none'


def make_presence_matrix(
    df,
    complex_col='pair',
    activity_col='relative_activity_percentile',
    source_col='source',
    flatten_interaction_types=False,
    specify_peptide_residues=True,
):
    required = {
        complex_col,
        activity_col,
        source_col,
        'receptor_residue_to_use',
        'Interaction_Type',
    }

    if specify_peptide_residues:
        required.add('peptide_residue_anonymized')

    missing = required.difference(df.columns)
    if missing:
        raise KeyError(f'Missing columns: {sorted(missing)}')

    df = df.copy()

    df['contact'] = df['receptor_residue_to_use'].astype(str)

    if flatten_interaction_types:
        df['contact'] += '|contact'
    else:
        df['contact'] += '|' + df['Interaction_Type'].astype(str)

    if specify_peptide_residues:
        df['contact'] += '|' + df['peptide_residue_anonymized'].astype(str)

    metadata_check = df.groupby(complex_col)[[activity_col, source_col]].nunique(dropna=False)
    if (metadata_check > 1).any().any():
        bad = metadata_check[(metadata_check > 1).any(axis=1)].index.tolist()
        raise ValueError(f'Complexes with inconsistent experimental metadata: {bad[:10]}')

    metadata = df.groupby(complex_col, as_index=False).agg(
        source=(source_col, 'first'),
        activity=(activity_col, 'first'),
    )

    presence = (
        df[[complex_col, 'contact']]
        .drop_duplicates()
        .assign(present=1)
        .pivot(index=complex_col, columns='contact', values='present')
        .fillna(0)
        .astype('uint8')
        .reset_index()
    )

    annotations = df.groupby('contact', as_index=False).agg(
        gpcr_residue=('gpcr_residue', 'first'),
        receptor_residue_to_use=('receptor_residue_to_use', 'first'),
        peptide_residue=('peptide_residue_anonymized', 'first'),
        interaction_types_observed=(
            'Interaction_Type',
            lambda values: '; '.join(sorted(set(values.dropna().astype(str)))),
        ),
    )

    if flatten_interaction_types:
        annotations['interaction_type'] = 'Any contact'
    else:
        annotations['interaction_type'] = annotations['interaction_types_observed']

    return metadata.merge(
        presence,
        on=complex_col,
        validate='one_to_one',
    ), annotations


def van_elteren_test(values, present, strata):
    numerator = 0.0
    variance = 0.0
    effect_sum = 0.0
    effect_weight_sum = 0.0
    source_effects = []

    for source in pd.unique(strata):
        mask = strata == source
        y = np.asarray(values[mask], dtype=float)
        g = np.asarray(present[mask], dtype=bool)
        n_present = int(g.sum())
        n_absent = int((~g).sum())
        n = n_present + n_absent

        if n_present == 0 or n_absent == 0:
            continue

        ranks = rankdata(y, method='average')
        rank_sum = ranks[g].sum()
        expected = n_present * (n + 1) / 2
        rank_ss = np.square(ranks - ranks.mean()).sum()
        rank_sum_var = n_present * n_absent * rank_ss / (n * (n - 1))

        if rank_sum_var <= 0:
            continue

        weight = 1 / (n + 1)
        numerator += weight * (rank_sum - expected)
        variance += weight ** 2 * rank_sum_var

        u = rank_sum - n_present * (n_present + 1) / 2
        rank_biserial = 2 * u / (n_present * n_absent) - 1
        effect_weight = weight * n_present * n_absent
        effect_sum += effect_weight * rank_biserial
        effect_weight_sum += effect_weight
        source_effects.append(f'{source}:{rank_biserial:.3f}({n_present}/{n_absent})')

    if variance == 0 or effect_weight_sum == 0:
        return None

    z = numerator / np.sqrt(variance)
    return {
        'z_value': z,
        'p_value': 2 * norm.sf(abs(z)),
        'effect': effect_sum / effect_weight_sum,
        'n_sources': len(source_effects),
        'source_effects': '; '.join(source_effects),
    }


def run_enrichment(data, min_n=5):
    metadata_cols = {'pair', 'source', 'activity'}
    contact_cols = [col for col in data.columns if col not in metadata_cols]
    results = []

    for contact in contact_cols:
        present = data[contact].to_numpy(dtype=bool)
        n_present = int(present.sum())
        n_absent = int((~present).sum())

        if n_present < min_n or n_absent < min_n:
            continue

        test = van_elteren_test(
            values=data['activity'].to_numpy(),
            present=present,
            strata=data['source'].to_numpy(),
        )
        if test is None:
            continue

        present_values = data.loc[present, 'activity']
        absent_values = data.loc[~present, 'activity']
        results.append({
            'contact': contact,
            'n_present': n_present,
            'n_absent': n_absent,
            'frequency': n_present / (n_present + n_absent),
            'median_present': present_values.median(),
            'median_absent': absent_values.median(),
            'median_difference': present_values.median() - absent_values.median(),
            **test,
        })

    results = pd.DataFrame(results)
    if results.empty:
        return results

    results['p_fdr_bh'] = multipletests(results['p_value'], method='fdr_bh')[1]
    results['p_bonferroni'] = multipletests(results['p_value'], method='bonferroni')[1]
    results['p_corrected'] = results['p_fdr_bh']
    return results.sort_values('p_corrected').reset_index(drop=True)



def safe_filename(value):
    return re.sub(r'[^A-Za-z0-9._-]+', '_', str(value)).strip('_')



def plot_enriched_contacts_boxplots(
    data,
    selected,
    output_dir,
    ):
    plot_dir = Path(output_dir) / 'significantly-enriched-contacts'
    plot_dir.mkdir(parents=True, exist_ok=True)

    for _, result in selected.iterrows():
        contact = result['contact']
        if contact not in data.columns:
            print(f'Skipping {contact}: contact column not found in presence matrix.')
            continue

        plot_df = data[['pair', 'source', 'activity', contact]].copy()
        plot_df['contact_status'] = np.where(plot_df[contact].astype(bool), 'Present', 'Absent')
        plot_df = plot_df.dropna(subset=['activity', 'source'])

        counts = plot_df['contact_status'].value_counts()
        if not {'Present', 'Absent'}.issubset(counts.index):
            print(f'Skipping {contact}: both present and absent groups are required.')
            continue

        order = ['Present', 'Absent']
        fig, ax = plt.subplots(figsize=(2.4, 2.6))

        sns.boxplot(
            data=plot_df,
            x='contact_status',
            y='activity',
            order=order,
            width=0.55,
            showfliers=False,
            color='white',
            linewidth=0.8,
            boxprops={'facecolor': 'none'},
            ax=ax,
        )

        sns.stripplot(
            data=plot_df,
            x='contact_status',
            y='activity',
            order=order,
            color='black',
            jitter=0.22,
            size=1.5,
            alpha=0.5,
            linewidth=0,
            rasterized=True,
            ax=ax,
        )

        ax.set_xlabel('')
        ax.set_ylabel('Within-source activity percentile')
        ax.set_xticks(
            [0, 1],
            labels=[
                f'Present\nn={int(counts["Present"])}',
                f'Absent\nn={int(counts["Absent"])}',
            ],
        )
        ax.set_ylim(-0.03, 1.03)
        ax.set_title(
            f'{contact}\n'
            f'Δ median activity percentile={result["median_difference"]:+.3f}\n'
            f'rank-biserial effect={result["effect"]:+.3f}, '
            f'adj. P={result["p_corrected"]:.2g}',
            fontsize=8,
            pad=5,
        )
        sns.despine(ax=ax)
        fig.tight_layout()

        filename = safe_filename(contact)
        fig.savefig(plot_dir / f'{filename}.png', dpi=300, bbox_inches='tight', facecolor='white')
        fig.savefig(plot_dir / f'{filename}.svg', dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)


def plot_significantly_enriched_contacts(
    data,
    results,
    output_dir,
    p_thresh=0.01,
    contacts_to_plot=None,
):
    required_data = {'pair', 'source', 'activity'}
    missing_data = required_data.difference(data.columns)
    if missing_data:
        raise KeyError(f'Missing data columns: {sorted(missing_data)}')

    required_results = {'contact', 'effect', 'p_corrected'}
    missing_results = required_results.difference(results.columns)
    if missing_results:
        raise KeyError(f'Missing result columns: {sorted(missing_results)}')

    selected = results.loc[results['p_corrected'] < p_thresh].copy()

    if contacts_to_plot is not None:
        if isinstance(contacts_to_plot, str):
            contacts_to_plot = [contacts_to_plot]
        selected = selected[selected['contact'].isin(contacts_to_plot)]

    if selected.empty:
        print(f'No contacts available for boxplots at adjusted P < {p_thresh}.')
        return

    plot_enriched_contacts_boxplots(
        data,
        selected,
        output_dir,
    )
    

def run_analysis(
    contacts,
    output_dir,
    flatten_interaction_types=False,
    specify_peptide_residues=True,
    significance_threshold=0.01,
    minimum_effect_size=0.15,
    minimum_sources=2,
):
    analysis_name = (
        'flattened_interactions'
        if flatten_interaction_types
        else 'by_interaction_type'
    )

    analysis_dir = output_dir / analysis_name
    analysis_dir.mkdir(parents=True, exist_ok=True)

    data, annotations = make_presence_matrix(
        contacts,
        flatten_interaction_types=flatten_interaction_types,
        specify_peptide_residues=specify_peptide_residues,
    )

    results = run_enrichment(data, min_n=5)
    results = results.merge(
        annotations,
        on='contact',
        how='left',
        validate='one_to_one',
    )

    results['passes_p_threshold'] = (
        results['p_corrected'] < significance_threshold
    )
    results['passes_effect_threshold'] = (
        results['effect'].abs() >= minimum_effect_size
    )
    results['passes_source_threshold'] = (
        results['n_sources'] >= minimum_sources
    )
    results['passes_primary_thresholds'] = (
        results['passes_p_threshold']
        & results['passes_effect_threshold']
        & results['passes_source_threshold']
    )
    # just p-val
    significant = results.loc[results['passes_p_threshold']].copy()
    # selected
    selected = results.loc[
        results['passes_primary_thresholds']
    ].copy()

    results.to_csv(
        analysis_dir / 'contact_enrichment_results.csv',
        index=False,
    )

    significant.to_csv(
        analysis_dir / 'contact_enrichment_results.significant.csv',
        index=False,
    )

    selected.to_csv(
        analysis_dir / 'contact_enrichment_results.significant.selected.csv',
        index=False,
    )

    plot_enriched_contacts_boxplots(
        data,
        selected,
        analysis_dir,
    )

    plot_freq_vs_pvalue_ncomms_wide(
        results,
        analysis_dir / 'contact_enrichment_effect.svg',
        hue_col='effect',
        p_thresh=significance_threshold,
        add_labels=False,
    )

    # plot labelled too
    plot_freq_vs_pvalue_ncomms_wide(
        results,
        analysis_dir / 'contact_enrichment_effect_labeled.svg',
        hue_col='effect',
        p_thresh=significance_threshold,
        add_labels=True,
    )

    if not flatten_interaction_types:
        plot_freq_vs_pvalue_ncomms_wide(
            results,
            analysis_dir / 'contact_enrichment_interaction_type.svg',
            hue_col='interaction_type',
            effect_col='effect',
            categorical_palette='tab10',
            p_thresh=significance_threshold,
            add_labels=False,
        )

    print(f'{analysis_name}: tested {len(results):,} contacts across {len(data):,} peptide complexes; {len(selected):,} passed adjusted P < {significance_threshold}, |effect| >= {minimum_effect_size}, and {len(selected):,} passed n_sources >= {minimum_sources}.')

    return results


def plot_freq_vs_pvalue_ncomms_wide(
    results,
    save_path=None,
    type='both',
    hue_col='effect',
    effect_col='effect',
    p_col='p_corrected',
    freq_col='frequency',
    p_thresh=0.05,
    figsize=(2.5, 1.75),
    cmap='coolwarm_r',
    categorical_palette='tab10',
    s=12,
    alpha=0.9,
    show_colorbar=True,
    add_labels=True,
):
    required_cols = {hue_col, effect_col, p_col, freq_col}
    missing_cols = required_cols - set(results.columns)
    if missing_cols:
        raise KeyError(f'Missing columns: {sorted(missing_cols)}')

    plot_df = results[list(required_cols) + ['contact']].replace([np.inf, -np.inf], np.nan).dropna().copy()

    if type == 'negative-only':
        plot_df = plot_df[plot_df[effect_col] < 0]
    elif type == 'positive-only':
        plot_df = plot_df[plot_df[effect_col] > 0]
    elif type != 'both':
        raise ValueError("type must be 'both', 'negative-only', or 'positive-only'")

    if plot_df.empty:
        raise ValueError('No enrichment results available to plot.')

    plot_df['minus_log10_p'] = -np.log10(
        plot_df[p_col].clip(lower=np.nextafter(0, 1))
    )

    numeric_hue = pd.api.types.is_numeric_dtype(plot_df[hue_col])

    fig, ax = plt.subplots(figsize=figsize)

    if numeric_hue:
        hue_abs_max = max(float(plot_df[hue_col].abs().max()), 0.1)
        hue_abs_max = np.ceil(hue_abs_max * 10) / 10
        norm = TwoSlopeNorm(vmin=-hue_abs_max, vcenter=0, vmax=hue_abs_max)

        scatter = ax.scatter(
            plot_df[freq_col],
            plot_df['minus_log10_p'],
            c=plot_df[hue_col],
            cmap=cmap,
            norm=norm,
            s=s,
            alpha=alpha,
            edgecolors='black',
            linewidths=0.2,
            rasterized=True,
        )

        if show_colorbar:
            cax = ax.inset_axes([0.70, 0.79, 0.26, 0.055])
            colorbar = fig.colorbar(scatter, cax=cax, orientation='horizontal')
            colorbar.set_ticks([-hue_abs_max, 0, hue_abs_max])
            colorbar.ax.set_title('Rank-biserial effect', fontsize=7, pad=1)
            colorbar.ax.tick_params(length=1.5, width=0.5, pad=1, labelsize=6.5)
            colorbar.outline.set_linewidth(0.5)

    else:
        categories = sorted(plot_df[hue_col].astype(str).unique())
        colors = sns.color_palette(categorical_palette, n_colors=len(categories))
        color_map = dict(zip(categories, colors))

        for category in categories:
            subset = plot_df[plot_df[hue_col].astype(str) == category]

            ax.scatter(
                subset[freq_col],
                subset['minus_log10_p'],
                color=color_map[category],
                label=category,
                s=s,
                alpha=alpha,
                edgecolors='black',
                linewidths=0.2,
                rasterized=True,
            )

        if show_colorbar:
            handles = [
                Line2D(
                    [0],
                    [0],
                    marker='o',
                    linestyle='',
                    markerfacecolor=color_map[category],
                    markeredgecolor='black',
                    markeredgewidth=0.2,
                    markersize=4,
                    label=category,
                )
                for category in categories
            ]

            ax.legend(
                handles=handles,
                title=hue_col.replace('_', ' ').title(),
                frameon=False,
                fontsize=6.5,
                title_fontsize=7,
                loc='upper right',
                bbox_to_anchor=(1, 1),
                handletextpad=0.4,
                labelspacing=0.3,
            )

    if add_labels:
        top_3_pos_eff = plot_df.loc[plot_df[effect_col] > 0].nsmallest(3, p_col)
        top_3_neg_eff = plot_df.loc[plot_df[effect_col] < 0].nsmallest(3, p_col)
        labeled = pd.concat([top_3_pos_eff, top_3_neg_eff], axis=0)
        print(labeled)
        
        for _, row in labeled.iterrows():
            ax.annotate(
                row['contact'],
                xy=(row[freq_col], row['minus_log10_p']),
                xytext=(5, 5),
                textcoords='offset points',
                fontsize=5,
                alpha=0.7,
            )

    threshold_y = -np.log10(p_thresh)
    ax.axhline(
        threshold_y,
        color='black',
        linestyle=(0, (3, 2)),
        linewidth=0.7,
        alpha=0.8,
    )

    ax.set_xlim(0, 1)
    ax.set_ylim(
        0,
        max(1, plot_df['minus_log10_p'].max(), threshold_y) * 1.08,
    )

    ax.set_xlabel('Contact frequency', labelpad=2)
    ax.set_ylabel(r'$-\log_{10}(P_{\mathrm{adj}})$', labelpad=2)
    ax.set_xticks([0, 0.5, 1])
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.tick_params(axis='both', length=2.5, width=0.6, pad=1.5)

    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(0.6)
    ax.set_axisbelow(True)
    ax.grid(axis='y', linestyle=':', linewidth=0.5, alpha=0.45)

    fig.subplots_adjust(left=0.13, right=0.98, bottom=0.25, top=0.95)

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, facecolor='white')
        fig.savefig(save_path.with_suffix('.svg'), dpi=300, facecolor='white')

    plt.close(fig)
    


def main():
    repo_root = Path(__file__).resolve().parents[2]
    input_path = repo_root / 'zenodo-data' / 'in-silico' / '04-af3-contact-enrichment' / 'tables' / 'contacts-with-activities.csv'
    output_dir = repo_root / 'output' / 'analysis' / 'af3-contact-enrichment'
    if not input_path.is_file():
        raise FileNotFoundError(f'Missing Zenodo contact table: {input_path}')

    contacts = pd.read_csv(input_path)

    run_analysis(
        contacts,
        output_dir,
        flatten_interaction_types=False,
        specify_peptide_residues=True,
        significance_threshold=0.01,
        minimum_effect_size=0.15,
        minimum_sources=2,
    )

    run_analysis(
        contacts,
        output_dir,
        flatten_interaction_types=True,
        specify_peptide_residues=True,
        significance_threshold=0.01,
        minimum_effect_size=0.15,
        minimum_sources=2,
        
    )

    
if __name__ == '__main__':
    main()
