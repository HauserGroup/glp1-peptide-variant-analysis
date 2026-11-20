#!/usr/bin/env python3
"""Create Figure 2d from peptide and receptor neighborhood effects.
Compare peptide binding shifts with loss in neighboring GLP1R Mini-Gs variants.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "data",
    )
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[2] / "output/figures/fig2d_receptor_neighborhood.svg")
    return parser.parse_args()


def main():
    args = parse_args()
    data = pd.read_csv(
        args.data_root
        / "receptor-neighborhood/"
        "peptide_receptor_neighborhoods.csv"
    )
    data = data[(data.cutoff_A == 6) & data.peptide_position.notna()]
    data = data.drop_duplicates("peptide_position")

    x = data.median_MiniGs_loss.to_numpy(float)
    y = data.binding_effect_delta_vs_WT.to_numpy(float)
    valid = np.isfinite(x) & np.isfinite(y)
    rho, _ = spearmanr(x[valid], y[valid])

    fig, ax = plt.subplots(figsize=(4.2, 3.5))

    for _, row in data.iterrows():
        if pd.isna(row.median_MiniGs_loss) or pd.isna(row.binding_effect_delta_vs_WT):
            continue

        color = "#168d83" if row.binding_effect_delta_vs_WT >= 0 else "#c03952"
        ax.scatter(
            row.median_MiniGs_loss,
            row.binding_effect_delta_vs_WT,
            c=color,
            s=24,
        )
        ax.annotate(
            row.Compound,
            (row.median_MiniGs_loss, row.binding_effect_delta_vs_WT),
            xytext=(3, 3),
            textcoords="offset points",
            fontsize=6,
        )

    ax.axhline(0, c="#999999", lw=0.8)
    ax.text(
        0.98,
        0.98,
        f"Spearman ρ = {rho:.2f}\nn = {valid.sum()}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
    )
    ax.set(
        xlabel="Median loss in neighboring GLP1R variants",
        ylabel="Peptide binding change vs WT (logit)",
    )
    ax.spines[["top", "right"]].set_visible(False)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
