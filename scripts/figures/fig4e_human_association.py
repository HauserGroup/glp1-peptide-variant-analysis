#!/usr/bin/env python3
"""Create Figure 4e from published aggregate association results.
Plot effect estimates and confidence intervals without requiring Regenie files.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "data",
    )
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[2] / "output/figures/fig4e_human_association.svg")
    return parser.parse_args()


def main():
    args = parse_args()
    data = pd.read_csv(
        args.data_root
        / "source-data/fig4e_published_aggregate_associations.csv"
    )

    fig, ax = plt.subplots(figsize=(7.3, 4))
    n_rows = len(data)

    for i, row in data.iterrows():
        y = n_rows - 1 - i
        xerr = [[row.beta - row.ci_lower], [row.ci_upper - row.beta]]
        ax.errorbar(
            row.beta,
            y,
            xerr=xerr,
            fmt="s",
            color="black",
            capsize=0,
            ms=5,
        )
        ax.text(
            1.02,
            y,
            f"{row.beta:+.3f} ({row.ci_lower:+.3f}, {row.ci_upper:+.3f})   "
            f"{row.n_carriers}   {row.p_acat_o:.3f}",
            va="center",
            fontsize=7,
        )

    ax.axvline(0, c="black", ls="--", lw=0.8)
    ax.set(
        yticks=range(n_rows),
        yticklabels=data.phenotype.iloc[::-1],
        xlabel="Effect size (β)",
        xlim=(-0.9, 1.0),
        ylim=(-0.7, n_rows - 0.3),
    )
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
