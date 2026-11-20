#!/usr/bin/env python3
"""Create Figure 3c from β-arrestin efficacy measurements.
Plot selected variant Emax values and SEM from Supplementary Table 7.
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
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[2] / "output/figures/fig3c_arrestin_efficacy.svg")
    return parser.parse_args()


def main():
    args = parse_args()
    data = pd.read_csv(
        args.data_root
        / "source-data/fig3_pharmacology_from_supplementary_table_7.csv"
    )
    data["arrestin_Emax"] = pd.to_numeric(data.arrestin_Emax, errors="coerce")
    data["arrestin_Emax_SEM"] = pd.to_numeric(
        data.arrestin_Emax_SEM,
        errors="coerce",
    )

    variants = ["Q23P", "S14T", "A24G", "F28L"]
    data = data.set_index("variant").loc[variants].reset_index()

    fig, ax = plt.subplots(figsize=(4, 2.7))
    x = range(len(data))

    ax.bar(x, data.arrestin_Emax, color="#d47c72")
    ax.errorbar(
        x,
        data.arrestin_Emax,
        yerr=data.arrestin_Emax_SEM,
        fmt="none",
        c="black",
        capsize=2,
        lw=0.8,
    )
    ax.axhline(100, ls="--", c="black", lw=0.8)
    ax.set(
        xticks=x,
        xticklabels=data.variant,
        ylabel="β-arrestin Emax (% WT)",
        ylim=(0, 115),
    )
    ax.spines[["top", "right"]].set_visible(False)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
