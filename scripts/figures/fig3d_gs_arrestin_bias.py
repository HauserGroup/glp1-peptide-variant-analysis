#!/usr/bin/env python3
"""Create Figure 3d from Gαs-versus-arrestin bias values.
Plot ΔΔlog(Emax/EC50) estimates and SEM from Supplementary Table 8.
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
        default=Path(__file__).resolve().parents[2] / "zenodo-data",
    )
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[2] / "output/figures/fig3d_gs_arrestin_bias.svg")
    return parser.parse_args()


def main():
    args = parse_args()
    data = pd.read_csv(
        args.data_root
        / "figure-source-data/fig3_bias_from_supplementary_table_8.csv"
    )
    data["bias"] = pd.to_numeric(data.bias, errors="coerce")
    data["bias_sem"] = pd.to_numeric(data.bias_sem, errors="coerce")
    data = data[data.bias.notna()]

    variants = ["A30S", "E9Q", "A30V", "S14T", "R36Q", "F28L", "Q23P", "A24G"]
    data = data.set_index("variant").loc[variants].reset_index()

    fig, ax = plt.subplots(figsize=(4, 3))
    x = range(len(data))

    ax.bar(x, data.bias, color="#678fb9")
    ax.errorbar(
        x,
        data.bias,
        yerr=data.bias_sem,
        fmt="none",
        c="black",
        capsize=2,
        lw=0.8,
    )
    ax.axhline(0, c="black", lw=0.8)
    ax.set(
        xticks=x,
        xticklabels=data.variant,
        ylabel="Gαs − β-arrestin bias",
    )
    ax.tick_params(axis="x", rotation=45)
    ax.spines[["top", "right"]].set_visible(False)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
