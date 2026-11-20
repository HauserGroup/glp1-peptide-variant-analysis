#!/usr/bin/env python3
"""Create Figure 1b from proglucagon GEMME scores.
Summarize substitutions by residue and plot the smoothed median profile.
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
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[2] / "output/figures/fig1b_gemme_profile.svg")
    return parser.parse_args()


def main():
    args = parse_args()
    data_path = (
        args.data_root
        / "figure-source-data/fig1b_original_gemme_scores.csv"
    )
    data = pd.read_csv(data_path)

    profile = data.groupby("position")["gemme_score"].median().sort_index()
    smooth = profile.rolling(7, center=True, min_periods=1).mean()

    fig, ax = plt.subplots(figsize=(8, 2.8))
    ax.scatter(
        data.position,
        data.gemme_score,
        s=4,
        c="#b8b8b8",
        alpha=0.5,
        lw=0,
    )
    ax.plot(smooth.index, smooth.values, c="black", lw=1.5)

    regions = [
        (53, 81, "#f0f0f0"),
        (98, 128, "#f8eadb"),
        (146, 179, "#e8ecf4"),
    ]
    for lower, upper, color in regions:
        ax.axvspan(lower, upper, color=color, alpha=0.5, zorder=-1)

    ax.set(xlim=(0, 181), xlabel="GCG residue", ylabel="GEMME score")
    ax.spines[["top", "right"]].set_visible(False)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
