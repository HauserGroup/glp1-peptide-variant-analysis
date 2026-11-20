#!/usr/bin/env python3
"""Create Figure 2a from single-concentration assay measurements.
Plot estimated marginal means with technical replicates for each assay.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "data",
    )
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[2] / "output/figures/fig2a_single_concentration.svg")
    return parser.parse_args()


def main():
    args = parse_args()
    source_dir = args.data_root / "source-data"
    emm = pd.read_csv(source_dir / "fig2a_all_assays_emm.csv")
    replicates = pd.read_csv(
        source_dir / "fig2a_technical_replicates_from_supplementary_table_4.csv"
    )

    order = (
        emm[emm.assay == "Binding"]
        .sort_values("pct", ascending=False)
        .Compound.tolist()
    )
    excluded = {"WT", "Buffer", "E3A", "F4L", "GLP1 1-37", "GLP1 7-36", "Exenatide"}
    order = ["WT"] + [variant for variant in order if variant not in excluded] + ["Buffer"]

    assays = [
        ("Binding", "Competition Binding", "#2e85b8"),
        ("Arrestin", "Arrestin Recruitment", "#d93535"),
        ("cAMP", "cAMP accumulation", "#888888"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(9, 7), sharey=True)
    rng = np.random.default_rng(4)

    for ax, (assay, raw_name, color) in zip(axes, assays):
        assay_emm = emm[emm.assay == assay].set_index("Compound")
        assay_replicates = replicates[replicates.assay == raw_name]

        for i, variant in enumerate(order):
            replicate_name = "GLP1 7-37" if variant == "WT" else variant
            points = (
                assay_replicates[assay_replicates.variant == replicate_name]
                .percent.dropna()
                .to_numpy()
            )
            points = np.clip(points, 0, 150)
            jitter = rng.uniform(-0.15, 0.15, len(points))

            ax.scatter(
                points,
                i + jitter,
                s=10,
                facecolors="white",
                edgecolors="#999999",
                lw=0.5,
                zorder=2,
            )

            if variant in assay_emm.index:
                row = assay_emm.loc[variant]
                xerr = [
                    [max(0, row["pct"] - row["pct.low"])],
                    [max(0, row["pct.high"] - row["pct"])],
                ]
                ax.errorbar(
                    row["pct"],
                    i,
                    xerr=xerr,
                    fmt="o",
                    ms=5,
                    color=color,
                    ecolor="black",
                    capsize=0,
                    mec="black",
                    zorder=4,
                )

        wt_response = emm[(emm.assay == assay) & (emm.Compound == "WT")].pct.iloc[0]
        ax.axvline(wt_response, ls="--", c="black", lw=0.8)
        ax.set(title=assay, xlim=(0, 150), xlabel="Response (%)")
        ax.spines[["top", "right"]].set_visible(False)

    axes[0].set_yticks(range(len(order)), order, fontsize=8)
    axes[0].invert_yaxis()
    fig.tight_layout()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
