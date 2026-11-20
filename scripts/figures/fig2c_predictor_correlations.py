#!/usr/bin/env python3
"""Create Figure 2c from predictor and experimental correlations.
Plot absolute Spearman correlations for the deposited 22-variant benchmark.
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
        default=Path(__file__).resolve().parents[2] / "zenodo-data",
    )
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[2] / "output/figures/fig2c_predictor_correlations.svg")
    return parser.parse_args()


def annotate_heatmap(ax, matrix, skip_upper_triangle=False):
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            if skip_upper_triangle and col >= row:
                continue

            value = matrix[row, col]
            text_color = "white" if value > 0.75 else "black"
            ax.text(
                col,
                row,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=7,
                color=text_color,
            )


def main():
    args = parse_args()
    data = pd.read_csv(
        args.data_root
        / "in-silico/02-predictor-benchmark/benchmark-scores-22-variants.csv"
    )

    predictors = {
        "Boltz-2": "boltz2_iptm",
        "Chai": "chai_iptm",
        "FEP ΔΔG": "BAR_ddG",
        "AF2": "AF2_iptm+ptm",
        "AF3": "af3_iptm",
        "MTNN potency": "predicted_potency_GLP-1R",
        "ProteoCast": "proteocast_score",
        "RIA ΔΔG": "dG_cross",
        "AlphaMissense": "alphamissense_pathogenicity_score",
    }
    experimental = {
        "Binding": "mean_binding",
        "Arrestin": "mean_arrestin",
        "cAMP": "mean_cAMP",
    }

    predictor_corr = np.array(
        [
            [
                abs(data[[experiment, predictor]].corr("spearman").iloc[0, 1])
                for predictor in predictors.values()
            ]
            for experiment in experimental.values()
        ]
    )
    experimental_corr = (
        data[list(experimental.values())].corr("spearman").abs().to_numpy()
    )

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(9, 2.8),
        gridspec_kw={"width_ratios": [3, 1]},
    )
    panels = [
        (axes[0], predictor_corr, list(predictors), list(experimental)),
        (axes[1], experimental_corr, list(experimental), list(experimental)),
    ]

    for ax, matrix, x_labels, y_labels in panels:
        ax.imshow(matrix, cmap="RdBu_r", vmin=0, vmax=1)
        ax.set_xticks(
            range(len(x_labels)),
            x_labels,
            rotation=45,
            ha="right",
            fontsize=7,
        )
        ax.set_yticks(range(len(y_labels)), y_labels, fontsize=8)
        ax.set_xticks(np.arange(-0.5, matrix.shape[1], 1), minor=True)
        ax.set_yticks(np.arange(-0.5, matrix.shape[0], 1), minor=True)
        ax.grid(which="minor", color="white", lw=0.5)
        ax.tick_params(which="minor", bottom=False, left=False)

    annotate_heatmap(axes[0], predictor_corr)
    annotate_heatmap(axes[1], experimental_corr, skip_upper_triangle=True)

    axes[0].set_title("Predictors")
    axes[1].set_title("Experimental metrics")
    fig.tight_layout()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
