#!/usr/bin/env python3
"""Create Figure 4d from glucose AUC and in-vitro potency data.
Regress cohort mean AUC against logEC50 for WT and seven variants.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd
from scipy.stats import linregress


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "zenodo-data",
    )
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[2] / "output/figures/fig4d_auc_potency_regression.svg")
    return parser.parse_args()


def calculate_auc(rows, times, names, group_indices):
    auc = {}

    for name in names:
        start = 1 + 12 * group_indices[name]
        data = np.array(
            [
                [
                    float(value) if isinstance(value, (int, float)) else np.nan
                    for value in row[start : start + 12]
                ]
                for row in rows[2:6]
            ]
        )
        values = np.trapezoid(data, times, axis=0)
        values = values[np.isfinite(values)]
        auc[name] = (
            np.mean(values),
            np.std(values, ddof=1) / np.sqrt(len(values)),
        )

    return auc


def main():
    args = parse_args()
    sheet = openpyxl.load_workbook(
        args.data_root / "in-vivo/in vivo data.xlsx",
        read_only=True,
        data_only=True,
    ).active
    rows = list(sheet.values)
    times = np.array([float(rows[i][0]) for i in range(2, 6)])

    names = ["WT", "E9Q", "S14T", "A24G", "F28L", "A30V", "A30S", "R36Q"]
    group_indices = {
        "WT": 1,
        "E9Q": 2,
        "S14T": 3,
        "A24G": 4,
        "F28L": 5,
        "A30V": 6,
        "A30S": 7,
        "R36Q": 8,
    }
    auc = calculate_auc(rows, times, names, group_indices)

    pharmacology = pd.read_csv(
        args.data_root
        / "figure-source-data/fig3_pharmacology_from_supplementary_table_7.csv"
    ).set_index("variant")
    assays = [
        ("Binding", "binding", "#ec8e52"),
        ("β-arrestin", "arrestin", "#e33e45"),
        ("Gαs", "Gs", "#4665aa"),
        ("cAMP", "cAMP", "#168a4a"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(6.5, 6))

    for ax, (title, assay, color) in zip(axes.flat, assays):
        data = []
        for name in names:
            lookup = "GLP1(7-37)" if name == "WT" else name
            log_ec50 = -float(pharmacology.loc[lookup, assay + "_pEC50"])
            data.append((name, auc[name][0], auc[name][1], log_ec50))

        x = np.array([record[1] for record in data])
        y = np.array([record[3] for record in data])
        fit = linregress(x, y)

        for name, x_value, x_sem, y_value in data:
            point_color = "black" if name == "WT" else color
            ax.errorbar(
                x_value,
                y_value,
                xerr=x_sem,
                fmt="o",
                ms=4,
                c=point_color,
                capsize=2,
            )

        x_line = np.array([min(x), max(x)])
        y_line = fit.intercept + fit.slope * x_line
        ax.plot(x_line, y_line, c=color, lw=1.4)
        ax.text(
            0.04,
            0.95,
            f"R² = {fit.rvalue**2:.2f}\nP = {fit.pvalue:.3g}",
            transform=ax.transAxes,
            va="top",
            fontsize=8,
        )
        ax.set(
            title=title,
            xlabel="Glucose AUC (mmol·min/L)",
            ylabel="logEC50 (M)",
        )
        ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
