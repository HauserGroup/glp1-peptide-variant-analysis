#!/usr/bin/env python3
"""Create Figure 4b from mouse blood-glucose measurements.
Plot cohort means and SEM values across the measured time course.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import openpyxl


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "zenodo-data",
    )
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[2] / "output/figures/fig4b_glucose_timecourse.svg")
    return parser.parse_args()


def main():
    args = parse_args()
    sheet = openpyxl.load_workbook(
        args.data_root / "in-vivo/in vivo data.xlsx",
        read_only=True,
        data_only=True,
    ).active
    rows = list(sheet.values)

    times = np.array([float(rows[i][0]) for i in range(2, 6)])
    names = ["Vehicle", "WT", "E9Q", "S14T", "A24G", "F28L", "A30V", "A30S", "R36Q"]
    colors = ["#666666", "#000000", "#ce84b7", "#b7d978", "#2585aa", "#7350a3", "#e77d85", "#87a876", "#234c90"]

    fig, ax = plt.subplots(figsize=(5, 3.3))

    for group_index, (name, color) in enumerate(zip(names, colors)):
        start = 1 + 12 * group_index
        data = np.array(
            [
                [
                    float(value) if isinstance(value, (int, float)) else np.nan
                    for value in row[start : start + 12]
                ]
                for row in rows[2:6]
            ]
        )
        mean = np.nanmean(data, axis=1)
        n = np.sum(np.isfinite(data), axis=1)
        sem = np.nanstd(data, axis=1, ddof=1) / np.sqrt(n)

        ax.errorbar(
            times,
            mean,
            yerr=sem,
            label=name,
            color=color,
            marker="o",
            ms=3,
            lw=1.3,
            capsize=2,
        )

    ax.set(
        xlabel="Time (min)",
        ylabel="Blood glucose (mmol/L)",
        xlim=(-2, 62),
        ylim=(0, 26),
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(
        fontsize=7,
        ncol=2,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
