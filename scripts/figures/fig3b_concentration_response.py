#!/usr/bin/env python3
"""Create Figure 3b from normalized concentration-response experiments.
Plot assay means, SEM values, and fitted pharmacology curves for each variant.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "data",
    )
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[2] / "output/figures/fig3b_concentration_response.svg")
    return parser.parse_args()


def read_response_records(rows, variants, start_column):
    records = {variant: {} for variant in variants}

    for i, row in enumerate(rows):
        if not str(row[0]).startswith("log[ligand]"):
            continue

        for response_row in rows[i + 1 :]:
            if response_row[0] is None:
                break
            if not isinstance(response_row[0], (int, float)):
                continue

            concentration = float(response_row[0])
            for j, variant in enumerate(variants):
                column = start_column + 3 * j
                if column >= len(response_row):
                    continue

                value = response_row[column]
                if isinstance(value, (int, float)) and np.isfinite(value):
                    records[variant].setdefault(concentration, []).append(float(value))

    return records


def main():
    args = parse_args()
    workbook = openpyxl.load_workbook(
        args.data_root / "in-vitro/Raw and normalized data CRCs compiled.xlsx",
        read_only=True,
        data_only=True,
    )
    pharmacology = pd.read_csv(
        args.data_root
        / "source-data/fig3_pharmacology_from_supplementary_table_7.csv"
    ).set_index("variant")

    variants = [
        "WT",
        "E9Q",
        "S14T",
        "Q23P",
        "A24G",
        "F28L",
        "I29N",
        "I29S",
        "A30V",
        "A30S",
        "R36Q",
    ]
    colors = dict(zip(variants, plt.get_cmap("tab20").colors[: len(variants)]))
    colors["WT"] = "black"

    panels = [
        ("β-arrestin", "Arrestin recruitment", 47),
        ("Binding", "Binding", 47),
        ("Gαs", "Gs", 69),
        ("cAMP", "cAMP", 69),
    ]
    assay_keys = {
        "Arrestin recruitment": "arrestin",
        "Binding": "binding",
        "Gs": "Gs",
        "cAMP": "cAMP",
    }

    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True, sharey=True)

    for ax, (title, sheet_name, start_column) in zip(axes.flat, panels):
        rows = list(workbook[sheet_name].values)
        records = read_response_records(rows, variants, start_column)

        for variant in variants:
            concentrations = sorted(records[variant])
            means = [np.mean(records[variant][x]) for x in concentrations]
            sem = [
                np.std(records[variant][x], ddof=1) / np.sqrt(len(records[variant][x]))
                if len(records[variant][x]) > 1
                else 0
                for x in concentrations
            ]

            if sheet_name == "Binding":
                means = [100 - value for value in means]

            ax.errorbar(
                concentrations,
                means,
                yerr=sem,
                color=colors[variant],
                marker="o",
                linestyle="none",
                ms=2,
                capsize=1,
                label=variant,
            )

            assay = assay_keys[sheet_name]
            lookup = "GLP1(7-37)" if variant == "WT" else variant
            pec50 = pd.to_numeric(
                pharmacology.loc[lookup, assay + "_pEC50"],
                errors="coerce",
            )
            emax = pd.to_numeric(
                pharmacology.loc[lookup, assay + "_Emax"],
                errors="coerce",
            )

            if np.isfinite(pec50) and np.isfinite(emax):
                x_fit = np.linspace(-11, -6, 150)
                y_fit = emax / (1 + 10 ** (-pec50 - x_fit))
                ax.plot(x_fit, y_fit, color=colors[variant], lw=1)

        ax.set(
            title=title,
            xlim=(-11.3, -5.7),
            ylim=(-10, 140),
            xlabel="log10[ligand] (M)",
            ylabel="Response (% WT)",
        )
        ax.spines[["top", "right"]].set_visible(False)

    axes[0, 1].legend(
        fontsize=7,
        ncol=3,
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
    )
    fig.tight_layout()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
