#!/usr/bin/env python3
"""Create Figure 3e from WT-referenced potency shifts.
Plot logEC50 changes across assays and mark censored lower-bound values.
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
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[2] / "output/figures/fig3e_potency_heatmap.svg")
    return parser.parse_args()


def main():
    args = parse_args()
    sheet = openpyxl.load_workbook(
        args.data_root / "in-vivo/CRC potency.xlsx",
        read_only=True,
        data_only=True,
    ).active
    rows = list(sheet.values)

    variants = ["I29N", "I29S", "Q23P", "A24G", "S14T", "F28L", "A30V", "A30S", "E9Q", "R36Q"]
    workbook_order = ["WT", "E9Q", "S14T", "Q23P", "A24G", "F28L", "I29N", "I29S", "A30V", "A30S", "R36Q"]
    columns = {variant: i for i, variant in enumerate(workbook_order, start=2)}
    assays = [("Binding", 3), ("cAMP", 6), ("Arrestin", 9), ("Gs", 12)]

    values = np.zeros((4, len(variants)))
    censored = np.zeros_like(values, dtype=bool)

    for assay_index, (_, row_index) in enumerate(assays):
        for variant_index, variant in enumerate(variants):
            raw_value = rows[row_index][columns[variant]]
            text = str(raw_value)
            censored[assay_index, variant_index] = text.startswith(">")
            values[assay_index, variant_index] = float(text.replace(">", ""))

    fig, ax = plt.subplots(figsize=(7, 3))
    image = ax.imshow(
        values,
        cmap="RdBu_r",
        vmin=-1,
        vmax=3.2,
        aspect="auto",
    )
    ax.set(
        xticks=range(len(variants)),
        xticklabels=variants,
        yticks=range(len(assays)),
        yticklabels=[name for name, _ in assays],
    )
    ax.tick_params(axis="x", rotation=45)

    for row in range(censored.shape[0]):
        for col in range(censored.shape[1]):
            if censored[row, col]:
                ax.text(col, row, "*", ha="center", va="center", fontsize=11)

    fig.colorbar(
        image,
        ax=ax,
        shrink=0.7,
        label="ΔlogEC50 relative to WT",
    )
    fig.tight_layout()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
