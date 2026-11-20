#!/usr/bin/env python3
"""Create the Figure 1d contact-frequency inset.
Compare contact frequencies between unobserved and observed variants.
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
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "output/figures/fig1d_contact_frequency.svg",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    source_dir = args.data_root / "figure-source-data"
    receptor_peptide_path = (
        args.data_root
        / "in-silico/01-af2-mutational-scan/tables/contact_pair_enrichment_fisher.csv"
    )

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(6, 3),
        sharex=True,
        sharey=True,
    )
    panels = [
        ("Intra-receptor", source_dir / "fig1d_intra_receptor_all_contacts.csv"),
        ("Receptor–peptide", receptor_peptide_path),
    ]

    for ax, (title, table_path) in zip(axes, panels):
        data = pd.read_csv(table_path)
        x_col = "frequency_group_a_variant_level"
        y_col = "frequency_group_b_variant_level"

        colors = ["#bdbdbd"] * len(data)
        top = data.nlargest(5, "delta_frequency").index
        bottom = data.nsmallest(5, "delta_frequency").index

        for index in top:
            colors[data.index.get_loc(index)] = "#d62728"
        for index in bottom:
            colors[data.index.get_loc(index)] = "#1f77b4"

        ax.scatter(data[x_col], data[y_col], c=colors, s=5, edgecolors="none")
        ax.plot([0, 1], [0, 1], "--", c="#999999", lw=0.8)
        ax.set(
            title=title,
            xlim=(0, 1),
            ylim=(0, 1),
            xlabel="Unobserved variants",
        )

    axes[0].set_ylabel("Observed variants")
    fig.tight_layout()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
