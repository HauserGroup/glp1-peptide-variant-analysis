#!/usr/bin/env python3
"""Create the plotted portions of Figure 2e.
Show AF3 contact enrichment and summarize effects by peptide position.
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
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[2] / "output/figures/fig2e_af3_contact_enrichment.svg")
    return parser.parse_args()


def main():
    args = parse_args()
    table_dir = args.data_root / "in-silico/04-af3-contact-enrichment/tables"
    contacts = pd.read_csv(table_dir / "contact-enrichment-all-tested.csv")
    summary = pd.read_csv(table_dir / "peptide-position-effect-summary.csv")

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(10, 3),
        gridspec_kw={"width_ratios": [1, 2]},
    )

    adjusted_p = contacts.p_corrected.clip(lower=np.finfo(float).tiny)
    colors = np.where(contacts.effect >= 0, "#4977aa", "#ed8952")

    ax = axes[0]
    ax.scatter(
        contacts.frequency,
        -np.log10(adjusted_p),
        c=colors,
        s=12,
        alpha=0.6,
        edgecolors="black",
        linewidths=0.2,
    )
    ax.axhline(-np.log10(0.05), ls="--", c="black", lw=0.7)
    ax.set(xlabel="Contact frequency", ylabel="−log10(adjusted P)")

    for _, row in contacts.nsmallest(5, "p_corrected").iterrows():
        label = str(row.contact).replace("|contact|", "–")
        p_value = max(row.p_corrected, np.finfo(float).tiny)
        ax.annotate(
            label,
            (row.frequency, -np.log10(p_value)),
            fontsize=6,
        )

    ax = axes[1]
    positions = contacts.peptide_residue.astype(str).str.extract(r"(\d+)")[0].astype(float)
    jitter = np.random.default_rng(2025).uniform(-0.17, 0.17, len(contacts))
    ax.scatter(
        positions + jitter,
        contacts.effect,
        s=8,
        c=colors,
        alpha=0.45,
        lw=0,
    )
    ax.fill_between(
        summary.peptide_position,
        summary.q1,
        summary.q3,
        color="#87ad70",
        alpha=0.3,
        lw=0,
    )
    ax.plot(
        summary.peptide_position,
        summary.median_effect,
        c="#2381bb",
        lw=1.5,
        marker="o",
        ms=2,
    )
    ax.axhline(0, ls="--", c="black", lw=0.7)
    ax.set(
        xlabel="GLP-1 peptide position",
        ylabel="Activity effect",
        xlim=(6, 38),
    )

    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
