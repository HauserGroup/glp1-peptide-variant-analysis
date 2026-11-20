#!/usr/bin/env python3
"""Create Figure 1c for GLP-1 amino-acid substitutions.
Recreate the grouped residue grid with gnomAD, SNV, and in-vitro overlays.
Keep the plotting code compact while matching the original visual style.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle


GLP1_SEQUENCE = "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGRG"

AA_GROUPS = {
    "positive_basic": ["R", "K", "H"],
    "negative_acidic": ["D", "E"],
    "polar_uncharged": ["N", "Q", "S", "T"],
    "special": ["G", "P", "C"],
    "hydrophobic_aliphatic": ["A", "V", "I", "L", "M"],
    "aromatic": ["F", "Y", "W"],
}

GROUP_COLORS = {
    "positive_basic": "#1f77b4",
    "negative_acidic": "#d62728",
    "polar_uncharged": "#2ca02c",
    "special": "#7f7f7f",
    "hydrophobic_aliphatic": "#ff7f0e",
    "aromatic": "#9467bd",
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "zenodo-data",
    )
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[2] / "output/figures/fig1c_gnomad_variant_grid.svg")
    return parser.parse_args()


def amino_acid_colors():
    return {
        amino_acid: GROUP_COLORS[group]
        for group, residues in AA_GROUPS.items()
        for amino_acid in residues
    }


def build_visual_rows(spacer_height=0.5):
    rows = []
    heights = []
    groups = list(AA_GROUPS.values())

    for group_index, residues in enumerate(groups):
        rows.extend(residues)
        heights.extend([1.0] * len(residues))

        if group_index < len(groups) - 1:
            rows.append(None)
            heights.append(spacer_height)

    return rows, np.asarray(heights, dtype=float)


def plot_variant_grid(data, output):
    aa_colors = amino_acid_colors()
    positions = list(range(7, 38))
    x_labels = [f"{aa}\n{pos}" for aa, pos in zip(GLP1_SEQUENCE, positions)]

    visual_rows, row_heights = build_visual_rows()
    y_edges = np.concatenate(([0.0], np.cumsum(row_heights)))
    x_edges = np.arange(len(positions) + 1)

    plot_values = np.full((len(visual_rows), len(positions)), np.nan)
    row_lookup = {aa: i for i, aa in enumerate(visual_rows) if aa is not None}

    snv_cells = []
    in_vitro_cells = []

    for row in data.itertuples(index=False):
        position = int(row.position)
        alt_aa = row.alt_aa

        if alt_aa not in row_lookup or not 7 <= position <= 37:
            continue

        x_index = position - 7
        y_index = row_lookup[alt_aa]

        if bool(row.in_gnomad):
            plot_values[y_index, x_index] = float(row.gnomAD_AF)

        if bool(row.snv_accessible):
            snv_cells.append((x_index, y_index))

        if bool(row.in_vitro):
            in_vitro_cells.append((x_index, y_index))

    positive = plot_values[np.isfinite(plot_values) & (plot_values > 0)]
    vmin = positive.min() if positive.size else 1e-6
    vmax = positive.max() if positive.size else 1.0

    reds = plt.get_cmap("Reds")
    cmap = LinearSegmentedColormap.from_list(
        "truncated_Reds",
        reds(np.linspace(0.15, 0.90, 256)),
    )
    cmap.set_bad("white")

    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    mesh = ax.pcolormesh(
        x_edges,
        y_edges,
        np.ma.masked_invalid(plot_values),
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        shading="flat",
        edgecolors="none",
    )

    for x_index, y_index in snv_cells:
        y0 = y_edges[y_index]
        ax.add_patch(
            Rectangle(
                (x_index, y0),
                1,
                row_heights[y_index],
                fill=False,
                edgecolor="black",
                linewidth=0.7,
                zorder=4,
            )
        )

    for x_index, y_index in in_vitro_cells:
        y_center = y_edges[y_index] + row_heights[y_index] / 2
        ax.plot(
            x_index + 0.5,
            y_center,
            "o",
            color="black",
            markersize=3,
            markeredgewidth=0,
            zorder=5,
        )

    for y_index, aa in enumerate(visual_rows):
        if aa is None:
            continue

        y0, y1 = y_edges[y_index], y_edges[y_index + 1]
        ax.vlines(
            np.arange(len(positions) + 1),
            y0,
            y1,
            colors="black",
            linestyles=":",
            linewidth=0.6,
            alpha=0.25,
            zorder=3,
        )
        ax.hlines(
            [y0, y1],
            0,
            len(positions),
            colors="black",
            linestyles=":",
            linewidth=0.6,
            alpha=0.25,
            zorder=3,
        )

    ax.set_xlim(0, len(positions))
    ax.set_ylim(y_edges[-1], 0)
    ax.set_aspect("equal")

    ax.set_xticks(np.arange(len(positions)) + 0.5)
    ax.set_xticklabels(x_labels, fontsize=7)

    y_positions = [
        y_edges[i] + row_heights[i] / 2
        for i, aa in enumerate(visual_rows)
        if aa is not None
    ]
    y_labels = [aa for aa in visual_rows if aa is not None]
    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=7)

    for label, wt_aa in zip(ax.get_xticklabels(), GLP1_SEQUENCE):
        label.set_color(aa_colors[wt_aa])
        label.set_fontweight("bold")

    for label, alt_aa in zip(ax.get_yticklabels(), y_labels):
        label.set_color(aa_colors[alt_aa])
        label.set_fontweight("bold")

    ax.set_xlabel("GLP-1 residue", fontsize=9)
    ax.set_ylabel("Substituted residue", fontsize=9)
    ax.tick_params(axis="both", length=0, pad=2)

    for spine in ax.spines.values():
        spine.set_visible(False)

    colorbar = fig.colorbar(
        mesh,
        ax=ax,
        shrink=0.40,
        pad=0.02,
        fraction=0.035,
    )
    colorbar.set_label("gnomAD allele frequency", fontsize=8)
    colorbar.outline.set_linewidth(0.6)
    colorbar.ax.tick_params(length=2, width=0.6, labelsize=7)

    ticks = [value for value in np.arange(1, 8) * 1e-6 if vmin <= value <= vmax]
    if ticks:
        colorbar.set_ticks(ticks)
        colorbar.set_ticklabels([f"{value:.0e}" for value in ticks])

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    data_path = (
        args.data_root
        / "in-silico/01-af2-mutational-scan/single-nucleotide-variant-accessibility/gnomAD_v4_gcg_missense_glp1_7-37.all_variants.csv"
    )
    data = pd.read_csv(data_path)
    plot_variant_grid(data, args.output)


if __name__ == "__main__":
    main()
