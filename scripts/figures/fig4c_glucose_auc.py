#!/usr/bin/env python3
"""Create Figure 4c from individual mouse glucose AUC values.
Plot group means, SEM, individual animals, and stored Dunnett comparisons.
"""

import argparse
import json
from pathlib import Path
from zipfile import ZipFile

import matplotlib.pyplot as plt
import numpy as np
import openpyxl


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "data",
    )
    parser.add_argument(
        "--zenodo-root",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "zenodo-data",
    )
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[2] / "output/figures/fig4c_glucose_auc.svg")
    return parser.parse_args()


def read_dunnett_pvalues(prism_path):
    with ZipFile(prism_path) as archive:
        for path in archive.namelist():
            if not path.endswith("/sheet.json") or not path.startswith("analyses/"):
                continue

            metadata = json.loads(archive.read(path))
            if metadata.get("title") != "Ordinary one-way ANOVA of AUC":
                continue

            results_path = path.replace("sheet.json", "results.json")
            results = json.loads(archive.read(results_path))
            tests = results["content"]["ordinaryAnovaModel"]["multipleComparisons"]["tests"]
            break
        else:
            raise ValueError("Dunnett AUC analysis not found in Prism file")

    comparisons = next(iter(tests.values()))["value"]
    return {
        record["id"]["col"]["title"]: record["value"]["pValueAdjusted"]["value"]
        for record in comparisons.values()
    }


def main():
    args = parse_args()
    sheet = openpyxl.load_workbook(
        args.zenodo_root / "in-vivo/in vivo data.xlsx",
        read_only=True,
        data_only=True,
    ).active
    rows = list(sheet.values)

    times = np.array([float(rows[i][0]) for i in range(2, 6)])
    names = ["Vehicle", "S14T", "F28L", "A24G", "WT", "R36Q", "A30V", "E9Q", "A30S"]
    offsets = {
        "Vehicle": 0,
        "WT": 1,
        "E9Q": 2,
        "S14T": 3,
        "A24G": 4,
        "F28L": 5,
        "A30V": 6,
        "A30S": 7,
        "R36Q": 8,
    }
    colors = ["#aaaaaa", "#bedf8a", "#7051a4", "#268aad", "#111111", "#244e8f", "#e77f86", "#cd91bd", "#88aa78"]

    fig, ax = plt.subplots(figsize=(6, 3.4))
    rng = np.random.default_rng(13)

    for i, (name, color) in enumerate(zip(names, colors)):
        start = 1 + 12 * offsets[name]
        data = np.array(
            [
                [
                    float(value) if isinstance(value, (int, float)) else np.nan
                    for value in row[start : start + 12]
                ]
                for row in rows[2:6]
            ]
        )
        auc = np.trapezoid(data, times, axis=0)
        auc = auc[np.isfinite(auc)]
        mean = np.mean(auc)
        sem = np.std(auc, ddof=1) / np.sqrt(len(auc))

        ax.bar(i, mean, color=color, edgecolor="black", width=0.75)
        ax.scatter(
            i + rng.uniform(-0.17, 0.17, len(auc)),
            auc,
            c="black",
            s=8,
            zorder=3,
        )
        ax.errorbar(i, mean, yerr=sem, fmt="none", c="black", capsize=2)

    pvalues = read_dunnett_pvalues(args.zenodo_root / "in-vivo/In vivo data.prism")
    for height, variant in [(1100, "S14T"), (960, "F28L"), (840, "A30S")]:
        pvalue = pvalues[variant]
        i = names.index(variant)
        j = names.index("WT")
        ax.plot(
            [i, i, j, j],
            [height - 25, height, height, height - 25],
            color="black",
            lw=0.8,
        )
        ax.text(
            (i + j) / 2,
            height + 10,
            "**" if pvalue < 0.01 else "*",
            ha="center",
            fontsize=9,
        )

    ax.set(
        xticks=range(len(names)),
        xticklabels=names,
        ylabel="Glucose AUC (mmol·min/L)",
        ylim=(0, 1220),
    )
    ax.tick_params(axis="x", rotation=45)
    ax.spines[["top", "right"]].set_visible(False)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
