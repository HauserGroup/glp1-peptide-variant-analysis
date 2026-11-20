#!/usr/bin/env python3
"""Create the structural layer used in Figure 2b.
Align AF2 models in PyMOL and color mutated side chains by binding response.
"""

import argparse
import csv
from pathlib import Path

from matplotlib import colormaps
from pymol import cmd, finish_launching


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
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "output/figures/fig2b_variant_structure.png",
    )
    return parser.parse_args()


def load_binding_values(path):
    binding = {}
    with path.open() as handle:
        for row in csv.DictReader(handle):
            if row["assay"] == "Binding":
                binding[row["Compound"]] = float(row["pct"])
    return binding


def binding_color(percent):
    """Map 0% binding to red and 100% to blue, matching Fig. 2b."""
    value = max(0.0, min(1.0, percent / 100))
    return colormaps["coolwarm_r"](value)[:3]


def main():
    args = parse_args()
    models_dir = (
        args.zenodo_root
        / "in-silico/02-predictor-benchmark/af2-relaxed-selected-structures"
    )
    binding_path = (
        args.data_root
        / "source-data/fig2a_all_assays_emm.csv"
    )
    binding = load_binding_values(binding_path)

    finish_launching(["pymol", "-cq"])
    cmd.reinitialize()
    cmd.bg_color("white")

    cmd.load(str(models_dir / "GLP1_7-37.pdb"), "wt")
    cmd.hide("everything", "all")
    cmd.show("cartoon", "wt and chain B")
    cmd.cartoon("tube", "wt and chain B")
    cmd.color("gray70", "wt and chain B")

    for model in sorted(models_dir.glob("*.pdb")):
        variant = model.stem
        if variant == "GLP1_7-37":
            continue
        if variant not in binding:
            raise ValueError(f"No binding value for {variant}")

        cmd.load(str(model), variant)
        cmd.align(f"{variant} and chain A", "wt and chain A")
        cmd.hide("everything", variant)

        residue = int(variant[1:-1]) - 6
        selection = f"{variant} and chain B and resi {residue}"
        if cmd.count_atoms(selection) == 0:
            raise ValueError(f"No mutated side-chain atoms for {variant} at residue {residue}")
        cmd.show("sticks", selection)

        color_name = f"bind_{variant}"
        cmd.set_color(color_name, binding_color(binding[variant]))
        cmd.color(color_name, selection)

    cmd.remove("hydro")
    cmd.orient("wt and chain B")
    cmd.turn("z", 90)
    cmd.zoom("wt and chain B", buffer=2)
    cmd.set("ray_opaque_background", 1)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    cmd.png(str(args.output), width=1000, height=1500, ray=1)


if __name__ == "__main__":
    main()
