#!/usr/bin/env python3
"""Recreate the supplied Figure 2b PyMOL scene from the Zenodo structures."""

import argparse
import sys
import re
from pathlib import Path

from pymol import cmd, finish_launching

# PyMOL `run` supplies PyMOL's own __file__; the code filename is this script.
REPO_ROOT = Path(sys._getframe().f_code.co_filename).resolve().parents[2]
MODEL_SUBDIR = Path("in-silico/02-predictor-benchmark/af2-relaxed-selected-structures")
WT_OBJECT = "GLP1_7_37"

# RGB values from the original working fig2b_variant_structure.pml scene.
VARIANT_COLORS = {
    "A24G": (199, 54, 53),
    "A30S": (102, 135, 237),
    "A30V": (242, 201, 180),
    "D15A": (207, 69, 60),
    "E9Q": (97, 128, 233),
    "F28L": (240, 138, 108),
    "G22V": (247, 177, 148),
    "G35C": (234, 212, 200),
    "G35S": (241, 204, 184),
    "G37E": (128, 163, 250),
    "I29N": (209, 73, 63),
    "I29S": (221, 95, 75),
    "K26N": (247, 175, 145),
    "Q23P": (180, 4, 38),
    "R36L": (158, 190, 255),
    "R36Q": (103, 136, 238),
    "S14T": (195, 46, 49),
    "S17R": (183, 13, 40),
    "T11I": (236, 129, 101),
    "T11N": (187, 27, 44),
    "V33L": (213, 219, 229),
    "W31C": (238, 208, 192),
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=REPO_ROOT / "zenodo-data")
    parser.add_argument(
        "--output", type=Path,
        default=REPO_ROOT / "output/figures/fig2b_variant_structure.png",
    )
    parser.add_argument("--session", type=Path, help="Optional PyMOL session (.pse) to save")
    return parser.parse_args()


def peptide_residue(variant):
    match = re.fullmatch(r"[A-Z](\d+)[A-Z]", variant)
    if match is None:
        raise ValueError(f"Invalid variant name: {variant}")
    return int(match.group(1)) - 6  # GLP-1 position 7 is chain-B residue 1.


def set_view():
    """ Orient ECD to top.
    """
    cmd.set_view(
        (-0.905282080, -0.155481130, 0.395332634,
         -0.043526996, -0.891766012, -0.450396538,
         0.422572792, -0.424943060, 0.800534189,
         0.000000000, 0.000000000, -147.697387695,
         1.924023628, -15.177440643, -2.490204811,
         -76.454101562, 371.848937988, 20.000000000)
    )


def build_scene(models_dir):
    if not models_dir.is_dir():
        raise FileNotFoundError(f"Zenodo model folder not found: {models_dir}")
    required = [models_dir / "GLP1_7-37.pdb"] + [
        models_dir / f"{variant}.pdb" for variant in VARIANT_COLORS
    ]
    missing = [path.name for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing {len(missing)} Zenodo PDBs in {models_dir}: {missing[:5]}")

    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("retain_order", 1)

    cmd.load(str(required[0]), WT_OBJECT)
    cmd.hide("everything", WT_OBJECT)
    cmd.show("cartoon", WT_OBJECT)
    cmd.color("white", f"{WT_OBJECT} and chain A")
    cmd.color("black", f"{WT_OBJECT} and chain B")

    for variant, rgb in VARIANT_COLORS.items():
        cmd.load(str(models_dir / f"{variant}.pdb"), variant)
        cmd.hide("everything", variant)
        cmd.show("cartoon", variant)
        cmd.super(
            f"{variant} and chain A and polymer.protein",
            f"{WT_OBJECT} and chain A and polymer.protein",
        )
        cmd.color("white", f"{variant} and chain A")
        cmd.color("black", f"{variant} and chain B")

        selection = f"{variant}_mutation"
        cmd.select(selection, f"{variant} and chain B and resi {peptide_residue(variant)}")
        if cmd.count_atoms(selection) == 0:
            raise ValueError(f"No mutated residue found for {variant}")
        cmd.show("sticks", selection)
        color_name = f"{variant}_color"
        cmd.set_color(color_name, list(rgb))
        cmd.color(color_name, selection)

    cmd.hide("cartoon", "chain A and resi 420-")
    cmd.hide("cartoon", "chain A and resi 1-20")
    cmd.set("cartoon_transparency", 0.5, "chain A")
    cmd.set("cartoon_transparency", 0.7, "chain B")
    cmd.set("cartoon_side_chain_helper", 1)
    cmd.set("cartoon_tube_radius", 0.4)
    cmd.set("orthoscopic", 1)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_shadows", 0)
    cmd.cartoon("loop")
    cmd.hide("everything", "elem H")
    cmd.set("cartoon_sampling", 10)
    cmd.set("stick_radius", 0.5)
    cmd.set("specular", 0)
    cmd.set("ray_trace_gain", 0.1)
    cmd.set("antialias", 2)


    cmd.disable("all")
    cmd.enable(WT_OBJECT)
    anchor = f"{WT_OBJECT} and chain B and resi 15 and name CA"
    cmd.origin(anchor)
    cmd.center(anchor)
    cmd.orient(WT_OBJECT)
    cmd.zoom(f"{WT_OBJECT} and chain B", buffer=2)

    set_view()

    # hide cartoons except for WT
    cmd.enable("all")
    for variant in VARIANT_COLORS:
        cmd.hide("cartoon", f"{variant} and chain A")
        cmd.hide("cartoon", f"{variant} and chain B")

    cmd.show("cartoon", WT_OBJECT)
    cmd.deselect()


def main():
    args = parse_args()
    finish_launching(["pymol", "-cq"])
    build_scene(args.data_root / MODEL_SUBDIR)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    cmd.png(str(args.output), width=1000, height=1500, ray=1)
    if args.session is not None:
        args.session.parent.mkdir(parents=True, exist_ok=True)
        cmd.save(str(args.session))
    print(f"Wrote {args.output}")
    if args.session is not None:
        print(f"Wrote {args.session}")


if __name__ == "__main__":
    main()
elif __name__ == "pymol":
    # `run path/to/this_script.py` builds the scene in the open PyMOL session.
    build_scene(REPO_ROOT / "zenodo-data" / MODEL_SUBDIR)
