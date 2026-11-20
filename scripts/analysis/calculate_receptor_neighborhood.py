#!/usr/bin/env python3
"""Integrate GLP-1 peptide binding variants with GLP1R variant phenotypes.

Inputs are the two supplied CSV files and PDB 6VCB in mmCIF format. The script
parses the structure without a specialist structural-biology dependency, maps
minimum heavy-atom distances between GLP-1 and GLP1R residues, performs
residue-level permutation tests, and creates processed tables and figures.

The receptor screen used exendin-4, whereas the peptide screen used GLP-1.
Accordingly, all cross-dataset analyses are framed as structural concordance,
not as paired ligand-receptor epistasis or causal variance decomposition.


Run with the repository analysis runner or this script directly.

"""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
# svg should be editable (svg fonttype none)
plt.rcParams["svg.fonttype"] = "none"

SEED = 20260820
DIRECT_CUTOFF = 4.5
PROXIMAL_CUTOFF = 6.0
SENSITIVITY_CUTOFFS = (4.0, 4.5, 5.0, 6.0, 8.0)
N_PERM = 100_000

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

COL_SURFACE = "Log Delta Surface expression"
COL_MINIGS = "Log Delta Mini-Gs"
COL_ENDO = "Log Delta GLP-1R endocytosis response"
COL_LOF_MINIGS = "LoF vs WT"
COL_SPECIFIC_MINIGS = "LoF, controlled for surface expression"
COL_LOF_ENDO = "LoF vs WT endocytosis"
COL_SPECIFIC_ENDO = "LoF, controlled for surface expression endocytosis"


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize whitespace and Greek delta while retaining readable labels."""
    out = df.copy()
    cols = []
    for col in out.columns:
        col = re.sub(r"\s+", " ", col.strip()).replace("Δ", "Delta")
        cols.append(col)
    out.columns = cols
    return out


def read_atom_site(cif_path: Path) -> pd.DataFrame:
    """Read standard one-row-per-atom PDBx/mmCIF atom_site records from 6VCB."""
    rows = []
    with cif_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(("ATOM ", "HETATM ")):
                continue
            fields = line.split()
            if len(fields) < 21:
                raise ValueError(f"Unexpected atom_site row: {line[:100]}")
            rows.append(
                {
                    "group": fields[0],
                    "element": fields[2],
                    "atom": fields[19],
                    "label_chain": fields[6],
                    "auth_resid": int(fields[16]),
                    "resname": fields[17],
                    "auth_chain": fields[18],
                    "x": float(fields[10]),
                    "y": float(fields[11]),
                    "z": float(fields[12]),
                    "model": int(fields[20]),
                }
            )
    atoms = pd.DataFrame(rows)
    if atoms.empty:
        raise ValueError("No atom_site records found in mmCIF file")
    return atoms


def residue_coordinates(atoms: pd.DataFrame, chain: str) -> dict[tuple[int, str], np.ndarray]:
    subset = atoms[
        (atoms["group"] == "ATOM")
        & (atoms["auth_chain"] == chain)
        & (atoms["model"] == 1)
        & (atoms["element"] != "H")
    ]
    return {
        (int(resid), str(resname)): group[["x", "y", "z"]].to_numpy(float)
        for (resid, resname), group in subset.groupby(["auth_resid", "resname"], sort=True)
    }


def contact_distance_table(atoms: pd.DataFrame) -> pd.DataFrame:
    receptor = residue_coordinates(atoms, "R")
    peptide = residue_coordinates(atoms, "P")
    rows = []
    for (pep_pos, pep_res3), pep_xyz in peptide.items():
        for (rec_pos, rec_res3), rec_xyz in receptor.items():
            delta = pep_xyz[:, None, :] - rec_xyz[None, :, :]
            distance = float(np.sqrt(np.sum(delta * delta, axis=2)).min())
            rows.append(
                {
                    "peptide_position": pep_pos,
                    "peptide_residue": AA3_TO_1.get(pep_res3, pep_res3),
                    "receptor_position": rec_pos,
                    "receptor_residue": AA3_TO_1.get(rec_res3, rec_res3),
                    "minimum_heavy_atom_distance_A": distance,
                }
            )
    return pd.DataFrame(rows).sort_values(
        ["peptide_position", "minimum_heavy_atom_distance_A", "receptor_position"]
    )


def parse_peptide_variants(path: Path) -> tuple[pd.DataFrame, float, float]:
    """Parse peptide data and define two WT-referenced effect outcomes.

    ``binding_effect_delta_vs_WT`` retains direction (positive, increased;
    negative, decreased). ``binding_effect_magnitude_vs_WT`` is its absolute
    magnitude and must not be described as loss because it includes gains.
    Both retain the original effect scale supplied in the peptide dataset.
    """
    peptide = pd.read_csv(path)
    required = {"Compound", "effect", "pct"}
    missing = required - set(peptide.columns)
    if missing:
        raise ValueError(f"Peptide table is missing columns: {sorted(missing)}")
    peptide["pct"] = pd.to_numeric(peptide["pct"], errors="raise")
    if not peptide["pct"].between(0, 100).all():
        raise ValueError("Peptide pct values must lie between 0 and 100")
    match = peptide["Compound"].str.extract(r"^([A-Z])(\d+)([A-Z])$")
    peptide["reference_aa"] = match[0]
    peptide["peptide_position"] = pd.to_numeric(match[1], errors="coerce")
    peptide["alternate_aa"] = match[2]
    wt = peptide.loc[peptide["Compound"] == "WT", "effect"]
    if len(wt) != 1:
        raise ValueError("Peptide table must contain exactly one WT row")
    wt_effect = float(wt.iloc[0])
    wt_pct = float(peptide.loc[peptide["Compound"] == "WT", "pct"].iloc[0])
    if wt_pct <= 0:
        raise ValueError("WT pct must be greater than zero")
    peptide["binding_effect_delta_vs_WT"] = peptide["effect"] - wt_effect
    peptide["binding_effect_magnitude_vs_WT"] = peptide[
        "binding_effect_delta_vs_WT"
    ].abs()
    peptide["binding_displacement_delta_pp_vs_WT"] = peptide["pct"] - wt_pct
    peptide["binding_displacement_retained_vs_WT"] = peptide["pct"] / wt_pct
    peptide["in_6VCB_active_peptide"] = peptide["peptide_position"].between(7, 37)
    peptide.loc[peptide["Compound"] == "WT", "in_6VCB_active_peptide"] = True
    return peptide, wt_effect, wt_pct


def parse_receptor_variants(path: Path) -> pd.DataFrame:
    receptor = normalize_columns(pd.read_csv(path))
    receptor = receptor.rename(
        columns={
            "SE.1": "SE Mini-Gs",
            "N.1": "N Mini-Gs",
            "SE.2": "SE endocytosis",
            "N.2": "N endocytosis",
            "LoF vs WT.1": COL_LOF_ENDO,
            "LoF, controlled for surface expression.1": COL_SPECIFIC_ENDO,
        }
    )
    # After whitespace normalization, the first occurrence keeps its original name.
    receptor = receptor.rename(
        columns={
            "Log Delta Surface expression": COL_SURFACE,
            "Log Delta Mini-Gs": COL_MINIGS,
            "Log Delta GLP-1R endocytosis response": COL_ENDO,
            "LoF vs WT": COL_LOF_MINIGS,
            "LoF, controlled for surface expression": COL_SPECIFIC_MINIGS,
        }
    )
    parsed = receptor["Protein consequence"].str.extract(r"^([A-Z])(\d+)([A-Z*]+)$")
    if parsed.isna().any(axis=None):
        bad = receptor.loc[parsed.isna().any(axis=1), "Protein consequence"].tolist()
        raise ValueError(f"Could not parse receptor consequences: {bad}")
    receptor["reference_aa"] = parsed[0]
    receptor["receptor_position"] = parsed[1].astype(int)
    receptor["alternate_aa"] = parsed[2]
    receptor["Mini-Gs efficiency proxy"] = receptor[COL_MINIGS] - receptor[COL_SURFACE]
    receptor["Endocytosis efficiency proxy"] = receptor[COL_ENDO] - receptor[COL_SURFACE]
    receptor["Mini-Gs loss magnitude"] = (-receptor[COL_MINIGS]).clip(lower=0)
    receptor["Endocytosis loss magnitude"] = (-receptor[COL_ENDO]).clip(lower=0)
    return receptor


def add_receptor_distances(receptor: pd.DataFrame, distances: pd.DataFrame) -> pd.DataFrame:
    minima = distances.groupby("receptor_position")["minimum_heavy_atom_distance_A"].min()
    receptor = receptor.copy()
    receptor["minimum_distance_to_GLP1_A"] = receptor["receptor_position"].map(minima)
    receptor["distance_class"] = pd.cut(
        receptor["minimum_distance_to_GLP1_A"],
        bins=[-np.inf, DIRECT_CUTOFF, PROXIMAL_CUTOFF, np.inf],
        labels=["Direct (<=4.5 A)", "Proximal (4.5-6 A)", "Distal (>6 A)"],
    ).astype(object)
    receptor.loc[receptor["minimum_distance_to_GLP1_A"].isna(), "distance_class"] = "Unresolved"
    return receptor


def aggregate_receptor_residues(receptor: pd.DataFrame) -> pd.DataFrame:
    grouped = receptor.groupby("receptor_position", as_index=False).agg(
        receptor_residue=("reference_aa", "first"),
        n_receptor_variants=("Protein consequence", "size"),
        minimum_distance_to_GLP1_A=("minimum_distance_to_GLP1_A", "first"),
        median_surface_effect=(COL_SURFACE, "median"),
        median_MiniGs_effect=(COL_MINIGS, "median"),
        median_endocytosis_effect=(COL_ENDO, "median"),
        median_abs_MiniGs_effect=(COL_MINIGS, lambda s: float(np.median(np.abs(s)))),
        median_abs_endocytosis_effect=(COL_ENDO, lambda s: float(np.median(np.abs(s)))),
        max_MiniGs_loss=("Mini-Gs loss magnitude", "max"),
        max_endocytosis_loss=("Endocytosis loss magnitude", "max"),
        fraction_MiniGs_LoF=(COL_LOF_MINIGS, "mean"),
        fraction_specific_MiniGs_LoF=(COL_SPECIFIC_MINIGS, "mean"),
        fraction_endocytosis_LoF=(COL_LOF_ENDO, "mean"),
        fraction_specific_endocytosis_LoF=(COL_SPECIFIC_ENDO, "mean"),
        median_MiniGs_loss=("Mini-Gs loss magnitude", "median"),
    )
    grouped["distance_class"] = pd.cut(
        grouped["minimum_distance_to_GLP1_A"],
        bins=[-np.inf, DIRECT_CUTOFF, PROXIMAL_CUTOFF, np.inf],
        labels=["Direct (<=4.5 A)", "Proximal (4.5-6 A)", "Distal (>6 A)"],
    ).astype(object)
    grouped.loc[grouped["minimum_distance_to_GLP1_A"].isna(), "distance_class"] = "Unresolved"
    return grouped


def permutation_mean_difference(
    values: np.ndarray, groups: np.ndarray, rng: np.random.Generator, n_perm: int = N_PERM
) -> tuple[float, float]:
    values = np.asarray(values, float)
    groups = np.asarray(groups, bool)
    keep = np.isfinite(values)
    values, groups = values[keep], groups[keep]
    observed = float(values[groups].mean() - values[~groups].mean())
    n_group = int(groups.sum())
    exceed = 0
    for _ in range(n_perm):
        selected = rng.choice(len(values), size=n_group, replace=False)
        perm_group = np.zeros(len(values), dtype=bool)
        perm_group[selected] = True
        diff = values[perm_group].mean() - values[~perm_group].mean()
        exceed += abs(diff) >= abs(observed)
    return observed, (exceed + 1) / (n_perm + 1)


def benjamini_hochberg(p_values: list[float]) -> np.ndarray:
    p = np.asarray(p_values, float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def receptor_proximity_tests(residue: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    resolved = residue.dropna(subset=["minimum_distance_to_GLP1_A"]).copy()
    rows = []
    metrics = [
        "median_abs_MiniGs_effect",
        "median_MiniGs_loss",
        "max_MiniGs_loss",
        "fraction_specific_MiniGs_LoF",
        "median_abs_endocytosis_effect",
        "max_endocytosis_loss",
        "fraction_specific_endocytosis_LoF",
    ]
    for background_name, subset in [
        ("all resolved receptor residues", resolved),
        ("local receptor residues <=12 A", resolved[resolved["minimum_distance_to_GLP1_A"] <= 12]),
    ]:
        direct = subset["minimum_distance_to_GLP1_A"] <= DIRECT_CUTOFF
        for metric in metrics:
            diff, p_value = permutation_mean_difference(subset[metric].to_numpy(), direct.to_numpy(), rng)
            rows.append(
                {
                    "background": background_name,
                    "metric": metric,
                    "n_direct_residues": int(direct.sum()),
                    "n_background_residues": int((~direct).sum()),
                    "direct_mean": float(subset.loc[direct, metric].mean()),
                    "background_mean": float(subset.loc[~direct, metric].mean()),
                    "mean_difference": diff,
                    "permutation_p": p_value,
                }
            )
    tests = pd.DataFrame(rows)
    tests["BH_adjusted_p"] = np.nan
    for background, idx in tests.groupby("background").groups.items():
        tests.loc[idx, "BH_adjusted_p"] = benjamini_hochberg(tests.loc[idx, "permutation_p"].tolist())
    return tests


def peptide_receptor_neighborhoods(
    peptide: pd.DataFrame, receptor: pd.DataFrame, distances: pd.DataFrame
) -> pd.DataFrame:
    variants = peptide[peptide["Compound"] != "WT"].copy()
    outputs = []
    for cutoff in SENSITIVITY_CUTOFFS:
        close = distances[distances["minimum_heavy_atom_distance_A"] <= cutoff]
        mapping = close.groupby("peptide_position")["receptor_position"].apply(set).to_dict()
        for row in variants.itertuples(index=False):
            position = getattr(row, "peptide_position")
            in_structure = bool(getattr(row, "in_6VCB_active_peptide"))
            rec_positions = mapping.get(int(position), set()) if in_structure and pd.notna(position) else set()
            linked = receptor[receptor["receptor_position"].isin(rec_positions)]
            outputs.append(
                {
                    "Compound": getattr(row, "Compound"),
                    "peptide_position": position,
                    "cutoff_A": cutoff,
                    "in_6VCB_active_peptide": in_structure,
                    "binding_effect_delta_vs_WT": getattr(
                        row, "binding_effect_delta_vs_WT"
                    ),
                    "binding_effect_magnitude_vs_WT": getattr(
                        row, "binding_effect_magnitude_vs_WT"
                    ),
                    "binding_displacement_delta_pp_vs_WT": getattr(
                        row, "binding_displacement_delta_pp_vs_WT"
                    ),
                    "binding_displacement_retained_vs_WT": getattr(
                        row, "binding_displacement_retained_vs_WT"
                    ),
                    "n_structural_receptor_residues": len(rec_positions),
                    "n_assayed_receptor_residues": int(linked["receptor_position"].nunique()),
                    "n_assayed_receptor_variants": int(len(linked)),
                    "median_abs_MiniGs_effect": float(np.median(np.abs(linked[COL_MINIGS]))) if len(linked) else np.nan,
                    "max_MiniGs_loss": float(linked["Mini-Gs loss magnitude"].max()) if len(linked) else np.nan,
                    "fraction_specific_MiniGs_LoF": float(linked[COL_SPECIFIC_MINIGS].mean()) if len(linked) else np.nan,
                    "median_abs_endocytosis_effect": float(np.median(np.abs(linked[COL_ENDO]))) if len(linked) else np.nan,
                    "max_endocytosis_loss": float(linked["Endocytosis loss magnitude"].max()) if len(linked) else np.nan,
                    "fraction_specific_endocytosis_LoF": float(linked[COL_SPECIFIC_ENDO].mean()) if len(linked) else np.nan,
                    "linked_receptor_positions": ";".join(map(str, sorted(rec_positions & set(receptor["receptor_position"])))),
                    "linked_receptor_variants": ";".join(sorted(linked["Protein consequence"].tolist())),
                    "median_MiniGs_loss": (float(linked["Mini-Gs loss magnitude"].median()) if len(linked) else np.nan),
                }
            )
    out = pd.DataFrame(outputs)
    out["coverage_class"] = np.select(
        [
            ~out["in_6VCB_active_peptide"],
            out["n_assayed_receptor_variants"] > 0,
            out["n_structural_receptor_residues"] > 0,
        ],
        ["Outside 6VCB peptide", "Assayed receptor neighbor", "Structural contacts, no assayed receptor variant"],
        default="No receptor neighbor at threshold",
    )
    return out


def peptide_neighborhood_correlations(neighborhoods: pd.DataFrame) -> pd.DataFrame:
    rows = []
    metrics = [
        "median_abs_MiniGs_effect",
        "median_MiniGs_loss",
        "max_MiniGs_loss",
        "fraction_specific_MiniGs_LoF",
        "median_abs_endocytosis_effect",
        "max_endocytosis_loss",
        "fraction_specific_endocytosis_LoF",
    ]
    for cutoff, subset in neighborhoods.groupby("cutoff_A"):
        for outcome in [
            "binding_effect_delta_vs_WT",
            "binding_effect_magnitude_vs_WT",
        ]:
            for metric in metrics:
                data = subset[[outcome, metric]].dropna()
                if len(data) >= 4 and data[metric].nunique() > 1 and data[outcome].nunique() > 1:
                    rho, p_value = stats.spearmanr(data[outcome], data[metric])
                else:
                    rho, p_value = np.nan, np.nan
                rows.append(
                    {
                        "cutoff_A": cutoff,
                        "outcome": outcome,
                        "receptor_neighborhood_metric": metric,
                        "n_peptide_variants": len(data),
                        "spearman_rho": rho,
                        "p_value": p_value,
                    }
                )
    result = pd.DataFrame(rows)
    result["BH_adjusted_p"] = np.nan
    for (cutoff, outcome), idx in result.groupby(["cutoff_A", "outcome"]).groups.items():
        valid = result.loc[idx, "p_value"].notna()
        valid_idx = result.loc[idx].index[valid]
        if len(valid_idx):
            result.loc[valid_idx, "BH_adjusted_p"] = benjamini_hochberg(result.loc[valid_idx, "p_value"].tolist())
    return result


def peptide_position_correlations(neighborhoods: pd.DataFrame) -> pd.DataFrame:
    """Repeat correlations after collapsing substitutions at the same peptide site."""
    rows = []
    metrics = [
        "median_abs_MiniGs_effect",
        "median_MiniGs_loss",
        "max_MiniGs_loss",
        "fraction_specific_MiniGs_LoF",
        "median_abs_endocytosis_effect",
        "max_endocytosis_loss",
        "fraction_specific_endocytosis_LoF",
    ]
    for cutoff, subset in neighborhoods.groupby("cutoff_A"):
        position = subset.groupby("peptide_position", as_index=False).agg(
            binding_effect_delta_vs_WT=(
                "binding_effect_delta_vs_WT", "median"
            ),
            binding_effect_magnitude_vs_WT=(
                "binding_effect_magnitude_vs_WT", "median"
            ),
            **{metric: (metric, "first") for metric in metrics},
        )
        for outcome in [
            "binding_effect_delta_vs_WT",
            "binding_effect_magnitude_vs_WT",
        ]:
            for metric in metrics:
                data = position[[outcome, metric]].dropna()
                if len(data) >= 4 and data[metric].nunique() > 1 and data[outcome].nunique() > 1:
                    rho, p_value = stats.spearmanr(data[outcome], data[metric])
                else:
                    rho, p_value = np.nan, np.nan
                rows.append(
                    {
                        "cutoff_A": cutoff,
                        "outcome": outcome,
                        "receptor_neighborhood_metric": metric,
                        "n_unique_peptide_positions": len(data),
                        "spearman_rho": rho,
                        "p_value": p_value,
                    }
                )
    result = pd.DataFrame(rows)
    result["BH_adjusted_p"] = np.nan
    for (cutoff, outcome), idx in result.groupby(["cutoff_A", "outcome"]).groups.items():
        valid_idx = result.loc[idx].index[result.loc[idx, "p_value"].notna()]
        if len(valid_idx):
            result.loc[valid_idx, "BH_adjusted_p"] = benjamini_hochberg(result.loc[valid_idx, "p_value"].tolist())
    return result


def position_predictive_models(neighborhoods: pd.DataFrame) -> pd.DataFrame:
    """Estimate in-sample R2 and leave-one-position-out Q2 for simple mappings."""
    rows = []
    metrics = [
        "median_abs_MiniGs_effect",
        "median_MiniGs_loss",
        "max_MiniGs_loss",
        "fraction_specific_MiniGs_LoF",
        "median_abs_endocytosis_effect",
        "max_endocytosis_loss",
        "fraction_specific_endocytosis_LoF",
    ]
    for cutoff, subset in neighborhoods.groupby("cutoff_A"):
        position = subset.groupby("peptide_position", as_index=False).agg(
            binding_effect_delta_vs_WT=("binding_effect_delta_vs_WT", "median"),
            binding_effect_magnitude_vs_WT=(
                "binding_effect_magnitude_vs_WT", "median"
            ),
            **{metric: (metric, "first") for metric in metrics},
        )
        for outcome in [
            "binding_effect_delta_vs_WT",
            "binding_effect_magnitude_vs_WT",
        ]:
            for metric in metrics:
                data = position[[outcome, metric]].dropna()
                y = data[outcome].to_numpy(float)
                x = data[metric].to_numpy(float)
                if len(data) < 4 or np.unique(x).size < 2 or np.unique(y).size < 2:
                    continue
                design = np.column_stack([np.ones(len(x)), x])
                coefficients = np.linalg.lstsq(design, y, rcond=None)[0]
                fitted = design @ coefficients
                denominator = np.sum((y - y.mean()) ** 2)
                r_squared = 1 - np.sum((y - fitted) ** 2) / denominator
                predictions = np.empty(len(y))
                for i in range(len(y)):
                    keep = np.arange(len(y)) != i
                    loo_coef = np.linalg.lstsq(design[keep], y[keep], rcond=None)[0]
                    predictions[i] = design[i] @ loo_coef
                q_squared = 1 - np.sum((y - predictions) ** 2) / denominator
                rows.append(
                    {
                        "cutoff_A": cutoff,
                        "outcome": outcome,
                        "receptor_neighborhood_metric": metric,
                        "n_unique_peptide_positions": len(data),
                        "intercept": coefficients[0],
                        "slope": coefficients[1],
                        "in_sample_R_squared": r_squared,
                        "leave_one_position_out_Q_squared": q_squared,
                    }
                )
    return pd.DataFrame(rows)


def sensitivity_coverage(neighborhoods: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cutoff, group in neighborhoods.groupby("cutoff_A"):
        in_structure = group[group["in_6VCB_active_peptide"]]
        covered = in_structure["n_assayed_receptor_variants"] > 0
        rows.append(
            {
                "cutoff_A": cutoff,
                "n_peptide_variants_in_structure": len(in_structure),
                "n_variants_with_assayed_receptor_neighbor": int(covered.sum()),
                "fraction_variants_covered": float(covered.mean()),
                "n_unique_peptide_positions_in_structure": int(in_structure["peptide_position"].nunique()),
                "n_unique_positions_covered": int(in_structure.loc[covered, "peptide_position"].nunique()),
            }
        )
    return pd.DataFrame(rows)


def configure_plotting() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.5,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    fig.savefig(output_dir / f"{stem}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(output_dir / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(output_dir / f"{stem}.svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def figure1_peptide_and_contacts(
    peptide: pd.DataFrame, receptor_residue: pd.DataFrame, distances: pd.DataFrame, output_dir: Path
) -> None:
    variants = peptide[peptide["Compound"] != "WT"].sort_values(
        ["peptide_position", "Compound"], na_position="last"
    ).reset_index(drop=True)
    fig = plt.figure(figsize=(7.2, 7.6), constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=[1.1, 1.25])

    ax = fig.add_subplot(grid[0])
    y = np.arange(len(variants))
    colors = np.where(
        variants["binding_effect_delta_vs_WT"] < -0.5,
        "#B23A48",
        np.where(
            variants["binding_effect_delta_vs_WT"] > 0.5,
            "#16857A",
            "#777777",
        ),
    )
    xerr = np.vstack([variants["pct"] - variants["pct.low"], variants["pct.high"] - variants["pct"]])
    for i in range(len(variants)):
        ax.errorbar(
            variants.loc[i, "pct"], y[i], xerr=xerr[:, i : i + 1], fmt="o",
            color=colors[i], ecolor=colors[i], ms=4, capsize=2, lw=1,
        )
    ax.axvline(float(peptide.loc[peptide["Compound"] == "WT", "pct"].iloc[0]), color="#222222", lw=1, ls="--")
    ax.set_yticks(y, variants["Compound"])
    ax.invert_yaxis()
    ax.set_xlim(-2, 103)
    ax.set_xlabel("Binding displacement (%)")
    ax.set_title("a  Experimental GLP-1 variant binding", loc="left", fontweight="bold")
    ax.text(0.58, 0.97, "Dashed line, WT estimate", transform=ax.transAxes, ha="left", va="top", color="#555555")
    outside = np.flatnonzero(~variants["in_6VCB_active_peptide"].to_numpy())
    for idx in outside:
        ax.text(101, idx, "outside 6VCB", va="center", ha="right", fontsize=6.5, color="#555555")

    ax2 = fig.add_subplot(grid[1])
    near = distances[distances["minimum_heavy_atom_distance_A"] <= PROXIMAL_CUTOFF].copy()
    variant_positions = receptor_residue.dropna(subset=["minimum_distance_to_GLP1_A"])
    variant_positions = variant_positions[variant_positions["minimum_distance_to_GLP1_A"] <= PROXIMAL_CUTOFF]
    rec_order = sorted(variant_positions["receptor_position"].unique())
    pep_order = list(range(7, 38))
    rec_index = {p: i for i, p in enumerate(rec_order)}
    pep_index = {p: i for i, p in enumerate(pep_order)}
    residue_lookup = receptor_residue.set_index("receptor_position")
    for row in near.itertuples(index=False):
        if row.receptor_position not in rec_index:
            continue
        effect = residue_lookup.loc[row.receptor_position, "median_MiniGs_effect"]
        color = mpl.colormaps["RdBu_r"](mpl.colors.Normalize(-0.8, 0.8, clip=True)(effect))
        size = 10 + 45 * (PROXIMAL_CUTOFF - row.minimum_heavy_atom_distance_A) / PROXIMAL_CUTOFF
        ax2.scatter(rec_index[row.receptor_position], pep_index[row.peptide_position], s=size, c=[color], edgecolor="white", lw=0.3)
    labels = [
        f"{residue_lookup.loc[p, 'receptor_residue']}{p}" for p in rec_order
    ]
    ax2.set_xticks(range(len(rec_order)), labels, rotation=55, ha="right")
    ax2.set_yticks(range(len(pep_order)), pep_order)
    ax2.invert_yaxis()
    ax2.set_xlabel("GLP1R residue with an assayed human variant")
    ax2.set_ylabel("GLP-1 residue")
    ax2.set_title("b  6VCB interface contacts linked to receptor-variant measurements", loc="left", fontweight="bold")
    ax2.set_xlim(-0.7, len(rec_order) - 0.3)
    ax2.set_ylim(len(pep_order) - 0.3, -0.7)
    sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-0.8, 0.8), cmap="RdBu_r")
    cbar = fig.colorbar(sm, ax=ax2, fraction=0.025, pad=0.02)
    cbar.set_label("Median log Delta Mini-Gs")
    save_figure(fig, output_dir, "Figure_1_peptide_binding_and_interface")


def figure2_receptor_proximity(
    receptor_residue: pd.DataFrame,
    tests: pd.DataFrame,
    output_dir: Path,
    panel_label: str = "b",
) -> None:
    """Plot residue-level GLP1R sensitivity by interface proximity."""

    classes = [
        "Direct (<=4.5 A)",
        "Proximal (4.5-6 A)",
        "Distal (>6 A)",
    ]
    class_labels = [
        "Direct",
        "Proximal",
        "Distal",
    ]

    # Each value represents the median loss across variants assayed
    # at one receptor residue.
    data = [
        receptor_residue.loc[
            receptor_residue["distance_class"].eq(distance_class),
            "median_MiniGs_loss",
        ]
        .dropna()
        .to_numpy(float)
        for distance_class in classes
    ]

    sample_counts = [len(values) for values in data]

    fig, ax = plt.subplots(
        figsize=(3.15, 2.75),
        constrained_layout=True,
    )

    boxes = ax.boxplot(
        data,
        patch_artist=True,
        widths=0.55,
        showfliers=False,
        medianprops={
            "color": "#222222",
            "linewidth": 1.1,
        },
        whiskerprops={
            "color": "#555555",
            "linewidth": 0.8,
        },
        capprops={
            "color": "#555555",
            "linewidth": 0.8,
        },
        boxprops={
            "edgecolor": "#555555",
            "linewidth": 0.8,
        },
    )

    box_colors = [
        "#D8727D",
        "#E7AF72",
        "#A8B3BE",
    ]

    for patch, color in zip(boxes["boxes"], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.80)

    # Overlay individual receptor residues.
    rng = np.random.default_rng(SEED)

    for position, values in enumerate(data, start=1):
        ax.scatter(
            rng.normal(position, 0.055, len(values)),
            values,
            s=13,
            color="#333333",
            alpha=0.60,
            edgecolor="white",
            linewidth=0.25,
            zorder=3,
        )

    tick_labels = [
        f"{label}\n$n$ = {count}"
        for label, count in zip(class_labels, sample_counts)
    ]

    ax.set_xticks(
        [1, 2, 3],
        tick_labels,
    )
    ax.set_xlabel(
        "Distance from GLP-1"
    )
    ax.set_ylabel(
        "Median Mini-Gs loss\nat GLP1R residue"
    )
    ax.set_title(
        f"{panel_label}  GLP1R sensitivity by interface proximity",
        loc="left",
        fontweight="bold",
        fontsize=9,
        pad=7,
    )

    test = tests[
        (tests["background"] == "all resolved receptor residues")
        & (tests["metric"] == "median_MiniGs_loss")
    ]

    if not test.empty:
        p_value = float(test["permutation_p"].iloc[0])

        ax.text(
            0.98,
            0.97,
            (
                "Direct versus non-direct\n"
                rf"permutation $P$ = {p_value:.3g}"
            ),
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=7,
        )

    ax.tick_params(axis="both", labelsize=7)
    ax.margins(y=0.10)

    save_figure(
        fig,
        output_dir,
        "Figure_2b_GLP1R_sensitivity_by_interface_proximity",
    )



def spearman_permutation_test(
    x: np.ndarray,
    y: np.ndarray,
    n_perm: int = N_PERM,
    seed: int = SEED,
) -> tuple[float, float]:
    """Two-sided permutation test for Spearman correlation."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    if len(x) < 4 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan, np.nan
    rank_x = stats.rankdata(x).astype(float)
    rank_y = stats.rankdata(y).astype(float)
    rank_x -= rank_x.mean()
    rank_y -= rank_y.mean()
    denominator = float(np.sqrt(np.sum(rank_x**2) * np.sum(rank_y**2)))
    observed = float(np.dot(rank_x, rank_y) / denominator)
    rng = np.random.default_rng(seed)
    exceed = 0
    for _ in range(n_perm):
        permuted = float(np.dot(rank_x, rng.permutation(rank_y)) / denominator)
        exceed += abs(permuted) >= abs(observed)
    return observed, (exceed + 1) / (n_perm + 1)


def figure3_paired_neighborhoods(
    neighborhoods: pd.DataFrame,
    output_dir: Path,
    absolute_binding_change: bool = False,
) -> None:
    """Compare receptor-neighborhood loss with signed or absolute peptide change.

    Peptide measurements are collapsed by position because all variants at one
    peptide position have the same structural receptor neighborhood. Absolute
    magnitude is calculated for each variant before position-level aggregation,
    so opposing substitutions cannot cancel. This is structural concordance,
    not an epistasis analysis.
    """
    metric = "median_MiniGs_loss"
    signed_outcome = "binding_effect_delta_vs_WT"
    magnitude_outcome = "binding_effect_magnitude_vs_WT"
    outcome = magnitude_outcome if absolute_binding_change else signed_outcome
    subset = neighborhoods[neighborhoods["cutoff_A"] == PROXIMAL_CUTOFF].copy()
    raw_covered = subset.dropna(
        subset=[metric, signed_outcome, magnitude_outcome]
    ).copy()

    covered = (
        raw_covered.groupby("peptide_position", as_index=False)
        .agg(
            position_label=("Compound", lambda s: "/".join(sorted(set(s)))),
            binding_effect_delta_vs_WT=(signed_outcome, "median"),
            binding_effect_magnitude_vs_WT=(magnitude_outcome, "median"),
            median_MiniGs_loss=(metric, "first"),
        )
        .sort_values("peptide_position")
    )

    x = covered[metric].to_numpy(float)
    y = covered[outcome].to_numpy(float)
    signed_y = covered[signed_outcome].to_numpy(float)
    rho, p_value = spearman_permutation_test(x, y)

    linked_receptor_variants = {
        variant
        for value in raw_covered["linked_receptor_variants"].fillna("")
        for variant in value.split(";")
        if variant
    }

    fig, axes = plt.subplots(1, 2, figsize=(6.0, 2.7))
    fig.subplots_adjust(left=0.11, right=0.99, top=0.89, bottom=0.28, wspace=0.42)

    axes[0].axhline(0, color="#888888", linewidth=0.7, zorder=0)
    axes[0].scatter(
        x,
        y,
        c=np.select(
            [signed_y < 0, signed_y > 0],
            ["#B23A48", "#16857A"],
            default="#777777",
        ),
        s=30,
        edgecolor="white",
        linewidth=0.5,
        zorder=2,
    )
    for row in covered.itertuples(index=False):
        axes[0].annotate(
            row.position_label,
            (row.median_MiniGs_loss, getattr(row, outcome)),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=6.2,
        )

    axes[0].set_xlabel("Median Mini-Gs loss among\nneighboring GLP1R variants")
    axes[0].set_ylabel(
        "Absolute GLP-1 binding change\nversus WT (original effect scale)"
        if absolute_binding_change
        else "GLP-1 binding change versus WT\n(original effect scale; signed)"
    )
    axes[0].set_title(
        (
            "a  Magnitude of interface-associated change (6 A)"
            if absolute_binding_change
            else "a  Direction of interface-associated change (6 A)"
        ),
        loc="left",
        fontweight="bold",
    )
    axes[0].text(
        0.03,
        0.97,
        (
            rf"Spearman $\rho$ = {rho:.2f}"
            "\n"
            rf"Permutation $P$ = {p_value:.3f}"
            "\n"
            rf"Peptide positions: {len(covered)}"
            "\n"
            rf"Unique receptor variants: {len(linked_receptor_variants)}"
        ),
        transform=axes[0].transAxes,
        ha="left",
        va="top",
        fontsize=7,
    )
    axes[0].text(
        0.97,
        0.03,
        "Red, decrease; green, increase",
        transform=axes[0].transAxes,
        ha="right",
        va="bottom",
        fontsize=6.2,
        color="#555555",
    )

    plot = subset.sort_values(outcome).reset_index(drop=True)
    colors = {
        "Assayed receptor neighbor": "#2E6F95",
        "Structural contacts, no assayed receptor variant": "#D99A35",
        "No receptor neighbor at threshold": "#8C8C8C",
        "Outside 6VCB peptide": "#D1D1D1",
    }
    axes[1].barh(
        np.arange(len(plot)),
        plot[outcome],
        color=[colors[value] for value in plot["coverage_class"]],
        height=0.72,
    )
    axes[1].axvline(0, color="#333333", linewidth=0.7)
    axes[1].set_yticks(np.arange(len(plot)), plot["Compound"], fontsize=6.2)
    axes[1].set_xlabel(
        "Absolute GLP-1 binding change\nversus WT (original effect scale)"
        if absolute_binding_change
        else "GLP-1 binding change versus WT\n(original effect scale; signed)"
    )
    axes[1].set_title(
        "b  Structural coverage of peptide effects",
        loc="left",
        fontweight="bold",
    )

    handles = [mpl.patches.Patch(color=color, label=label) for label, color in colors.items()]
    fig.legend(
        handles=handles,
        frameon=False,
        loc="lower center",
        ncol=2,
        fontsize=6.2,
        bbox_to_anchor=(0.5, 0.015),
        columnspacing=1.0,
    )
    fig.suptitle(
        (
            "Cross-dataset concordance: absolute peptide-effect magnitude"
            if absolute_binding_change
            else "Cross-dataset concordance: signed peptide effects"
        ),
        fontsize=10,
        fontweight="bold",
        y=0.99,
    )
    stem = (
        "Figure_3_cross_dataset_concordance_absolute"
        if absolute_binding_change
        else "Figure_3_cross_dataset_concordance_signed"
    )
    save_figure(fig, output_dir, stem)



def build_summary(
    peptide: pd.DataFrame,
    receptor: pd.DataFrame,
    residue: pd.DataFrame,
    neighborhoods: pd.DataFrame,
    tests: pd.DataFrame,
    position_correlations: pd.DataFrame,
    predictive_models: pd.DataFrame,
    coverage: pd.DataFrame,
) -> dict:
    direct_residues = residue[residue["distance_class"] == "Direct (<=4.5 A)"]
    prox = neighborhoods[neighborhoods["cutoff_A"] == PROXIMAL_CUTOFF]
    in_structure = prox[prox["in_6VCB_active_peptide"]]
    primary_positions = (
        prox.dropna(
            subset=[
                "binding_effect_delta_vs_WT",
                "binding_effect_magnitude_vs_WT",
                "median_MiniGs_loss",
            ]
        )
        .groupby("peptide_position", as_index=False)
        .agg(
            binding_effect_delta_vs_WT=("binding_effect_delta_vs_WT", "median"),
            binding_effect_magnitude_vs_WT=(
                "binding_effect_magnitude_vs_WT", "median"
            ),
            median_MiniGs_loss=("median_MiniGs_loss", "first"),
        )
    )
    signed_rho, signed_p = spearman_permutation_test(
        primary_positions["median_MiniGs_loss"].to_numpy(float),
        primary_positions["binding_effect_delta_vs_WT"].to_numpy(float),
    )
    absolute_rho, absolute_p = spearman_permutation_test(
        primary_positions["median_MiniGs_loss"].to_numpy(float),
        primary_positions["binding_effect_magnitude_vs_WT"].to_numpy(float),
    )
    signed_model = predictive_models[
        (predictive_models["cutoff_A"] == PROXIMAL_CUTOFF)
        & (predictive_models["outcome"] == "binding_effect_delta_vs_WT")
        & (predictive_models["receptor_neighborhood_metric"] == "median_MiniGs_loss")
    ].iloc[0]
    absolute_model = predictive_models[
        (predictive_models["cutoff_A"] == PROXIMAL_CUTOFF)
        & (predictive_models["outcome"] == "binding_effect_magnitude_vs_WT")
        & (predictive_models["receptor_neighborhood_metric"] == "median_MiniGs_loss")
    ].iloc[0]
    primary_test = tests[
        (tests["background"] == "all resolved receptor residues")
        & (tests["metric"] == "median_MiniGs_loss")
    ].iloc[0]
    return {
        "n_peptide_rows_including_WT": int(len(peptide)),
        "n_peptide_variants": int((peptide["Compound"] != "WT").sum()),
        "n_peptide_variants_in_6VCB_active_segment": int(in_structure.shape[0]),
        "peptide_variants_outside_6VCB": peptide.loc[
            (peptide["Compound"] != "WT") & (~peptide["in_6VCB_active_peptide"]), "Compound"
        ].tolist(),
        "n_receptor_variants": int(len(receptor)),
        "n_unique_receptor_residues": int(receptor["receptor_position"].nunique()),
        "n_receptor_variants_resolved_in_6VCB": int(receptor["minimum_distance_to_GLP1_A"].notna().sum()),
        "n_unique_assayed_receptor_residues_direct": int(len(direct_residues)),
        "n_peptide_variants_with_assayed_receptor_neighbor_at_6A": int(
            (in_structure["n_assayed_receptor_variants"] > 0).sum()
        ),
        "fraction_peptide_variants_covered_at_6A": float(
            (in_structure["n_assayed_receptor_variants"] > 0).mean()
        ),
        "signed_position_level_spearman_rho": float(signed_rho),
        "signed_position_level_permutation_p": float(signed_p),
        "absolute_position_level_spearman_rho": float(absolute_rho),
        "absolute_position_level_permutation_p": float(absolute_p),
        "signed_position_level_in_sample_R_squared": float(
            signed_model["in_sample_R_squared"]
        ),
        "signed_position_level_leave_one_out_Q_squared": float(
            signed_model["leave_one_position_out_Q_squared"]
        ),
        "absolute_position_level_in_sample_R_squared": float(
            absolute_model["in_sample_R_squared"]
        ),
        "absolute_position_level_leave_one_out_Q_squared": float(
            absolute_model["leave_one_position_out_Q_squared"]
        ),
        "direct_vs_other_residue_median_MiniGs_loss_mean_difference": float(
            primary_test["mean_difference"]
        ),
        "direct_vs_other_residue_median_MiniGs_loss_permutation_p": float(
            primary_test["permutation_p"]
        ),
        "coverage_sensitivity": coverage.to_dict(orient="records"),
    }


def main() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    inputs = repo_root / "data" / "receptor-neighborhood" / "inputs"
    parser = argparse.ArgumentParser()
    parser.add_argument("--peptide", type=Path, default=inputs / "glp1_peptide_vars.csv")
    parser.add_argument("--receptor", type=Path, default=inputs / "glp1_receptor_vars.csv")
    parser.add_argument("--structure", type=Path, default=inputs / "6VCB.cif")
    parser.add_argument("--output-dir", type=Path, default=repo_root / "output" / "analysis" / "receptor-neighborhood")
    args = parser.parse_args()
    for input_path in (args.peptide, args.receptor, args.structure):
        if not input_path.is_file():
            raise FileNotFoundError(f"Missing receptor-neighborhood input: {input_path}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir = args.output_dir / "figures"
    table_dir = args.output_dir / "tables"
    figure_dir.mkdir(exist_ok=True)
    table_dir.mkdir(exist_ok=True)

    peptide, wt_effect, wt_pct = parse_peptide_variants(args.peptide)
    receptor = parse_receptor_variants(args.receptor)
    atoms = read_atom_site(args.structure)
    distances = contact_distance_table(atoms)
    receptor = add_receptor_distances(receptor, distances)
    receptor_residue = aggregate_receptor_residues(receptor)
    tests = receptor_proximity_tests(receptor_residue)
    neighborhoods = peptide_receptor_neighborhoods(peptide, receptor, distances)
    correlations = peptide_neighborhood_correlations(neighborhoods)
    position_correlations = peptide_position_correlations(neighborhoods)
    predictive_models = position_predictive_models(neighborhoods)
    coverage = sensitivity_coverage(neighborhoods)

    direct_map = neighborhoods[neighborhoods["cutoff_A"] == DIRECT_CUTOFF][
        ["Compound", "n_assayed_receptor_variants", "n_assayed_receptor_residues"]
    ].rename(
        columns={
            "n_assayed_receptor_variants": "n_direct_assayed_receptor_variants",
            "n_assayed_receptor_residues": "n_direct_assayed_receptor_residues",
        }
    )
    proximal_map = neighborhoods[neighborhoods["cutoff_A"] == PROXIMAL_CUTOFF][
        ["Compound", "coverage_class", "linked_receptor_positions", "linked_receptor_variants"]
    ]
    peptide = peptide.merge(direct_map, on="Compound", how="left").merge(proximal_map, on="Compound", how="left")

    peptide.to_csv(table_dir / "peptide_variants_annotated.csv", index=False)
    receptor.to_csv(table_dir / "receptor_variants_annotated.csv", index=False)
    receptor_residue.to_csv(table_dir / "receptor_residue_summary.csv", index=False)
    distances.to_csv(table_dir / "all_residue_pair_distances.csv", index=False)
    distances[distances["minimum_heavy_atom_distance_A"] <= PROXIMAL_CUTOFF].to_csv(
        table_dir / "interface_contacts_within_6A.csv", index=False
    )
    neighborhoods.to_csv(table_dir / "peptide_receptor_neighborhoods.csv", index=False)
    tests.to_csv(table_dir / "receptor_proximity_permutation_tests.csv", index=False)
    correlations.to_csv(table_dir / "peptide_neighborhood_correlations.csv", index=False)
    position_correlations.to_csv(table_dir / "peptide_position_correlations.csv", index=False)
    predictive_models.to_csv(table_dir / "peptide_position_predictive_models.csv", index=False)
    coverage.to_csv(table_dir / "coverage_sensitivity.csv", index=False)

    configure_plotting()
    figure1_peptide_and_contacts(peptide, receptor_residue, distances, figure_dir)
    figure2_receptor_proximity(receptor_residue, tests, figure_dir)
    figure3_paired_neighborhoods(neighborhoods, figure_dir)
    figure3_paired_neighborhoods(
        neighborhoods, figure_dir, absolute_binding_change=True
    )

    summary = build_summary(
        peptide,
        receptor,
        receptor_residue,
        neighborhoods,
        tests,
        position_correlations,
        predictive_models,
        coverage,
    )
    summary["WT_binding_effect"] = wt_effect
    summary["WT_binding_displacement_pct"] = wt_pct
    with (args.output_dir / "analysis_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, allow_nan=False)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    """

    """
    main()
    
