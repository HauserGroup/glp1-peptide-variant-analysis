"""
Parse the rosetta relax scores from the relaxed PDBs of the GLP-1 variants, and compute features like:
- whole-pose total score change
- mutated-residue total score change
- per-term score changes at the mutated residue (fa_atr, fa_rep, fa_sol, fa_elec, fa_dun)
- hydrogen-bond summary change at the mutated residue (hbond_sr_bb + hbond_lr_bb + hbond_bb_sc + hbond_sc)

"""

import pandas as pd
import pathlib
import re

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
ZENODO_SCAN = REPO_ROOT / 'zenodo-data' / 'in-silico' / '01-af2-mutational-scan'
ANALYSIS_OUTPUT = REPO_ROOT / 'output' / 'analysis' / 'af2-mutational-scan'


ROSETTA_COLS = """
label fa_atr fa_rep fa_sol fa_intra_rep fa_intra_sol_xover4 lk_ball_wtd
fa_elec pro_close hbond_sr_bb hbond_lr_bb hbond_bb_sc hbond_sc
dslf_fa13 omega fa_dun p_aa_pp yhh_planarity ref rama_prepro total
""".split()

NUMERIC_COLS = [c for c in ROSETTA_COLS if c != "label"]


def _load_variant_pdbs():
    pdbs = ZENODO_SCAN / "relaxed-variant-structures"
    return sorted(pdbs.glob("*.pdb"))


def _load_wt_pdb():
    return ZENODO_SCAN / "reference-structure-and-contacts" / "GLP1_7-37.pdb"


def _parse_mutation_from_filename(pdb_file):
    """
    Extract WT AA, position, MUT AA from filenames like:
    L14R.pdb, Y13R_relaxed.pdb, etc.
    """
    m = re.search(r"([A-Z])(\d+)([A-Z])", pdb_file.stem)
    if not m:
        raise ValueError(f"Could not parse mutation from filename: {pdb_file.name}")
    wt_aa, pos, mut_aa = m.group(1), int(m.group(2)), m.group(3)
    return wt_aa, pos, mut_aa


def _get_rosetta_relax_scores_from_pdb(pdb_file):
    start_line = "#BEGIN_POSE_ENERGIES_TABLE"
    end_line = "#END_POSE_ENERGIES_TABLE"

    with open(pdb_file, "r") as f:
        lines = f.readlines()

    in_table = False
    table_lines = []
    for line in lines:
        if line.startswith(start_line):
            in_table = True
            continue
        if line.startswith(end_line):
            break
        if in_table:
            table_lines.append(line.strip())

    df = pd.DataFrame([line.split() for line in table_lines], columns=ROSETTA_COLS)
    df = df[df["label"] != "label"].copy()

    # Pandas 3 infers string dtype here; object dtype permits numeric assignment below.
    weights = df[df["label"] == "weights"].iloc[0].astype(object).copy()
    pose = df[df["label"] == "pose"].iloc[0].astype(object).copy()
    df = df[~df["label"].isin(["weights", "pose"])].copy()

    # convert numeric columns
    for col in NUMERIC_COLS:
        df[col] = pd.to_numeric(df[col], errors="coerce")
        if col != "total" or pose[col] != "NA":
            try:
                pose[col] = float(pose[col])
            except Exception:
                pass
            try:
                weights[col] = float(weights[col])
            except Exception:
                pass

    return df, weights, pose


def _extract_position_from_label(label):
    """
    Examples:
    ARG:NtermProteinFull_1 -> 1
    PRO_2 -> 2
    GLN_3 -> 3
    """
    m = re.search(r"_(\d+)$", label)
    if not m:
        raise ValueError(f"Could not extract residue position from label: {label}")
    return int(m.group(1))


def _add_position_column(df):
    df = df.copy()
    df["position"] = df["label"].apply(_extract_position_from_label)
    return df


def _get_row_for_position(df, position):
    hits = df[df["position"] == position]
    if len(hits) != 1:
        raise ValueError(f"Expected exactly one row for position {position}, found {len(hits)}")
    return hits.iloc[0]


def _sum_hbond_terms(row):
    return (
        row["hbond_sr_bb"]
        + row["hbond_lr_bb"]
        + row["hbond_bb_sc"]
        + row["hbond_sc"]
    )


def _build_variant_feature_row(variant_pdb, wt_df, wt_pose):
    wt_aa, position, mut_aa = _parse_mutation_from_filename(variant_pdb)

    mut_df, _, mut_pose = _get_rosetta_relax_scores_from_pdb(variant_pdb)
    mut_df = _add_position_column(mut_df)

    wt_row = _get_row_for_position(wt_df, position)
    mut_row = _get_row_for_position(mut_df, position)

    out = {
        "variant": f"{wt_aa}{position}{mut_aa}",
        "position": position,
        "wt_aa": wt_aa,
        "mut_aa": mut_aa,
        "pdb_file": variant_pdb.name,

        # whole-pose change
        "wt_pose_total": wt_pose["total"],
        "mut_pose_total": mut_pose["total"],
        "delta_pose_total": mut_pose["total"] - wt_pose["total"],

        # mutated-residue change
        "wt_res_total": wt_row["total"],
        "mut_res_total": mut_row["total"],
        "delta_res_total": mut_row["total"] - wt_row["total"],

        # per-term deltas at mutated position
        "delta_fa_atr": mut_row["fa_atr"] - wt_row["fa_atr"],
        "delta_fa_rep": mut_row["fa_rep"] - wt_row["fa_rep"],
        "delta_fa_sol": mut_row["fa_sol"] - wt_row["fa_sol"],
        "delta_fa_elec": mut_row["fa_elec"] - wt_row["fa_elec"],
        "delta_fa_dun": mut_row["fa_dun"] - wt_row["fa_dun"],

        # hydrogen-bond summary
        "wt_hbond_sum": _sum_hbond_terms(wt_row),
        "mut_hbond_sum": _sum_hbond_terms(mut_row),
        "delta_hbond_sum": _sum_hbond_terms(mut_row) - _sum_hbond_terms(wt_row),

        # optional extra interpretable terms
        "delta_pro_close": mut_row["pro_close"] - wt_row["pro_close"],
        "delta_rama_prepro": mut_row["rama_prepro"] - wt_row["rama_prepro"],
        "delta_ref": mut_row["ref"] - wt_row["ref"],
    }

    return out


def classify_structural_impact(df, score_col="delta_res_total"):
    df = df.copy()
    q90 = df[score_col].quantile(0.90)
    q75 = df[score_col].quantile(0.75)

    def _label(x):
        if x >= q90:
            return "strongly_disruptive"
        elif x >= q75:
            return "disruptive"
        else:
            return "permissive"

    df["impact_class"] = df[score_col].apply(_label)
    return df



def main():
    variants = _load_variant_pdbs()
    wt_pdb = _load_wt_pdb()
    if not variants or not wt_pdb.is_file():
        raise FileNotFoundError(f"Missing Zenodo structural inputs under {ZENODO_SCAN}")

    wt_df, wt_weights, wt_pose = _get_rosetta_relax_scores_from_pdb(wt_pdb)
    wt_df = _add_position_column(wt_df)

    rows = []
    for p in variants:
        try:
            row = _build_variant_feature_row(p, wt_df, wt_pose)
            rows.append(row)
        except Exception as exc:
            raise RuntimeError(f"Failed to calculate Rosetta features for {p.name}") from exc

    out_df = pd.DataFrame(rows).sort_values(["position", "variant"])
    out_df = classify_structural_impact(out_df, score_col="delta_res_total")
    print(out_df.head())
    print(out_df["impact_class"].value_counts())
    ANALYSIS_OUTPUT.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(ANALYSIS_OUTPUT / "d1_mutation_rosetta_delta_features.csv", index=False)


if __name__ == "__main__":
    main()