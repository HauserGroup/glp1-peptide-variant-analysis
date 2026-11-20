import pathlib
import re
import subprocess
import pandas as pd

from calculate_af2_rosetta_features import _load_variant_pdbs, _load_wt_pdb, ZENODO_SCAN, ANALYSIS_OUTPUT


def get_position_from_variant(variant):
    m = re.match(r"([A-Z])(\d+)([A-Z])$", variant)
    if not m:
        raise ValueError(f"Invalid variant: {variant}")
    return int(m.group(2))


def run_freesasa(pdb_path):
    result = subprocess.run(
        ["freesasa", "--format=seq", str(pdb_path)],
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout

def get_residue_sasa_from_seq_output(seq_output, chain, resid):
    for line in seq_output.splitlines():
        parts = line.split()
        if len(parts) < 4:
            continue
        if parts[0] != "SEQ":
            continue

        line_chain = parts[1]
        try:
            line_resid = int(parts[2])
        except ValueError:
            continue

        if line_chain != chain or line_resid != resid:
            continue

        # take last column as SASA (robust)
        try:
            sasa = float(parts[-1])
            return sasa
        except ValueError:
            continue

    return None


def classify_sasa(sasa):
    if sasa is None:
        return None
    if sasa < 20:
        return "buried"
    elif sasa < 50:
        return "intermediate"
    else:
        return "exposed"


def calculate_variant_sasa(wt_pdb, variant_pdb, peptide_chain="B"):
    variant = variant_pdb.stem
    pos = get_position_from_variant(variant)
    peptide_resid = pos - 6  # GLP-1 7-37 -> chain-B residue numbering

    wt_out = run_freesasa(wt_pdb)
    var_out = run_freesasa(variant_pdb)

    wt_sasa = get_residue_sasa_from_seq_output(wt_out, peptide_chain, peptide_resid)
    var_sasa = get_residue_sasa_from_seq_output(var_out, peptide_chain, peptide_resid)

    return {
        "variant": variant,
        "position": pos,
        "wt_sasa": wt_sasa,
        "var_sasa": var_sasa,
        "delta_sasa": None if wt_sasa is None or var_sasa is None else var_sasa - wt_sasa,
        "wt_burial_class": classify_sasa(wt_sasa),
        "var_burial_class": classify_sasa(var_sasa),
    }


def main():
    wt_pdb = _load_wt_pdb()
    variant_pdbs = _load_variant_pdbs()
    if not variant_pdbs or not wt_pdb.is_file():
        raise FileNotFoundError(f"Missing Zenodo structural inputs under {ZENODO_SCAN}")

    rows = []
    for i, p in enumerate(variant_pdbs):
        print(f"Processing variant {i+1}/{len(variant_pdbs)}: {p.name}")
        rows.append(calculate_variant_sasa(wt_pdb, p))

    out_df = pd.DataFrame(rows)
    ANALYSIS_OUTPUT.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(ANALYSIS_OUTPUT / "d3_variant_freesasa.csv", index=False)


if __name__ == "__main__":
    main()
