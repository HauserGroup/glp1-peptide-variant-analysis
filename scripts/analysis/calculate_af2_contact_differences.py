"""
Exposed        ●  (OR > 1)
Buried         ●  (OR < 1)
Interface      ●  (OR < 1)
Destabilizing  ●  (OR < 1)
"""
import contextlib
import logging
import os
import pathlib
import re
import pandas as pd
import tempfile

from calculate_af2_rosetta_features import _load_variant_pdbs, _load_wt_pdb, ZENODO_SCAN, ANALYSIS_OUTPUT


def get_position_from_variant(variant):
    """
    H7K -> position 7
    """
    valid_aa = "ACDEFGHIKLMNPQRSTVWY"
    if len(variant) < 3:
        raise ValueError(f"Invalid variant format: {variant}")
    wt_aa = variant[0]
    mut_aa = variant[-1]
    if wt_aa not in valid_aa:
        raise ValueError(f"Invalid wild-type amino acid: {wt_aa} in variant: {variant}")
    if mut_aa not in valid_aa:
        raise ValueError(f"Invalid mutant amino acid: {mut_aa} in variant: {variant}")
    position_str = variant[1:-1]
    if not position_str.isdigit():
        raise ValueError(f"Invalid position: {position_str} in variant: {variant}")
    position = int(position_str)
    return position


def _analyze_one_pdb(
    pdb_path: str,
    out_csv_path: str,
    chains: tuple[str, ...] = ("A", "B"),
) -> str:
    """
    Worker function run in a separate process.
    Writes CSV to out_csv_path and returns pdb_path on success.
    """

    with open(os.devnull, "w") as devnull, contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
        # Import inside worker so each process initializes independently
        from pandaprot import PandaProt

        # Initialize with PDB file
        analyzer = PandaProt(pdb_path, chains=chains)
        # Map all interactions
        interactions = analyzer.map_interactions(
            include_intrachain=False,
            distance_cutoff=4.5,
            )

        # Generate detailed report
        ## Update 13feb2026, had to enable van-der-waals in generate_report, by default they are excluded
        report = analyzer.generate_report(out_csv_path)
        return pdb_path
    
def get_resnum_from_pandaprot_residue(r : str):
    """
    Leu 9 -> 9
    """
    parts = r.split()
    if len(parts) != 2:
        raise ValueError(f"Invalid residue format: {r}")
    resnum_str = parts[1]
    if not resnum_str.isdigit():
        raise ValueError(f"Invalid residue number: {resnum_str} in residue: {r}")
    return int(resnum_str)


def _update_pandaprot_residue_number(res_num, change):
    if change == 0:
        return res_num
    else:
        parts = res_num.split()
        if len(parts) != 2:
            raise ValueError(f"Invalid residue format: {res_num}")
        resname = parts[0]
        resnum_str = parts[1]
        if not resnum_str.isdigit():
            raise ValueError(f"Invalid residue number: {resnum_str} in residue: {res_num}")
        new_resnum = int(resnum_str) + change
        return f"{resname} {new_resnum}"

def filter_pandaprot_report_for_position(report_csv_path, position, 
                                         chain_pair,
                                         peptide_resnum_col,
                                         chain_pos_correction=None):
    """
    in the example: chain_pair is "A-B" for GLP-1R (chain A) and GLP-1 (chain B).
    peptide_resnum_col is "resnum1" if the peptide is chain1 (A)
    and "resnum2" if the peptide is chain2 (B).

    we filter contacts based on the peptide position.
    """
    df = pd.read_csv(report_csv_path)
    chain1, chain2 = chain_pair.split("-")
    chain1_correction = chain_pos_correction[0] if chain_pos_correction else 0
    chain2_correction = chain_pos_correction[1] if chain_pos_correction else 0


    df = df[(df['Chain1'] == chain1) & (df['Chain2'] == chain2)]

    # add this residue change
    df['Residue1'] = df['Residue1'].apply(lambda x: _update_pandaprot_residue_number(x, chain1_correction))
    df['Residue2'] = df['Residue2'].apply(lambda x: _update_pandaprot_residue_number(x, chain2_correction))

    df['resnum1'] = df['Residue1'].apply(get_resnum_from_pandaprot_residue)
    df['resnum2'] = df['Residue2'].apply(get_resnum_from_pandaprot_residue)
    # make sure Chain1==chain1 and Chain2==chain2
    # filter for resnum1 == position
    df = df[df[peptide_resnum_col] == position]
    # add a 'contact pair' which is
    df['contact_pair'] = df['Residue1'] + "|" + df['Residue2'] + "|" + df['Interaction_Type']

    return df

def calculate_contact_count_difference(wt_pdb, variant_pdb):
    """
    Compare WT vs variant contacts for the mutated GLP-1 residue (chain B).
    Adds:
        - new receptor contacts
        - lost receptor contacts
    """
    variant = variant_pdb.stem
    pos = get_position_from_variant(variant_pdb.stem)
    pos = pos - 6  # convert GLP-1 numbering (7-37) to chain-B residue numbering

    # with tempfile.NamedTemporaryFile(suffix=".csv") as wt_csv, tempfile.NamedTemporaryFile(suffix=".csv") as var_csv:
    #     # # # _analyze_one_pdb(wt_pdb, wt_csv.name)
    #     # # # _analyze_one_pdb(variant_pdb, var_csv.name)

    # pre-caculated reports for speed, stored in reference/ and pdbs/ directories
    wt_csv = ZENODO_SCAN / "reference-structure-and-contacts" / "GLP1_7-37_pandaprot" / "pandaprot_results_report.csv"
    var_csv = ZENODO_SCAN / "relaxed-variant-contacts" / f"{variant}_pandaprot" / "pandaprot_results_report.csv"
    if wt_csv.exists() and var_csv.exists():
        wt_df = filter_pandaprot_report_for_position(
            wt_csv,
            pos,
            chain_pair="A-B",
            peptide_resnum_col="resnum2",
            chain_pos_correction=[-23, 0], # correct GPCR residue number.
        )
        var_df = filter_pandaprot_report_for_position(
            var_csv,
            pos,
            chain_pair="A-B",
            peptide_resnum_col="resnum2",
            chain_pos_correction=[0, 0],
        )
    else:
        raise FileNotFoundError(f"Missing pre-calculated pandaprot report CSVs for {variant}. Expected at: {wt_csv} and {var_csv}")

    # deduplicate contact pairs
    wt_pairs = set(wt_df["contact_pair"])
    var_pairs = set(var_df["contact_pair"])

    # unique receptor residues
    wt_rec_residues = set(wt_df["Residue1"])
    var_rec_residues = set(var_df["Residue1"])

    # new / lost receptor residues (cleaner, recommended for analysis)
    lost_residues = wt_rec_residues - var_rec_residues
    if len(lost_residues) < 1:
        lost_residues
    new_residues = var_rec_residues - wt_rec_residues

    return {
        "wt_unique_rec_residues": len(wt_rec_residues),
        "wt_total_unique_contacts": len(wt_pairs),
        "var_unique_rec_residues": len(var_rec_residues),
        "var_total_unique_contacts": len(var_pairs),

        "delta_unique_rec_residues": len(var_rec_residues) - len(wt_rec_residues),
        "delta_total_unique_contacts": len(var_pairs) - len(wt_pairs),

        # new useful metrics
        "n_lost_receptor_contacts": len(lost_residues),
        "n_new_receptor_contacts": len(new_residues),
        # optional detailed lists
        "lost_receptor_contacts": ";".join(sorted(lost_residues)),
        "new_receptor_contacts": ";".join(sorted(new_residues)),
        # keep your original lists
        "wt_contact_list": ";".join(sorted(wt_pairs)),
        "var_contact_list": ";".join(sorted(var_pairs)),
    }



def main():
    ""
    wt_pdb = _load_wt_pdb()
    variant_pdbs = _load_variant_pdbs()
    print(f"Found {len(variant_pdbs)} variant PDBs")
    if not variant_pdbs or not wt_pdb.is_file():
        raise FileNotFoundError(f"Missing Zenodo structural inputs under {ZENODO_SCAN}")
    rows = []
    for i, p in enumerate(variant_pdbs):
        print(f"Processing {p.name} ({i+1}/{len(variant_pdbs)})")
        contact_counts = calculate_contact_count_difference(wt_pdb, p)
        row = {
            'variant': p.stem,
            **contact_counts
        }
        rows.append(row)
    out_df = pd.DataFrame(rows)

    ANALYSIS_OUTPUT.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(ANALYSIS_OUTPUT / "d2_variant_contact_count_differences.csv", index=False)
       


if __name__ == "__main__":
    main()
