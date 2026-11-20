#!/usr/bin/env python3
"""Enumerate GLP-1(7-37) amino-acid variants accessible by one RNA substitution.

Adapted from the script deposited with the study's in-silico data. Uses the
GCG transcript FASTA supplied in zenodo-data and writes a CSV to output/analysis/.
"""

import argparse
from itertools import product
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

BASES = "ACGU"
REPO_ROOT = Path(__file__).resolve().parents[2]
EXPECTED_PROHORMONE = (
    "MKSIYFVAGLFVMLVQGSWQRSLQDTEEKSRSFSASQADPLSDPDQMNEDKRHSQGTFTSDYSKYLDSRRAQDFVQWLMNTKRNRNNIAKRHDEFERHAEGTFTSDVSSYLEGQAAKEFIAWLVKGRGRRDFPEEVAIVEELGRRHADGSFSDEMNTILDNLAARDFINWLIQTKITDRK"
)
EXPECTED_PEPTIDE = "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGRG"


def translate(rna: str) -> str:
    return str(Seq(rna).translate(to_stop=True))


def single_nt_aa_variants(rna: str, first_aa_position: int = 7) -> pd.DataFrame:
    rows = []
    for aa_index in range(len(rna) // 3):
        reference = rna[3 * aa_index:3 * aa_index + 3]
        reference_aa = translate(reference)
        for nt_index, alternate in product(range(3), BASES):
            if alternate == reference[nt_index]:
                continue
            mutated = reference[:nt_index] + alternate + reference[nt_index + 1:]
            alternate_aa = translate(mutated)
            if alternate_aa in ("", reference_aa):
                continue
            position = first_aa_position + aa_index
            rows.append({
                "aa_position": position, "ref_codon": reference,
                "mut_codon": mutated, "ref_aa": reference_aa,
                "alt_aa": alternate_aa, "nt_position_in_codon": nt_index + 1,
                "variant_name": f"{reference_aa}{position}{alternate_aa}",
            })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, default=REPO_ROOT / "zenodo-data/in-silico/01-af2-mutational-scan/single-nucleotide-variant-accessibility/rna.fna")
    parser.add_argument("--output", type=Path,
                        default=REPO_ROOT / "output/analysis/glp1_7_37_single_nucleotide_variants.csv")
    args = parser.parse_args()
    records = list(SeqIO.parse(args.fasta, "fasta"))
    if len(records) != 1:
        parser.error(f"Expected one transcript in {args.fasta}; found {len(records)}")
    cds = str(records[0].seq).upper()[99:642]
    if translate(cds) != EXPECTED_PROHORMONE:
        parser.error("The GCG coding sequence does not match the expected prohormone")
    peptide_rna = cds[97 * 3:128 * 3]
    if translate(peptide_rna) != EXPECTED_PEPTIDE:
        parser.error("The GLP-1(7-37) sequence does not match the expected peptide")
    variants = single_nt_aa_variants(peptide_rna)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    variants.to_csv(args.output, index=False)
    print(f"Wrote {len(variants)} substitutions ({variants.variant_name.nunique()} unique variants) to {args.output}")


if __name__ == "__main__":
    main()
