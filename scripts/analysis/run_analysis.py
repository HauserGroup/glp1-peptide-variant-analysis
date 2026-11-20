"""Run the available calculation scripts from the repository root."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = {
    "nucleotide-variants": "calculate_nucleotide_accessibility.py",
    "af2-rosetta-features": "calculate_af2_rosetta_features.py",
    "af2-contact-differences": "calculate_af2_contact_differences.py",
    "af2-sasa": "calculate_af2_sasa.py",
    "predictor-benchmark": "merge_predictor_scores.py",
    "af3-contact-enrichment": "calculate_af3_contact_enrichment.py",
    "receptor-neighborhood": "calculate_receptor_neighborhood.py",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", choices=SCRIPTS, help="Run one calculation")
    args = parser.parse_args()
    os.environ.setdefault("MPLBACKEND", "Agg")
    selected = {args.only: SCRIPTS[args.only]} if args.only else SCRIPTS
    for name, relpath in selected.items():
        script = Path(__file__).resolve().parent / relpath
        print(f"Running {name}: {script.relative_to(REPO_ROOT)}", flush=True)
        subprocess.run([sys.executable, str(script)], cwd=REPO_ROOT, check=True)
    print("Analysis complete. Results are in output/analysis/.", flush=True)


if __name__ == "__main__":
    main()
