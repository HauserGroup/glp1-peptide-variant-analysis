#!/usr/bin/env python3
"""Run the manuscript figure scripts against the Zenodo data folder."""

import argparse
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
STRUCTURE_SCRIPT = "fig2b_variant_structure.py"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=REPO_ROOT / "zenodo-data")
    parser.add_argument("--output-dir", type=Path, default=REPO_ROOT / "output/figures")
    parser.add_argument(
        "--pymol-python", type=Path,
        help="Python executable whose environment includes the PyMOL module; enables Figure 2b",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    scripts = sorted(SCRIPT_DIR.glob("fig*.py"))
    if not scripts:
        parser.error("No figure scripts found")
    for script in scripts:
        if script.name == STRUCTURE_SCRIPT and args.pymol_python is None:
            print(f"SKIP {script.name}: pass --pymol-python to render the structural panel", flush=True)
            continue
        python = args.pymol_python if script.name == STRUCTURE_SCRIPT else Path(sys.executable)
        suffix = ".png" if script.name == STRUCTURE_SCRIPT else ".svg"
        destination = args.output_dir / (script.stem + suffix)
        start = time.monotonic()
        command = [str(python), str(script), "--data-root", str(args.data_root),
                   "--output", str(destination)]
        subprocess.run(command, check=True)
        print(f"OK   {destination.name} ({time.monotonic() - start:.1f} s)", flush=True)
    print(f"Figure outputs: {args.output_dir.resolve()}", flush=True)


if __name__ == "__main__":
    main()
