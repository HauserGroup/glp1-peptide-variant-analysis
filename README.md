# GLP-1 peptide variant analysis

Code and selected inputs for the plotted panels in “Population variation in GLP-1 shapes incretin receptor signaling and metabolic function.” The larger supporting dataset is the [Zenodo record](https://doi.org/10.5281/zenodo.23036732).

## Set up

1. Download and extract the Zenodo data into `zenodo-data/` at the repository root. The extracted folder should contain `zenodo-data/in-silico/`, `zenodo-data/in-vivo/`, and the other top-level data folders from the record. This local folder is ignored by Git. The download is assumed complete for the commands below.
2. From the repository root, create a Python environment and install the dependencies:

   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip
   python -m pip install -r requirements.txt
   ```

   Tested with Python 3.11.3 and 3.12.10 on macOS. The requirements install NumPy, pandas, Matplotlib, SciPy, openpyxl, Biopython, seaborn, and statsmodels. The SASA calculation also needs the FreeSASA command-line executable (`brew install freesasa` on macOS). Figure 2b needs PyMOL's Python module, installed separately; the chosen PyMOL Python executable must also have Matplotlib available.

## Run analysis and figures

Run the available calculations first:

```bash
python scripts/analysis/generate_tables.py
```

This writes the nucleotide accessibility table, three AF2 structural tables, the predictor benchmark merge, AF3 contact enrichment, and receptor-neighborhood results under `output/analysis/`. The scripts use the larger Zenodo archive plus small direct inputs in `data/`. See [scripts/analysis/README.md](scripts/analysis/README.md) for each script's direct inputs, output, and verification status. The FEP extractor is excluded until the PMX replica outputs are available.

Then generate the plotted figure panels:

```bash
python scripts/figures/generate_panels.py
```

This writes 15 SVG panels to `output/figures/`. To include the Figure 2b structure PNG, use a Python executable that can import PyMOL:

```bash
python scripts/figures/generate_panels.py --pymol-python /path/to/pymol-python
```

On the tested Mac, that executable is `/Applications/PyMOL.app/Contents/bin/python`. The full run writes 15 SVGs and one PNG. `output/` is ignored by Git. Review the files there before moving approved results into `expected-output/`. The existing Figure 2b reference PNG predates the binding-color and camera corrections in the new output.

The figure scripts read the small, selected inputs in `data/`, which are provided with this repository. Their input provenance is listed in [data/README.md](data/README.md). The local Zenodo copy supplies the FASTA for the analysis command and contains the larger supporting data. Individual figure scripts accept `--data-root` and `--output`; for example:

```bash
python scripts/figures/fig1b_gemme_profile.py --output output/figures/fig1b_gemme_profile.svg
```

## Scope

The calculation scripts in `scripts/analysis/` recreate the deposited AF2 mutation features, contact differences, SASA values, 22-variant benchmark merge, AF3 contact enrichment, and receptor-neighborhood table from their direct inputs. The verified tables matched the Zenodo CSVs exactly. The 16 plotting scripts in `scripts/figures/` read curated tables in `data/` and produce manuscript panels. Figure 2b uses Figure 2a binding values to color mutant side chains; its PNG is a structural layer for the assembled publication panel.

The benchmark merge uses intermediate metric CSVs included in `data/benchmark/metric-inputs/`. The upstream predictor and assay parsing steps require additional raw sources and are not reproduced by this command. AF3 contact enrichment starts from the deposited contact/activity table rather than individual AF3 model files. The AF2 contact calculation uses deposited PandaProt reports. FEP extraction and PMX simulations are outside the tested batch. This repository therefore recreates the named deposited tables from their direct inputs, but not every table from primary instrument or simulation output.

These scripts create plotted panels, not the final assembled publication layouts. Figures 1a and 4a are workflow illustrations; Figure 3a and structural portions of Figures 1d and 2e were rendered or assembled manually. Individual-level All of Us data and unfiltered Regenie output are not included; Figure 4e uses a released aggregate table.

## License and citation

Original code is released under the [MIT License](LICENSE). Data files may have different rights and provenance; consult [data/README.md](data/README.md) and the Zenodo record before redistributing them. If you use this code, cite the manuscript and the repository version used.
