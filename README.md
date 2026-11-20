# GLP-1 peptide variant analysis

Code for the plotted panels in “Population variation in GLP-1 shapes incretin receptor signaling and metabolic function.” The larger supporting dataset is the [Zenodo record](https://zenodo.org/records/23036732?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6Ijk4MWM3OTkyLWNkMjktNDY4NS05NTQ1LWYyNDMxNTc4ZjhiYSIsImRhdGEiOnt9LCJyYW5kb20iOiIwOWE4NmY0OGQ2NzMyNmEwMmE3MWRlNmQyZGQyYzE5MiJ9.-_RtzewdC_hv0NHv2Dg38yNXHvpxjdgmrHkU7hYwUwUW0-rEUqVAf2G235Vk8k9ixBIAI3lPlI2e8UvFjH92hg).

## Set up

1. Download and extract the Zenodo data into `zenodo-data/` at the repository root. The extracted folder should contain `zenodo-data/in-silico/`, `zenodo-data/in-vivo/`, `zenodo-data/in-vitro/`, and `zenodo-data/figure-source-data/`. This local folder is ignored by Git.

2. From the repository root, create a Python environment and install the dependencies:

   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip
   python -m pip install -r requirements.txt
   ```

   Tested with Python 3.11.3 and 3.12.10 on macOS. The requirements install NumPy, pandas, Matplotlib, SciPy, openpyxl, Biopython, seaborn, and statsmodels. The SASA calculation also needs the FreeSASA command-line executable (`brew install freesasa` on macOS). Figure 2b needs PyMOL's Python module, installed separately; it does not need Matplotlib in the PyMOL Python environment.

## Run analysis and figures

Run the available calculations first:

```bash
python scripts/analysis/run_analysis.py
```

This writes the nucleotide accessibility table, three AF2 structural tables, the predictor benchmark merge, AF3 contact enrichment, and receptor-neighborhood results under `output/analysis/`. All direct inputs are read from `zenodo-data/`. On the tested laptops runtime was less than 2 minutes.


Then generate the plotted figure panels:

```bash
python scripts/figures/generate_panels.py
```
This writes 15 SVG panels to `output/figures/`. On the tested laptops runtime was less than 1 minute. To include the Figure 2b structure PNG, use a Python executable that can import PyMOL:

```bash
python scripts/figures/generate_panels.py --pymol-python /path/to/pymol-python
```

To render only Figure 2b with the PyMOL Python executable:

```bash
/path/to/pymol-python scripts/figures/fig2b_variant_structure.py
```
Rendering of variants with the Pymol took less than 3 minutes on the tested laptops.
# for example, on mac: /Applications/PyMOL.app/Contents/bin/python scripts/figures/fig2b_variant_structure.py


## Scope
The calculation scripts in `scripts/analysis/` recreate the deposited AF2 mutation features, contact differences, SASA values, 22-variant benchmark merge, AF3 contact enrichment, and receptor-neighborhood table from their direct inputs. The generated tables were checked against Zenodo. With pandas 3.1, the benchmark and receptor tables differ only by floating-point roundoff (at most 1.5 × 10⁻¹⁴); the other four tested tables match exactly. The 16 plotting scripts in `scripts/figures/` read tables from `zenodo-data/` and produce manuscript panels. Figure 2b reproduces the colors, alignments, and camera of the supplied PyMOL scene using PDBs from Zenodo; its PNG is a structural layer for the assembled publication panel.

The benchmark merge uses intermediate metric CSVs in `zenodo-data/in-silico/02-predictor-benchmark/metric-inputs/`. The upstream predictor and assay parsing steps require additional raw sources and are not reproduced by this command. AF3 contact enrichment starts from the deposited contact/activity table rather than individual AF3 model files. The AF2 contact calculation uses deposited PandaProt reports. FEP extraction and PMX simulations are outside the tested batch. This repository therefore recreates the named deposited tables from their direct inputs, but not every table from primary instrument or simulation output.

These scripts create plotted panels, not the final assembled publication layouts. Figures 1a and 4a are workflow illustrations; Figure 3a and structural portions of Figures 1d and 2e were rendered or assembled manually. Individual-level All of Us data and unfiltered Regenie output are not included; Figure 4e uses a released aggregate table.


## License and citation

Original code is released under the [MIT License](LICENSE). Data files may have different rights and provenance; consult the Zenodo data README and record before redistributing them. If you use this code, cite the manuscript and the repository version used.
