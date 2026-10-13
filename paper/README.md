# Si/graphite degradation model: code for the paper

This folder reproduces the paper's simulations and figures with the modified
PyBaMM in this repository (`packages/pybamm`, branch
`paper/si-gr-degradation`, tag `si-gr-degradation-v1.0`). One shared
parameter set covers four test conditions:

| Condition | Cell | Temperature | SoC window | Per-condition settings |
|---|---|---|---|---|
| `25deg_full_soc` | CELL064 | 25 °C | full | F0 0.7, width 5e-3 |
| `45deg_full_soc` | CELL017 | 45 °C | full | F0 0.8, width 3e-4 |
| `25deg_15-95_soc` | CELL009 | 25 °C | 15–95% | F0 0.7, width 1e-3 |
| `25deg_20-80_soc` | CELL026 | 25 °C | 20–80% | F0 0.7, width 3e-4 |

- **F0** is the cell's measured first-RPT expansion scale k.
- **width** is the pore-buffering transition width, the one fitted
  per-condition value.

## Contents

| Path | What |
|---|---|
| `reproduce_paper.ipynb` | the notebook: settings, simulation or loading, scores, figures |
| `sigr/conditions.py` | the four conditions and their cycling protocols |
| `sigr/parameters.py` | model options and the shared parameter set (`V6`) |
| `sigr/simulate.py` | formation, aligned C/20 RPTs, C/3 ageing, solver retries, result reduction |
| `sigr/data.py` | experimental-data loading (EFC since each cell's RPT1) |
| `sigr/analysis.py` | scores: SoH gap, voltage RMSE, expansion, k, knee, LLI/LAM, crack state |
| `sigr/plots.py` | figures |
| `tools/run_condition.py` | command-line run of one condition |
| `tools/build_data_bundle.py` | assembles `data/` from the study's files (author use) |
| `environment.yml` | pinned environment |

The model changes to PyBaMM are documented in
`packages/pybamm/si_gr_expansion_precursor/CHANGES.md`, and the degradation
equations in
`packages/pybamm/si_gr_expansion_precursor/pouch_cell_064/conditions_test_matrix/study_si_deg/si_volume_cracks/MODEL_EQUATIONS.md`.

## Install

```bash
conda env create -f environment.yml
conda activate sigr-paper
```

This installs the modified PyBaMM from the tagged commit:
`pip install "git+https://github.com/Dharshannan/PyBaMM.git@si-gr-degradation-v1.0#subdirectory=packages/pybamm"`.

- **Name clash:** it installs as `pybamm`, so use a fresh environment.
- **Solver:** the IDAKLU solver (`pybammsolvers`) is the standard PyPI wheel.

## Data

Download the data record (Zenodo DOI: *to be added*) and unpack it here:

- `data/`: per cell, the capacity fade, LLI, LAM, reversible expansion and
  expansion scale k, plus the C/20 RPT discharge curves. Also the
  beginning-of-life OCPs (`data/ocp/`). EFC = throughput / (2 × 2.5 Ah).
- `outputs/` (optional): the precomputed model outputs
  (`<condition>_results.npz`, `<condition>_rpt_curves.csv`).

## Run

```bash
jupyter lab reproduce_paper.ipynb
```

**In the notebook:**
- `RUN_SIMULATIONS = False` (default) loads `outputs/` and regenerates every
  score and figure.
- `RUN_SIMULATIONS = True` re-simulates: about 1–1.5 h per condition, run in
  parallel processes.

**From the command line,** one condition at a time:

```bash
python tools/run_condition.py 25deg_full_soc --data data --out outputs
```

**C-rate and EFC basis.** The paper runs use the model's capacity,
2.5947 Ah, for both (C/3 = 0.865 A). The tests use the 2.5 Ah nominal
capacity (C/3 = 0.833 A). Set `RATE_BASIS_AH` and `EFC_BASIS_AH` (or
`--rate-basis-ah` / `--efc-basis-ah`) to 2.5 for the test-protocol basis.

## Reproducibility

- **Validated against the study code.** The parameter set matches the
  study's fitting code key for key (201 parameters, relative difference
  < 1e-12) in all four conditions.
- **Tolerance:** small numerical differences across operating systems and
  BLAS builds are expected, of order 0.1 pp in SoH.

## Licence

- **The modified PyBaMM:** BSD-3-Clause. It is based on PyBaMM
  (https://github.com/pybamm-team/PyBaMM) at upstream commit `d63699ddf`;
  see `LICENSE.txt` at the repository root, whose PyBaMM copyright notice applies.
- **This folder:** BSD-3-Clause.
- **The data record:** its own licence (see the record).
