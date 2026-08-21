---
name: project-conda-env
description: Which Python/conda environment to use when running scripts against this local pybamm checkout
metadata:
  type: project
---

Scripts in this repo (e.g. `si_gr_expansion_precursor/**/*.py`) must be run with the `PyBaMM_Pressure` conda environment's Python, not the base anaconda Python on PATH.

**Why:** The default `python` on PATH resolves to base anaconda, which has a separately pip-installed `pybamm` package (missing all local-repo changes, e.g. custom parameter sets like `si_gr_expansion`, the pore-buffering submodel, etc.). The `PyBaMM_Pressure` env has pybamm installed in editable mode pointing at `C:\Shannan_PhD_Stuff_Local\PyBaMM_Pressure\packages\pybamm\src\pybamm`, so it picks up local source edits immediately.

**How to apply:** Always invoke scripts as `"/c/Users/ds3420/AppData/Local/anaconda3/envs/PyBaMM_Pressure/python" <script.py>` (Bash/Git-Bash path form) rather than bare `python`. If a bare `python` run raises `ValueError: '<param set>' is not a valid parameter set` or otherwise seems to be missing recent local changes, this environment mismatch is the first thing to check.
