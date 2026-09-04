# Pouch cell 064 — pipeline plan

Status: **not yet started** — this folder is currently an empty skeleton
(`conditions_test_matrix/`, `degradation_test_matrix/`, `expansion_pattern_test/`,
`test_pore_buffering/`, each with a placeholder `.gitkeep`). This document is the plan for
what goes in each, written before the work starts so the procedure is settled up front.

## What this folder is for

`../test_model/` develops and validates every submodel (composite Si/graphite electrode,
SEI/plating/LAM degradation, pore-buffering volume partition, electrolyte dry-out) against
a **synthetic** ageing recipe of this project's own design (`si_gr_expansion`, hand-tuned
multipliers chosen to produce a plausible-looking knee) — there is no real target to fit
to, so "success" there means self-consistency and the mechanisms behaving physically.

This folder repeats the same pipeline on a **real, experimentally-parameterised pouch
cell** (BoL parameters already fitted, from a separate project — to be ported in). The
target changes fundamentally: every submodel's degradation-rate constants get fit
against the cell's actual measured **capacity fade** and **expansion** data, not judged
against a self-designed recipe. Experimental data location: **to be added**.

## What has to change from `test_model/`'s scripts

Noted while surveying `test_model/` before writing this plan — every script there
hardcodes `pybamm.ParameterValues("si_gr_expansion")`; none take the parameter set as a
variable. Porting means either:
- a new parameter set module (analogous to `si_gr_expansion.py`) built from the
  already-fitted BoL parameters, or
- parameterising the existing scripts to accept a chemistry/parameter-set argument,

plus **re-deriving every tuned recipe constant** (`SI_MULT`, `GR_DIV`, `SI_CRIT_STRESS`,
`GR_CRIT_STRESS`, pore-buffering's `eps_min_transfer`/`eps_transfer_width`/
`f_transmit_min`, etc.) — these are fit to the toy chemistry's specific behaviour, not
physically general, and won't transfer as-is to a different cell/electrode design.

## Pipeline stages

Ordered so each stage's output is what the next stage needs; `degradation_test_matrix/`
is the integrative stage everything else feeds into (mirrors how `test_model/
degradation_test_matrix_plan.md` itself builds on the other three).

### 1. `conditions_test_matrix/` — match the experimental protocol

Define the cycling/RPT protocol (C-rates, temperature, RPT schedule, voltage cutoffs) to
match what was actually run on the real cell, so every later simulation is compared
against experimental data collected under the *same* conditions, not a convenient
default.

### 2. `test_pore_buffering/` — expansion submodel, structural validation

Port the pore-buffering (volume-partition) submodel with the ported BoL parameters and
confirm it runs and produces physically sensible within-cycle "breathing" behaviour
before any degradation is added — same two-simulation pattern as `test_model/`'s own
version (`"pore buffering": "false"` reference vs `"true"`), tracking both the internal
`"Negative electrode transfer ratio k"` diagnostic and the dilatometry-style observable
`k = δ_cell/δ_particle`.

### 3. `expansion_pattern_test/` — fit the expansion-precursor transition

Sweep the "physical" transition's `f0` (BoL transmitted-fraction plateau) and
`eps_transfer_width` — as `test_model/expansion_pattern_test/physical_f0_width_grid_test.py`
does — but now against the **experimental** expansion curve: find the parameter values
that reproduce the real cell's own early-life expansion-amplitude behaviour, including
(if present) an amplitude-slope precursor signal ahead of the capacity knee.

### 4. `degradation_test_matrix/` — fit degradation to experimental capacity fade

The integrative stage. Replace `test_model`'s hand-picked recipe-tuning with a genuine
fit: adjust the degradation submodel rate constants (SEI, plating, LAM/cracking, and
dry-out coupling if relevant) until the simulated capacity-fade curve **and** its knee
location match the experimental one. Two sub-pipelines carry over directly:

- **DMA fitting** (`test_DMA/`, ported from `dma_ocp_fit_final.py`/`dma_final_plot.py`):
  the established whole-domain OCV+DV+dV/dQ cost, cross-RPT constraints, and per-RPT
  search strategy (`"restart"` for well-behaved RPTs, `"global"` ensemble+RMSE-selection
  for RPTs whose cost landscape has multiple basins) — this becomes the **production**
  method here, not a comparison against ground truth, since there is no Method A
  (PyBaMM's own state) for a real cell. Fit it directly to the experimental RPT
  discharge (pOCV) curves to extract LLI/LAM_Gr/LAM_Si/LAM_pos as an independent
  cross-check on what the capacity-fade fit attributes degradation to.
- **Dry-out coupling** (`test_dry_out/`, using the (still top-level, unmoved)
  `../ec_dryout/` wrapper): only relevant if the experimental data shows dry-out-related
  behaviour (e.g. resistance growth inconsistent with LAM/LLI alone).

### 5. Cross-validation

Check that one fitted parameter set reproduces **both** the capacity-fade curve and the
expansion curve simultaneously — real degradation mechanisms drive both signals together
(LAM changes both capacity and how much each electrode swells), so fitting them
independently and inconsistently would mean the mechanism attribution is wrong even if
each curve looks fit in isolation.

## Known risk to watch for

`../submodel_stability_audit.md` flags `eps_solid`'s missing `bounds=` in
`loss_active_material.py` as the most likely single point of solver failure near the
capacity knee (plus a couple of other candidate singularities). Real experimental data
is noisier than the synthetic recipe and may expose this harder — worth checking early
rather than after a long fit attempt stalls.

## Inputs still needed

- [ ] BoL-fitted parameter set (from the other project) — to be ported in
- [ ] Experimental capacity-fade data — location to be added
- [ ] Experimental expansion/dilatometry data — location to be added
- [ ] Experimental RPT discharge (pOCV) curves, if available, for the DMA fit
