# Pouch cell 064 — pipeline plan

**See [`CELL064_final_parameterisation.md`](CELL064_final_parameterisation.md) for the
current, standing-reference parameter set and fit quality** (this README's status
narrative below is the detailed history of how it was reached, item by item).

Status: **stage 4 (degradation) has a capacity-fade fit whose knee SoH level, timing, LLI
magnitude, AND full post-knee trajectory (out to and past the real cell's final data
point) all match closely, out of plan order — with one important caveat found since (a
remaining numerical singularity at the porosity floor causes run-to-run non-determinism
in exactly how far a given config reaches, and the two most obvious candidate fixes have
now both been tested and ruled out — see below and items 14-15).** BoL parameters
(`parameters/cell064_parameters.py`) and the experimental aging data
(`degradation_test_matrix/experimental_data/`) are ported in; the degradation-rate
fitting script (`degradation_test_matrix/cell064_degradation_fit.py`) can reach a
NATURAL completion (the SoH-termination floor, not a crash) as far as EFC~232, with its
SoH curve passing almost exactly through both real RPT4 (EFC171, 57.1%) and RPT5
(EFC198, 50.0%) — but this is not a guaranteed outcome of a given config, see items 14-15 —
see that script's module docstring "TUNING STATUS" (items 9-17) for the full history.
**Stage 3 (expansion fitting) has now started too (item 17)** — see the dedicated section
below.
**v2 update (item 16):** the default is now `"stress and reaction-driven"` LAM for Si
(a pre-existing, published PyBaMM option, not new code) with a calibrated
`beta_LAM_sei=15e-6` — this couples Si's active-material loss to its own SEI-on-cracks
current density, which (unlike particle stress, confirmed flat across life) genuinely
jumps ~7x right at the knee. A new `score_rpt_gap()` metric (interpolates the model's
own C/20 RPT trajectory onto each real RPT's EFC) is now the standing fit-quality check
instead of eyeballing the plot; v2's calibrated default gets a 2.0pp mean RPT-vs-RPT gap
(down from ~13-20pp) and LAM_Si=54.2% (up from v1's ~12-33% ceiling), with the pre-knee
fit unchanged.
Four key findings got the fit itself there: (1) CELL064's BoL negative-electrode porosity
(0.35, `cell064_parameters.py`) was defined with no binder/conductive-additive volume
subtracted out, and reducing it to 0.25 (a physically-defensible correction, verified not
to disturb the BoL C/20 fit via `_porosity_bol_check.py`) is what let the knee's SoH
level move at all — every degradation-rate constant tried on its own or jointly (SEI
rate, LAM rate, critical stress, crack rate) had left it stuck at ~85%; (2) two bugs in
the run loop were masking how far the simulation could actually get — a corrupted-data
artifact after a solver error was falsely triggering early termination (fixed via
`MIN_VALID_CAP_AH` filtering), and a genuinely separate unrecoverable-exception failure
mode is now retried at 4 alternate solver tolerances before giving up; (3) the negative-
electrode porosity floor (now 0.021-0.022, still being narrowed) sets how deep the
post-knee collapse can go via the Bruggeman effective-transport relation; (4) the actual
root cause of one recurring post-knee solver crash mode: `src/pybamm/models/submodels/
interface/sei/sei_growth.py`'s "ec reaction limited" SEI flux had an uncapped exponential
term that could numerically overflow near the porosity-floor singularity (a real,
now-fixed core-library bug — the same smooth cap already existed for the sibling
"reaction limited" option, just was never wired into this branch); a second candidate
(`eps_solid` integrating through zero unbounded, per `submodel_stability_audit.md`) was
checked directly via a diagnostic and ruled out — active material fraction never drops
below ~65% of its BoL value even at the crash point. **Found since (item 14) and tested
(item 15, negative result):** a direct diagnostic confirms porosity does pin exactly at
the floor value for many consecutive batches before any crash, but the crash itself is
not deterministic — re-running an identical config gave three different completed-cycle
outcomes (362, 336, 174) — and final-cycles-reached vs. `EXPONENT_MAX_SEI` is
non-monotonic at a fixed floor. The natural next suspect, the porosity floor's own
softplus clamp (`reaction_driven_porosity.py`, sharpness k_eps=100) being too sharp and
creating a stiff Jacobian right at the pin point, was tested directly by temporarily
widening the transition (k_eps=20, k_eps=5) — and does not explain or fix it either:
k_eps=20 still crashed with the same failure signature, just at a different point within
the same noisy range; k_eps=5 was smooth enough that porosity never even reached the
floor within a physically relevant cycle count, which would silently break the knee's
timing rather than fix anything. That temporary override was reverted (core library is
unchanged from upstream). The non-determinism is a confirmed, currently unexplained
property of the marginal/stiff numerical regime at the pinned floor, surviving both of
the two most obvious candidate causes being ruled out. Practical implication:
single-run cycles-reached comparisons between configs are a noisy signal (treat a result
as one draw from a ~170-360+ cycle range for this config), not a clean one. Known open
items: the porosity-floor singularity itself (items 14-15) is not yet fixed, only
characterized, with two candidate causes now ruled out; DMA cross-check still not wired;
graphite/positive-electrode LAM not yet separately retuned. Item 16 (v2) substantially
closed both the RPT-vs-RPT gap and the LAM_Si undershoot that were previously listed
here as open — see item 16 for the calibration sweep and the option to push
`beta_LAM_sei` slightly higher (~17-20e-6) later for an even tighter RPT match if wanted.
Stages 1-3 below are still skeletons
(`.gitkeep` placeholders) — degradation was tackled first per explicit instruction,
ahead of the original stage order, since pore buffering/expansion do not feed back into
the electrochemistry (see stage 4's own note on this) and so can be retuned
independently, after this fit, without redoing it.

## What this folder is for

`../test_model/` develops and validates every submodel (composite Si/graphite electrode,
SEI/plating/LAM degradation, pore-buffering volume partition, electrolyte dry-out) against
a **synthetic** ageing recipe of this project's own design (`si_gr_expansion`, hand-tuned
multipliers chosen to produce a plausible-looking knee) — there is no real target to fit
to, so "success" there means self-consistency and the mechanisms behaving physically.

This folder repeats the same pipeline on a **real, experimentally-parameterised pouch
cell** (BoL parameters already fitted, from a separate project,
`Si_Gr_Expansion_Precursor/Pouch_Data/model_cell064/` — ported into
`parameters/cell064_parameters.py`). The target changes fundamentally: every submodel's
degradation-rate constants get fit against the cell's actual measured **capacity fade**
and **expansion** data, not judged against a self-designed recipe. Experimental aging
data is ported into `degradation_test_matrix/experimental_data/` (see that folder's
`NOTES.md` for provenance and column meaning); the source project's raw discharge/
expansion time series and BoL-fitting scripts stay in
`Si_Gr_Expansion_Precursor/Pouch_Data/` (a separate repo) and are not duplicated here.

## What has to change from `test_model/`'s scripts

Noted while surveying `test_model/` before writing this plan — every script there
hardcodes `pybamm.ParameterValues("si_gr_expansion")`; none take the parameter set as a
variable. **Resolved**: rather than a registered parameter set or a chemistry-argument
refactor of every `test_model/` script, `parameters/cell064_parameters.py` is a plain
module (not a PyBaMM entry point) whose `get_parameter_values()` starts from
`si_gr_expansion`'s own dict and layers CELL064's real BoL geometry/composition/OCPs/
capacity/voltage-window on top via `dict.update()` — every degradation submodel default
(SEI-on-cracks, plating, stress-driven LAM, the project's stability-fix parameters) is
inherited unchanged, matching the instruction not to define a whole new parameter set.
A C/20 formation-discharge smoke test reproduces CELL064's real RPT001 capacity to
within ~0.7% (2.435 Ah simulated vs. 2.452 Ah measured), confirming the port is sound.

Still true and still needed: **re-deriving every tuned recipe rate constant** (`SI_MULT`,
`GR_DIV`, `SI_LAM_PROP`, `GR_LAM_PROP`, `SI_CRIT_STRESS`, `GR_CRIT_STRESS`, pore-
buffering's closure porosity/transition width/plateau, etc.) — these are fit to
`si_gr_expansion`'s own illustrative composition (20% Si by anode volume) and capacity
(5.0 A.h), not physically general, and do not transfer as-is to CELL064's real
composition (8.9% Si by anode volume, 2.59 A.h) or electrode-area-driven geometry. This
is exactly what `degradation_test_matrix/cell064_degradation_fit.py` is for.

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

**Started (item 17 in `degradation_test_matrix/cell064_degradation_fit.py`'s TUNING
STATUS).** Three new diagnostic figures (`cell064_degradation_fit_result.png` — SoH/LAM/
LLI/RPT-voltage, `cell064_expansion_fit_result.png` — reversible expansion + k,
`cell064_overall_summary_result.png` — everything combined) and two new standing
fit-quality metrics (`score_expansion_shape`, `score_k_gap`, alongside item 12's
`score_rpt_gap`) are in place. Key findings so far:

- **Absolute expansion magnitude is not yet resolved.** `"Number of electrodes connected
  in parallel to make a cell"` (inherited un-overridden at 1.0) was the obvious first
  candidate for scaling the model's single-representative-electrode-pair thickness
  change up to the real pouch cell's whole-stack dilatometry reading (~8-10x too small
  as-is) — **rejected before being tested**, on tracing it: it's actually the
  current-collector area term (`A_cc = L_y*L_z*n_electrodes_parallel` in
  `geometric_parameters.py`), which drives capacity and every C-rate/current-density
  term throughout the electrochemistry, not a thickness-only knob — scaling it up would
  silently invalidate the whole stage-4 degradation fit. No physical layer count is
  documented anywhere in the source project's data either. Adopted instead: compare
  **normalised** reversible expansion (each series divided by its own early-life value)
  so shape can be fit now without resolving magnitude.
- **`f0`/`width` grid swept** (`f0` in {0.7, 0.55, 0.4} x `width` in {0.01, 0.05, 0.12})
  against the normalised expansion shape and `k`: `f0=0.7` (the inherited default)
  dominates every width tried; within it, `width=0.05` (adopted, was 0.01) gives the
  best shape match (RMSE 0.096, vs. 0.151 at the old default) and a close second-best `k`
  match. Pre-knee normalised expansion now tracks the real near-flat trend closely, and
  the knee-peak + post-knee decline both track well.
- **`k` has a structural ceiling the current pore-buffering formulation cannot exceed:**
  `k = dv_thickness/dv_solid <= 1` by construction whenever the electrode is net
  swelling, but real `k` rises to ~1.5 post-knee. Model `k` matches real closely at the
  first three RPTs (0.705/0.72/0.945 vs. real 0.705/0.756/0.949) then hits this ceiling.
  Not chased further per instruction — would need extending the pore-buffering submodel
  itself (genuine new physics), not a parameter fit.

**Phased joint sweep (`sweeps/`, item 19):** a multi-phase campaign over the six
degradation/expansion knobs — BoL negative-electrode porosity (`C64_NEG_POROSITY_BOL`),
porosity floor (`C64_NEG_POROSITY_FLOOR`), SEI exponent cap (`C64_EXPONENT_MAX_SEI`),
SEI kinetic-rate multiplier (`C64_K_SEI_MULT`), joint stress-driven-LAM proportional-rate
multiplier (`C64_LAM_PROP_MULT`, new), and reaction-driven-LAM factor
(`C64_SI_BETA_LAM_SEI`) — to improve the C/20 RPT capacity-fade fit, the reversible-
expansion fit, **and the knee position** together (several knobs, BoL porosity and the
SEI rate especially, move the knee, so they can't be fit independently of the fade
curve). Layout: shared machinery in `sweeps/cell064_sweep_lib.py`; one driver script per
phase (`cell064_joint_sweep_phase{1,2,3}.py`); results, logs and figures per phase in
`sweeps/phase{1,2,3}/`; colour matrix via `cell064_sweep_heatmap.py <phaseN>`. Every run
is a headless (`C64_NO_PLOT=1`) subprocess of `cell064_degradation_fit.py`,
`C64_MAX_CYCLES` capped at 250 (nothing to fit past EFC ~200). New standing metric
`score_knee` (knee-EFC + onset-SoH error vs the ~101/~91% anchors) and a pure knee-timing
knob `C64_TIMESCALE_MULT` (scales `K_SEI_MULT`/`LAM_PROP_MULT`/`SI_BETA_LAM_SEI` together)
were added for this.
- **Phase 1** — 2^(6-1) res-VI fractional factorial (32 runs + baseline). `K_SEI_MULT`
  is the dominant knob (~4–6× the rest); BoL porosity 0.28 is toxic unless paired with
  high `K_SEI_MULT`; no combo beats the baseline on all of RPT+expansion+k. Best "earlier
  knee + expansion" corner: **run31** (poro 0.28 / floor 0.028 / exp_max 100 / `k_sei`
  1.24e-3 / `lam` 0.6 / `beta` 7e-6) — best expansion RMSE in the sweep (0.069), knee at
  EFC ~98 (anchor ~101).
- **Phase 2** — 3-level grid around run31 (`k_sei` × BoL porosity × floor). Best: **run01**
  (`k_sei` 1.0e-3 / BoL 0.25 / floor 0.024) — expansion RMSE 0.052, knee EFC ~108. Floor
  is now the dominant knob (0.024 much better than 0.032).
- **Phase 3** — finer brackets (floor to 0.020, `k_sei` 0.9–1.6e-3). Confirmed floor
  0.023–0.024 as the bracket and BoL 0.25 as best — **but** revealed that every high-`k_sei`
  top combo *crashes* at batch ~5, so its RPT/k scores are computed against a single
  pre-knee RPT and are near-vacuous (the post-knee capacity fit is untested and actually
  too fast).
- **Phase 4** (tie-break) — `k_sei` {8e-4, 9e-4, 1.0e-3} × floor {0.023–0.025} with
  `MAX_CYCLES` raised to 350 so RPT4/RPT5 are back in range. `k_sei` sets the knee EFC in
  discrete steps (7.7e-4→133, 0.9e-3→118, 1.0e-3→108); getting the knee to ~101 *and* the
  best expansion RMSE *and* stability is not achievable with these 6 knobs.

**Adopted (item 19 final):** the **run07** recipe — `C64_K_SEI_MULT` **1.0e-3**,
`C64_NEG_POROSITY_FLOOR` **0.023**, `C64_EXPONENT_MAX_SEI` **100**, `C64_LAM_PROP_MULT`
**0.6**, `C64_SI_BETA_LAM_SEI` **7e-6** (BoL porosity stays 0.25). These five are now the
defaults in `cell064_degradation_fit.py`. Stable natural completion, validated through
RPT4 (gap −0.8pp), RPT mean gap 1.7pp (= old default's 1.6), expansion RMSE 0.102 (old
0.096), and the **knee moves from EFC ~133 to ~108** (error vs. the ~101 anchor: +32 →
+6.3) — a much better-aligned knee at essentially no cost. The residual ~6-EFC lateness
and the expansion-RMSE-vs-stability tradeoff need the deferred porosity-floor softplus
change (item 15), not more knob tuning. All sweep artefacts: `sweeps/phase{1,2,3,4}/`.

`f0`/`width` are NOT in this sweep — settled in item 17.

Not yet done: resolving absolute expansion magnitude (electrode-count question still
open); an amplitude-slope precursor-signal check ahead of the capacity knee.

**Next steps queued for a future session (item 18 — documented only, not yet tried; the
current fit is already judged good enough for publication as-is):** push the capacity-
fade knee a bit earlier and make it a bit smoother (currently a fairly sharp bend), via
some joint combination of the degradation-rate constants (`K_SEI_MULT`, `SI_LAM_PROP`/
`GR_LAM_PROP`, `SI`/`GR_CRIT_STRESS`, `SI`/`GR_LAM_EXP`), `NEG_POROSITY_FLOOR`, and
`EXPONENT_MAX_SEI` — see item 18 in the script's docstring for the full reasoning
(earlier timing and a smoother bend are somewhat in tension with each other, so this
will likely need at least 2 knobs tuned jointly, not a single-parameter fix).

### 4. `degradation_test_matrix/` — fit degradation to experimental capacity fade

**Started.** `cell064_degradation_fit.py` reuses the exact submodel coupling validated in
`../test_model/test_pore_buffering/pore_buffering_degradation_test_1000cyc.py` (SEI
reaction-limited + on-cracks + porosity change, stress-driven LAM, particle swelling on
graphite / swelling-and-cracking on silicon) on top of `cell064_parameters.py`, cycled
under CELL064's own real protocol (continuous ~C/3 ageing, C/20 RPT every ~50 EFC, real
voltage window). Pore buffering is deliberately left off for this stage — it only
repartitions expansion for the thickness-change diagnostic and has no feedback into the
reaction/porosity equations, so it cannot affect capacity/LLI/LAM (confirmed by reading
the submodel's `get_coupled_variables`, not assumed) — meaning stage 2/3 can be done
after this stage without invalidating it. Replace `test_model`'s hand-picked recipe-
tuning with a genuine fit: adjust the degradation submodel rate constants (SEI, LAM/
cracking, and dry-out coupling if relevant) until the simulated capacity-fade curve
**and** its knee location (EFC ~152, from the expansion-peak anchor) match the
experimental one — not yet converged, see the script's docstring. Two sub-pipelines
carry over directly:

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

- [x] BoL-fitted parameter set — ported (`parameters/cell064_parameters.py`)
- [x] Experimental capacity-fade data — ported (`degradation_test_matrix/experimental_data/`)
- [x] Experimental expansion/dilatometry data — ported (same folder; not yet used, see stage 3)
- [ ] Experimental RPT discharge (pOCV) curves — exist in the source project
  (`Pouch_Data/cell064_data/data_timeseries/low_rate_c20/`) but not ported; not needed for
  the DMA cross-check as currently planned (DMA runs on the *model's own* simulated RPT
  curves, not the real ones — see stage 4), so only worth porting if a direct real-vs-
  simulated RPT voltage-curve comparison is wanted later.
