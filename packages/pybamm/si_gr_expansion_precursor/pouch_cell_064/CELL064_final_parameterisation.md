# CELL064 — final parameterisation (stages 3 & 4)

**Status: final, for now.** This document is the standing reference for the CELL064
capacity-fade (stage 4) and expansion (stage 3) parameterisation of
`degradation_test_matrix/cell064_degradation_fit.py`. It supersedes the item-by-item
"TUNING STATUS" narrative in that script's docstring as the place to look for *the
current answer*; the docstring remains the detailed lab notebook of how it was reached.

## Adopted parameter set

All values below are the current defaults in `cell064_degradation_fit.py` (each
overridable via its `C64_*` env var for further experimentation):

| parameter | value | role |
|---|---|---|
| `NEG_POROSITY_BOL` | 0.25 | BoL negative-electrode porosity (sets knee onset SoH) |
| `NEG_POROSITY_FLOOR` | 0.023 | porosity floor (sets post-knee collapse depth/steepness) |
| `EXPONENT_MAX_SEI` | 100 | ec-reaction-limited SEI exponent cap |
| `K_SEI_MULT` | 1.0e-3 | SEI kinetic-rate multiplier — the primary **knee-timing** lever |
| `LAM_PROP_MULT` | 0.6 | joint multiplier on Si+Gr stress-driven LAM proportional rate |
| `SI_BETA_LAM_SEI` | 7e-6 | reaction-driven LAM factor for Si |
| `SI_LAM_OPTION` | `"stress and reaction-driven"` | Si's LAM mechanism (built-in PyBaMM option, Ai2019/Reniers2019) |
| `F0_BASELINE` | 0.7 | pore-buffering transmitted-fraction plateau |
| `WIDTH_BASELINE` | 0.05 | pore-buffering transition width |

Everything else not listed here is unchanged from `parameters/cell064_parameters.py`'s
ported BoL fit.

## Fit quality achieved

Run to a natural stop (no solver crash), the model's own EFC now extends to ~230,
past the real dataset's final point (RPT5, EFC ~248 absolute / ~197 shifted):

- **Capacity fade (RPT-vs-RPT):** mean |gap| = **1.9pp** across all 3 validated
  in-range real RPTs (RPT2, RPT4, RPT5) — RPT4 gap **−0.8pp**, RPT5 gap **+2.3pp**.
- **Knee position:** model knee (steepest SoH descent) at **EFC ~108**, against the
  real expansion-peak anchor of **~101** (error +6.3 EFC) — down from ~133 before this
  parameterisation.
- **Reversible expansion shape (normalised):** RMSE **~0.10** against the real
  per-cycle expansion trend, with the pre-knee rise and the knee-peak location both
  well matched (see plot below).
- **Expansion scale k:** matches real closely pre-knee (e.g. model 0.74 vs. real 0.75
  at EFC 46, model 0.88 vs. real 0.95 at EFC 98) before hitting a **structural ceiling
  at k=1** post-knee — see the key finding below.

![CELL064 final fit — degradation and expansion overview](degradation_test_matrix/cell064_overall_summary_result.png)

## Key finding: post-knee expansion is underpredicted because k is capped at ≤1

The single biggest remaining mismatch is that the model's reversible-expansion and
`k` curves flatten out post-knee while the real cell's keep climbing (real `k` reaches
~1.5 by end of life; the model cannot exceed 1.0). This is **not a parameter-fitting
shortfall** — it is a structural property of the current pore-buffering formulation in
`src/pybamm/models/submodels/porosity/reaction_driven_porosity.py`:

```
dv_thickness = dv_solid - dv_buffered
k = x_average(dv_thickness) / x_average(dv_solid)
```

Since `dv_buffered >= 0` whenever the electrode is net swelling (`dv_solid > 0`),
`dv_thickness <= dv_solid` always holds, so `k <= 1` by construction. No combination of
`f0`, `width`, or any of the degradation-rate knobs above can push `k` past 1 — it was
confirmed directly in stage-3 tuning (TUNING STATUS item 17) and re-confirmed across
the entire 4-phase sweep campaign below (every single run's `k` plateaus at exactly
1.0 post-knee, regardless of knob values).

**What real `k` > 1 physically implies:** post-knee, the real cell's *measured* cell
thickness change exceeds its own particle-level swelling — i.e. some additional
expansion mechanism (e.g. plating, crack-opening, or stack-pressure redistribution) is
amplifying the observable expansion beyond what pore-buffering alone can produce.
Capturing that would require **extending the pore-buffering submodel itself** (a
genuine model-physics change), not further parameter tuning. This is left as a clearly
scoped follow-up (see `submodel_stability_audit.md`'s note and TUNING STATUS item 15
for the related, still-deferred porosity-floor softplus-sharpness work, which sits in
the same submodel).

## Model vs. real C/20 RPT discharge curves

As a direct visual check (not just the single-number SoH/voltage comparisons above),
the model's own C/20 RPT discharge curves (voltage vs. discharge capacity) are overlaid
against the real RPT1/RPT2/RPT4/RPT5 discharge curves, each pair matched by nearest EFC
(the model's own RPT cadence drifts from the real RPT schedule post-knee, since EFC
accumulates roughly proportionally to SoH once capacity has faded, not 1:1 with cycle
count — so an exact RPT-index match isn't meaningful late in life):

![CELL064 model vs. real C/20 RPT discharge curves](degradation_test_matrix/cell064_rpt_discharge_overlay.png)

All four pairs matched within ≤4.5 EFC (RPT1: Δ0.5, RPT2: Δ3.5, RPT4: Δ1.0, RPT5:
Δ4.5). **Pre-knee (RPT1, RPT2):** model and real overlap almost exactly, both in
capacity and full voltage-curve shape — confirms the BoL fit and the early-life
degradation-rate calibration are sound. **Post-knee (RPT4, RPT5):** the discharge-
capacity endpoints still line up reasonably closely (within ~0.05 Ah), consistent with
the RPT-vs-RPT gaps above, but the **model's voltage sags noticeably lower through the
mid-discharge plateau (~3.5–3.8 V)** than the real curve, which holds a flatter plateau
before its own drop. So the model is capturing *how much* capacity is lost post-knee
well, but not the full *shape* of how that loss is expressed in the voltage curve —
plausibly a LAM-partition (Gr vs. Si) or SEI-resistance-growth detail rather than the
capacity-fade rate itself. Not chased further here; noted as an open item.

## How this was reached: the 4-phase joint sweep campaign

Six degradation/expansion knobs (BoL porosity, porosity floor, `EXPONENT_MAX_SEI`,
`K_SEI_MULT`, `LAM_PROP_MULT`, `SI_BETA_LAM_SEI`) move the knee, the expansion fit, and
numerical stability together, so they were tuned in one joint campaign rather than
independently:

1. **Phase 1** (33 runs) — 2^(6-1) resolution-VI fractional factorial. Found
   `K_SEI_MULT` dominant and a strong interaction with BoL porosity (high porosity is
   only survivable with high `K_SEI_MULT`).
2. **Phase 2** (28 runs) — 3-level grid around the phase-1 winner. Found the porosity
   floor is the next-most-important knob, with 0.024 clearly better than 0.032.
3. **Phase 3** (26 runs) — finer floor/`k_sei` grid. Revealed that the highest-scoring
   combos were all *crashing* mid-run, before their own RPTs reached the real RPT4/RPT5
   EFCs — meaning their apparently-excellent RPT/k scores were computed against a
   single pre-knee RPT and were effectively meaningless. This was the key methodological
   correction of the campaign.
4. **Phase 4** (10 runs, tie-break) — re-ran the leading `k_sei` values with the cycle
   budget extended so RPT4/RPT5 are back in range. Found `K_SEI_MULT` sets the knee EFC
   in discrete steps (7.7e-4→~133, 0.9e-3→~118, 1.0e-3→~108) and that landing the knee
   at ~101 **and** the best expansion RMSE **and** numerical stability simultaneously is
   not achievable with these six knobs — floor 0.024 gives a better expansion RMSE
   (0.052) but is solver-unstable past EFC~200, while floor 0.023 (adopted here) is
   stable and validated but leaves ~6 EFC of residual knee lateness.

All sweep code and results: `degradation_test_matrix/sweeps/` — shared driver
(`cell064_sweep_lib.py`), one script per phase (`cell064_joint_sweep_phase{1..4}.py`),
a colour-matrix heatmap tool (`cell064_sweep_heatmap.py`), and each phase's results
CSV/heatmap/logs/figures under `sweeps/phase{1,2,3,4}/`.

## Open items

- The k≤1 ceiling above (needs a pore-buffering submodel extension).
- The post-knee RPT discharge-curve **shape** mismatch noted above (capacity endpoint
  matches, mid-plateau voltage doesn't) — plausibly a LAM Gr/Si partition or SEI-
  resistance detail.
- Absolute expansion magnitude is still unresolved (compared via normalised shape only
  — see TUNING STATUS item 17's rejection of `n_electrodes_parallel` as a scaling fix).
- LAM_Si magnitude still undershoots the crude (non-DMA) experimental estimate;
  deprioritised per explicit instruction.
- DMA cross-check on the model's own simulated RPT curves is still not wired.
- Graphite/positive-electrode LAM not yet separately retuned.
