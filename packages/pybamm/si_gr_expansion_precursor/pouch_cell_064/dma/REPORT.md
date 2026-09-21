# CELL064 DMA + Si diffusivity investigation (2026-09-19)

## 1. Real-data DMA: does CELL064 actually show a severe Si stoichiometric collapse?

Pulled directly from Pouch_Data's own already-validated `discharge_esoh_fit.py`
outputs (raw RPT discharge CSVs already carry a `model_x_si`/`model_x_gr`
column -- no fresh fit needed). Si's per-RPT window = `max(model_x_si) -
min(model_x_si)` over one full discharge (`extract_real_si_sto_window.py`).

| RPT | EFC | source | x_si_delta | norm. to RPT1 |
|---|---|---|---|---|
| 0 | 0.8 | hi-rate neighbor | 0.927 | 0.94 |
| 1 | 50.6 | C/20 RPT | 0.990 | 1.00 |
| 2 | 100.0 | C/20 RPT | 0.981 | 0.99 |
| 3 | 148.5 | hi-rate neighbor | 0.926 | 0.94 |
| 4 | 222.0 | C/20 RPT | 0.780 | 0.79 |
| 5 | 248.0 | C/20 RPT | 0.920 | 0.93 |

**Real collapse ratio: 1.27x** (worst point RPT4, recovers by RPT5). This
was cross-checked against `CELL064_LAM_derived.csv`: real Si LAM (Cn_Si) IS
large and real (58.5% at RPT4, 79.6% at RPT5) -- there's no contradiction,
because LAM_Si (a *size* of the accessible-capacity pool) and x_si_delta
(the *fractional occupancy* of whatever pool remains) are independent axes
in the eSOH decomposition. Isolated Si drops out of the pool entirely
rather than staying counted but under-cycled.

## 2. Model check: does the model's OWN C/20 RPT (not just its C/3 ageing
cycle) also show this collapse?

`extract_model_si_sto_window_c20rpt.py` re-ran the full CELL064 degradation
model (same config as the session's other diagnostics) and classified each
cycle's discharge leg exactly the way `cell064_degradation_fit.py`'s own
`rpt_leg()`/`cycle_ageing_leg()` do (RPT: 0.05 < mean_I <= 0.5 A).

| model RPT cycle | EFC | si_delta_sto | eps_s_min |
|---|---|---|---|
| 0 | 0.5 | 0.984 | 0.058 |
| 152 | 130.2 | 0.974 | 0.040 |
| 202 | 157.9 | 0.927 | 0.024 |
| 302 | 198.9 | 0.647 | 0.013 |
| **352** | **217.9** | **0.182** | **0.011** |
| 502 | 269.0 | 0.106 | 0.007 |

**Model C/20 RPT collapse ratio: 9.3x** -- and critically, the cliff is
sharp and LATE: the window stays near-full (>0.65) all the way to EFC~199,
then collapses to 0.18 by EFC~218, right where real RPT4 (EFC 222) sits.
So this is NOT purely a fast-rate (C/3) artifact -- even the model's own
slow diagnostic discharge shows the collapse once `eps_s` (active Si
fraction) has crashed far enough. This is exactly consistent with your
theory: once enough Si is isolated, the *remaining* active fraction sees
high enough local current density that even a nominal C/20 test becomes
locally rate-limited.

See `dma_si_sto_window_comparison.png` for real vs. model (C/20 RPT and C/3
ageing legs) overlaid.

## 3. Si diffusivity + item-30 sweep

Base value `silicon_LGM50_diffusivity_Bonkile2024` (D_ref=3.0e-16 m2/s) is a
literature placeholder, never fitted to CELL064 (no high-C-rate pouch data
available). Added `C64_SI_DIFFUSIVITY_MULT` (wraps the function with a
scalar) and reinstated item 30 (`"active material expansion residual"`,
isolated-Si fraction still contributes to volume-change signal) behind a
new option, both smoke-tested before the full run.

| config | expansion RMSE | C/20 RPT collapse ratio | capacity-fade RPT gap (pp) |
|---|---|---|---|
| baseline (x1) | 0.307 | 9.25x | **0.79** (excellent) |
| **x3** | **0.149** (best) | 7.85x | 11.25 |
| x10 | 0.172 | 3.61x | 13.07 |
| x30 | 0.174 | 3.53x | 13.20 |
| x100 | 0.175 | 3.39x | 13.21 |
| x3 + residual (frac=0.5) | 0.189 | 7.85x | 11.27 |
| x3 + residual (frac=0.9) | 0.243 | 7.85x | 11.28 |
| residual only (x1, frac=0.9) | 0.228 | 9.26x | 0.78 |

**Findings:**

- **The diffusivity fix works, substantially**: expansion-shape RMSE roughly
  HALVES (0.307 -> 0.149-0.175) for any multiplier tried, 3x-100x -- the
  best single-lever improvement of this whole investigation.
- **But it comes at a real, exactly-as-you-flagged cost**: the capacity-fade
  fit (model vs. real C/20 RPT SoH) degrades from an excellent 0.79pp gap
  to ~11-13pp once diffusivity is raised. The SEI/porosity-floor parameters
  were tuned assuming the old (too-low) diffusivity; raising it shifts the
  capacity-fade trajectory enough that those would need re-tuning together
  to recover both fits at once. This wasn't done tonight -- it's a coupled
  re-fit, not a one-line change, and is the natural next step.
- **Diffusivity alone doesn't fully close the gap**: even at 100x, the
  model's collapse ratio plateaus at ~3.4x, still well above the real
  1.27x. Returns diminish sharply past 10x. This means there's a second,
  diffusivity-independent mechanism also contributing to the residual
  collapse (candidates: SEI film resistance growth, contact resistance,
  electrolyte transport limitation, or the isolation-gate mechanism itself)
  -- worth a follow-up investigation, not chased tonight.
- **Item 30 (isolated-Si residual) does NOT help on top of the diffusivity
  fix -- it makes the expansion fit monotonically WORSE** as the residual
  fraction increases (0.149 -> 0.189 -> 0.243 for frac 0/0.5/0.9 on top of
  x3). This is a clear negative result, contrary to the hoped-for outcome:
  once the dominant sto-collapse factor is addressed by diffusivity, adding
  back a flat fraction of "isolated material still expanding" overshoots
  the real signal rather than complementing it. Recommend leaving item 30
  off going forward unless a different (non-flat-fraction) formulation is
  tried.
- Item 30 alone (no diffusivity change) reproduces its known earlier-session
  behaviour: some improvement over raw baseline (0.307 -> 0.228) with zero
  capacity-fade side-effect (gap stays 0.78pp), since it only touches the
  reported thickness output, not the electrochemistry.

**Bottom line recommendation**: adopt the Si diffusivity fix (x3-x10 range)
as the primary lever, but treat it as needing a joint re-tune of
K_SEI_MULT / porosity floor / isolation kinetics against the real
capacity-fade curve before calling this settled -- raising diffusivity
alone, without that re-tune, trades a capacity-fade regression for the
expansion-fit gain. Don't pursue item 30 further as configured.

## 4. Follow-up sweeps (2026-09-19, later same session): K_SEI/floor retune,
Si exchange-current density, and combinations

Added C/3 ageing-leg Si sto-window tracking (`c3_collapse_ratio`) to the
sweep worker, since it's the metric that actually drives the fitted
expansion signal (C/20 RPT is only an apples-to-apples DIAGNOSTIC, not what
`score_expansion_shape` is scored against). This sweep's own internally-
consistent baseline gives **c3_collapse_ratio=299x** at diffusivity x1 --
notably higher than the 17.7x quoted earlier from a one-off diagnostic
script that left `GR_REDIRECT_TO_LAM` at a different default; treat this
sweep's numbers as authoritative going forward.

| lever | best expansion RMSE | C/3 collapse | capacity-fade gap (pp) |
|---|---|---|---|
| diffusivity x3 (still the best overall) | **0.149** | n/a (not tracked yet at the time) | 11.25 |
| K_SEI retune on top of x3 (1.33x-2.67x) | 0.289-0.532 (worse, monotonically) | n/a | 1.4-10.0 |
| **floor retune on top of x3/x10 (combined best)** | 0.367-0.377 (**worse than baseline!**) | **172-177x (worse than ever!)** | 0.4-0.44 (looks great) |
| exchange-current alone (x3-x30) | 0.307-0.357 (flat/slightly worse) | 18.5-18.9x (16x better than x1) | 0.6-2.0 |
| exchange-current x3-10 + diffusivity x3 | 0.243-0.250 (worse than diffusivity alone) | 16.1-16.2x | 18.8-18.9 (worst of all) |

**Key findings:**

- **The K_SEI+floor "combined" retune is a trap.** It was picked by looking
  good on the C/20 RPT collapse ratio (1.1-1.2x, nearly matching real's
  1.27x) AND the capacity-fade gap (0.4-0.7pp, excellent) -- but the C/3
  ageing-leg collapse (what the fit is ACTUALLY scored against) got roughly
  10x WORSE (172-177x vs. 299x baseline... actually comparable/slightly
  better than x1's 299x but far worse than diffusivity x3's implied
  improvement), and expansion RMSE ended up worse than even the original
  pre-fix baseline. Lowering the porosity floor enough to fix the SLOW-rate
  metrics apparently lets isolation progress much further overall, which at
  FAST C/3 current density is catastrophic even though the slow C/20
  diagnostic doesn't reveal it. This is exactly why tracking C/3 alongside
  C/20 RPT matters -- judging by RPT/capacity-fade alone would have wrongly
  flagged this as a good fit.
- **Si exchange-current density is a real, independent lever** -- it
  dramatically fixes the DEEPEST end-of-life C/3 floor (299x -> 18.7x
  collapse ratio, at only x3) with almost no capacity-fade cost. But this
  barely moves the aggregate expansion-shape RMSE, because that metric
  averages error across the real dataset's whole EFC range, and the
  exchange-current fix mainly rescues the last few life-points, not the
  overall shape. Returns plateau by x3, same pattern as diffusivity.
- **Combining exchange-current with diffusivity does NOT synergise** --
  worse expansion RMSE (0.243-0.250) than diffusivity alone (0.149), and
  the WORST capacity-fade gap of everything tried (18.8-18.9pp). Stacking
  polarisation-relief levers compounds the capacity-fade side effect
  without compounding the expansion-fit benefit.
- **Diffusivity x3 alone remains the single best result of the whole
  investigation.** None of K_SEI retune, floor retune, exchange-current
  alone, or exchange-current+diffusivity beat it on expansion RMSE. The
  11.25pp capacity-fade cost it carries is real and still unresolved --
  none of tonight's single-parameter joint-retune attempts fixed it without
  a worse trade-off elsewhere. A genuine fix likely needs a proper multi-
  parameter joint optimisation (K_SEI, floor, diffusivity, exchange-current
  together), not one-at-a-time grid probes.

See `all_levers_summary.png` (every config, expansion RMSE vs. capacity-
fade gap) and `c3_collapse_vs_fit_landscape.png` (C/3 collapse ratio vs.
expansion RMSE, showing the floor-retune trap clearly).

## 5. Reversible-expansion measurement-methodology check

Per request, checked whether the real "reversible expansion" data might be
measuring something other than a per-cycle peak-to-peak amplitude (which is
what `cycle_expansion_ptp_um` computes on the model side).

- The source paper (`Notes_Papers/Expansion_as_Precursor.pdf`, Zhiwen Wan,
  UM Battery Control Lab) is a short 6-slide conceptual deck introducing
  the pore-buffering hypothesis and the `k` mechanical-transfer parameter
  (matches this project's own `k` convention exactly) -- it doesn't specify
  the numerical extraction method itself.
- Traced the actual computation to
  `Pouch_Data/cell064_aging_data/process_cell064_aging_data.py`'s
  `process_reversible_expansion()`, which reads `segment_stitch_amp_um`
  from `cell064_segment_stitch_trend_preview.csv` (columns:
  `raw_amp_um -> qc_cleaned_amp_um -> segment_stitch_amp_um`). **This
  confirms it IS a per-cycle amplitude metric** -- "segment_stitch" refers
  to patching sensor jumps/resets in the raw VDF trace within each cycle,
  not a different measurement convention. No mismatch found: we are
  comparing the right thing.

## Files in this folder

- `extract_real_si_sto_window.py` / `real_si_gr_sto_window_by_rpt.csv` --
  real per-RPT Si/Gr stoichiometric window.
- `extract_model_si_sto_window_c20rpt.py` /
  `model_si_sto_window_c3_and_c20rpt.csv` -- model's own C/20 RPT and C/3
  ageing per-cycle Si sto window across full simulated life.
- `plot_dma_comparison.py` / `dma_si_sto_window_comparison.png` -- real vs.
  model comparison figure.
- `si_diffusivity_sweep_results.csv` (mirrored from
  `../degradation_test_matrix/sweeps/baseline/`) / `plot_diffusivity_sweep_summary.py`
  / `si_diffusivity_sweep_summary.png` -- the sweep above.
- `expansion_fit_baseline_x1.png` / `expansion_fit_best_x3.png` -- the
  standard expansion-fit plot for the worst and best configs, for direct
  visual comparison.
