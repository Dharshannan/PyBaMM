# CELL064 experimental aging data

Copied from `Si_Gr_Expansion_Precursor/Pouch_Data/cell064_aging_data/data/` (a
separate project/repo -- see that project's `process_cell064_aging_data.py` /
`reconstruct_capacity_fade_knee.py` for how each file was derived; not
duplicated here). Only 6 RPT checkpoints exist (RPT0-5, EFC 0.8-248), so
every series below is sparse -- treat shapes between points as interpolation,
not measurement.

- **`CELL064_capacity_fade.csv`** -- the primary fitting target. 4 real C/20
  RPT discharges (RPT1, 2, 4, 5) + 2 high-rate-neighbor substitutes (RPT0,
  3, whose own C/20 discharge failed QC) filled in and rate-averaged over
  their before/after events, flagged via `source`. RPT1 (EFC 50.6) is the
  100%-retention reference, matching `cell064_parameters.py`'s own BoL
  anchor (RPT001) -- not RPT0, since RPT0's discharge capacity comes from a
  different (high) rate and would bias the retention percentage.
- **`CELL064_capacity_fade_highrate_all.csv`** -- all 10 individual
  high-rate-neighbor events (denser, but a different, non-C/20 rate) --
  context only, not part of the rate-consistent fitting series above.
- **`CELL064_capacity_fade_knee_reconstruction.csv`** -- a hinge-regression
  reconstruction (NOT a measurement) through the 6 capacity points, with the
  knee EFC taken from the expansion peak below. Useful as a visual guide to
  the likely underlying shape, not a fit target in its own right.
- **`CELL064_reversible_expansion.csv`** -- per-cycle QC-cleaned reversible
  expansion amplitude, EFC 1.9-247 (much denser than capacity, real-time
  cycling data). Its smoothed peak (EFC ~152) is the knee-location anchor
  used throughout this project. Not yet used for degradation-rate fitting
  (only capacity/LLI/LAM so far, per the current task) -- reserved for
  `../../expansion_pattern_test/` once pore buffering is retuned for
  CELL064.
- **`CELL064_k_expansion_scale.csv`** -- expansion scale
  k = delta_rev,cell / delta_rev,particle, same use as above (reserved for
  the expansion-fitting stage).
- **`CELL064_LLI.csv`** -- `charge_LLI`/`joule_LLI` are the source package's
  own fitted REMAINING cyclable-lithium-inventory values (a capacity-like
  state, decreasing from BOL), not cumulative LLI. `*_loss` columns
  (`reference RPT1 value - value`) are the increasing-from-zero quantity
  "LLI" conventionally means, derived here for direct comparison with
  PyBaMM's own LLI output. Only 4 points (RPT1, 2, 4, 5); no RPT0/RPT3/
  high-rate-resolution LLI exists in the source package.
- **`CELL064_LAM_derived.csv`** -- **NOT independently re-fit LAM** -- simple
  capacity-retention ratios (`1 - Cn_Gr(t)/Cn_Gr(RPT1)` etc.) computed from
  the source package's own already-fitted `discharge_Cn_Si/Cn_Gr/Cp`
  columns. Per the task instruction, prefer running the project's own DMA
  method (`../../test_model/degradation_test_matrix/test_DMA/`) on this
  model's simulated RPT discharge curves over trusting these numbers
  directly when judging a candidate degradation-rate fit -- these retention
  ratios carry whatever eSOH-fit error the source package's own Cn_Si/Cn_Gr/
  Cp identification has, are only 4 points, and (per that package's
  `process_lam()` docstring) are explicitly flagged as a simplification, not
  a rigorous DMA. Useful for order-of-magnitude sanity-checking only:
  LAM_NE_silicon dominates the trend (0% at EFC100 -> 58.5% at EFC222 ->
  79.6% at EFC248), LAM_NE_graphite stays modest (<10% throughout), LAM_PE
  is small and roughly flat after an early ~8% step -- i.e. a Si-fatigue-
  dominated knee, consistent with CELL064's high silicon anode fraction
  (48% of anode capacity, see `cell064_parameters.py`'s docstring).
- **`CELL064_LAM_derived_highrate.csv`**, **`CELL064_k_expansion_scale_highrate_all.csv`**
  -- denser high-rate-neighbor-resolution versions of the above, referenced
  to their own first point (different rate/fit conditions) -- context only.
