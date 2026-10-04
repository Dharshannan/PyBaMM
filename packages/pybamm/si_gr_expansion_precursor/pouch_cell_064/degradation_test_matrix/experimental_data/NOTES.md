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

# CELL017 (45 degC) experimental set — `CELL017_45C/`

CELL017 is at 45 degC and 103 kPa / 15 psi, the same pressure as CELL064. It
was built from the IC_PYBAMM_FIG1 data package (v1, 2026-09-27) by
`build_cell017_45C_data.py`, in the same file formats as the CELL064 files
above:
- `CELL017_capacity_fade.csv`
- `CELL017_LLI.csv`
- `CELL017_LAM_derived.csv`
- `CELL017_k_expansion_scale.csv`
- `CELL017_reversible_expansion.csv`
- `low_rate_c20/` (RPT0-7 C/20 discharges)

The high_temp_45C fork overlays and scores all of them at 318.15 K
(`load_cell017()`). Differences forced by what the package contains:

- **Reference RPT.** RPT0 (EFC 6.6) is a 25 degC pre-characterisation, so
  RPT1 (EFC 8.72) is the reference and zero, like CELL064's RPT1. All 8 RPTs
  have valid C/20 fits; there are no high-rate substitutes.
- **LLI.** This is the package's common `LLI_Ah` (x_start*Cn_total + y_start*Cp,
  from the discharge fit). CELL064's `charge_LLI` came from a charge-curve fit
  that isn't shipped; the two definitions agree within ~1-5% on CELL064.
- **Reversible expansion.** Per-cycle discharge-swing amplitude from the
  lifetime pickle, via `extract_lifetime_reversible_expansion.py`. Validated on
  CELL064 against its own pipeline output: median -2.9%, correlation 0.92.
- **Knee.** The manifest capacity knee is at EFC 318.9 raw (~310 since RPT1);
  the expansion transition is at 297.1.
- `CELL017_45C_capacity_fade.csv` (top level) is the user's earlier quick
  capacity/2.5 Ah overlay, superseded by the folder above.

## Partial-SoC-window cells: CELL009 (15-95%) and CELL026 (20-80%), 25 degC, 103 kPa

Built by `build_partial_soc_data.py` from the same package, into
`CELL009_25C_15-95/` and `CELL026_25C_20-80/`, in the CELL017 formats plus
`CELLxxx_cycling_protocol.json`.
- **Which cells.** The manifest's partial-window cells are CELL009 and CELL054
  (15-95%) and CELL026 (20-80%). CELL025 is a 45 degC / 103 kPa cell, a twin
  of CELL017. CELL054 is skipped because its data is incomplete.
- **Reference RPT.** RPT1 is the reference and zero. k at RPT1 is 0.705
  (CELL009) and 0.708 (CELL026), so F0 = 0.70, the same as CELL064.
- **Capacity-only RPTs.** CELL009 RPT2/4 and CELL026 RPT2/3/5 have no eSoH
  fit: they are marked invalid only for missing raw expansion. Their C/20
  capacity (QC `v_q_span_Ah`) is kept in `capacity_fade` (`esoh_fit` =
  False), but they have no LLI/LAM/k rows and no curve files.
- **Ageing control** (detected from the lifetime data, constant over life):
  - **CELL009:** C/3 CC-CV charge to 4.150 V (CV to ~14 mA), then C/3 CC
    discharge to 3.140 V. Voltage-limited at both ends; the Ah per discharge
    falls from 1.94 to ~0.94 Ah with ageing.
  - **CELL026:** C/3 CC-CV charge to 3.969 V (CV to ~14 mA), then a C/3
    discharge of a fixed 1.500 Ah (coulomb counting, 60% of 2.5 Ah) with a
    2.5 V safety floor. The discharge end voltage drifts from 3.18 V down to
    the 2.5 V floor after the knee, and the Ah per cycle then drops (0.87,
    then 0.68 Ah).
- **Knees.** The manifest capacity knees are at EFC 298.2 (CELL009) and
  301.4 (CELL026), about twice CELL064's.
- **eSoH noise.** CELL026's graphite LAM is negative (−8 to −11%) at RPT0/4/6.
  Its RPT1 Cn_Gr (1.368 Ah) is low compared with the others (1.48-1.52), so
  treat CELL026's Gr LAM as noisy.
- **EFC basis.** Every package EFC, including CELL064's and CELL017's, is
  throughput / (2 x 2.5 Ah). The model's `efc_from_throughput` divides by
  2 x 2.5947 Ah, so at equal EFC the model has ~3.8% more throughput. That
  convention is used unchanged throughout the fits so far.
