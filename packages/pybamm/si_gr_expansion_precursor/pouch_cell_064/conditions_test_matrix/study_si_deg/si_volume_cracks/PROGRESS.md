# si_volume_cracks: progress log (branch test_volume_strain_cracking)

Working log for the strain-driven cracking study. The plan is in
`../PROPOSAL_strain_driven_cracking.md`. Runs use `run.sh` (aligned RPTs,
outputs in this folder).

## Stage 1: cathode LAM saturation (power-law damping on the positive)

- **New fork knob:** `C64_POS_STRESS_LAM_DAMP_EXP=n` sets the positive
  electrode's `"stress-driven LAM damping"` to `"power"`, i.e. stress ×
  (eps_s/eps_s0)^n. It works with the existing `C64_POS_LAM_MULT`
  (β multiplier).
- **Saturation estimate.** The rate goes as (1−L)^(2n) for m_p = 2. Matching
  the v3 25 °C LAM_PE (~6.5% at EFC 200) needs:

  | n | β multiplier |
  |---|---|
  | 25 | ~8 |
  | 50 | ~120 |
  | 75 | ~2300 |

  Larger n gives a flatter plateau.

## Stage 2 (built early, while stage 1 runs): strain-fatigue cracking in core

- **Core option** `"particle cracking growth"`: `"Paris"` (default) /
  `"Paris + strain fatigue"` (B1) / `"Paris + strain fatigue (contraction)"`
  (B2). See CHANGES.md item 15.
- **Fork knobs:** `C64_SI_STRAIN_CRACK_KV`, `C64_SI_STRAIN_CRACK_MODE`
  (`contraction` default | `abs`), `C64_SI_STRAIN_CRACK_EA`.
- **Two fixes found in smoke tests:**
  1. Literature `t_change` is ΔV/V0, so the volume ratio is 1 + t_change.
  2. sto is clamped to [1e-6, 1]. At exactly 0, the Jacobian of s^n_eff is
     NaN and IDA stalls at the end of the formation discharge.
- **Calibration** (v3, k_V = 1e-6 smoke tests, B2):

  | Window | Cumulative Si strain E per EFC | Paris crack growth |
  |---|---|---|
  | Full window, 25 °C | ~1.50 (≈ ln 4) | 34 nm by EFC 50 |
  | 15–95% | ~1.31 (0.88×) | 0.5 nm by EFC 80 |

  So route A switches off in the partial window while route B keeps ~90%.
- **k_V estimate.** The CELL009 knee needs l/l0 ≈ 2–2.4 at EFC ~290 (E ≈ 380),
  giving k_V(B2) ≈ 0.8/380 ≈ 2e-3. Stage-3 grid: 1.5e-3 to 3e-3.

## Stage 1 results: full windows (aligned RPTs)

| Run | Cathode LAM at 25 °C, EFC 49 / 171 / 197 (real 7.9 / 6.2 / 6.9) | 25 °C mean SoH gap | Cathode LAM at 45 °C, EFC 51 → 400 (real 4.3, 3.3, 4.3, 3.8, 5.1, 3.0) | 45 °C mean SoH gap |
|---|---|---|---|---|
| v3 reference (linear) | 2.2 / 5.8 / 6.5 | 3.3pp | 0.8 → 5.3 | 0.8pp |
| S1a n=25, β×8 | 4.6 / 6.5 / 6.7 | 2.8pp | 2.8 → 6.2 | 1.2pp |
| S1b n=50, β×120 | 5.6 / 6.5 / 6.6 | **2.2pp** | 4.5 → 6.4 | 1.4pp |
| S1c n=75, β×2300 | 5.9 / 6.5 / 6.6 | 2.3pp | 5.2 → 6.4 | 1.5pp |

- **25 °C:** saturation front-loads the cathode LAM, hits the same end LAM,
  and improves SoH. S1b is the best.
- **45 °C:** every saturating law plateaus at ~6.2–6.4%, against a real ~4%.
  The plateau depends only logarithmically on the stress rate, and running
  2× the EFC cancels the lower 45 °C rate.
- **Real cathode LAM ordering:** 25 °C full ~6.9% > 25 °C 15–95% ~4.5% ≈
  45 °C ~4%. That ordering follows the cathode stress. A single β then needs
  a ~4× lower 45 °C rate.
- **Testing:** an isothermal cathode-LAM activation energy
  (`C64_POS_LAM_EAC`, exactly 1 at 25 °C) of −40 and −55 kJ/mol on S1b at
  45 °C (runs S1e).

## Stage 1 conclusion: cathode LAM, adopted `stage1_v4a.env`

| Variant | 25 °C end cathode LAM (real 6.9) | 45 °C cathode LAM (real 3.0–5.1) | CELL009 at EFC 402 (real 4.8) | Mean SoH gap 25 / 45 / CELL009 |
|---|---|---|---|---|
| v3 linear | 6.5 | 0.8 → 5.3 | 18.2 (needed the ×0.3 patch) | 3.3 / 0.8 / – |
| **n=50, β×120, Ea −55 kJ/mol** | **6.6** | **3.2 → 5.0** | 7.5 | **2.2 / 1.1 / 2.2** |
| n=50, β×30, no Ea | 5.3 | 3.2 → 5.1 | 6.1 | 2.5 / 1.1 / 2.2 |
| n=50, β×120, Ea −40 kJ/mol | (6.6) | 3.5 → 5.4 | – | – / 1.2 / – |

- **Adopted: `stage1_v4a.env`** (shared 25/45 °C) and
  **`stage1_v4a_CELL009.env`** (adds crack 1.1e4 and width 1e-3). The CELL009
  cathode ×0.3 patch is gone.
- **Known gap.** CELL009 cathode LAM still overshoots (+2.7pp). A
  damping-type saturation plateaus at ~the same level in every condition;
  the real plateau tracks the cathode stress. A stress-threshold
  (weakest-particle) law is proposed for review in PROPOSAL §10.1. It is not
  implemented, since it is a new mechanism.
- **The −55 kJ/mol cathode activation energy** is a shared parameter of a
  temperature-dependent process, but its sign (slower cathode fracture at
  high T) needs your judgement. The alternative, β×30 with no Ea, matches
  45 °C equally well but undershoots 25 °C (5.3%).
- **Cleanup:** non-final stage-1 runs deleted. Kept S1b (25 °C), S1e
  −55 kJ/mol (45 °C), S1n n50/β120 (CELL009).

## Stage 3: route B (strain-fatigue cracking) sweeps

### S3a: CELL009, B2, shared Paris 1e3, width 1e-3, cathode n50/β×30 (pre-fix code)

| k_V | SoH gap at EFC 79 / 145 / 184 / 260 / 329 / 382 / 402 (pp) | Mean SoH gap | Voltage RMSE |
|---|---|---|---|
| 1.5e-3 | … / +13.6 at 329 | 8.0pp | – |
| 2e-3 | −1.6 / −3.2 / −2.6 / −1.8 / +4.3 / +4.1 / +5.8 | 3.4pp | 0.077 V |
| 3e-3 | −1.8 / −3.6 / −3.3 / −4.8 / −19.4 / −7.1 / −2.6 | 6.1pp | 0.119 V |

- **Knee timing scales steeply with k_V.** Best k_V is ~2.3e-3.
- **Route B dominates crack growth:** Paris adds only ~2 nm by EFC 260 in the
  15–95% window.

### Bug found and fixed: the expansion hump disappeared with route B on

- **Cause.** In the first version, the strain term looked up the particle rhs
  inside `get_coupled_variables`, after the mechanics had already written its
  thickness change. The KeyError deferred the submodel; its retry ran after
  pore buffering and overwrote the buffered thickness with the unbuffered
  one.
- **Effect.** The CELL009 swing was 11.9 µm/cycle instead of 7.4. Expansion
  output only: the k-ratio, LAM and SoH curves matched S1n to within noise.
- **Fix.** The strain term is now added in `CrackPropagation.set_rhs`, plus a
  pore-buffering deferral guard in `reaction_driven_porosity.py` (CHANGES
  item 15). The smoke test gives a 7.3 µm swing with route B on and off.
- **Invalid outputs.** S3a and S3c **expansion** panels predate the fix and
  are invalid. Their degradation results are valid.

### S3d: CELL009 rerun on the fixed code with the shared v4a cathode law

Setup: `stage3_base.env`, width 1e-3, shared Paris 1e3, B2.

| k_V | SoH gap at EFC 329 / 382 / 402 (pp) | Mean SoH gap | Voltage RMSE | Expansion RMSE |
|---|---|---|---|---|
| 2e-3 | +4.6 / +4.1 / +5.7 | 3.5pp | 0.080 V | 0.100 |
| **2.3e-3** | −4.4 / −0.7 / +2.2 | **2.6pp** | 0.095 V | 0.177 |
| 2.6e-3 | −12.1 / −4.0 / −0.3 | 4.0pp | 0.110 V | 0.254 |

- **Pre-knee offset.** The −2 to −3.5pp before the knee is mostly the
  cathode-LAM overshoot (7.5% vs 4.8%; stage-1 known gap).
- **Comparison.** S1n, which needed the ×11 crack patch, gave 2.2pp.

### S3c / S3e: 25 °C (CELL064) with route B

| Run | SoH gap at EFC 49 / 171 / 197 (pp) | Mean SoH gap | 91% crossing (real ~101) |
|---|---|---|---|
| S3c k_V 2e-3, Paris 1e3 (pre-fix) | −2.7 / −5.6 / −2.8 | 3.7pp | 98 |
| S3e k_V 2.3e-3, Paris 0.7e3 | −2.3 / +1.7 / +0.8 | 1.6pp | 114 |
| **S3e k_V 2.3e-3, Paris 0.75e3** | −2.4 / −0.4 / −0.2 | **1.0pp** | 110 |
| S3e k_V 2.3e-3, Paris 0.85e3 | −2.5 / −3.2 / −1.6 | 2.5pp | 104 |

Route B is ~20% of crack growth at 25 °C (45 vs 194 nm by EFC 198), so
Paris is lowered ×0.75 to compensate.

### S3c / S3f: 45 °C (CELL017) needs route B slowed when hot

| Run | Mean SoH gap | Crack growth by EFC 400, Paris / strain (nm) |
|---|---|---|
| S3c k_V 2e-3, no E_V, `SI_CRACK_EAC` 0 | 21.9pp | 120 / 102 |
| S3c k_V 2e-3, no E_V, `SI_CRACK_EAC` 90k | 14.2pp | 13 / 76 |
| S3f k_V 2.3e-3, E_V +80k | 3.5pp | – |
| S3f k_V 2.3e-3, E_V +100k | 2.3pp | 3.8 / 2.7 |
| **S3f k_V 2.3e-3, E_V +120k** | **1.6pp** (all RPTs −1.0 to −2.2) | 3.4 / 1.6 |

E_V values are in the `SI_CRACK_EAC` convention (positive = slower when hot).

- **Why the 45 °C knee is so sensitive.** It is a pore-clogging knee (S1e:
  l/l0 only 1.15 by EFC 400). Any extra crack SEI drives the porosity to its
  floor early.
- **Consequence.** Route B must be essentially off at 45 °C (k_V ÷21), the
  same direction and physics as route A's +90k: lithiated Si is more ductile
  when hot.
- **Open question.** Whether ÷21 is physically plausible.

## Stage 3 shared set: `stage3_v5.env`

- **Contents:** v3 + v4a cathode law + Paris ×0.75e3 + B2 k_V 2.3e-3 with
  E_V +120 kJ/mol.
- **Per-condition values:** only T, F0 and WIDTH (CELL009 1e-3). The CELL009
  ×11 crack patch is gone.
- **25 °C:** identical to S3e Paris 0.75e3, since E_V has no effect at
  25 °C.
- **Running:** V5_45degC, V5_CELL009_15-95, and the CELL026 blind test at
  width 1e-3 and 5e-3.

### Sign convention for E_V (changed 2026-10-05)

- **New convention.** `C64_SI_STRAIN_CRACK_EA` now uses the SAME convention
  as `SI_CRACK_EAC` (Ai2020), exp[Ea/R (1/T − 1/298.15)], so positive means
  slower when hot (lithiated Si more ductile). Route A's +90 kJ/mol is kept.
- **Old S3f tags.** The S3f runs (tagged `EV-80000` / `-100000` / `-120000`)
  were launched with the old standard-sign code, so they are E_V = +80 /
  +100 / +120 kJ/mol in the new convention.
- **What the 45 °C runs show.** Without temperature dependence (S3c), route B
  knees ~100 EFC early at 45 °C. Route B needs a positive E_V, the same
  direction as route A.

### v5 results (aligned RPTs)

| Condition | Mean SoH gap | Gap after knee (pp) | Voltage RMSE | Expansion RMSE | Stage-1 / v3 reference |
|---|---|---|---|---|---|
| 25 °C CELL064 | **1.0pp** | −0.4 / −0.2 at EFC 171 / 197 | 0.098 V | 0.109 | S1b 2.2pp |
| 45 °C CELL017 | **1.1pp** | −0.6 / −1.1 at EFC 334 / 400 | 0.067 V | 0.120 | S1e 1.1pp |
| 15–95% CELL009 | **2.3pp** | −2.8 / +0.1 / +2.8 at EFC 329 / 382 / 402 | 0.092 V | 0.164 | S1n 2.2pp, which needed the ×11 crack patch |
| 20–80% CELL026 (blind, width 1e-3) | 5.9pp | −9.2 / −11.2 / +5.5 at EFC 293 / 324 / 383 | 0.111 V | 0.441 | – |
| 20–80% CELL026 (blind, width 5e-3) | 6.7pp | −12.5 / −13.6 / +5.4 | 0.113 V | 0.509 | – |

- **25 °C row.** V5_25degC, rerun from `stage3_v5.env`. It is identical to
  S3e Paris 0.75e3, since E_V has no effect at 25 °C.
- **45 °C with E_V +150k (V5b).** Same mean gap (1.1pp); after the knee
  +0.9 / +0.6pp vs −0.6 / −1.1pp at 120k. The 45 °C data can't tell 120k
  from 150k, so v5 keeps 120k.
- **Cleanup (2026-10-05).** S3a–S3f, the smoke runs and V5b were deleted.
  Kept: stage-1 finals (S1b, S1e, S1n) and V5_* (25 °C, 45 °C, CELL009,
  CELL026 ×2).

**Why the CELL026 blind test misses:**

- **Route B over-predicts the 20–80% strain per EFC.**

  | Cell | Cumulative Si contraction strain | Strain per EFC | l/l0 |
  |---|---|---|---|
  | CELL026 (20–80%) | 524 by EFC 293 | 1.79 | 3.2 at 293 |
  | CELL009 (15–95%) | 469 by EFC 330 | 1.42 | 2.8 at 330 |

  - The per-cycle Si swing is about the same (~1.0 vs ~1.1). Removing the top
    of the window removes mostly graphite capacity.
  - But 20–80% packs ~1.4× more cycles into each EFC, so its Si cracks grow
    faster.
  - Real data: the two cells knee at about the same EFC (~293 vs ~290).
  - So pure cumulative log-strain is not the controlling variable across
    window widths. This is the discriminating test flagged in PROPOSAL §4.
- **Pre-knee drift.** Model LLI runs ahead (7.9% vs 3.9% at EFC 167), and so
  does the cathode LAM (6.4–7.0% vs ~0–4%).
- **Possible fixes, for review; not implemented, since they change the law:**
  1. **Coffin–Manson amplitude exponent p > 1** (PROPOSAL §3.4). It
     penalises large swings, but per-cycle swings are nearly equal here, so
     it would need p ≈ 3. Weak lever.
  2. **Throughput driver**, D ∝ |j_Si| (PROPOSAL §4 caveat). It weights Si
     use by charge rather than by log volume. The log weighting
     over-emphasises low-sto Si, where the 20–80% window sits.
  3. **Accept CELL026's window-specific behaviour.** Its knee is sharper,
     which a pore-clogging threshold gives more readily than gradual
     cracking.

## v6: earlier 25 °C knee (user request, 2026-10-05) → `stage3_v6.env` (CURRENT)

**Goal:** the user wants the 25 °C knee earlier, accepting an undershoot at
the last two RPTs, and the 20–80% knee later if possible.

**Paris rate:**

| Paris | 25 °C: 91% crossing (real ~101) / knee-fit EFC / mean gap | 45 °C, E_V 150k: mean gap (after knee) | CELL009: mean gap (at EFC 329) |
|---|---|---|---|
| 0.75e3 (v5) | 110 / 141 / 1.0pp | 1.1pp (+0.9 / +0.6) | 2.3pp (−2.8) |
| **1.0e3 (v6)** | **97 / 122 / 3.9pp** (−6.0 / −3.0 at 171 / 197) | **1.0pp** (−0.2 / −0.6) | **2.6pp** (−4.4) |
| 1.15e3 | 90 / 114 / 4.8pp (−7.8 / −3.7) | 1.2pp (−0.9 / −1.3) | 2.7pp (−5.2) |

**CELL026 width** (Paris 1.0e3):

| Width | Mean gap | Gap at EFC 293 / 324 (pp) |
|---|---|---|
| 1e-3 (v5, Paris 0.75e3) | 5.9pp | −9.2 / −11.2 |
| 5e-4 | 6.0pp | −9.2 / −11.9 |
| **3e-4** | **5.6pp** | −7.9 / −10.6 |

Width is a weak lever for the 20–80% knee. The −2.5 to −5.9pp pre-knee
offset (lithium loss and cathode LAM running ahead) and the route-B strain
per EFC (§ v5) are what's left.

**v6 = v5 with Paris back to v3's 1e3 and E_V +150k:**

| Condition | Mean SoH gap | Voltage RMSE | Expansion RMSE |
|---|---|---|---|
| 25 °C | 3.9pp (knee on time, late RPTs undershoot) | 0.105 V | 0.142 |
| 45 °C | 1.0pp | 0.066 V | 0.118 |
| CELL009 | 2.6pp | 0.095 V | 0.177 |
| CELL026 (width 3e-4) | 5.6pp | 0.109 V | 0.398 |

Kept runs: V6_25degC_crack1.0e3, V6_45degC_crack1.0e3_EV150k,
V6_CELL009_crack1.0e3 and V6_CELL026_crack1.0e3_w3e-4, plus the stage-1
finals. The V5 runs and the other V6 variants were deleted; v5's plots are
in commit df68b5301.

## Protocol-basis fix and final plots (2026-10-08)

- **Currents.** The real tests use C-rates of the 2.5 Ah nominal capacity:
  ageing 0.833 A, C/20 RPT 0.125 A. The model used its own 2.5947 Ah, so
  its currents were ~3.8% high. Fork `RATE_BASIS_AH` = 2.5 now sets
  "Nominal cell capacity [A.h]".
- **CELL026 discharge.** The fixed 1.5 Ah is now 6480 s at 0.833 A (was
  6243 s at 0.865 A). Still exactly 1.5 Ah and never updated, as in the
  test.
- **EFC basis.** Fork `EFC_BASIS_AH` = 2.5 now matches the data's
  throughput / (2 × 2.5 Ah); the model's EFC was ~3.8% low at equal
  throughput. The aligned RPTs therefore land at ~3.8% less model
  throughput than before.
- **Reruns.** v6 (`stage3_v6.env`) rerun for all 4 conditions into
  `final_plots/`. Plot names are `<condition>_<plot>.png` (25deg_full_soc,
  45deg_full_soc, 25deg_15-95_soc, 25deg_20-80_soc) via
  `C64_FIG_NAMES=condition`. Figure titles carry the cell
  (`CELL_LABEL`, set by rpt_soc_plots.py).

**Final results** (v6, corrected currents and EFC basis; `final_plots/`, with
data in `final_data/`):

| Condition | Mean SoH gap (before fix) | Voltage RMSE | Expansion RMSE | Notes |
|---|---|---|---|---|
| 25deg_full_soc (CELL064) | 3.3pp (3.9) | 0.104 V | 0.124 | 91% crossing at EFC 98 (real ~101) |
| 45deg_full_soc (CELL017) | 1.2pp (1.0) | 0.059 V | 0.112 | after knee +0.8 / +1.5pp |
| 25deg_15-95_soc (CELL009) | 2.4pp (2.6) | 0.084 V | 0.133 | after knee +0.4 / +2.1 / +4.3pp |
| 25deg_20-80_soc (CELL026) | 4.5pp (5.6) | 0.105 V | 0.334 | knee −4.9 / −5.8pp at EFC 293 / 324 (were −7.9 / −10.6) |

**`final_plots_v2/` and `final_data_v2/` (2026-10-09).** The same v6 set
and the same titles and names, but with the PRE-fix basis
(`C64_RATE_BASIS_AH` = `C64_EFC_BASIS_AH` = 2.5947, so C/3 = 0.865 A).
They reproduce the earlier v6 scores exactly:

| Condition | Mean SoH gap (v2, pre-fix) | Mean SoH gap (fixed) |
|---|---|---|
| 25 °C | 3.9pp | 3.3pp |
| 45 °C | 1.0pp | 1.2pp |
| 15–95% | 2.6pp | 2.4pp |
| 20–80% | 5.6pp | 4.5pp |

The corrected runs' logs are now in `final_data/`.

## Isolation-gate dose form (branch `test_isolation_reform`, 2026-10-13)

**Setup:**
- Option `"porosity isolation gate": "dose"` (CHANGES item 16): dD/dt = g*,
  LAM = k_iso·D·a·j_SEI/F with k_iso = Y(T)/τ, and no (1 − g) cut on SEI
  growth.
- Fork knob: `C64_ISO_GATE_FORM=dose`.
- Same v6 set and corrected 2.5 Ah basis as `final_plots/`.
- Outputs in `final_plots_v3/` and `final_data_v3/`.

**Results:**

| Condition | Mean SoH gap, dose (v3) | Mean SoH gap, lagged (final) | Largest per-RPT difference |
|---|---|---|---|
| 25 °C full | 3.3pp | 3.3pp | 0.1pp |
| 45 °C full | 1.2pp | 1.2pp | 0.1pp |
| 15–95% | 2.4pp | 2.4pp | 0.1pp |
| 20–80% | 4.5pp | 4.5pp | 0.1pp |

- Voltage and expansion RMSE agree to ≤ 0.0004.
- This confirms MODEL_EQUATIONS §1.3 and §1.8: the lagged gate is an
  integrator, and the SEI redirect is numerically negligible. The dose form
  is the same model with one constant fewer (k_iso,0 = Y0/τ).

**`final_plots_v4/` and `final_data_v4/`.** Dose-form gate (as v3) with the
OLD 2.5947 Ah current/EFC basis (as v2). Identical to v2 within 0.1pp:

| Condition | Mean SoH gap, v4 | Mean SoH gap, v2 |
|---|---|---|
| 25 °C | 3.9pp | 3.9pp |
| 45 °C | 1.0pp | 1.0pp |
| 15–95% | 2.6pp | 2.6pp |
| 20–80% | 5.6pp | 5.6pp |

**The four final sets:**

| Folder | Gate | Current basis |
|---|---|---|
| `final_plots` | lagged | fixed |
| `final_plots_v2` | lagged | old |
| `final_plots_v3` | dose | fixed |
| `final_plots_v4` | dose | old |
