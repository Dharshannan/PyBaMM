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
