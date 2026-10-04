# study_si_deg — cracking-dominated Si SEI growth, 25 degC refit

**Status: in progress, started 2026-09-28. Current 25 °C defaults = `25degC_baseline1`; cracking-driven 25/45 °C path under development (see the last three sections).** This is a deliberately
ISOLATED investigation — nothing here touches the old baseline, new
baseline, or high_temp_45C/narrow_voltage_window forks. Everything lives in
this folder.

## Hypothesis

The 45 degC knee-delay mechanism work (see `../CONDITIONS.md`'s
"Alternative theory" section) proposed that higher temperature reduces Si
particle cracking (lower overpotential + faster diffusivity → smaller
concentration gradients → less driving stress for crack propagation), which
reduces "SEI on cracks" growth, which reduces total porosity closure —
explaining the real high-T knee delay without needing an exotic negative
apparent activation energy on bulk SEI.

This study asks the next question: **if cracking-driven SEI growth is the
dominant physical mechanism, does it also explain the 25 degC baseline
itself?** I.e., instead of treating "ec reaction limited" bulk SEI as the
main knee-timing driver (as the current new-baseline recipe does), can a
**"solvent-diffusion limited" SEI law for Si, with most of the growth
coming from cracking (not bulk) pre-knee**, refit the real 25 degC data as
well as or better than the new baseline — with the post-knee
porosity-isolation LAM mechanism left unchanged?

## Plan

1. **Diagnostic**: quantify the bulk-vs-cracks SEI split in the CURRENT
   high_temp_45C fork's latest recipe (`check_crack_vs_bulk_sei.py`), as a
   reference point for "how much crack-area amplification is already
   happening" before deliberately pushing it further here.
2. **New model**: `cell064_degradation_fit_crack_dominated.py` — a full
   standalone fork of the shared baseline script, starting from the
   new-baseline's own parameter values (isolation ON, radius/diffusivity
   retune, etc.), with Si's SEI option switched to `"solvent-diffusion
   limited"` and `SI_CRACK_MULT` raised well above the "genuine" baseline
   value so cracks dominate pre-knee porosity closure over bulk SEI.
3. **Fit target**: match or beat the new baseline's own fit quality against
   real 25 degC CELL064 data, using the SAME metrics/scoring functions
   already in the shared script:
   - `score_rpt_gap` — SoH via C/20 RPT-vs-RPT gap (new-baseline's own
     benchmark not yet re-measured this session — see "Benchmark" below)
   - `score_voltage_shape` / `score_voltage_shape_soc` — RPT discharge V(Q)
     shape / depth-of-discharge
   - `score_expansion_shape` — reversible expansion RMSE
   - `score_k_gap` — expansion scale k
   - LAM (Si/Gr/PE split) and LLI — visual comparison against
     `load_experimental_lam()`/`load_experimental_lli()` (no dedicated
     score_* function exists for these; compare against
     `dma/FINDINGS_2026-09-20.md`'s crude LAM_Si estimate of 58-80% and the
     LLI trajectory shown in the new-baseline plots)
4. Tune ALL degradation-only parameters (D_sol/SOLVENT_MULT, SI_CRACK_MULT,
   SI_BETA_LAM_SEI, SI_CRIT_STRESS/SI_LAM_PROP, NEG_POROSITY_BOL/FLOOR,
   EXPONENT_MAX_SEI, F0/WIDTH, SI_REDIRECT_LAM_YIELD/SI_LAM_ISO_EXPONENT/
   SI_TAU_LAM_ISO) — NOT temperature, NOT voltage window, NOT anything
   related to the other two condition forks.

## Benchmark to match or beat

`CELL064_final_parameterisation.md` documents the OLD baseline's fit
(mean |gap|=1.9pp, knee EFC~108 vs anchor 101, expansion RMSE~0.10) — this
is NOT the right benchmark (confirmed 2026-09-27: that document predates
the new-baseline's isolation-mechanism + radius/diffusivity retune). The
new baseline's own quantitative scores have not been explicitly written
down anywhere found so far — **first step of the tuning work should be
running the new baseline once (via `degradation_test_matrix/
cell064_degradation_fit.py` with the exact `build_new_baseline_curves.py`
env-var overrides) and recording its score_rpt_gap/score_voltage_shape/
score_expansion_shape/score_k_gap output as the actual target**, rather
than assuming the old baseline's numbers apply.

## Starting configuration (2026-09-28, untested as a whole yet)

All new-baseline values inherited unchanged except:

| Parameter | New baseline | study_si_deg starting guess |
|---|---|---|
| SEI option (Si) | `"ec reaction limited"` | `"solvent-diffusion limited"` |
| `SOLVENT_MULT` (new) | n/a | `34.0` (anchor from the earlier, non-cracking-dominated knee-delay work: gave EFC~101.8 at SI_CRACK_MULT=1) |
| `SI_CRACK_MULT` | `0.1` | `20.0` (cranked up so crack area dominates — untested guess) |

**Expected issue, not yet resolved**: increasing crack area at the SAME
`D_sol` will likely close porosity FASTER (earlier knee) than the mult=34
anchor assumed, since that anchor was calibrated at `SI_CRACK_MULT=1`, not
20. `SOLVENT_MULT` will likely need to come DOWN from 34 to compensate and
land back on EFC~101 — this compensation has NOT been done yet; the
starting config above is an initial, deliberately-untuned first guess.

## Diagnostic result (2026-09-28)

Ran `check_crack_vs_bulk_sei.py` against the CURRENT high_temp_45C fork
(unmodified, its own latest recipe: SOLVENT_MULT=5, SI_CRACK_MULT=1.0 —
i.e. the "genuine", non-inflated crack rate). Result: **cracks already
contribute 65.6% of total porosity-closing SEI thickness vs bulk's 34.4%**,
stable across the whole run (roughness plateaus at ~2.91, i.e. crack area
~1.9x the base surface area). No artificial inflation of `SI_CRACK_MULT`
needed to get a cracking-dominated split — the "genuine" value already
gives one. This ratio should be ~independent of `SOLVENT_MULT` (both bulk
and cracks scale by the same D_sol; only the AREA ratio, set by roughness,
determines the split) — corrected the model's starting `SI_CRACK_MULT`
from an initial (wrong) guess of 20.0 down to 1.0 accordingly.

## iter00: first working config (2026-09-28)

Starting config: new-baseline values throughout, SEI=solvent-diffusion for
Si (SOLVENT_MULT=34, SI_CRACK_MULT=1.0) — i.e. literally the same
(SOLVENT_MULT, SI_CRACK_MULT) pair as the earlier knee-delay
investigation's own 25 degC anchor point, now checked for the FIRST time
against the full real-data metric suite (that investigation had real-data
overlay deliberately disabled).

**Bug found and fixed en route**: this fork's copied `GR_K_SEI_MULT`
default was left at the OLD baseline's `_K_SEI_MULT_BASE` (1.0e-3) instead
of the new baseline's `6e-6` — a ~167x too-fast graphite SEI rate, since
graphite still uses "ec reaction limited" and its k_sei is NOT inert. This
caused a completely broken first attempt (knee at EFC~28, SoH crashing to
31% by EFC~228). Fixed by hardcoding `GR_K_SEI_MULT`'s default to 6e-6
(and `SI_K_SEI_MULT` to 3e-4 for consistency, though the latter is inert
for Si now).

**Result after the fix** (`iter00_starting_config`):

| Metric | Result | Reference (old baseline, `CELL064_final_parameterisation.md`) |
|---|---|---|
| Knee EFC | **101.8** (err +0.4 vs anchor ~101) | ~108 (err +6.3) |
| RPT-vs-RPT mean \|gap\| | 8.2pp (RPT2 -3.3, RPT4 +9.0, RPT5 +12.2) | 1.9pp |
| Expansion shape RMSE | 0.175 | ~0.10 |
| Voltage shape RMSE (mean) | 0.079 V | not directly comparable (not reported in that doc) |
| k-gap mean \|gap\| | 0.306 (dominated by the same known k<=1 structural ceiling post-knee documented in the old-baseline doc — real k reaches 1.5+, model cannot exceed 1.0) | — |

**Reading (corrected — the first pass of this note had the sign backwards)**:
knee TIMING is essentially solved on the first attempt (better than the old
baseline even). `score_rpt_gap`'s convention is `model - real`, `+ve = model
retains more`. RPT4/RPT5 gaps are **+9.0pp/+12.2pp** (positive) — the model
RETAINS MORE capacity than real post-knee, i.e. it UNDER-degrades /
degrades too SLOWLY after the knee (not "too fast/too deep" as first
written). RPT2 (pre-knee, EFC=49.4) gap is -3.3pp — model degrades very
slightly faster than real pre-knee, a small effect. So the fix needed is:
INCREASE the post-knee degradation rate (isolation-redirect yield and/or
reaction-driven LAM), not decrease it. Likely cause: switching Si from
"ec reaction limited" to "solvent-diffusion limited" changes bulk SEI
current's (`a_j_sei`, the isolation-redirect's driving variable) post-knee
trajectory shape (continuously 1/L_sei-decaying vs ec-reaction-limited's
eta_SEI-sustained form) — the redirect yield calibrated for the old
mechanism's magnitude may now under-drive LAM_Si post-knee.

## New baseline's own benchmark (measured 2026-09-28, `new_baseline_benchmark_run`)

Ran `degradation_test_matrix/cell064_degradation_fit.py` with
`build_new_baseline_curves.py`'s exact env-var overrides. **This is the
real target**, not the old-baseline numbers above:

| Metric | New baseline (real target) | Old baseline (for reference only) | iter00 (this study) |
|---|---|---|---|
| Knee EFC | 116.1 (err +14.7) | ~108 (err +6.3) | **101.8 (err +0.4) — already better** |
| RPT-vs-RPT mean \|gap\| | **1.5pp** (RPT2 -0.9, RPT4 -2.5, RPT5 -1.2) | 1.9pp | 8.2pp — needs work |
| Expansion shape RMSE | 0.157 | ~0.10 | 0.175 — close, slightly worse |
| Voltage shape RMSE (mean) | 0.097 V | n/a | **0.079 V — already better** |
| k-gap mean \|gap\| | 0.291 | n/a | 0.306 — close, slightly worse |

Interesting: the new baseline's OWN knee (116.1) is actually further from
the real anchor (101) than the OLD baseline's (108) — the isolation
mechanism + radius/diffusivity retune improved other things (RPT gap 1.9→
1.5pp) at some cost to knee timing precision. study_si_deg's iter00 already
beats the new baseline on knee timing and voltage shape RMSE on the very
first attempt. The remaining real gap is the RPT-vs-RPT gap (8.2 vs 1.5pp)
and, to a lesser extent, expansion RMSE and k-gap (both close, slightly
worse). Given the corrected sign reading above, **closing the RPT gap means
increasing post-knee degradation rate** (model currently under-degrades
post-knee, e.g. +9.0pp/+12.2pp at RPT4/RPT5 relative to real, meaning model
retains too much capacity) — likely via `SI_REDIRECT_LAM_YIELD` and/or
`SI_BETA_LAM_SEI`, since bulk `a_j_sei` (the isolation-redirect's driving
variable) likely decays faster post-knee under "solvent-diffusion limited"
than it did under "ec reaction limited".

## Tuning campaign (2026-09-28, iter01-iter11)

Priority per instruction: close the RPT-vs-RPT gap first (8.2pp → target
1.5pp) by increasing post-knee degradation rate, while watching knee timing
(101.8, already better than new-baseline's 116.1) doesn't regress.

| iter | Change from iter00 | Knee EFC (err) | RPT gap mean\|gap\| | Expansion RMSE | k-gap mean | Voltage RMSE |
|---|---|---|---|---|---|---|
| 00 | baseline (SOLVENT_MULT=34, SI_CRACK_MULT=1.0, yield=0.02, exp=3) | 101.8 (+0.4) | 8.2pp | 0.175 | 0.306 | 0.079 V |
| 01 | yield 0.02→0.04 | 102.3 (+0.9) | 6.1pp | 0.186 | 0.305 | 0.073 V |
| 02 | yield→0.08 | 88.1 (-13.3) — knee regressed | 4.9pp | 0.218 | 0.304 | 0.071 V |
| 03 | yield=0.08 + exponent 3→6 (sharpens gate, recovers knee) | 101.7 (+0.3) | 4.3pp | 0.231 | 0.306 | 0.070 V |
| 04 | yield→0.12, exp=6 | 100.4 (-1.0) | 3.5pp | 0.254 | 0.306 | 0.073 V |
| 05 | yield=0.08, exp=6, F0 0.7→0.8 | 100.7 (-0.7) | 4.3pp | **0.254 (worse)** | 0.316 (worse) | 0.071 V |
| 06 | yield=0.08, exp=6, F0→0.6 | 102.6 (+1.2) | 4.3pp | **0.320 (worse)** | 0.293 (best k so far) | 0.070 V |
| 07 | yield→0.20, exp=6 | 98.9 (-2.5) | 2.9pp | 0.280 | 0.305 | 0.073 V |
| 08 | yield=0.20, exp 6→10 | **101.1 (-0.3)** | **2.6pp** | 0.281 | 0.306 | 0.074 V |
| 09 | yield→0.30, exp=10 | 99.3 (-2.1) — no further RPT-gap gain | 2.7pp | 0.307 | 0.306 | 0.074 V |
| 10 | yield=0.20, exp=10, SI_BETA_LAM_SEI 1e-7→5e-7 | 101.1 (-0.3) | 2.6pp — **no effect** | 0.281 | 0.306 | 0.074 V |
| 11 | yield=0.20, exp=10, SI_TAU_LAM_ISO 2e8→5e7 | 89.4 (-12.1) — knee regressed badly | 3.0pp — **no real gain** | 0.354 (worse) | 0.306 | 0.072 V |

**Findings from this sweep**:
- `SI_REDIRECT_LAM_YIELD` is the correct primary lever for the RPT gap —
  confirmed monotonic improvement 0.02→0.20 (8.2→2.6pp), but pushing yield
  alone regresses the knee; must be paired with a sharper
  `SI_LAM_ISO_EXPONENT` (3→6→10) to keep the gate's pre-knee leakage
  suppressed and hold the knee near 101.
- Diminishing/negative returns past yield=0.20: yield=0.30 (iter09) gave no
  further RPT-gap improvement and cost 2 EFC of knee accuracy. **iter08
  (yield=0.20, exponent=10, everything else at iter00's values) is the best
  point found on the RPT-gap/knee trade-off curve.**
- `SI_BETA_LAM_SEI` (reaction-driven LAM, a channel independent of the
  isolation gate) has **no measurable effect** even at 5x its default —
  confirms this session's earlier finding that bulk `a_j_sei` (what this
  term scales) is tiny under "solvent-diffusion limited" Si SEI, so this is
  not a useful lever here.
- `F0` moves the WRONG direction from the 45 degC work: raising it (0.7→0.8)
  worsens both expansion RMSE and k-gap here, because real `k` needs to rise
  ABOVE 1.0 post-knee (reaches ~1.5) — a higher `F0` shrinks the achievable
  range toward the model's 1.0 ceiling instead of widening it. Lowering it
  (0.7→0.6) helped k-gap marginally (0.293, close to the new-baseline's own
  0.291) but worsened expansion RMSE further (0.320). **F0=0.7 (unchanged
  from new baseline) is the best value found for expansion RMSE specifically**
  — the k-gap is dominated by the structural k≤1 ceiling regardless (see
  below), so it's not worth trading expansion RMSE for a small k-gap gain.
- `SI_TAU_LAM_ISO` reduction (making the gate saturate faster) was tried to
  address the persistent RPT5 (latest EFC) lag specifically — it badly hurt
  knee timing (89.4) for no real RPT-gap gain (3.0pp, worse than iter08's
  2.6pp). Not a useful lever at this setting; reverted.
- **Residual RPT5 lag**: even at the best config (iter08), RPT4 is matched
  almost exactly (-0.1pp) but RPT5 (the latest, EFC~197) still lags by
  +4.6pp (model retains too much capacity at the very end of the tested
  range). None of SI_BETA_LAM_SEI, SI_TAU_LAM_ISO, or further yield pushes
  closed this specific residual without costing other metrics. Not yet
  resolved — candidates not yet tried: NEG_POROSITY_FLOOR (deeper floor →
  more sustained post-knee decline), EXPONENT_MAX_SEI, or accepting this as
  the residual gap (2.6pp mean is still close to new-baseline's 1.5pp and
  far better than the un-tuned 8.2pp starting point).
- **Expansion RMSE has regressed from iter00 across every yield-increasing
  change** (0.175→0.28-ish) and has not been recovered by F0 tuning. This is
  the main OPEN issue: closing the RPT gap via more post-knee LAM_Si
  currently costs expansion-shape accuracy. Not yet resolved.

## Adopted config (2026-09-28) — iter08

**Baked into `cell064_degradation_fit_crack_dominated.py` as the new
defaults** (both changes documented inline at their definitions):
`SI_REDIRECT_LAM_YIELD`: 0.02 → **0.20**; `SI_LAM_ISO_EXPONENT`: 3.0 →
**10.0**. Everything else unchanged from iter00 (new-baseline values +
Si SEI = solvent-diffusion-limited, `SOLVENT_MULT=34`, `SI_CRACK_MULT=1.0`).
Running the script with no env var overrides now reproduces iter08.

### Final scores vs. the new-baseline benchmark

| Metric | New baseline (target) | iter08 (adopted) | Verdict |
|---|---|---|---|
| Knee EFC | 116.1 (err +14.7) | **101.1 (err -0.3)** | **Beats target** |
| RPT-vs-RPT mean \|gap\| | **1.5pp** | 2.6pp | Close, not quite matched |
| Expansion shape RMSE | **0.157** | 0.281 | Worse — open issue |
| Voltage shape RMSE (mean) | 0.097 V | **0.074 V** | **Beats target** |
| k-gap mean \|gap\| | 0.291 | 0.306 | Close, not quite matched |
| LAM/LLI visual | crude real LAM_Si estimate 58-80% (dma/FINDINGS_2026-09-20.md) | not re-checked visually this session — see caveat below | unverified |

**Honest assessment**: iter08 beats the new baseline on 2 of 5 quantitative
metrics (knee timing, voltage shape) and comes close on 2 others (RPT gap,
k-gap), but expansion shape RMSE is clearly worse (0.281 vs 0.157) and is
the main unresolved gap. This is NOT simply an untuned parameter — every
attempt to push the RPT-gap improvement further (higher
`SI_REDIRECT_LAM_YIELD`, sharper exponent) made expansion RMSE WORSE, and
`F0` tuning (the lever that fixed an analogous amplitude issue in the 45
degC work) moves the WRONG direction here (real k needs to exceed 1.0
post-knee, so raising F0 shrinks headroom rather than adding it). This
looks like a genuine structural tension in this specific config, not yet
resolved:

- **Overall, "as good as or better than the new baseline" is NOT yet fully
  achieved** — call it a partial success: the cracking-dominated mechanism
  gets knee timing and RPT voltage shape genuinely right (better than the
  new baseline even), and gets capacity-fade (RPT gap) and k-gap within
  reasonable range of the target, but has NOT matched the new baseline's
  expansion-shape fit.
- **Root-cause hypothesis for the expansion regression** (not verified):
  increasing `SI_REDIRECT_LAM_YIELD` increases post-knee LAM_Si (active
  material loss), which per `base_mechanics.py`'s `v_change = eps_s_eff *
  (t_change_now - t_change_init)` formula directly reduces the modeled
  particle-level expansion amplitude post-knee — on top of the pre-existing
  k≤1 structural ceiling (a real cell's k keeps rising past 1.0 post-knee;
  this pore-buffering formulation cannot represent that, a known,
  documented, pre-existing limitation per `CELL064_final_parameterisation.md`,
  not something introduced by this study or fixable by parameter tuning).
  So closing the RPT-gap by increasing post-knee LAM makes the ALREADY
  present k-ceiling mismatch worse in absolute normalised-RMSE terms, since
  it further suppresses the model's own post-knee expansion signal while
  the real signal keeps climbing.
- **`NEG_POROSITY_FLOOR` tried (iter12, 0.035→0.025, on top of iter08)**:
  not a clean win. RPT5's lag is essentially fixed (+4.6pp→+0.4pp), but
  RPT4 now overshoots the other way (-0.1pp→-4.4pp) — the floor change
  shifts WHERE the mismatch sits rather than removing it (mean gap
  unchanged at 2.6pp). It also made voltage-shape RMSE worse (0.074→0.090V)
  and expansion RMSE worse still (0.281→0.330), and cost 2 EFC of knee
  accuracy (101.1→99.0). Not adopted — reverted to iter08's 0.035.
- **Not yet tried, candidates for future work**: decoupling the LAM-rate increase from the expansion/thickness-change
  calculation somehow, or accepting the residual expansion-RMSE gap as a
  genuine structural cost of the cracking-dominated + isolation-redirect
  combination and prioritizing the (already good) SoH/knee/voltage-shape
  fit instead, depending on which metric the user weights most.
- **LAM/LLI visual comparison — DONE, and this is the real headline
  finding**: opened `cell064_overall_summary_result_iter08_yield0p20_exp10.png`.
  LLI matches real data closely (model ~40-47% vs real diamonds 42%/49% at
  EFC 171/197 — good). But **LAM_Si dramatically overshoots**: model climbs
  past 90% by EFC~230, well above the crude real ceiling of 58-80%
  (dma/FINDINGS_2026-09-20.md), and its sigmoidal rise appears to start
  earlier (~EFC 75-100) than the two real LAM_Si points (~0% at EFC~50, 58%
  at EFC~170) suggest it should. **This, not just the pre-existing k≤1
  ceiling, is almost certainly the real root cause of the expansion-RMSE
  regression** — pushing `SI_REDIRECT_LAM_YIELD` to 0.20 to close the
  aggregate SoH (RPT-gap) mismatch has overshot the LAM_Si split itself.
  The reversible-expansion plot confirms this directly: model peaks early
  (~EFC 75-100, real peaks ~EFC 100-125) and then crashes much harder than
  real (down to ~0.4 normalised by EFC~230 vs real staying ~0.6-1.2 over
  the same range) — consistent with too much active material being removed
  too fast, not just the k-ceiling alone.
  **Correction after checking iter00's and iter04's plots too**: the
  LAM_Si overshoot and the expansion peak-then-crash shape mismatch are
  BOTH already visible at iter00 (yield=0.02, before ANY of this session's
  tuning) and iter04 (yield=0.12) looks nearly identical to iter08
  (yield=0.20) in both panels. So this is NOT primarily something the
  yield/exponent tuning introduced or worsened much — it is a STRUCTURAL
  feature of the cracking-dominated solvent-diffusion-limited SEI +
  isolation-redirect mechanism itself, present from the very first working
  config. The yield tuning genuinely improved aggregate SoH (RPT-gap
  8.2→2.6pp) essentially "for free" without meaningfully worsening an
  already-present LAM_Si/expansion problem (expansion RMSE did drift
  0.175→0.281 across the sweep, but the qualitative shape mismatch was
  already there at 0.175, not created by the tuning).
  **Implication**: retuning yield down will NOT fix the LAM_Si/expansion
  mismatch — that requires addressing the underlying mechanism (why does
  cracking-dominated SEI growth push Si's own LAM curve up faster/higher
  and the expansion curve's peak earlier/sharper than the "ec reaction
  limited" new baseline does?). Two live hypotheses, neither tested this
  session: (1) the crack area's own roughness growth (SI_CRACK_MULT) is
  itself accelerating over life in a way that compounds with the isolation
  redirect — try REDUCING `SI_CRACK_MULT` below 1.0 and see if LAM_Si comes
  down without needing to touch yield; (2) `SOLVENT_MULT=34` may need
  re-lowering once the crack contribution's own scaling is better
  understood, since it was calibrated (this session, and in the earlier
  knee-delay work) purely against knee EFC, never against the LAM_Si
  trajectory shape. **This is the most important unresolved item for
  anyone continuing this study** — start by sweeping `SI_CRACK_MULT` down
  from 1.0 (e.g. 0.3, 0.5) at iter08's yield/exponent settings and check the
  LAM_Si plot directly, not just the RPT-gap number.

### iter13: tested the SI_CRACK_MULT hypothesis — surprising, important negative result

Ran `SI_CRACK_MULT` 1.0→0.3 (a 3.3x reduction) on top of iter08. Result:
**every metric is essentially unchanged** (RPT-gap 2.6pp, knee EFC=101.1,
expansion RMSE 0.284 vs 0.281 — all within noise) and the overall-summary
plot is visually indistinguishable from iter08's, including the LAM_Si
curve (still overshoots to >85%) and the expansion peak-then-crash shape.

**This is a bigger deal than a simple "lever doesn't work" result.** A
3.3x reduction in crack rate should, if cracking genuinely dominates
pre-knee porosity closure as the diagnostic found (65.6% at the 45 degC
recipe), have delayed the knee noticeably — it didn't move at all. Two
possible explanations, NEITHER confirmed:
1. Graphite's own SEI growth (`GR_K_SEI_MULT=6e-6`, "ec reaction limited",
   completely unaffected by `SI_CRACK_MULT`, which only touches Si) may be
   the actual dominant contributor to overall electrode porosity closure
   in THIS specific recipe, with Si's cracking-dominated mechanism playing
   a smaller role than assumed. If true, **this would undercut the study's
   own central hypothesis** — the knee timing here might be set mostly by
   an unremarkable, unchanged graphite mechanism, not by the deliberately-
   designed Si cracking mechanism at all.
2. Alternatively, the 65.6% split measured in the diagnostic (run against
   the 45 degC recipe, `SOLVENT_MULT=5`) may not transfer to this 25 degC
   recipe's `SOLVENT_MULT=34` — worth re-running `check_crack_vs_bulk_sei.py`-
   style diagnostics directly against THIS script's config before concluding
   graphite dominates.

**Whoever continues this study should resolve this before doing any more
Si-mechanism tuning** — if graphite's SEI is actually what's setting the
knee, the "cracking-dominated" story for the 25 degC fit needs rethinking
from the ground up, independent of the LAM_Si/expansion-RMSE issue above.

### Follow-up diagnostic (`check_graphite_vs_si_porosity.py`) — hypothesis #1 REFUTED

Built a direct diagnostic comparing graphite's own porosity-closing SEI
thickness against Si's combined (bulk+cracks) contribution, at iter08's
adopted config (crack_mult=1.0, the default). Result:
**graphite = 9.6%, Si (bulk+cracks) = 90.4%** of total porosity-closing SEI
thickness at end of run (EFC~229) — Si genuinely dominates, confirmed
directly, not assumed. **Hypothesis #1 (graphite secretly drives the knee)
is REFUTED.** Within Si's own total, the split matches the earlier 45 degC
diagnostic almost exactly (bulk 34.4%, cracks 65.6%, from the raw
per-EFC numbers: e.g. at EFC=228.7, L_si_bulk=849.4nm vs
L_si_cracks_contrib=(838.9nm raw)×(roughness-1=1.932)=1620.6nm, giving
849.4/(849.4+1620.6)=34.4% bulk).

**So the real open question is narrower than first thought**: Si (with
cracks genuinely dominant within it, 65.6%) sets ~90% of the porosity
closure — yet cutting `SI_CRACK_MULT` by 3.3x (iter13) still didn't move
the knee or any other metric. This diagnostic was run at the SCRIPT'S
DEFAULT `SI_CRACK_MULT=1.0`, not re-run at 0.3, so the exact mechanism
by which crack rate stops mattering isn't yet nailed down — plausible
explanation not yet checked: `SOLVENT_MULT=34` may be high enough that
BULK SEI alone (34% of the total, unaffected by crack rate) is already
sufficient to close porosity to the isolation-gate's trigger point within
the tested EFC range, making the ADDITIONAL cracks contribution "extra"
rather than rate-limiting for KNEE TIMING specifically (even though it's
the majority of total SEI mass) — i.e. the gate may saturate off of
whichever contribution reaches the threshold FIRST, and bulk alone might
already be enough, with cracks only affecting how far EXCESS porosity
closes beyond that point (relevant to post-knee behavior, not knee onset).
**Next step for whoever continues**: re-run this same diagnostic at
`C64_SI_CRACK_MULT=0.3` and `=3.0` to see how much the ABSOLUTE Si
contribution and its timing actually shift, and separately check
`Negative electrode porosity` directly (via `C64_DIAG_POROSITY=1`) against
`NEG_POROSITY_FLOOR` at each RPT to see when the gate's headroom actually
approaches zero, rather than reasoning about it indirectly through SEI
thickness alone.

## Tuning campaign 2 (2026-09-28, iter17–26) and adopted 25 °C baselines

Goal: keep iter16's knee/LAM_Si match, add graphite isolation LAM so that
Gr LAM climbs from the knee to the DMA values (~3.5% at EFC 171, ~9% at
EFC 197), and fix the post-knee shape. Scored on the same metrics as
before. "Onset" below = SoH-91% crossing (where the visible bend starts);
the printed "knee EFC" is the steepest-descent point, which sits 20–30 EFC
later.

### Two core PyBaMM bugs found and fixed along the way (uncommitted)

1. **Bulk OCP pushed onto the electrode mesh** (`base_ocp.py`,
   `_apply_ocp_aging_deformation`). The bulk LAM fraction divided the
   x-averaged eps_solid by the x-dependent `epsilon_s` FunctionParameter,
   promoting the bulk OCP to the "negative electrode" domain. This was
   harmless for Si, because its bulk OCP is diagnostic only. Once graphite
   (the primary phase) was in scope, it broke `ocv_bulk` / `"Local ECM
   resistance [Ohm]"` with a (5,1) vs (11,1) ShapeError. Fix: `x_average`
   the BOL reference.
2. **Aging-deformation options applied to every phase with porosity LAM**
   (`base_mechanics.py`, `base_ocp.py`, `base_battery_model.py`). Both
   `"open-circuit potential aging deformation"` and `"volume change aging
   deformation"` were global options, scoped only by the LAM option.
   Turning on `GR_REDIRECT_TO_LAM` therefore switched them on for graphite
   too. The volume-change one REPLACES the phase's own `t_change()` with Si's
   `1 + 3*sto^n`. That is never a no-op: for graphite it is ~30x its real
   volume change, which corrupted iter18's expansion/k (knee/SoH/LAM were
   unaffected, since `t_change` only feeds `Cell thickness change`). Fix:
   both options now accept per-phase tuples. The script passes
   `(("false", "true"), "false")`, i.e. Si only.

### Campaign summary

| Iter | Change | Onset | RPT gap | Exp. RMSE | Note |
|---|---|---|---|---|---|
| 18 | mult 35, Si yield 0.20, Gr LAM on (0.03) | ~150 | 6.1pp | (bugged) | knee late |
| 19a/b | mult 45 / 52 | 79 / 70 | 4.4 / 5.7 | 0.28 / 0.34 | knee early; Si LAM overshoots |
| 20 | mult 40, Si 0.15, Gr 0.015 | 89 | 2.3 | 0.19 | first good one |
| 21b | + Si tau 2e8→4e8 | 93 | 2.0 | 0.125 | tau spreads post-knee |
| 22a | mult 38, Si 0.17, Gr 0.02 | 96 | 1.8 | 0.117 | |
| 22c | tau 8e8 | 96 | 3.7 | 0.078 | best expansion; post-knee SoH too high |
| 23c | mult 36, Si 0.22, tau 8e8 | 104 | 3.4 | 0.074 | onset later, per user |
| 24c | Si 0.26, tau 6e8 | 101 | 2.2 | 0.097 | Si LAM 80% at EFC 197 |
| **25b** | **Si 0.32** | **99.6** | **1.7** | **0.115** | **→ 25degC_baseline1** |
| 26 | tau 5e8 | 98.4 | 1.8 | 0.131 | **→ 25degC_baseline2** (Si LAM ~82.5%, worse V fit) |

Findings from the campaign:
- **Knee timing vs SOLVENT_MULT is steep and non-linear:** mult 35 → 45
  moves onset ~150 → ~79. Raising mult also raises the SEI current at the
  moment the gate trips, so the post-knee Si LAM rise gets steeper too.
  Timing and post-knee slope therefore can't be tuned independently with
  mult alone.
- **The Si gate time constant (`SI_TAU_LAM_ISO`) is the post-knee width
  lever:** it barely moves the onset.
- **Si redirect yield has diminishing returns above ~0.25**, because the
  redirect term self-limits (`remaining_frac` and area both shrink).
- **Post-knee LLI stays ~4pp short of the real 49% at EFC 197 whatever the
  yield.** A direct particle-inventory check (`Total lithium lost from
  particles`) confirmed there is no hidden lithium loss behind the plotted
  LLI: the true inventory loss is ~2–3pp *lower* than side+LAM-trapped.
- **BoL porosity** (0.17→0.14 on iter17's config, "iter17b") moves the knee
  a lot (onset ~150 → ~87) but also sharpens the collapse.

### Adopted defaults: 25degC_baseline1 (= iter25b)

Now the script's defaults (a run with no `C64_*` env vars reproduces it):
`SOLVENT_MULT=36`, `SI_REDIRECT_LAM_YIELD=0.32`, `SI_TAU_LAM_ISO=6e8`,
`SI_LAM_ISO_EXPONENT=10`, `F0=0.7`, `WIDTH=3e-4`, `SI_CRIT_STRESS=3e7`,
`GR_REDIRECT_TO_LAM=1`, `GR_REDIRECT_LAM_YIELD=0.02`,
`GR_LAM_ISO_EXPONENT=10`, and both Si deformations on (Si only).

| Metric | baseline1 | baseline2 (tau 5e8) |
|---|---|---|
| Onset (91% crossing) / steepest knee | 99.6 / 128 | 98.4 / 119 |
| RPT gap at EFC 49 / 171 / 197 | −2.8 / 0.0 / +2.3pp (mean 1.7) | −2.8 / −0.6 / +1.9 (1.8) |
| Expansion RMSE / k gap | 0.115 / 0.305 | 0.131 / 0.305 |
| Voltage-shape RMSE | 0.075 V | 0.086 V |
| Si LAM at EFC 171 / 197 (real 58 / 80%) | ~72 / ~81 | ~75 / ~82.5 |
| Gr LAM at EFC 197 (real 9%) | ~9.5 | ~9 |
| LLI at EFC 197 (real 49%) | 45.2% | 45.6% |

Plots and logs: `*_25degC_baseline1.*`, `*_25degC_baseline2.*`.

Known gaps:
- Si LAM is ~14pp high at EFC 171. The model's Si LAM decelerates after the
  knee, whereas the real data accelerates from 58% to 80% between RPT4 and
  RPT5.
- RPT5 is +2.3pp.
- k is capped at 1.0 by the pore-buffering formulation (real ~1.5).

## The iter13 puzzle resolved: Si cracks are static in this recipe (2026-09-29)

A 20-cycle simulation that tracked the crack variables directly gave:

| | 25 °C | 45 °C |
|---|---|---|
| Si crack length | 2.000e-8 → 2.001e-8 m | no change |
| Roughness | ~2.908 at both temperatures | ~2.908 |
| Peak Si tangential stress | 4.5e7 Pa | 2.0e7 Pa |

- Roughness is set entirely by the INITIAL crack length (2e-8 m) × crack
  density (3.18e15 m⁻²).
- Crack growth (`k_cr = 3.9e-20`, Paris `m = 2.2`) is negligible at either
  stress level.

That is why iter13's `SI_CRACK_MULT` 1.0 → 0.3 changed nothing. It is also
why `SI_CRACK_EAC`, including the high_temp_45C fork's own
`SI_CRACK_EAC=38000`, is inert. "SEI on cracks" still supplies ~65% of Si
pore closure, but only on that fixed initial crack area.

`SI_LAM_EAC` is also inert here, because Si stress stays below the 3e7 Pa
critical stress.

**45 °C from baseline1, unchanged** (`C64_T_AMBIENT_K=318.15`, new switch):
- Onset ~48, i.e. EARLIER than 25 °C, not the real ~300.
- Reached the 45% floor at cycle 539 (EFC ~293).
- SEI Arrhenius (E_sei = 38 kJ/mol, ×2.6) wins, and nothing crack-related
  can respond.
- Adding `SI_CRACK_EAC` = 38k or 80k plus `SI_LAM_EAC` = 40k gave SoH
  identical to 0.01pp.

## Crack-growth sweep: making the reduced-cracking mechanism real (2026-09-29)

Hypothesis (user): higher T → faster diffusion and lower overpotential →
lower stress → less crack growth → less crack SEI → slower pore closure →
later knee. For that to act, cracks must actually GROW at 25 °C.

**Coupling:** cracking doesn't change the SEI law per unit area. It adds
fresh surface: "SEI on cracks" grows on `a·(roughness − 1)` with its own
thickness. New crack area thins that film, and diffusion-limited j ~ 1/L,
so crack growth raises the SEI growth rate per electrode volume. That feeds
both pore closure (`L_sei_cr·(roughness − 1)`) and LLI.

**Setup:** `crack_growth_sweep.py` runs the high_temp_45C fork unmodified
apart from env vars. It keeps the fork's own `SOLVENT_MULT=5` and prints a
per-RPT crack table: `l_cr/l0`, roughness, peak stress, crack share of pore
closure, and bulk vs crack LLI. It also saves `crack_diagnostics_<tag>.png`.
Two opt-in switches were added to the fork: `C64_T_AMBIENT_K` (default
318.15) and `C64_DIAG_EXTRA_VARS`. At 298.15 K the fork now overlays and
scores against the real 25 °C data.

**25 °C, `SOLVENT_MULT=5`:**

| `SI_CRACK_MULT` | l/l0 at EFC ~300 | Crack share of pore closure | Crack LLI | Onset |
|---|---|---|---|---|
| 1 | 1.005 | 65% | 3.3% | none (91% SoH at EFC 309) |
| 1e2 | 1.6 | 72% | 4.3% | none (~91% at 309) |
| 1e3 | 14 | 94% | 24% | 130 |
| **1.5e3** | 14.6 | 94% | 27% | **106** |
| 2e3 | 14.8 | 95% | 29% | 92 |
| 3e3 | 14.9 | 95% | 31% | 76 |
| 1e4 | 15.1 (saturated by EFC ~50) | 96% | 34% | 48 |

Cracks saturate at ~300 nm (~15× l0), probably a geometric cap tied to the
retuned small Si radius. Shapes: a genuine knee and a dip-then-hump
expansion, the shape the fork targeted at 45 °C. Stress feedback: after the
knee, the remaining active Si carries more current and stress rises to
~110 MPa.

**45 °C, same parameters (`SI_CRACK_MULT=1.5e3`):**

| `SI_CRACK_EAC` | Pre-knee peak stress | l/l0 at EFC ~270 | Onset | Steepest knee |
|---|---|---|---|---|
| 25 °C reference | ~43 MPa | 14 | 106 | 93 |
| 0 (stress only) | ~16–23 MPa | 2.6 | 213 | 263 |
| 38000 | ~16–20 MPa | 1.4 | 263 | **310** |

**Result:**
- One parameter set gives a 25 °C onset ~106 and a 45 °C knee ~310. The
  delay comes from genuinely less cracking at high T.
- Lower stress alone doubles the onset EFC, with no activation energy on
  cracking. A physically plausible `SI_CRACK_EAC` adds the rest.

**Caveats:**
- This is the fork's recipe at 25 °C, not a fitted 25 °C baseline. Its
  post-knee decline is too gentle (~61% SoH at EFC 320 against the real
  ~50% at EFC 197).
- Its LAM/LLI split has not been checked against DMA. That run
  (`fork25C_crackx1.5e3_vsData`) was in progress at the time of writing.
- Deferred (user, 2026-09-29): tune 45 °C against an onset-based target (~300
  onset would need `SI_CRACK_EAC` ~50000).

## Fitting the cracking-driven config to the real 25 °C data (2026-09-29)

Same base as the crack-growth sweep: the high_temp_45C fork run through
`crack_growth_sweep.py` with `C64_T_AMBIENT_K=298.15`, `SOLVENT_MULT=5`
(fork default) and `SI_CRACK_MULT=1.5e3`. At 298.15 K the fork now overlays
and scores the real data. Per-phase, Si-only deformation scoping was ported
into the fork so graphite isolation LAM can be switched on safely.

| Si yield | Gr isolation (yield, exp 10) | F0 | Onset | RPT gap (49/171/197) | Mean | Exp. RMSE | V RMSE |
|---|---|---|---|---|---|---|---|
| 0.02 (fork) | off | 0.85 | 105.7 | −1.0/+19.3/+22.2 | 14.1 | 0.085 | 0.078 |
| 0.08 | 0.02 | 0.85 | 99.7 | −1.0/+3.3/+5.1 | 3.1 | 0.159 | 0.073 |
| 0.15 | 0.02 | 0.85 | 98.3 | −1.0/+0.2/+2.7 | 1.3 | 0.223 | 0.069 |
| 0.25 | 0.02 | 0.85 | 97.1 | −1.0/−1.8/+1.2 | 1.4 | 0.271 | 0.079 |
| 0.10 | off | 0.85 | 100.8 | −1.0/+9.2/+12.4 | 7.5 | 0.163 | 0.075 |
| 0.15 | 0.005 | 0.85 | 99.2 | −1.0/+4.7/+7.9 | 4.5 | 0.210 | 0.071 |
| 0.15 | 0.01 | 0.85 | 98.8 | −1.0/+2.9/+5.8 | 3.2 | 0.215 | 0.070 |
| 0.15 | 0.01 | **0.7** | 97.9 | −1.0/+2.8/+5.8 | 3.2 | **0.118** | 0.069 |
| 0.15 | 0.013 | 0.77 | 98.3 | −1.0/+1.9/+4.8 | 2.6 | 0.160 | 0.069 |
| 0.18 | 0.01 | 0.75 | 97.8 | −1.0/+2.0/+5.3 | 2.8 | 0.163 | 0.069 |
| **0.15** | **0.013** | **0.7** | **97.8** | **−1.0/+1.9/+4.8** | **2.5** | **0.121** | **0.069** |

Findings:
- **Crack SEI supplies the late LLI** that `25degC_baseline1` could never
  reach. With low Si yield, LLI matched both DMA points exactly.
- **Graphite yield is a strong SoH lever** once Si is ~80% lost, because
  graphite capacity becomes limiting.
- **Graphite yield needs to be ~2× lower here than in baseline1.** The
  fork's graphite SEI current is larger, so at 0.02 Gr LAM overshoots to
  ~17% and LAM-trapped Li lifts LLI ~5pp.
- **`F0=0.7` fixes the post-knee expansion tail and k at BoL (0.70).** The
  cost is a peak overshoot (~1.47 against ~1.23). Higher `F0` lowers the
  peak but worsens the RMSE.

**Candidate "crack baseline"** (last row): fork + `SI_CRACK_MULT=1.5e3`,
`SI_REDIRECT_LAM_YIELD=0.15`, `GR_REDIRECT_TO_LAM=1`,
`GR_REDIRECT_LAM_YIELD=0.013`, `GR_LAM_ISO_EXPONENT=10`, `F0=0.7`.
- **Remaining gaps at 25 °C:** RPT5 +4.8pp, Si LAM ~14pp high at EFC 171,
  Gr LAM ~10.5/12.5% against 3.5/9%, and the expansion peak overshoot.
- **At 45 °C** the Gr-0.01 version with `SI_CRACK_EAC=38000` gives onset
  248 and steepest-descent knee 284, with a sharp visible knee at ~265 and
  a dip-then-hump expansion (peak ~1.27). It reaches the 45% floor at cycle
  882. Moving the knee to ~300 is deferred 45 °C tuning (`SI_CRACK_EAC`
  ~45000).

## crack_baseline_v1 — one parameter set for 25 °C and 45 °C (2026-09-29)

Saved as `crack_baseline_v1.env` (run it the same way as v2, below). It uses the high_temp_45C fork, via
`crack_growth_sweep.py`, with the fork's defaults for everything not listed,
including `SOLVENT_MULT=5` and `SI_ESEI=GR_ESEI=38000`. **The only difference
between the two conditions is `C64_T_AMBIENT_K`.** `SI_CRACK_EAC` and
`SI_LAM_EAC` are exactly inert at 298.15 K.

| Parameter | Value | Role |
|---|---|---|
| `SI_CRACK_MULT` | 1.3e3 | makes Si cracks actually grow (static at 1), setting knee timing at both T |
| `SI_TAU_LAM_ISO` | 4e8 | knee width |
| `SI_REDIRECT_LAM_YIELD` | 0.30 | post-knee Si LAM |
| `GR_REDIRECT_TO_LAM` / yield / exponent | 1 / 0.016 / 10 | graphite isolation LAM |
| `F0` | 0.7 | expansion tail + k at BoL |
| `SI_CRACK_EAC` | 70000 | Si cracking ÷5.9 at 45 °C (on top of ~2× lower stress); positions the 45 °C knee |

The 45 °C real data (CELL017, 15 psi) is in
`degradation_test_matrix/experimental_data/CELL017_45C_capacity_fade.csv`.
It uses capacity / 2.5 Ah, EFC since its first RPT, and is overlaid and
scored by the fork at 318.15 K.

| | 25 °C (vs CELL064 DMA) | 45 °C (vs CELL017 SoH) |
|---|---|---|
| Onset (91% crossing) / steepest knee | 105.5 / 117.5 (real ~101) | 272.6 / 310.9 (real ~256) |
| RPT gaps | −0.9 / +2.8 / +5.0pp at EFC 49 / 171 / 197 (mean 2.9) | −0.5 / +0.3 / +0.4 / +0.2 / **−5.3** / **+4.7**pp at EFC 51 / 148 / 197 / 245 / 334 / 400 (mean 1.9) |
| Other | expansion RMSE 0.105, V RMSE 0.070 V | model-only otherwise |

The 25 °C result is `opt25C_siy0.30_gry0.016` (Ea is inert there), copied as
`*_crack_baseline_v1_25degC*`. The 45 °C row is from
`crack_baseline_v1_45degC`, run with exactly this set. The nearest earlier
run, `opt45C_siy0.22_gry0.016_Ea70000` (Si yield 0.22), gave onset 278, gaps
−0.5/+0.3/+0.4/+0.4/**−3.5/+6.1** and a mean of 1.9pp.

Known gaps, the same at both T:
- Post-knee decline is too gradual.
- At 25 °C the true particle-inventory LLI is 34.6/38.5% against the DMA's
  41.7/49.1% at EFC 171/197.
- At 45 °C the real 334→400 collapse (−26pp) is steeper than the model's
  (~−16pp).

Two leads that were ruled out:
- A smaller initial crack length (`C64_SI_CRACK_L0`) is only a timing knob.
  At a matched onset it gives identical LLI, because cracks hit the
  particle-radius cap either way.
- Si yield barely matters at 45 °C.

Next under test: the opt-in `"stress-driven LAM damping": "power"`
(`C64_SI_STRESS_LAM_DAMP_EXP`), aimed at a late, accelerating Si LAM.

## crack_baseline_v2 — refit on the fixed core (2026-09-30, retune R1–R3)

**Why a refit.** `crack_baseline_v1` and everything before it ran with the
composite SEI-on-cracks bug (CHANGES.md item 14): the Si crack-SEI current
entered the charge balance with area $a$ instead of $a(\rho-1)$. Film growth
and the LLI counter used the crack area, so they did not match the lithium
that actually left the particles. That caused the graphite over-lithiation at
top of charge and the low late-life 25 °C C/20 curves (the "overpotential").
After the fix, `lithium_budget.py` closes to ≈0.

Saved as `crack_baseline_v2.env`. To run it (the shell drivers were removed
in the 2026-09-30 tidy-up):

```bash
set -a; source crack_baseline_v2.env; set +a
C64_T_AMBIENT_K=298.15 C64_F0=0.7 C64_OUT_TAG=v2_25degC python crack_growth_sweep.py   # 45 degC: 318.15 / 0.8
```
It is **one parameter set**; only $T$ and $F_0$ (0.70 / 0.80, the first-RPT
k) differ between the two cells.

| Parameter | v1 | **v2 (R3c)** |
|---|---|---|
| `SI_CRACK_MULT` | 1.3e3 | **1.0e3** |
| `SI_CRACK_EAC` | 70000 | **90000** |
| `SI_ESEI` | 38000 (fork) | **35000** |
| `SI_TAU_LAM_ISO` / `SI_LAM_ISO_EXPONENT` | 4e8 / 15 | **3.5e8 / 20** |
| `SI_REDIRECT_LAM_YIELD` / `SI_YIELD_EAC` | 0.30 / – | **0.25 / 35000** (0.61 at 45 °C) |
| `GR_REDIRECT_LAM_YIELD` / `GR_YIELD_EAC` / exp | 0.016 / – / 10 | **0.018 / 43000 / 10** (0.054 at 45 °C) |
| `SI_STRESS_LAM_DAMP_EXP` | (linear) | **0.25** (`"power"`) |

**Retune iterations.** Every run uses the fixed core. Gaps are SoH model − real,
in pp.

| Run | Change vs. previous best | 25 °C mean gap | 45 °C mean gap |
|---|---|---|---|
| fixed_* | v1 + fix, no retune | 2.2 (onset 84, too early) | 3.0 |
| R1a–d | crack mult 0.9–1.0e3, yields | 3.1–4.7 | 2.4–3.3 |
| R2a | crack 1.0e3, Si 0.28, Gr 0.018, τ 5e8, Si-yield Ea 35k | 2.5 | 1.6 |
| R2b | τ 4e8, Si 0.25 | 2.3 | 1.4 |
| R2c (45 only) | Si-yield Ea 50k, τ 5e8 | – | 1.2 |
| R3a | R2b + isolation exponent 25 | 2.3 | 1.7 |
| R3b | R2b + Si 0.28, Si-yield Ea 30k | 2.1 | 1.4 |
| **R3c** | **R2b + τ 3.5e8** | **2.0** | **1.2** |

**R3c in detail.**

| | 25 °C (CELL064) | 45 °C (CELL017) |
|---|---|---|
| SoH gaps | −2.3 / +1.1 / +2.7 at EFC 49 / 171 / 197 | −1.0 / −0.8 / −0.9 / −0.8 / +0.1 / +3.8 at EFC 51 / 148 / 197 / 245 / 334 / 400 |
| Onset (91%) / steepest knee | 96.6 / 129.6 (real ~101) | 257.8 / 346.6 (real knee ~310) |
| Gr LAM, last RPT | 8.3% (real 9.1%) | 10.2% (real 10.5%) |
| Si LAM | 64.6 / 80.4% (real 58.5 / 79.6) at EFC 171 / 197 | 22.0 / 67.5% (real 8.0 / 58.0) at EFC 334 / 400 |
| Particle LLI | 38.1 / 42.8% (real 41.7 / 49.1) | 15.8 / 37.5% (real 15.8 / 41.6) |
| Expansion RMSE / V RMSE | 0.113 / 0.071 V | 0.116 / 0.054 V |

What changed and what's left:
- **Fixed:** the late 25 °C RPT curves now overlay the real RPT4/5; the
  offset was the bug, not resistance or plating.
- **Knee width:** a shorter τ (3.5e8) sharpens the post-knee drop at both
  temperatures without moving the onset. The isolation exponent (R3a) and
  the Si-yield Ea (R3b) did less.
- **Remaining 45 °C gap:** Si LAM at EFC 334 is too high (22% vs 8%) while
  SoH at 400 is still +3.8pp. The model's Si LAM starts too early and ends
  too shallow; the real Si collapse is later and steeper.
- **Remaining 25 °C gap:** RPT5 is +2.7pp, and particle LLI at EFC 197 is 6pp
  short of the DMA.
- **Pre-knee:** the model sits 1–2pp below the data (−2.3pp at 25 °C
  EFC 49). That's BoL capacity or early SEI, not the knee.
- **k:** still monotonic, so the dip is not reproduced (deferred, part 3).

**Width follow-up (R4, 2026-09-30).** R3c with only `C64_WIDTH`, the pore
buffering transition width $w$, changed:

| Run | $w$ | SoH mean gap | Expansion RMSE | k at EFC 98 (real 0.949) |
|---|---|---|---|---|
| R3c 25 °C | 3e-4 | 2.0 | 0.113 | 0.757 |
| R4c 25 °C | 5e-3 | 1.5 | 0.110 | 0.871 |
| R4d 25 °C | 5e-2 | 1.4 | 0.143 (k rises too early, 0.81 at EFC 49) | 0.974 |
| R3c 45 °C | 3e-4 | 1.2 | 0.116 | – |
| R4 45 °C | 2e-4 / 1.5e-4 | 1.7 / 2.0 | 0.104 / 0.099 | – |

- **Width is not expansion-only.** Buffered swelling occupies pore volume,
  $\varepsilon = \varepsilon_{\mathrm{struct}} - \Delta v_{\mathrm{buf}}$
  (`reaction_driven_porosity.py`). That feeds the isolation gate and
  transport, so $w$ also shifts SoH and LAM.
- **At 45 °C,** a narrower $w$ improves the expansion shape but makes the
  post-knee SoH worse and moves the knee later.
- **Kept R3c ($w$ = 3e-4 at both T) for now.** R4c is a promising 25 °C
  expansion option, at the cost of more Si LAM at EFC 171 (68% vs 58.5%).
- **If adopted later:** a T-dependent $w$ (5e-3 at 25 °C, ≤ 3e-4 at 45 °C)
  implies $E_w \approx -110$ kJ/mol.
- **Limit of any $w$:** k is capped at 1, but the real 25 °C k reaches ~1.5
  after the knee.

**Measured-k test (option 2, 2026-09-30, `k_compliance_test.py`).** R3c at
45 °C with degradation unchanged. The cell thickness was rebuilt afterwards
with the measured CELL017 k, linearly interpolated between RPTs, in place of
the model's compliance k.

| Variant | Expansion RMSE |
|---|---|
| Model's own compliance k | 0.116 |
| A: measured k on negative-electrode swelling only | 0.128 |
| B: measured k on the whole electrode stack (neg + pos) | **0.094** |

- **B reproduces the shape.** It gets the pre-knee dip and the hump height.
  Its tail at EFC 400 is too high: 0.745 against a real 0.644.
- **A overdoes both.** The dip is too deep and the post-knee overshoot too
  large. The positive electrode is 39% of the model's breathing, so which
  basis is right depends on how rev_particle was defined in the data.
- **Particle-level swelling** (model vs real rev_cell/k_exp) has RMSE
  0.095:
  - Pre-knee, the model's swelling falls ~5% while the real one stays flat.
    That matches the model's early Si LAM of 2–4% against a real ~0–1.5%.
  - At EFC 400 the model's swelling is too high (0.44 vs 0.37). That
    matches too little late Si LAM, and SoH +3.8pp.
  - Between EFC 250 and 305 the real curve drops earlier than the model's.
    This is at least partly the linear k interpolant between RPT5 and RPT6.
- **Conclusion:** most of the 45 °C expansion gap is k, both the pre-knee
  dip and the rise above 1. The rest is the Si LAM timing: too early before
  the knee, too little at EFC 400.

**25 °C late-RPT voltage offset: polarisation vs OCV (2026-09-30).** The
per-temperature configs (25 °C = R4c, 45 °C = R3c) are plotted with
`rpt_soc_plots.py`: RPT voltage and dV/d(SOC) vs depth of discharge. At
45 °C, RPT1–6 overlay the data. At 25 °C, RPT4/5 sit about 0.1 V low
across the whole discharge.

- **Overpotential breakdown** (`overpotential_breakdown.py`, R4c). The
  late-life growth is almost all electrolyte transport:
  - Electrolyte ohmic + concentration losses rise from ~6 to ~45 mV once
    the negative porosity reaches its 0.035 floor at the knee (EFC 129).
  - SEI film adds ~6 mV.
  - The isolation gate reads the same porosity and floor as transport
    (`loss_active_material.py`).
- **Parameter-only levers tried at 25 °C** (the 45 °C versions were stopped):
  - SEI resistivity ×0.1: no voltage change.
  - Electrolyte conductivity and diffusivity ×2: recovers ~20 mV, but SoH
    after the knee rises to +9/+11pp and the onset moves 7 EFC earlier.
- **C/50 test** (`c50_rpt_test.py`). R4c is unchanged; a C/50 discharge
  branches off the same charged state at every RPT. Mean OCV − V at C/20
  vs C/50 is 20 → 8 mV at BoL and 70 → 37 mV late. At RPT4/5 the C/50
  curve recovers about half of the offset.
  - **The rest is OCV.** Even the model OCV sits 20–40 mV below the real
    C/20 curve mid-discharge, and its end-of-discharge drop comes earlier.
    The dV/d(SOC) shape is the same at C/20 and C/50, so the shape mismatch
    is thermodynamic: electrode balance or LLI.
  - **Conclusion:** consistent with the lumped model over-predicting pore
    polarisation, which is about half of the offset. The other half is an
    OCV / balance deficit.
  - **Caveat:** polarisation only halves for a 2.5× lower current, so part
    of it doesn't scale with rate.

**θ isolation-lithium-trapping re-test and aligned RPTs (2026-09-30).** θ
(CHANGES.md item 13) was re-added and re-tested on the fixed core. The runs
are R6a (θ = 0.5) and R6b (θ = 1.0), each on top of the per-temperature
final config. The knee and Si LAM are unchanged at both temperatures,
because θ adds no SEI film and fills no pores.

| | 25 °C θ = 0 / 0.5 / 1.0 | 45 °C θ = 0 / 0.5 / 1.0 | Real (25 / 45 °C) |
|---|---|---|---|
| LLI at last RPT | ~43 / 43.7 / 47.4% | 37.5 / 38.9 / 40.4% | 49.1 / 41.6% |
| Gr LAM at last RPT | 8.8 / 7.6 / 6.6% | 10.2 / 7.5 / 5.4% | 9.1 / 10.5% |
| Mean SoH gap | 1.5 / 1.5 / 2.7pp | 1.2 / 1.2 / 1.1pp | – |

- **At 45 °C,** θ = 1.0 brings the EFC 400 gap from +3.8 to +1.2pp.
- **Graphite LAM falls roughly in proportion to θ** at both temperatures,
  because graphite LAM comes from the graphite SEI redirect. Adopting θ would
  need the graphite redirect yield re-raised.

**Aligned RPTs** (`aligned_rpt_cycling.py`, `rpt_soc_plots.py --aligned`).
The model RPTs are placed within 0.3 EFC of the real ones. The 50-cycle
RPT schedule had put them 10–25 EFC off after the knee, and at 25 °C the
voltage gap changes ~50 mV per model RPT. Mean V gap (model − real, 5–95%
DoD) at the aligned RPT4 / RPT5:
- R4c: −84 / −94 mV.
- R6b θ = 1.0: −68 / −75 mV.

So θ = 1.0 recovers ~16–19 mV of the ~90 mV gap. With the RPTs aligned,
R4c's SoH and LAM are essentially unchanged: mean gap 1.5pp, steepest knee
126.0.

## crack_baseline_v3 — adopted baseline (2026-09-30, R7 / R7a)

Saved as `crack_baseline_v3.env`. It is crack_baseline_v2 plus:
- **θ = 1.0** (`C64_SI_ISO_LI_TRAP`): isolation lithium trapping, CHANGES
  item 13.
- **Graphite redirect yield 0.018 → 0.025,** restoring the Gr LAM that θ
  had lowered.
- **`GR_YIELD_EAC` 43 → 56.5 kJ/mol,** because θ lowered Gr LAM more at
  45 °C.

Per temperature, only T, F0 and WIDTH differ: 25 °C 0.70 / 5e-3; 45 °C
0.80 / 3e-4. To run it, with RPTs aligned to the real EFCs:

```bash
python rpt_soc_plots.py 25 --aligned     # and 45
```

| Aligned RPTs | 25 °C (R7) | 45 °C (R7a) | Real (25 / 45 °C) |
|---|---|---|---|
| SoH gaps | −2.2 / −5.0 / −2.6pp at EFC 49 / 171 / 197 | −0.9 / −0.8 / −0.9 / −0.7 / +0.8 / +0.5pp at EFC 51 / 148 / 197 / 245 / 334 / 400 | – |
| Mean SoH gap | 3.3pp | **0.8pp** | – |
| LLI at last RPT | 48.7% | 41.3% | 49.1 / 41.6% |
| Gr LAM at last RPT | 8.4% | 8.8% | 9.1 / 10.5% |
| Si LAM at last RPT | 83.6% | 62.3% | 79.6 / 58.0% |
| Onset / steepest knee | 98.7 / 126.6 | 257.6 / 361.1 | ~101 / ~310 |
| Expansion RMSE | 0.127 | 0.111 | – |
| Voltage, SOC-normalised | mean gap −66 / −72 mV at RPT4 / RPT5 | RMSE 8–16 mV (RPT1–5), 30 / 36 mV (RPT6 / 7) | – |

Chosen over v2 (no θ) for its more physical LLI at both temperatures, the
better 25 °C voltage fit, and the best 45 °C fit. On SoH alone, v2 is more
balanced (25 °C 1.5pp, 45 °C 1.2pp, with RPTs not aligned at 45 °C).

Known gaps:
- **Si LAM runs ahead of the DMA** at 25 °C EFC 171–197 (+4–9pp) and at
  45 °C EFC 334 (18.6% vs 8%). Gr LAM at 45 °C EFC 334 is too low (0.9% vs
  6.8%).
  - This compresses the Si-dominated end of discharge, so the late-RPT
    dV/d(SOC) drops early at both temperatures.
  - Candidate fix, untested: lower `SI_LAM_EAC` (acts at 45 °C only) plus a
    slightly lower Si redirect yield.
- **25 °C SoH undershoots after the knee:** −5.0pp at EFC 171.
- **About half of the remaining 25 °C voltage offset is lumped-porosity
  polarisation** (C/50 test).
- **The k dip and k > 1 are not modelled.** The measured-k post-hoc test
  shows that is most of the 45 °C expansion gap.

## Stage 3: partial-SoC-window cycling at 25 °C (2026-09-30, crack_baseline_v3)

**Cells.** CELL009 (15–95% SoC) and CELL026 (20–80% SoC), both 25 °C and
103 kPa. The data is built by
`experimental_data/build_partial_soc_data.py` into `CELL009_25C_15-95/` and
`CELL026_25C_20-80/`. CELL025, which I had assumed was the second
partial-window cell, is actually a 45 °C twin of CELL017. k at RPT1 is 0.70
for both, so F0 = 0.70.

**Real ageing control** (detected from the lifetime data, constant over life):
- **CELL009:** voltage-limited at both ends. C/3 CC-CV to 4.150 V (CV to
  14 mA), then C/3 CC discharge to 3.140 V.
- **CELL026:** C/3 CC-CV to 3.969 V (CV to 14 mA), then a coulomb-counted
  1.500 Ah C/3 discharge with a 2.5 V safety floor.

The model reproduces both exactly (`rpt_soc_plots.py 25p1595|25p2080 --aligned`,
`PARTIAL_AGEING`). RPTs are aligned to every real RPT, including the
capacity-only ones. Each RPT is preceded by a 10 min rest and a charge to
4.195 V: CC-CV for CELL026, CV only for CELL009.

| | CELL009 15–95% | CELL026 20–80% |
|---|---|---|
| SoH gap, pre-knee | −2.8 to −4.7pp (EFC 79–260) | −3.1 to −5.7pp (EFC 62–231) |
| SoH gap, post-knee | +12.6 / +34.8 / +41.3pp at EFC 329 / 382 / 402 | +11.8 / +40.4pp at EFC 324 / 383 |
| LLI at last RPT | 12.8% (real 53.8%) | 12.0% (real 51.0%) |
| Si LAM at last RPT | 0.4% (real 73.4%) | 0.3% (real 79.3%) |
| Gr LAM at last RPT | 2.4% (real 11.3%) | 2.2% (real 4.6%) |
| C/20 voltage-shape RMSE | 0.093 V | 0.082 V |

**Result: v3 predicts no knee for either window.** Both real cells knee
at about EFC 290–300, only ~2× CELL064's (full window). The model fades
slightly too fast before the knee (LLI ~6% vs ~4%) and then barely degrades.

**Working hypothesis.** The 25 °C knee in v3 runs through Si crack growth.
Paris law with m = 2.2 on the stress swing leads to crack-SEI, then pore
closure, then the isolation gate. A narrower window cuts the Si stress swing,
so the cracks never grow enough to close the pores. The real knee barely
depends on window size (CELL009 293 vs CELL026 298 EFC), so the real trigger
is much less stress-swing-sensitive than the model's.

**Next.** Confirm with crack length, roughness and porosity vs EFC for full
vs partial windows. Then look for a trigger that scales with throughput or
time rather than with stress swing.

**Workflow notes.**
- **SoH normalisation.** For partial windows the fork now normalises the
  model RPT SoH by the model's own RPT1 (`_rpt_soh_ref`). The first C/3
  ageing discharge is only the window's capacity.
- **Solver bug.** A step that ends instantly makes pybamm's
  `Solution.__add__` leave the stitched solution's t one sample shorter than
  its variables, when output variables are restricted. The RPT pre-charge
  avoids it; the core fix was spawned as a separate task.

### CELL009 (15–95%) partial-window tuning: P1 sweep (2026-10-01)

These are **CELL009-specific deviations from the shared 25/45 °C set
(crack_baseline_v3).** Everything else is v3 at 25 °C (F0 0.70, WIDTH 5e-3).

- **`C64_POS_LAM_MULT=0.3`** (new fork knob; default 1 = unchanged). It
  scales the cathode stress-driven LAM rate (2.78e-7 s⁻¹). The real
  cathode LAM jumps early and then plateaus at ~4–5% (CELL064 ~7%, CELL009
  ~4–5% from EFC 260–402). The model's law grows roughly linearly: it gives
  6.5% at CELL064's EFC 200, but 18.2% at CELL009's EFC 402.
  - ×0.3 gives 5.5% at EFC 402.
  - It also removes ~1–2pp of the pre-knee over-fade.
- **`C64_SI_CRACK_MULT=1e4`** (shared value 1e3). This is a stand-in for Si
  cracking driven by volume change rather than by the lithiation-gradient
  stress the Paris law uses. A partial window shrinks the model's stress
  swing, so cracks barely grow: 1.12× by EFC 402 at 1e3. The real cell
  still knees at about EFC 290.

| Si crack mult (cathode LAM ×0.3) | 1e3 (P1a) | 3e3 (P1b) | **1e4 (P1c)** | 3e4 (P1d) |
|---|---|---|---|---|
| Mean SoH gap | 15.3pp | 14.9pp | **3.4pp** | 18.6pp |
| Gap at EFC 329 / 382 / 402 | +15.5 / +38.0 / +44.6 | +15.0 / +37.1 / +42.8 | **+2.6 / +3.5 / +5.1** | −37.6 / −18.1 / −12.2 |
| Si LAM at EFC 382 (real 56.0%) | 0.3% | ~1% | **56.8%** | 87.0% |
| Crack length l/l0 at EFC 260 | 1.07 | 1.23 | 1.97 | 10.2 |
| Porosity reaches the 0.035 floor | no (0.060) | no | at ~EFC 330 | before EFC 260 |
| C/20 voltage-shape RMSE | 0.088 V | 0.086 V | **0.071 V** | 0.151 V |

**P1c is the CELL009 candidate.**
- Pre-knee gap: −2.3 to −4.0pp.
- LLI at EFC 329 / 382 / 402: 19.6 / 41.6 / 46.0%, against a real 23.0 /
  46.3 / 53.8%.
- Gr LAM: about half the DMA value.
- Knee timing is very sensitive between 1e4 and 3e4. P2 refinement, both
  with cathode LAM ×0.3:
  - P2a (1.2e4): gap at EFC 260 / 329 / 382 / 402 of −3.5 / −8.9 / −2.8 /
    +0.4pp (slightly early), mean 3.7pp, V RMSE 0.088 V.
  - P2b (1.5e4): −6.7 / −20.8 / −8.5 / −4.2pp (early), mean 7.3pp,
    V RMSE 0.116 V.
- **P1c (1e4) is kept.** The optimum is around 1.05–1.1e4, within
  cell-to-cell scatter.

The knee goes through the same route as at full window: crack growth,
crack SEI, porosity reaching the floor, then the isolation gate. Crack SEI
is most of the Si LLI (7.9% vs 2.1% from bulk SEI at EFC 402). A physical
volume-change crack driver would be a new mechanism; propose it before
implementing.

### CELL009 P3–P5: initial crack length, finer crack rate, width (2026-10-01)

All runs use cathode LAM ×0.3.

| Initial crack length | Crack mult | Pre-knee gap (EFC 79–260) | Gap at EFC 329 | Mean | V RMSE |
|---|---|---|---|---|---|
| 2e-8 m (default) | 1e4 (P1c) | −2.3 to −4.0 | +2.6 | 3.4 | 0.071 |
| 2e-8 m | **1.1e4 (P3a)** | −2.4 to −4.1 | **−3.4** | **2.7** | 0.079 |
| 2e-8 m | 1.2e4 (P2a) | −2.4 to −4.2 | −8.9 | 3.7 | 0.088 |
| 1e-8 m | 1.1 / 1.3e4 (P3b/c) | 0.0 to −2.6 | +16 (no knee) | ~15 | ~0.084 |
| 1e-8 m | 2.5e4 (P4a) | −1.6 to −3.2 | −7.5 | 3.3 | 0.091 |
| 1e-8 m | 3.5e4 (P4b) | knee at ~260 | −30.7 | 12.2 | 0.139 |
| 5e-9 m | 1.2e4 (P3d) | +1.3 to −1.7 | +18 (no knee) | 15.8 | 0.087 |
| 5e-9 m | 4e4 / 6e4 (P4c/d) | knee before 329 | −17 / −38 | 6.0 / 17.1 | 0.11 / 0.15 |

**Shorter initial cracks only help pre-knee while there is no knee.**
Growth goes as ~l^1.1, nearly exponential, so the doubling time barely
depends on l0. Restoring the knee at EFC ~330 needs a much higher rate, and
that gives similar mid-life crack area and the same pre-knee fade.
**P3a (1.1e4, default l0) is the best CELL009 setting.**

**Width test on P3a** (P5a 1e-2, P5b 2e-2):

| Width | C/3 step: one-cycle jump, at EFC | k at EFC 260 (real 0.72) | Expansion RMSE | Mean gap |
|---|---|---|---|---|
| 5e-3 | +9.2 mAh at 222 | 0.94 | 0.214 | 2.7pp |
| 1e-2 | +9.8 mAh at 209 | 0.97 | 0.254 | 2.9pp |
| 2e-2 | +9.7 mAh at 195 | 0.98 | 0.285 | 3.0pp |

A larger width moves the step earlier, but doesn't smooth it. The step is
therefore not the shape of the k(eps) transition. It is most likely the hard
cap `dv_buffered = min((1-k) dv_solid, headroom)` in
`reaction_driven_porosity.py` engaging as the structural porosity nears the
0.08 closure porosity. Width stays at 5e-3.

**The step only shows here:** in a mid-window voltage-limited discharge
(3.14 V), small polarisation changes become visible Ah. In full-window runs
the same transition coincides with the knee and is hidden.

**RPT1 reference artefact.** Model RPT1 is the formation C/20 from
initial_soc = 1.0. A C/20 RPT straight after one standard C/3 CC-CV (C/50)
recharge gives 2.402 vs 2.443 Ah (−1.7%). About 0.4% of that is one cycle
of SEI, so **~1.3pp of every CELL009 SoH gap is reference, not
degradation.** Fixed (see the next section).

### CELL009 with the conditioned RPT1 reference, and a smaller width (2026-10-01)

**RPT1 reference fix (implemented).** For partial-window runs,
`aligned_rpt_cycling.py` now runs a conditioned C/20 RPT straight after
formation (a C/20 discharge after the standard C/3 CC-CV recharge) and uses
it as the 100% reference (`RPT_SOH_REF_INDEX = 1` in the fork). Full-window
runs are unchanged.

| CELL009 (cathode LAM ×0.3, crack 1.1e4) | Gap at EFC 79 / 145 / 184 / 260 | Gap at 329 / 382 / 402 | Mean | Expansion RMSE | Hump peak EFC (real ~330) | V RMSE |
|---|---|---|---|---|---|---|
| P3a, old reference, w 5e-3 | −2.4 / −4.1 / −3.6 / −3.1 | −3.4 / −0.1 / +2.5 | 2.7 | 0.214 | – | 0.079 |
| P3a_ref, w 5e-3 | −0.9 / −2.7 / −2.2 / −2.0 | −7.0 / −1.4 / +1.6 | 2.5 | 0.259 | ~279 | 0.087 |
| P3a_ref, w 2e-3 | −0.9 / −2.7 / −2.2 / −1.9 | −4.9 / −0.5 / +2.3 | 2.2 | 0.207 | ~289 | 0.083 |
| **P3a_ref, w 1e-3** | −0.9 / −2.7 / −2.2 / −1.8 | **−2.8 / +0.3 / +2.9** | **1.9** | **0.172** | ~298 | **0.080** |

- **The reference fix closes ~1.3pp of the pre-knee gap.** Its extra
  full-depth C/20 cycle at BoL (high Si stress, ~35 MPa) also moves the
  knee ~10 EFC earlier. That matches the real start better, since the real
  cell also had RPT0 and RPT1.
- **A smaller width delays the k transition.** At EFC 329, k is 0.84 for
  1e-3 against 0.96 for 5e-3 (real 0.95). That moves the expansion hump
  later and also slightly delays the knee.
- **ADOPTED for CELL009 (user, 2026-10-01): P3a_ref with width 1e-3**,
  saved as `crack_baseline_v3_CELL009.env`. CELL009-specific
  deviations from the v3 25 °C values: cathode LAM ×0.3, crack mult 1.1e4,
  width 1e-3 (25 °C shared value 5e-3).

### Reference: plain v3 on CELL009, extended to EFC 752 (2026-10-01)

Run `v3_CELL009_15-95_extended800_25degC`: shared v3 parameters, none of the
CELL009 deviations, conditioned RPT1 reference, model-only RPTs every 50 EFC
after EFC 402 (`RPT_EXTEND_TO_EFC`).

| EFC | 260 | 330 | 403 | 453 | 503 | 553 | 602 | 652 | 752 |
|---|---|---|---|---|---|---|---|---|---|
| Model C/20 SoH | 91.1 | 89.2 | 87.2 | 85.1 | 78.0 | 62.4 | 49.2 | 41.2 | 32.6 |
| Si LAM % | – | – | 0.4 | 1.7 | 11.4 | 36.1 | 57.5 | 70.3 | 82.8 |
| Crack l/l0 | 1.07 | 1.09 | 1.11 | 1.13 | 1.15 | 1.18 | 1.26 | 1.38 | 1.71 |
| Porosity | 0.090 | 0.076 | 0.062 | 0.052 | 0.044 | 0.040 | 0.038 | 0.037 | 0.037 |

- **v3's knee is at EFC ~500–550, against ~290 real (≈1.8× later).**
- **The knee comes through slow pore clogging, not crack growth.** Bulk SEI
  plus SEI on the static initial cracks closes the pores; cracks have grown
  only 1.15× by the knee and accelerate only afterwards.
- **Cathode LAM reaches 28% by EFC 752** (real plateau ~4–5%).

So v3 does predict a delayed knee for the partial window. The missing ~250
EFC is consistent with a missing cracking route, which is what
PROPOSAL_strain_driven_cracking.md addresses.

## Log

- 2026-09-28: folder created, diagnostic + new model script set up.
  Diagnostic run complete (cracks dominate at 65.6%, genuine rate). First
  working config found (iter00, knee EFC=101.8) after fixing a
  GR_K_SEI_MULT copy-over bug. Post-knee fit (8.2pp RPT gap, 0.175
  expansion RMSE) not yet as good as even the old baseline's 1.9pp/0.10.
- 2026-09-28: measured the new baseline's own real benchmark (1.5pp RPT
  gap, 0.157 expansion RMSE, 0.291 k-gap, knee EFC=116.1) — this, not the
  old baseline's numbers, is the actual target.
- 2026-09-28: ran iter01-iter11 tuning `SI_REDIRECT_LAM_YIELD`,
  `SI_LAM_ISO_EXPONENT`, `F0`, `SI_BETA_LAM_SEI`, `SI_TAU_LAM_ISO`.
  Adopted iter08 (yield=0.20, exponent=10) as the new default in the
  script: beats the new baseline on knee timing and voltage-shape RMSE,
  close on RPT-gap and k-gap, but expansion RMSE regressed (0.281 vs
  target 0.157) as a side effect of the post-knee LAM increase needed to
  close the RPT gap — see "Adopted config" section above for full detail,
  root-cause hypothesis, and untried next steps (NEG_POROSITY_FLOOR
  tuning, queued but blocked by a persistent tool error this session).
  **Study paused here, not fully complete** — expansion-shape fit remains
  an open gap against the new-baseline benchmark.
- 2026-09-28: tuning campaign 2 (iter17–26). Fixed two core bugs:
  bulk-OCP domain, and per-phase deformation scoping. Graphite isolation LAM
  added. Adopted iter25b as `25degC_baseline1` (the script defaults) and kept
  iter26 as `25degC_baseline2`. Deleted all other plots and logs.
- 2026-09-29: found that Si cracks are static, which resolves iter13; 45 °C
  from baseline1 knees too early (~48). Crack-growth sweep on the
  high_temp_45C fork: `SI_CRACK_MULT=1.5e3` gives a 25 °C onset ~106 and a
  45 °C knee ~310 (`SI_CRACK_EAC=38000`). The reduced-cracking mechanism now
  works. Next: fit the 25 °C post-knee/LAM/LLI on this config.
- 2026-09-29: fitted the cracking config to the real 25 °C data. Candidate
  crack baseline: Si 0.15, Gr 0.013, F0 0.7. 25 °C RPT gap 2.5pp, expansion
  0.121, V 0.069; 45 °C knee ~265–285 from the same set.
- 2026-09-29: added the 45 °C CELL017 SoH overlay. Saved `crack_baseline_v1`
  (one set for 25 and 45 °C). Found that the initial crack length is only a
  timing knob. Added the opt-in stress-LAM damper option for testing.
- 2026-09-30: found and fixed the composite SEI-on-cracks charge-balance bug
  (CHANGES item 14). Tried θ-trapping and SEI suppression, both reverted.
  Retuned R1–R3 on the fixed core and adopted R3c as `crack_baseline_v2`
  (25 °C mean gap 2.0pp, 45 °C 1.2pp, one shared set).
- 2026-09-30: R4 width test (kept R3c). Tidied the folder: removed all logs,
  shell drivers and non-baseline plots.
- 2026-09-30: C/50 polarisation test; θ re-added and re-tested; aligned-RPT
  driver added (`aligned_rpt_cycling.py`). Adopted R7 / R7a as
  `crack_baseline_v3` (θ = 1, Gr yield 0.025, `GR_YIELD_EAC` 56.5 kJ/mol).
- 2026-09-30: stage 3. Built the CELL009/CELL026 partial-window data; ran v3
  with their real protocols. No model knee for either window.
- 2026-10-01: CELL009 tuning (P1-P5, plus the reference fix). Adopted
  P3a_ref_w1e-3: cathode LAM x0.3, crack mult 1.1e4, width 1e-3, conditioned
  RPT1 reference. Proposal for strain-driven cracking written
  (`PROPOSAL_strain_driven_cracking.md`).
