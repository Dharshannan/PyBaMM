# CELL064 conditions test matrix

## Baseline parameter reference (read this first — do not confuse these three)

**Correction, 2026-09-26**: earlier text in this file and in-session analysis
briefly compared a 45 degC run against `degradation_test_matrix/
cell064_overall_summary_result.png` (the untagged plot) and called it "the
baseline". That file is the **old baseline** (see below) — a pre-LAM/LLI-
split-tuning config, NOT the dma-folder's best-fit-to-experiment reference.
The correct comparison target for shape/timing questions is the **new
baseline** (`dma/new_baseline_*` files). This section exists so that
mistake doesn't happen again — three genuinely distinct configs are in play
in this project, and they must not be conflated:

### 1. Old baseline ("true old baseline", 2026-09-14 config)

The canonical pre-session config, from BEFORE the whole Si-LAM/porosity-
isolation/redirect-to-LAM investigation started. Reconstructed exactly by
`dma/build_true_old_baseline.py` (runs `cell064_degradation_fit.py` with
**zero env var overrides** — i.e. every value below is the shared script's
own hardcoded default). Produces `degradation_test_matrix/
cell064_overall_summary_result.png` (untagged) when run with no
`C64_OUT_TAG`, or the `*_true_old_baseline*` files when run via the build
script itself (which also bumps `NEG_POROSITY_FLOOR` 0.023→0.028 and
`MAX_CYCLES` 250→350 purely for solver stability/reach, not a recipe
change). **No porosity-isolation mechanism at all** (`SI_BETA_LAM_ISO=0`,
`SI_REDIRECT_TO_LAM=0`) — LAM_Si is stress-driven only, no knee-coupled
mechanism. This is the run whose reversible-expansion panel is flat-then-
sharp-peak (matches real data reasonably) but whose SoH/LAM fit is the
one the whole later investigation was motivated to improve on.

| Parameter | Value |
|---|---|
| SEI option | `"ec reaction limited"` (both Gr and Si) |
| SI_LAM_OPTION | `"stress and reaction-driven"` (no porosity isolation) |
| SI_BETA_LAM_ISO / SI_REDIRECT_TO_LAM | `0.0` / `0` (isolation mechanism OFF) |
| SI_LAM_ISO_EXPONENT / SI_TAU_LAM_ISO | `4.0` / `1.0` (irrelevant, mechanism off) |
| SI_OCP_AGING_DEFORM / SI_VOLUME_CHANGE_AGING_DEFORM | `0` / `0` (both OFF) |
| SI_CRIT_STRESS / GR_CRIT_STRESS | `2.2e8` / `3.0e7` Pa |
| SI_LAM_EXP (m_LAM) | `2.5` |
| NEG_POROSITY_FLOOR | `0.023` (0.028 in the stability-patched rebuild) |
| EXPONENT_MAX_SEI | `100.0` |
| TIMESCALE_MULT / LAM_PROP_MULT | `1.0` / `0.6` |
| SI_LAM_PROP / GR_LAM_PROP | `3.078e-8` / `4.634e-7` s⁻¹ |
| SI_CRACK_MULT / GR_CRACK_MULT | `0.1` / `0.1` |
| SI_BETA_LAM_SEI | `7e-6` |
| F0 / WIDTH | `0.7` / `0.05` |
| NEG_POROSITY_BOL | `0.25` |
| K_SEI_MULT (both Gr/Si base) | `1.0e-3` |
| SI_DIFFUSIVITY_MULT / SI_PARTICLE_RADIUS_MULT / SI_EXCHANGE_CURRENT_MULT | `1.0` / `1.0` / `1.0` (no retune) |
| MAX_CYCLES / SOH_FLOOR / EFC_LIMIT | `250` / `45.0` / `inf` (build script: `350`/unchanged/unchanged) |

### 2. New baseline (dma-folder reference, best fit to experiment)

The current adopted, best-fit-to-real-data config, per `dma/
build_new_baseline_curves.py` (the authoritative source — read directly,
not reconstructed from memory) and `dma/FINDINGS_2026-09-20.md`. This is
what `../degradation_test_matrix/cell064_degradation_fit.py` should be run
WITH these overrides to reproduce (plots: `dma/new_baseline_*`,
`degradation_test_matrix/cell064_overall_summary_result_new_baseline_*.png`).
**Porosity isolation is ON** here (unlike old baseline) — this is the
config `CELL064_final_parameterisation.md` documents as current-best.

| Parameter | Value | Changed from old baseline? |
|---|---|---|
| SEI option | `"ec reaction limited"` (both) | No (unchanged) |
| SI_REDIRECT_TO_LAM | `1` | **Yes** — isolation mechanism ON |
| SI_BETA_LAM_ISO | `1.0` | **Yes** (was 0.0) |
| SI_TAU_LAM_ISO | `2e8` | **Yes** (was 1.0) |
| SI_REDIRECT_LAM_YIELD | `0.02` | new knob (n/a in old baseline) |
| SI_LAM_ISO_EXPONENT | `3.0` | new knob (n/a in old baseline) |
| SI_OCP_AGING_DEFORM | `1` | **Yes** (was 0) |
| SI_VOLUME_CHANGE_AGING_DEFORM | `1` | **Yes** (was 0) |
| GR_REDIRECT_TO_LAM / GR_REDIRECT_LAM_YIELD | `0` / `0.0` | unchanged (graphite never redirects) |
| NEG_POROSITY_FLOOR | `0.035` | Yes (was 0.023) |
| EXPONENT_MAX_SEI | `70` | Yes (was 100.0) |
| SI_BETA_LAM_SEI | `1e-7` | Yes (was 7e-6) |
| SI_PARTICLE_RADIUS_MULT | `0.2` | **Yes** (was 1.0) — the radius/diffusivity retune |
| SI_DIFFUSIVITY_MULT | `3` | **Yes** (was 1.0) |
| NEG_POROSITY_BOL | `0.17` | Yes (was 0.25) |
| SI_K_SEI_MULT | `3e-4` | Yes (replaces the old `K_SEI_MULT=1e-3` scheme) |
| GR_K_SEI_MULT | `6e-6` | Yes |
| F0 | `0.7` | unchanged |
| WIDTH | `0.03` | Yes (was 0.05) — narrowed to preserve roughly the same width/porosity-margin ratio after NEG_POROSITY_BOL/FLOOR shrank the margin (0.25−0.023=0.227 old vs. 0.17−0.035=0.135 new; 0.05/0.227≈0.22 ≈ 0.03/0.135≈0.22) |
| MAX_CYCLES / SOH_FLOOR / EFC_LIMIT | `700` / `30` / `205` | Yes (extended reach + explicit EFC cap) |
| All else (SI_CRIT_STRESS, SI_LAM_PROP, GR_LAM_PROP, SI_CRACK_MULT, etc.) | same as old baseline | No |

### 3. Conditions-test-matrix high-temperature config (this investigation, current state)

`high_temp_45C/cell064_degradation_fit_45C.py` — a **standalone fork**
starting from the new baseline recipe above (same isolation/retune values
inherited as its own hardcoded defaults), with the operating temperature
changed to 318.15 K AND, as of 2026-09-25/26, a deliberate **mechanism
change** for Si's SEI growth (needed to reproduce the real high-T knee
delay — see "Knee-timing investigation" section below for the full
derivation). This is NOT just a condition change on top of the new
baseline — the SEI option itself differs for Si.

| Parameter | New baseline value | High-T (45 degC) fork value | Why changed |
|---|---|---|---|
| Ambient/initial temperature | 298.15 K | **318.15 K** | the condition being tested |
| SEI option | `"ec reaction limited"` (both) | `(("ec reaction limited", "solvent-diffusion limited"), "none")` — **Si only** switched | "ec reaction limited" cannot reproduce a delayed (not accelerated) knee at higher T under any rate-constant tuning tried (see investigation section) |
| Secondary: SEI solvent diffusivity | n/a (option unused) | `2.5e-22 * SOLVENT_MULT`, `SOLVENT_MULT=34` | calibrated so this fork's OWN 25 degC prediction matches the real anchor knee EFC~101 |
| Secondary: SEI growth activation energy | `38000` J/mol (shared default, both phases) | `-37500` J/mol (Si only; graphite stays `38000`) | calibrated so 45 degC knee lands at the real target EFC~300 (deliberately negative — see physical caveat in investigation section) |
| SI_CRACK_MULT | `0.1` | `1.0` | restored genuine (unsuppressed) cracking so LAM_Si gets a gradual, cracking-driven contribution from BOL rather than only the isolation step |
| SI_LAM_ISO_EXPONENT | `3.0` | `15.0` (as of 2026-09-26) | 3.0 gave a knee EFC number near the 45 degC target but NO visible knee shape (smooth curve); 15.0 restores a genuine visible bend at 45 degC while leaving the 25 degC shape/timing essentially unchanged |
| Everything else | — | inherited unchanged from new baseline | not part of this investigation |

**Update 2026-09-26 (later same day) — expansion shape mostly resolved**:
the "front-loaded vs. back-loaded porosity decline" theory above was based
on comparing against the WRONG reference plot (the old baseline, not the
new baseline — see the correction at the top of this file). Against the
correct new-baseline reference (`dma/new_baseline_*`, which itself shows a
gradual rise from EFC=0, peaking essentially AT the knee, ratio
peak-EFC/knee-EFC≈1.0-1.14), the actual finding is simpler: `WIDTH=0.03`
made the k-transition (hence the reversible-expansion peak) arrive too
EARLY relative to knee timing (~60-80% of knee-EFC, not ~100-114%).
Narrowing `WIDTH` to **0.001** (a further 30x reduction) fixes this — k now
stays flat near `f0` through ~180 EFC and peaks at ~90-100% of knee-EFC,
matching the reference's timing. This is now the fork's default (see its
own code comment for the derivation).

**Still open: hump AMPLITUDE.** With `WIDTH=0.001`, the peak magnitude is
~1.49x (normalised) — larger than the new baseline's own peak (~1.32
model / ~1.2 real). Peak amplitude turned out to be governed almost
entirely by `F0` (`k` runs from `f0` at BOL to `1` at the floor, so a
smaller `f0` gives a bigger swing) and is essentially independent of
`WIDTH` (iter06 at width=0.03 and iter18 at width=0.001 both peaked at
~1.49-1.50, confirming amplitude and timing are governed by separate
knobs). Tested `F0=0.85`: peak drops to ~1.15 (now undershooting the
~1.2-1.32 target). True value is likely in the **0.75-0.82** range — not
yet pinned down; a run at `F0=0.79` was queued but stopped before
completion per instruction, session paused here. `F0_BASELINE` in the fork
is left at the new-baseline's own `0.7` until this is resolved (do not
change without re-checking hump amplitude against `dma/new_baseline_ocp`
as the reference).

### High-temperature (45 degC) fork — FINAL applied parameters (2026-09-28)

Confirmed working and memory-safe (see engineering note below). This
supersedes every earlier "latest"/"converged" table above — those are kept
for derivation history only.

| Parameter | Value | Note |
|---|---|---|
| SEI option | `(("ec reaction limited", "solvent-diffusion limited"), "none")` | Si only switched |
| `SOLVENT_MULT` | `5.0` | sets `Secondary: SEI solvent diffusivity = 2.5e-22*5` |
| `SI_ESEI` (Secondary: SEI growth activation energy) | `38000.0` J/mol | plain positive/standard — the earlier negative-Ea mechanism was abandoned |
| `GR_ESEI` (Primary: SEI growth activation energy) | `38000.0` J/mol | unchanged from default |
| `SI_CRIT_STRESS` | `6e6` Pa | down from literature default 2.2e8 — real Si stress here is only ~1.3e7 Pa (confirmed via `C64_DIAG_STRESS`), so the default threshold left stress-driven LAM_Si permanently negligible |
| `SI_LAM_EAC` (new) | `40000.0` J/mol | positive Arrhenius boost on `beta_LAM` (stress-driven LAM_Si rate constant), which has no T-dependence by default in this recipe |
| `SI_CRACK_MULT` | `1.0` | genuine (unsuppressed) cracking, was 0.1 |
| `SI_LAM_ISO_EXPONENT` | `15.0` | was 3.0 (new baseline) — sharpens the isolation-gate transition so a genuine knee SHAPE (not just a matching EFC number) emerges |
| `WIDTH` (pore buffering transition width) | `0.0003` | was 0.03 (new baseline) — compresses the reversible-expansion peak's timing to arrive at ~knee EFC instead of ~60-80% of the way there |
| `F0` (transmitted fraction plateau) | `0.85` | was 0.7 (new baseline) — brings the expansion peak AMPLITUDE down toward the real target (~1.06-1.15 achieved vs ~1.2-1.32 target; still slightly under, not further tuned) |

**Result**: knee SoH-91% crossing at EFC~310.7 (target ~300-310, "looks ok" per
explicit confirmation), genuine visible knee bend (not just a matching
number), reversible expansion shows a real initial dip (1.0→~0.965 by
EFC~200) then a sharp peak (~1.09) right near the knee (~EFC 330) — same
qualitative shape as the real 45 degC/new-baseline reference, though the
peak amplitude is still somewhat undershooting the ~1.2-1.32 target.

**Engineering note**: this fork had NO `output_variables` memory
optimization at all (unlike the shared script) until 2026-09-28 — a plain
400-500 cycle verification run was OOM-killed. Ported the same
`_CORE_OUTPUT_VARIABLES`/`_RESTRICT_OUTPUT_VARIABLES` mechanism from the
shared script into this fork (applied to both the main and retry solvers),
plus a new `C64_DIAG_CRACK_SPLIT=1` flag for the study_si_deg
crack-vs-bulk-SEI diagnostic. Verified fixed: a 450-EFC run (550 cycles)
now completes without incident.

**Not yet resolved** (open, lower priority): expansion peak amplitude
still ~10-15% below the real target; `F0` between 0.85 and 0.7 (e.g.
~0.78-0.80) would likely close this gap but wasn't pinned down further.

### Correction to the FINAL table above (2026-09-29): crack terms are inert, cracks are static

**Omission in the table:** the fork's code also sets `SI_CRACK_EAC=38000`
(Si cracking rate ÷2.6 at 45 °C, `(1/T − 1/298.15)` convention), which
the table doesn't list.

**Neither crack-rate term does anything in this recipe.** Si (and graphite)
crack length stays at its initial 2e-8 m and roughness at ~2.908 at both 25
and 45 °C (direct crack-variable check). Crack GROWTH is negligible, so
`SI_CRACK_MULT` and `SI_CRACK_EAC` are both inert. `SI_LAM_EAC` is
essentially inert too, since stress-driven LAM barely activates. The FINAL
recipe's 45 °C knee (~310) therefore comes from `SOLVENT_MULT=5`, a rate
recalibration applied at every temperature, not from reduced high-T
cracking.

**Making the "Alternative theory" mechanism real** (full detail in
`study_si_deg/FINDINGS.md`, "Crack-growth sweep"): keep this fork's recipe
and `SOLVENT_MULT=5`, and raise `SI_CRACK_MULT` to **1.5e3** so cracks
actually grow (×~15 by the knee at 25 °C):

| Condition | Onset (SoH-91% crossing) | Steepest-descent knee |
|---|---|---|
| 25 °C | 106 (real ~101) | 93 |
| 45 °C, `SI_CRACK_EAC=0` | 213 | 263 |
| 45 °C, `SI_CRACK_EAC=38000` | 263 | **310** (real ~300) |

Lower stress at 45 °C (~16–20 vs ~43 MPa pre-knee, from faster Si
diffusion; Paris m = 2.2) alone doubles the knee EFC. This is the first
config in this investigation that gets both temperature anchors from one
parameter set through a physical cracking mechanism. The 25 °C post-knee,
LAM and LLI fit on this config is still to be done.

New fork switches, both opt-in with defaults unchanged:
- `C64_T_AMBIENT_K` sets the operating temperature. At 298.15 K the fork
  overlays and scores against the real 25 °C data.
- `C64_DIAG_EXTRA_VARS` is a `|`-separated list of extra output variables.

### Narrow-window (partial SoC) real data anchors (2026-09-27)

Per instruction, real experimental data exists for this condition (SoH only
— no full DMA split yet): **knee at EFC~298**, and a reversible-expansion
pattern qualitatively similar to the 45 degC condition (rises to a hump
near the knee) **but without an initial dip** — expansion stays roughly
flat/constant pre-knee, then rises to the hump, unlike 45 degC's
dip-then-hump.

**First attempt result (2026-09-27, this session, NOT yet calibrated
against this target)**: ran the narrow-window fork at 25 degC using the
latest 45 degC-tuned parameters as-is (`SOLVENT_MULT=5`, `SI_ESEI=38000`
positive, `SI_CRIT_STRESS=6e6`, `SI_LAM_EAC=40000`, `WIDTH=0.0003`,
`F0=0.85`) up to EFC~427 (637 cycles) — **no knee emerged at all** in this
range (LAM_Si stayed ~0%, LAM_positive dominated at ~21%, expansion
declined monotonically with no dip, no flat plateau, no hump). This
confirms the current parameterization is miscalibrated for this condition,
not merely "not run far enough": since the real target (298) is LOWER than
where this run still shows nothing, `SOLVENT_MULT` likely needs to
INCREASE (faster porosity closure) for this specific window, unlike the
wide-window case. Also, since the real data shows NO initial dip here, the
stress-driven LAM bump (`SI_CRIT_STRESS`/`SI_LAM_EAC`) that was added
specifically to produce the 45 degC dip probably should NOT be part of
this condition's recipe at all — just the plain porosity/isolation
mechanism, with `WIDTH`/`F0` retuned for a flat-then-hump (not
dip-then-hump) shape.

**Not yet done** (paused here per session length): actually retuning
`SOLVENT_MULT` (and dropping the stress-LAM bump) for this condition to
hit EFC~298 with the correct flat-then-hump expansion shape.

**Engineering note for whoever resumes this**: the narrow-voltage-window
fork (`conditions_test_matrix/narrow_voltage_window/
cell064_degradation_fit_narrow_window.py`) needs ~4-13x more CYCLES than
the wide-window fork to reach a comparable EFC range (narrower window moves
less capacity per cycle) — a single-process run at that cycle count
produced an unmanageable ~55GB pickled solution and OOM'd even a
chunked/resumable retry. The fix that worked: `run_lean_narrow25C.py` (in
`degradation_test_matrix/sweeps/porosity_isolation_calib/`) splits each
50-cycle batch into two solve calls — 49 ageing-only cycles with
`store_first_last=True` (first/last sample per step only, ~50x less data)
and the 1 RPT cycle (always last in the batch, since
`RPT_INTERVAL==BATCH_SIZE==50`) with full within-step resolution (needed
for the RPT voltage-curve panel). Verified to reproduce identical SoH
numbers to a plain single-solve-per-batch run at matched cycle counts, and
completed 637 cycles in a single process without OOM. Reuse this pattern
(or port `store_first_last`/`_CORE_OUTPUT_VARIABLES` properly into the
narrow-window fork itself, which currently has neither) for any further
narrow-window work. The shared script (`cell064_degradation_fit.py`) also
now has an opt-in `C64_STORE_FIRST_LAST=1` flag wired in (default off) for
the same purpose, though it applies uniformly to every step (including
RPTs) — fine for a pure SoH/EFC check, not for anything needing the
voltage-curve panel.

### Alternative theory (2026-09-26, documented for future work — NOT implemented)

**Hypothesis**: the real high-T knee delay may come from reduced Si
**cracking** at higher temperature, not from an exotic negative apparent
activation energy on bulk SEI growth. Physical chain: higher T -> lower
overpotential and faster (genuinely Arrhenius) particle diffusivity ->
smaller concentration gradients within Si particles -> less
diffusion-induced/electrochemical stress -> less driving force for crack
propagation (also consistent with lithiated Si becoming more ductile at
modestly elevated T — a real phenomenon, already noted in this fork's own
`SI_CRACK_EAC` comment from the 2026-09-25 crack-Arrhenius test) -> less
crack SURFACE AREA -> less "SEI on cracks" growth -> **less total porosity
closure** at high T -> knee delayed. At low T, the reverse: more cracking,
more SEI-on-cracks growth, more porosity closure, earlier knee (matching
the real ~101 EFC anchor).

**Why the earlier crack-Arrhenius test's "zero effect" finding does NOT
rule this out**: that test (2026-09-25, `SI_CRACK_EAC`) checked whether
crack-rate T-dependence affects **LAM_Si** via `a_j_sei` in
`loss_active_material.py` — correctly zero, since neither the reaction-
driven nor porosity-isolation LAM term ever reads the "SEI on cracks"
current (see "silicon cracking is NOT coupled to LAM_Si" elsewhere in this
file). But cracking has a **second, separate** pathway that test never
isolated: verified directly in `reaction_driven_porosity.py` (line ~164),
overall electrode porosity decline is
`L_sei + L_pl + L_dead + L_sei_cr*(roughness-1)` — i.e. SEI-on-cracks
thickness, scaled by crack-driven roughness, IS a genuine additive
contributor to **total porosity closure** (hence to knee timing via the
isolation gate's `headroom_frac`), completely independently of the LAM
pathway. A genuine crack-rate Arrhenius effect on knee TIMING was never
actually tested in isolation — the 2026-09-25 test's conclusion ("zero
effect") was specifically about LAM_Si, and got generalised too broadly at
the time.

**Two implementation paths for the future** (neither started):
1. **Quick/approximate**: keep treating this as an effective reduced bulk
   SEI rate at higher T (what the current negative-`E_sei` mechanism
   already does numerically) — cheap, already working, but not tied to a
   real mechanism.
2. **Proper**: give the EXISTING `silicon_cracking_rate_Ai2020`-style
   function (already wired via `SI_CRACK_EAC`, currently defaulted to 0 in
   the shared baseline) a genuine positive `Eac_cr` (crack rate decreases
   with T, per the convention documented in this fork's own code), keep
   Si's SEI as `"ec reaction limited"` (no need for the speculative
   negative `E_sei` at all), and **re-fit both temperatures**: 25 degC
   against real data (as the new baseline already is, but now also
   checking the crack-rate/SEI-on-cracks contribution reproduces the real
   LAM split), and 45 degC to confirm the delay emerges from genuinely
   reduced cracking rather than being imposed via an ad-hoc rate constant.
   This would let the whole 45 degC fork drop the negative-`E_sei`
   mechanism entirely, replacing it with something with a real physical
   handle (crack propagation's known stress/ductility temperature
   dependence) instead of a back-calculated, phenomenological number.

Per instruction, no changes or sweeps for this were made now — this is a
documented direction for a future session.

**Summary of what's confirmed vs. still open, end of 2026-09-26 session:**
- Knee timing at both 25 degC (~101.8, matches anchor) and 45 degC (~352,
  a deliberate ~15% higher than the original ~300 target — user confirmed
  "the knee EFC looks ok here", so not re-tuned back down) — CONFIRMED,
  with genuine visible knee SHAPE at both temperatures (not just a matching
  number — this was itself a correction of an earlier mistake this
  session, see the mechanism-search writeup above).
- Reversible-expansion peak TIMING relative to knee (flat pre-knee, sharp
  peak at/near the knee) — CONFIRMED via `WIDTH=0.001`.
- Reversible-expansion peak AMPLITUDE (~1.2-1.32 target) — OPEN, `F0`
  somewhere in 0.75-0.82 likely closes this, not yet found.
- Silicon cracking's coupling to LAM_Si — investigated and found to be
  architecturally ABSENT (crack current never reaches either LAM term, see
  below) — not fixed (would need a core-code change), deliberately not
  pursued further after stress-driven LAM_Si was shown to be the more
  promising (and ultimately sufficient, via WIDTH/F0) lever for the
  expansion-shape questions instead.

Tracks operating-condition variants of the CELL064 degradation fit: same
tuned degradation recipe as the baseline
(`../degradation_test_matrix/cell064_degradation_fit.py`, current adopted
config per `../CELL064_final_parameterisation.md`), with only the stated
operating condition changed. The point is to see how the *same* mechanism
recipe extrapolates outside the conditions it was actually fitted to — this
is exploratory, not a fit, since no real CELL064 data exists at these
altered conditions.

Each condition is a **standalone forked copy** of the full baseline script
(not an import + env-var override), so it can be tuned/run independently
without touching the shared baseline. See each fork's own docstring header
for exactly what was changed and why. If a genuine bug fix or recipe change
is made in one copy, it needs to be ported back to the baseline and the
other fork(s) manually — nothing here shares state automatically.

**Degradation recipe active in both forks** (identical `MODEL_OPTIONS_BASE`
in both, unmodified from the baseline except where noted): SEI = "ec
reaction limited" (Gr + Si, each own kinetic rate constant), SEI porosity
change on, SEI on cracks on, SEI film resistance "distributed", particle
mechanics "swelling and cracking" (Gr + Si) / "swelling only" (positive),
LAM: graphite stress-driven only, silicon stress- **and** reaction-driven
**and, as of 2026-09-25, porosity isolation** (see below), pore buffering
on ("physical"), thermal isothermal. OCP aging deformation, volume-change
aging deformation, and the active-material expansion residual (item 30) are
all off. All particle diffusivities (Gr/Si/NMC) and both electrolyte
transport parameters (diffusivity, conductivity) carry genuine Arrhenius
temperature dependence, verified directly from
`src/pybamm/input/parameters/lithium_ion/si_gr_expansion.py` — real
activation energies (17–48 kJ/mol range), not zeroed out, and the
multiplier wrapper functions in `cell064_degradation_fit.py` pass `T`
through to those functions unmodified (current multipliers are 1.0, i.e.
no-ops, in both forks).

**2026-09-25 update — Si porosity-gated isolation enabled:** both forks
were changed from `SI_BETA_LAM_ISO=0.0` (off) to `SI_BETA_LAM_ISO=1.0` with
`SI_REDIRECT_TO_LAM=1` (selects the "redirect-to-LAM" branch, the form
`dma/FINDINGS_2026-09-20.md`'s Finding 2 identifies as "the one used for
CELL064's Si"). This makes `SI_LAM_OPTION = "stress and reaction-driven and
porosity isolation"`. See "Findings" below — this changed the results
drastically, and not in a way that looks like a genuine per-condition
effect (see the caveat there before reading too much into it).

Each fork's `plot_all()` was also edited to not require real experimental
data (see "No real-data overlay" below) — that's the other intentional
difference from the baseline script, beyond the one condition variable each
is named for.

## Conditions

| # | Folder | Condition | Temperature | Voltage window | Everything else | Status |
|---|---|---|---|---|---|---|
| 0 | `../degradation_test_matrix/` | Baseline (reference) | 298.15 K (25 degC), isothermal | 2.5–4.195 V | — | Fitted, see `../CELL064_final_parameterisation.md` |
| 1 | `high_temp_45C/` | High temperature | **318.15 K (45 degC), isothermal** | 2.5–4.195 V (unchanged) | Identical to baseline | Rerun complete (2026-09-25, isolation ON): hit 45% SoH floor at EFC~37, 49/50 cycles |
| 2 | `narrow_voltage_window/` | Narrow (reduced-SoC) voltage window | 298.15 K (unchanged) | **3.15–4.2 V** | Identical to baseline | Rerun complete (2026-09-25, isolation ON): hit 45% SoH floor at EFC~49, 98/100 cycles |

Earlier (isolation-OFF) run numbers: 45 degC reached 206/250 cycles, final
SoH 51.5% @ EFC~159; narrow window reached 245/250 cycles, final SoH 53.3%
@ EFC~136. See "Findings" for both sets side by side.

## Running

Each folder's script is self-contained and runnable directly:

```
python conditions_test_matrix/high_temp_45C/cell064_degradation_fit_45C.py
python conditions_test_matrix/narrow_voltage_window/cell064_degradation_fit_narrow_window.py
```

Plots save into the same folder as the script (`SCRIPT_DIR`-relative, same
convention as the baseline), e.g. `cell064_degradation_fit_result.png`,
`cell064_overall_summary_result.png`, etc.

**No real-data overlay (yet):** collaborators are running the physical
experiment for each of these conditions and will provide real data in a few
days. Until then, both forks' `plot_all()` uses empty (but correctly
columned) placeholder DataFrames instead of `load_experimental_*()`, so
every plot renders model-only curves — no real-data scatter/dashed lines,
and the `score_*()` fit-quality functions that need real data are skipped
(only `score_knee()`, a pure model-self diagnostic, still prints). Wiring
the real CSVs in later is a matter of swapping the empty-DataFrame block in
each fork's `plot_all()` back to the `load_experimental_*()` calls (see the
comment left there) — the file paths/column schemas are unchanged from the
baseline's `experimental_data/` convention, whenever the new CSVs are
placed there or a condition-specific `experimental_data/` folder.

Both forks currently inherit the baseline's full cycle budget
(`MAX_TOTAL_CYCLES=250`, `SOH_TERMINATION_PERCENT=45`, etc., all still
env-var overridable via `C64_*` exactly as in the baseline script) — expect
a run to take a similar amount of time/memory as a baseline run.

## Findings

Both runs completed cleanly (`reached MAX_TOTAL_CYCLES`, no solver crash) and
produced all diagnostic plots (`cell064_degradation_fit_result.png`,
`cell064_expansion_fit_result.png`, `cell064_overall_summary_result.png`,
`stoichiometry_evolution.png`, `cell064_thermal_diagnostics_result.png` in
each folder). These are **initial, model-only observations** — exploratory,
not validated against real data, and the runs stopped at the 250-cycle
budget well short of the EFC range (~170–200) the baseline was actually
fitted/validated against, so some of this may simply be "not run far enough
past the knee yet" rather than a genuine condition-driven effect. Treat
accordingly; revisit once collaborators' real data lands.

### Cross-condition sanity check: temperature dependence is real

BOL internal resistance differs exactly as Arrhenius kinetics would predict:
**~80 mOhm at 45 degC** vs **~140–150 mOhm at 25 degC** (narrow-window run,
same temperature as baseline). This confirms CELL064's own particle/
electrolyte transport parameters carry genuine temperature dependence (unlike
the older `si_gr_expansion` precursor model, where the conditions-test-matrix
plan for *that* model found cracking-rate/SEI-kinetics activation energies
hard-set to zero) — a real, physically meaningful T-effect is entering both
runs, not just an isothermal label change.

### High temperature (45 degC)

- **Knee timing barely moves**: knee EFC = 112.1 (vs. the real 25 degC
  anchor of ~101, i.e. only +10.7 EFC later) despite the 20 degC jump — the
  onset of the knee is nearly insensitive to temperature in this recipe.
- **But post-knee collapse is much steeper**: SoH falls from ~90% to 51.5%
  in only ~50 EFC after the knee (EFC 110→160), a faster post-knee collapse
  than the baseline shows over a comparable EFC span.
- **LLI totally dominates LAM_Si** in the degradation-mode split (~50% LLI
  vs. ~16% LAM_Si by EFC 160) — the opposite balance from the real baseline
  (LAM_Si dominant, ~58–80% by EFC 171–197). Caveat above applies: this run
  only reaches ~50 EFC past its own knee, well short of the real RPT4/RPT5
  comparison points, so this may partly reflect "too early post-knee" rather
  than a genuine temperature-driven shift toward LLI. Worth re-running with
  a larger `C64_MAX_CYCLES` before drawing a firm conclusion.
- **Reversible expansion overshoots then collapses below its BOL value**:
  peaks at ~1.43x early-life value right at the knee (EFC~110), then falls
  to ~0.80x by EFC~160 — i.e. ends up *below* where it started. Same pattern
  seen in the narrow-window run (see below), so this looks like a shared
  structural behavior (LAM-driven collapse of the particle-swelling signal
  once the electrode is sufficiently isolated) rather than something
  specific to high temperature.
- **Si stoichiometric window stays reasonably healthy** post-knee (narrows
  to ~0.68–0.86 by EFC~159, not the severe ~9-300x collapse the pre-fix
  baseline used to show) — the diffusivity/particle-radius retune that fixed
  this for the 25 degC baseline generalizes to 45 degC without a new
  pathology.
- **Real post-knee self-heating would likely exceed what this isothermal run
  captures**: the (unfed-back) quasi-steady temperature estimate reaches
  ~15 K above the 45 degC ambient by EFC~150, driven by internal resistance
  climbing from ~80 to ~350 mOhm. A real 45 degC cell would likely run
  hotter than 45 degC post-knee, self-accelerating further — not modelled
  here (would need `THERMAL_OPTION="lumped"`, a separate, bigger change).

### Narrow voltage window (3.15–4.2 V)

- **Knee arrives earlier in EFC terms**: knee EFC = 93.7 (vs. anchor ~101,
  i.e. −7.7 EFC). **Important caveat**: EFC is defined as throughput
  capacity / (2 × nominal capacity), and a narrower voltage window inherently
  moves less capacity per cycle — so the same wall-clock/cycle-count aging
  maps onto a *smaller* EFC value purely as a bookkeeping consequence of the
  normalisation, not necessarily because the cell degrades "faster" in any
  absolute sense. The SoH-91% crossing at EFC=40.6 (vs. anchor ~101) should
  be read with the same caveat.
- **A small non-monotonic wiggle appears in the SoH curve** around EFC
  75–90, just before the knee (a slight plateau/tick-up, visible in
  `cell064_overall_summary_result.png`'s SoH panel) — not obviously
  physical; worth a closer look (possible batch-boundary or solver-retry
  artifact) before reading too much into the exact pre-knee shape.
- **Same LLI-dominant split** as the 45 degC run (~40% LLI vs. ~12% LAM_Si
  by EFC 136) — same caveat as above about not yet being far enough past
  the knee for a fair comparison to real RPT4/RPT5.
- **Same expansion overshoot-then-collapse pattern**: peaks ~1.34x around
  EFC~90 (at the knee), falls to ~0.78x by EFC~136 — reinforces that this is
  a shared model behavior, not condition-specific.
- **BOL Si stoichiometric window is naturally narrower** (~0.77–0.98 vs. the
  full 0–1 range), as physically expected from cycling a smaller voltage/SoC
  swing — degrades further post-knee but not catastrophically.
- Internal resistance and estimated self-heating rise on the same order as
  the 45 degC run (~140→365 mOhm; ~16 K estimated rise above ambient by
  EFC~135), despite this condition being at the baseline's own 25 degC —
  consistent with the same post-knee resistance-growth mechanism, just
  triggered slightly earlier (EFC~93 vs ~110).

### Rerun with Si porosity-gated isolation enabled (2026-09-25)

Both conditions were rerun with `SI_BETA_LAM_ISO=1.0`/`SI_REDIRECT_TO_LAM=1`
(previously 0.0/off — see the recipe note above). Result: **both conditions
collapse to the 45% SoH floor almost immediately and almost identically**,
regardless of temperature or voltage window:

| | Knee EFC (vs. anchor ~101) | SoH at knee | 45%-floor hit at | Cycles before floor |
|---|---|---|---|---|
| 45 degC | 31.2 (−70.2) | 52.5% | EFC~37 | 49/50 |
| Narrow window | 29.9 (−71.5) | 69.4% | EFC~49 | 98/100 |

In both, LAM_Si races to ~100% by EFC~35–40 (vs. topping out at ~16% and
~12% respectively in the isolation-off runs), and reversible expansion
crashes to ~0.2–0.4x its BOL value.

**This is almost certainly not a genuine temperature or voltage-window
effect** — it's the same catastrophic collapse, at essentially the same EFC,
under two very different conditions (45 degC isothermal vs. 25 degC narrow
window). That pattern points to the porosity-isolation mechanism's own
strength being the dominant factor, not the condition. `SI_BETA_LAM_ISO=1.0`
was carried over from `dma/FINDINGS_2026-09-20.md`'s "Finding 3's old
baseline" — a config built around a *different* BoL porosity (0.13) and
`SI_K_SEI_MULT` (6e-4) than the currently-adopted baseline (BoL 0.25,
`K_SEI_MULT` 1.0e-3). It was never recalibrated against the recipe these
two forks actually use, so `1.0` is likely simply too aggressive here, not
a meaningful physical result. **Recommendation:** before drawing any
conclusion about the isolation mechanism's condition-sensitivity, sweep
`SI_BETA_LAM_ISO` down (e.g. 0.01–0.3) against the *baseline* condition
first, find a value that reproduces something like the un-isolated
baseline's own knee timing (~EFC 100–110), and only then compare across
conditions at that calibrated value — otherwise both these reruns just
show "the mechanism is on and very strong," not anything about 45 degC or
the narrow window specifically.

### Open follow-ups

- Calibrate `SI_BETA_LAM_ISO` (and/or `SI_REDIRECT_LAM_YIELD`,
  `SI_LAM_ISO_EXPONENT`, `SI_TAU_LAM_ISO`) against the baseline condition
  before trusting any isolation-enabled cross-condition comparison (see
  above).
- Re-run both (isolation on, once calibrated, and/or isolation off) with a
  larger `C64_MAX_CYCLES`/lower `C64_SOH_FLOOR` to reach an EFC range
  comparable to real RPT4/RPT5 (~170–200) before drawing firm conclusions
  about the LLI/LAM_Si split shift.
- Investigate the narrow-window SoH wiggle at EFC~75–90 seen in the
  isolation-off run (numerical artifact vs. genuine).
- Once collaborators' real data arrives for either condition, wire it back
  into each fork's `plot_all()` (see the "No real-data overlay" note above).
- The solvent-diffusion mechanism below (2026-09-25/26) has only been
  checked against the knee-EFC metric, at exactly two temperatures. It has
  NOT been checked against the fuller multi-panel fit (LAM split, LLI,
  expansion shape, stoichiometry) the baseline's "ec reaction limited"
  recipe was tuned against — do that before treating it as more than a
  working hypothesis for the knee-timing puzzle specifically.
- Narrow-voltage-window condition has not yet been re-tested with this new
  mechanism (only `high_temp_45C/` has).

## Knee-timing investigation: why does the real 45 degC knee arrive at
## EFC~300 instead of earlier? (2026-09-25/26)

**Motivating observation** (real experimental data, not yet in this repo):
at 45 degC the cell's knee arrives at EFC~300 — *later*, not earlier, than
the 25 degC baseline's own anchor knee at EFC~101. This is the opposite of
the "obvious" expectation that higher T accelerates SEI growth (hence
faster porosity closure, hence an *earlier* knee). Below is the mechanism
search that eventually reproduced this, including a mistake made and
caught along the way — kept in for anyone re-deriving this later.

### Mechanisms ruled out (kept Si's SEI as "ec reaction limited")

All of these were tested at 45 degC against the baseline's "ec reaction
limited" SEI option for Si (unchanged growth law, just perturbing rate
parameters), all with **zero or insufficient effect on the knee's onset**:

- **D_ec magnitude scaling** and **D_ec decoupled activation energy**
  (temporarily required core `lithium_ion_parameters.py`/`sei_growth.py`
  edits, since reverted): D_ec only enters once SEI has already thickened
  past onset — reshapes the post-onset tail, never moves the onset itself,
  regardless of parameterisation.
- **Crack-rate Arrhenius** (`SI_CRACK_EAC`, still wired into this fork's
  `_silicon_cracking_rate_scaled`, default kept at the physically-motivated
  but untested 38000 J/mol): zero effect, because the isolation-redirect
  term (`a_j_sei`) that drives the knee excludes the "SEI on cracks" current
  entirely — confirmed via direct read of `loss_active_material.py`.
- **eta_SEI's natural partial compensation**: real (measured net per-cycle
  ratio 1.16x vs. a naive 2.62x Arrhenius boost with no compensation at
  all), but nowhere near enough — net effect is still acceleration, not the
  ~3x delay needed.

**Conclusion**: with Si's SEI kept as "ec reaction limited", *any* positive
`E_sei` (needed since "SEI grows faster at higher temperature" is genuine
physics, confirmed directly from `sei_growth.py`) can only ever accelerate
the knee at higher T, not delay it — the `eta_SEI` overpotential term is
the only compensating channel available, and it's insufficient by itself.
A different growth law was needed, not just different rate constants.

### The (initially mis-evaluated) pivot: "solvent-diffusion limited" for Si

Switched Si's SEI option to `"solvent-diffusion limited"` (graphite
unchanged): `j_sei = -D_sol*c_sol*F/L_sei`, no `eta_SEI` term at all. Scaled
the (unfitted, literature-default) `D_sol` up via a flat multiplier until
the 45 degC knee landed at EFC~300 (`SOLVENT_MULT=5.2` did this: knee
EFC=300.2–300.9 by both metrics).

**This was reported as a success — incorrectly.** The mistake: not checking
whether the *same* configuration's 25 degC prediction was self-consistent
with the real 25 degC anchor before declaring the 45 degC number a win.
Running the identical config at 25 degC gave knee EFC≈504–539 — over 5x
later than the real anchor of ~101. Going from this (uncalibrated) 25 degC
value to 45 degC, the knee moved *earlier* (504→300) — ordinary Arrhenius
acceleration, the exact opposite of the delay the real data shows. The
EFC=300 "hit" at 45 degC was a numerical coincidence of an arbitrarily-sized
rate constant, not evidence the mechanism reproduces the real physics.

### The corrected, jointly-calibrated mechanism

Key structural fact: the universal Arrhenius factor applied in
`sei_growth.py`, `exp(E_sei/R*(1/T_ref - 1/T))`, is **identically 1 at
T=T_ref=298.15K regardless of E_sei's value or sign**. This decouples the
two calibration targets cleanly:

- **`D_sol` magnitude** (`SOLVENT_MULT`) alone sets the **25 degC absolute**
  knee EFC (E_sei is irrelevant there by construction).
- **Si's own `"Secondary: SEI growth activation energy [J.mol-1]"`**
  (independently settable from graphite's, which stays at the default
  +38000 J/mol) alone sets the **45-vs-25 degC ratio**. A ratio >1 (delayed
  knee at higher T, as the real data needs) requires this to be
  **negative** — an effectively negative apparent activation energy for
  this specific pathway.

Two-step calibration (using the established `EFC_knee ~ mult^-0.85` scaling
from the earlier 45 degC-only sweep, applied as a starting estimate then
corrected against actual runs):

| Step | `C64_SOLVENT_MULT` | `C64_SI_ESEI` [J/mol] | T | Result (RPT-based / continuous knee EFC) | Target |
|---|---|---|---|---|---|
| 1 | 34 | 38000 (default, irrelevant at T_ref) | 25 degC | 88.9 / **101.8** | ~101 |
| 2 | 34 | -25000 (first estimate) | 45 degC | 218.0 / 218.7 | ~300 (undershoot) |
| 3 | 34 | **-37500** (corrected) | 45 degC | 301.3 / **302.1** | ~300 |

**Converged config**: `SOLVENT_MULT=34` (i.e. `D_sol = 8.5e-21 m2.s-1` at
T_ref), `Secondary: SEI growth activation energy = -37500 J/mol`,
`SI_CRACK_MULT=1.0` (genuine, unsuppressed Si cracking restored — gives
LAM_Si a gradual, cracking-driven contribution from BOL rather than a sharp
isolation-only step, per instruction). This single combination, run at both
temperatures with nothing else changed, gives:

- 25 degC: knee EFC ≈ 101.8 (continuous) / 88.9 (RPT) — matches the real
  anchor (~101).
- 45 degC: knee EFC ≈ 302.1 (continuous) / 301.3 (RPT) — matches the real
  high-T target (~300).

This is now applied as the default recipe in `high_temp_45C/
cell064_degradation_fit_45C.py` (SEI option changed to
`(("ec reaction limited", "solvent-diffusion limited"), "none")`,
`SI_CRACK_MULT` default 0.1→1.0, two new env-var-overridable constants
`SOLVENT_MULT`/`SI_ESEI` added with the converged values as defaults — see
that file's own comments for the full derivation, duplicated there for
anyone reading the fork in isolation). All exploratory tuning-run plots
(iter01 through iter03, plus the initial mis-evaluated 5.2-multiplier runs)
are saved under
`../degradation_test_matrix/sweeps/porosity_isolation_calib/` with
`_iterNN_...` suffixes — kept per instruction, not cleaned up.

**Physical caveat, not to be understated**: a genuinely negative apparent
activation energy for SEI growth is unusual. It has *some* grounding in
passivation-film literature (a competing, higher-T-favoured pathway can
form a more compact/protective interphase, lowering the net observed rate
despite faster elementary kinetics) — but here it is a back-calculated,
phenomenological number fit to reproduce two knee-EFC values, not derived
from a literature source or real high-T CELL064 data (none exists yet).
It has only been validated against the single knee-EFC metric at exactly
two temperatures — not against the fuller diagnostic fit (LAM split, LLI
split, expansion shape, stoichiometric window) the baseline recipe was
actually tuned against. Treat as a working hypothesis that closes the
knee-timing puzzle quantitatively, pending real high-T data and a broader
diagnostic check, not as a validated replacement mechanism.

### Second correction: matching the EFC number is not matching the knee

The "converged config" above was *also* reported prematurely. The knee-EFC
*number* matched (~300), but the underlying SoH-vs-EFC curve at 45 degC had
**no visible knee at all** — a smooth, continuously-curving decline from
EFC 0 to 400+, no bend anywhere near 300. This was caught by actually
looking at the plot, not by trusting the printed metric — the same failure
mode as the earlier single-temperature-fit mistake, just less obvious
(this time the *number* was right, only the *shape* was wrong). Checked
retroactively: the original mis-evaluated `SOLVENT_MULT=5.2` run (the very
first "success") never had a genuine knee shape either, at any point in
this investigation — every solvent-diffusion config tried so far, at any
D_sol/E_sei decomposition, gave a smooth curve at 45 degC.

**Why**: at a fixed operating temperature, only the *absolute* effective
`D_sol` (i.e. `SOLVENT_MULT * arrhenius(T, E_sei)`, not how it's split
between the two factors) enters the growth law. The 25 degC config's
genuine sharp knee runs at `D_sol_eff = 8.5e-21` (T_ref, arrhenius=1). Every
45 degC config tried (mult=5.2 default-E, mult=34/E=-25000, mult=34/
E=-37500, mult=34/both-electrodes-E=-37500) works out to `D_sol_eff` in the
`3.3–6.6e-21` range — all well below that — and all of them, without
exception, gave a smooth non-knee curve. Trying to hit both a *later*
timing and *this same absolute magnitude* by scaling `D_sol`/`E_sei` alone
is self-contradicting: a large-enough `D_sol_eff` to trigger a genuine
sharp isolation cascade inherently reaches the porosity floor early, and
slowing it down (via mult or E_sei, same effect either way) to delay the
timing also drops it below the threshold needed for the cascade to runaway
sharply within a normal cycle budget — it just fizzles into a smooth,
gradually-decelerating-then-reaccelerating curve instead.

**Fix**: sharpen the isolation gate itself
(`"Secondary: Negative electrode porosity-isolation LAM exponent"`,
`SI_LAM_ISO_EXPONENT`) — a structural nonlinearity-shape parameter,
independent of the SEI growth rate that sets the *timing*. Swept 3.0 (this
fork's prior default) → 8.0 → 15.0 at the 45 degC config
(`SOLVENT_MULT=34`, `Secondary SEI growth activation energy=-37500`,
graphite's left unmodified at +38000 — reverting the both-electrodes test
above, which didn't help since it didn't address the real cause): 8.0 gave
a visibly sharper but still gradual bend (knee EFC~309–310); **15.0 gives a
genuine, clearly visible knee** (flat-ish decline from EFC 0–~280, then a
distinctly steeper decline from ~280 onward — same qualitative shape as the
25 degC curve), with RPT-based knee EFC~311.1 / SoH-91% crossing EFC~284.3,
both close to the ~300 target.

**Checked for self-consistency** (this parameter is not itself
temperature-dependent, so it must be shared, not re-tuned per condition):
re-ran the 25 degC config with the same `SI_LAM_ISO_EXPONENT=15.0` —
its own genuine sharp knee is preserved, right at the real anchor (SoH-91%
crossing EFC=97.6, RPT-based knee EFC=89.6, vs. anchor ~101), if anything
slightly *improved* over the exponent=3.0 baseline (was 84.6).

**Final converged config** (supersedes the "converged config" above):
`SOLVENT_MULT=34`, `Secondary: SEI growth activation energy=-37500 J/mol`,
`SI_CRACK_MULT=1.0`, **`SI_LAM_ISO_EXPONENT=15.0`** (was 3.0). Applied as
the new defaults in `high_temp_45C/cell064_degradation_fit_45C.py`. Both
temperatures now show a genuine, visually-verified knee at the right
location — not just a matching number. The physical caveat above about the
negative apparent activation energy still applies in full; this correction
doesn't change that, only the shape-verification failure it was compounding.

**Process lesson, worth restating**: a printed "knee EFC" number is not
suffient evidence of a genuine knee. Always look at the actual SoH-vs-EFC
plot (and ideally the LAM_Si panel, which shows the same transition more
starkly) before reporting a tuning target as met — the steepest-single-
interval metric can and does land on smooth, non-knee-shaped curves.
