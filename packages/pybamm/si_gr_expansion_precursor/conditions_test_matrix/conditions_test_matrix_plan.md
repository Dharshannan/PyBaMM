# Conditions Test Matrix: Plan

Status: **planning only** — no implementation yet, per request. Review first, then kick off.

## 0. Purpose and scope

Test whether the **same** degradation recipe (the tuned `si_gr_expansion` baseline used
throughout `test_pore_buffering/` and `expansion_pattern_test/` — `SI_MULT=18`,
`LAM_MULT=1.0`, unchanged) reproduces the *qualitative* pre-hump/hump/post-knee pattern
across the 5 real operating conditions in `../fig2_population_lead_time_v2.png`
(documented in `../expansion_precursor_test_plan.md` Part A.0g), by varying only
**operating-condition inputs** (temperature, C-rate, voltage window) and the
**pore-buffering shape** (`f0`, `eps_transfer_width` — per A.0l/A.0m, standing in for
pressure) — not the degradation-kinetics recipe itself.

**Explicit scope decision, per instruction**: low pressure is represented purely as a
pore-buffering *parameter fit* (A.0l's `f0=0.4, width=0.17`), not a physical
`K_stack(P_stack)` coupling in the core model equations. This plan does **not** revisit
that decision — see `expansion_precursor_test_plan.md` Part C.4 for why the literal
pressure-coupling route was tried and set aside (A.0h).

This is a different axis from `../degradation_test_matrix/` (which varies *mechanisms* at
one fixed operating condition): here the mechanism recipe is held fixed and only the
5 conditions vary, mirroring the real data's own experimental design.

## 1. The 5 conditions

| # | Condition | Temperature | Ageing C-rate | Ageing voltage window | Pore buffering `(f0, width)` | Real reference (Lead EFC → Knee EFC) |
|---|---|---|---|---|---|---|
| 1 | Baseline | 298.15 K (25 °C) | C/3 | 2.6–4.2 V | `(0.7, 0.01)` — existing default | 121 → 145 |
| 2 | **Low pressure** (proxy) | 298.15 K | C/3 | 2.6–4.2 V | **`(0.4, 0.17)`** — A.0l's validated point | 119 → 158 |
| 3 | High C-rate | 298.15 K | **2C** | 2.6–4.2 V | `(0.7, 0.01)` | 95 → 143 |
| 4 | High temperature | **318.15 K (45 °C)** | C/3 | 2.6–4.2 V | `(0.7, 0.01)` | 296 → 319 |
| 5 | Narrow SoC window | 298.15 K | C/3 | **~3.0–4.15 V** (approximating real 3.15–4.12 V) | `(0.7, 0.01)` | 255 → 298 |

Everything else — `SI_MULT`, `LAM_MULT`, `SI_CRIT_STRESS`, `SI_LAM_PROP_BASELINE`,
`GR_DIV`, `GR_CRIT_STRESS`, `GR_LAM_PROP_BASELINE`, crack-rate multipliers, porosity
floor, SEI exponent caps, composition (`GR_EPS`/`SI_EPS`/`SI_MAX_CONC_NEEDED`) — stays
**exactly** the tuned baseline recipe from `test_pore_buffering/
pore_buffering_degradation_test_1000cyc.py` / `expansion_pattern_test/
physical_f0_width_grid_test.py`, for all 5 conditions. Formation stays fixed (full-range,
gentle-rate) for all 5 — only the *ageing* protocol/parameters vary per condition, matching
how every earlier sweep script in this project has treated formation as a fixed BOL
reference point.

**Known minor mismatch worth a decision**: existing scripts' ageing cycle uses "Discharge
at C/3 until **2.5 V**" (not 2.6 V, the real baseline's actual cutoff). Recommend aligning
to 2.6 V for this specific test matrix, since reproducing the real conditions closely is
the whole point here — flagged as a one-line change, not a re-tune, since the tuned
recipe's knee timing was calibrated against the 2.5 V convention and shifting to 2.6 V by
0.1 V is unlikely to move it meaningfully, but worth a quick single-condition check (Phase
1) before assuming that.

## 2. Phase 0 — Temperature-dependence audit (do this first; changes what Condition 4 means)

Requested explicitly: check that temperature actually propagates correctly through the
SEI/cracking models and the Arrhenius terms on particle/electrolyte diffusion, **before**
running the high-temperature condition and interpreting its result. Read-only investigation
this session already turned up concrete findings, not just "worth checking":

| Mechanism | Status | Evidence |
|---|---|---|
| Graphite particle diffusivity | **Genuine Arrhenius** | `graphite_LGM50_diffusivity_Chen2020` (`si_gr_expansion.py:101`), `E_D_s=3.03e4 J/mol` |
| Silicon particle diffusivity | **Genuine Arrhenius** | `silicon_LGM50_diffusivity_Bonkile2024` (`:342`), `E_D_s=4.82e4 J/mol` |
| Electrolyte diffusivity/conductivity | **Genuine Arrhenius** | `electrolyte_diffusivity_Nyman2008_arrhenius` / `electrolyte_conductivity_Nyman2008_arrhenius` (`:1019-1022`), not overridden by any script's `build_param_updates` |
| Graphite exchange-current density | Indirect only | `graphite_LGM50_electrolyte_exchange_current_density_Chen2020` has its own Arrhenius (`E_r=35000`, `:127`), but this is the *intercalation* reaction, not SEI |
| **Graphite/silicon particle cracking rate** | **Not implemented — confirmed, not assumed** | Both `graphite_cracking_rate_Ai2020` (`:214`) and `silicon_cracking_rate_Ai2020` (`:476`) contain `Eac_cr = 0  # to be implemented` in the source — the Arrhenius scaffolding exists but is hard-set to zero activation energy, i.e. **no temperature dependence at all** by design, even before any script touches it. On top of that, every sweep script in this project (`build_param_updates`) overrides the cracking-rate parameter with a **plain scalar** (`BASE_CRACK_RATE * SI_CRACK_RATE_MULT`), which would strip out any `Eac_cr` even if it *were* nonzero. |
| **SEI reaction exchange current density (both phases)** | **Weak/indirect only** | `"Primary/Secondary: SEI reaction exchange current density [A.m-2]"` (`:813`, `:830`) are plain scalar constants in the *default* parameter set, not Arrhenius `FunctionParameter`s. The `"reaction limited"` SEI submodel's own rate expression still has *some* T-sensitivity through its Butler-Volmer-style overpotential term, but the pre-exponential rate constant itself doesn't scale with T — real SEI growth is well known to be strongly Arrhenius-activated, so this recipe likely **under-predicts** temperature acceleration of LLI specifically. |
| LAM proportional/exponential rate constants | Indirect only | Plain constants (`SI_LAM_PROP_BASELINE`, `GR_LAM_PROP_BASELINE`), no explicit T term; picks up *some* T-sensitivity indirectly via T-dependent diffusivity feeding into concentration-gradient-driven stress. |

**Decision needed before Condition 4 can be run meaningfully** (surface for review, not
resolved here):
- **(a) Run as-is, flag the limitation.** Cheapest: Condition 4 will show whatever
  acceleration comes from particle/electrolyte diffusivity's genuine Arrhenius terms alone
  (which do feed into concentration-polarisation-driven stress and hence LAM/cracking
  indirectly, and into cell resistance/knee timing) — but cracking rate and SEI-driven LLI
  will show **no direct** rate acceleration. Given A.0's own high-temperature hypothesis
  (A.2) already speculated the *opposite*-than-naive timing (later knee at higher T) might
  come from reduced diffusion-driven stress outweighing faster SEI kinetics — with SEI
  kinetics not actually T-accelerated in this recipe at all, that hypothesis's SEI-kinetics
  side is moot by construction, which changes what a Condition-4 result would even be
  testing. Worth stating explicitly in whatever writeup follows this test.
- **(b) Restore cracking-rate Arrhenius behaviour.** Cheap, no new physics: stop
  overriding the crack-rate parameter with a flat constant; instead scale the *function*
  (e.g. `lambda T_dim: silicon_cracking_rate_Ai2020(T_dim) * SI_CRACK_RATE_MULT`) so the
  already-present (if currently zeroed) Arrhenius scaffolding is at least wired through —
  though `Eac_cr=0` in the source means this alone still does nothing unless a nonzero
  activation energy is also chosen (see next point). Setting a literature-plausible
  `Eac_cr` (e.g. a value in the range typically reported for stress-corrosion/fatigue
  crack growth acceleration, on the order of 1–5×10⁴ J/mol, matching the *diffusivity*
  activation energies already in this same parameter set as a starting-point analogy) is a
  genuinely new modelling decision, not a "restore," and should be treated as its own
  small calibration exercise, not a rubber-stamp default.
- **(c) Add genuine Arrhenius SEI kinetics.** Moderate lift: wrap
  `"Primary/Secondary: SEI reaction exchange current density [A.m-2]"` as a
  `FunctionParameter(T)` with a chosen activation energy, mirroring the pattern already
  used for particle/electrolyte diffusivity in this same file. Also a new modelling
  decision (which `Ea` to pick), not a bug fix.

**Recommendation for the first pass**: (a), explicitly documented as a known limitation,
to get a genuine first read on what the *existing* Arrhenius pathways (diffusivity →
stress → LAM/cracking, and diffusivity/conductivity → resistance → knee timing) alone
produce, before deciding whether (b)/(c) are worth the extra calibration effort. This
matches the project's established "test the cheap thing first, escalate only if it doesn't
close the gap" pattern (e.g. A.0 → A.0a → A.0i → A.0j → A.0l's escalation path).

## 3. Phase 1 — Single-condition smoke tests

Run each of the 5 conditions individually (not yet the full comparison), to 50% SoH
(matching the now-established convention from `expansion_pattern_test/
physical_f0_width_grid_test.py`), confirming each:
- Solves without solver failure across its full formation + ageing-to-50%-SoH run.
- Produces a plausible SoH trajectory (monotonically declining, no numerical artefacts).
- Produces a plausible reversible-amplitude trajectory (right order of magnitude, no
  NaNs/discontinuities).

Specifically worth checking during this phase, one at a time:
- **Condition 3 (high C-rate)**: does the 2C ageing cycle need adjusted current/voltage
  step definitions (e.g. `"Discharge at 2C until 2.6 V"` vs. the existing `"Discharge at
  C/3 until ..."` string), and does the IDAKLU solver's tolerances (currently
  `root_tol=1e-06, atol=1e-06, rtol=1e-06`, unchanged across all prior scripts) still
  converge cleanly at the faster rate — 2C is a substantially stiffer transient than C/3.
- **Condition 4 (high temperature)**: confirm setting `"Ambient temperature [K]"` and
  `"Initial temperature [K]"` to 318.15 (rather than enabling a full non-isothermal
  thermal submodel) is sufficient — every script in this project runs isothermal by
  default (no `"thermal"` option set), and Arrhenius terms only need a single `T` value,
  so a fixed elevated ambient/initial temperature should be enough without opening up
  thermal-submodel scope creep. Confirm this assumption explicitly rather than silently
  relying on it.
- **Condition 5 (narrow SoC window)**: confirm the ageing cycle's cutoffs (approximating
  3.15–4.12 V) don't collide with the existing formation protocol's own targets (formation
  ends at a full-range 4.2 V hold — ageing then needs to *start* its first discharge from
  that full-range state down to the narrow window's lower bound, which is fine, but the
  narrow window's upper cutoff, ~4.15 V, sits *below* formation's 4.2 V, so the very first
  ageing charge step needs to target the narrower 4.15 V, not 4.2 V, from cycle 1 onward —
  a detail easy to get wrong by copy-pasting the baseline ageing tuple).
- **The 2.5 V → 2.6 V voltage-cutoff realignment** (§1's "known minor mismatch"): run
  baseline at both cutoffs briefly and confirm knee timing doesn't shift meaningfully
  before committing to 2.6 V for the whole matrix.

## 4. Phase 2 — Knee/Lead-EFC detection methodology

The real data's own figure overlays visible piecewise-linear fits (a flat/gently-sloped
segment, then a steeper post-knee segment) — the detection method should mirror this
rather than inventing an unrelated one, so model and real "Lead EFC"/"Knee EFC" numbers are
comparable on the same basis.

**Proposed method**: two-segment (or three-segment, for the "Lead" point too) piecewise
linear regression on capacity-vs-EFC, via `scipy.optimize` (e.g. a simple grid search over
candidate breakpoints minimising total squared residual of two independent line fits, or
`scipy.optimize.curve_fit` on a continuous piecewise-linear model with the breakpoint as a
free parameter) — a standard, reproducible, already-visually-validated-against-the-source-
figure approach, not a novel heuristic. Concretely:
- **Knee EFC**: the breakpoint of a 2-segment piecewise-linear fit to discharge capacity
  vs. EFC over the condition's full run.
- **Lead EFC**: per the source figure's own definition (visible as a *second*, earlier
  dashed line in panel a, and explicitly as the x-axis in panel b, `EFC − EFC_knee`) —
  the point where the *reversible expansion* trace first detectably departs from its own
  early-life trend, ahead of the knee. Propose the same piecewise-linear-breakpoint method
  applied to the expansion-vs-EFC trace instead of capacity, as the most direct analogue.

**EFC conversion**: all of this project's scripts so far report throughput on the
x-axis as raw `"Throughput capacity [A.h]"`, not EFC. Standard conversion:
`EFC = Throughput capacity [A.h] / (2 * initial measured capacity [A.h])` (one full cycle
= one full discharge + one full charge = `2×` capacity of throughput). Needed to make any
of this comparable to the real data's EFC-based x-axis — not yet applied in any script to
date; this test matrix would be the first place it's needed.

**Validation**: run the breakpoint method on baseline first and sanity-check the resulting
knee EFC against the already-known throughput-based knee location from
`physical_f0_width_grid_test.py`'s baseline run (~850–900 A·h, ~140 cycles) before trusting
it on the other 4 conditions.

## 5. Phase 3 — Full comparison script and outputs

One script, looping over the 5 conditions (parameterised the same way
`physical_f0_width_grid_test.py` parameterises its 6-point grid), each run to 50% SoH,
producing:
- A 5-condition comparison figure mirroring `fig2_population_lead_time_v2.png`'s panel-a
  layout as closely as practical: discharge capacity and reversible expansion vs. EFC
  (using the Phase 2 conversion), 5 side-by-side (or overlaid) panels, with detected
  Lead/Knee EFC annotated the same way.
- An EFC-aligned overlay panel (mirroring panel b): all 5 conditions' reversible expansion
  vs. `EFC − EFC_knee`, on one axis.
- A Lead-EFC-by-category-style summary (mirroring panel c, though with only 5 single-run
  points rather than the real data's repeated/multi-cell statistics — report as single
  points, not means ± spread, and say so explicitly on the plot to avoid implying
  statistics that aren't there).

## 6. Phase 4 — Cross-check against real numbers (A.0m-style)

For each condition, report the same comparison A.0m already did for low pressure alone:
BoL amplitude, peak amplitude, peak/BoL ratio, Lead EFC, Knee EFC, model vs. real,
side by side in a table. Explicitly expect (per A.0g/A.0m's own framing): a qualitative
shape match is the goal, not a quantitative fit — the absolute magnitude scale is expected
to differ (as it already does for low pressure alone, ~6×, since this diagnostic recipe
isn't calibrated to the real paper's specific cell). Report the *ordering* across
conditions (e.g. does the model's Lead-EFC ranking across the 5 conditions match the real
data's — C-rate longest, temperature shortest, per panel c) as the more meaningful check
than any single condition's absolute numbers.

## 7. Open decisions to confirm before implementation

1. Phase 0's (a)/(b)/(c) choice for temperature (recommend (a) first).
2. The 2.5 V → 2.6 V cutoff realignment (recommend yes, after a quick Phase-1 check).
3. Whether "narrow SoC window" should use the paper's exact 3.15–4.12 V or the
   approximated ~3.0–4.15 V suggested above (pick whichever keeps the ageing cycle
   syntactically simple against the existing formation protocol — exact values TBD during
   Phase 1's smoke test).
4. Whether Phase 3's script should reuse `physical_f0_width_grid_test.py` as a direct
   template (recommended — same recipe, same SoH-termination helper, same plotting
   conventions) or be written fresh.
