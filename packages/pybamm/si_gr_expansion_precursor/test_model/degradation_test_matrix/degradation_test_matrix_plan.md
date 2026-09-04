# Degradation Mechanism Test Matrix: Plan

Status: **the full Phase-2 matrix is planning only** — no implementation yet, per request,
matching `../conditions_test_matrix/`'s "plan first" convention. **§1a's wrapper extension
(`sei_ec_coupling="ec_reaction"`) has since been implemented and smoke-tested**, per an
explicit follow-up instruction to do that specific piece first — see
`../ec_dryout/implementation_plan.md` §7 for the full change record. **§1a's base-recipe
retune is also now done**: `k_sei`×0.0017 (short recipe) and its `v2` 1000-cycle
`TIMESCALE_STRETCH` stretch (`k_sei`÷7, LAM proportional terms÷4.7) both closely match
`"reaction limited"`'s knee/full-lifetime/LLI/LAM targets — see the "Joint calibration"
entry below. Everything else (1b onward, and all of Phase 2) remains plan-only.

## 0. Recap: what this is testing, and how it relates to the main plan

This elaborates Part B of `../expansion_precursor_test_plan.md` (page 6's 5-row mechanism
matrix: Si cracking / SEI porosity loss / dry-out, each on or off, pore buffering as the
constant "lens"), at a **single**, fixed operating condition (the tuned baseline recipe —
this is the axis orthogonal to `../conditions_test_matrix/`, which fixes mechanisms and
varies operating conditions instead). B.1/B.2/B.5 already establish the row definitions and
row-by-row PyBaMM option mapping; this document adds the concrete parameter values, the
dry-out/composite integration work needed to actually run rows 4–5, and the full script
plan.

## 1. Phase 1 — Couple the dry-out wrapper to the composite base model (do this first)

**Why this is genuinely new work, not just "turn on an existing flag."**
`ec_dryout_wrapper.py` already supports `composite=True` as a first-class, *default*
parameter (`run_ec_dryout_degradation(..., composite=True)`, `ECDryoutLedger.__init__`,
`sei_capacity_loss_Ah`, `extract_batch_series` all branch on it correctly — confirmed by
reading the source, not assumed) — but it has never actually been *run* against
`si_gr_expansion`. Every existing dry-out validation (`ruihe_dryout_validation.py`) uses
`OKane2022` (`composite=False`) instead, per this session's earlier explicit pivot away
from `si_gr_expansion` for that specific validation. So the composite code path is written
and self-consistent (e.g. it already correctly updates *both* `"Primary:"` and
`"Secondary:"` bulk-solvent-concentration parameters from a single shared EC-pool ledger
value each batch — verified in `run_ec_dryout_degradation`, not just assumed from the
docstring) — what's missing is exercising it against the actual tuned recipe and pore
buffering together for the first time.

**1a. The SEI submodel mismatch (B.3) — resolved by extending the wrapper, not the recipe.**
Superseded update: rather than recalibrating the tuned recipe under `"solvent-diffusion
limited"` (option (a) originally proposed here), `ec_dryout_wrapper.py` has been extended
instead — **done, not just planned** — to support `"ec reaction limited"` (Yang et al.
2017), a genuine reaction-kinetic SEI law (Butler-Volmer-style, can accelerate into a knee
— unlike `"solvent-diffusion limited"`'s purely self-limiting diffusion law, which cannot)
that is *also* genuinely coupled to an EC-concentration parameter (`c_ec_0`) the wrapper
can drive dry-out through. This sits functionally between `"reaction limited"` (no
concentration term at all — nothing to attach dry-out to) and `"solvent-diffusion
limited"` (concentration-only, no reaction kinetics) — exactly the missing middle ground
option (b) above was looking for, except it turns out to already exist upstream in
vanilla PyBaMM, needing **zero `src/pybamm` changes**, only wrapper-side parameter-key
selection. Full rationale, exact parameter-key mapping, and a passing end-to-end smoke
test (composite `si_gr_expansion`, pore buffering on, `sei_ec_coupling="ec_reaction"`, 10
cycles, `c_EC` genuinely evolving 4541 → 3300.0 → 2643.2 mol/m³) are documented in
`../ec_dryout/implementation_plan.md` §7. New API: `run_ec_dryout_degradation(...,
sei_ec_coupling="ec_reaction")` (default remains `"solvent_diffusion"`, exactly preserving
every previously-validated behaviour — confirmed via regression smoke test, output
unchanged).

**Still needed before Phase 1c (unchanged from the original plan, just re-scoped to
`"ec reaction limited"` instead of `"solvent-diffusion limited"`)**:
1. Switch `"SEI": "reaction limited"` → `"ec reaction limited"` in the base-case recipe
   (`test_pore_buffering/pore_buffering_degradation_test_1000cyc.py`'s options), **keep
   every other tuned parameter as-is** (`SI_MULT`, `SI_CRIT_STRESS`, `SI_LAM_PROP_BASELINE`,
   `GR_DIV`, `GR_CRIT_STRESS`, `GR_LAM_PROP_BASELINE`, crack-rate multipliers), and run it
   once in isolation (no dry-out yet) to see how far off the knee timing/LLI/LAM balance
   lands from the already-validated `"reaction limited"` baseline. `"ec reaction limited"`'s
   own kinetic parameters (`k_sei`, `D_ec`, `c_ec_0` — already present in
   `si_gr_expansion.py:824-826`/`:841-843`, likely unused/untuned defaults) will very
   plausibly need adjusting to land anywhere near the existing knee timing, since they've
   never been calibrated against this recipe's targets — expect this step to reveal how far
   off the defaults are, not to just work immediately.

   **Done — result: far off, as expected, cleanly and without solver issues.**
   `ec_reaction_limited_swap_test.py` ran `"reaction limited"` (already-validated baseline)
   side by side with `"ec reaction limited"` at its untouched defaults, same recipe
   otherwise, both to 50% SoH:

   | SEI option | Cycles to 50% SoH | Throughput to 50% SoH | Final LLI | Final LAM (neg) |
   |---|---|---|---|---|
   | `reaction limited` (baseline) | 160 | ~1,100 A·h | 33.0% | 23.4% |
   | `ec reaction limited` (untuned defaults) | **20** | **~100 A·h** | 40.1% | **2.4%** |

   `"ec reaction limited"`'s defaults are roughly **10–15× too fast** — the cell collapses
   almost entirely from LLI alone (LAM barely engages, 2.4% vs. 23.4%) within a single
   20-cycle batch, and `k` saturates to 1.0 (buffering transition fully engaged) almost
   immediately rather than tracking a genuine multi-hundred-cycle porosity decline. The
   expansion trace shows a compressed rise-then-fall over just ~5 data points — technically
   present, but far too fast to call a resolved "hump" in any meaningful sense. This ran
   cleanly with no solver failures (confirms the "far off" result is a genuine kinetics
   mismatch, not the exponent-cap/solver-robustness issue flagged separately below).

   **Step 2/3 retune, done — `k_sei` alone saturates; `D_ec` has real leverage but reveals
   a genuine timing/balance coupling.** `ec_reaction_limited_ksei_retune_test.py` swept
   `k_sei` down 2 orders of magnitude (×1 → ×0.01, i.e. `1e-12` → `1e-14` m/s): knee
   throughput only moved 29 → 122 A·h, clearly sub-linear/saturating on a log-log plot —
   confirms `j_sei = -F·c_0·k_exp/(1+L_over_D·k_exp)` is deep in the
   *diffusion*-dominated regime here (`L_over_D·k_exp ≫ 1`), where `k_sei` has
   diminishing leverage by construction. `ec_reaction_limited_dec_retune_test.py` then
   swept `D_ec` down 4 orders of magnitude at `k_sei` fixed at ×0.01:

   | `D_ec` mult | Knee throughput | Final LLI | Final LAM(neg) |
   |---|---|---|---|
   | ×1 | ~122 A·h | 35.4% | 7.3% |
   | ×0.1 | ~218 | 30.5% | 19.1% |
   | ×0.01 | ~574 | 23.4% | 56.6% |
   | ×0.001 | ~866 | 18.1% | 58.6% |
   | ×0.0001 | **~999** | 16.2% | 59.0% |
   | **Target (`reaction limited`)** | **~1,100** | **33.0%** | **23.4%** |

   `D_ec` has genuine, non-saturating leverage — knee timing reaches 91% of the target
   (999 vs. 1,100 A·h) and is still climbing at the sweep's lower end, unlike `k_sei`'s
   clear plateau. The SoH curves for `×0.01`–`×0.0001` also visibly develop a proper
   flat-plateau-then-knee shape (not the `×1`/`×0.1` cases' near-monotonic collapse),
   qualitatively much closer to the target.

   **But timing and LLI/LAM balance are coupled in a way that punishes a single-lever
   fix, and this is worth understanding, not just noting.** As `D_ec` drops, LAM's final
   share rises sharply (7% → 59%) while LLI's falls (35% → 16%) — the balance has now
   *inverted* relative to target (33% LLI / 23% LAM), overshooting past a good match
   somewhere in this sweep rather than landing on one. Mechanistically: `SI_LAM_PROP_
   BASELINE`/`GR_LAM_PROP_BASELINE` (LAM's own rate constants) are unchanged throughout —
   but slowing SEI/LLI down means the cell survives far more cycles before reaching a
   given SoH, and LAM accumulates roughly per-cycle (stress-driven), so it keeps
   accruing over all those extra cycles regardless of how slow LLI itself is. The tuned
   `"reaction limited"` recipe's LLI/LAM balance is therefore not just a property of its
   own rate constants in isolation — it's tied to *how many cycles* that recipe takes to
   reach its knee. Matching both axes at once likely needs a **joint** retune (LAM-side
   rate constants scaled down too, not just SEI/EC-side), not a sequential one-lever-at-
   a-time search. Not yet attempted — a genuine decision point on how much further
   calibration effort to invest before treating this as "good enough" for a first
   Phase-2 matrix pass (with the balance mismatch documented as a known limitation) versus
   continuing to a fuller joint calibration.

   **Separately, still open**: extending `"ec reaction limited"`'s `k_exp` term with the
   same smooth exponent cap `"reaction limited"` already has (`exponent_max_sei`,
   `sei_growth.py`) was requested but not yet done, and needs care — `exponent_max_sei`
   (`"{pref}SEI reaction exponent cap"`) is currently only defined by `si_gr_expansion.py`
   and `Mayur2024.py`; confirmed via a direct test that `"reaction limited"` itself
   already fails against a standard parameter set lacking it (`Chen2020` + `"reaction
   limited"` raises `KeyError: 'SEI reaction exponent cap' not found`) — i.e. this
   narrow-compatibility trade-off is a pre-existing, already-accepted property of this
   project's fork, not something the cap addition would newly introduce. Also confirmed
   `"ec reaction limited"` has its own existing PyBaMM unit-test coverage
   (`base_lithium_ion_tests.py`) against parameter sets that do *not* define
   `exponent_max_sei` — adding an unconditional reference would break those tests, so this
   still needs the exponent cap gated in a way that doesn't affect callers who never define
   the parameter (e.g. only referenced when `si_gr_expansion`/`Mayur2024`-style parameter
   sets are in use) — not yet resolved, flagged for a decision before implementing.

   **Joint calibration, done — a single-point match found, then stretched to ~1000
   cycles.** Rather than continuing the `D_ec`-only sweep into its inverted-balance
   regime, `ec_reaction_limited_ksei_fine_sweep_test.py` fine-swept `k_sei` alone (`D_ec`
   left at `si_gr_expansion.py`'s default) between the two earlier bracketing points
   (×0.003/×0.001), fixing an earlier knee-vs-full-lifetime metric-comparison error along
   the way (both are now compared like-for-like). **`k_sei` × 0.0017 is an excellent
   single-point match on every axis at once**:

   | | Knee (SoH<90%) | Full lifetime (50% SoH) | Cycles | LLI | LAM(neg) |
   |---|---|---|---|---|---|
   | `reaction limited` target | 457.1 A·h | 1,115.6 A·h | 160 | 33.0% | 23.4% |
   | `ec reaction limited`, `k_sei`×0.0017 (anchor) | 408.3 A·h | 1,061.6 A·h | 160 | 33.2% | 21.7% |

   `ec_reaction_limited_1000cyc_stretch_test.py` then applied the project's established
   `TIMESCALE_STRETCH` pattern (`pore_buffering_degradation_test_1000cyc.py:126-139`) to
   this anchor, dividing every RATE constant (`k_sei`, both LAM proportional terms, both
   cracking rates) by 7 uniformly — the same factor that stretched `reaction limited`'s own
   ~160-cycle recipe to ~1000 cycles. Result: cycle count landed at 840 (84% of the ~1000
   target, reasonably close) and LLI tracked well (32.6% vs. 33.0%), but **LAM_neg
   undershot** (16.1% vs. 21.7% anchor / 23.4% target) rather than being preserved by the
   uniform stretch. Diagnosis: the run landed short of a clean 7× cycle-count
   extrapolation (840 vs. ideal 1,120 = 160×7, a ~0.75× shortfall) — LLI still tracks
   target because SEI growth follows `k_sei` directly, but LAM_neg's shortfall matches that
   same ~0.75× factor almost exactly (21.68% × 0.75 ≈ 16.3% ≈ the observed 16.10%), i.e.
   LAM wasn't behaving differently under stretch, it just had fewer effective cycles to
   accumulate over than a clean 7× extrapolation implies.

   **Fix — decouple the two rate channels — done, and it works.**
   `ec_reaction_limited_1000cyc_stretch_v2_test.py` keeps `k_sei`'s divisor at 7 (unchanged,
   preserves the good LLI/knee-timing match) but gives the Si/Gr LAM proportional terms
   their own, smaller divisor `LAM_STRETCH=4.7` (crack rate stays tied to ÷7, since it
   feeds SEI-on-cracks — an LLI pathway — not the stress-driven LAM term directly):

   | | Knee (SoH<90%) | Full lifetime (50% SoH) | Cycles | LLI | LAM(neg) |
   |---|---|---|---|---|---|
   | `reaction limited` target | 457.1 A·h | 1,115.6 A·h | 160 | 33.0% | 23.4% |
   | v1 uniform stretch (×7 on everything) | 2,477.4 A·h | 5,990.6 A·h | 840 | 32.6% | 16.1% |
   | **v2 LAM-boosted (`k_sei`÷7, LAM÷4.7)** | 2,202.5 A·h | 5,683.9 A·h | **800** | **32.4%** | **23.1%** |

   LAM_neg now matches target almost exactly (23.1% vs. 23.4%), LLI stays close (32.4% vs.
   33.0%), and the SoH/expansion-amplitude curves both show a cleanly resolved knee/hump at
   the ~1000-cycle scale (hump: 10.1 → 13.4 µm). Cycle count (800) dipped slightly from v1
   (840) since faster LAM growth reaches 50% SoH a touch sooner — an expected trade-off, not
   a regression. **This is treated as the working 1000-cycle `"ec reaction limited"`
   recipe for Phase 1b/1c and the Phase 2 matrix**, pending confirmation. `v2` also adds
   RPT (C/10, once per batch) tracking, plotted as orange-triangle scatter over the C/3
   ageing curve — matching `pore_buffering_degradation_test_1000cyc.py`'s own RPT
   convention — on both the SoH and expansion-amplitude panels.

2. If close (same order of magnitude knee EFC, no qualitative shape change in expansion),
   treat the existing `k_sei`/`D_ec`/`c_ec_0` defaults as "good enough" for a first coupled
   dry-out run and defer a full retune.
3. If far off, retune `k_sei`/`c_ec_0` (the LLI-rate-setting side, analogous to how
   `j0_sei` was scaled by `SI_MULT`/`GR_DIV` for `"reaction limited"`) rather than
   re-deriving the whole recipe — `SI_CRIT_STRESS`/LAM constants govern the *mechanical*
   (LAM) side and shouldn't need to move just because the SEI rate law changed.

**1b. Re-derive dry-out kinetics calibration for the composite recipe.**
`ruihe_dryout_validation.py`'s calibration (`SEI_DIFFUSIVITY_MULT=100.0`,
`SEI_MOLAR_VOLUME_DIV=5.0`, scaling OKane2022's `"SEI solvent diffusivity"`/`"SEI partial
molar volume"` under `"solvent-diffusion limited"`) doesn't directly transfer to `"ec
reaction limited"`'s different parameter set (`D_ec`/`c_ec_0` govern EC transport there,
not `"SEI solvent diffusivity"`/`"SEI partial molar volume"`). Once 1a lands on a working
`"ec reaction limited"` composite parameterisation, run the analogous calibration check:
confirm dry-out (at `r_eres=0%`, no reservoir) actually engages within a reasonable cycle
count, watching for the same self-limiting-decay trap `ruihe_dryout_validation.py` found
(naive rate scaling alone wasn't enough there; both `dn_EC` and porosity decrease decayed
together as `c_EC` dropped) — expect `D_ec` and/or the EC "molar volume"-equivalent lever
to need similar joint tuning, not just one constant.

**1c. Validation sweep — the actual "does it work as intended" check requested.**
Once 1a/1b land, reproduce `ruihe_dryout_validation.py`'s own validation structure
(`r_eres ∈ {0%, 6%, 9%}`, `UPDATE_EVERY_N_CYCLES=20`), but on the composite + pore-buffering
base model, and — the genuinely new angle this phase adds — track **reversible expansion
amplitude and the hump**, not just capacity/SoH (the only thing the earlier OKane2022
validation looked at). Two things to check, both currently unknown:
- Does dry-out severity (via `r_eres`) shift the **knee timing** the same clean,
  monotonic way it did for OKane2022 (`ruihe_dryout_validation_fig3.png`'s validated
  result), now under composite + pore buffering?
- Does dry-out severity change the **expansion hump's** shape or timing at all, or does it
  only affect capacity/LLI-side behaviour while expansion stays governed purely by the
  pore-buffering partition (a real, falsifiable, currently-unknown question — dry-out
  changes electrolyte/EC availability, which affects SEI growth rate and hence `ε_struct`'s
  decline rate, which *should* shift when the buffering transition engages, but by how
  much relative to the capacity-side shift is untested)?

## 2. Phase 2 — Full mechanism test matrix

### 2.1 Row definitions (elaborating B.2's table with concrete values)

| Row | Si cracking | SEI porosity loss | Dry-out | Concrete change from tuned baseline |
|---|---|---|---|---|
| Healthy baseline | Off | **On, at normal/untuned rate** | Off | See §2.2 — redefined from B.2's original "Off" |
| Silicon fatigue only | **On** (Si crack rate at tuned value, Gr crack rate → 0) | **On, at normal/untuned rate** | Off | `"Secondary: Negative electrode cracking rate"` at tuned value, `"Primary: Negative electrode cracking rate"` → 0 |
| Pore clogging, Si intact | Off | **On, at tuned (accelerated) rate** | Off | Both phases' cracking rates → 0; SEI kinetics at the full tuned `SI_MULT`/`GR_DIV` |
| Dry-out, Si intact | Off | **On, at tuned rate** | **On** | As above, plus Phase 1's dry-out wrapper loop engaged |
| Fully coupled | **On** | **On, at tuned rate** | **On** | Everything at tuned values — the full recipe as currently validated |

### 2.2 Resolving the "pore buffering true for all" request against the `OptionError` guard

**Flagged conflict, needs your confirmation before implementation.** `"pore buffering":
"true"` requires `"SEI porosity change": "true"` (`base_battery_model.py:766-775`,
enforced as a hard `OptionError`) — it **cannot** be set for a row where SEI porosity
change is off, full stop, regardless of intent. B.2's original table had "Healthy baseline"
and "Silicon fatigue only" with SEI porosity loss **Off** specifically to isolate
"no porosity-driven degradation at all" — under that original definition, pore buffering
literally cannot be turned on for those two rows.

**Proposed resolution** (a redefinition of those two rows, not a workaround): keep SEI
porosity change **on** for every row, including "Healthy baseline," but at *normal*
(un-accelerated) kinetics rather than the tuned recipe's `SI_MULT=18`/boosted rates. This
matches how a real cell actually behaves — even a "healthy," non-pathological cell has
*some* ordinary calendar/cycling SEI growth and porosity decline; the matrix's real purpose
per B.1 is isolating the *accelerated/pathological* mechanisms (Si cracking, dry-out) on
top of that normal baseline, not testing a literally-zero-porosity-change idealisation.
Under this redefinition, pore buffering is meaningfully "true" (not a guaranteed no-op) for
**all 5 rows**, satisfying the request, while "Healthy baseline"/"Silicon fatigue only"
still isolate what they're meant to (no *accelerated* SEI/porosity loss, no dry-out) —
just no longer isolating literally-zero SEI activity, which B.2 itself already flagged as
"the expected, physically correct outcome... not a limitation" for the old definition, so
this redefinition doesn't lose anything B.2 was relying on.

**If this resolution isn't wanted**: the fallback is B.2's original definition, with pore
buffering explicitly `"false"` (not "true", and not omitted) for rows 1–2 only, and a note
on any comparison plot that those two rows are structurally buffering-independent — a
smaller change to this plan but reintroduces the guard conflict if "true for all" is a hard
requirement.

### 2.3 Buffered vs. unbuffered reference, per row

Per instruction: pore buffering `"true"` (`"physical"` transition) for the primary run of
every row (using §2.2's resolution), **plus** a `"false"` reference run per row — mirroring
the established `test_pore_buffering/` buffered-vs-unbuffered convention throughout this
project. That's **10 runs** total (5 rows × 2 buffering states), not 5.

### 2.4 Stopping criterion: 50% SoH, with an explicit "no knee" fallback

Run to 50% SoH per row (reusing `physical_f0_width_grid_test.py`'s `current_soh()` helper
and its bug-fixed termination pattern — **do not** re-introduce the direct
`sol.cycles[0]`/`[-1]` indexing bug found and fixed there), with a generous cycle-count
safety cap (start from that script's `MAX_TOTAL_CYCLES=400`, adjust per row once Phase-1-
style smoke tests show actual timescales). **Explicitly expected, not a failure mode**: per
B.5's own success criteria, "Healthy baseline" and "Silicon fatigue only" may never
reach 50% SoH within any reasonable cycle count once SEI kinetics are at normal
(un-accelerated) rates (§2.2) — LAM-only fade (silicon fatigue row) can be very slow, and
normal-rate SEI/LLI fade alone may plateau well above 50% SoH within a practical cycle
budget. For rows that hit the safety cap without reaching 50% SoH, report final SoH reached
and stop there — a genuine, reportable result (matching B.5's own prediction that these two
rows show "flat capacity, flat/no-hump expansion" and "capacity fade from LAM... but no
hump", not necessarily a full 50%-SoH trajectory).

### 2.5 Output plan

One script (or one parameterised runner + one plotting pass), producing:
- Per-row: SoH vs. throughput/EFC, cell-level reversible-amplitude vs. throughput/EFC,
  internal transfer ratio `k`, buffered vs. unbuffered overlaid — same panel style as
  `physical_f0_width_grid_test.py`.
- One 5-row (or 10-line, buffered+unbuffered) comparison figure, matching the established
  `_comparison_*.png` convention from `test_pore_buffering/`.
- Explicit per-row annotation of whether a knee/hump was reached within the SoH target or
  safety cap, so "no knee" rows are visually distinguishable from "ran out of budget."
- Reuse Phase 2's EFC conversion from `../conditions_test_matrix/` (`EFC = Throughput / (2
  × initial capacity)`) for axis consistency across both test-matrix efforts, if practical.

### 2.6 Success criteria (per B.5, restated concretely)

- *Healthy baseline*: flat-ish capacity fade (normal-rate SEI only), flat/no-hump
  expansion — buffering should show as a near-constant partition (` k` pinned near `f0`)
  since normal-rate porosity decline is slow enough that the transition band is unlikely to
  be reached within the safety cap.
- *Silicon fatigue only*: capacity fade now has a LAM contribution on top of normal-rate
  LLI; expansion should still stay close to flat/monotonic — **no hump expected**, since
  without *accelerated* SEI/porosity loss there's still no fast `ε_struct` decline to
  trigger the buffering transition early. A directly falsifiable prediction, worth
  checking explicitly rather than assumed.
- *Pore clogging, Si intact*: the first row where a hump is expected — driven by
  accelerated porosity decline alone, no cracking contribution. Compare its knee
  timing/expansion shape against the fully-tuned baseline (`test_pore_buffering/`'s own
  validated result) as a sanity check that removing cracking alone doesn't collapse the
  expansion signal.
- *Dry-out, Si intact*: as above, plus whatever Phase 1c's validation sweep found about
  dry-out's effect on knee timing and expansion shape — this row is where that finding
  gets exercised inside the full matrix context (alongside the other rows) rather than in
  isolation.
- *Fully coupled*: the closest analogue to the real experimental cell — main point of
  comparison against `../conditions_test_matrix/`'s baseline condition and, ultimately,
  the real data.

## 3. Open decisions to confirm before implementation

1. §1a: accept the "switch to `\"ec reaction limited\"`, keep other constants, check how
   far off" cheap first pass, or commit directly to a full `k_sei`/`D_ec`/`c_ec_0`
   recalibration? (The wrapper-level extension enabling this choice is already done —
   this decision is now purely about the *recipe*-side calibration approach.)
2. §2.2: confirm the "Healthy baseline"/"Silicon fatigue only" row redefinition
   (normal-rate SEI on, not literally zero) — this changes what those two rows mean
   relative to B.2's original table.
3. §2.4: confirm the 50%-SoH-or-safety-cap approach and that a "did not reach 50%" result
   is an acceptable, reportable outcome for some rows rather than something to force by
   extending the cap indefinitely.
4. Whether Phase 1's dry-out/composite validation sweep (1c) should be written up as its
   own standalone result (mirroring `ruihe_dryout_validation_fig3.png`/`_ratios.png`) before
   folding dry-out into the full Phase-2 matrix, or done in one pass — recommend
   standalone first, same reasoning as `conditions_test_matrix_plan.md`'s phased approach:
   validate the new integration in isolation before composing it with everything else.
