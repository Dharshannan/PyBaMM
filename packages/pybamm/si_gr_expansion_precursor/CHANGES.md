# Core PyBaMM Changes Log — Si/Gr Expansion Precursor Project

Tracks every change made to `packages/pybamm/src/pybamm/**` (the actual PyBaMM
library, not the analysis scripts in `si_gr_expansion_precursor/`) over the
course of this project, in chronological order. Update this file whenever a
core-library change lands — it's the single place to check "what did we
actually modify in PyBaMM itself, and why" without having to reconstruct it
from conversation history or `git log`.

Each entry: **what** changed, **why**, **files touched**, and whether it's a
new `Parameter` (so it needs a default in whichever parameter set(s) use it).

---

## 1. SEI reaction-exponent soft cap (`exponent_max_sei`)

**Why**: the "reaction limited" SEI growth law's exponent term
(`exp(-alpha_SEI * F/RT * eta_SEI)`) could blow up under large overpotential
swings (e.g. a big current step), causing `IDA_ERR_FAIL`/`IDA_CONV_FAIL`
solver failures. A smooth cap fixes the numerical failure and, as a side
effect, turned out to control how sharp the degradation "knee" is: 10
recovers/flattens post-knee behaviour, 13 overshoots into a near-vertical
cliff.

**What**: `capped_exponent = exponent_max * tanh(exponent / exponent_max)`
replaces the raw exponent in the `j_sei` formula. For `|exponent| <<
exponent_max` this is ~identical to the uncapped law (`tanh(x) ≈ x`), so
normal SEI growth is unaffected until `eta_SEI` swings extreme.

**New parameter**: `exponent_max_sei` (`"{Primary/Secondary}: SEI reaction
exponent cap"`), phase-level, added to `ParticleLithiumIonParameters`.
Default 11.0 in most sets; this project uses 10.0 (see the winner scripts).

**Files**:
- [src/pybamm/models/submodels/interface/sei/sei_growth.py](../src/pybamm/models/submodels/interface/sei/sei_growth.py) — the `tanh` cap, in the `"reaction limited"` branch of `get_coupled_variables`.
- [src/pybamm/parameters/lithium_ion_parameters.py](../src/pybamm/parameters/lithium_ion_parameters.py) — `ParticleLithiumIonParameters`, `self.exponent_max_sei = pybamm.Parameter(f"{pref}SEI reaction exponent cap")`.

---

## 2. Anti-Paris-law stress relief in crack propagation

**Why**: the original Paris-law crack growth rate keeps accelerating crack
growth even as a crack gets very large relative to the particle, which is
physically wrong (a crack approaching the particle radius should see reduced
driving stress, not runaway growth) and contributed to solver instability at
deep cycling.

**What**: added a `stress_relief = max(1 - l_cr/R_typ, 0)` factor multiplying
`stress_t_surf` before it enters the stress-intensity-factor (`dK_SIF`)
calculation, so the effective driving stress relaxes toward zero as the crack
length approaches the typical particle radius. Original unbounded Paris law
left commented-out in place (marked `TODO: Anti-Paris law update`) for easy
reference/rollback.

**New parameters**: none (reuses existing `R_typ`).

**Files**:
- [src/pybamm/models/submodels/particle_mechanics/crack_propagation.py](../src/pybamm/models/submodels/particle_mechanics/crack_propagation.py) — `get_coupled_variables`, just before `dK_SIF`.

---

## 3. Porosity floor as a real parameter (`epsilon_min`)

**Why**: `reaction_driven_porosity.py` had a hardcoded `eps_min =
pybamm.Scalar(0.08)` softplus floor with a `TODO` comment noting it
suppresses the genuine resistance-increase/knee-sharpening that a truly low
porosity would cause via the Bruggeman transport relations. Making it a real
parameter let it be swept (0.005–0.08) to find the stability/sharpness
trade-off — landed on 0.01 for the accepted recipe.

**What**: `eps_min = pybamm.Scalar(0.08)` → `eps_min = domain_param.epsilon_min`,
softplus-smoothed floor unchanged otherwise.

**New parameter**: `epsilon_min` (`"{Domain} electrode porosity floor"`),
domain-level, in `DomainLithiumIonParameters`. Defaulted to 0.08 (preserving
old behaviour) in **all 15** lithium-ion parameter sets (Mayur2024, Chen2020,
Chen2020_composite, Ai2020, Ecker2015, Ecker2015_graphite_halfcell,
Marquis2019, Mohtat2020, MSMR_example_set, OKane2022, NCA_Kim2011,
OKane2022_graphite_SiOx_halfcell, Prada2013, Xu2019, ORegan2022,
Ramadass2004) via a bulk regex-insertion script, so every existing parameter
set keeps working unchanged. This project's scripts override it to 0.01.

**Files**:
- [src/pybamm/parameters/lithium_ion_parameters.py](../src/pybamm/parameters/lithium_ion_parameters.py) — `DomainLithiumIonParameters._set_parameters()`.
- [src/pybamm/models/submodels/porosity/reaction_driven_porosity.py](../src/pybamm/models/submodels/porosity/reaction_driven_porosity.py) — softplus floor now reads the parameter.
- All 15 `src/pybamm/input/parameters/lithium_ion/*.py` files listed above — `"{Domain} electrode porosity floor": 0.08` added.

---

## 4. `si_gr_expansion` parameter set (new)

**Why**: needed a parameter set that reports a physically realistic
particle/electrode expansion signal for the Si secondary phase, using a
*measured* V/V₀-vs-stoichiometry curve instead of the analytic Ai2020/
Bonkile2024 linear model that Mayur2024 uses.

**What**: exact copy of `Mayur2024.py`, differing only in
`"Secondary: Negative electrode volume change"`: now
`silicon_volume_change_JiaGuo`, an interpolant built from a measured CSV
(`process_1D_data` + `pybamm.Interpolant`), instead of
`silicon_volume_change_Ai2020`. `"Primary: Negative electrode volume
change"` (graphite, Ai2020) is unchanged. Verified to reproduce **identical**
degradation numbers (knee, LAM, LLI) to Mayur2024 at the same recipe — the
volume-change function only feeds the reported thickness-change diagnostic
(`t_change`), not the stress equations (which use the separate `Omega`
partial-molar-volume parameter), so it doesn't touch degradation physics.

**New parameter data, no new Parameter definitions**: adds the BoL
pore-buffering defaults from item 6 below (this parameter set only, per
explicit instruction — not bulk-added to the other 14 sets).

**Files**:
- [src/pybamm/input/parameters/lithium_ion/si_gr_expansion.py](../src/pybamm/input/parameters/lithium_ion/si_gr_expansion.py) (new file).
- [src/pybamm/input/parameters/lithium_ion/data/volume_vs_sto_lithiation.csv](../src/pybamm/input/parameters/lithium_ion/data/volume_vs_sto_lithiation.csv) (new data file, ~250 rows, `Sto,Volume` columns, V/V₀ ratio peaking ~3.454 near sto≈0.963).
- `pyproject.toml` — `si_gr_expansion` entry under `[project.entry-points."pybamm_parameter_sets"]`.
- Installed package `entry_points.txt` — same entry patched directly (editable install needs both to pick up a new parameter set without reinstalling).

---

## 5. Pore-buffering (volume partition) submodel

**Why**: implements the "pore buffering" mechanism from
`Si_Gr_Expansion_Precursor/Notes_Papers/si-gr_expansion_precursor_notes.PDF`
("Making the Pore-Buffering Model Physical") — partitions active-material
swelling between pore-volume buffering and electrode thickness change,
instead of routing 100% of swelling to thickness (the previous behaviour).
Full design rationale in [pore_buffering_implementation_plan.md](pore_buffering_implementation_plan.md).
Scope: only the volume-partition mechanism (plan §5/§6, Changes 1–4 and 6).
Explicitly **excludes** the pressure-coupling phase (plan §3/§4, Change 5 —
`Stack pressure [Pa]` feeding the particle-stress BC, Yin et al. `K_stack`/
`C_pore(ε)` stack model) — deferred to a later phase.

**What**, four pieces:

1. **New model option** `"pore buffering"`: `"false"` (default, original
   behaviour) / `"true"`. Requires `"SEI porosity change"` or `"lithium
   plating porosity change"` to be `"true"` (validated, raises
   `pybamm.OptionError` otherwise) since the partition lives inside the
   reaction-driven porosity submodel, which is only built in that case.
2. **Three new domain-level parameters**, distinct from the existing
   `epsilon_min` (numerical/SEI-growth solver floor) — these describe the
   *mechanical* pore-network percolation/closure behaviour:
   - `eps_min_transfer` (`"{Domain} electrode pore buffering closure
     porosity"`) — the porosity at which the pore network's throat
     percolation closes (all further swelling forced to thickness).
   - `eps_max_transfer` (`"{Domain} electrode pore buffering upper
     porosity"`) — the upper "knee" porosity above which the transmitted
     fraction plateaus at `f0`.
   - `f_transmit_min` (`"{Domain} electrode transmitted fraction
     plateau"`) — `f0`, the fraction of swelling always transmitted to
     thickness even with abundant pore space.
3. **Build-order change**: `set_porosity_submodel()` moved to *after*
   `set_crack_submodel()` in `set_submodels()`, so the porosity submodel can
   read the per-phase electrode thickness-change variables particle
   mechanics just computed. Safe unconditionally (nothing between the old
   and new position consumes porosity).
4. **The partition itself**, in `ReactionDriven.get_coupled_variables`,
   scoped to the **negative electrode only** (the Si/Gr precursor physics is
   anode-specific — matches the notes' and Yin et al.'s framing, and avoids
   needing positive-electrode buffering parameters that don't exist yet):
   - Exposes `"{Domain} electrode structural porosity"` unconditionally
     (irreversible-only porosity: SEI/plating/dead-Li/crack terms, no
     buffering) — a harmless diagnostic even when buffering is off (equals
     total porosity in that case).
   - When `"pore buffering" == "true"` and the negative electrode has a
     particle-mechanics thickness-change output to partition: computes the
     three-regime transmitted fraction `f(eps_struct)` (Eq. 18 of the
     notes), applies the hard geometric ceiling (`min((1-f)*dv_solid,
     headroom)`, Eq. 19), and overwrites `"Negative electrode thickness
     change [m]"` with the *transmitted* portion only. Exposes
     `"Negative electrode buffered volume change"`, `"Negative electrode
     solid volume change"`, and `"Negative electrode transfer ratio k"`
     (read-off diagnostic only, Eq. 23 — never fed back into any equation).
   - Falls through to the original unbuffered `eps_struct` when the option
     is off, or when there's no thickness-change variable to partition
     (e.g. `"particle mechanics": "none"`).

**BoL parameterisation used** (per explicit instruction — assumed for now,
pending real porosimetry/dilatometry characterization): `k_BoL = 0.7`,
`eps_min_transfer = 0.12`. Since `eps_max_transfer` is set equal to the
electrode's own known BoL porosity (`ε_init = 0.25`), the three-regime
formula collapses algebraically at BoL to `f0 = k_BoL` — no fitting needed.
Literature check for the `eps_min_transfer` assumption (0.08–0.12 range) is
in the plan doc §3.3; defaults added to `si_gr_expansion.py` **only** (not
bulk-added to other parameter sets, unlike item 3's porosity floor):
```
"Negative electrode pore buffering closure porosity": 0.12,   # eps_min_transfer (assumed)
"Negative electrode pore buffering upper porosity": 0.25,     # eps_max_transfer = eps_init
"Negative electrode transmitted fraction plateau": 0.7,       # f0 = k_BoL
```

**Bugfix found during the degradation-cycle test** (not in the smoke test,
which is too short to reach it): `eps_struct` is the *raw* (unfloored)
irreversible porosity, which can go negative under heavy late-life
degradation (previously always hidden from every consumer by the numerical
`epsilon_min` softplus floor, applied only at the very end). Once `eps_struct
< eps_min_transfer`, `headroom = eps_struct - eps_min_transfer` went
negative, and `pybamm.minimum((1-f)*dv_solid, headroom)` then picked the
negative headroom even when `f` was already correctly clamped to 1 (making
`(1-f)*dv_solid` exactly 0) — injecting a spurious large-negative
`dv_buffered` that corrupted `k` (observed reaching -4.98) and the reported
thickness change late in life. Fixed by flooring `headroom` at 0 via
`pybamm.maximum(eps_struct_avg - eps_min_transfer, 0)` before the `minimum`
call. After the fix, `k` correctly saturates at exactly 1.0 and stays there
for the rest of life. The `eps_struct` diagnostic itself still legitimately
goes negative late in life (it's intentionally unfloored, matching what the
model was always doing internally before this project exposed it) — that's
expected, not a bug.

**Result** (252-cycle run to 50% SoH, same recipe as the accepted winner):
`k` rises from ~0.71 at BoL and saturates at exactly 1.0 around throughput
≈580–650 A.h — well before the capacity knee (65.3% SoH at throughput
1583.5 A.h, sharpness 2.87) — a genuine early-precursor signal, matching the
notes' central claim. The negative-electrode within-cycle expansion
amplitude shows a clear hump (14.5 → peaks ~17.5 µm right where `k`
saturates → declines to ~7.5 µm by 50% SoH), even though `k` itself just
saturates and stays flat rather than rising-then-falling — the hump in the
*reported thickness* comes from LAM shrinking the underlying swelling
capacity (`dv_solid`) after `k` has already saturated, not from `k`'s own
shape. **Caveat**: a true rise-then-fall hump in `k` itself needs the
deferred mechanistic-`f`/pressure-coupling phase (plan §5.6/§3/§4) — the
reduced three-regime `f(ε)` implemented here can only saturate and hold once
`ε_struct` passes `eps_min_transfer`, it has no mechanism to un-saturate.
Also note the knee moved much earlier than the non-buffered baseline
(65.3% vs ~82%) — expected, since buffering changes the porosity trajectory
feeding transport, not a regression, but the original recipe was tuned to
hit 82–93% *without* buffering, so it would need re-tuning if matching that
band under the new physics is desired.

**Second fix (design correction, per explicit instruction)**: the first cut
of `eps_k` also clamped the *reported/transport* porosity at
`eps_min_transfer` (`eps_k = max(eps_struct - dv_buffered,
eps_min_transfer)`), matching the notes' literal Eq. 21. In practice this
meant that once buffering saturates late in life (`dv_buffered -> 0`),
`eps_min_transfer` (0.12, the literature-informed BoL placeholder) silently
became the effective porosity floor for transport, **overriding** the
separately-tuned numerical floor (`"Negative electrode porosity floor"`,
0.01) this recipe's knee/sharpness was originally tuned against. Per this
session's earlier floor-sweep investigation, a higher floor gives an
earlier/less-sharp knee (0.08 → 67.0%, ~0.01 → 82.5%) — exactly matching the
65.3% seen with the un-corrected version. Corrected: `eps_min_transfer` now
only gates the *partition* (via `f` and `headroom`) — it is not a floor on
`eps_k`. `eps_k = eps_struct - dv_buffered` (unclamped at the
`eps_min_transfer` level), so actual porosity keeps decreasing with
`eps_struct` once buffering saturates, bounded only by the existing
numerical softplus floor. This deliberately departs from the notes' literal
Eq. 21 for this implementation.

**Result after the design correction** (same 252-ish-cycle run, defaults
unchanged — `eps_min_transfer=0.12` etc., no parameter retuning needed):
knee restored to **82.29% SoH** at throughput 805.5 A.h, sharpness 6.76 —
matching the non-buffered baseline (82.48%) closely. The cell-level
within-cycle expansion amplitude now shows the intended dip-then-rise hump
shape (starts ~1.08 µm → dips to ~0.4–0.5 µm through mid-life → rises
sharply to a peak ~1.9–2.0 µm right at/past the knee → settles to ~0.8–0.9
µm by 50% SoH), closely matching the LAM-driven hump observed in the
earlier non-buffered `test_expansion` runs (~1.1 → 0.42 → 1.95 → 0.9–1.2
µm) — now reinforced by `k` saturating at throughput ≈580 A.h, just before
the knee, so both mechanisms compound around the same point in life instead
of appearing as two disjoint features.

**Verified** (smoke test, `si_gr_expansion_precursor/test_pore_buffering/smoke_test.py`):
- `"pore buffering": "false"` reproduces the original behaviour exactly
  (structural porosity == total porosity, no buffering diagnostics present).
- `"pore buffering": "true"` solves; `k` starts at exactly 0.7000 at t=0
  (matches the BoL algebra) and stays ≤ 1; buffered negative-electrode
  thickness change is ~70.5% of unbuffered (matches `f ≈ 0.70–0.71`
  expected from the formula); within-cycle cell-level expansion amplitude
  (charge-end minus same-cycle discharge-end, same convention as
  `test_expansion`) stays positive and is damped by buffering too, though
  much more subtly at cell level (~1.7% reduction) than at the negative
  electrode alone (~30%) — because `"Cell thickness change [m]"` combines
  negative + positive electrode + thermal, and buffering only touches the
  negative electrode's own contribution.

**Files**:
- [src/pybamm/models/full_battery_models/base_battery_model.py](../src/pybamm/models/full_battery_models/base_battery_model.py) — option definition, docstring, validation.
- [src/pybamm/models/full_battery_models/lithium_ion/base_lithium_ion_model.py](../src/pybamm/models/full_battery_models/lithium_ion/base_lithium_ion_model.py) — `set_submodels()` reorder.
- [src/pybamm/parameters/lithium_ion_parameters.py](../src/pybamm/parameters/lithium_ion_parameters.py) — `eps_min_transfer`, `eps_max_transfer`, `f_transmit_min`.
- [src/pybamm/models/submodels/porosity/reaction_driven_porosity.py](../src/pybamm/models/submodels/porosity/reaction_driven_porosity.py) — the partition, `_transmitted_fraction`, k diagnostic.
- [src/pybamm/input/parameters/lithium_ion/si_gr_expansion.py](../src/pybamm/input/parameters/lithium_ion/si_gr_expansion.py) — BoL defaults.
- [si_gr_expansion_precursor/test_pore_buffering/](test_pore_buffering/) — smoke test + short degradation cycling validation. `pore_buffering_degradation_test.py` now runs **two** full degradation simulations (buffering off as a reference, then on) and produces three figures: the main SoH/k/amplitude diagnostic, a volume-vs-capacity figure (cell-level expansion amplitude, per-electrode contributions [negative buffered vs. positive unbuffered — they largely cancel, which is why the cell-level signal is ~10x smaller than either electrode alone], and capacity fade, all sharing the throughput x-axis), and a porosity "breathing" diagnostic (full time-resolved negative-electrode porosity within four representative ageing cycles at different life stages, showing the within-cycle reversible oscillation shrinking to exactly zero once `k` saturates at 1 — direct visual confirmation the partition's reversible/breathing component behaves as intended).
  - **Two distinct `k` definitions are now plotted together, per explicit clarification — these are NOT the same quantity**: (1) **internal/mechanistic k** — the partition's own read-off diagnostic (`"Negative electrode transfer ratio k"`, Eq. 23), computed entirely within the negative electrode from `dv_thickness/dv_solid`; rises smoothly from ~0.71 and saturates at exactly 1.0 once the buffer is exhausted. (2) **observable k = δ_cell/δ_particle** (the notes' original Eq. 1-2, with the "rev" qualifier read as within-cycle amplitude) — the buffered run's cell-level within-cycle amplitude divided by the unbuffered reference run's cell-level amplitude at the same throughput; this is what you'd actually compute from two dilatometry campaigns (or here, two simulations), not a model-internal state.

  - **Third bugfix, found while investigating why the observable k started near 1 instead of the expected ~0.7**: `"Cell thickness change [m]"` is built by `base_mechanics.py`'s `_aggregate_cell_thickness_change`, which runs *during the particle-mechanics submodels' own `get_coupled_variables` call* — i.e. **before** porosity/buffering runs (even after the build-order reorder, mechanics still finishes before porosity starts). It reads whatever `"Negative electrode thickness change [m]"` value exists in `variables` at that moment — the *unbuffered* one. PyBaMM expressions are immutable: overwriting `variables["Negative electrode thickness change [m]"]` later inside `ReactionDriven` does **not** retroactively propagate into the already-built `"Cell thickness change [m]"` expression, which stays permanently wired to the stale unbuffered value. Diagnosed by computing the "implied thermal amplitude" (`cell_amplitude - neg_amplitude - pos_amplitude`, which should be ~0 for an isothermal run) directly: it came out exactly 0.0000 µm for the unbuffered run but a spurious 6.22 µm for the buffered run — proof `"Cell thickness change"` wasn't using the buffered negative-electrode value. **Fixed** by having the buffering branch explicitly correct `"Cell thickness change [m]"` with a delta (`− thickness_change_unbuffered + thickness_change_buffered`) rather than re-deriving the full `neg + pos + thermal` formula (avoids duplicating/drifting from that formula).
    - **This changes the volume-vs-capacity and k plots materially.** Before the fix, buffered and unbuffered cell-level amplitude were nearly identical (~1.2 µm vs ~1.18 µm, ratio ~0.98) purely because the buffered run's cell signal was accidentally still using the unbuffered value. After the fix: unbuffered cell amplitude ≈ +1.2 µm, buffered ≈ **−5.0 µm** at BoL — the sign flips. This is because `δ_cell` is a near-total cancellation between the negative electrode (~+20.7 µm unbuffered) and positive electrode (~−19.5 µm) contributions, leaving only a ~+1.2 µm residual; damping just the negative electrode by the BoL transfer ratio (~0.70, to ~+14.5 µm) is enough to flip the sign of that already-small residual. This is a genuine, physically real property of this composition (anode and cathode swelling nearly cancelling at cell level) revealed by the fix — not a new bug — but it does mean the "observable" `k = δ_cell/δ_particle` is a numerically ill-conditioned diagnostic early in life for this specific composition.
    - **The full corrected (but still pre-item-6) trajectory was coherent**: `k_observable` started at ≈−4.5 at BoL, rose steadily, crossed zero around throughput ≈380–400 A.h, and converged to track the internal k (≈1) right around where buffering saturates (~580 A.h) — because once buffering saturates, the negative electrode's own contribution simply equals its unbuffered value again, so the two runs' cell-level signals converge for the rest of life.
  - **Fourth fix — `δ_particle` was defined wrong** (caught by inspecting a reference dilatometry/`k`-vs-EFC slide from the actual experiment this project models): re-reading the notes' own §4.4 definition, `δ_particle` is the **unbuffered negative electrode's own** within-cycle swelling (`n·L·∫Ω dc` for the anode, "cathode neglected" per Yin et al.) — **not** a second whole-cell simulation's cell-level amplitude. Using a second simulation's cell-level amplitude (as the third-bugfix version above did) wrongly pulled the positive electrode's own, unrelated, always-unbuffered swelling into the denominator. Corrected: `δ_particle` is now computed from the SAME buffered run, using `"Negative electrode solid volume change"` (`dv_solid`, the pre-partition unbuffered swelling in volume-fraction units) scaled to thickness — no second simulation needed for `k` (the unbuffered reference run is now kept only as an independent SoH/knee cross-check). This narrowed `k_observable`'s range substantially (from ≈[−4.5, 1.1] to ≈[−0.24, 0.13]) but didn't fully fix it, because the *numerator* (`δ_cell`) was still corrupted by item 6's cathode over-sizing.
  - **See item 6 below** for the cathode-calibration fix that, combined with this `δ_particle` correction, produced the final sane result: `k_observable` in a well-behaved 0.59–0.89 range, tracking the same shape as the internal k, and both `δ_cell`/electrode-contribution plots now qualitatively matching the reference experimental data (gradual rise-then-decline dominated by the anode, cathode contribution small and near-flat).

---

## 6. Cathode (`volume_change_Ai2020`) recalibration

**Why**: found while investigating why the corrected `k = δ_cell/δ_particle`
(item 5) still didn't look physically sane. `"Positive electrode partial
molar volume [m3.mol-1]"` (inherited unchanged from Mayur2024, itself from
Ai2020/Ai2019) gives `volume_change_Ai2020(sto) = Ω·c_s_max·sto` a value of
**0.79 at sto=1** — i.e. the model's cathode particle swells by ~79% of its
own volume over the full stoichiometry window. A reference experimental
slide (particle-expansion decomposition + `k`-vs-EFC panels from the actual
cell this project models) shows the real cathode (NMC) contributing an
essentially flat, near-**zero** particle expansion — `δ_rev,particle` in
that data is completely dominated by the anode (Si+graphite), matching the
notes' own "cathode neglected" framing (§4.1, Yin et al.). The same
reference data also confirms this project's other targets independently:
`k` starts in a ~0.85–0.9 band (not far from the `k_BoL=0.7` assumption used
for this project) and — only near end-of-life — **exceeds 1** ("potential
pore closure"), which the current volume-partition-only implementation
cannot reproduce by construction (`k ≤ 1` always, Eq. 23) since that late-
life overshoot needs additional expansion sources (lithium plating swelling,
pressure-driven cracking feedback) that belong to the deferred pressure-
coupling phase, not a bug in what's implemented today.

**What**: reduced `"Positive electrode partial molar volume [m3.mol-1]"`
from `1.25e-05` to `1.584686e-06` (≈7.9× smaller), targeting a ~10% max
cathode particle volume change (per explicit instruction, based on the
reference data's observed ceiling) instead of ~79%. Scoped to
`si_gr_expansion.py` only, per explicit instruction — `Mayur2024.py` and all
other parameter sets keep the original (inherited, seemingly-never-
validated-for-this-diagnostic) value.

**Side-effect risk**: this parameter is shared with the positive electrode's
*stress* equations (`Omega` in `base_mechanics.py`'s `stress_t_surf`/
`disp_surf`), not just the thickness-change diagnostic, and
`"loss of active material": "stress-driven"` applies to the positive
electrode too — so this could in principle shift the cathode's own
stress-driven LAM rate. Verified this doesn't meaningfully disturb the
already-tuned degradation recipe: knee moved from 82.29% to **81.98%** SoH
(throughput 805.5 → 790.7 A.h, sharpness 6.76 → 6.49) — well within the
82–93% target band and closely tracking the unbuffered reference run's own
knee (82.03%), so the composition/degradation recipe did not need
re-tuning.

**Result**: positive-electrode within-cycle amplitude dropped from
~−19.5 µm to **~−1.2 to −2.5 µm** (matching the reference data's near-
negligible cathode contribution); cell-level amplitude is now **cleanly
positive throughout** (12.0 → 9.4 µm, peak 15.3 µm, no more sign flips) and
closely tracks the negative electrode's own shape (just scaled down);
`k_observable` now sits in a **0.59–0.89** range, rising with the same
timing/shape as the internal k rather than diverging wildly from it. The
electrode-contributions and volume-vs-capacity plots now qualitatively match
the reference data's own `Expansion [µm]` panel shape (gradual rise then
decline, anode-dominated).

**Files**:
- [src/pybamm/input/parameters/lithium_ion/si_gr_expansion.py](../src/pybamm/input/parameters/lithium_ion/si_gr_expansion.py) — `"Positive electrode partial molar volume [m3.mol-1]"`.
- [si_gr_expansion_precursor/test_pore_buffering/pore_buffering_degradation_test.py](test_pore_buffering/pore_buffering_degradation_test.py) — `δ_particle` redefinition (item 5's fourth fix) lives here, not in core PyBaMM.

---

## 7. Knee-timed band retune + smooth (tanh) transmitted fraction

**Why (part 1 — band retuning)**: with the item 5 BoL-derived defaults
(`eps_min_transfer=0.12`, `eps_max_transfer=eps_init=0.25`), the transmitted
fraction started ramping immediately from BoL and saturated around
throughput ≈580 A.h — well before the knee (~800 A.h) — giving the
`"Cell thickness change"` diagnostic a smooth rise starting from BoL rather
than the reference dilatometry pattern (flat/stable expansion through most
of life, then a **sharp rise right at the knee**, then decline — see the
`k`-vs-EFC and `Expansion [µm]` panels in the reference slide). Swept
`(eps_min_transfer, eps_max_transfer)` as a narrow (~0.01 wide) band shifted
progressively lower via `eps_max_transfer_sweep.py`; found the transition
timing (and, via the porosity→transport feedback loop, the knee itself)
both shift with the band position. `[0.08, 0.09]` gave the closest
alignment: buffering-saturation transition ≈59 A.h before the knee, knee at
81.17–81.98% SoH (essentially unchanged from the pre-retune 81.98–82.29%
range), matching the reference pattern's qualitative shape closely for the
first time.

**Why (part 2 — smoothing)**: once tuned to a narrow band, the original
three-regime piecewise-linear `f(eps_struct)` (Eq. 18 of the notes) —
continuous in value but with kinks (discontinuous derivative) at
`eps_min`/`eps_max` — produced a visually and numerically near-step-function
transition (flagged directly: "this seems like there is a jump/
discontinuity"). Compared three candidate forms analytically first (no
simulation needed — a quick standalone numpy/matplotlib script, not kept,
producing [f_form_comparison.png](test_pore_buffering/f_form_comparison.png)):
the original piecewise-linear (visible slope-discontinuity, confirmed via
its derivative plot), a tanh/sigmoid blend (smooth everywhere, symmetric,
tightly bounded to the tuned band), and a compliance-ratio form derived from
the notes' §5.6 "mechanistic f" (`f = 1/(1 + K·C_pore(ε))`, physically
motivated but asymmetric with influence extending well outside the intended
band, and needing an extra smooth floor to avoid its own kink). Picked the
tanh blend: same three parameters (`eps_min_transfer`/`eps_max_transfer` now
interpreted as the ~90%/10% points of the transition rather than hard
cutoffs, `f_transmit_min` unchanged), no min/max clamping needed (tanh is
already bounded), C-infinity smooth everywhere.

**What**:
```python
eps_mid = (eps_min + eps_max) / 2
steepness = 2.1972 / (eps_max - eps_min)  # 2*arctanh(0.8): band edges sit at the 90%/10% points
blend = 0.5 * (1 - pybamm.tanh(steepness * (eps_struct - eps_mid)))
f = f0 + (1 - f0) * blend
```
replaces the old `f_shared = 1 - (1-f0)*(eps-eps_min)/(eps_max-eps_min)` +
`pybamm.minimum(pybamm.maximum(f_shared, f0), 1)` clamp.

**New defaults in `si_gr_expansion.py`**:
```
"Negative electrode pore buffering closure porosity": 0.08,   # was 0.12
"Negative electrode pore buffering upper porosity": 0.09,     # was 0.25 (=eps_init)
"Negative electrode transmitted fraction plateau": 0.7,       # unchanged
```
`f0 = k_BoL = 0.7` is still exact at BoL regardless of this change:
`eps_struct(BoL) = eps_init (0.245) >> eps_max_transfer`, so the blend is
saturated at its `f0` asymptote at BoL either way.

**Result** (same paired buffered/unbuffered run as item 6, now with the
retuned band + smooth transition): knee at **81.16% SoH**, throughput 830.1
A.h, sharpness 6.75 (vs the unbuffered reference's 82.03%/6.37 — still
closely tracking). Cell-level expansion amplitude: flat/declining
11.9→9.7 µm through throughput ≈0–700 A.h, smooth rise to a peak of
**14.6 µm right at/just after the knee** (throughput ≈800 A.h), decline to
9.3 µm by 50% SoH — a clean match to the reference `Expansion [µm]` panel's
shape. Both `k` definitions now show a smooth S-curve (internal k:
0.70→1.00; observable `k = δ_cell/δ_particle`: 0.58→0.88) instead of the
sharp near-vertical jump from the un-smoothed version.

**Caveat, still open**: this is now three empirically-tuned parameters
(`eps_min_transfer`, `eps_max_transfer`, and implicitly the tanh steepness
convention) fit to match the *shape* of one reference dataset — not a
measured or first-principles value. If real porosimetry/dilatometry data
for this specific electrode becomes available, re-derive
`eps_min_transfer` from that per the plan doc §3.4 rather than treating this
tuned value as physical.

**Files**:
- [src/pybamm/models/submodels/porosity/reaction_driven_porosity.py](../src/pybamm/models/submodels/porosity/reaction_driven_porosity.py) — `_transmitted_fraction`, now the tanh blend.
- [src/pybamm/input/parameters/lithium_ion/si_gr_expansion.py](../src/pybamm/input/parameters/lithium_ion/si_gr_expansion.py) — retuned `eps_min_transfer`/`eps_max_transfer` defaults.
- [si_gr_expansion_precursor/test_pore_buffering/eps_max_transfer_sweep.py](test_pore_buffering/eps_max_transfer_sweep.py) — the band-position sweep script (kept for reference; overwritten in place between v1/v2 passes, so it reflects the final v2 grid).

---

## 8. "pore buffering transition" model option: tanh vs. physical

**Why**: per explicit request, made the choice between the two
`_transmitted_fraction` forms (item 7) a real, switchable model option
instead of picking one — the tanh blend is semi-empirical (fits the shape,
doesn't derive it), while a properly-worked-out compliance-ratio form is
more defensible as "the physical version." Wanted to actually run and
compare both, not just implement blind.

**What**: new option `"pore buffering transition"`: `"tanh"` (default,
item 7's blend, renamed `_transmitted_fraction_tanh`) or `"physical"` (new
`_transmitted_fraction_physical`). Only read when `"pore buffering"` is
`"true"`. The physical form:
```python
headroom = softplus(eps_struct - eps_min_transfer, 0, 100)   # smooth max(., 0)
C_pore = 1 - exp(-headroom / (eps_max_transfer - eps_min_transfer))
K_Cmax = (1 - f0) / f0   # solved so f(eps -> infinity) = f0 exactly
f = 1 / (1 + K_Cmax * C_pore)
```
derived from the notes' §5.6 stiffness-balance `f = 1/(1 + K·C_pore(ε))`,
with `C_pore(ε)` modelled as a smooth function of headroom above closure
that saturates at a finite value (rather than growing unboundedly) — that
saturation is what gives the `f0` plateau as a true limit rather than
decaying to 0. Reuses the same three parameters
(`eps_min_transfer`/`eps_max_transfer`/`f_transmit_min`) as the tanh
version, no new parameters needed. `si_gr_expansion_precursor/
test_pore_buffering/pore_buffering_degradation_test.py` gained a
module-level `PORE_BUFFERING_TRANSITION` switch (`"tanh"`/`"physical"`) and
transition-suffixed output filenames so both can be run and compared
side by side without overwriting each other's results.

**Result** (same paired run as item 7, `PORE_BUFFERING_TRANSITION =
"physical"`): knee at 80.72% SoH (throughput 849.7 A.h, sharpness 6.82) —
close to the tanh version's 81.16%/830.1/6.75. `k`, expansion amplitude, and
electrode contributions all land in the same order of magnitude as the tanh
run (internal k 0.70→1.00, observable k 0.58→0.88, cell-level amplitude
peak 14.2 µm vs tanh's 14.6 µm). The qualitative difference between the two
forms: the physical/compliance-ratio version's rise is **more gradual and
asymmetric** — it starts deviating from the flat baseline earlier
(≈throughput 550–650 vs tanh's ≈700) and rises as a smoother, rounder bell
curve rather than tanh's sharper, more localized S-curve — consistent with
the asymmetric shape seen in the standalone analytical comparison
([f_form_comparison.png](test_pore_buffering/f_form_comparison.png)). Both
forms peak close to (physical: ~50 A.h after; tanh: ~-30 to -60 A.h before)
the knee and are fully smooth (no kinks) — genuinely a shape choice at this
point, not a correctness question. Not yet decided which better matches the
reference dilatometry data's peak shape (looks visually rounder/bell-like
in the reference slide, which might favor the physical form, but this
hasn't been checked quantitatively).

**Files**:
- [src/pybamm/models/full_battery_models/base_battery_model.py](../src/pybamm/models/full_battery_models/base_battery_model.py) — option definition/docstring/defaults.
- [src/pybamm/models/submodels/porosity/reaction_driven_porosity.py](../src/pybamm/models/submodels/porosity/reaction_driven_porosity.py) — `_transmitted_fraction` dispatcher, `_transmitted_fraction_tanh`, `_transmitted_fraction_physical`.
- [si_gr_expansion_precursor/test_pore_buffering/pore_buffering_degradation_test.py](test_pore_buffering/pore_buffering_degradation_test.py) — `PORE_BUFFERING_TRANSITION` switch, transition-suffixed output filenames.

---

## 9. Reparameterisation: `eps_max_transfer` → `eps_transfer_width`

**Why**: caught by direct question — with the item 7 retuned band
(`eps_min_transfer=0.08`, `eps_max_transfer=0.09`) and `eps_init=0.245`,
`eps_init` sits nearly 3× above `eps_max_transfer`. The notes' original
Eq. 18 three-regime picture assumes the cell genuinely occupies three
life-stages — starting in the "shared" regime at BoL (`eps_init < eps_max`)
and only later reaching the plateau/closure regimes. That assumption no
longer holds here: the cell starts (and stays, for most of life) deep
inside the "maximal buffering" plateau, only reaching the transition region
very close to the knee. The three regimes have collapsed into two lived
stages (flat plateau, then late transition-to-closure) — which is exactly
what we *wanted* for the knee-aligned shape, but it means
`eps_max_transfer` no longer represents "a porosity threshold the cell
passes through as a distinct regime change." In every implementation
(piecewise, tanh, physical) its only remaining role had become "define a
width together with `eps_min_transfer`" — so it's renamed to say that
directly.

**What**: `eps_max_transfer` parameter removed; replaced with
`eps_transfer_width` (`"{Domain} electrode pore buffering transition
width"`), used directly wherever the old code computed
`eps_max_transfer - eps_min_transfer`:
- Tanh: `eps_mid = eps_min_transfer + width/2`, `steepness = 2.1972/width`.
- Physical: `eps_transfer_width` is directly the compliance decay length
  (`C_pore(eps) = C_max*(1 - exp(-headroom/width))`).

`eps_min_transfer` is unchanged (still a genuine physical reference point —
the closure/percolation porosity) and `f_transmit_min` is unchanged. This
is a pure rename/reparameterisation with no numerical change: `si_gr_
expansion.py`'s default became `"Negative electrode pore buffering
transition width": 0.01` (was `eps_max_transfer=0.09` with
`eps_min_transfer=0.08`, same difference). Verified via the smoke test and
a direct tanh/physical solve check: identical `k(t=0)=0.7000` and all other
numbers unchanged from before the rename.

**Also fixed while doing this**: `smoke_test.py` had a stale hardcoded
`0.12` threshold in one assertion (`eps_struct_on.min() >= 0.12`), left over
from before the item 7 retune to 0.08 — corrected to 0.08. It hadn't been
failing (a short 2-cycle run stays well above either value), but it wasn't
testing the actual current default either. `eps_max_transfer_sweep.py` and
`peak_height_sweep.py` (which actively set the old parameter name via
`param.update()`) were updated to convert their band-endpoint sweep
variables to a width before calling `param.update()` — without this fix,
re-running either script would have silently done nothing (the old
parameter name is no longer read by anything), making the sweep
meaningless while appearing to run fine.

**Files**:
- [src/pybamm/parameters/lithium_ion_parameters.py](../src/pybamm/parameters/lithium_ion_parameters.py) — `eps_transfer_width` parameter definition.
- [src/pybamm/models/submodels/porosity/reaction_driven_porosity.py](../src/pybamm/models/submodels/porosity/reaction_driven_porosity.py) — both `_transmitted_fraction_*` methods updated.
- [src/pybamm/models/full_battery_models/base_battery_model.py](../src/pybamm/models/full_battery_models/base_battery_model.py) — `"pore buffering transition"` option docstring updated.
- [src/pybamm/input/parameters/lithium_ion/si_gr_expansion.py](../src/pybamm/input/parameters/lithium_ion/si_gr_expansion.py) — default renamed/updated.
- [si_gr_expansion_precursor/test_pore_buffering/](test_pore_buffering/) — `pore_buffering_degradation_test.py`, `smoke_test.py` (docstrings + the stale `0.12` fix), `eps_max_transfer_sweep.py`, `peak_height_sweep.py` (actual `param.update()` calls fixed).

---

## Open items / not yet done

- Pore-buffering pressure-coupling phase (plan §3/§4/Change 5) — deferred.
  This is also what's needed to reproduce the reference data's late-life
  `k > 1` overshoot ("potential pore closure") — see item 6.
- `eps_min_transfer` is an assumed placeholder (0.12), not measured — needs
  BoL porosimetry/tomography on the actual electrode per the plan doc §3.4.
- Cathode `"Positive electrode partial molar volume"` (item 6) is now a
  ~10% max-particle-volume-change placeholder based on eyeballing the
  reference data's ceiling, not a literature/measured value — worth
  replacing with an actual NMC (or whatever this cathode chemistry is)
  dilatometry/XRD-derived value if one becomes available.
- Porosity-floor/`exponent_max_sei` defaults (item 3) were bulk-added to all
  15 parameter sets; pore-buffering defaults and the cathode recalibration
  (items 5, 6) deliberately were not — only `si_gr_expansion.py` has them,
  since this is specific to this project's Si/Gr composite, not a
  general-purpose default.
