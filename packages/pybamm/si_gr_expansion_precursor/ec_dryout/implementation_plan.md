# EC Dry-Out / Solvent-Consumption Wrapper — Implementation Plan

Source paper: Ruihe Li, Simon O'Kane, Monica Marinescu, Gregory J Offer,
**"Modelling Solvent Consumption from SEI Layer Growth in Lithium-Ion
Batteries"**, *J. Electrochem. Soc.* 2022 169 060516.
(`C:\Users\ds3420\Downloads\solvent_consumption.pdf`)

**Status: implemented and validated.** `ec_dryout_wrapper.py` in this folder
is a complete, working script — every mechanism it relies on
(`starting_solution=` continuation across changed parameters, the R_Li
state-rescale, pore-buffering compatibility) was empirically verified
against the live PyBaMM install in this repo, not just designed on paper.
`spike_li_rescale.py` is the standalone regression test for the trickiest
piece (R_Li); run it after any PyBaMM upgrade.
`ruihe_dryout_validation.py` reproduces (a reduced-cycle-count version of)
the paper's own comparison: no-dry-out baseline vs. 0%/6%/9% extra
electrolyte, plotted paper-Fig.-3/A-3-style.

**Important fix applied during validation**: continuing each batch via
`solution.last_state` (needed for the R_Li rescale, §4.2) means each
batch's own returned `sol` object only contains *that batch's* cycle
history, not the full run's (`sol.cycles`/`sol["..."].entries` reset per
batch — verified directly, see §4.5). This doesn't affect the wrapper's
own physics (every ledger update reads end-of-batch values, which stay
correctly cumulative regardless), but it silently breaks any caller trying
to read a full-run trajectory off the final `sol`. Fixed by having
`run_ec_dryout_degradation` accumulate a `trajectory` dict batch-by-batch
as it runs and return it as a 4th return value — see §4.5 and the current
script listing in §5.

> **OCR caveat**: the PDF's text layer garbled some of the subscript/
> superscript-heavy equations, particularly Eqs. 20–32 (the electrolyte
> mixing/concentration-ratio derivations). Section 2 below re-derives those
> relationships independently from the paper's *prose* description (which
> extracted cleanly) and cross-checks the result against the garbled
> equation text where possible — they agree, and the re-derived form is
> what's implemented. Eqs. 1–19 (SEI growth, volume bookkeeping) extracted
> cleanly and are transcribed directly.

---

## 1. Why a wrapper, not a submodel

The paper's own implementation (their §"Model implementation", Fig. 2)
deliberately avoids adding new PyBaMM state variables/ODEs for the
electrolyte-volume bookkeeping. Instead, every *n* cycles they pause the
simulation, compute three scalar ratios from the just-completed solution,
and restart with three cell **parameters** patched: electrode area, bulk
Li⁺ concentration, and bulk EC concentration. They justify this ("The
proposed update method should give similar results as adding new
differential equations... since electrolyte consumption is a slow,
long-term process") and confirm it empirically in their Appendix by
sweeping the update frequency and finding negligible sensitivity.

This maps directly onto the batch-restart loop this project already uses
throughout `si_gr_expansion_precursor/test_pore_buffering/` (e.g.
`run_full_degradation`'s batched ageing loop) — the wrapper just also
patches parameters and (for one variable, §4.2) rescales part of the
carried-over state at each batch boundary, instead of continuing unchanged.

**No PyBaMM core changes are needed**, and (contrary to what an earlier
pass at this plan assumed) **no manual model-building/`set_initial_
conditions_from` gymnastics are needed either** — see §3.4/§4.2:
`sim.solve(starting_solution=prev_solution)`, the same continuation
mechanism this project's other scripts already use, was confirmed (by
direct testing, not just reading the source) to tolerate a **different**
`ParameterValues` object on the new `Simulation` than the one that produced
`prev_solution`. That single fact is what makes the whole wrapper simple.

---

## 2. Physical model (re-derived, clean notation)

### 2.1 SEI growth and its EC dependence (Eqs. 1–11, high-confidence transcription)

SEI-forming reaction: 2 Li⁺ + 2 EC + 2 e⁻ → (CH₂OCO₂Li)₂ + C₂H₄ — i.e. a
**1:1:1 stoichiometric ratio of Li⁺ : EC : electrons** (2 of each, 1 SEI
product).

Diffusion-limited SEI current density (solvent transport through the SEI
layer is rate-limiting):

```
j_SEI = -F * D_EC * c_EC / L_SEI
```

where `c_EC` is the (bulk) EC concentration and `L_SEI` the SEI thickness
the solvent has to diffuse through. SEI thickness grows as

```
dL_SEI/dt = -j_SEI * M_SEI / (2 * a_n * F * rho_SEI)
```

and this directly drives the electrode's porosity down (a native PyBaMM
behaviour, §3.3):

```
d(eps_n)/dt = -a_n * dL_SEI/dt
```

**This is precisely PyBaMM's built-in `"SEI": "solvent-diffusion limited"`
submodel** — the paper isn't proposing a new SEI kinetics model, only
making its `c_EC` input dynamic instead of constant.

### 2.2 Electrolyte volume bookkeeping (Eqs. 12–21, high-confidence)

Two volumes tracked outside PyBaMM:
- `V_eJR` — electrolyte volume actually inside the jelly-roll (JR) pores.
- `V_eres` — electrolyte volume in a reservoir outside the JR but inside
  the cell package (extra electrolyte added at manufacture).

At BOL: `V_eJR,0 = A_cell,0 * (L_n*eps_n0 + L_sep*eps_s0 + L_pos*eps_p0)`
(fully-wetted pores), `V_eres,0 = R_eres * V_eJR,0` for a user-chosen
"percent extra electrolyte" `R_eres` (the paper sweeps 0%, 6%, 9%).

Over an interval where `dn_EC` moles of EC are consumed by the SEI
reaction:
- Electrolyte volume lost to consumption: `dV_EC = dn_EC * M_EC / rho_EC`.
- Pore volume lost to SEI solid formation: **read directly from PyBaMM's
  own porosity output** (§2.1's last equation is already implemented
  natively) rather than re-derived from `dn_EC` and a separate SEI
  density/molar-volume constant this repo doesn't parameterise anyway
  (§3.2).

The paper argues that consuming 1 mol EC shrinks electrolyte volume by
more than the SEI solid that forms grows the pore-filling solid phase, so
electrolyte volume normally shrinks *faster* than pore volume, and the JR
would go under-filled without a reservoir. **This turned out NOT to hold
for this repo's `si_gr_expansion` parameter set** — see §6's validation
run: pore volume shrinks faster than electrolyte volume for this specific
chemistry's SEI parameters, so no shortfall (hence no dry-out) appeared
within 160 cycles at 0% reservoir. The wrapper's bookkeeping is correct
regardless of which direction wins — it's a property of the underlying
SEI molar-volume/density parameterisation, not a bug (§6 discusses this
further and what to check if you want to see actual dry-out with this
model).

### 2.3 The three per-update-step ratios (re-derived from prose; Eqs. 26/30/32 in the paper)

Let `t` be the start and `t+dt` the end of one update interval (a batch of
*n* cycles).

**Step 1 — EC consumed this interval**, from the SEI capacity-loss output
(§3.2):

```
dn_EC = delta_Q_SEI[A.h] * 3600 / F      # mol, via 1:1 Li:EC stoichiometry
```

**Step 2 — pore & electrolyte volumes**:

```
V_pore(t+dt)  = A_cell(t) * (L_n * eps_n_avg(t+dt) + L_sep*eps_s0 + L_pos*eps_p0)
dV_EC         = dn_EC * M_EC / rho_EC
V_eJR(t+dt-)  = V_eJR(t) - dV_EC          # before any reservoir refill
shortfall     = max(V_pore(t+dt) - V_eJR(t+dt-), 0)
```

**Step 3 — reservoir refill (if any)**:

```
dV_add        = min(shortfall, V_eres(t))
V_eres(t+dt)  = V_eres(t) - dV_add
V_eJR(t+dt)   = V_eJR(t+dt-) + dV_add
```

**Step 4 — dry-out ratio** (paper Eq. 26/27 — shrinks electrode area):

```
R_dry = clip(V_eJR(t+dt) / V_pore(t+dt), 0, 1)
A_cell(t+dt) = R_dry * A_cell(t)          # multiplicative, compounds over time
```

**Step 5 — bulk EC concentration update** (paper Eq. 31/32 — a pure mass
balance in the *bulk*, since `c_EC` isn't spatially resolved in PyBaMM,
§3.2):

```
n_EC(t)       = c_EC(t) * V_eJR(t)
n_EC(t+dt)    = n_EC(t) - dn_EC + c_EC0 * dV_add   # reservoir liquid carries BOL c_EC
c_EC(t+dt)    = n_EC(t+dt) / V_eJR(t+dt)
```

**Step 6 — Li⁺ concentration ratio** (paper Eq. 28–30 — a uniform rescale
of the *existing spatial profile*, applied only where reservoir liquid
mixed in):

```
R_Li = [c_Li_avg(t+dt-) * V_eJR(t) + c_Li0 * dV_add] / [c_Li_avg(t+dt-) * V_eJR(t+dt)]
```

using the X-averaged electrolyte Li⁺ concentration at the end of the batch
(`c_Li_avg(t+dt-)`, i.e. before mixing) as the reference. **If `dV_add = 0`
(no reservoir, or reservoir already empty), `R_Li = 1` identically** — the
"plain" concentration effect of shrinking pore volume around a fixed
amount of Li⁺ is **already** captured natively by PyBaMM's own electrolyte
mass conservation PDE (coupled to `eps_n(x,t)` already), so `R_Li`
represents *only* the extra correction from reservoir mixing.

### 2.4 What actually feeds back into the next batch

| Ratio / quantity | Applied to |
|---|---|
| `A_cell(t+dt)` (via `R_dry`) | `"Number of electrodes connected in parallel to make a cell"` |
| `c_EC(t+dt)` | `"{Primary,Secondary}: Bulk solvent concentration [mol.m-3]"` |
| `R_Li` | Rescales the state carried over via `starting_solution=` (§4.2) — implemented and verified |
| everything else (particle concentrations, SEI thickness/porosity state, mechanics state, temperature, pore-buffering's own outputs, …) | passed through unchanged automatically by `starting_solution=` |

---

## 3. Paper quantities → PyBaMM API (this repo, verified)

### 3.1 SEI kinetics

- Option: `options["SEI"] = "solvent-diffusion limited"` — implements
  exactly `j_SEI = -F * D_sol * c_sol / L_SEI`
  (`src/pybamm/models/submodels/interface/sei/sei_growth.py:207-209`).
  This is a **different** option from `"reaction limited"`, which is what
  every other script in this repo (`si_gr_expansion`, pore-buffering, etc.)
  uses — the EC-dry-out wrapper is a genuinely separate model
  configuration. `si_gr_expansion.py` already defines both
  `Primary:`/`Secondary:` prefixed solvent-diffusion-limited SEI parameters
  (verified: `"Primary: SEI solvent diffusivity [m2.s-1]": 2.5e-22`,
  `"Primary: Bulk solvent concentration [mol.m-3]": 2636.0`, and the
  `Secondary:` equivalents), so it runs against this project's own
  composite Si/Gr chemistry with no extra parameterisation needed.
- **Not** `"ec reaction limited"` — that's a different (Yang et al. 2017)
  mechanism with its own coupled EC-surface-concentration ODE.
- Caveat: this repo's `SEIGrowth`/`SEIThickness` classes track **one
  lumped SEI thickness**, not the paper's separate inner/outer layers with
  a split coefficient β (Eqs. 5–6). Treat the whole lumped `L_SEI` as the
  paper's outer (diffusion-limiting) layer — implicitly β=0. Does not
  affect the solvent-consumption bookkeeping, which only cares about
  `dn_EC` (via capacity loss, §3.2), not the thickness split.

### 3.2 Parameters and outputs

| Quantity | PyBaMM name | Notes |
|---|---|---|
| Bulk EC concentration `c_EC` | `"{Primary,Secondary}: Bulk solvent concentration [mol.m-3]"` | Fixed scalar `Parameter` — updating between batches is a plain `param.update(...)`. |
| EC diffusivity through SEI `D_EC` | `"{Primary,Secondary}: SEI solvent diffusivity [m2.s-1]"` | Constant, not updated by the wrapper. |
| Cumulative SEI capacity loss (→ `dn_EC`) | `"Loss of capacity to negative primary SEI [A.h]"` + `"...secondary SEI [A.h]"` (summed, composite negative electrode) | Cumulative from t=0 of the *whole run* (since `starting_solution=` carries the underlying state forward) — take `end_value − start_value` of the current batch. |
| Negative electrode porosity | `"X-averaged negative electrode porosity"` | |
| Electrolyte Li⁺ concentration | `"X-averaged electrolyte concentration [mol.m-3]"` (bulk ratio calc) | |
| Initial electrolyte Li⁺ concentration | `"Initial concentration in electrolyte [mol.m-3]"` | Used as the reservoir liquid's (BOL) concentration in the `R_Li` mass balance. |
| Electrode area | `"Number of electrodes connected in parallel to make a cell"`, `"Electrode width [m]"`, `"Electrode height [m]"` | See §4.3 — the wrapper scales the parallel-electrode count, not width/height. |
| Faraday constant | `pybamm.constants.F` (evaluate to a float) | **Not** `param["Faraday constant [C.mol-1]"]` — that key is deprecated/removed in this PyBaMM version (confirmed by running the wrapper: `KeyError: "Accessing 'Faraday constant [C.mol-1]' from ParameterValues is deprecated. Use pybamm.constants.F instead."`). |

There is no separate "SEI density"/"SEI molar mass" parameter in this repo
(only `"SEI partial molar volume [m3.mol-1]"`) — not needed anyway since
pore volume is read from PyBaMM's own porosity output (§2.2).

### 3.3 Porosity coupling already native to PyBaMM

`"SEI porosity change": "true"` (already used throughout this project)
implements exactly the paper's Eq. 8 (`d(eps_n)/dt = -a_n * dL_SEI/dt`)
internally. The wrapper reads `V_pore(t)` from the solved porosity output
rather than re-deriving it from `dn_EC`.

### 3.4 Cross-batch state continuation — corrected finding

An earlier pass at this plan assumed `pybamm.BaseModel.set_initial_
conditions_from` plus a manually-built `Simulation` would be required,
based on a (accurate but incomplete) reading of `tests/unit/
test_simulation.py`. **Direct testing overturned the "identical
parameters only" assumption behind that**: `Simulation.solve(starting_
solution=prev_solution)` — the exact mechanism this project's other
scripts already use for continuing between ageing batches — works
correctly even when the new `Simulation` has a **different**
`ParameterValues` object (different bulk EC concentration, different
electrode area), including with an `Experiment` attached (multi-step
cycling), and including with `"pore buffering": "true"`. This was verified
directly, not inferred:

```python
sol2 = sim2.solve(starting_solution=sol1)   # sim2's ParameterValues != sim1's
```

ran cleanly and produced physically correct results (verified: bulk EC
concentration in the continued solve reflected `sim2`'s new value, not
`sim1`'s). This makes the wrapper dramatically simpler than the original
plan's `set_initial_conditions_from`/manual-`build()` approach, which (also
verified directly) hits an internal `symbol_processor`/state-mapper
inconsistency specifically when combined with `Experiment`-based solving —
that path is **not used** by the final implementation; §4.2 explains what
replaced it.

---

## 4. Design decisions (final, implemented)

### 4.1 Everything except R_Li: plain parameter patch + `starting_solution=`

`R_dry` (via electrode area) and `R_EC` (via bulk EC concentration) are
both scalar `Parameter` updates. Between batches:

```python
param.update({
    "Primary: Bulk solvent concentration [mol.m-3]": c_EC_new,
    "Secondary: Bulk solvent concentration [mol.m-3]": c_EC_new,
    "Number of electrodes connected in parallel to make a cell": n_parallel_new,
})
model = pybamm.lithium_ion.DFN(options)
sim = pybamm.Simulation(model, parameter_values=param, experiment=batch_experiment)
sol = sim.solve(starting_solution=prev_solution.last_state)
```

No manual `build()`, no `set_initial_conditions_from` — `starting_
solution=` handles everything (particle concentrations, SEI state,
porosity's underlying drivers, cracking, active-material volume fractions,
pore-buffering's algebraic outputs) automatically and correctly.

### 4.2 R_Li: validated, exact, in-place state rescale

`set_initial_conditions_from` was abandoned for this (§3.4's caveat) in
favour of a lower-level but fully-validated mechanism: mutate the raw
solved state vector directly, on the SAME `Solution.last_state` object
that gets passed to `starting_solution=`, before passing it in.

The raw ODE/DAE state PyBaMM actually tracks for the electrolyte is
**`"Porosity times concentration [mol.m-3]"`** (a `Concatenation` over
negative electrode/separator/positive electrode), *not* `"Electrolyte
concentration [mol.m-3]"` directly — a finite-volume-conservative
formulation. Its slice in the discretised state vector is found via
`solution.last_state.all_models[-1].y_slices`, and — critically, verified
empirically (§6, `spike_li_rescale.py`) — this state has **zero reference
offset**, so a pure multiplicative rescale of the raw slice by `R_Li`
exactly reproduces the same rescale on the physical concentration:

```python
def apply_li_rescale(last_state, R_Li):
    if R_Li == 1.0:
        return last_state
    lm = last_state.all_models[-1]
    sl = next(
        slices[0] for var, slices in lm.y_slices.items()
        if var.name == "Porosity times concentration [mol.m-3]"
    )
    last_state.y[sl] = last_state.y[sl] * R_Li
    return last_state
```

Verified end-to-end (`spike_li_rescale.py`, passing): after `apply_li_
rescale(last_state, 1.10)` and `sim2.solve(starting_solution=last_state)`,
the new solution's `"X-averaged electrolyte concentration [mol.m-3]"` at
t=0 is exactly `1.10×` the pre-rescale value (ratio `1.100000` to solver
precision), while porosity and every other state variable (checked:
particle concentration) pass through **completely unaffected** — confirming
the rescale is isolated to exactly the electrolyte Li⁺ concentration, with
no cross-talk into porosity or any other state.

**Run `spike_li_rescale.py` after any PyBaMM upgrade** — this mechanism
relies on an internal state-naming/`y_slices` convention that isn't part
of PyBaMM's public API contract, so it's worth a quick regression check if
the installed PyBaMM version changes.

### 4.3 Electrode-area parameter to scale

`"Number of electrodes connected in parallel to make a cell"`, not
`"Electrode width [m]"`/`"Electrode height [m]"`:
- No mesh-extent role under any model configuration (those two double as
  literal 2D/2+1D current-collector mesh domain extents in
  `src/pybamm/geometry/battery_geometry.py`; this project's models don't
  use that mode, but the parallel-electrode-count parameter is safe
  regardless of future model-option changes).
- Pure multiplicative term in `A_cc = L_y * L_z * n_parallel`
  (`geometric_parameters.py`), so shrinking it by `R_dry` produces an
  identical area reduction.
- Track the **cumulative product** of all `R_dry` factors applied so far
  (`n_parallel(t) = n_parallel,0 * prod(R_dry_i)`), matching the paper's
  Eq. 27 compounding behaviour. Implemented as `ECDryoutLedger.area_scale`.

### 4.4 Pore-buffering compatibility — why no special-casing was needed

The user's requirement was that this wrapper work with `"pore buffering":
"true"` enabled. It does, with **zero** pore-buffering-specific code in
the wrapper, for a structural reason confirmed by re-reading
`reaction_driven_porosity.py`: the pore-buffering submodel introduces
**no new PyBaMM state variables** (`pybamm.Variable`s) of its own. It only
recomputes the *outputs* `"Negative electrode thickness change [m]"`,
`"Negative electrode transfer ratio k"`, and the reported porosity `eps_k`
as **algebraic functions of the SAME underlying state** (SEI
thickness/active-material-volume-fraction state, structural porosity
`eps_struct = epsilon_init + delta_eps_k`) that this wrapper already
carries forward via `starting_solution=`. Neither electrode area nor bulk
EC concentration enter the pore-buffering partition math (`_transmitted_
fraction`, `dv_buffered`, `dv_thickness`) at all. So:
- `R_dry`'s area shrink and `R_EC`'s concentration update are completely
  orthogonal to pore buffering's own logic — nothing to reconcile.
- `R_Li`'s rescale touches only the electrolyte concentration state, which
  pore buffering never reads.
- Every batch's fresh `pybamm.lithium_ion.DFN(options)` with `options["pore
  buffering"] = "true"` just re-derives the buffered porosity/thickness
  from whatever SEI/mechanics state was carried over — identical in kind
  to how it already works within a single (non-restarted) run.

This was directly confirmed, not just argued: the integration test in §6
ran the full wrapper (area shrink + EC concentration update + Li⁺ rescale,
all three simultaneously) with `"pore buffering": "true"`, and `"Negative
electrode transfer ratio k"` read back correctly (`k = 0.7000...` at BOL,
matching `f_transmit_min`) at the start of the continued batch.

`MODEL_OPTIONS` in `ec_dryout_wrapper.py` has pore buffering **on** by
default (`"pore buffering transition": "physical"`); set `options["pore
buffering"] = "false"` to compare against an unbuffered run using the same
wrapper, same pattern as `test_pore_buffering/pore_buffering_degradation_
test.py`'s `data_off`/`data_on` comparison.

### 4.5 Full-run trajectories: don't read them off the final `sol`

Continuing each batch via `solution.last_state` (§4.2's rescale mechanism
needs a `.last_state`-reduced Solution to find `y_slices` on) has a side
effect worth calling out explicitly, because it's easy to get wrong
silently: **the returned `sol` after N batches only contains the LAST
batch's own cycle history**, not the full run's. Verified directly:

```python
sol1 = sim1.solve(...)                       # 3-cycle experiment -> len(sol1.cycles) == 9
sol2 = sim2.solve(starting_solution=sol1.last_state)   # another 3 cycles
len(sol2.cycles)                             # == 10, NOT ~18
sol2["Throughput capacity [A.h]"].entries[0] # == 29.26, i.e. starts where batch 1 left off,
                                              # not from t=0 -- batch 1's data points are gone
```

(For comparison: `sim3.solve(starting_solution=sol1)` — the FULL solution,
not `.last_state` — correctly gives `len(sol3.cycles) == 18`. This project's
other scripts, e.g. `test_pore_buffering/pore_buffering_degradation_test.py`,
all pass the full solution to `starting_solution=`, which is why they don't
hit this — they don't need the `.last_state`-specific R_Li rescale hook.)

This does **not** corrupt the wrapper's own physics: `ECDryoutLedger.update()`
only ever reads `.entries[-1]` (the end-of-batch value) from each batch's
own `sol`, and that value is still correctly the true cumulative physical
quantity (the underlying ODE/DAE *state* — SEI amount, particle
concentrations, etc. — is carried forward correctly regardless; only the
*trajectory/history* of a single Solution object is truncated). But it
would silently corrupt any caller trying to plot a full-run trajectory by
reading `final_sol.cycles` or `final_sol["..."].entries` after the loop —
they'd only see the last batch.

**Fixed** by having `run_ec_dryout_degradation` call `extract_batch_series`
on every batch's `sol` right after it's produced (while that batch's own
local history is still intact) and accumulate the results into a
`trajectory` dict, returned as a 4th value. Use `trajectory`, not
`final_solution.cycles`, for any full-run analysis — see
`ruihe_dryout_validation.py` (§8) for the pattern.

---

## 5. `ec_dryout_wrapper.py` (full listing)

This is the actual file at `si_gr_expansion_precursor/ec_dryout/
ec_dryout_wrapper.py` — implemented and validated (section 6), not a sketch.

```python
"""
ec_dryout_wrapper.py -- standalone electrolyte dry-out / solvent-consumption
wrapper around PyBaMM, implementing Li et al. 2022 (J. Electrochem. Soc.
169 060516) "Modelling Solvent Consumption from SEI Layer Growth in
Lithium-Ion Batteries" as a batch-restart post-processing loop. See
si_gr_expansion_precursor/ec_dryout/implementation_plan.md for the full
derivation, API research and validation history behind this script.

No PyBaMM core changes. Every UPDATE_EVERY_N_CYCLES cycles, the just-
completed batch's solution is used to compute three ratios (R_dry, R_EC,
R_Li) per the paper's Eqs. 26/30/32, applied as:
  - R_dry -> shrinks "Number of electrodes connected in parallel to make
    a cell" (a pure area multiplier, safe for any current-collector mode --
    see plan doc section 3.5) for the NEXT batch's ParameterValues.
  - R_EC  -> sets "{Primary,Secondary}: Bulk solvent concentration
    [mol.m-3]" directly (it is a plain scalar Parameter in PyBaMM's
    "solvent-diffusion limited" SEI submodel, not a PDE state).
  - R_Li  -> rescales the transplanted electrolyte Li+ concentration
    state ("Porosity times concentration [mol.m-3]") in place on
    `solution.last_state.y`, exact and validated (spike_li_rescale.py;
    the state has zero reference offset, so a pure multiplicative
    rescale is correct).
Everything else (particle concentrations, SEI amount, active-material
volume fractions, cracking, pore-buffering's own algebraic outputs, ...)
is carried over unmodified via `sim.solve(starting_solution=...)`, which
this project already uses elsewhere and which was confirmed (this
implementation) to tolerate a DIFFERENT ParameterValues object between
batches -- it is not restricted to identical-parameter continuation as
originally assumed.

Compatible with "pore buffering": "true" with NO special-casing needed:
the pore-buffering submodel (reaction_driven_porosity.py) introduces no
new state variables of its own -- it only recomputes porosity/thickness
outputs from the SAME underlying SEI-thickness state this wrapper already
transplants, and neither electrode area nor bulk EC concentration enter
its partition equations. See implementation_plan.md section 4.4.
"""
import numpy as np
import pybamm

# ---------------------------------------------------------------------------
# EC (ethylene carbonate) physical constants -- NOT PyBaMM Parameters (the
# repo has no "solvent density"/"solvent molar mass" parameter; these are
# fixed physical constants of the specific solvent species, kept as plain
# module constants, matching how the paper's own Table A-II treats them as
# external to the DFN parameter table).
# ---------------------------------------------------------------------------
M_EC = 88.062e-3     # EC molar mass [kg/mol]
RHO_EC = 1321.0      # EC density [kg/m3] (~1.32 g/cm3)

MODEL_OPTIONS = {
    "particle phases": ("2", "1"),
    "open-circuit potential": (("single", "current sigmoid"), "single"),
    "SEI": "solvent-diffusion limited",   # required: matches the paper's
                                           # j_SEI = -F*D_EC*c_EC/L_SEI exactly
    "SEI porosity change": "true",        # gives V_pore(t) via porosity --
                                           # no need to re-derive from dn_EC
    "particle mechanics": ("swelling and cracking", "swelling only"),
    "SEI on cracks": "true",
    "loss of active material": "stress-driven",
    "pore buffering": "true",
    "pore buffering transition": "physical",
}

AREA_PARAM = "Number of electrodes connected in parallel to make a cell"
POROSITY_TIMES_CONC_NAME = "Porosity times concentration [mol.m-3]"


# ---------------------------------------------------------------------------
# R_Li state rescale -- validated mechanism, see spike_li_rescale.py
# ---------------------------------------------------------------------------
def _find_y_slice(solution_last_state, var_name):
    lm = solution_last_state.all_models[-1]
    for var, slices in lm.y_slices.items():
        if getattr(var, "name", None) == var_name:
            return slices[0]
    raise KeyError(
        f"{var_name!r} not found in y_slices -- PyBaMM's internal state "
        f"representation may have changed; re-run spike_li_rescale.py"
    )


def apply_li_rescale(last_state, R_Li):
    """In-place rescale of the transplanted electrolyte Li+ concentration
    state by R_Li (a no-op if R_Li == 1, i.e. whenever no reservoir liquid
    mixed in this interval). Mutates and returns last_state."""
    if R_Li == 1.0:
        return last_state
    sl = _find_y_slice(last_state, POROSITY_TIMES_CONC_NAME)
    last_state.y[sl] = last_state.y[sl] * R_Li
    return last_state


# ---------------------------------------------------------------------------
# External ledger: everything the paper treats as independent of PyBaMM's
# own state (JR/reservoir electrolyte volumes, bulk EC concentration,
# cumulative electrode-area shrink factor). Plain Python floats, updated
# once per batch from that batch's solved output.
# ---------------------------------------------------------------------------
class ECDryoutLedger:
    def __init__(self, param, r_eres, area_param_name=AREA_PARAM, composite=True):
        self.param = param
        self.area_param_name = area_param_name
        self.composite = composite
        self.area_scale = 1.0            # cumulative product of R_dry
        self.n_parallel_0 = param[area_param_name]
        c_ec_key = (
            "Primary: Bulk solvent concentration [mol.m-3]" if composite
            else "Bulk solvent concentration [mol.m-3]"
        )
        self.c_EC = param[c_ec_key]
        self.c_EC0 = self.c_EC   # BOL EC concentration (reservoir liquid's
                                  # concentration, per the paper's assumption
                                  # that reservoir electrolyte matches the
                                  # initial JR concentration -- Eq. 31)
        self.r_eres = r_eres
        self.V_pore = None
        self.V_eJR = None
        self.V_eres = None

    def initialise_from_formation_solution(self, sol0):
        """Call once, right after the formation solve, before any ageing
        batches. Reads BOL geometry/porosity from the SOLVED state
        (robust to porosity being a FunctionParameter, unlike trying to
        hand-compute eps0 from raw Parameters)."""
        p = self.param
        L_n = p["Negative electrode thickness [m]"]
        L_s = p["Separator thickness [m]"]
        L_p = p["Positive electrode thickness [m]"]
        eps_n0 = sol0["X-averaged negative electrode porosity"].entries[0]
        eps_s0 = p["Separator porosity"]
        eps_p0 = p["Positive electrode porosity"]
        A_cell0 = (
            p["Electrode width [m]"] * p["Electrode height [m]"] * self.n_parallel_0
        )
        self.V_pore = A_cell0 * (L_n * eps_n0 + L_s * eps_s0 + L_p * eps_p0)
        self.V_eJR = self.V_pore          # fully wetted at BOL
        self.V_eres = self.r_eres * self.V_eJR

    def update(self, delta_Q_SEI_Ah, eps_n_avg_new, c_Li_avg_end, F):
        """Advance the ledger by one batch. Returns a dict of the new
        absolute values ready to feed into the next batch's ParameterValues
        / initial-condition rescale.

        delta_Q_SEI_Ah : this batch's (end - start) cumulative
            "Loss of capacity to negative {primary,secondary} SEI [A.h]"
            (summed over phases -- both phases share the SEI reaction).
        eps_n_avg_new : "X-averaged negative electrode porosity" at the
            end of this batch.
        c_Li_avg_end : "X-averaged electrolyte concentration [mol.m-3]"
            at the end of this batch (pre-mixing reference for R_Li).
        F : Faraday constant [C/mol].
        """
        p = self.param
        L_n = p["Negative electrode thickness [m]"]
        L_s = p["Separator thickness [m]"]
        L_p = p["Positive electrode thickness [m]"]
        eps_s0 = p["Separator porosity"]
        eps_p0 = p["Positive electrode porosity"]
        A_cell = (
            p["Electrode width [m]"] * p["Electrode height [m]"]
            * self.n_parallel_0 * self.area_scale
        )

        # Step 1: EC consumed this interval (1:1:1 Li:EC:e- stoichiometry,
        # Eq. 1: 2Li+ + 2EC + 2e- -> SEI + gas)
        dn_EC = delta_Q_SEI_Ah * 3600.0 / F

        # Step 2: pore & electrolyte volumes
        V_pore_new = A_cell * (L_n * eps_n_avg_new + L_s * eps_s0 + L_p * eps_p0)
        dV_EC = dn_EC * M_EC / RHO_EC
        V_eJR_predryout = self.V_eJR - dV_EC
        shortfall = max(V_pore_new - V_eJR_predryout, 0.0)

        # Step 3: reservoir refill, if any
        dV_add = min(shortfall, self.V_eres)
        V_eres_new = self.V_eres - dV_add
        V_eJR_new = V_eJR_predryout + dV_add

        # Step 4: dry-out ratio -> cumulative area shrink
        R_dry = float(np.clip(V_eJR_new / V_pore_new, 0.0, 1.0))
        self.area_scale *= R_dry

        # Step 5: bulk EC concentration update
        n_EC_old = self.c_EC * self.V_eJR
        n_EC_new = n_EC_old - dn_EC + self.c_EC0 * dV_add
        c_EC_new = max(n_EC_new / V_eJR_new, 0.0)

        # Step 6: Li+ ratio (only non-trivial when dV_add > 0)
        if dV_add > 0:
            c_Li0 = p["Initial concentration in electrolyte [mol.m-3]"]
            R_Li = (
                c_Li_avg_end * self.V_eJR + c_Li0 * dV_add
            ) / (c_Li_avg_end * V_eJR_new)
        else:
            R_Li = 1.0

        # commit
        self.V_pore, self.V_eJR, self.V_eres, self.c_EC = (
            V_pore_new, V_eJR_new, V_eres_new, c_EC_new,
        )

        return dict(
            R_dry=R_dry, R_Li=R_Li, c_EC_new=c_EC_new,
            n_parallel_new=self.n_parallel_0 * self.area_scale,
            V_pore=V_pore_new, V_eJR=V_eJR_new, V_eres=V_eres_new,
            dn_EC=dn_EC, dV_add=dV_add,
        )


# ---------------------------------------------------------------------------
# batch-restart driver
# ---------------------------------------------------------------------------
def sei_capacity_loss_Ah(sol, composite=True):
    """Cumulative "Loss of capacity to negative SEI [A.h]" -- summed over
    both phases if composite=True (both phases share the SEI reaction,
    Eq. 1, and therefore both consume EC), or the single-phase variant
    otherwise (e.g. OKane2022, a plain single-electrode chemistry)."""
    if composite:
        return (
            sol["Loss of capacity to negative primary SEI [A.h]"].entries[-1]
            + sol["Loss of capacity to negative secondary SEI [A.h]"].entries[-1]
        )
    return sol["Loss of capacity to negative SEI [A.h]"].entries[-1]


def _cycle_ageing_leg_capacity(cyc):
    """(cap, discharge_end_thr) for the ageing-rate (>1 A mean current)
    discharge leg of this cycle, or (0.0, None) if not found."""
    for step in getattr(cyc, "steps", []):
        try:
            I = step["Current [A]"].entries
            q = step["Discharge capacity [A.h]"].entries
            thr_step = step["Throughput capacity [A.h]"].entries
        except (KeyError, TypeError, AttributeError):
            continue
        if I.size >= 2 and np.mean(I) > 1.0 and q.size >= 2:
            return float(q[-1] - q[0]), float(thr_step[-1])
    return 0.0, None


def extract_batch_series(sol, composite=True):
    """Per-cycle discharge capacity + degradation-mode values, for the
    cycles contained in THIS SINGLE BATCH's `sol` alone.

    IMPORTANT: because run_ec_dryout_degradation continues each batch via
    `solution.last_state` (a single-timestep reduction, needed for the
    R_Li rescale -- see apply_li_rescale), each batch's own `sol` does
    NOT retain earlier batches' cycle history the way this project's other
    scripts' `starting_solution=full_solution` pattern does (verified:
    passing `.last_state` gives a `sol.cycles` covering only that batch's
    own cycles, with `sol["Throughput capacity [A.h]"]` starting from
    that batch's own t=0, not the run's). The underlying STATE is still
    correctly carried forward (confirmed: end-of-batch cumulative values
    like SEI capacity loss are physically correct), so this only matters
    for reconstructing a FULL-RUN trajectory -- which is exactly why
    run_ec_dryout_degradation calls this function once per batch and
    accumulates the results into its returned `trajectory` dict, rather
    than callers trying to read a full history off the final returned
    `sol` (which would silently only reflect the last batch)."""
    Qt_full = sol["Throughput capacity [A.h]"].entries
    if composite:
        LLI_full = (
            sol["Loss of capacity to negative primary SEI [A.h]"].entries
            + sol["Loss of capacity to negative secondary SEI [A.h]"].entries
        )
    else:
        LLI_full = sol["Loss of capacity to negative SEI [A.h]"].entries
    # LAM variables only exist if "loss of active material" is enabled --
    # minimal test configs (e.g. the OKane2022 dry-out-only setup) turn
    # this off entirely, so degrade gracefully to all-zero rather than
    # raising a KeyError.
    try:
        LAM_neg_full = sol["Loss of active material in negative electrode [%]"].entries
        LAM_pos_full = sol["Loss of active material in positive electrode [%]"].entries
    except KeyError:
        LAM_neg_full = np.zeros_like(Qt_full)
        LAM_pos_full = np.zeros_like(Qt_full)

    thr_list, cap_list = [], []
    for cyc in sol.cycles:
        cap, thr = _cycle_ageing_leg_capacity(cyc)
        if cap > 0 and thr is not None:
            thr_list.append(thr)
            cap_list.append(cap)
    thr_arr = np.array(thr_list)
    cap_arr = np.array(cap_list)
    if thr_arr.size:
        LLI_arr = np.interp(thr_arr, Qt_full, LLI_full)
        LAM_neg_arr = np.interp(thr_arr, Qt_full, LAM_neg_full)
        LAM_pos_arr = np.interp(thr_arr, Qt_full, LAM_pos_full)
    else:
        LLI_arr = LAM_neg_arr = LAM_pos_arr = np.array([])
    return dict(
        thr=thr_arr, cap=cap_arr, LLI=LLI_arr,
        LAM_neg=LAM_neg_arr, LAM_pos=LAM_pos_arr,
    )


def run_ec_dryout_degradation(
    base_param,
    formation_experiment,
    ageing_cycle_steps,
    update_every_n_cycles,
    max_total_cycles,
    r_eres,
    options=None,
    apply_li_correction=True,
    composite=True,
):
    """Top-level batch-restart loop.

    base_param : pybamm.ParameterValues, NOT yet containing the dry-out
        wrapper's per-batch overrides (those get patched in fresh each
        batch as the loop runs).
    formation_experiment : pybamm.Experiment for the initial (BOL) solve.
    ageing_cycle_steps : the per-cycle experiment step tuple, e.g.
        ("Discharge at 1C until 2.5 V", "Charge at 0.3C until 4.2 V",
         "Hold at 4.2 V until C/100").
    update_every_n_cycles : batch size -- how many ageing cycles run
        between each EC-dryout ledger update (paper's own protocol uses
        ~16; this project's other scripts default to 20-50).
    r_eres : ratio of initial extra (reservoir) electrolyte to initial
        JR electrolyte volume (paper's Eq. 33; sweep 0.0, 0.06, 0.09 to
        reproduce the paper's three cases).
    options : model options dict; defaults to MODEL_OPTIONS (pore
        buffering ON, "physical" transition) if not given.
    composite : True for a composite (Primary:/Secondary:-prefixed)
        negative electrode chemistry (e.g. si_gr_expansion); False for a
        plain single-phase chemistry (e.g. OKane2022) -- controls which
        parameter/variable names are used for bulk solvent concentration
        and SEI capacity loss.

    Returns (final_solution, history, ledger, trajectory).
    `history` is a list of the dict returned by ECDryoutLedger.update()
    for each batch, in order. `trajectory` is a dict of numpy arrays
    (keys: thr, cap, LLI, LAM_neg, LAM_pos) giving the FULL-RUN, cycle-by-
    cycle degradation trajectory, accumulated batch-by-batch as the loop
    runs -- use this for full-run plotting, NOT `final_solution.cycles`
    (which only covers the last batch; see extract_batch_series's
    docstring for why).
    """
    options = dict(options or MODEL_OPTIONS)
    param = base_param.copy()
    F = float(pybamm.constants.F.evaluate())

    sim0 = pybamm.Simulation(
        pybamm.lithium_ion.DFN(options), parameter_values=param,
        experiment=formation_experiment,
    )
    sol = sim0.solve(initial_soc=1.0)
    print("[ec_dryout] formation cycle solved.")

    ledger = ECDryoutLedger(param, r_eres, composite=composite)
    ledger.initialise_from_formation_solution(sol)
    print(
        f"[ec_dryout] BOL: V_pore={ledger.V_pore:.4e} m3  V_eJR={ledger.V_eJR:.4e} m3  "
        f"V_eres={ledger.V_eres:.4e} m3 (r_eres={r_eres:.2%})"
    )

    Q_SEI_prev = sei_capacity_loss_Ah(sol, composite=composite)
    history = []
    traj_parts = dict(thr=[], cap=[], LLI=[], LAM_neg=[], LAM_pos=[])
    total_cycles = 0
    while total_cycles < max_total_cycles:
        batch_exp = pybamm.Experiment(
            [ageing_cycle_steps for _ in range(update_every_n_cycles)]
        )
        model = pybamm.lithium_ion.DFN(options)
        sim = pybamm.Simulation(model, parameter_values=param, experiment=batch_exp)

        last_state = sol.last_state
        if apply_li_correction and history:
            apply_li_rescale(last_state, history[-1]["R_Li"])

        sol = sim.solve(starting_solution=last_state)
        total_cycles += update_every_n_cycles

        batch_series = extract_batch_series(sol, composite=composite)
        for key in traj_parts:
            traj_parts[key].append(batch_series[key])

        Q_SEI_new = sei_capacity_loss_Ah(sol, composite=composite)
        delta_Q_SEI = Q_SEI_new - Q_SEI_prev
        Q_SEI_prev = Q_SEI_new

        eps_n_avg_new = sol["X-averaged negative electrode porosity"].entries[-1]
        c_Li_avg_end = sol["X-averaged electrolyte concentration [mol.m-3]"].entries[-1]
        step = ledger.update(delta_Q_SEI, eps_n_avg_new, c_Li_avg_end, F)
        history.append(step)

        if composite:
            ec_updates = {
                "Primary: Bulk solvent concentration [mol.m-3]": step["c_EC_new"],
                "Secondary: Bulk solvent concentration [mol.m-3]": step["c_EC_new"],
            }
        else:
            ec_updates = {"Bulk solvent concentration [mol.m-3]": step["c_EC_new"]}
        param.update(
            {**ec_updates, AREA_PARAM: step["n_parallel_new"]},
            check_already_exists=False,
        )

        print(
            f"[ec_dryout] batch {len(history)}: {total_cycles} cycles  "
            f"R_dry={step['R_dry']:.4f}  R_Li={step['R_Li']:.4f}  "
            f"c_EC={step['c_EC_new']:.1f} mol/m3  V_eres={step['V_eres']:.3e} m3  "
            f"dn_EC={step['dn_EC']:.3e} mol"
        )

    trajectory = {
        key: (np.concatenate(parts) if parts else np.array([]))
        for key, parts in traj_parts.items()
    }
    return sol, history, ledger, trajectory


if __name__ == "__main__":
    param = pybamm.ParameterValues("si_gr_expansion")
    formation_exp = pybamm.Experiment(
        [
            "Discharge at 0.1C until 2.5 V",
            "Charge at C/3 until 4.2 V",
            "Hold at 4.2 V until C/100",
        ]
    )
    ageing_cycle = (
        "Discharge at 1C until 2.5 V",
        "Charge at 0.3C until 4.2 V",
        "Hold at 4.2 V until C/100",
    )
    sol, history, ledger, trajectory = run_ec_dryout_degradation(
        param, formation_exp, ageing_cycle,
        update_every_n_cycles=20, max_total_cycles=60,
        r_eres=0.06,   # 6% extra electrolyte, matching one of the paper's cases
    )
    print(f"\nFinal: {len(history)} batches, {sum(h['dn_EC'] for h in history):.3e} "
          f"mol EC consumed total, final c_EC={history[-1]['c_EC_new']:.1f} mol/m3, "
          f"cumulative area scale={ledger.area_scale:.4f}")
    print(f"Full-run trajectory: {trajectory['thr'].size} cycles, "
          f"capacity {trajectory['cap'][0]:.4f} -> {trajectory['cap'][-1]:.4f} A.h")
```

## 6. Validation performed

Run against the live PyBaMM install in this repo (`PyBaMM_Pressure` conda
env — see the memory note on this; the base anaconda `python` has a
separately-installed pybamm missing all local changes and will fail with
`'si_gr_expansion' is not a valid parameter set`).

1. **`spike_li_rescale.py`** (standalone regression test, passing):
   confirms `starting_solution=` tolerates a simultaneously-changed bulk
   EC concentration *and* electrode area, confirms the R_Li rescale is
   exact (`ratio: 1.100000`, expected `1.1`), and confirms porosity is
   completely unaffected by the rescale.
2. **Full wrapper smoke test** (`ec_dryout_wrapper.py`'s own `__main__`,
   3 batches × 20 cycles, `r_eres=0.06`, pore buffering ON): ran cleanly
   end-to-end. `c_EC` decreased monotonically each batch (2636.0 → 2616.5
   → 2601.2 → 2588.3 mol/m³) — the paper's headline negative-feedback
   direction. `R_dry`/`R_Li` stayed at 1.0000 throughout, i.e. no dry-out
   yet at only 60 cycles with a 6% reservoir — consistent with the paper's
   own ~390-cycle deviation point for a comparable (6%) case; this run
   simply didn't go far enough to reach it.
3. **0%-reservoir edge case, extended to 160 cycles**: `R_dry` stayed
   exactly `1.0000000000` the entire run — i.e. **no shortfall ever
   appeared** for this chemistry (`V_eJR` stayed above `V_pore` the whole
   time: e.g. at cycle 160, `V_eJR=5.3268e-6 m3` vs `V_pore=5.0232e-6
   m3`). This means, for `si_gr_expansion`'s specific SEI/porosity
   parameterisation, **pore volume shrinks faster than electrolyte volume
   consumed** — the opposite of the paper's general Eq. 19 argument (which
   was derived for their own SEI density/molar-volume assumptions, not
   this repo's). The wrapper's bookkeeping is verified correct (§6.1-6.2);
   whether *this specific chemistry* exhibits dry-out is a separate,
   legitimate physics question. **If reproducing the paper's dry-out
   phenomenon specifically is the goal**, the next step is checking this
   repo's implicit SEI molar-volume/density (via `"SEI partial molar
   volume [m3.mol-1]"` and the porosity-change relation in
   `reaction_driven_porosity.py`) against the paper's own Table A-I/A-II
   values, or simply running far more cycles — not a wrapper bug to fix.
4. **Pore-buffering compatibility**: the integrated test (area shrink + EC
   concentration + Li⁺ rescale, all three together, `"pore buffering":
   "true"`) solved cleanly across a batch restart, and `"Negative electrode
   transfer ratio k"` read back correctly (`k=0.7000` at the start of the
   continued batch, matching `f_transmit_min`'s BoL value) — confirming
   the pore-buffering submodel's own state/outputs survive the
   `starting_solution=` continuation with updated parameters intact.

5. **First `ruihe_dryout_validation.py` run, si_gr_expansion chemistry**
   (baseline + `r_eres` ∈ {0%, 6%, 9%}, `UPDATE_EVERY_N_CYCLES=20`,
   `MAX_TOTAL_CYCLES=160`): completed cleanly. Result: **all four
   conditions tracked each other almost exactly** — capacity retention
   94.19% (baseline) vs 94.20% (all three dry-out conditions, identical to
   2 d.p.), because `R_dry`/`R_Li` stayed at exactly `1.0000` throughout
   every `r_eres` case (consistent with point 3: this chemistry never
   triggers a shortfall in any practical cycle budget). The tiny
   baseline-vs-dry-out gap **was** the paper's secondary negative-feedback
   mechanism showing through (bulk EC concentration monotonically
   decreasing, 2616.5 → 2541.1 mol/m³), just without the headline
   `r_eres`-dependent dry-out effect itself.

6. **Superseded by an OKane2022-based version** (§6a below) at the user's
   request, specifically to demonstrate `R_dry` actually dropping below 1
   — `ruihe_dryout_validation.py` now contains the OKane2022 version, same
   filename, same output plot filenames (overwritten). The si_gr_expansion
   run above is kept here as a documented data point (it's genuinely
   useful: it's what first established the "ratio issue" finding in point
   3), but is no longer reproducible from a script in this folder as-is —
   see git history / this section if you need to reconstruct it.

### 6a. OKane2022-based dry-out demonstration (current `ruihe_dryout_validation.py`)

Switched to `pybamm.ParameterValues("OKane2022")` (a plain, single/
non-composite negative electrode) with `MODEL_OPTIONS = {"SEI":
"solvent-diffusion limited", "SEI porosity change": "true", "pore
buffering": "false"}` — i.e. everything except SEI growth and its
porosity coupling turned off, per the user's request, specifically to
isolate and demonstrate the dry-out mechanism itself (as opposed to the
earlier si_gr_expansion run, which already separately validated pore-
buffering compatibility, §4.4, but never showed genuine dry-out).

**Calibration finding (important, not just "increase the rate")**: the
natural first attempt — boosting `"SEI solvent diffusivity [m2.s-1]"`
alone — was tried up to a 3000× multiplier over a full 200 cycles, and
`R_dry` stayed at *exactly* `1.0000` the entire time regardless. This
turned out to be structural, not a rate problem: `j_SEI ∝ c_EC`, so both
the EC-consumption rate (shrinking `V_eJR` via the wrapper's external
`M_EC`/`RHO_EC`) and PyBaMM's own porosity-decrease rate (shrinking
`V_pore` via `"SEI partial molar volume [m3.mol-1]"`) are driven by the
same `j_SEI` and decay together self-limitingly as `c_EC` drops — boosting
the overall rate just makes both decay faster in lockstep, never changing
their *ratio*. That ratio is set by OKane2022's `"SEI partial molar
volume"` (9.585e-05 m³/mol) being large relative to the wrapper's EC molar
volume (`M_EC/RHO_EC` = 6.67e-05 m³/mol) — this parameterisation implies
*more* solid volume forms per mole of EC reacted than EC's own molecular
volume, so `V_pore` structurally shrinks at least as fast as `V_eJR` for
*any* rate.

**Fix**: reduce `"SEI partial molar volume [m3.mol-1]"` (denser/thinner
simulated SEI per mole reacted) *alongside* boosting the diffusivity.
Calibrated by direct testing:

| Diffusivity × | Molar volume ÷ | Result (200 cycles, r_eres=0%) |
|---|---|---|
| 3000× | 1× (unchanged) | `R_dry` stayed exactly `1.0000` the whole run |
| 1000× | 10× | `R_dry` hit `0.9075` within just 10 cycles, 20% capacity fade by cycle 40 — too aggressive |
| **100×** | **5×** | **chosen** — `R_dry` settles into a sustained ~0.973→0.996 range across the full 200 cycles, 86.3% capacity retention at cycle 200 (`r_eres=0%`) — gradual, clearly visible, not catastrophic |

This is a transparent, documented calibration choice for *demonstrating*
the mechanism (both multipliers are declared as named constants at the
top of the script with this rationale in the docstring) — not a claim
about OKane2022's real physical SEI density.

**Full 4-condition run result** (baseline + `r_eres` ∈ {0%, 6%, 9%},
`UPDATE_EVERY_N_CYCLES=20`, `MAX_TOTAL_CYCLES=200`): see
`ruihe_dryout_validation_fig3.png` (capacity retention / LLI, cf. paper
Fig. 3(a-b) — LAM omitted, disabled in `MODEL_OPTIONS`) and
`ruihe_dryout_validation_ratios.png` (R_dry/R_Li/bulk EC concentration vs
cycle, cf. paper Fig. A-3(d-f)/A-6) in this folder.

---

## 7. Folder layout (current)

```
si_gr_expansion_precursor/ec_dryout/
    implementation_plan.md               <- this file
    ec_dryout_wrapper.py                  <- implemented, validated (section 6)
    spike_li_rescale.py                   <- passing regression test for the R_Li mechanism;
                                              re-run after any PyBaMM upgrade
    ruihe_dryout_validation.py            <- baseline vs. 0%/6%/9% r_eres comparison,
                                              reduced-cycle-count replica of paper Fig. 3/A-3
    ruihe_dryout_validation_fig3.png      <- output of the above (cf. paper Fig. 3(a-c))
    ruihe_dryout_validation_ratios.png    <- output of the above (cf. paper Fig. A-3(d-f)/A-6)
```

Future work (not needed for the wrapper itself): re-run
`ruihe_dryout_validation.py` with a much larger `MAX_TOTAL_CYCLES` (or a
smaller `NEG_POROSITY_FLOOR`-style SEI-density adjustment, see point 3) to
actually reach the dry-out regime and see the `r_eres`-dependent capacity
divergence the paper reports.
