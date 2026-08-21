# Pore-Buffering (Volume Partition) Model — Implementation Plan

Source: `Si_Gr_Expansion_Precursor/Notes_Papers/si-gr_expansion_precursor_notes.PDF`
("Making the Pore-Buffering Model Physical", 2026-07-31).

**Scope of this plan**: implement the volume-partition mechanism (PDF §5, Eqs. 17–19,
the three-regime transmitted fraction, and the read-off diagnostic `k`), i.e. Changes
1, 2, 3, 4, 6 from PDF §6. **Explicitly deferred**: PDF §3 (particle-stress pressure BC),
§4 (Yin et al. stack-pressure / `K_stack`, `C_pore(ε)`, `P(t)` module constraint), and
§6.7/Change 5 (feeding `Stack pressure [Pa]` into `σr(R)`). Those require a stack-mechanics
submodel and a pressure state that don't exist yet — separate follow-on phase.

**Naming note (per your correction):** the closure floor used by the buffering partition
(`ε_min` in Eq. 18/19) is a **new, separate parameter**, not the existing
`"{Domain} electrode porosity floor"` I added earlier this session (which is a numerical/
SEI-growth stability floor for the softplus in `reaction_driven_porosity.py`). Physically
they represent different things — one is "porosity low enough that the solver needs
protecting", the other is "porosity low enough that the pore network percolation-closes
and stops absorbing swelling." They may end up close numerically after fitting, but they
are fit and reasoned about independently. Below I call it `eps_min_transfer`
(Parameter string: `"{Domain} electrode pore buffering closure porosity"`).

---

## 1) Scripts/folders to access

**Core submodels (the physics):**
- [src/pybamm/models/submodels/porosity/reaction_driven_porosity.py](src/pybamm/models/submodels/porosity/reaction_driven_porosity.py) — `ReactionDriven.get_coupled_variables`. Currently computes `eps_k = epsilon_init + delta_eps_k` (irreversible reaction-film loss only) then softplus-floors it at the existing `epsilon_min` (numerical floor). This is exactly the PDF's `ε_struct` (Eq. 20) *before* it gets a name.
- [src/pybamm/models/submodels/porosity/base_porosity.py](src/pybamm/models/submodels/porosity/base_porosity.py) — `_get_standard_porosity_variables`, builds the `"Porosity"` concatenation and per-domain `"{Domain} porosity"`/`"X-averaged {domain} porosity"` variables from an `eps_dict`. Whatever final `ε_total` value we compute must flow through this same helper so the global `"Porosity"` variable stays consistent.
- [src/pybamm/models/submodels/particle_mechanics/base_mechanics.py](src/pybamm/models/submodels/particle_mechanics/base_mechanics.py) — `_extract_mechanical_variables` (computes per-phase `v_change` = the PDF's per-phase `dv_solid` contribution, before scaling by `n_electrodes_parallel * L`), `_aggregate_phase_thickness_changes` (sums primary+secondary → `"{Domain} electrode thickness change [m]"`), `_aggregate_cell_thickness_change` (sums neg+pos+thermal → `"Cell thickness change [m]"`). This is PDF §6.1's existing `d_t` machinery — it currently routes 100% of swelling to thickness (Gap 1 in PDF §6.2).
- [src/pybamm/models/submodels/particle_mechanics/swelling_only.py](src/pybamm/models/submodels/particle_mechanics/swelling_only.py) and [crack_propagation.py](src/pybamm/models/submodels/particle_mechanics/crack_propagation.py) — both call `self._get_mechanical_results(variables)` from `get_coupled_variables`; this is where per-phase thickness change gets computed and injected into `variables`. `crack_propagation.py` additionally reads `stress_t_surf` for the crack-growth law — relevant later for the pressure phase, not this one.

**Model assembly (build order — this is the part that needs care):**
- [src/pybamm/models/full_battery_models/lithium_ion/base_lithium_ion_model.py](src/pybamm/models/full_battery_models/lithium_ion/base_lithium_ion_model.py) — `set_submodels()` (lines 43–64) calls, in order: `set_porosity_submodel()` (line 45) → `set_interface_utilisation_submodel()` (46) → `set_crack_submodel()` (47, registers `{domain} {phase}particle mechanics`) → `set_active_material_submodel()` (48) → `set_transport_efficiency_submodels()` (49) → `set_convection_submodel()` (50) → ... → electrolyte submodels (55–56). PyBaMM calls each submodel's `get_coupled_variables` in this same insertion order, accumulating into one shared `variables` dict. **Porosity is currently built before particle mechanics**, so today `ReactionDriven` cannot see per-phase swelling — it doesn't need to, because there's no buffering yet. `set_porosity_submodel()` itself lives at line 461, branching `Constant` vs `ReactionDriven` on `options["SEI porosity change"]`.

**Parameters:**
- [src/pybamm/parameters/lithium_ion_parameters.py](src/pybamm/parameters/lithium_ion_parameters.py) — `DomainLithiumIonParameters._set_parameters()`. `epsilon_init` (~line 285) and `epsilon_min` (line 293, `pybamm.Parameter(f"{Domain} electrode porosity floor")`, added this session) live here. `PhaseParameters.t_change` (line 897) and the existing `Omega`, `E` accessors live in the same file (phase-level, not domain-level) — these already exist and don't need touching.
- [src/pybamm/input/parameters/lithium_ion/si_gr_expansion.py](src/pybamm/input/parameters/lithium_ion/si_gr_expansion.py) — the parameter set to extend with the new buffering parameters (`eps_min_transfer`, `eps_max`, `f0`) as defaults. Already has per-phase `t_change` (graphite primary from Ai2020, silicon secondary from the CSV) and `epsilon_min = 0.08` default (inherited from the Mayur2024-derived bulk edit) — this stays untouched; the new params are additive.
- `pyproject.toml` / installed `entry_points.txt` — no change needed, `si_gr_expansion` is already registered.

**Empirical validation (this project's existing pipeline, not generic pytest):**
- [si_gr_expansion_precursor/test_expansion/](si_gr_expansion_precursor/test_expansion/) — `winner_deepdive_si20pct_expansion.py` / `..._1000cyc.py`. These already compute `"Cell thickness change [m]"`, the within-cycle expansion amplitude, and interpolate onto a per-cycle throughput grid — this is the harness to extend for `k(EFC)` once buffering exists, and it's also the source of the *unbuffered baseline* run needed for parameterisation (§3 below).
- Note: `tests/` in this checkout is minimal (no existing unit tests reference `ReactionDriven` or the mechanics classes by name — the checked-out test tree only has `conftest.py`/`shared.py`/`test_scripts.py`). Formal unit-test coverage for the new submodel should follow the pattern of neighboring PyBaMM submodel tests when this lands upstream, but for this project the primary validation loop is the deep-dive scripts above, not pytest.

---

## 2) Modifications

### 2.1 New parameters (`lithium_ion_parameters.py`, `DomainLithiumIonParameters._set_parameters()`)

Add three new domain-level parameters, right after the existing `epsilon_min` line:

```python
self.epsilon_min = pybamm.Parameter(f"{Domain} electrode porosity floor")          # existing (numerical/SEI stability)
self.eps_min_transfer = pybamm.Parameter(f"{Domain} electrode pore buffering closure porosity")  # NEW — mechanical pore-closure floor
self.eps_max_transfer = pybamm.Parameter(f"{Domain} electrode pore buffering upper porosity")    # NEW — upper "knee" porosity
self.f_transmit_min = pybamm.Parameter(f"{Domain} electrode transmitted fraction plateau")        # NEW — f0, dimensionless (0, 1]
```

Default values go into `si_gr_expansion.py` (and, if you want the change available to
every parameter set the way the numerical floor is, the same bulk-edit treatment applied
to the other 14 lithium-ion sets earlier this session — but I'd hold off on that until
the fit in §3 gives real numbers, otherwise the "defaults" are just placeholders).

### 2.2 `reaction_driven_porosity.py` — expose the structural porosity separately

Split the single `eps_k` computation into an irreversible-only `eps_struct` (exposed as a
new variable) and keep the existing softplus-floored value as a fallback for domains/phases
where buffering isn't active (e.g. positive electrode, or when the option is off):

```python
domain_param = self.param.domain_params[domain.split()[0]]
eps_struct = domain_param.epsilon_init + delta_eps_k   # irreversible only: SEI, plating, dead Li, cracks

if domain != "separator":
    variables[f"{Domain} electrode structural porosity"] = eps_struct
    eps_k = pybamm.softplus(eps_struct, domain_param.epsilon_min, k_eps)  # unchanged numerical floor
else:
    eps_k = eps_struct

eps_dict[domain] = eps_k
```

This is a pure addition (new variable, existing behavior unchanged) — safe to land on its
own, independent of the rest, and it's also the first thing needed for the parameterisation
step in §3 (you can extract `ε_struct(EFC)` from a baseline run before the partition exists).

### 2.3 Build-order change (`base_lithium_ion_model.py`)

The partition needs per-phase `dv_solid` (from particle mechanics, built at step 47)
*before* it can produce the final porosity (built at step 45) that transport submodels
(steps 49+) consume. Reorder:

```python
def set_submodels(self, build):
    self.set_external_circuit_submodel()
    self.set_interface_utilisation_submodel()
    self.set_crack_submodel()              # moved up: particle mechanics now runs before porosity
    self.set_porosity_submodel()           # moved down: can now read per-phase thickness-change vars
    self.set_active_material_submodel()
    self.set_transport_efficiency_submodels()
    ...
```

This reorder is safe to make unconditionally: `set_interface_utilisation_submodel()` and
`set_crack_submodel()` don't consume porosity, and everything that *does* consume porosity
(`set_transport_efficiency_submodels()`, `set_convection_submodel()`, electrolyte
submodels) still runs strictly after the new porosity position. Worth a quick grep for any
other consumer of `"{Domain} electrode porosity"` between the old and new position before
committing to this — I didn't find one, but the check is cheap.

### 2.4 `ReactionDriven.get_coupled_variables` — the partition itself (Eqs. 17–19)

Gated behind a new model option so existing parameter sets/behavior are untouched by
default (mirrors how `options["SEI porosity change"]` already branches `Constant` vs
`ReactionDriven` in `set_porosity_submodel`). Proposed option: `options["pore buffering"]`
= `"false"` (default) / `"true"`.

When `"true"`, after computing `eps_struct` (§2.2) and before the softplus floor:

```python
if self.options["pore buffering"] == "true" and domain != "separator":
    # 1. total solid swelling this step, summed over phases (shared pore space).
    #    Recovered from the thickness-change variables particle mechanics already
    #    computed (built earlier now, per the §2.3 reorder) rather than recomputing
    #    t_change() here — divide back out the n_electrodes_parallel * L scaling.
    thickness_change_unbuffered = variables[f"{Domain} electrode thickness change [m]"]
    dv_solid = thickness_change_unbuffered / (
        self.param.n_electrodes_parallel * domain_param.L
    )

    # 2. mechanical preference f(eps_struct): three-regime transmitted fraction (Eq. 18)
    f = self._transmitted_fraction(eps_struct, domain_param)

    # 3. hard geometric ceiling: can't buffer more void than exists above closure (Eq. 19)
    headroom = eps_struct - domain_param.eps_min_transfer
    dv_buffered = pybamm.minimum((1 - f) * dv_solid, headroom)  # binds only on net swelling
    dv_thickness = dv_solid - dv_buffered

    # 4. corrected thickness change (overwrites the unbuffered value from mechanics)
    variables[f"{Domain} electrode thickness change [m]"] = (
        self.param.n_electrodes_parallel * dv_thickness * domain_param.L
    )
    variables[f"{Domain} electrode buffered volume change"] = dv_buffered
    variables[f"{Domain} electrode solid volume change"] = dv_solid

    # 5. total porosity seen by transport = structural minus reversible buffering
    eps_struct = pybamm.maximum(eps_struct - dv_buffered, domain_param.eps_min_transfer)
    # (falls through to the existing softplus-floor line using the numerical eps_min)
```

with the three-regime fraction as its own small method:

```python
def _transmitted_fraction(self, eps_struct, domain_param):
    eps_min = domain_param.eps_min_transfer
    eps_max = domain_param.eps_max_transfer
    f0 = domain_param.f_transmit_min
    f_shared = 1 - (1 - f0) * (eps_struct - eps_min) / (eps_max - eps_min)
    return pybamm.minimum(pybamm.maximum(f_shared, f0), 1)
```

Use `pybamm.smooth_min`/a tanh-blend at the two piecewise corners if the IDA solver
struggles with the kinks (same pattern as the `exponent_max_sei` soft-cap already in
`sei_growth.py`, and the porosity floor's softplus) — flag this as a likely-needed
follow-up rather than building it in from the start, since it adds tuning surface
(smoothing sharpness) you don't want before the base partition is validated.

**Where `"Cell thickness change [m]"` recomputes**: `_aggregate_cell_thickness_change`
in `base_mechanics.py` sums `"Negative/Positive electrode thickness change [m]"` — since
those variables are now the *corrected* (post-buffering) values by the time anything
downstream reads `"Cell thickness change [m]"`, no change is needed there, **provided**
the negative electrode's mechanics submodel and the porosity submodel both write into the
same shared `variables` dict before that aggregation is read — true once the §2.3 reorder
lands, since mechanics (47) → porosity (45→ new position after 47) → nothing recomputes
`"Cell thickness change [m]"` again until it's actually consumed downstream (thermal,
diagnostics). Double-check this by grepping for read sites of `"{Domain} electrode
thickness change [m]"` between the new porosity position and wherever cell-level
aggregation is finalized, to make sure nothing reads the stale unbuffered value in between.

### 2.5 `k` as a read-off diagnostic (Eq. 23) — pure addition, wire into nothing

Small addition at the end of the same `get_coupled_variables`, or in a tiny new method:

```python
dv_th = pybamm.x_average(dv_thickness)
dv_sol = pybamm.x_average(dv_solid)
variables[f"{Domain} electrode transfer ratio k"] = dv_th / (dv_sol + 1e-30)
```

Per the PDF (§5.5), this must never feed back into any equation — it's output-only. If a
future measured trajectory shows `k > 1`, that's a sign of missing *additive* reversible
thickness terms (Li plating swelling, cathode lattice strain) outside this partition —
not a reason to let this `k` exceed 1.

### 2.6 What NOT to touch in this phase

- `base_mechanics.py`'s `_compute_stress_and_displacement` (`stress_r_surf = pybamm.Scalar(0)`) — this is PDF §6.7/Change 5, deferred.
- No `"Stack pressure [Pa]"` variable, no `K_stack`/`C_pore(ε)`/`E_s(ϕ)` — PDF §3/§4, deferred.
- `crack_propagation.py`'s stress-coupled cracking law is untouched — it already reads `stress_t_surf`, which isn't modified by this phase.

---

## 3) Parameterisation procedure (assuming experiments are already collected)

**Correction (per your note):** `eps_min_transfer`, `eps_max_transfer`, `f_transmit_min`
(`f0`) are **BoL constants**, determined only from the very first cycle (plus, for
`eps_min_transfer`, complementary BoL microstructural characterization — see §3.3). They
are held fixed for the entire simulation; they are never re-fit against the full-life
trajectory. What evolves over life is `ε_struct(t)`, driven entirely by the
already-calibrated SEI/plating/crack degradation physics (fit earlier this session against
capacity fade, independent of buffering). `f(ε_struct(t))` — and hence the model's
predicted `k(t)` — then evolves as a pure *consequence* of that, not because the three
buffering constants themselves change. A full-life measured `k(EFC)` trajectory, if you
have one, is used only to **validate** this prediction (§3.4 step 5), never to calibrate
it. This replaces the least-squares-over-EFC approach from the previous draft of this plan.

Two data types feed this, matching the PDF's particle-vs-cell distinction (§4.4's `k =
δ_cell/δ_particle`).

### 3.1 Particle-level expansion → `t_change(sto)` per phase (mostly already done)

This calibrates PDF Eq. 10's `Ω(c) dc` (via `t_change`), not the buffering parameters
directly — it's the "unbuffered particle reference" input to everything else.

- Source data: operando/ex-situ XRD lattice-volume vs. stoichiometry, or half-cell
  dilatometry, for each active material separately (graphite, Si/SiOx).
- Format: 2-column CSV (`Sto`, `Volume` as V/V₀), same convention as
  `volume_vs_sto_lithiation.csv` already in `src/pybamm/input/parameters/lithium_ion/data/`.
- Procedure: drop the new CSV into that `data/` folder, load with
  `pybamm.parameters.process_1D_data(filename, path=path)`, wrap in a
  `pybamm.Interpolant(x, y, sto, interpolator="linear")`-returning function exactly as
  `silicon_volume_change_JiaGuo` does in `si_gr_expansion.py`, and point
  `"Primary/Secondary: {Domain} electrode volume change"` at it.
- If the *actual* test cell's electrode blend differs from the literature Ai2020
  graphite curve currently used for the primary phase, get a matched half-cell
  measurement for it too rather than mixing a literature graphite curve with a
  cell-specific Si curve — inconsistent particle-level references will bias the `k` fit
  in §3.3 in a way that's hard to distinguish from genuine buffering.

### 3.2 Cell-level expansion → BoL measurement only

This is the new data needed for `eps_max_transfer`/`f0`. Unlike the previous draft of
this plan, you do **not** need a full-life expansion trajectory to calibrate anything —
only the very first (formation/characterization) cycle.

- Source data: stack/cell thickness (dilatometry, fixture displacement sensor, or
  equivalent), measured over **cycle 1 only**. The same within-cycle amplitude convention
  already used in `test_expansion`'s `expansion_amplitude` (charge-end minus
  discharge-end) is exactly what's needed — a single scalar `δ_cell,BoL` for that first
  cycle, not a life-long series.
- A full-life series is still useful *if you happen to have it*, but only as an optional
  validation check in §3.4 step 5 — it plays no role in determining the three parameters.

### 3.3 Why BoL data alone pins `f0` and `eps_max_transfer`, but not `eps_min_transfer`

During cycle 1, `ε_struct ≈ ε_init` — the SEI/plating/crack terms that drive `ε_struct`
away from `ε_init` are ~0 this early. Eq. 18's `f(ε_struct)` depends only on `ε_struct`
(not on SOC, not on the reversible/breathing porosity), so it is **constant across the
entire first cycle** — even a SOC-resolved BoL dilatometry curve collapses to one number:

```
k_BoL = δ_cell,BoL / δ_particle,BoL
```

(`δ_particle,BoL` computed the same way as before — the model's `dv_solid` over the same
cycle-1 window, using the already-parameterised per-phase `t_change` functions and BoL
`eps_s`.)

One measured number cannot pin three unknowns. Two physically-motivated choices close the
gap without inventing free parameters:

- **Set `eps_max_transfer = ε_init`** (the electrode's own known/measured initial
  porosity — from electrode density and mass loading, or the same characterization pass
  used for `eps_min_transfer` below). This matches the PDF's own framing that a Si/Gr
  cell "starts with `ε_init < ε_max`" and "operates in the shared regime from BoL" —
  pinning `ε_max` at `ε_init` puts BoL exactly at the top of the shared-regime line.
- With `eps_max_transfer = ε_init`, Eq. 18's shared-regime formula evaluated at
  `ε_struct = ε_init` collapses algebraically: the `(ε_struct − eps_min)/(eps_max − eps_min)`
  ratio becomes exactly 1, so `k_BoL = f0` directly — **no fitting needed for `f0`**, just
  read it off the single BoL measurement.
- **`eps_min_transfer` cannot be constrained by BoL electrochemical data at all** — it
  only matters once `ε_struct` drops below `ε_init` as irreversible degradation
  accumulates, which by definition hasn't happened yet at BoL. It needs an independent
  **BoL microstructural characterization** of the pristine electrode's pore network —
  mercury intrusion porosimetry (MIP), BET, or tomography — to identify the porosity at
  which the pore network's throat percolation closes. This is still "BoL data" in your
  sense (measured once, on the fresh electrode, not fit against a life trajectory) but
  it's a structural measurement, not an electrochemical/dilatometry one.

#### Literature check for the assumed `eps_min_transfer ≈ 0.08–0.12` (until that
characterization exists)

No paper found reports a directly-measured percolation/closure porosity for a *cycled,
swollen* Si/graphite composite pore network — the specific quantity Eq. 18/19 needs. The
closest support:

- **Strongest anchor** — "Modeling Volume and Porosity Change: Transitioning from Graphite
  Anodes to Silicon-Dominant Anodes," *J. Electrochem. Soc.* 2026,
  [10.1149/1945-7111/ae5c3d](https://iopscience.iop.org/article/10.1149/1945-7111/ae5c3d)
  (author list not confirmed from search — check the DOI page directly before citing).
  A DFN model coupled to solid mechanics for a Si-dominant anode / NCA cell, with porosity
  updated from active-material + inactive-material + domain-compression volume changes —
  structurally the same idea as this plan's `ε_struct`/partition, just without an explicit
  buffered-vs-transmitted split. At **40% initial anode porosity**, the in-operation
  minimum porosity during charging reaches **~8%**, described as "very close to the
  absolute limit." At a lower **35% initial porosity**, porosity reaches **0%** during C/10
  charging and the simulation terminates in solver errors — i.e. the same failure mode the
  numerical `epsilon_min` floor in `reaction_driven_porosity.py` was added to prevent this
  session, at a design point close to real cells. This is a simulated operational minimum
  under one specific design/rate, not a characterized percolation threshold, but it's the
  best available quantitative match — and it lands almost exactly on the **low end** of
  your assumed range.
- **Qualitative confirmation the mechanism is real** — Pfeifer et al. (contrast-matched
  SANS on Si-graphite anodes), *J. Electrochem. Soc.* 2019,
  [10.1149/2.0781906jes](https://iopscience.iop.org/article/10.1149/2.0781906jes): pores
  become essentially completely filled with SEI product during cycling (pore clogging is
  directly observed), consistent with the notes' `ϕ_min`/closure mechanism, but no
  porosity percentage is reported for the closure point.
- **Directional support, no threshold number** — "Influence of Initial Porosity on the
  Expansion Behavior of Electrodes in Lithium-Ion Batteries," *J. Electrochem. Soc.* 2023,
  [10.1149/1945-7111/acd2fe](https://iopscience.iop.org/article/10.1149/1945-7111/acd2fe):
  dilatometry across graphite (30/40/50%) and Si-graphite (45–65%) initial porosities
  confirms lower initial porosity → more irreversible thickness change transmitted, i.e.
  the same direction as `f(ε)` rising as `ε` falls, but doesn't fit a closure porosity.
- **Don't conflate with a similarly-named but different quantity**: FIB/SEM tomography
  studies of graphite anode material (e.g. "Optimization of Pore Characteristics of
  Graphite-Based Anode for Li-Ion Batteries," *Materials* 2023) report *closed porosity*
  of 21–26% — but that's **intraparticle** closed porosity inside pristine graphite
  particles (pores sealed off within the particle itself, present before any cycling),
  not the **interparticle** electrode-level pore network percolation/closure floor this
  model needs. Same term, different physical quantity — a common trap when literature
  searching this topic.
- General percolation theory (transport-property divergence near a percolation threshold
  `ε_p`) confirms qualitatively that such a threshold exists and that tortuosity/effective
  diffusivity diverge near it, but the threshold value is strongly geometry/coordination-
  number dependent with no universal number that transfers to a specific composite battery
  electrode — which is why the direct simulation match above is more useful here than
  generic percolation theory.

**Recommendation**: keep `eps_min_transfer` in the 0.08–0.12 range as assumed for now,
weighted toward the **0.08–0.10** end given the direct simulation anchor, until the BoL
porosimetry/tomography characterization in §3.4 step 1 replaces it with a measured value
for the actual electrode.

### 3.4 Procedure, concretely

1. **Microstructural characterization** (pristine/BoL electrode): porosimetry or
   tomography on a fresh coupon → `eps_min_transfer` (percolation/closure porosity).
   Independently confirm `ε_init` (electrode density/mass loading, or the same
   porosimetry pass) → `eps_max_transfer = ε_init`.
2. **First-cycle electrochemical + dilatometry**: measure `δ_cell,BoL` over the
   formation/characterization cycle; compute `δ_particle,BoL` over the same cycle from
   the model's already-parameterised `t_change` functions and BoL `eps_s`; set
   `f0 = δ_cell,BoL / δ_particle,BoL`.
3. **Hold `(eps_min_transfer, eps_max_transfer, f0)` fixed** as model constants for the
   rest of the simulation — they are not revisited.
4. **Run the full life simulation** with buffering on, using these fixed constants and
   the already-independently-calibrated degradation parameters from earlier this session
   (composition v7, `SI_MULT`/`GR_DIV`/etc. — untouched by this fit). `ε_struct(t)`
   evolves via the existing SEI/LAM/crack physics; `f(ε_struct(t))` — and hence the
   model's predicted `k(t)` and buffered `"Cell thickness change [m]"` — evolves as pure
   **output**, not because the three constants change.
5. **Validation only, if full-life expansion data exists.** Compare the model-predicted
   `k(EFC)` trajectory against measured full-life data. A mismatch — specifically in
   *when* the sharp late-life rise happens relative to the capacity knee, which is the
   whole point of the precursor claim — is a signal to revisit the `eps_min_transfer`
   *microstructural* estimate (e.g. percolation closure under electrochemical cycling may
   differ from the as-manufactured pristine measurement), not to back-fit it against the
   life trajectory. If it still doesn't match after that, that's the trigger to move to
   the mechanistic `f` (PDF Eq. 24, `1/K_stack / (1/K_stack + C_pore(ε))`) in the deferred
   pressure-coupling phase, rather than continuing to force the reduced three-regime
   version to fit.

### 3.5 What this procedure does *not* need

- No full-life cell-expansion trajectory to calibrate anything (only cycle 1, plus one
  BoL microstructural characterization pass).
- No pressure sensor / stack-force data — deferred phase only.
- No re-fit of the existing LLI/LAM/SEI degradation parameters — those were calibrated
  against capacity fade independently and `ε_struct(t)` is a downstream *output* of them,
  not something this procedure touches.
