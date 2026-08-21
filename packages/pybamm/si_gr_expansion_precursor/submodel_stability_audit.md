# Task C — submodel stability audit (read-only)

Composite Si/Gr DFN degradation model, Mayur2024 parameter set. Goal: identify
submodel expressions that can go singular/stiff as the cell fades toward and
past an ~88% SoH knee. No submodel files were edited — this is a candidate
list only, per the plan's guardrails.

| File | Expression (line#) | Failure mode | Candidate softening |
|---|---|---|---|
| `active_material/loss_active_material.py` | `eps_solid = pybamm.Variable(...)` (36-46) — no `bounds=`, unlike `c_s` in `fickian_diffusion.py` | Stress-driven `deps_solid_dt` (115-119) has no floor; under sustained tensile stress it can integrate through 0. Then `a = 3*eps_solid/R` (`base_active_material.py:103`) → 0 or negative, and since `a*j = i_applied` locally, `j = i/a → ∞` — a current-focusing singularity. **This is the most likely single point of failure near the knee** (see Task A diagnostic notes). | Add `bounds=(eps_s_min, eps_s_max)` to the Variable (mirrors `c_s`); brake `deps_solid_dt` multiplicatively as `eps_solid→eps_s_floor` (same softplus idea as the porosity floor). |
| `active_material/loss_active_material.py` | `j_stress_LAM = -beta_LAM * ((stress_h_surf_tensile - stress_h_surf_min)/stress_critical)**m_LAM` (114-119) | `m_LAM` typically fractional; derivative blows up near `stress=0`, adding stiffness once combined with the unfloored integration above. | Floor `eps_solid` (above); optionally smooth the `stress>0` mask with a softplus instead of a hard comparison. |
| `particle/fickian_diffusion.py` + `base_particle.py` | `bounds=(0, c_max)` on `c_s` (fickian_diffusion.py:46,55,71,97) | `bounds` only trips a solver *event* — it does not clamp. Mid-step residual evaluations can transiently see `c_s<0` or `c_s>c_max`, feeding straight into `sqrt(c_max - c_surf)` in kinetics → NaN. Silicon's `D(c,T)` and stress factor are also non-monotonic near `x→0/1`. | Add a smooth clamp on `c_s_surf` before it's used in kinetics/OCP, not just on the state variable itself. |
| `interface/kinetics/base_kinetics.py`, `butler_volmer.py`, `Mayur2024.py` | exchange-current-density `... * (c_s_max - c_s_surf)**0.5` (Mayur2024.py:161, ~895) | As `x→1`, `j0→0` and `d(j0)/dc_s_surf → -∞` (square-root cusp) — step-size collapse right where near-surface saturation is common in Si. | Floor the sqrt argument: `sqrt(max(c_s_max - c_s_surf, eps*c_max))`. |
| `open_circuit_potential/` (silicon OCP fns in `Mayur2024.py`) | `silicon_ocp_lithiation_Mark2016` / `nmc_LGM50_ocp_Chen2020`: `+ 1e-4*(1/sto + 1/(sto-1))` (lines 288, 518). Note `silicon_ocp_delithiation_Mark2016` has **no** such term. | Regularizer diverges (not saturates) as `sto→0/1`. The "current sigmoid" hysteresis blend is itself smooth, but it linearly combines lithiation/delithiation branches, so the singular lithiation branch dominates regardless of current direction. **Tested directly (Task A): disabling this hysteresis did not fix the observed crash**, so it's a latent risk, not the primary driver here. | Clip `sto` to `[eps, 1-eps]` before evaluating the OCP polynomial, or replace with a bounded penalty. |
| `porosity/reaction_driven_porosity.py` | `eps_min = pybamm.Scalar(0.08)`, `k=100` softplus (85-87) | Floor is active and correctly wired in. **Mismatch found:** `knee_sweep.py`'s `BASE_OVERRIDES` sets `"Negative/Positive electrode minimum porosity": 0.10`, but this file never reads that parameter — the real floor is the hardcoded `0.08`, silently overriding the intended `0.10`. | Either read `eps_min` from a `pybamm.Parameter(...)` here, or drop the dead override from `knee_sweep.py` so it doesn't imply a floor that isn't active. |
| `particle_mechanics/crack_propagation.py` | anti-Paris `stress_relief = max(1-l_cr/R_typ, 0)` (118-122); event `1-max(l_cr)/R_typ` (196) | Confirmed working as intended: caps `l_cr < R_typ`, `dl_cr→0` before saturation. The radius event should no longer fire in normal operation. | None needed — patch looks complete. Optionally soften the hard `(stress_eff>=0)` mask for a smoother Jacobian. |
| `particle_mechanics/base_mechanics.py` | stress BC (`sigma_r(R)=0`) fine; but `eps_s` (unfloored, see row 1) feeds thickness-change/roughness/`a_cr` aggregation used by SEI-on-cracks | A negative `eps_solid` propagates unphysical sign flips into thickness/roughness/crack-area terms. | Resolved transitively by flooring `eps_solid`. |
| `interface/sei/sei_growth.py` | reaction-limited `j_sei` (no self-limiting concentration term, unlike EC-diffusion-limited SEI) | On a growing crack area, `a_cr`-scaled `j_sei` is unbounded except indirectly via crack-length saturation and SEI-resistance feedback. Plausibly the real driver of local Li consumption/porosity stress near the knee — more likely than LAM given LAM is deliberately tame in this parameter set. | Cap `a_cr`-scaled growth with a soft saturation on SEI resistance growth rate, or use `ec reaction limited` for the cracked-area contribution specifically. |
| `electrolyte_diffusion/`, `electrolyte_conductivity/` (Bruggeman) | `tor = eps**b` reads the same floored `"...porosity"` variable populated by `reaction_driven_porosity.py` | Confirmed: no unfloored bypass found for this lithium-ion DFN path (the unfloored `ReactionDrivenODE` variant is lead-acid only). | None needed for this model; latent trap only if a future change swaps in the unfloored porosity ODE. |

## Summary

The most probable single point of failure is the **unbounded `eps_solid`** in
`active_material/loss_active_material.py`. It's the only state variable in
this model's degradation stack without a floor/bounds, and the mechanism it
enables — `a = 3·eps_solid/R → 0` forcing `j = i/a → ∞` — is consistent with
everything observed in the Task A diagnostic: a reproducible, tolerance-
independent, stoichiometry-independent, hysteresis-independent solver
collapse (`IDA_ERR_FAIL`/`IDA_CONV_FAIL`) right in the aggressive corner of
the parameter space the knee-tuning sweep needs to explore.
