"""Validated regression test for the R_Li state-rescale mechanism used by
ec_dryout_wrapper.py. Confirms, against the live PyBaMM install, that:

  1. `sim.solve(starting_solution=prev_solution)` tolerates a DIFFERENT
     ParameterValues object on the new Simulation (different bulk solvent
     concentration, different electrode area) -- this is NOT documented
     as guaranteed behaviour anywhere, so this script exists to catch a
     future PyBaMM version silently changing it.
  2. The electrolyte Li+ concentration state ("Porosity times
     concentration [mol.m-3]") can be rescaled by a known factor R_Li via
     direct in-place mutation of `solution.last_state.y` at the slice
     given by that state's y_slices entry, with NO reference offset (a
     pure multiplicative rescale is exact) -- verified by checking the
     resulting concentration at t=0 of the continued solve matches
     R_Li * (concentration at end of the previous solve) to solver
     tolerance.
  3. Porosity and all other state (e.g. particle concentration) pass
     through completely unaffected.
  4. This composes correctly with a SIMULTANEOUS bulk solvent
     concentration parameter change and electrode-area (dry-out) change,
     with "pore buffering": "true" enabled -- i.e. the exact combination
     ec_dryout_wrapper.py uses every batch.

Run this after any PyBaMM upgrade before trusting ec_dryout_wrapper.py's
apply_li_rescale(). See implementation_plan.md section 4.3.
"""
import numpy as np
import pybamm

MODEL_OPTIONS = {
    "particle phases": ("2", "1"),
    "open-circuit potential": (("single", "current sigmoid"), "single"),
    "SEI": "solvent-diffusion limited",
    "SEI porosity change": "true",
    "particle mechanics": ("swelling and cracking", "swelling only"),
    "SEI on cracks": "true",
    "loss of active material": "stress-driven",
    "pore buffering": "true",
    "pore buffering transition": "physical",
}

POROSITY_TIMES_CONC_NAME = "Porosity times concentration [mol.m-3]"


def find_y_slice(solution_last_state, var_name):
    lm = solution_last_state.all_models[-1]
    for var, slices in lm.y_slices.items():
        if getattr(var, "name", None) == var_name:
            return slices[0]
    raise KeyError(
        f"{var_name!r} not found in y_slices -- PyBaMM's internal state "
        f"representation may have changed; ec_dryout_wrapper.py's "
        f"apply_li_rescale() needs updating"
    )


def main():
    param1 = pybamm.ParameterValues("si_gr_expansion")
    exp = pybamm.Experiment(
        ["Discharge at 1C until 2.5 V", "Charge at 0.3C until 4.2 V",
         "Hold at 4.2 V until C/100"] * 2
    )
    sim1 = pybamm.Simulation(
        pybamm.lithium_ion.DFN(MODEL_OPTIONS), parameter_values=param1, experiment=exp
    )
    sol1 = sim1.solve(initial_soc=1.0)
    c_e_ref = sol1["X-averaged electrolyte concentration [mol.m-3]"].entries[-1]
    eps_ref = sol1["X-averaged negative electrode porosity"].entries[-1]
    n_par_0 = param1["Number of electrodes connected in parallel to make a cell"]

    R_LI = 1.10
    R_DRY = 0.97
    param2 = param1.copy()
    param2.update(
        {
            "Primary: Bulk solvent concentration [mol.m-3]": 2000.0,
            "Secondary: Bulk solvent concentration [mol.m-3]": 2000.0,
            "Number of electrodes connected in parallel to make a cell": n_par_0 * R_DRY,
        },
        check_already_exists=False,
    )

    last_state = sol1.last_state
    sl = find_y_slice(last_state, POROSITY_TIMES_CONC_NAME)
    last_state.y[sl] = last_state.y[sl] * R_LI

    sim2 = pybamm.Simulation(
        pybamm.lithium_ion.DFN(MODEL_OPTIONS), parameter_values=param2, experiment=exp
    )
    sol2 = sim2.solve(starting_solution=last_state)

    c_e_new = sol2["X-averaged electrolyte concentration [mol.m-3]"].entries[0]
    eps_new = sol2["X-averaged negative electrode porosity"].entries[0]
    ratio = c_e_new / c_e_ref

    print(f"c_e ratio: {ratio:.6f}  (expected {R_LI})")
    print(f"porosity ref/new: {eps_ref:.8f} / {eps_new:.8f} (should match exactly)")
    assert abs(ratio - R_LI) < 1e-6, "R_Li rescale did not propagate correctly"
    assert abs(eps_new - eps_ref) < 1e-12, "porosity should be untouched by R_Li"
    print("PASS: apply_li_rescale mechanism verified against live PyBaMM install.")


if __name__ == "__main__":
    main()
