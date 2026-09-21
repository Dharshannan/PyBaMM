"""Quick standalone check: does reducing 'Negative electrode porosity' alone
(leaving Primary/Secondary active-material volume fractions untouched) shift
the BoL C/20 discharge capacity/voltage/dV-dQ curve? Porosity only enters the
electrolyte-transport equations, whose effect should be small at C/20's very
low rate -- this verifies that empirically before trusting any degradation
run at a reduced porosity. Uses the exact same model options/VAR_PTS/
degradation-recipe PARAM_UPDATES as cell064_degradation_fit.py (only the
porosity override differs per run), so this is a like-for-like BoL check.

Usage:
    "<PyBaMM_Pressure conda env>/python" _porosity_bol_check.py
"""
import os
import sys

import numpy as np
import pybamm

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
import cell064_degradation_fit as c64  # noqa: E402

UPPER_CUTOFF_V = c64.UPPER_CUTOFF_V
LOWER_CUTOFF_V = c64.LOWER_CUTOFF_V

c20_exp = pybamm.Experiment(
    [
        f"Charge at C/3 until {UPPER_CUTOFF_V} V",
        f"Hold at {UPPER_CUTOFF_V} V until C/50",
        f"Discharge at C/20 until {LOWER_CUTOFF_V} V",
    ]
)

results = {}
for label, porosity in [("baseline (0.35)", None), ("reduced (0.25)", 0.25), ("reduced (0.15)", 0.15)]:
    model = pybamm.lithium_ion.DFN(dict(c64.MODEL_OPTIONS_BASE))
    param = c64.build_parameter_values()
    if porosity is not None:
        param.update({"Negative electrode porosity": porosity}, check_already_exists=False)
    sim = pybamm.Simulation(model, experiment=c20_exp, parameter_values=param,
                             var_pts=c64.VAR_PTS, solver=pybamm.IDAKLUSolver(root_tol=1e-06, atol=1e-06, rtol=1e-06))
    sol = sim.solve()
    disch = sol.cycles[-1].steps[-1]
    cap = disch["Discharge capacity [A.h]"].entries
    v = disch["Terminal voltage [V]"].entries
    q_final = cap[-1] - cap[0]
    results[label] = (cap - cap[0], v)
    print(f"{label}: C/20 discharge capacity = {q_final:.4f} Ah", flush=True)

base_cap, base_v = results["baseline (0.35)"]
for label in ["reduced (0.25)", "reduced (0.15)"]:
    cap, v = results[label]
    v_interp = np.interp(base_cap, cap, v)
    max_dv = np.max(np.abs(v_interp - base_v))
    print(f"{label} vs baseline: max |dV| over full discharge = {max_dv * 1000:.2f} mV, "
          f"capacity diff = {(cap[-1] - base_cap[-1]) * 1000:.2f} mAh", flush=True)
