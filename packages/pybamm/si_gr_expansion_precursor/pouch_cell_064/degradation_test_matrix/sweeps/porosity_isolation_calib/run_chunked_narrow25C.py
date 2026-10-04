# -*- coding: utf-8 -*-
"""Chunked/resumable runner for the narrow-voltage-window 25 degC check
(check_solvent_diffusion_sei_narrow.py's config), to avoid the long single-
process runs that have twice been killed by the system's low-memory
background-task reaper on this machine.

Runs C64_CHUNK_CYCLES ageing cycles (default 300) per invocation, then saves
the FULL cumulative PyBaMM Solution (via pybamm's own pickle-based
Solution.save/pybamm.load, confirmed to be the supported mechanism -- see
solvers/solution.py) plus a little bookkeeping (total cycles so far, stop
reason) to a SINGLE pickle file at C64_RESUME_PKL. Re-running this same
script with the same C64_RESUME_PKL path picks up exactly where the last
invocation left off (loads the saved solution, continues the ageing-batch
loop with `sim.solve(starting_solution=last_sol)`), and overwrites that same
pickle file with the new, larger cumulative state -- so there is always
exactly one pkl file holding "everything solved so far", not one file per
chunk. Stops early (without waiting for the full C64_CHUNK_CYCLES) if
SOH_TERMINATION_PERCENT or EFC_LIMIT is reached first, exactly like
run_degradation()'s own loop.

Usage: run repeatedly (e.g. from a shell loop) until it prints
"OVERALL STOP" -- each call processes one chunk and returns.
"""
import os
import pickle
import sys

sys.path.insert(
    0,
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))),
        "conditions_test_matrix",
        "narrow_voltage_window",
    ),
)
os.environ.setdefault("C64_SI_BETA_LAM_ISO", "1.0")
os.environ.setdefault("C64_SI_REDIRECT_TO_LAM", "1")
os.environ.setdefault("C64_SI_REDIRECT_LAM_YIELD", "0.02")
os.environ.setdefault("C64_SI_TAU_LAM_ISO", "2e8")
os.environ.setdefault("C64_SI_LAM_ISO_EXPONENT", "3.0")
os.environ.setdefault("C64_SI_OCP_AGING_DEFORM", "1")
os.environ.setdefault("C64_SI_VOLUME_CHANGE_AGING_DEFORM", "1")
os.environ.setdefault("C64_NEG_POROSITY_FLOOR", "0.035")
os.environ.setdefault("C64_EXPONENT_MAX_SEI", "70")
os.environ.setdefault("C64_SI_BETA_LAM_SEI", "1e-7")
os.environ.setdefault("C64_SI_PARTICLE_RADIUS_MULT", "0.2")
os.environ.setdefault("C64_SI_DIFFUSIVITY_MULT", "3")
os.environ.setdefault("C64_NEG_POROSITY_BOL", "0.17")
os.environ.setdefault("C64_SI_K_SEI_MULT", "3e-4")
os.environ.setdefault("C64_GR_K_SEI_MULT", "6e-6")
os.environ.setdefault("C64_SI_CRACK_MULT", os.environ.get("C64_SI_CRACK_MULT_TEST", "1.0"))
os.environ.setdefault("C64_F0", "0.7")
os.environ.setdefault("C64_WIDTH", "0.03")
os.environ.setdefault("C64_SOH_FLOOR", "30")
os.environ.setdefault("C64_NO_PLOT", "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pybamm  # noqa: E402
import cell064_degradation_fit_narrow_window as m  # noqa: E402

m.load_experimental_capacity_fade = lambda: pd.DataFrame(
    columns=["rpt", "efc", "discharge_capacity_Ah", "source", "capacity_retention_pct"])
m.load_experimental_lam = lambda: pd.DataFrame(
    columns=["rpt", "efc", "discharge_Cn_Si", "discharge_Cn_Gr", "discharge_Cp",
             "LAM_NE_graphite_pct", "LAM_NE_silicon_pct", "LAM_PE_pct"])
m.load_experimental_lli = lambda: pd.DataFrame(
    columns=["rpt", "efc", "charge_LLI", "joule_LLI", "charge_LLI_loss", "joule_LLI_loss"])
m.load_experimental_reversible_expansion = lambda: pd.DataFrame(
    columns=["cyc_seq", "efc", "reversible_expansion_um"])
m.load_experimental_k_expansion = lambda: pd.DataFrame(columns=["rpt", "efc", "k", "source"])
m.plot_si_ocp_aging = lambda *a, **k: None

m.MODEL_OPTIONS_BASE["SEI"] = (("ec reaction limited", "solvent-diffusion limited"), "none")

SOLVENT_MULT = float(os.environ.get("C64_SOLVENT_MULT", 1.0))
m.PARAM_UPDATES["Secondary: SEI solvent diffusivity [m2.s-1]"] = 2.5e-22 * SOLVENT_MULT
AMBIENT_T = float(os.environ.get("C64_AMBIENT_T", 298.15))
m.PARAM_UPDATES["Ambient temperature [K]"] = AMBIENT_T
m.PARAM_UPDATES["Initial temperature [K]"] = AMBIENT_T
SI_ESEI = float(os.environ.get("C64_SI_ESEI", 38000.0))
m.PARAM_UPDATES["Secondary: SEI growth activation energy [J.mol-1]"] = SI_ESEI
GR_ESEI = float(os.environ.get("C64_GR_ESEI", 38000.0))
m.PARAM_UPDATES["Primary: SEI growth activation energy [J.mol-1]"] = GR_ESEI

SI_LAM_PROP_BASE = float(os.environ.get("C64_SI_LAM_PROP", 3.078e-8))
SI_LAM_EAC = float(os.environ.get("C64_SI_LAM_EAC", 0.0))


def _silicon_beta_lam_scaled(T_dim):
    arrhenius = np.exp(SI_LAM_EAC / 8.314 * (1 / 298.15 - 1 / T_dim))
    return SI_LAM_PROP_BASE * arrhenius


if SI_LAM_EAC != 0.0:
    m.PARAM_UPDATES["Secondary: Negative electrode LAM constant proportional term [s-1]"] = (
        _silicon_beta_lam_scaled
    )

print(f"SOLVENT_MULT={SOLVENT_MULT}, SI_ESEI={SI_ESEI}, GR_ESEI={GR_ESEI}, "
      f"SI_LAM_EAC={SI_LAM_EAC}, AMBIENT_T={AMBIENT_T}, "
      f"voltage window {m.LOWER_CUTOFF_V}-{m.UPPER_CUTOFF_V} V", flush=True)

CHUNK_CYCLES = int(os.environ.get("C64_CHUNK_CYCLES", 300))
RESUME_PKL = os.environ.get("C64_RESUME_PKL")
if not RESUME_PKL:
    raise SystemExit("C64_RESUME_PKL must be set to a file path")
assert CHUNK_CYCLES % m.BATCH_SIZE == 0, "C64_CHUNK_CYCLES must be a multiple of BATCH_SIZE"


def make_sim(experiment, this_solver, model, param):
    return pybamm.Simulation(
        model, parameter_values=param, experiment=experiment,
        solver=this_solver, var_pts=m.VAR_PTS,
    )


if os.path.exists(RESUME_PKL):
    with open(RESUME_PKL, "rb") as f:
        state = pickle.load(f)
    last_sol = state["sol"]
    total_cycles_requested = state["total_cycles_requested"]
    batch_num = state["batch_num"]
    print(f"Resumed from {RESUME_PKL}: {total_cycles_requested} cycles already "
          f"completed ({batch_num} batches).", flush=True)
else:
    options = dict(m.MODEL_OPTIONS_BASE)
    model = pybamm.lithium_ion.DFN(options)
    param = m.build_parameter_values()
    sim = make_sim(m.formation_exp, m.solver, model, param)
    last_sol = sim.solve(initial_soc=1.0)
    total_cycles_requested = 0
    batch_num = 0
    print("No existing state -- solved formation cycle fresh.", flush=True)

options = dict(m.MODEL_OPTIONS_BASE)
model = pybamm.lithium_ion.DFN(options)
param = m.build_parameter_values()

stop_reason = None
chunk_start_cycles = total_cycles_requested
while total_cycles_requested - chunk_start_cycles < CHUNK_CYCLES and total_cycles_requested < m.MAX_TOTAL_CYCLES:
    batch_num += 1
    new_sol = None
    last_exc = None
    for attempt, tol_mult in enumerate([1.0] + m.RETRY_TOL_MULTIPLIERS):
        try_tol = m.SOLVER_TOL * tol_mult
        this_solver = (m.solver if tol_mult == 1.0
                       else pybamm.IDAKLUSolver(
                           root_tol=try_tol, atol=try_tol, rtol=try_tol,
                           output_variables=m._CORE_OUTPUT_VARIABLES if m._RESTRICT_OUTPUT_VARIABLES else None,
                       ))
        sim = make_sim(m.ageing_batch_exp, this_solver, model, param)
        try:
            new_sol = sim.solve(starting_solution=last_sol)
            if attempt > 0:
                print(f"[BATCH {batch_num}] recovered on retry {attempt} (tol={try_tol:.0e})", flush=True)
            break
        except Exception as exc:  # noqa: BLE001
            last_exc = exc
            continue
    if new_sol is None:
        print(f"[BATCH {batch_num}] solver failure after {total_cycles_requested} "
              f"successfully completed ageing cycles: {type(last_exc).__name__}: {str(last_exc)[:200]}",
              flush=True)
        stop_reason = f"solver_failure_batch_{batch_num}"
        break

    last_sol = new_sol
    total_cycles_requested += m.BATCH_SIZE

    age_cap = []
    for cyc in last_sol.cycles:
        cp, rate, _, _ = m.cycle_ageing_leg(cyc)
        if cp > 0 and rate > 0.5:
            age_cap.append(cp)
    if age_cap:
        soh = 100.0 * age_cap[-1] / age_cap[0]
        print(f"[BATCH {batch_num}] {total_cycles_requested} ageing cycles requested, "
              f"{len(age_cap)} completed, SoH={soh:.2f}%", flush=True)
        if soh <= m.SOH_TERMINATION_PERCENT:
            stop_reason = f"reached_{m.SOH_TERMINATION_PERCENT:.0f}pct_floor"
            print("Reached SoH floor -- stopping.", flush=True)
            break
    else:
        print(f"[BATCH {batch_num}] {total_cycles_requested} ageing cycles requested, "
              "no completed ageing legs found yet", flush=True)

    if np.isfinite(m.EFC_LIMIT):
        current_efc = float(m.efc_from_throughput(last_sol["Throughput capacity [A.h]"].entries[-1]))
        if current_efc >= m.EFC_LIMIT:
            stop_reason = f"reached_efc_limit_{m.EFC_LIMIT:.0f}"
            print(f"[BATCH {batch_num}] EFC={current_efc:.1f} >= C64_EFC_LIMIT={m.EFC_LIMIT:.1f} "
                  "-- stopping.", flush=True)
            break

with open(RESUME_PKL, "wb") as f:
    pickle.dump({"sol": last_sol, "total_cycles_requested": total_cycles_requested,
                 "batch_num": batch_num}, f, protocol=pickle.HIGHEST_PROTOCOL)
print(f"Saved cumulative state to {RESUME_PKL}: {total_cycles_requested} cycles, "
      f"{batch_num} batches.", flush=True)

if stop_reason is None and total_cycles_requested >= m.MAX_TOTAL_CYCLES:
    stop_reason = "reached MAX_TOTAL_CYCLES"

if stop_reason is not None:
    print(f"OVERALL STOP: {stop_reason}", flush=True)
else:
    print(f"CHUNK DONE (not overall stop) -- run again with the same C64_RESUME_PKL "
          f"to continue from {total_cycles_requested} cycles.", flush=True)

try:
    current_efc = float(m.efc_from_throughput(last_sol["Throughput capacity [A.h]"].entries[-1]))
    print(f"Cumulative EFC so far: {current_efc:.1f}", flush=True)
except Exception as exc:  # noqa: BLE001
    print(f"(EFC readout skipped: {exc})", flush=True)
