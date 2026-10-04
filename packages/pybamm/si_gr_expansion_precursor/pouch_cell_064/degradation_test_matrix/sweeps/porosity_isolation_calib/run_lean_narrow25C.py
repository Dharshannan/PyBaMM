# -*- coding: utf-8 -*-
"""Single-process, memory-lean runner for the narrow-voltage-window 25 degC
verification (no chunking/pickling needed -- 2026-09-27, after the 55GB
pickle problem ruled that out). Splits each batch into TWO solve calls
instead of one: 49 ageing-only cycles with a LEAN solver
(store_first_last=True -- only first/last sample per step, not every
internal timestep), then the 1 RPT cycle (always the last position in a
batch, since RPT_INTERVAL==BATCH_SIZE==50 here) with the FULL-RESOLUTION
solver (needed for rpt_discharge_curve()'s V(Q) shape). Both chained via
starting_solution, so the final solution has full within-step resolution
only for RPT cycles and first/last-only for the other 49/50 of cycles --
should cut memory roughly ~50x for the bulk of a long run while keeping
every RPT's voltage curve intact.
"""
import os
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
os.environ.setdefault("C64_NO_PLOT", "0")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pybamm  # noqa: E402
import cell064_degradation_fit_narrow_window as m  # noqa: E402

# No real data overlay (exploratory mechanism test).
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

assert m.RPT_INTERVAL == m.BATCH_SIZE, (
    "This lean runner assumes exactly one RPT at the END of each batch "
    "(RPT_INTERVAL == BATCH_SIZE); adjust if that assumption changes."
)

# The narrow-voltage-window fork predates the output_variables memory
# optimization added to the shared script this session -- it has no
# _CORE_OUTPUT_VARIABLES/_RESTRICT_OUTPUT_VARIABLES of its own, so this
# runner defines its own copy rather than depending on fork internals that
# don't exist there.
_CORE_OUTPUT_VARIABLES = [
    "Current [A]",
    "Discharge capacity [A.h]",
    "Throughput capacity [A.h]",
    "Terminal voltage [V]",
    "Cell thickness change [m]",
    "Negative electrode transfer ratio k",
    "Volume-averaged cell temperature [K]",
    "Local ECM resistance [Ohm]",
    "X-averaged negative primary particle surface stoichiometry",
    "X-averaged negative secondary particle surface stoichiometry",
    "Loss of active material in negative electrode [%]",
    "Loss of active material in positive electrode [%]",
    "Loss of active material in primary phase in negative electrode [%]",
    "Loss of active material in secondary phase in negative electrode [%]",
    "Total capacity lost to side reactions [A.h]",
    "Loss of lithium due to loss of primary active material in negative electrode [mol]",
    "Loss of lithium due to loss of secondary active material in negative electrode [mol]",
    "Loss of lithium due to loss of active material in positive electrode [mol]",
]
lean_solver = pybamm.IDAKLUSolver(
    root_tol=m.SOLVER_TOL, atol=m.SOLVER_TOL, rtol=m.SOLVER_TOL,
    output_variables=_CORE_OUTPUT_VARIABLES,
    store_first_last=True,
)
full_res_solver = pybamm.IDAKLUSolver(
    root_tol=m.SOLVER_TOL, atol=m.SOLVER_TOL, rtol=m.SOLVER_TOL,
    output_variables=_CORE_OUTPUT_VARIABLES,
)

ageing_only_exp = pybamm.Experiment(
    [m._batch_cycle_steps(i) for i in range(1, m.BATCH_SIZE)]  # 1..BATCH_SIZE-1: all plain ageing legs
)
rpt_only_exp = pybamm.Experiment(
    [m._batch_cycle_steps(m.BATCH_SIZE)]  # position BATCH_SIZE: the RPT leg
)


def make_sim(experiment, this_solver, model, param):
    return pybamm.Simulation(
        model, parameter_values=param, experiment=experiment,
        solver=this_solver, var_pts=m.VAR_PTS,
    )


options = dict(m.MODEL_OPTIONS_BASE)
model = pybamm.lithium_ion.DFN(options)
param = m.build_parameter_values()

sim = make_sim(m.formation_exp, full_res_solver, model, param)
last_sol = sim.solve(initial_soc=1.0)
print("Formation cycle solved.", flush=True)

EFC_LIMIT = m.EFC_LIMIT
MAX_TOTAL_CYCLES = m.MAX_TOTAL_CYCLES
SOH_TERMINATION_PERCENT = m.SOH_TERMINATION_PERCENT

total_cycles_requested = 0
stop_reason = "reached MAX_TOTAL_CYCLES"
batch_num = 0
while total_cycles_requested < MAX_TOTAL_CYCLES:
    batch_num += 1

    # Part A: 49 ageing-only cycles, lean solver.
    sim = make_sim(ageing_only_exp, lean_solver, model, param)
    try:
        mid_sol = sim.solve(starting_solution=last_sol)
    except Exception as exc:  # noqa: BLE001
        print(f"[BATCH {batch_num}] ageing-leg solver failure: {type(exc).__name__}: {str(exc)[:200]}",
              flush=True)
        stop_reason = f"solver_failure_batch_{batch_num}_ageing"
        break

    # Part B: 1 RPT cycle, full-resolution solver.
    sim = make_sim(rpt_only_exp, full_res_solver, model, param)
    try:
        new_sol = sim.solve(starting_solution=mid_sol)
    except Exception as exc:  # noqa: BLE001
        print(f"[BATCH {batch_num}] RPT-leg solver failure: {type(exc).__name__}: {str(exc)[:200]}",
              flush=True)
        stop_reason = f"solver_failure_batch_{batch_num}_rpt"
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
        if soh <= SOH_TERMINATION_PERCENT:
            stop_reason = f"reached_{SOH_TERMINATION_PERCENT:.0f}pct_floor"
            print("Reached SoH floor -- stopping.", flush=True)
            break
    else:
        print(f"[BATCH {batch_num}] {total_cycles_requested} ageing cycles requested, "
              "no completed ageing legs found yet", flush=True)

    if np.isfinite(EFC_LIMIT):
        current_efc = float(m.efc_from_throughput(last_sol["Throughput capacity [A.h]"].entries[-1]))
        if current_efc >= EFC_LIMIT:
            stop_reason = f"reached_efc_limit_{EFC_LIMIT:.0f}"
            print(f"[BATCH {batch_num}] EFC={current_efc:.1f} >= C64_EFC_LIMIT={EFC_LIMIT:.1f} "
                  "-- stopping.", flush=True)
            break

print(f"Stop reason: {stop_reason}", flush=True)

results = m.extract_results(last_sol)
m.plot_all(results, last_sol)


def rpt_based_knee(results):
    rpt_cap = results["rpt_cap_arr"]
    rpt_thr = results["rpt_thr"]
    if rpt_cap.size < 2:
        print("RPT-based knee: unavailable (fewer than 2 RPTs)")
        return None
    rpt_efc = m.efc_from_throughput(rpt_thr)
    soh = 100.0 * rpt_cap / rpt_cap[0]
    slopes = np.diff(soh) / np.diff(rpt_efc)
    idx = int(np.argmin(slopes))
    knee_efc = rpt_efc[idx]
    print(f"RPT-based knee: steepest interval is RPT#{idx}->RPT#{idx + 1} "
          f"(EFC {rpt_efc[idx]:.1f}->{rpt_efc[idx + 1]:.1f}, "
          f"SoH {soh[idx]:.1f}%->{soh[idx + 1]:.1f}%, "
          f"slope={slopes[idx]:.3f} pp/EFC) -> knee EFC~{knee_efc:.1f}")
    return knee_efc


rpt_based_knee(results)
