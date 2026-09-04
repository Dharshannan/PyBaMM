"""
dma_baseline_run.py -- runs the project's tuned "ec reaction limited" anchor recipe
(k_sei x0.0017, f0=0.7, width=0.01 -- identical to ../ec_reaction_limited_swap_test.py /
../sweeps/ec_reaction_limited_ksei_sweep_test.py's best-match point) to 50% SoH, with a
C/20 pseudo-OCV RPT inserted every 10 ageing cycles (mirroring test_pore_buffering/
pore_buffering_degradation_test_1000cyc.py's RPT convention, at the slower rate this
project's DMA pipeline needs -- see dma_ocp_fit.py's module docstring for why).

Saves ALL of the following to a single .npz, for reuse by both dma_method_a_plot.py
(direct model-state-variable degradation modes) and dma_method_b_plot.py (composite-DMA
curve-fit degradation modes, from the RPT pOCV curves saved here) -- the expensive
simulation only needs to run once:
  - age_thr, age_cap, age_LLI, age_LAM_gr, age_LAM_si, age_LAM_pos, age_amplitude, age_k:
    dense, per-(C/3-ageing-)cycle series (used for the smooth "knee reconstruction" line
    and the amplitude/k panel, identically in both plots).
  - rpt_thr, rpt_cap: per-RPT (C/20) discharge capacity and its throughput.
  - rpt_LLI, rpt_LAM_gr, rpt_LAM_si, rpt_LAM_pos: per-RPT direct model-state values
    (Method A), interpolated onto each RPT's own throughput for an apples-to-apples
    comparison against Method B (which can only ever be evaluated AT RPT points).
  - rpt_Q_list, rpt_V_list: the actual (discharge capacity, voltage) arrays for each RPT's
    C/20 discharge leg -- the raw pseudo-OCV curves dma_method_b_plot.py fits.
  - nominal_cap, EFC = throughput / (2*nominal_cap) conversion is done downstream (kept as
    raw throughput here so both plotting scripts share one conversion point).
"""
import os

import numpy as np
import pybamm

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

pybamm.set_logging_level("NOTICE")

MODEL_OPTIONS_BASE = {
    "particle phases": ("2", "1"),
    "open-circuit potential": (("single", "current sigmoid"), "single"),
    "SEI": "ec reaction limited",
    "SEI porosity change": "true",
    "particle mechanics": ("swelling and cracking", "swelling only"),
    "SEI on cracks": "true",
    "loss of active material": "stress-driven",
    "pore buffering": "true",
    "pore buffering transition": "physical",
}

CAPACITY_PROBE_OPTIONS = {
    "particle phases": ("2", "1"),
    "open-circuit potential": (("single", "current sigmoid"), "single"),
}

VAR_PTS = {
    "x_n": 5, "x_s": 5, "x_p": 5,
    "r_n_prim": 20, "r_n_sec": 20, "r_n": 20, "r_p": 20,
}

SI_MAX_CONC_DEFAULT = 278000.0
NOMINAL_CAP_AH = 5.0


def scaled_si_initial_conc(si_max_conc):
    return 275220.0 * (si_max_conc / SI_MAX_CONC_DEFAULT)


def measure_capacity_constants(si_max_conc):
    param = pybamm.ParameterValues("si_gr_expansion")
    param.update(
        {
            "Secondary: Maximum concentration in negative electrode [mol.m-3]": si_max_conc,
            "Secondary: Initial concentration in negative electrode [mol.m-3]": scaled_si_initial_conc(si_max_conc),
        },
        check_already_exists=False,
    )
    EPS_PROBE = 0.1
    p = param.copy()
    p.update(
        {
            "Primary: Negative electrode active material volume fraction": EPS_PROBE,
            "Secondary: Negative electrode active material volume fraction": EPS_PROBE,
        }
    )
    model = pybamm.lithium_ion.DFN(CAPACITY_PROBE_OPTIONS)
    sim = pybamm.Simulation(
        model, parameter_values=p, var_pts=VAR_PTS,
        experiment=pybamm.Experiment(["Rest for 1 second"]),
    )
    sol = sim.solve()
    cap_gr = float(np.squeeze(sol["Negative electrode primary phase capacity [A.h]"].entries[0]))
    cap_si = float(np.squeeze(sol["Negative electrode secondary phase capacity [A.h]"].entries[0]))
    return cap_gr / EPS_PROBE, cap_si / EPS_PROBE


def solve_capacity_split_at_fixed_volume(si_vol_frac, cap_gr_frac, k_gr, target_cap):
    gr_eps = target_cap * cap_gr_frac / k_gr
    si_eps = (si_vol_frac / (1 - si_vol_frac)) * gr_eps
    cap_si_target = target_cap * (1 - cap_gr_frac)
    k_si_needed = cap_si_target / si_eps
    si_max_conc_needed = SI_MAX_CONC_DEFAULT * (k_si_needed / k_si_orig_global)
    return gr_eps, si_eps, si_max_conc_needed


print("Measuring capacity constants at DEFAULT Si max concentration (pure Si)...")
k_gr_orig, k_si_orig = measure_capacity_constants(SI_MAX_CONC_DEFAULT)
k_si_orig_global = k_si_orig

CAP_GR_FRAC = 0.55
SI_VOL_FRAC = 0.20
GR_EPS, SI_EPS, SI_MAX_CONC_NEEDED = solve_capacity_split_at_fixed_volume(
    SI_VOL_FRAC, CAP_GR_FRAC, k_gr_orig, NOMINAL_CAP_AH)
print(f"Solved composition: Gr_eps={GR_EPS:.5f}  Si_eps={SI_EPS:.5f}  "
      f"Si_max_conc={SI_MAX_CONC_NEEDED:.1f} mol/m3")

SI_CRIT_STRESS = 2.2e8
SI_LAM_PROP_BASELINE = 7.5e-7
SI_LAM_EXP = 2.5
GR_CRIT_STRESS = 3.0e7
GR_LAM_PROP_BASELINE = 4.5e-6
SI_CRACK_RATE_MULT = 0.1
GR_CRACK_RATE_MULT = 0.1
BASE_CRACK_RATE = 3.9e-20
NEG_POROSITY_FLOOR = 0.01
EXPONENT_MAX_SEI = 10.0
F0_BASELINE = 0.7
WIDTH_BASELINE = 0.01
K_SEI_DEFAULT = 1e-12
K_SEI_ANCHOR_MULT = 0.0017

BATCH_SIZE = 10           # = RPT_INTERVAL -- one C/20 RPT per batch
RPT_INTERVAL = 10
RPT_RATE = "C/20"
MAX_TOTAL_CYCLES = 400
SOH_TARGET = 0.5

solver = pybamm.IDAKLUSolver(root_tol=1e-06, atol=1e-06, rtol=1e-06)

formation_exp = pybamm.Experiment(
    [
        "Discharge at 0.1C until 2.5 V",
        "Charge at C/3 until 4.2 V",
        "Hold at 4.2 V until C/100",
    ]
)


def _batch_cycle_steps(local_position):
    if local_position % RPT_INTERVAL == 0:
        return (
            f"Discharge at {RPT_RATE} until 2.5 V",
            "Charge at C/3 until 4.2 V",
            "Hold at 4.2 V until C/100",
        )
    return (
        "Discharge at C/3 until 2.5 V",
        "Charge at C/3 until 4.2 V",
        "Hold at 4.2 V until C/100",
    )


ageing_batch_exp = pybamm.Experiment(
    [_batch_cycle_steps(i) for i in range(1, BATCH_SIZE + 1)]
)


def build_param():
    param = pybamm.ParameterValues("si_gr_expansion")
    param.update(
        {
            "Primary: Negative electrode active material volume fraction": GR_EPS,
            "Secondary: Negative electrode active material volume fraction": SI_EPS,
            "Secondary: Maximum concentration in negative electrode [mol.m-3]": SI_MAX_CONC_NEEDED,
            "Secondary: Initial concentration in negative electrode [mol.m-3]": scaled_si_initial_conc(SI_MAX_CONC_NEEDED),
            "Nominal cell capacity [A.h]": NOMINAL_CAP_AH,
            "Secondary: Negative electrode Young's modulus [Pa]": 1.0e10,
            "Secondary: Negative electrode LAM constant proportional term [s-1]": SI_LAM_PROP_BASELINE,
            "Secondary: Negative electrode LAM constant exponential term": SI_LAM_EXP,
            "Secondary: Negative electrode critical stress [Pa]": SI_CRIT_STRESS,
            "Secondary: Negative electrode cracking rate": BASE_CRACK_RATE * SI_CRACK_RATE_MULT,
            "Primary: Negative electrode LAM constant exponential term": 2.0,
            "Primary: Negative electrode critical stress [Pa]": GR_CRIT_STRESS,
            "Primary: Negative electrode LAM constant proportional term [s-1]": GR_LAM_PROP_BASELINE,
            "Primary: Negative electrode cracking rate": BASE_CRACK_RATE * GR_CRACK_RATE_MULT,
            "Negative electrode porosity floor": NEG_POROSITY_FLOOR,
            "Primary: SEI reaction exponent cap": EXPONENT_MAX_SEI,
            "Secondary: SEI reaction exponent cap": EXPONENT_MAX_SEI,
            "Negative electrode pore buffering transition width": WIDTH_BASELINE,
            "Negative electrode transmitted fraction plateau": F0_BASELINE,
            "Primary: SEI kinetic rate constant [m.s-1]": K_SEI_DEFAULT * K_SEI_ANCHOR_MULT,
            "Secondary: SEI kinetic rate constant [m.s-1]": K_SEI_DEFAULT * K_SEI_ANCHOR_MULT,
        },
        check_already_exists=False,
    )
    return param


def ageing_leg(cyc):
    """(cap, discharge_end_thr, charge_end_thr) for the C/3 ageing-rate discharge leg."""
    cap, discharge_end_thr, charge_end_thr = 0.0, None, None
    for step in getattr(cyc, "steps", []):
        try:
            I = step["Current [A]"].entries
            q = step["Discharge capacity [A.h]"].entries
            thr_step = step["Throughput capacity [A.h]"].entries
        except (KeyError, TypeError, AttributeError):
            continue
        if I.size < 2:
            continue
        mean_I = np.mean(I)
        if mean_I > 1.0 and q.size >= 2:
            c = float(q[-1] - q[0])
            if c > cap:
                cap, discharge_end_thr = c, float(thr_step[-1])
        elif mean_I < -1e-3:
            charge_end_thr = float(thr_step[-1])
    return cap, discharge_end_thr, charge_end_thr


def rpt_leg(cyc):
    """(cap, thr, Q_array, V_array) for the C/20 RPT discharge leg, or (0, None, None,
    None) if this cycle has no RPT leg. Q_array is discharge capacity FROM THE START of
    the RPT discharge (0 at 4.2 V) -- the raw pseudo-OCV curve."""
    for step in getattr(cyc, "steps", []):
        try:
            I = step["Current [A]"].entries
            q = step["Discharge capacity [A.h]"].entries
            V = step["Terminal voltage [V]"].entries
            thr_step = step["Throughput capacity [A.h]"].entries
        except (KeyError, TypeError, AttributeError):
            continue
        if I.size < 2:
            continue
        mean_I = np.mean(I)
        if 0.01 < mean_I <= 1.0 and q.size >= 2:
            cap = float(q[-1] - q[0])
            Q_arr = q - q[0]
            return cap, float(thr_step[-1]), Q_arr, np.array(V)
    return 0.0, None, None, None


def current_soh(sol):
    caps = []
    for cyc in sol.cycles:
        cap, discharge_end_thr, _ = ageing_leg(cyc)
        if cap > 0 and discharge_end_thr is not None:
            caps.append(cap)
    if len(caps) < 2:
        return None
    return caps[-1] / caps[0]


print(f"\n=== DMA baseline run: 'ec reaction limited' anchor, RPT={RPT_RATE} every "
      f"{RPT_INTERVAL} cycles ===")
param = build_param()
model = pybamm.lithium_ion.DFN(MODEL_OPTIONS_BASE)


def make_sim(experiment):
    return pybamm.Simulation(
        model, parameter_values=param, experiment=experiment,
        solver=solver, var_pts=VAR_PTS,
    )


sim = make_sim(formation_exp)
sol = sim.solve(initial_soc=1.0)
print("formation solved.")

total_cycles = 0
stop_reason = "reached MAX_TOTAL_CYCLES safety cap"
while total_cycles < MAX_TOTAL_CYCLES:
    sim = make_sim(ageing_batch_exp)
    try:
        sol = sim.solve(starting_solution=sol)
    except Exception as exc:  # noqa: BLE001
        print(f"solver failure after {total_cycles} cycles: {type(exc).__name__}: {str(exc)[:300]}")
        stop_reason = "solver_failure"
        break
    total_cycles += BATCH_SIZE
    soh_now = current_soh(sol)
    print(f"{total_cycles} cycles solved. SoH={'n/a' if soh_now is None else f'{soh_now:.3f}'}")
    if soh_now is not None and soh_now <= SOH_TARGET:
        stop_reason = f"reached SOH_TARGET ({SOH_TARGET:.0%})"
        break

print(f"stop reason: {stop_reason}")

# ---------------------------------------------------------------------------
# extract full-run series
# ---------------------------------------------------------------------------
Qt_full = sol["Throughput capacity [A.h]"].entries
tc_cell_full = sol["Cell thickness change [m]"].entries
k_full = sol["Negative electrode transfer ratio k"].entries
LLI_full = sol["Loss of lithium inventory [%]"].entries
LAM_gr_full = sol["Loss of active material in primary phase in negative electrode [%]"].entries
LAM_si_full = sol["Loss of active material in secondary phase in negative electrode [%]"].entries
LAM_pos_full = sol["Loss of active material in positive electrode [%]"].entries

age_thr, age_cap, age_discharge_end, age_charge_end = [], [], [], []
rpt_thr, rpt_cap, rpt_Q_list, rpt_V_list = [], [], [], []
for cyc in sol.cycles:
    cap, discharge_end_thr, charge_end_thr = ageing_leg(cyc)
    if cap > 0 and discharge_end_thr is not None and charge_end_thr is not None:
        age_thr.append(discharge_end_thr)
        age_cap.append(cap)
        age_discharge_end.append(discharge_end_thr)
        age_charge_end.append(charge_end_thr)
        continue
    r_cap, r_thr, Q_arr, V_arr = rpt_leg(cyc)
    if r_cap > 0 and r_thr is not None:
        rpt_thr.append(r_thr)
        rpt_cap.append(r_cap)
        rpt_Q_list.append(Q_arr)
        rpt_V_list.append(V_arr)

age_thr = np.array(age_thr)
age_cap = np.array(age_cap)
age_discharge_end = np.array(age_discharge_end)
age_charge_end = np.array(age_charge_end)
rpt_thr = np.array(rpt_thr)
rpt_cap = np.array(rpt_cap)

tc_charge_end = np.interp(age_charge_end, Qt_full, tc_cell_full)
tc_discharge_end = np.interp(age_discharge_end, Qt_full, tc_cell_full)
age_amplitude = tc_charge_end - tc_discharge_end
age_k = np.interp(age_charge_end, Qt_full, k_full)
age_LLI = np.interp(age_discharge_end, Qt_full, LLI_full)
age_LAM_gr = np.interp(age_discharge_end, Qt_full, LAM_gr_full)
age_LAM_si = np.interp(age_discharge_end, Qt_full, LAM_si_full)
age_LAM_pos = np.interp(age_discharge_end, Qt_full, LAM_pos_full)

rpt_LLI = np.interp(rpt_thr, Qt_full, LLI_full)
rpt_LAM_gr = np.interp(rpt_thr, Qt_full, LAM_gr_full)
rpt_LAM_si = np.interp(rpt_thr, Qt_full, LAM_si_full)
rpt_LAM_pos = np.interp(rpt_thr, Qt_full, LAM_pos_full)

nominal_cap = param["Nominal cell capacity [A.h]"]

print(f"\nAgeing points: {age_thr.size}   RPT points: {rpt_thr.size}")
print(f"Final: LLI={age_LLI[-1]:.2f}%  LAM_gr={age_LAM_gr[-1]:.2f}%  "
      f"LAM_si={age_LAM_si[-1]:.2f}%  LAM_pos={age_LAM_pos[-1]:.2f}%  "
      f"final SoH={age_cap[-1]/age_cap[0]:.3f}")

outpath = os.path.join(SCRIPT_DIR, "dma_baseline_run_data.npz")
np.savez(
    outpath,
    nominal_cap=nominal_cap,
    age_thr=age_thr, age_cap=age_cap, age_amplitude=age_amplitude, age_k=age_k,
    age_LLI=age_LLI, age_LAM_gr=age_LAM_gr, age_LAM_si=age_LAM_si, age_LAM_pos=age_LAM_pos,
    rpt_thr=rpt_thr, rpt_cap=rpt_cap,
    rpt_LLI=rpt_LLI, rpt_LAM_gr=rpt_LAM_gr, rpt_LAM_si=rpt_LAM_si, rpt_LAM_pos=rpt_LAM_pos,
    rpt_Q_list=np.array(rpt_Q_list, dtype=object),
    rpt_V_list=np.array(rpt_V_list, dtype=object),
)
print(f"\nSaved: {outpath}")
