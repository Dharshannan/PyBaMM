import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import numpy as np
import pybamm

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

"""
1000-CYCLE (TIMESCALE-STRETCHED) VERSION of pore_buffering_degradation_test.py.
Same composition, degradation recipe and pore-buffering submodel as the short
(~144-cycle) test, but with all RATE constants scaled by
1/TIMESCALE_STRETCH (TIMESCALE_STRETCH=7.0, same treatment as
test_expansion/winner_deepdive_si20pct_expansion_1000cyc.py) to stretch the
same knee shape out over ~1000 cycles instead of ~130-150. Thresholds/
exponents (crit_stress, LAM exponents, porosity floor, exponent_max_sei) are
left unchanged since they set the knee's SHAPE, not its timing.

Runs the SAME two full degradation simulations as the short test (pore
buffering off/on, transition=PORE_BUFFERING_TRANSITION), then produces, for
EACH run separately:
  - the LAM/LLI/porosity-evolution diagnostic grid (ported from
    winner_deepdive_si20pct_expansion_1000cyc.py's 2x3 grid) -- ONE figure
    per run (unbuffered, buffered), not combined;
  - the electrode-contribution (NE/PE/cell thickness-change vs capacity)
    3-panel plot -- ONE figure per run this time (previously only the
    buffered run had this).

And, shared across both runs:
  - the SoH/k(both definitions)/expansion-amplitude 3-panel (buffered run
    only -- k is only meaningful with buffering on);
  - the porosity "breathing" diagnostic 2x2 (buffered run only -- shows the
    buffering reversible component, meaningless without it);
  - the side-by-side unbuffered-vs-buffered comparison 2x2 (the "combined"
    plot the user explicitly asked to keep).

See CHANGES.md items 5-10 for the pore-buffering submodel's full history.
"""

pybamm.set_logging_level("NOTICE")

MODEL_OPTIONS_BASE = {
    "particle phases": ("2", "1"),
    "open-circuit potential": (("single", "current sigmoid"), "single"),
    "SEI": "reaction limited",
    "SEI porosity change": "true",
    "particle mechanics": ("swelling and cracking", "swelling only"),
    "SEI on cracks": "true",
    "loss of active material": "stress-driven",
}

CAPACITY_PROBE_OPTIONS = {
    "particle phases": ("2", "1"),
    "open-circuit potential": (("single", "current sigmoid"), "single"),
}

PORE_BUFFERING_TRANSITION = "physical"

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

# ---------------------------------------------------------------------------
# degradation recipe: same as the short test EXCEPT rate constants scaled by
# 1/TIMESCALE_STRETCH (thresholds/exponents/floor unchanged) -- see
# test_expansion/winner_deepdive_si20pct_expansion_1000cyc.py, same pattern.
# ---------------------------------------------------------------------------
TIMESCALE_STRETCH = 7.0

SI_CRIT_STRESS = 2.2e8       # unchanged (threshold)
SI_MULT = 18.0 / TIMESCALE_STRETCH
SI_LAM_PROP = 7.5e-7 / TIMESCALE_STRETCH
SI_LAM_EXP = 2.5              # unchanged (shape exponent)
GR_DIV = 4.0 * TIMESCALE_STRETCH
GR_CRIT_STRESS = 3.0e7        # unchanged (threshold)
GR_LAM_PROP = 4.5e-6 / TIMESCALE_STRETCH
SI_CRACK_RATE_MULT = 0.1 / TIMESCALE_STRETCH
GR_CRACK_RATE_MULT = 0.1 / TIMESCALE_STRETCH
BASE_CRACK_RATE = 3.9e-20
NEG_POROSITY_FLOOR = 0.01     # unchanged (shape)
EXPONENT_MAX_SEI = 10.0       # unchanged (shape)

TARGET_BAND = (82.0, 93.0)
SOH_TERMINATION_PERCENT = 50.0

# Env-var overrides so a fast smoke test can be run without editing the
# file: PB1000_BATCH_SIZE / PB1000_RPT_INTERVAL / PB1000_MAX_CYCLES.
BATCH_SIZE = int(os.environ.get("PB1000_BATCH_SIZE", 50))
MAX_TOTAL_CYCLES = int(os.environ.get("PB1000_MAX_CYCLES", 1300))
RPT_INTERVAL = int(os.environ.get("PB1000_RPT_INTERVAL", 50))
RPT_RATE = "C/10"
assert BATCH_SIZE % RPT_INTERVAL == 0, "BATCH_SIZE must be a multiple of RPT_INTERVAL"

PARAM_UPDATES = {
    "Primary: Negative electrode active material volume fraction": GR_EPS,
    "Secondary: Negative electrode active material volume fraction": SI_EPS,
    "Secondary: Maximum concentration in negative electrode [mol.m-3]": SI_MAX_CONC_NEEDED,
    "Secondary: Initial concentration in negative electrode [mol.m-3]": scaled_si_initial_conc(SI_MAX_CONC_NEEDED),
    "Nominal cell capacity [A.h]": NOMINAL_CAP_AH,
    "Secondary: Negative electrode Young's modulus [Pa]": 1.0e10,
    "Secondary: SEI reaction exchange current density [A.m-2]": 1.5e-07 * SI_MULT,
    "Secondary: Negative electrode LAM constant proportional term [s-1]": SI_LAM_PROP,
    "Secondary: Negative electrode LAM constant exponential term": SI_LAM_EXP,
    "Secondary: Negative electrode critical stress [Pa]": SI_CRIT_STRESS,
    "Secondary: Negative electrode cracking rate": BASE_CRACK_RATE * SI_CRACK_RATE_MULT,
    "Primary: SEI reaction exchange current density [A.m-2]": 1.5e-07 / GR_DIV,
    "Primary: Negative electrode LAM constant exponential term": 2.0,
    "Primary: Negative electrode critical stress [Pa]": GR_CRIT_STRESS,
    "Primary: Negative electrode LAM constant proportional term [s-1]": GR_LAM_PROP,
    "Primary: Negative electrode cracking rate": BASE_CRACK_RATE * GR_CRACK_RATE_MULT,
    "Negative electrode porosity floor": NEG_POROSITY_FLOOR,
    "Primary: SEI reaction exponent cap": EXPONENT_MAX_SEI,
    "Secondary: SEI reaction exponent cap": EXPONENT_MAX_SEI,
}

solver = pybamm.IDAKLUSolver(root_tol=1e-06, atol=1e-06, rtol=1e-06)


def cycle_ageing_leg(cyc):
    cap, rate = 0.0, 0.0
    discharge_end_thr, charge_end_thr = None, None
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
                cap, rate = c, float(mean_I)
                discharge_end_thr = float(thr_step[-1])
        elif mean_I < -1e-3:
            charge_end_thr = float(thr_step[-1])
    return cap, rate, discharge_end_thr, charge_end_thr


def rpt_leg(cyc):
    cap, rate = 0.0, 0.0
    discharge_end_thr, charge_end_thr = None, None
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
        if 0.1 < mean_I <= 1.0 and q.size >= 2:
            c = float(q[-1] - q[0])
            if c > cap:
                cap, rate = c, mean_I
                discharge_end_thr = float(thr_step[-1])
        elif mean_I < -1e-3:
            charge_end_thr = float(thr_step[-1])
    return cap, rate, discharge_end_thr, charge_end_thr


def find_discharge_step(cyc):
    """Return the raw C/3 discharge Step object for this cycle (or None) --
    used for the voltage-curve diagnostic panel."""
    for step in getattr(cyc, "steps", []):
        try:
            I = step["Current [A]"].entries
            q = step["Discharge capacity [A.h]"].entries
        except (KeyError, TypeError, AttributeError):
            continue
        if q.size >= 2 and np.mean(I) > 1.0:
            return step
    return None


def soh_from_solution(sol_, cap0=None):
    cyc, cap = [], []
    for i, c in enumerate(sol_.cycles):
        cp, rate, _, _ = cycle_ageing_leg(c)
        if cp > 0 and rate > 1.0:
            cyc.append(i)
            cap.append(cp)
    cyc = np.array(cyc)
    cap = np.array(cap)
    if cap.size == 0:
        return cyc, np.array([])
    c0 = cap[0] if cap0 is None else cap0
    return cyc, 100.0 * cap / c0


def find_knee_bisector(thr, soh_, pre_frac=0.3, post_frac=0.15):
    n = len(soh_)
    if n < 10:
        return None
    n_pre = max(5, int(n * pre_frac))
    n_post = max(5, int(n * post_frac))
    m1, b1 = np.polyfit(thr[:n_pre], soh_[:n_pre], 1)
    m2, b2 = np.polyfit(thr[-n_post:], soh_[-n_post:], 1)
    if abs(m1 - m2) < 1e-12:
        return None
    x_int = (b2 - b1) / (m1 - m2)
    y_int = m1 * x_int + b1
    return x_int, y_int, m1, m2


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


def run_full_degradation(pore_buffering, transition=PORE_BUFFERING_TRANSITION):
    options = dict(MODEL_OPTIONS_BASE)
    options["pore buffering"] = pore_buffering
    if pore_buffering == "true":
        options["pore buffering transition"] = transition
    model = pybamm.lithium_ion.DFN(options)
    param = pybamm.ParameterValues("si_gr_expansion")
    param.update(PARAM_UPDATES, check_already_exists=False)

    def make_sim(experiment):
        return pybamm.Simulation(
            model, parameter_values=param, experiment=experiment,
            solver=solver, var_pts=VAR_PTS,
        )

    tag = f"[pore buffering={pore_buffering}]"
    sim = make_sim(formation_exp)
    last_sol = sim.solve(initial_soc=1.0)
    print(f"{tag} Formation cycle solved.")

    total_cycles_requested = 0
    stop_reason = "reached MAX_TOTAL_CYCLES"
    batch_num = 0
    while total_cycles_requested < MAX_TOTAL_CYCLES:
        batch_num += 1
        sim = make_sim(ageing_batch_exp)
        try:
            new_sol = sim.solve(starting_solution=last_sol)
        except Exception as exc:  # noqa: BLE001
            print(f"{tag} [BATCH {batch_num}] solver failure after "
                  f"{total_cycles_requested} successfully completed ageing cycles: "
                  f"{type(exc).__name__}: {str(exc)[:200]}")
            stop_reason = f"solver_failure_batch_{batch_num}"
            break

        last_sol = new_sol
        total_cycles_requested += BATCH_SIZE

        cyc, soh = soh_from_solution(last_sol)
        if soh.size:
            print(f"{tag} [BATCH {batch_num}] {total_cycles_requested} ageing cycles "
                  f"requested, {soh.size} completed, SoH={soh[-1]:.2f}%")
            if soh[-1] <= SOH_TERMINATION_PERCENT:
                stop_reason = f"reached_{SOH_TERMINATION_PERCENT:.0f}pct_floor"
                print(f"{tag} Reached {SOH_TERMINATION_PERCENT:.0f}% SoH floor -- stopping.")
                break
        else:
            print(f"{tag} [BATCH {batch_num}] {total_cycles_requested} ageing cycles "
                  f"requested, no completed ageing legs found yet")

    print(f"{tag} Stop reason: {stop_reason}")
    print(f"{tag} sol.termination:", getattr(last_sol, "termination", "<unknown>"))
    return last_sol


def extract_cycle_data(sol):
    """Per-cycle SoH, within-cycle cell-level expansion amplitude, RPT
    points and knee detection -- everything meaningful for BOTH the
    buffered and unbuffered runs. Also returns the raw charge/discharge-end
    throughput arrays (including RPT ones) so per-electrode amplitude
    extraction (extract_electrode_thickness, below) can reuse them."""
    Qt_full = sol["Throughput capacity [A.h]"].entries
    tc_cell_full = sol["Cell thickness change [m]"].entries

    age_cyc, age_cap, age_thr = [], [], []
    age_discharge_end_thr, age_charge_end_thr = [], []
    rpt_cyc, rpt_cap_list, rpt_thr = [], [], []
    rpt_discharge_end_thr, rpt_charge_end_thr = [], []
    for i, cyc in enumerate(sol.cycles):
        cap, rate, discharge_end_thr, charge_end_thr = cycle_ageing_leg(cyc)
        if cap > 0 and discharge_end_thr is not None and charge_end_thr is not None:
            thr = float(cyc["Throughput capacity [A.h]"].entries[-1])
            age_cyc.append(i)
            age_cap.append(cap)
            age_thr.append(thr)
            age_discharge_end_thr.append(discharge_end_thr)
            age_charge_end_thr.append(charge_end_thr)
            continue
        rpt_cap, rpt_rate, rpt_discharge_end, rpt_charge_end = rpt_leg(cyc)
        if rpt_cap > 0 and rpt_discharge_end is not None and rpt_charge_end is not None:
            rpt_cyc.append(i)
            rpt_cap_list.append(rpt_cap)
            rpt_thr.append(float(cyc["Throughput capacity [A.h]"].entries[-1]))
            rpt_discharge_end_thr.append(rpt_discharge_end)
            rpt_charge_end_thr.append(rpt_charge_end)

    age_cyc = np.array(age_cyc)
    age_cap = np.array(age_cap)
    age_thr = np.array(age_thr)
    age_discharge_end_thr = np.array(age_discharge_end_thr)
    age_charge_end_thr = np.array(age_charge_end_thr)
    rpt_cyc = np.array(rpt_cyc)
    rpt_cap_arr = np.array(rpt_cap_list)
    rpt_thr = np.array(rpt_thr)
    rpt_discharge_end_thr = np.array(rpt_discharge_end_thr)
    rpt_charge_end_thr = np.array(rpt_charge_end_thr)

    soh_full = 100 * age_cap / age_cap[0]
    below = np.where(soh_full <= SOH_TERMINATION_PERCENT)[0]
    end = (below[0] + 1) if below.size else len(soh_full)
    age_cyc, age_cap, age_thr, soh = age_cyc[:end], age_cap[:end], age_thr[:end], soh_full[:end]
    age_discharge_end_thr = age_discharge_end_thr[:end]
    age_charge_end_thr = age_charge_end_thr[:end]

    rpt_keep = rpt_thr <= age_thr[-1]
    rpt_cyc, rpt_cap_arr, rpt_thr = rpt_cyc[rpt_keep], rpt_cap_arr[rpt_keep], rpt_thr[rpt_keep]
    rpt_discharge_end_thr = rpt_discharge_end_thr[rpt_keep]
    rpt_charge_end_thr = rpt_charge_end_thr[rpt_keep]
    rpt_soh = 100 * rpt_cap_arr / age_cap[0] if rpt_cap_arr.size else np.array([])

    tc_cell_charge_end = np.interp(age_charge_end_thr, Qt_full, tc_cell_full)
    tc_cell_discharge_end = np.interp(age_discharge_end_thr, Qt_full, tc_cell_full)
    cell_amplitude = tc_cell_charge_end - tc_cell_discharge_end

    rpt_cell_amplitude = np.array([])
    if rpt_charge_end_thr.size:
        rpt_tc_charge_end = np.interp(rpt_charge_end_thr, Qt_full, tc_cell_full)
        rpt_tc_discharge_end = np.interp(rpt_discharge_end_thr, Qt_full, tc_cell_full)
        rpt_cell_amplitude = rpt_tc_charge_end - rpt_tc_discharge_end

    knee = find_knee_bisector(age_thr, soh)
    knee_thr = knee_soh = sharpness = None
    if knee is not None:
        knee_thr, knee_soh, m1, m2 = knee
        sharpness = abs(m2) / (abs(m1) + 1e-12)

    return dict(
        sol=sol, Qt_full=Qt_full, tc_cell_full=tc_cell_full,
        age_cyc=age_cyc, age_cap=age_cap, age_thr=age_thr, soh=soh,
        age_discharge_end_thr=age_discharge_end_thr, age_charge_end_thr=age_charge_end_thr,
        rpt_cyc=rpt_cyc, rpt_cap_arr=rpt_cap_arr, rpt_thr=rpt_thr, rpt_soh=rpt_soh,
        rpt_discharge_end_thr=rpt_discharge_end_thr, rpt_charge_end_thr=rpt_charge_end_thr,
        cell_amplitude=cell_amplitude, rpt_cell_amplitude=rpt_cell_amplitude,
        knee_thr=knee_thr, knee_soh=knee_soh, sharpness=sharpness,
    )


def extract_electrode_thickness(sol, data):
    """NE/PE thickness-change within-cycle amplitude -- these variables
    exist regardless of pore buffering (thickness change is a core
    particle-mechanics output), so this works identically for the
    unbuffered and buffered runs."""
    Qt_full = data["Qt_full"]
    tc_neg_full = sol["Negative electrode thickness change [m]"].entries
    tc_pos_full = sol["Positive electrode thickness change [m]"].entries

    tc_neg_charge_end = np.interp(data["age_charge_end_thr"], Qt_full, tc_neg_full)
    tc_neg_discharge_end = np.interp(data["age_discharge_end_thr"], Qt_full, tc_neg_full)
    neg_amplitude = tc_neg_charge_end - tc_neg_discharge_end

    tc_pos_charge_end = np.interp(data["age_charge_end_thr"], Qt_full, tc_pos_full)
    tc_pos_discharge_end = np.interp(data["age_discharge_end_thr"], Qt_full, tc_pos_full)
    pos_amplitude = tc_pos_charge_end - tc_pos_discharge_end

    return dict(
        tc_neg_full=tc_neg_full, tc_pos_full=tc_pos_full,
        neg_amplitude=neg_amplitude, pos_amplitude=pos_amplitude,
    )


def extract_degradation_diagnostics(sol, data):
    """LAM/LLI/porosity-evolution diagnostics, ported from
    winner_deepdive_si20pct_expansion_1000cyc.py, generalised to run on
    either the unbuffered or buffered solution."""
    Qt_full = data["Qt_full"]
    age_thr = data["age_thr"]
    cutoff_time_idx = np.searchsorted(Qt_full, age_thr[-1], side="right")

    def full_res_var(name):
        return sol[name].entries[:cutoff_time_idx]

    Qt = Qt_full[:cutoff_time_idx]

    LLI_builtin = full_res_var("Loss of lithium inventory [%]")
    Q_side_reactions_total = full_res_var("Total capacity lost to side reactions [A.h]")
    LLI_reconstructed = 100 * Q_side_reactions_total / NOMINAL_CAP_AH
    pct_diff = 100 * abs(LLI_builtin[-1] - LLI_reconstructed[-1]) / max(LLI_builtin[-1], 1e-9)
    if pct_diff < 2.0:
        LLI_to_plot = LLI_builtin
        LLI_label = "LLI (built-in, verified to include SEI-on-cracks)"
    else:
        LLI_to_plot = LLI_reconstructed
        LLI_label = "LLI (reconstructed: SEI + SEI-on-cracks, all domains)"

    LAM_neg = full_res_var("Loss of active material in negative electrode [%]")
    LAM_pos = full_res_var("Loss of active material in positive electrode [%]")
    LAM_gr = full_res_var("Loss of active material in primary phase in negative electrode [%]")
    LAM_si = full_res_var("Loss of active material in secondary phase in negative electrode [%]")

    Q_neg_prim_sei = full_res_var("Loss of capacity to negative primary SEI [A.h]")
    Q_neg_sec_sei = full_res_var("Loss of capacity to negative secondary SEI [A.h]")
    Q_neg_prim_sei_cr = full_res_var("Loss of capacity to negative primary SEI on cracks [A.h]")
    Q_neg_sec_sei_cr = full_res_var("Loss of capacity to negative secondary SEI on cracks [A.h]")
    LLI_neg_prim_sei = 100 * Q_neg_prim_sei / NOMINAL_CAP_AH
    LLI_neg_sec_sei = 100 * Q_neg_sec_sei / NOMINAL_CAP_AH
    LLI_neg_prim_sei_cr = 100 * Q_neg_prim_sei_cr / NOMINAL_CAP_AH
    LLI_neg_sec_sei_cr = 100 * Q_neg_sec_sei_cr / NOMINAL_CAP_AH

    eps_n_full = sol["Negative electrode porosity"].entries[:, :cutoff_time_idx]
    eps_n_avg = full_res_var("X-averaged negative electrode porosity")
    eps_n_cc = eps_n_full[0, :]
    eps_n_sep_side = eps_n_full[-1, :]

    cycle_voltage_traces = {}
    for i, cyc in enumerate(sol.cycles):
        step = find_discharge_step(cyc)
        if step is None:
            continue
        try:
            q = step["Discharge capacity [A.h]"].entries
            V = step["Voltage [V]"].entries
        except (KeyError, TypeError, AttributeError):
            continue
        if q.size >= 2:
            cycle_voltage_traces[i] = (q - q[0], V)

    return dict(
        Qt=Qt, LLI_to_plot=LLI_to_plot, LLI_label=LLI_label,
        LAM_neg=LAM_neg, LAM_pos=LAM_pos, LAM_gr=LAM_gr, LAM_si=LAM_si,
        LLI_neg_prim_sei=LLI_neg_prim_sei, LLI_neg_sec_sei=LLI_neg_sec_sei,
        LLI_neg_prim_sei_cr=LLI_neg_prim_sei_cr, LLI_neg_sec_sei_cr=LLI_neg_sec_sei_cr,
        eps_n_cc=eps_n_cc, eps_n_avg=eps_n_avg, eps_n_sep_side=eps_n_sep_side,
        cycle_voltage_traces=cycle_voltage_traces,
    )


def plot_diagnostic_grid(data, ddata, label, tag, outname):
    """2x3 LAM/LLI/porosity-evolution diagnostic grid for ONE run (either
    unbuffered or buffered) -- ported from
    winner_deepdive_si20pct_expansion_1000cyc.py."""
    age_thr, soh = data["age_thr"], data["soh"]
    rpt_thr, rpt_soh = data["rpt_thr"], data["rpt_soh"]
    knee_thr, knee_soh = data["knee_thr"], data["knee_soh"]
    Qt = ddata["Qt"]

    fig, ax = plt.subplots(2, 3, figsize=(19, 10))

    ax[0, 0].plot(age_thr, soh, "-", color="steelblue", lw=1.5, label="SoH")
    ax[0, 0].axhspan(*TARGET_BAND, color="orange", alpha=0.15,
                      label=f"target knee band ({TARGET_BAND[0]:.0f}-{TARGET_BAND[1]:.0f}%)")
    if knee_thr is not None:
        ax[0, 0].axvline(knee_thr, color="k", ls=":", lw=1)
        ax[0, 0].plot(knee_thr, knee_soh, "ko", ms=6, label=f"knee: {knee_soh:.1f}%")
    ax[0, 0].axhline(SOH_TERMINATION_PERCENT, color="red", ls="--", lw=0.8,
                      label=f"{SOH_TERMINATION_PERCENT:.0f}% floor")
    ax[0, 0].set_ylabel("SoH [%]")
    ax[0, 0].set_title(f"{label}: SoH vs throughput")
    ax[0, 0].legend(loc="lower left", fontsize=8)
    ax[0, 0].set_xlabel("Throughput capacity [A.h]")

    ax[0, 1].plot(Qt, ddata["LLI_to_plot"], label=ddata["LLI_label"], color="crimson", lw=1.5)
    ax[0, 1].plot(Qt, ddata["LAM_neg"], label="LAM negative electrode", color="darkorange", lw=1.5)
    ax[0, 1].plot(Qt, ddata["LAM_pos"], label="LAM positive electrode", color="seagreen", lw=1.5)
    ax[0, 1].set_ylabel("Degradation mode [%]")
    ax[0, 1].set_xlabel("Throughput capacity [A.h]")
    ax[0, 1].set_title("Overview: LLI and LAM")
    ax[0, 1].legend(fontsize=8)

    ax[0, 2].plot(Qt, ddata["LAM_gr"], label="LAM Gr (primary)", color="saddlebrown", lw=1.3)
    ax[0, 2].plot(Qt, ddata["LAM_si"], label="LAM Si (secondary)", color="darkorange", lw=1.3)
    ax[0, 2].plot(Qt, ddata["LLI_neg_prim_sei"], label="Gr interfacial SEI (LLI contrib.)",
                  color="teal", lw=1.3, ls="--")
    ax[0, 2].plot(Qt, ddata["LLI_neg_sec_sei"], label="Si interfacial SEI (LLI contrib.)",
                  color="purple", lw=1.3, ls="--")
    ax[0, 2].plot(Qt, ddata["LLI_neg_sec_sei_cr"], label="Si SEI-on-cracks (LLI contrib.)",
                  color="crimson", lw=1.5, ls="-")
    ax[0, 2].plot(Qt, ddata["LLI_neg_prim_sei_cr"], label="Gr SEI-on-cracks (rate suppressed 10x)",
                  color="gray", lw=1.0, ls=":")
    ax[0, 2].set_ylabel("Contribution [%]")
    ax[0, 2].set_xlabel("Throughput capacity [A.h]")
    ax[0, 2].set_title("Negative electrode degradation contributions")
    ax[0, 2].legend(fontsize=7, ncol=2)

    all_v_cycles = sorted(ddata["cycle_voltage_traces"].keys())
    stride = max(1, len(all_v_cycles) // 40)
    selected = all_v_cycles[::stride]
    if selected:
        cmap = cm.viridis
        norm = Normalize(vmin=min(selected), vmax=max(selected))
        for i in selected:
            q_rel, V = ddata["cycle_voltage_traces"][i]
            ax[1, 0].plot(q_rel, V, color=cmap(norm(i)), lw=0.9)
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax[1, 0])
        cbar.set_label("Cycle number")
    ax[1, 0].set_xlabel("Discharge capacity [A.h]")
    ax[1, 0].set_ylabel("Voltage [V]")
    ax[1, 0].set_title(f"Discharge voltage curves ({len(selected)} shown)")

    ax[1, 1].scatter(age_thr, soh, s=14, color="steelblue", label="C/3 ageing SoH", zorder=3)
    if rpt_soh.size:
        ax[1, 1].scatter(rpt_thr, rpt_soh, s=28, color="crimson", marker="^",
                          label=f"RPT ({RPT_RATE}) SoH", zorder=4)
    ax[1, 1].axhspan(*TARGET_BAND, color="orange", alpha=0.12)
    ax[1, 1].axhline(SOH_TERMINATION_PERCENT, color="red", ls="--", lw=0.8)
    ax[1, 1].set_xlabel("Throughput capacity [A.h]")
    ax[1, 1].set_ylabel("SoH [%]")
    ax[1, 1].set_title("SoH: C/3 ageing vs RPT capacity (scatter)")
    ax[1, 1].legend(fontsize=8)

    ax[1, 2].plot(Qt, ddata["eps_n_cc"], label="CC-side (x=0)", color="steelblue", lw=1.3)
    ax[1, 2].plot(Qt, ddata["eps_n_avg"], label="X-averaged", color="darkorange", lw=1.3)
    ax[1, 2].plot(Qt, ddata["eps_n_sep_side"], label="Separator-side", color="crimson", lw=1.3)
    ax[1, 2].axhline(NEG_POROSITY_FLOOR, color="gray", ls=":", lw=1, label=f"floor={NEG_POROSITY_FLOOR}")
    ax[1, 2].set_xlabel("Throughput capacity [A.h]")
    ax[1, 2].set_ylabel("Negative electrode porosity")
    ax[1, 2].set_title("Porosity evolution (the knee's underlying mechanism)")
    ax[1, 2].legend(fontsize=8)

    for row in ax:
        for a in row:
            a.grid(alpha=0.3)
    fig.suptitle(f"{label} -- LAM/LLI/porosity diagnostic grid (si_gr_expansion, {tag})")
    plt.tight_layout()
    outpath = os.path.join(SCRIPT_DIR, outname)
    plt.savefig(outpath, dpi=140, bbox_inches="tight")
    print(f"Saved: {outpath}")


def plot_volume_vs_capacity(data, edata, label, tag, outname, ne_label, pe_label, cell_label, cell_color):
    """3-panel cell-level amplitude / NE+PE+cell contributions / capacity
    fade, for ONE run (unbuffered or buffered)."""
    age_thr = data["age_thr"]
    knee_thr, knee_soh = data["knee_thr"], data["knee_soh"]
    rpt_thr = data["rpt_thr"]

    fig_vc, ax_vc = plt.subplots(3, 1, figsize=(9, 12), sharex=True)

    ax_vc[0].plot(age_thr, data["cell_amplitude"] * 1e6, "o-", ms=3, color=cell_color, label="C/3 ageing")
    if data["rpt_cell_amplitude"].size:
        ax_vc[0].scatter(rpt_thr, data["rpt_cell_amplitude"] * 1e6, marker="^", color="tab:orange",
                          s=30, label=f"RPT ({RPT_RATE})", zorder=5)
    if knee_thr is not None:
        ax_vc[0].axvline(knee_thr, color="red", ls="--", alpha=0.6, label=f"knee ({knee_soh:.1f}%)")
    ax_vc[0].set_ylabel("Cell-level within-cycle\nexpansion amplitude [um]")
    ax_vc[0].set_title(f"{label}: electrode contributions, cell thickness change and "
                        f"capacity fade\n(si_gr_expansion, {tag})")
    ax_vc[0].legend(fontsize=8)
    ax_vc[0].grid(alpha=0.3)

    ax_vc[1].plot(age_thr, edata["neg_amplitude"] * 1e6, "s-", ms=3, color="tab:brown", label=ne_label)
    ax_vc[1].plot(age_thr, edata["pos_amplitude"] * 1e6, "d-", ms=3, color="tab:cyan", label=pe_label)
    ax_vc[1].plot(age_thr, data["cell_amplitude"] * 1e6, "o-", ms=3, color=cell_color, label=cell_label)
    if knee_thr is not None:
        ax_vc[1].axvline(knee_thr, color="red", ls="--", alpha=0.6)
    ax_vc[1].axhline(0, color="gray", lw=0.8, alpha=0.6)
    ax_vc[1].set_ylabel("Within-cycle expansion\namplitude by electrode [um]")
    ax_vc[1].legend(fontsize=8)
    ax_vc[1].grid(alpha=0.3)

    ax_vc[2].plot(age_thr, data["age_cap"], "o-", ms=3, color="tab:blue", label="C/3 ageing")
    if data["rpt_cap_arr"].size:
        ax_vc[2].scatter(rpt_thr, data["rpt_cap_arr"], marker="^", color="tab:orange",
                          s=30, label=f"RPT ({RPT_RATE})", zorder=5)
    if knee_thr is not None:
        ax_vc[2].axvline(knee_thr, color="red", ls="--", alpha=0.6)
    ax_vc[2].set_ylabel("Discharge capacity [A.h]")
    ax_vc[2].set_xlabel("Throughput capacity [A.h]")
    ax_vc[2].legend(fontsize=8)
    ax_vc[2].grid(alpha=0.3)

    plt.tight_layout()
    outpath = os.path.join(SCRIPT_DIR, outname)
    plt.savefig(outpath, dpi=150)
    print(f"Saved: {outpath}")


# ---------------------------------------------------------------------------
# reference run: pore buffering OFF
# ---------------------------------------------------------------------------
print("\n=== Reference run: pore buffering OFF (1000-cycle stretched) ===")
sol_off = run_full_degradation("false")
data_off = extract_cycle_data(sol_off)
edata_off = extract_electrode_thickness(sol_off, data_off)
ddata_off = extract_degradation_diagnostics(sol_off, data_off)
print(f"Reference (unbuffered) run: {data_off['age_cap'].size} ageing cycles, "
      f"final SoH={data_off['soh'][-1]:.2f}%, "
      f"knee={data_off['knee_soh']}, sharpness={data_off['sharpness']}")

# ---------------------------------------------------------------------------
# main run: pore buffering ON
# ---------------------------------------------------------------------------
print(f"\n=== Main run: pore buffering ON, transition={PORE_BUFFERING_TRANSITION} (1000-cycle stretched) ===")
sol = run_full_degradation("true", PORE_BUFFERING_TRANSITION)
data_on = extract_cycle_data(sol)
edata_on = extract_electrode_thickness(sol, data_on)
ddata_on = extract_degradation_diagnostics(sol, data_on)

Qt_full = data_on["Qt_full"]
age_thr, soh = data_on["age_thr"], data_on["soh"]
age_discharge_end_thr, age_charge_end_thr = data_on["age_discharge_end_thr"], data_on["age_charge_end_thr"]
rpt_thr, rpt_soh = data_on["rpt_thr"], data_on["rpt_soh"]
cell_amplitude, rpt_cell_amplitude = data_on["cell_amplitude"], data_on["rpt_cell_amplitude"]
knee_thr, knee_soh, sharpness = data_on["knee_thr"], data_on["knee_soh"], data_on["sharpness"]
neg_amplitude, pos_amplitude = edata_on["neg_amplitude"], edata_on["pos_amplitude"]

print(f"\nCompleted {data_on['age_cap'].size} ageing (C/3) cycles total, {data_on['rpt_cyc'].size} RPT cycles")
print(f"Final SoH = {soh[-1]:.2f}%  (final throughput {age_thr[-1]:.1f} A.h)")
if knee_thr is not None:
    print(f"Knee (bisector): throughput={knee_thr:.1f} A.h, SoH={knee_soh:.2f}%, "
          f"sharpness={sharpness:.2f}")
else:
    print("Knee detection failed (not enough points)")

# buffered-run-only diagnostics: internal k, structural porosity, "observable" k
k_full = sol["Negative electrode transfer ratio k"].entries
eps_struct_full = sol["Negative electrode structural porosity"].entries.mean(axis=0)
eps_full = sol["Negative electrode porosity"].entries.mean(axis=0)

k_at_charge_end = np.interp(age_charge_end_thr, Qt_full, k_full)
eps_struct_at_charge_end = np.interp(age_charge_end_thr, Qt_full, eps_struct_full)

dv_solid_full = sol["Negative electrode solid volume change"].entries
_ref_param = pybamm.ParameterValues("si_gr_expansion")
n_electrodes_parallel = _ref_param["Number of electrodes connected in parallel to make a cell"]
L_neg = _ref_param["Negative electrode thickness [m]"]
tc_particle_unbuffered_full = dv_solid_full * n_electrodes_parallel * L_neg

tc_particle_charge_end = np.interp(age_charge_end_thr, Qt_full, tc_particle_unbuffered_full)
tc_particle_discharge_end = np.interp(age_discharge_end_thr, Qt_full, tc_particle_unbuffered_full)
particle_amplitude = tc_particle_charge_end - tc_particle_discharge_end

k_cell_measurable = cell_amplitude / particle_amplitude

print(f"\nInternal k (negative electrode, mechanistic) at charge-end: "
      f"t=0 -> {k_at_charge_end[0]:.4f}, final -> {k_at_charge_end[-1]:.4f} "
      f"(max over run: {k_full.max():.4f})")
print(f"Observable k = delta_cell/delta_particle at charge-end: t=0 -> "
      f"{k_cell_measurable[0]:.4f}, final -> {k_cell_measurable[-1]:.4f} "
      f"(min={k_cell_measurable.min():.4f}, max={k_cell_measurable.max():.4f})")
print(f"Cell-level within-cycle amplitude: {cell_amplitude[0]*1e6:.3f} um -> "
      f"{cell_amplitude[-1]*1e6:.3f} um (peak {np.abs(cell_amplitude).max()*1e6:.3f} um)")

# ---------------------------------------------------------------------------
# Figure 1: SoH / k (both definitions) / expansion amplitude -- buffered
# run only (k is only meaningful with buffering on)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(3, 1, figsize=(9, 11), sharex=True)

ax[0].plot(age_thr, soh, "o-", ms=3, color="tab:blue", label="C/3 ageing")
if rpt_thr.size:
    ax[0].scatter(rpt_thr, rpt_soh, marker="^", color="tab:orange", s=30, label=f"RPT ({RPT_RATE})", zorder=5)
if knee_thr is not None:
    ax[0].axvline(knee_thr, color="red", ls="--", alpha=0.6, label=f"knee ({knee_soh:.1f}%)")
ax[0].axhspan(*TARGET_BAND, color="green", alpha=0.08, label="target band")
ax[0].set_ylabel("SoH [%]")
ax[0].set_title(f"Pore-buffering 1000-cycle degradation test (si_gr_expansion, transition={PORE_BUFFERING_TRANSITION})")
ax[0].legend(fontsize=8)
ax[0].grid(alpha=0.3)

ax[1].plot(age_thr, k_at_charge_end, "o-", ms=3, color="tab:purple",
           label="Internal k (negative electrode, mechanistic)")
ax[1].plot(age_thr, k_cell_measurable, "s-", ms=3, color="tab:red",
           label="Observable k = delta_cell / delta_particle (cell-level)")
ax[1].axhline(0.7, color="gray", ls=":", alpha=0.7, label="BoL f0 = k_BoL = 0.7")
if knee_thr is not None:
    ax[1].axvline(knee_thr, color="red", ls="--", alpha=0.6)
ax[1].set_ylabel("Transfer ratio k [-]")
ax[1].legend(fontsize=8)
ax[1].grid(alpha=0.3)

ax[2].plot(age_thr, cell_amplitude * 1e6, "o-", ms=3, color="tab:green", label="Cell-level")
ax[2].plot(age_thr, neg_amplitude * 1e6, "s-", ms=3, color="tab:brown", label="Negative electrode-level")
if knee_thr is not None:
    ax[2].axvline(knee_thr, color="red", ls="--", alpha=0.6)
ax[2].set_ylabel("Within-cycle expansion\namplitude [um]")
ax[2].set_xlabel("Throughput capacity [A.h]")
ax[2].legend(fontsize=8)
ax[2].grid(alpha=0.3)

plt.tight_layout()
outpath = os.path.join(SCRIPT_DIR, f"pore_buffering_degradation_test_1000cyc_result_{PORE_BUFFERING_TRANSITION}.png")
plt.savefig(outpath, dpi=150)
print(f"\nSaved: {outpath}")

# ---------------------------------------------------------------------------
# Figure 2a/2b: electrode-contribution (NE/PE/cell) plots, ONE per run
# ---------------------------------------------------------------------------
plot_volume_vs_capacity(
    data_on, edata_on, "Buffered", f"transition={PORE_BUFFERING_TRANSITION}",
    f"pore_buffering_degradation_test_1000cyc_volume_vs_capacity_{PORE_BUFFERING_TRANSITION}.png",
    ne_label="Negative electrode (buffered)", pe_label="Positive electrode (unbuffered)",
    cell_label="Cell (negative + positive + thermal)", cell_color="tab:green",
)
plot_volume_vs_capacity(
    data_off, edata_off, "Unbuffered", "pore buffering=false",
    "pore_buffering_degradation_test_1000cyc_volume_vs_capacity_unbuffered.png",
    ne_label="Negative electrode (unbuffered)", pe_label="Positive electrode (unbuffered)",
    cell_label="Cell (negative + positive + thermal)", cell_color="tab:gray",
)

# ---------------------------------------------------------------------------
# Figure 3: porosity "breathing" diagnostic -- buffered run only. Target
# throughput points chosen as FRACTIONS of the run's own final throughput
# (rather than fixed A.h values) so this generalises across knee-position
# shifts between the short and 1000-cycle-stretched runs.
# ---------------------------------------------------------------------------
BREATH_TARGET_FRACS = [0.05, 0.5, 0.8, 0.95]
BREATH_TARGET_THR = [f * age_thr[-1] for f in BREATH_TARGET_FRACS]
breath_cyc_indices = []
for t in BREATH_TARGET_THR:
    idx = int(np.argmin(np.abs(age_thr - t)))
    breath_cyc_indices.append(int(data_on["age_cyc"][idx]))

fig_br, ax_br = plt.subplots(2, 2, figsize=(11, 8), sharey=True)
for ax_i, cyc_idx, target_t in zip(ax_br.ravel(), breath_cyc_indices, BREATH_TARGET_THR):
    cyc_sol = sol.cycles[cyc_idx]
    cyc_thr = cyc_sol["Throughput capacity [A.h]"].entries
    cyc_thr_rel = cyc_thr - cyc_thr[0]
    cyc_eps_total = cyc_sol["Negative electrode porosity"].entries.mean(axis=0)
    cyc_eps_struct = cyc_sol["Negative electrode structural porosity"].entries.mean(axis=0)
    k_here = float(np.interp(cyc_thr[-1], Qt_full, k_full))

    ax_i.plot(cyc_thr_rel, cyc_eps_total, "-", color="tab:blue", label="Total porosity")
    ax_i.plot(cyc_thr_rel, cyc_eps_struct, "--", color="gray", label="Structural porosity")
    ax_i.set_title(f"~{target_t:.0f} A.h (cycle {cyc_idx}, k={k_here:.3f})", fontsize=10)
    ax_i.grid(alpha=0.3)
    ax_i.legend(fontsize=7)

for ax_i in ax_br[-1, :]:
    ax_i.set_xlabel("Throughput within cycle [A.h]")
for ax_i in ax_br[:, 0]:
    ax_i.set_ylabel("Negative electrode\nporosity [-]")

fig_br.suptitle("Pore-buffering \"breathing\" diagnostic (1000-cycle): within-cycle porosity "
                 "oscillation vs. structural (irreversible) porosity, at different life stages")
plt.tight_layout()
outpath_br = os.path.join(SCRIPT_DIR, f"pore_buffering_degradation_test_1000cyc_breathing_{PORE_BUFFERING_TRANSITION}.png")
plt.savefig(outpath_br, dpi=150)
print(f"Saved: {outpath_br}")

# ---------------------------------------------------------------------------
# Figure 4: combined side-by-side comparison -- unbuffered (left) vs.
# buffered (right), expansion amplitude (top) and capacity fade (bottom)
# ---------------------------------------------------------------------------
fig_cmp, ax_cmp = plt.subplots(2, 2, figsize=(13, 9), sharey="row")

for col, (label, data, color) in enumerate([
    ("Unbuffered (pore buffering=false)", data_off, "tab:gray"),
    (f"Buffered (transition={PORE_BUFFERING_TRANSITION})", data_on, "tab:green"),
]):
    ax_exp = ax_cmp[0, col]
    ax_cap = ax_cmp[1, col]

    ax_exp.plot(data["age_thr"], data["cell_amplitude"] * 1e6, "o-", ms=3, color=color, label="C/3 ageing")
    if data["rpt_cell_amplitude"].size:
        ax_exp.scatter(data["rpt_thr"], data["rpt_cell_amplitude"] * 1e6, marker="^",
                        color="tab:orange", s=30, label=f"RPT ({RPT_RATE})", zorder=5)
    if data["knee_thr"] is not None:
        ax_exp.axvline(data["knee_thr"], color="red", ls="--", alpha=0.6,
                        label=f"knee ({data['knee_soh']:.1f}%)")
    ax_exp.set_title(label)
    ax_exp.legend(fontsize=8)
    ax_exp.grid(alpha=0.3)

    ax_cap.plot(data["age_thr"], data["age_cap"], "o-", ms=3, color=color, label="C/3 ageing")
    if data["rpt_cap_arr"].size:
        ax_cap.scatter(data["rpt_thr"], data["rpt_cap_arr"], marker="^",
                        color="tab:orange", s=30, label=f"RPT ({RPT_RATE})", zorder=5)
    if data["knee_thr"] is not None:
        ax_cap.axvline(data["knee_thr"], color="red", ls="--", alpha=0.6)
    ax_cap.set_xlabel("Throughput capacity [A.h]")
    ax_cap.legend(fontsize=8)
    ax_cap.grid(alpha=0.3)

ax_cmp[0, 0].set_ylabel("Cell-level within-cycle\nexpansion amplitude [um]")
ax_cmp[1, 0].set_ylabel("Discharge capacity [A.h]")
fig_cmp.suptitle("Unbuffered vs. buffered (1000-cycle): cell-level expansion and capacity fade "
                  f"(si_gr_expansion, transition={PORE_BUFFERING_TRANSITION})")
plt.tight_layout()
outpath_cmp = os.path.join(SCRIPT_DIR, f"pore_buffering_degradation_test_1000cyc_comparison_{PORE_BUFFERING_TRANSITION}.png")
plt.savefig(outpath_cmp, dpi=150)
print(f"Saved: {outpath_cmp}")

# ---------------------------------------------------------------------------
# Figure 5a/5b: LAM/LLI/porosity-evolution diagnostic grid, ONE per run
# ---------------------------------------------------------------------------
plot_diagnostic_grid(
    data_off, ddata_off, "Unbuffered", "pore buffering=false",
    "pore_buffering_degradation_test_1000cyc_diagnostics_unbuffered.png",
)
plot_diagnostic_grid(
    data_on, ddata_on, "Buffered", f"transition={PORE_BUFFERING_TRANSITION}",
    f"pore_buffering_degradation_test_1000cyc_diagnostics_{PORE_BUFFERING_TRANSITION}.png",
)

print("\nAll plots saved.")
