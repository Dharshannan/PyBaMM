import os
import matplotlib.pyplot as plt
import numpy as np
import pybamm

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

"""
SHORT DEGRADATION-CYCLE TEST for the pore-buffering (volume partition)
submodel, per si_gr_expansion_precursor/pore_buffering_implementation_plan.md
and CHANGES.md item 5. Same accepted winning composition and degradation
recipe as test_expansion/winner_deepdive_si20pct_expansion.py (si_gr_expansion
parameter set, 20% Si by AM volume, SiOx composition), run to the SAME 50%
SoH floor used throughout this project (NOT the 1000-cycle timescale-
stretched protocol) -- this is a short, ~144-cycle sanity check that pore
buffering ("pore buffering": "true") solves cleanly through a full
degradation run and produces sane k / volume-change trajectories as ageing
proceeds, not a full deep-dive.

BoL parameterisation: eps_min_transfer, eps_transfer_width, f_transmit_min
(f0). Originally set from a BoL-only placeholder (plan doc section 3), then
retuned via eps_max_transfer_sweep.py so the buffering-saturation
transition times against the capacity knee -- current defaults (0.08, 0.01,
0.7) live in si_gr_expansion.py; see CHANGES.md items 7-9 for the full
retuning/reparameterisation history (the parameter was originally named
"eps_max_transfer" -- item 9 explains why that was dropped in favour of an
explicit transition-width parameter).

Runs TWO full degradation simulations, same composition/degradation recipe,
differing only in options["pore buffering"]:
  - "false" (reference/unbuffered) -- needed for the "observable" k =
    delta_cell / delta_particle definition below.
  - "true" (main run) -- everything else (SoH, electrode contributions,
    breathing diagnostic) uses this one.

Tracks, per ageing cycle: SoH (C/3 + RPT), TWO different k definitions
(see below), and the within-cycle expansion amplitude (charge-end minus
same-cycle discharge-end, same convention as test_expansion -- cancels the
fixed BoL reference for a genuinely positive swelling signal) at cell,
negative-electrode and positive-electrode level.

Two k definitions, NOT the same quantity:
  - "Negative electrode transfer ratio k" (internal/mechanistic): the
    partition's own read-off diagnostic (Eq. 23 of the pore-buffering
    notes), k = dv_thickness/dv_solid, computed ENTIRELY within the
    negative electrode, from the buffered run alone.
  - "k = delta_cell/delta_particle" (observable, per the notes' original
    Eq. 1-2 with the "rev" = reversible/amplitude qualifier, delta_particle
    per Sec. 4.4): the CELL-level (buffered) within-cycle expansion
    amplitude, "Cell thickness change [m]", divided by delta_particle = the
    UNBUFFERED NEGATIVE ELECTRODE's own within-cycle amplitude ("cathode
    neglected", per Yin et al.) -- i.e. what you'd actually compute from a
    dilatometry measurement (delta_cell) plus the anode's own particle-level
    reference (delta_particle), NOT a second whole-cell simulation. An
    earlier version of this script used a second unbuffered simulation's
    cell-level amplitude as delta_particle, which wrongly pulled the
    positive electrode's own (unrelated, always-unbuffered, and -- per
    CHANGES.md -- implausibly large, ~79% particle volume change) swelling
    into the denominator, producing a wildly unstable ratio driven mostly by
    cathode/anode cancellation rather than by pore buffering. The unbuffered
    reference run (data_off below) is now used only for the SoH/knee cross-
    check, not for k.
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

# Switch between the two "pore buffering transition" options -- "tanh"
# (semi-empirical smooth sigmoid, the default, see CHANGES.md item 7) or
# "physical" (compliance-ratio form derived from a pore-network/stack
# stiffness balance, CHANGES.md item 8).
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

# degradation recipe -- unchanged from the accepted winner
SI_CRIT_STRESS = 2.2e8
SI_MULT = 18.0
SI_LAM_PROP = 7.5e-7
SI_LAM_EXP = 2.5
GR_DIV = 4.0
GR_CRIT_STRESS = 3.0e7
GR_LAM_PROP = 4.5e-6
SI_CRACK_RATE_MULT = 0.1
GR_CRACK_RATE_MULT = 0.1
BASE_CRACK_RATE = 3.9e-20
NEG_POROSITY_FLOOR = 0.01
EXPONENT_MAX_SEI = 10.0

TARGET_BAND = (82.0, 93.0)
SOH_TERMINATION_PERCENT = 50.0
BATCH_SIZE = 20
MAX_TOTAL_CYCLES = 300
RPT_INTERVAL = 10
RPT_RATE = "C/10"
assert BATCH_SIZE % RPT_INTERVAL == 0

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
    """Return (cap, rate, discharge_end_thr, charge_end_thr) for the C/3
    ageing legs of this cycle, or zeros/None if not found. *_end_thr are
    throughput capacity [A.h] at the end of the discharge/charge steps
    specifically, used to sample thickness-change/k at matching points for
    the within-cycle expansion amplitude."""
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
    """Same convention as cycle_ageing_leg but for the slow-rate RPT
    discharge -- also returns discharge_end_thr/charge_end_thr so RPT-cycle
    within-cycle expansion amplitude can be computed the same way."""
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
    """Build a fresh model/param pair with the given "pore buffering" (and,
    when it's "true", "pore buffering transition") option and run formation
    + batched ageing to the SOH_TERMINATION_PERCENT floor (or
    MAX_TOTAL_CYCLES). Returns the final Solution."""
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
    points and knee detection -- everything that's meaningful for BOTH the
    buffered and unbuffered runs (k / structural-porosity / electrode-level
    amplitude are buffered-run-only and computed separately in the main
    script body)."""
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
        cell_amplitude=cell_amplitude, rpt_cell_amplitude=rpt_cell_amplitude,
        knee_thr=knee_thr, knee_soh=knee_soh, sharpness=sharpness,
    )


# ---------------------------------------------------------------------------
# reference run: pore buffering OFF -- needed for the "observable"
# k = delta_cell/delta_particle definition (delta_particle here = what the
# cell-level amplitude would have been with no buffering at all)
# ---------------------------------------------------------------------------
print("\n=== Reference run: pore buffering OFF ===")
sol_off = run_full_degradation("false")
data_off = extract_cycle_data(sol_off)
print(f"Reference (unbuffered) run: {data_off['age_cap'].size} ageing cycles, "
      f"final SoH={data_off['soh'][-1]:.2f}%, "
      f"knee={data_off['knee_soh']}, sharpness={data_off['sharpness']}")

# ---------------------------------------------------------------------------
# main run: pore buffering ON
# ---------------------------------------------------------------------------
print(f"\n=== Main run: pore buffering ON, transition={PORE_BUFFERING_TRANSITION} ===")
sol = run_full_degradation("true", PORE_BUFFERING_TRANSITION)
data_on = extract_cycle_data(sol)

Qt_full = data_on["Qt_full"]
tc_cell_full = data_on["tc_cell_full"]
age_cyc, age_cap, age_thr, soh = data_on["age_cyc"], data_on["age_cap"], data_on["age_thr"], data_on["soh"]
age_discharge_end_thr, age_charge_end_thr = data_on["age_discharge_end_thr"], data_on["age_charge_end_thr"]
rpt_cyc, rpt_cap_arr, rpt_thr, rpt_soh = data_on["rpt_cyc"], data_on["rpt_cap_arr"], data_on["rpt_thr"], data_on["rpt_soh"]
cell_amplitude, rpt_cell_amplitude = data_on["cell_amplitude"], data_on["rpt_cell_amplitude"]
knee_thr, knee_soh, sharpness = data_on["knee_thr"], data_on["knee_soh"], data_on["sharpness"]

print(f"\nCompleted {age_cap.size} ageing (C/3) cycles total, {rpt_cyc.size} RPT cycles")
print(f"Final SoH = {soh[-1]:.2f}%  (final throughput {age_thr[-1]:.1f} A.h)")
if knee_thr is not None:
    print(f"Knee (bisector): throughput={knee_thr:.1f} A.h, SoH={knee_soh:.2f}%, "
          f"sharpness={sharpness:.2f}")
else:
    print("Knee detection failed (not enough points)")

# buffered-run-only diagnostics: internal k, structural porosity, per-
# electrode amplitude
tc_neg_full = sol["Negative electrode thickness change [m]"].entries
tc_pos_full = sol["Positive electrode thickness change [m]"].entries
k_full = sol["Negative electrode transfer ratio k"].entries
eps_struct_full = sol["Negative electrode structural porosity"].entries.mean(axis=0)
eps_full = sol["Negative electrode porosity"].entries.mean(axis=0)

tc_neg_charge_end = np.interp(age_charge_end_thr, Qt_full, tc_neg_full)
tc_neg_discharge_end = np.interp(age_discharge_end_thr, Qt_full, tc_neg_full)
neg_amplitude = tc_neg_charge_end - tc_neg_discharge_end

tc_pos_charge_end = np.interp(age_charge_end_thr, Qt_full, tc_pos_full)
tc_pos_discharge_end = np.interp(age_discharge_end_thr, Qt_full, tc_pos_full)
pos_amplitude = tc_pos_charge_end - tc_pos_discharge_end

k_at_charge_end = np.interp(age_charge_end_thr, Qt_full, k_full)
eps_struct_at_charge_end = np.interp(age_charge_end_thr, Qt_full, eps_struct_full)
eps_at_charge_end = np.interp(age_charge_end_thr, Qt_full, eps_full)

# "observable" k = delta_cell/delta_particle, per the notes' own definition
# (Sec. 4.4): delta_particle is the UNBUFFERED NEGATIVE ELECTRODE's own
# swelling (n*L*integral(Omega dc), "cathode neglected" -- Yin et al.), NOT
# a whole separate cell-level simulation. Using a second simulation's
# cell-level amplitude (as an earlier version of this script did) wrongly
# drags the positive electrode's own, unrelated, always-unbuffered swelling
# into the denominator -- with this composition's ~79% cathode particle
# volume change (see CHANGES.md), that produced a wildly unstable ratio
# driven mostly by cathode-vs-anode cancellation, not by buffering. The
# correct delta_particle is available from THIS SAME run: "Negative
# electrode solid volume change" (dv_solid) is exactly the pre-partition
# unbuffered swelling in volume-fraction units -- scale to thickness the
# same way the submodel does (n_electrodes_parallel * L) and take the same
# charge-end-minus-discharge-end amplitude as everything else.
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
print(f"Observable k = delta_cell/delta_particle (delta_particle = unbuffered "
      f"negative electrode's own swelling) at charge-end: t=0 -> "
      f"{k_cell_measurable[0]:.4f}, final -> {k_cell_measurable[-1]:.4f} "
      f"(min={k_cell_measurable.min():.4f}, max={k_cell_measurable.max():.4f})")
print(f"Negative electrode unbuffered particle amplitude (delta_particle): "
      f"{particle_amplitude[0]*1e6:.3f} um -> {particle_amplitude[-1]*1e6:.3f} um")
print(f"Negative electrode structural porosity: {eps_struct_at_charge_end[0]:.4f} -> "
      f"{eps_struct_at_charge_end[-1]:.4f}")
print(f"Cell-level within-cycle amplitude: {cell_amplitude[0]*1e6:.3f} um -> "
      f"{cell_amplitude[-1]*1e6:.3f} um (peak {np.abs(cell_amplitude).max()*1e6:.3f} um)")
print(f"Negative-electrode within-cycle amplitude: {neg_amplitude[0]*1e6:.3f} um -> "
      f"{neg_amplitude[-1]*1e6:.3f} um (peak {np.abs(neg_amplitude).max()*1e6:.3f} um)")
print(f"Positive-electrode within-cycle amplitude: {pos_amplitude[0]*1e6:.3f} um -> "
      f"{pos_amplitude[-1]*1e6:.3f} um (peak {np.abs(pos_amplitude).max()*1e6:.3f} um)")

# ---------------------------------------------------------------------------
# plot: SoH / k (both definitions) / expansion amplitude, shared x-axis
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(3, 1, figsize=(9, 11), sharex=True)

ax[0].plot(age_thr, soh, "o-", ms=3, color="tab:blue", label="C/3 ageing")
if rpt_thr.size:
    ax[0].scatter(rpt_thr, rpt_soh, marker="^", color="tab:orange", s=30, label=f"RPT ({RPT_RATE})", zorder=5)
if knee_thr is not None:
    ax[0].axvline(knee_thr, color="red", ls="--", alpha=0.6, label=f"knee ({knee_soh:.1f}%)")
ax[0].axhspan(*TARGET_BAND, color="green", alpha=0.08, label="target band")
ax[0].set_ylabel("SoH [%]")
ax[0].set_title(f"Pore-buffering short degradation test (si_gr_expansion, transition={PORE_BUFFERING_TRANSITION})")
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
outpath = os.path.join(SCRIPT_DIR, f"pore_buffering_degradation_test_result_{PORE_BUFFERING_TRANSITION}.png")
plt.savefig(outpath, dpi=150)
print(f"\nSaved: {outpath}")

# ---------------------------------------------------------------------------
# volume change vs. capacity: cell-level expansion + electrode contributions
# + capacity fade, shared x-axis, same convention as
# test_expansion/winner_deepdive_si20pct_expansion.py
# ---------------------------------------------------------------------------
fig_vc, ax_vc = plt.subplots(3, 1, figsize=(9, 12), sharex=True)

ax_vc[0].plot(age_thr, cell_amplitude * 1e6, "o-", ms=3, color="tab:green", label="C/3 ageing")
if rpt_cell_amplitude.size:
    ax_vc[0].scatter(rpt_thr, rpt_cell_amplitude * 1e6, marker="^", color="tab:orange",
                      s=30, label=f"RPT ({RPT_RATE})", zorder=5)
if knee_thr is not None:
    ax_vc[0].axvline(knee_thr, color="red", ls="--", alpha=0.6, label=f"knee ({knee_soh:.1f}%)")
ax_vc[0].set_ylabel("Cell-level within-cycle\nexpansion amplitude [um]")
ax_vc[0].set_title("Pore-buffering: electrode contributions, cell thickness change "
                    f"and capacity fade (si_gr_expansion, transition={PORE_BUFFERING_TRANSITION})")
ax_vc[0].legend(fontsize=8)
ax_vc[0].grid(alpha=0.3)

ax_vc[1].plot(age_thr, neg_amplitude * 1e6, "s-", ms=3, color="tab:brown", label="Negative electrode (buffered)")
ax_vc[1].plot(age_thr, pos_amplitude * 1e6, "d-", ms=3, color="tab:cyan", label="Positive electrode (unbuffered)")
ax_vc[1].plot(age_thr, cell_amplitude * 1e6, "o-", ms=3, color="tab:green", label="Cell (negative + positive + thermal)")
if knee_thr is not None:
    ax_vc[1].axvline(knee_thr, color="red", ls="--", alpha=0.6)
ax_vc[1].axhline(0, color="gray", lw=0.8, alpha=0.6)
ax_vc[1].set_ylabel("Within-cycle expansion\namplitude by electrode [um]")
ax_vc[1].legend(fontsize=8)
ax_vc[1].grid(alpha=0.3)

ax_vc[2].plot(age_thr, age_cap, "o-", ms=3, color="tab:blue", label="C/3 ageing")
if rpt_cap_arr.size:
    ax_vc[2].scatter(rpt_thr, rpt_cap_arr, marker="^", color="tab:orange",
                      s=30, label=f"RPT ({RPT_RATE})", zorder=5)
if knee_thr is not None:
    ax_vc[2].axvline(knee_thr, color="red", ls="--", alpha=0.6)
ax_vc[2].set_ylabel("Discharge capacity [A.h]")
ax_vc[2].set_xlabel("Throughput capacity [A.h]")
ax_vc[2].legend(fontsize=8)
ax_vc[2].grid(alpha=0.3)

plt.tight_layout()
outpath_vc = os.path.join(SCRIPT_DIR, f"pore_buffering_degradation_test_volume_vs_capacity_{PORE_BUFFERING_TRANSITION}.png")
plt.savefig(outpath_vc, dpi=150)
print(f"Saved: {outpath_vc}")

# ---------------------------------------------------------------------------
# porosity "breathing" diagnostic: full time-resolved (not per-cycle-
# sampled) negative electrode porosity within a handful of representative
# ageing cycles at different life stages, to directly show the pore-
# buffering reversible/breathing component. dv_buffered tracks instantaneous
# stoichiometry (Eq. 17-19 of the pore-buffering notes) ONLY while the
# transmitted fraction f < 1; once k saturates (dv_buffered -> 0 identically,
# both on charge and discharge), porosity should stop oscillating within a
# cycle and just track the slow structural-porosity trend. Comparing cycles
# before/at/after the k-saturation throughput (~580 A.h in this run) should
# show that oscillation shrinking to nothing.
# ---------------------------------------------------------------------------
BREATH_TARGET_THR = [50, 550, 850, 1050]
breath_cyc_indices = []
for t in BREATH_TARGET_THR:
    idx = int(np.argmin(np.abs(age_thr - t)))
    breath_cyc_indices.append(int(age_cyc[idx]))

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
    ax_i.set_title(f"~{target_t} A.h (cycle {cyc_idx}, k={k_here:.3f})", fontsize=10)
    ax_i.grid(alpha=0.3)
    ax_i.legend(fontsize=7)

for ax_i in ax_br[-1, :]:
    ax_i.set_xlabel("Throughput within cycle [A.h]")
for ax_i in ax_br[:, 0]:
    ax_i.set_ylabel("Negative electrode\nporosity [-]")

fig_br.suptitle("Pore-buffering \"breathing\" diagnostic: within-cycle porosity oscillation "
                 "vs. structural (irreversible) porosity, at different life stages")
plt.tight_layout()
outpath_br = os.path.join(SCRIPT_DIR, f"pore_buffering_degradation_test_breathing_{PORE_BUFFERING_TRANSITION}.png")
plt.savefig(outpath_br, dpi=150)
print(f"Saved: {outpath_br}")

# ---------------------------------------------------------------------------
# side-by-side comparison: unbuffered (left column) vs. buffered (right
# column), expansion amplitude (top row) and capacity fade (bottom row),
# same throughput x-axis, sharing y-axis per row for a direct visual
# comparison of the buffering effect.
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
fig_cmp.suptitle("Unbuffered vs. buffered: cell-level expansion and capacity fade "
                  f"(si_gr_expansion, transition={PORE_BUFFERING_TRANSITION})")
plt.tight_layout()
outpath_cmp = os.path.join(SCRIPT_DIR, f"pore_buffering_degradation_test_comparison_{PORE_BUFFERING_TRANSITION}.png")
plt.savefig(outpath_cmp, dpi=150)
print(f"Saved: {outpath_cmp}")
