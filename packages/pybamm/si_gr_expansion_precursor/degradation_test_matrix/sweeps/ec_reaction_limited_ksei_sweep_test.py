"""
ec_reaction_limited_ksei_sweep_test.py -- single consolidated k_sei sweep, replacing the
five earlier one-off sweep scripts (ksei_retune, dec_retune, ksei_fine_sweep, width_f0,
ksei_refine_validated_width -- all deleted). Reproduces ec_reaction_limited_ksei_retune_
test.py's original 4-panel plot layout (SoH-vs-throughput, expansion amplitude, knee-vs-
k_sei, LLI/LAM balance), but over a single sweep spanning the full range explored across
those earlier scripts: K_SEI_MULT from 0.02 down to 1e-05, plus 0.0017 (the single-point
match found by the now-deleted ksei_fine_sweep_test.py) inserted explicitly since it does
not sit on the coarse log grid.

D_ec/c_ec_0 are left at si_gr_expansion.py's defaults throughout (k_sei is the only lever
swept here) -- see degradation_test_matrix_plan.md for the fuller multi-lever history
(D_ec sweep, width/f0 sweep) that is no longer reproduced by a standalone script but
remains documented there.

Per instruction, the "reaction limited" full-lifetime throughput reference line/legend
entry (~1,100-1,116 A.h, shown in the original ksei_retune_test_result.png) has been
removed from both the SoH and knee-timing panels -- the LLI/LAM% target lines in the
balance panel are retained.
"""
import os
import matplotlib.pyplot as plt
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

# Base-case recipe -- identical to ec_reaction_limited_swap_test.py / the validated
# "reaction limited" baseline. Held fixed throughout this sweep -- only K_SEI_MULT varies.
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
K_SEI_DEFAULT = 1e-12  # si_gr_expansion.py's untouched default, both phases

BATCH_SIZE = 20
MAX_TOTAL_CYCLES = 400
SOH_TARGET = 0.5

# Reference LLI/LAM balance from ec_reaction_limited_swap_test.py's "reaction limited" run
# (throughput target intentionally not included -- removed from the plot per instruction).
REACTION_LIMITED_REF = dict(cycles=160, LLI=33.0, LAM_neg=23.4)

# Full sweep range: 0.02 down to 1e-05, plus 0.0017 (ksei_fine_sweep_test.py's best single-
# point match) inserted in sorted position.
K_SEI_MULT_SWEEP = [0.02, 0.01, 0.003, 0.0017, 0.001, 0.0003, 0.0001, 0.00001]

solver = pybamm.IDAKLUSolver(root_tol=1e-06, atol=1e-06, rtol=1e-06)

formation_exp = pybamm.Experiment(
    [
        "Discharge at 0.1C until 2.5 V",
        "Charge at C/3 until 4.2 V",
        "Hold at 4.2 V until C/100",
    ]
)
ageing_cycle = (
    "Discharge at C/3 until 2.5 V",
    "Charge at C/3 until 4.2 V",
    "Hold at 4.2 V until C/100",
)
ageing_batch_exp = pybamm.Experiment([ageing_cycle for _ in range(BATCH_SIZE)])


def build_param_updates(k_sei_mult):
    return {
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
        "Primary: SEI kinetic rate constant [m.s-1]": K_SEI_DEFAULT * k_sei_mult,
        "Secondary: SEI kinetic rate constant [m.s-1]": K_SEI_DEFAULT * k_sei_mult,
        # D_ec/c_ec_0 left at si_gr_expansion.py's defaults -- k_sei is the only lever swept.
    }


def cycle_ageing_leg(cyc):
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
                cap = c
                discharge_end_thr = float(thr_step[-1])
        elif mean_I < -1e-3:
            charge_end_thr = float(thr_step[-1])
    return cap, discharge_end_thr, charge_end_thr


def current_soh(sol):
    caps = []
    for cyc in sol.cycles:
        cap, discharge_end_thr, charge_end_thr = cycle_ageing_leg(cyc)
        if cap > 0 and discharge_end_thr is not None and charge_end_thr is not None:
            caps.append(cap)
    if len(caps) < 2:
        return None
    return caps[-1] / caps[0]


def run_case(label, k_sei_mult):
    param = pybamm.ParameterValues("si_gr_expansion")
    param.update(build_param_updates(k_sei_mult), check_already_exists=False)
    model = pybamm.lithium_ion.DFN(MODEL_OPTIONS_BASE)

    def make_sim(experiment):
        return pybamm.Simulation(
            model, parameter_values=param, experiment=experiment,
            solver=solver, var_pts=VAR_PTS,
        )

    tag = f"[{label}]"
    sim = make_sim(formation_exp)
    sol = sim.solve(initial_soc=1.0)
    print(f"{tag} formation solved.")

    total_cycles = 0
    stop_reason = "reached MAX_TOTAL_CYCLES safety cap"
    while total_cycles < MAX_TOTAL_CYCLES:
        sim = make_sim(ageing_batch_exp)
        try:
            sol = sim.solve(starting_solution=sol)
        except Exception as exc:  # noqa: BLE001
            print(f"{tag} solver failure after {total_cycles} cycles: "
                  f"{type(exc).__name__}: {str(exc)[:300]}")
            stop_reason = "solver_failure"
            break
        total_cycles += BATCH_SIZE
        soh_now = current_soh(sol)
        print(f"{tag} {total_cycles} cycles solved. "
              f"SoH={'n/a' if soh_now is None else f'{soh_now:.3f}'}")
        if soh_now is not None and soh_now <= SOH_TARGET:
            stop_reason = f"reached SOH_TARGET ({SOH_TARGET:.0%})"
            break

    print(f"{tag} stop reason: {stop_reason}")

    Qt_full = sol["Throughput capacity [A.h]"].entries
    tc_cell_full = sol["Cell thickness change [m]"].entries
    k_full = sol["Negative electrode transfer ratio k"].entries

    thr_list, cap_list, discharge_end_list, charge_end_list = [], [], [], []
    for cyc in sol.cycles:
        cap, discharge_end_thr, charge_end_thr = cycle_ageing_leg(cyc)
        if cap > 0 and discharge_end_thr is not None and charge_end_thr is not None:
            thr_list.append(float(cyc["Throughput capacity [A.h]"].entries[-1]))
            cap_list.append(cap)
            discharge_end_list.append(discharge_end_thr)
            charge_end_list.append(charge_end_thr)

    thr_arr = np.array(thr_list)
    cap_arr = np.array(cap_list)
    soh_arr = cap_arr / cap_arr[0] if len(cap_arr) else cap_arr
    discharge_end_arr = np.array(discharge_end_list)
    charge_end_arr = np.array(charge_end_list)

    if thr_arr.size:
        tc_charge_end = np.interp(charge_end_arr, Qt_full, tc_cell_full)
        tc_discharge_end = np.interp(discharge_end_arr, Qt_full, tc_cell_full)
        cell_amplitude = tc_charge_end - tc_discharge_end
        k_at_charge_end = np.interp(charge_end_arr, Qt_full, k_full)
    else:
        cell_amplitude = np.array([])
        k_at_charge_end = np.array([])

    LLI_final = float(sol["Loss of lithium inventory [%]"].entries[-1])
    LAM_neg_final = float(sol["Loss of active material in negative electrode [%]"].entries[-1])

    below_90 = np.where(soh_arr < 0.9)[0]
    knee_thr = float(thr_arr[below_90[0]]) if below_90.size else float("nan")
    full_thr = float(thr_arr[-1]) if thr_arr.size else float("nan")

    print(f"{tag} final: LLI={LLI_final:.2f}%  LAM_neg={LAM_neg_final:.2f}%  "
          f"final SoH={soh_arr[-1] if soh_arr.size else float('nan'):.3f}  "
          f"cycles={total_cycles}  knee(SoH<90%)~{knee_thr:.1f} A.h  "
          f"full_throughput~{full_thr:.1f} A.h  stop_reason={stop_reason}")

    return dict(
        label=label, k_sei_mult=k_sei_mult,
        thr=thr_arr, cap=cap_arr, soh=soh_arr, cell_amplitude=cell_amplitude,
        k=k_at_charge_end, LLI_final=LLI_final, LAM_neg_final=LAM_neg_final,
        knee_thr=knee_thr, full_thr=full_thr, total_cycles=total_cycles,
        stop_reason=stop_reason,
    )


results = []
for k_sei_mult in K_SEI_MULT_SWEEP:
    print(f"\n=== k_sei_mult={k_sei_mult:g} (k_sei={K_SEI_DEFAULT*k_sei_mult:.2e} m/s) ===")
    results.append(run_case(f"x{k_sei_mult:g}", k_sei_mult))

# ---------------------------------------------------------------------------
# plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(2, 2, figsize=(13, 10))
cmap = plt.cm.plasma(np.linspace(0.05, 0.85, len(results)))

for r, c in zip(results, cmap):
    if r["thr"].size:
        ax[0, 0].plot(r["thr"], r["soh"] * 100, "o-", ms=3, color=c,
                      label=f"{r['label']} (knee~{r['knee_thr']:.0f} A.h)")
ax[0, 0].axhline(50, color="gray", ls="--", lw=1)
ax[0, 0].set_xlabel("Throughput capacity [A.h]")
ax[0, 0].set_ylabel("SoH (discharge capacity) [%]")
ax[0, 0].set_title("State of Health vs. throughput")
ax[0, 0].legend(fontsize=7)
ax[0, 0].grid(alpha=0.3)

for r, c in zip(results, cmap):
    if r["thr"].size:
        ax[0, 1].plot(r["thr"], r["cell_amplitude"] * 1e6, "o-", ms=3, color=c, label=r["label"])
ax[0, 1].set_xlabel("Throughput capacity [A.h]")
ax[0, 1].set_ylabel("Cell-level within-cycle amplitude [um]")
ax[0, 1].set_title("Reversible expansion amplitude")
ax[0, 1].legend(fontsize=7)
ax[0, 1].grid(alpha=0.3)

k_sei_vals = [r["k_sei_mult"] for r in results]
knee_vals = [r["knee_thr"] for r in results]
ax[1, 0].plot(k_sei_vals, knee_vals, "o-", color="tab:purple")
ax[1, 0].set_xscale("log")
ax[1, 0].set_xlabel("k_sei multiplier (log scale)")
ax[1, 0].set_ylabel("Knee throughput (SoH<90%) [A.h]")
ax[1, 0].set_title("Knee timing vs. k_sei")
ax[1, 0].grid(alpha=0.3)

x = np.arange(len(results))
width = 0.35
ax[1, 1].bar(x - width / 2, [r["LLI_final"] for r in results], width, label="LLI [%]", color="crimson")
ax[1, 1].bar(x + width / 2, [r["LAM_neg_final"] for r in results], width, label="LAM negative [%]", color="darkorange")
ax[1, 1].axhline(REACTION_LIMITED_REF["LLI"], color="crimson", ls=":", lw=1)
ax[1, 1].axhline(REACTION_LIMITED_REF["LAM_neg"], color="darkorange", ls=":", lw=1)
ax[1, 1].set_xticks(x)
ax[1, 1].set_xticklabels([r["label"] for r in results])
ax[1, 1].set_ylabel("Final cumulative [%]")
ax[1, 1].set_title("LLI vs. LAM balance (dotted = reaction limited target)")
ax[1, 1].legend(fontsize=8)
ax[1, 1].grid(alpha=0.3)

fig.suptitle(
    "'ec reaction limited' k_sei sweep (0.02 -> 1e-05, plus 0.0017 fine-sweep point)"
)
plt.tight_layout()
os.makedirs(os.path.join(SCRIPT_DIR, "..", "pics"), exist_ok=True)
outpath = os.path.join(SCRIPT_DIR, "..", "pics", "ec_reaction_limited_ksei_retune_test_result.png")
plt.savefig(outpath, dpi=150)
print(f"\nSaved: {outpath}")

print("\n--- summary ---")
print(f"Reference LLI/LAM balance (reaction limited): cycles={REACTION_LIMITED_REF['cycles']}  "
      f"LLI={REACTION_LIMITED_REF['LLI']:.1f}%  LAM_neg={REACTION_LIMITED_REF['LAM_neg']:.1f}%")
for r in results:
    print(f"{r['label']}: cycles={r['total_cycles']}  "
          f"final SoH={r['soh'][-1] if r['soh'].size else float('nan'):.3f}  "
          f"knee(SoH<90%)~{r['knee_thr']:.1f} A.h  full_throughput~{r['full_thr']:.1f} A.h  "
          f"LLI={r['LLI_final']:.2f}%  LAM_neg={r['LAM_neg_final']:.2f}%  "
          f"stop_reason={r['stop_reason']}")
