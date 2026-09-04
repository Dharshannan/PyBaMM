"""
ec_reaction_limited_1000cyc_stretch_test.py -- root-causes why EVERY stretch attempt so
far (v1 uniform, v2 LAM-boosted, v3 intermediate) produced a hump noticeably sharper than
the unstretched anchor's, even v1, which preserved the LAM:k_sei:crack-rate ratio exactly.
That ruled out "the LAM/SEI ratio was broken" as the (sole) explanation -- something about
the stretch is wrong even when every OUTER rate constant is stretched uniformly.

Root cause, found by checking whether TIMESCALE_STRETCH is actually an exact time-
reparametrization of "ec reaction limited"'s own kinetics, not just of the LAM/crack
submodels: j_sei = -F*c0*k_exp/(1 + (L_sei/D_ec)*k_exp), k_exp = k_sei*exp(-a*F/RT*eta).
For a stretched system (k_sei' = k_sei/S) to reproduce the SAME j_sei(L_sei) curve as the
anchor, just S times slower (i.e. dL_sei/d(tau) = (1/S) dL_sei/dt at every L_sei -- a clean
slow-motion replay), solving j_sei(L; k_sei/S, D_ec') = j_sei(L; k_sei, D_ec)/S for all L
requires BOTH k_sei' = k_sei/S AND D_ec' = D_ec/S -- the diffusion-length parameter has to
shrink by the SAME factor, not stay at its default. v1/v2/v3 all left D_ec untouched at
si_gr_expansion.py's default (2e-18) while dividing only k_sei by S -- this makes
(L_sei/D_ec)*k_exp come out S times SMALLER than the self-similar target at any given
L_sei, pushing the stretched run further into the kinetics-dominated regime than the
anchor was at the same accumulated SEI thickness. Since L_sei grows largest late in life
(exactly where the buffering transition/hump sits), this distorts precisely the region
that looked "sharp" -- consistent with why even v1 (correct LAM:crack:k_sei ratio) still
came out sharper than the anchor: the SEI-growth ODE itself was never made self-similar,
only the LAM/crack submodels were.

Fix: divide D_ec by the SAME TIMESCALE_STRETCH as k_sei (in addition to LAM proportional
terms and crack rate, all uniformly, exactly as v1 did for everything else). If this
restores genuine self-similarity, cycle count should land close to the clean 7x
extrapolation (~1120, not v1's 840-cycle shortfall), which would ALSO resolve v1's
LAM_neg undershoot as a side effect -- without needing v2/v3's ad-hoc LAM_STRETCH
decoupling at all. f0/width stay at the ORIGINAL (0.7, 0.01) -- not the A.0l/A.0m
validated point, per instruction.
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

K_SEI_ANCHOR_MULT = 0.0017
K_SEI_DEFAULT = 1e-12
D_EC_DEFAULT = 2e-18  # si_gr_expansion.py default -- NEVER divided by stretch until now

# Original (never-changed) pore-buffering parameters.
SI_CRIT_STRESS = 2.2e8
SI_LAM_EXP = 2.5
GR_CRIT_STRESS = 3.0e7
NEG_POROSITY_FLOOR = 0.01
EXPONENT_MAX_SEI = 10.0
F0_BASELINE = 0.7
WIDTH_BASELINE = 0.01

SI_LAM_PROP_BASELINE = 7.5e-7
GR_LAM_PROP_BASELINE = 4.5e-6
SI_CRACK_RATE_MULT_BASELINE = 0.1
GR_CRACK_RATE_MULT_BASELINE = 0.1
BASE_CRACK_RATE = 3.9e-20

# Single, fully-uniform stretch factor -- now ALSO applied to D_ec, restoring the SEI
# growth ODE's own self-similarity (see module docstring). No separate LAM_STRETCH.
TIMESCALE_STRETCH = 7.0

BATCH_SIZE = 20
MAX_TOTAL_CYCLES = 1500
SOH_TARGET = 0.5

RPT_INTERVAL = BATCH_SIZE
RPT_RATE = "C/10"
assert BATCH_SIZE % RPT_INTERVAL == 0, "BATCH_SIZE must be a multiple of RPT_INTERVAL"

REACTION_LIMITED_REF = dict(knee=457.1, full=1115.6, cycles=160, LLI=33.0, LAM_neg=23.4)
ANCHOR_REF = dict(knee=408.3, full=1061.6, cycles=160, LLI=33.22, LAM_neg=21.68)
V1_UNIFORM_STRETCH_REF = dict(knee=2477.4, full=5990.6, cycles=840, LLI=32.64, LAM_neg=16.10)
V2_LAM47_REF = dict(knee=2202.5, full=5683.9, cycles=800, LLI=32.36, LAM_neg=23.10)

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


def build_param_updates(stretch):
    return {
        "Primary: Negative electrode active material volume fraction": GR_EPS,
        "Secondary: Negative electrode active material volume fraction": SI_EPS,
        "Secondary: Maximum concentration in negative electrode [mol.m-3]": SI_MAX_CONC_NEEDED,
        "Secondary: Initial concentration in negative electrode [mol.m-3]": scaled_si_initial_conc(SI_MAX_CONC_NEEDED),
        "Nominal cell capacity [A.h]": NOMINAL_CAP_AH,
        "Secondary: Negative electrode Young's modulus [Pa]": 1.0e10,
        "Secondary: Negative electrode LAM constant proportional term [s-1]": SI_LAM_PROP_BASELINE / stretch,
        "Secondary: Negative electrode LAM constant exponential term": SI_LAM_EXP,
        "Secondary: Negative electrode critical stress [Pa]": SI_CRIT_STRESS,
        "Secondary: Negative electrode cracking rate": BASE_CRACK_RATE * SI_CRACK_RATE_MULT_BASELINE / stretch,
        "Primary: Negative electrode LAM constant exponential term": 2.0,
        "Primary: Negative electrode critical stress [Pa]": GR_CRIT_STRESS,
        "Primary: Negative electrode LAM constant proportional term [s-1]": GR_LAM_PROP_BASELINE / stretch,
        "Primary: Negative electrode cracking rate": BASE_CRACK_RATE * GR_CRACK_RATE_MULT_BASELINE / stretch,
        "Negative electrode porosity floor": NEG_POROSITY_FLOOR,
        "Primary: SEI reaction exponent cap": EXPONENT_MAX_SEI,
        "Secondary: SEI reaction exponent cap": EXPONENT_MAX_SEI,
        "Negative electrode pore buffering transition width": WIDTH_BASELINE,
        "Negative electrode transmitted fraction plateau": F0_BASELINE,
        "Primary: SEI kinetic rate constant [m.s-1]": K_SEI_DEFAULT * K_SEI_ANCHOR_MULT / stretch,
        "Secondary: SEI kinetic rate constant [m.s-1]": K_SEI_DEFAULT * K_SEI_ANCHOR_MULT / stretch,
        # NEW vs. v1/v2/v3: D_ec now divided by the SAME stretch factor as k_sei, to keep
        # the "ec reaction limited" j_sei formula's own kinetics-vs-diffusion balance
        # (L_sei/D_ec * k_exp) self-similar across the stretch, not just the outer
        # LAM/crack rates.
        "Primary: EC diffusivity [m2.s-1]": D_EC_DEFAULT / stretch,
        "Secondary: EC diffusivity [m2.s-1]": D_EC_DEFAULT / stretch,
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


def rpt_leg(cyc):
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
        if 0.1 < mean_I <= 1.0 and q.size >= 2:
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


def run_case(label, stretch):
    param = pybamm.ParameterValues("si_gr_expansion")
    param.update(build_param_updates(stretch), check_already_exists=False)
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
    rpt_thr_list, rpt_cap_list, rpt_discharge_end_list, rpt_charge_end_list = [], [], [], []
    for cyc in sol.cycles:
        cap, discharge_end_thr, charge_end_thr = cycle_ageing_leg(cyc)
        if cap > 0 and discharge_end_thr is not None and charge_end_thr is not None:
            thr_list.append(float(cyc["Throughput capacity [A.h]"].entries[-1]))
            cap_list.append(cap)
            discharge_end_list.append(discharge_end_thr)
            charge_end_list.append(charge_end_thr)
            continue
        rpt_cap, rpt_discharge_end, rpt_charge_end = rpt_leg(cyc)
        if rpt_cap > 0 and rpt_discharge_end is not None and rpt_charge_end is not None:
            rpt_thr_list.append(float(cyc["Throughput capacity [A.h]"].entries[-1]))
            rpt_cap_list.append(rpt_cap)
            rpt_discharge_end_list.append(rpt_discharge_end)
            rpt_charge_end_list.append(rpt_charge_end)

    thr_arr = np.array(thr_list)
    cap_arr = np.array(cap_list)
    soh_arr = cap_arr / cap_arr[0] if len(cap_arr) else cap_arr
    discharge_end_arr = np.array(discharge_end_list)
    charge_end_arr = np.array(charge_end_list)

    rpt_thr_arr = np.array(rpt_thr_list)
    rpt_cap_arr = np.array(rpt_cap_list)
    rpt_discharge_end_arr = np.array(rpt_discharge_end_list)
    rpt_charge_end_arr = np.array(rpt_charge_end_list)
    rpt_soh_arr = rpt_cap_arr / rpt_cap_arr[0] if rpt_cap_arr.size else rpt_cap_arr

    if thr_arr.size:
        tc_charge_end = np.interp(charge_end_arr, Qt_full, tc_cell_full)
        tc_discharge_end = np.interp(discharge_end_arr, Qt_full, tc_cell_full)
        cell_amplitude = tc_charge_end - tc_discharge_end
        k_at_charge_end = np.interp(charge_end_arr, Qt_full, k_full)
    else:
        cell_amplitude = np.array([])
        k_at_charge_end = np.array([])

    rpt_cell_amplitude = np.array([])
    if rpt_thr_arr.size:
        rpt_tc_charge_end = np.interp(rpt_charge_end_arr, Qt_full, tc_cell_full)
        rpt_tc_discharge_end = np.interp(rpt_discharge_end_arr, Qt_full, tc_cell_full)
        rpt_cell_amplitude = rpt_tc_charge_end - rpt_tc_discharge_end

    LLI_final = float(sol["Loss of lithium inventory [%]"].entries[-1])
    LAM_neg_final = float(sol["Loss of active material in negative electrode [%]"].entries[-1])

    below_90 = np.where(soh_arr < 0.9)[0]
    knee_thr = float(thr_arr[below_90[0]]) if below_90.size else float("nan")
    full_thr = float(thr_arr[-1]) if thr_arr.size else float("nan")
    peak_idx = int(np.argmax(cell_amplitude)) if cell_amplitude.size else None
    peak_thr = float(thr_arr[peak_idx]) if peak_idx is not None else float("nan")
    post_peak_frac = (100.0 * (full_thr - peak_thr) / full_thr) if full_thr else float("nan")

    print(f"{tag} final: LLI={LLI_final:.2f}%  LAM_neg={LAM_neg_final:.2f}%  "
          f"final SoH={soh_arr[-1] if soh_arr.size else float('nan'):.3f}  "
          f"cycles={total_cycles}  knee(SoH<90%)~{knee_thr:.1f} A.h  "
          f"full_throughput~{full_thr:.1f} A.h  peak_thr~{peak_thr:.1f} A.h  "
          f"post_peak_frac~{post_peak_frac:.1f}%  n_rpt={rpt_thr_arr.size}  "
          f"stop_reason={stop_reason}")

    return dict(
        label=label, stretch=stretch,
        thr=thr_arr, cap=cap_arr, soh=soh_arr, cell_amplitude=cell_amplitude,
        k=k_at_charge_end, LLI_final=LLI_final, LAM_neg_final=LAM_neg_final,
        knee_thr=knee_thr, full_thr=full_thr, total_cycles=total_cycles,
        peak_thr=peak_thr, post_peak_frac=post_peak_frac,
        stop_reason=stop_reason,
        rpt_thr=rpt_thr_arr, rpt_soh=rpt_soh_arr, rpt_cell_amplitude=rpt_cell_amplitude,
    )


print(f"\n=== TIMESCALE_STRETCH={TIMESCALE_STRETCH:g} applied UNIFORMLY to k_sei, D_ec, "
      f"LAM proportional terms, AND crack rate (v4: D_ec now included) ===")
result = run_case(f"stretch={TIMESCALE_STRETCH:g}_with_Dec", TIMESCALE_STRETCH)

# ---------------------------------------------------------------------------
# plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(2, 2, figsize=(13, 10))

ax[0, 0].plot(result["thr"], result["soh"] * 100, "o-", ms=3, color="tab:blue",
              label=f"C/3 ageing (full~{result['full_thr']:.0f} A.h, {result['total_cycles']} cyc)")
if result["rpt_thr"].size:
    ax[0, 0].scatter(result["rpt_thr"], result["rpt_soh"] * 100, marker="^", color="tab:orange",
                      s=30, label=f"RPT ({RPT_RATE})", zorder=5)
ax[0, 0].axhline(50, color="gray", ls="--", lw=1)
ax[0, 0].set_xlabel("Throughput capacity [A.h]")
ax[0, 0].set_ylabel("SoH (discharge capacity) [%]")
ax[0, 0].set_title("State of Health vs. throughput (v4: D_ec also stretched)")
ax[0, 0].legend(fontsize=8)
ax[0, 0].grid(alpha=0.3)

ax[0, 1].plot(result["thr"], result["cell_amplitude"] * 1e6, "o-", ms=3, color="tab:blue",
              label="C/3 ageing")
if result["rpt_thr"].size:
    ax[0, 1].scatter(result["rpt_thr"], result["rpt_cell_amplitude"] * 1e6, marker="^",
                      color="tab:orange", s=30, label=f"RPT ({RPT_RATE})", zorder=5)
ax[0, 1].axvline(result["peak_thr"], color="gray", ls=":", lw=1,
                  label=f"peak (post-peak span {result['post_peak_frac']:.0f}% of life)")
ax[0, 1].set_xlabel("Throughput capacity [A.h]")
ax[0, 1].set_ylabel("Cell-level within-cycle amplitude [um]")
ax[0, 1].set_title("Reversible expansion amplitude (v4: D_ec also stretched)")
ax[0, 1].legend(fontsize=8)
ax[0, 1].grid(alpha=0.3)

ax[1, 0].plot(result["thr"], result["k"], "o-", ms=3, color="tab:blue")
ax[1, 0].axhline(0.7, color="gray", ls=":", lw=1, label="BoL f0=0.7")
ax[1, 0].set_xlabel("Throughput capacity [A.h]")
ax[1, 0].set_ylabel("k [-]")
ax[1, 0].set_title("Internal transfer ratio k")
ax[1, 0].legend(fontsize=8)
ax[1, 0].grid(alpha=0.3)

labels = ["reaction limited\n(target)", "anchor\n(unstretched, 160cyc)",
          "v1: uniform x7,\nD_ec NOT stretched\n(840cyc)", "v2: LAM/4.7,\nD_ec NOT stretched\n(800cyc)",
          f"v4: uniform x7,\nD_ec stretched too\n({result['total_cycles']}cyc)"]
lli_vals = [REACTION_LIMITED_REF["LLI"], ANCHOR_REF["LLI"], V1_UNIFORM_STRETCH_REF["LLI"],
            V2_LAM47_REF["LLI"], result["LLI_final"]]
lam_vals = [REACTION_LIMITED_REF["LAM_neg"], ANCHOR_REF["LAM_neg"], V1_UNIFORM_STRETCH_REF["LAM_neg"],
            V2_LAM47_REF["LAM_neg"], result["LAM_neg_final"]]
xpos = np.arange(len(labels))
width_plot = 0.35
ax2 = ax[1, 1]
ax2.bar(xpos - width_plot / 2, lli_vals, width_plot, label="LLI [%]", color="crimson")
ax2.bar(xpos + width_plot / 2, lam_vals, width_plot, label="LAM negative [%]", color="darkorange")
ax2.axhline(REACTION_LIMITED_REF["LLI"], color="crimson", ls=":", lw=1)
ax2.axhline(REACTION_LIMITED_REF["LAM_neg"], color="darkorange", ls=":", lw=1)
ax2.set_xticks(xpos)
ax2.set_xticklabels(labels, fontsize=7)
ax2.set_ylabel("Final cumulative [%]")
ax2.set_title("LLI vs. LAM balance across all variants (dotted = target)")
ax2.legend(fontsize=8)
ax2.grid(alpha=0.3)

fig.suptitle(
    f"v4: f0={F0_BASELINE:g}/width={WIDTH_BASELINE:g} (ORIGINAL), TIMESCALE_STRETCH={TIMESCALE_STRETCH:g} "
    "applied to k_sei, D_ec, LAM, AND crack rate -- testing whether D_ec self-similarity fixes the sharpness"
)
plt.tight_layout()
os.makedirs(os.path.join(SCRIPT_DIR, "pics"), exist_ok=True)
outpath = os.path.join(SCRIPT_DIR, "pics", "ec_reaction_limited_1000cyc_stretch_test_result.png")
plt.savefig(outpath, dpi=150)
print(f"\nSaved: {outpath}")

print("\n--- summary ---")
print(f"Reaction limited target: knee={REACTION_LIMITED_REF['knee']:.1f} A.h  "
      f"full={REACTION_LIMITED_REF['full']:.1f} A.h  cycles={REACTION_LIMITED_REF['cycles']}  "
      f"LLI={REACTION_LIMITED_REF['LLI']:.1f}%  LAM_neg={REACTION_LIMITED_REF['LAM_neg']:.1f}%")
print(f"Anchor (unstretched): knee={ANCHOR_REF['knee']:.1f} A.h  full={ANCHOR_REF['full']:.1f} A.h  "
      f"cycles={ANCHOR_REF['cycles']}  LLI={ANCHOR_REF['LLI']:.1f}%  LAM_neg={ANCHOR_REF['LAM_neg']:.1f}%")
print(f"v1 (uniform x7, D_ec NOT stretched): full~{V1_UNIFORM_STRETCH_REF['full']:.1f} A.h  "
      f"cycles={V1_UNIFORM_STRETCH_REF['cycles']}  LLI={V1_UNIFORM_STRETCH_REF['LLI']:.2f}%  "
      f"LAM_neg={V1_UNIFORM_STRETCH_REF['LAM_neg']:.2f}%")
print(f"v4 (uniform x7, D_ec stretched too): full~{result['full_thr']:.1f} A.h  "
      f"cycles={result['total_cycles']}  LLI={result['LLI_final']:.2f}%  "
      f"LAM_neg={result['LAM_neg_final']:.2f}%  peak_thr~{result['peak_thr']:.1f} A.h  "
      f"post_peak_frac~{result['post_peak_frac']:.1f}% (anchor's own post_peak_frac ~27-29%)  "
      f"stop_reason={result['stop_reason']}")
