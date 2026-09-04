"""
dryout_molarvol_sweep_test.py -- proper sweep of SEI_MOLAR_VOLUME_DIV, replacing the two
one-off single-point attempts (ec_reaction_limited_dryout_molarvol_test.py at /5.0 then
/1.5) once both turned out to sit on the wrong sides of the trade-off:
  - /1.5: R_dry STILL pinned at exactly 1.0000 (zero dry-out signal) -- yet the tuned
    recipe was already distorted (225 cycles to 50% SoH vs. ~150-165 unmodified baseline,
    LAM_neg 29-34% vs. ~20-23% target).
  - /5.0: R_dry genuinely moved (min 0.9863) -- but the distortion is severe (300 cycles,
    LAM_neg ballooned to ~55%, the sharp tuned knee replaced by a smooth featureless decline).

This sweeps the divisor across both known points plus the gap between them, tracking BOTH
axes of the trade-off at once (min R_dry achieved, AND how much the divisor alone --
independent of the dry-out wrapper, via each point's own "no solvent consumption" baseline
-- distorts cycles-to-50%-SoH / final LLI / final LAM_neg away from the tuned recipe's own
150-165 cycle / ~18% / ~20-23% character), to see whether ANY divisor in this range gives a
genuine R_dry signal without unacceptable distortion, or whether the two are fundamentally
in conflict for this recipe (in which case Phase 1c's dry-out validation may need to lean on
the already-confirmed EC-concentration/self-limiting-SEI channel instead of area/R_dry
dry-out as the representative effect).

Everything else (ec_consumption_scale=0.1, r_eres=0%, soh_target=0.5,
UPDATE_EVERY_N_CYCLES=15, matched solver/var_pts) matches the two prior single-point
scripts, so all three divisor points are directly comparable.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pybamm

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(SCRIPT_DIR, "..", "..", "ec_dryout"))
from ec_dryout_wrapper import (  # noqa: E402
    extract_batch_series,
    run_ec_dryout_degradation,
)

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

UPDATE_EVERY_N_CYCLES = 15
MAX_TOTAL_CYCLES = 400
SOH_TARGET = 0.5
R_ERES = 0.0
EC_CONSUMPTION_SCALE = 0.1
SEI_MOLAR_VOLUME_DEFAULT = 9.585e-05

# The two already-run single-point results, reused as reference (not re-run):
#   /1.5: R_dry pinned at 1.0000, cycles=225, LLI=17.12%, LAM_neg=34.48% (dry-out run)
#   /5.0: R_dry min 0.9863,      cycles=300, LLI=13.03%, LAM_neg=56.34% (dry-out run)
KNOWN_POINTS = {
    1.5: dict(min_R_dry=1.0000, cycles=225, LLI=17.12, LAM_neg=34.48),
    5.0: dict(min_R_dry=0.9863, cycles=300, LLI=13.03, LAM_neg=56.34),
}
# Fill the gap plus one point beyond 5.0, to see the trend clearly.
DIV_SWEEP_NEW = [2.0, 2.5, 3.0, 4.0]

# Tuned-recipe reference (SEI partial molar volume untouched, from
# ec_reaction_limited_dryout_coupling_test.py's robustness-pass run):
UNMODIFIED_REF = dict(cycles=165, LLI=17.98, LAM_neg=22.72)

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


def build_param(molar_vol_div):
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
            "Primary: SEI partial molar volume [m3.mol-1]": SEI_MOLAR_VOLUME_DEFAULT / molar_vol_div,
            "Secondary: SEI partial molar volume [m3.mol-1]": SEI_MOLAR_VOLUME_DEFAULT / molar_vol_div,
        },
        check_already_exists=False,
    )
    return param


def run_case(molar_vol_div):
    param = build_param(molar_vol_div)
    tag = f"[div={molar_vol_div:g}]"
    sol_dryout, history, ledger, traj_dryout = run_ec_dryout_degradation(
        param, formation_exp, ageing_cycle,
        update_every_n_cycles=UPDATE_EVERY_N_CYCLES, max_total_cycles=MAX_TOTAL_CYCLES,
        r_eres=R_ERES, options=MODEL_OPTIONS_BASE, composite=True,
        sei_ec_coupling="ec_reaction", ec_consumption_scale=EC_CONSUMPTION_SCALE,
        soh_target=SOH_TARGET, solver=solver, var_pts=VAR_PTS,
    )
    min_R_dry = min(h["R_dry"] for h in history)
    cycles = UPDATE_EVERY_N_CYCLES * len(history)
    nominal_cap = param["Nominal cell capacity [A.h]"]
    LLI_final = 100 * traj_dryout["LLI"][-1] / nominal_cap if traj_dryout["LLI"].size else float("nan")
    LAM_neg_final = traj_dryout["LAM_neg"][-1] if traj_dryout["LAM_neg"].size else float("nan")
    print(f"{tag} cycles={cycles}  min_R_dry={min_R_dry:.6f}  LLI={LLI_final:.2f}%  "
          f"LAM_neg={LAM_neg_final:.2f}%")
    return dict(
        div=molar_vol_div, traj=traj_dryout, history=history,
        min_R_dry=min_R_dry, cycles=cycles, LLI=LLI_final, LAM_neg=LAM_neg_final,
    )


results = []
for div in DIV_SWEEP_NEW:
    print(f"\n=== SEI_MOLAR_VOLUME_DIV={div:g} ===")
    results.append(run_case(div))

# Merge in the known points (not re-run) for the summary plots
all_points = [dict(div=d, **v) for d, v in KNOWN_POINTS.items()] + [
    dict(div=r["div"], min_R_dry=r["min_R_dry"], cycles=r["cycles"], LLI=r["LLI"], LAM_neg=r["LAM_neg"])
    for r in results
]
all_points.sort(key=lambda d: d["div"])

# ---------------------------------------------------------------------------
# plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(2, 2, figsize=(13, 10))
cmap = plt.cm.viridis(np.linspace(0.1, 0.85, len(results)))

for r, c in zip(results, cmap):
    if r["traj"]["thr"].size:
        cap0 = r["traj"]["cap"][0]
        ax[0, 0].plot(r["traj"]["thr"], 100 * r["traj"]["cap"] / cap0, "o-", ms=3,
                      color=c, label=f"div={r['div']:g}")
ax[0, 0].axhline(50, color="gray", ls="--", lw=1)
ax[0, 0].set_xlabel("Throughput capacity [A.h]")
ax[0, 0].set_ylabel("SoH (discharge capacity) [%]")
ax[0, 0].set_title("SoH vs. throughput (newly-run divisors, dry-out branch)")
ax[0, 0].legend(fontsize=8)
ax[0, 0].grid(alpha=0.3)

divs = [p["div"] for p in all_points]
min_rdry_vals = [p["min_R_dry"] for p in all_points]
ax[0, 1].plot(divs, min_rdry_vals, "o-", color="tab:red")
ax[0, 1].axhline(1.0, color="gray", ls=":", lw=1, label="no dry-out (R_dry=1.0)")
ax[0, 1].set_xlabel("SEI_MOLAR_VOLUME_DIV")
ax[0, 1].set_ylabel("min R_dry achieved")
ax[0, 1].set_title("Dry-out signal strength vs. divisor")
ax[0, 1].legend(fontsize=8)
ax[0, 1].grid(alpha=0.3)

cycles_vals = [p["cycles"] for p in all_points]
ax[1, 0].plot(divs, cycles_vals, "o-", color="tab:purple", label="dry-out run")
ax[1, 0].axhline(UNMODIFIED_REF["cycles"], color="black", ls=":", lw=1.5,
                  label=f"unmodified recipe (~{UNMODIFIED_REF['cycles']} cyc)")
ax[1, 0].set_xlabel("SEI_MOLAR_VOLUME_DIV")
ax[1, 0].set_ylabel("Cycles to 50% SoH")
ax[1, 0].set_title("Recipe distortion: cycle-life shift vs. divisor")
ax[1, 0].legend(fontsize=8)
ax[1, 0].grid(alpha=0.3)

lam_vals = [p["LAM_neg"] for p in all_points]
lli_vals = [p["LLI"] for p in all_points]
ax[1, 1].plot(divs, lli_vals, "o-", color="crimson", label="LLI [%]")
ax[1, 1].plot(divs, lam_vals, "o-", color="darkorange", label="LAM negative [%]")
ax[1, 1].axhline(UNMODIFIED_REF["LLI"], color="crimson", ls=":", lw=1)
ax[1, 1].axhline(UNMODIFIED_REF["LAM_neg"], color="darkorange", ls=":", lw=1)
ax[1, 1].set_xlabel("SEI_MOLAR_VOLUME_DIV")
ax[1, 1].set_ylabel("Final cumulative [%]")
ax[1, 1].set_title("Recipe distortion: LLI/LAM balance vs. divisor (dotted = unmodified)")
ax[1, 1].legend(fontsize=8)
ax[1, 1].grid(alpha=0.3)

fig.suptitle(
    "SEI_MOLAR_VOLUME_DIV sweep -- does ANY divisor give genuine R_dry movement "
    "without unacceptable recipe distortion?"
)
plt.tight_layout()
os.makedirs(os.path.join(SCRIPT_DIR, "pics"), exist_ok=True)
outpath = os.path.join(SCRIPT_DIR, "pics", "dryout_molarvol_sweep_test_result.png")
plt.savefig(outpath, dpi=150)
print(f"\nSaved: {outpath}")

print("\n--- summary (all points, sorted by divisor) ---")
print(f"Unmodified recipe reference: cycles~{UNMODIFIED_REF['cycles']}  "
      f"LLI~{UNMODIFIED_REF['LLI']:.1f}%  LAM_neg~{UNMODIFIED_REF['LAM_neg']:.1f}%")
for p in all_points:
    print(f"div={p['div']:g}: min_R_dry={p['min_R_dry']:.6f}  cycles={p['cycles']}  "
          f"LLI={p['LLI']:.2f}%  LAM_neg={p['LAM_neg']:.2f}%")
