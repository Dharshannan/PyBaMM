"""
dryout_reres_sweep_test.py -- reproduces Li et al. 2022's own headline validation
structure (Fig. 3: baseline "no solvent consumption" vs. r_eres in {0%, 6%, 9%}) against
the "ec reaction limited" recipe, per the paper passage: "cells without solvent
consumption will overestimate the LLI due to SEI but underestimate the LAM... the cell
with 9% initial extra electrolyte can retain 7% more capacity... compared with the cell
with 0% initial extra electrolyte."

UPDATE: switched from the "SEI partial molar volume" divisor to the new
ec_molar_volume_scale lever (ec_dryout_wrapper.py), per explicit follow-up instruction --
the recipe should keep its own tuned KNEE shape (like the /1.5 divisor's relatively
preserved shape), not the smooth featureless decline /5.0 produced. dryout_molarvol_sweep_
test.py's own data explains why "SEI partial molar volume" can't give both: the knee was
already eroding by divisor~2 (cycles 165->285), well before R_dry moved meaningfully
(still pinned at EXACTLY 1.0 through divisor 2.0, only 0.9996 at 2.5) -- and worse, at
divisor<=2.0 (which includes 1.5), R_dry stays pinned at EXACTLY 1.0 because
ECDryoutLedger.update's reservoir-refill step (`dV_add = min(shortfall, V_eres)`) only ever
draws on V_eres when `shortfall = max(V_pore_new - V_eJR_predryout, 0) > 0` -- with
shortfall=0 throughout, r_eres=0%/6%/9% would produce IDENTICAL trajectories regardless of
how much reservoir electrolyte is specified, making a divisor=1.5 r_eres sweep completely
uninformative, not just "gentle."

ec_molar_volume_scale (new, wrapper-only) fixes this cleanly: it scales dV_EC (electrolyte
volume consumed per mole EC reacted) directly in the ledger's OWN bookkeeping, leaving
"SEI partial molar volume" -- and therefore PyBaMM's own simulated porosity-decline/LAM/
knee physics -- completely untouched. SEI_MOLAR_VOLUME_DIV is back to 1.0 (i.e. not
applied at all) here; EC_MOLAR_VOLUME_SCALE=2.0 is a first candidate, chosen to clear the
~1.44x crossover (SEI partial molar volume's default 9.585e-05 m3/mol / EC's own molar
volume 6.667e-05 m3/mol) with some margin, without yet knowing how much R_dry movement (or
recipe distortion, if any) it produces -- that's what this run itself will show.

Produces two figures, both extending ruihe_dryout_validation.py's own structure with a LAM
panel (that script's OKane2022 config had LAM disabled; this recipe's is very much on):
  - dryout_reres_sweep_test_trajectories.png: SoH / LLI / LAM(neg) vs. throughput, baseline
    + all three r_eres conditions overlaid -- cf. paper Fig. 3(a-c).
  - dryout_reres_sweep_test_diagnostics.png: R_dry / R_Li / c_EC / V_eres / dn_EC /
    n_parallel vs. cycle, one colour per r_eres condition -- cf. paper Fig. A-3(d-f)/A-6.
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
# "SEI partial molar volume" is left at si_gr_expansion.py's own default (NOT divided) --
# preserves the tuned recipe's exact knee shape. R_dry is instead driven by the new
# ec_molar_volume_scale wrapper lever below (see module docstring).

UPDATE_EVERY_N_CYCLES = 1  # was 15 -- per-cycle granularity, avoids the batch-boundary
                            # SoH-overshoot artifact that confounded the r_eres comparison
                            # at N=15 (steep late-stage collapse landing well past the 50%
                            # check before it could be caught)
MAX_TOTAL_CYCLES = 400
SOH_TARGET = 0.5
EC_CONSUMPTION_SCALE = 0.1
EC_MOLAR_VOLUME_SCALE = 5.0  # scale=2.0: no effect (R_dry=1.0). scale=10.0: total collapse
                              # (R_dry crashed to 0.0, c_EC spiked to 15773 then crashed --
                              # a genuine runaway feedback). scale=5.0: min R_dry=0.7854, a
                              # real, substantial, but gradual/well-behaved dip -- the
                              # bisection's landing point.

# Paper's own three r_eres cases.
R_ERES_SWEEP = [0.0, 0.06, 0.09]

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
            # "SEI partial molar volume" deliberately NOT overridden -- stays at
            # si_gr_expansion.py's own default, preserving the tuned recipe's knee shape.
        },
        check_already_exists=False,
    )
    return param


def run_baseline_no_dryout(base_param):
    param = base_param.copy()
    sim0 = pybamm.Simulation(
        pybamm.lithium_ion.DFN(MODEL_OPTIONS_BASE), parameter_values=param,
        experiment=formation_exp, solver=solver, var_pts=VAR_PTS,
    )
    sol = sim0.solve(initial_soc=1.0)
    print("[baseline] formation cycle solved.")

    traj_parts = dict(thr=[], cap=[], LLI=[], LAM_neg=[], LAM_pos=[], cell_amplitude=[])
    total_cycles = 0
    while total_cycles < MAX_TOTAL_CYCLES:
        batch_exp = pybamm.Experiment([ageing_cycle for _ in range(UPDATE_EVERY_N_CYCLES)])
        model = pybamm.lithium_ion.DFN(MODEL_OPTIONS_BASE)
        sim = pybamm.Simulation(
            model, parameter_values=param, experiment=batch_exp,
            solver=solver, var_pts=VAR_PTS,
        )
        sol = sim.solve(starting_solution=sol.last_state)
        total_cycles += UPDATE_EVERY_N_CYCLES

        batch_series = extract_batch_series(sol, composite=True)
        for key in traj_parts:
            traj_parts[key].append(batch_series[key])

        soh_now = (
            batch_series["cap"][-1] / traj_parts["cap"][0][0]
            if batch_series["cap"].size else float("nan")
        )
        print(f"[baseline] {total_cycles} cycles  SoH={100*soh_now:.2f}%")
        if not np.isnan(soh_now) and soh_now <= SOH_TARGET:
            print(f"[baseline] reached SOH_TARGET ({SOH_TARGET:.0%}) at {total_cycles} cycles -- stopping.")
            break

    trajectory = {
        key: (np.concatenate(parts) if parts else np.array([]))
        for key, parts in traj_parts.items()
    }
    return sol, trajectory


print(f"\n=== Baseline: no solvent consumption (SEI partial molar volume UNTOUCHED) ===")
base_param = build_param()
_, traj_base = run_baseline_no_dryout(base_param)
nominal_cap = base_param["Nominal cell capacity [A.h]"]
cap0 = traj_base["cap"][0]

conditions = []
for r_eres in R_ERES_SWEEP:
    print(f"\n=== Dry-out coupled: r_eres={r_eres:.0%} ===")
    dryout_param = build_param()
    _, history, ledger, traj = run_ec_dryout_degradation(
        dryout_param, formation_exp, ageing_cycle,
        update_every_n_cycles=UPDATE_EVERY_N_CYCLES, max_total_cycles=MAX_TOTAL_CYCLES,
        r_eres=r_eres, options=MODEL_OPTIONS_BASE, composite=True,
        sei_ec_coupling="ec_reaction", ec_consumption_scale=EC_CONSUMPTION_SCALE,
        ec_molar_volume_scale=EC_MOLAR_VOLUME_SCALE,
        soh_target=SOH_TARGET, solver=solver, var_pts=VAR_PTS,
    )
    conditions.append((r_eres, traj, history))

colors = {0.0: "tab:red", 0.06: "tab:orange", 0.09: "tab:green"}

# ---------------------------------------------------------------------------
# Figure 1: SoH / (LLI+LAM combined) / reversible expansion -- cf. paper Fig. 3(a-c),
# extended with the expansion panel this project tracks throughout.
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(17, 5))

ax[0].plot(traj_base["thr"], 100 * traj_base["cap"] / cap0, "^-", ms=3,
           color="tab:gray", label="no solvent consumption")
for r_eres, traj, _ in conditions:
    ax[0].plot(traj["thr"], 100 * traj["cap"] / cap0, "o-", ms=3,
               color=colors[r_eres], label=f"{r_eres:.0%} extra electrolyte")
ax[0].axhline(50, color="gray", ls="--", lw=1)
ax[0].set_xlabel("Throughput capacity [A.h]")
ax[0].set_ylabel("SoH (discharge capacity) [%]")
ax[0].set_title("(a) Capacity retention")
ax[0].legend(fontsize=8)
ax[0].grid(alpha=0.3)

ax[1].plot(traj_base["thr"], 100 * traj_base["LLI"] / nominal_cap, "-", lw=2,
           color="tab:gray", label="LLI, no solvent consumption")
ax[1].plot(traj_base["thr"], traj_base["LAM_neg"], "--", lw=2,
           color="tab:gray", label="LAM, no solvent consumption")
for r_eres, traj, _ in conditions:
    ax[1].plot(traj["thr"], 100 * traj["LLI"] / nominal_cap, "-", lw=2,
               color=colors[r_eres], label=f"LLI, {r_eres:.0%}")
    ax[1].plot(traj["thr"], traj["LAM_neg"], "--", lw=2,
               color=colors[r_eres], label=f"LAM, {r_eres:.0%}")
ax[1].set_xlabel("Throughput capacity [A.h]")
ax[1].set_ylabel("Final cumulative [%]  (solid=LLI, dashed=LAM)")
ax[1].set_title("(b) LLI + LAM (negative electrode)")
ax[1].legend(fontsize=6, ncol=2)
ax[1].grid(alpha=0.3)

ax[2].plot(traj_base["thr"], traj_base["cell_amplitude"] * 1e6, "^-", ms=3,
           color="tab:gray", label="no solvent consumption")
for r_eres, traj, _ in conditions:
    ax[2].plot(traj["thr"], traj["cell_amplitude"] * 1e6, "o-", ms=3,
               color=colors[r_eres], label=f"{r_eres:.0%}")
ax[2].set_xlabel("Throughput capacity [A.h]")
ax[2].set_ylabel("Cell-level within-cycle amplitude [um]")
ax[2].set_title("(c) Reversible expansion amplitude")
ax[2].legend(fontsize=8)
ax[2].grid(alpha=0.3)

fig.suptitle(
    f"'ec reaction limited' dry-out r_eres sweep -- ec_molar_volume_scale={EC_MOLAR_VOLUME_SCALE:g} "
    "(SEI partial molar volume untouched), "
    f"ec_consumption_scale={EC_CONSUMPTION_SCALE:g}, cut off at {SOH_TARGET:.0%} SoH -- cf. paper Fig. 3(a-c)"
)
plt.tight_layout()
os.makedirs(os.path.join(SCRIPT_DIR, "pics"), exist_ok=True)
out1 = os.path.join(SCRIPT_DIR, "pics", "dryout_reres_sweep_test_trajectories.png")
plt.savefig(out1, dpi=150)
print(f"\nSaved: {out1}")

# ---------------------------------------------------------------------------
# Figure 2: every ledger diagnostic vs. cycle -- cf. paper Fig. A-3(d-f)/A-6
# ---------------------------------------------------------------------------
fig2, ax2 = plt.subplots(3, 2, figsize=(14, 12), sharex=True)
for r_eres, traj, history in conditions:
    cyc_x = UPDATE_EVERY_N_CYCLES * np.arange(1, len(history) + 1)
    ax2[0, 0].plot(cyc_x, [h["R_dry"] for h in history], "o-", ms=4,
                   color=colors[r_eres], label=f"{r_eres:.0%}")
    ax2[0, 1].plot(cyc_x, [h["R_Li"] for h in history], "o-", ms=4, color=colors[r_eres])
    ax2[1, 0].plot(cyc_x, [h["c_EC_new"] for h in history], "o-", ms=4, color=colors[r_eres])
    ax2[1, 1].plot(cyc_x, [h["V_eres"] for h in history], "o-", ms=4, color=colors[r_eres])
    ax2[2, 0].plot(cyc_x, [h["dn_EC"] for h in history], "o-", ms=4, color=colors[r_eres])
    ax2[2, 1].plot(cyc_x, [h["n_parallel_new"] for h in history], "o-", ms=4, color=colors[r_eres])

ax2[0, 0].set_ylabel("R_dry\n(V_eJR / V_pore)")
ax2[0, 0].set_title("Electrolyte dry-out ratio")
ax2[0, 0].legend(fontsize=8)
ax2[0, 0].grid(alpha=0.3)
ax2[0, 1].set_ylabel("R_Li\n(Li+ mixing ratio)")
ax2[0, 1].set_title("Lithium-ion concentration change ratio")
ax2[0, 1].grid(alpha=0.3)
ax2[1, 0].set_ylabel("Bulk EC concentration\n[mol/m3]")
ax2[1, 0].set_title("EC concentration")
ax2[1, 0].grid(alpha=0.3)
ax2[1, 1].set_ylabel("Reservoir electrolyte\nvolume V_eres [m3]")
ax2[1, 1].set_title("Reservoir depletion")
ax2[1, 1].grid(alpha=0.3)
ax2[2, 0].set_ylabel("EC consumed this batch\ndn_EC [mol]")
ax2[2, 0].set_xlabel("Ageing cycle number")
ax2[2, 0].set_title("Per-batch EC consumption")
ax2[2, 0].grid(alpha=0.3)
ax2[2, 1].set_ylabel("Effective parallel electrodes\n(cumulative area scale)")
ax2[2, 1].set_xlabel("Ageing cycle number")
ax2[2, 1].set_title("Cumulative electrode-area shrink")
ax2[2, 1].grid(alpha=0.3)

fig2.suptitle(
    f"Dry-out ledger diagnostics -- 'ec reaction limited', ec_molar_volume_scale={EC_MOLAR_VOLUME_SCALE:g}, "
    f"ec_consumption_scale={EC_CONSUMPTION_SCALE:g}"
)
plt.tight_layout()
out2 = os.path.join(SCRIPT_DIR, "pics", "dryout_reres_sweep_test_diagnostics.png")
plt.savefig(out2, dpi=150)
print(f"Saved: {out2}")

print("\n--- summary ---")
print(f"Baseline (no solvent consumption): final capacity retention = "
      f"{100 * traj_base['cap'][-1] / cap0:.2f}%  final LLI={100*traj_base['LLI'][-1]/nominal_cap:.2f}%  "
      f"final LAM_neg={traj_base['LAM_neg'][-1]:.2f}%")
for r_eres, traj, history in conditions:
    total_dnEC = sum(h["dn_EC"] for h in history)
    print(f"r_eres={r_eres:.0%}: final capacity retention = "
          f"{100 * traj['cap'][-1] / cap0:.2f}%  final LLI={100*traj['LLI'][-1]/nominal_cap:.2f}%  "
          f"final LAM_neg={traj['LAM_neg'][-1]:.2f}%  cycles={UPDATE_EVERY_N_CYCLES*len(history)}  "
          f"min R_dry={min(h['R_dry'] for h in history):.6f}  total EC consumed={total_dnEC:.3e} mol  "
          f"final c_EC={history[-1]['c_EC_new']:.1f} mol/m3")
