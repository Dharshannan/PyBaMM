"""
ec_reaction_limited_dryout_molarvol_test.py -- run IN PARALLEL with
ec_reaction_limited_dryout_coupling_test.py, testing whether reducing "SEI partial molar
volume" gets R_dry actually moving.

ec_reaction_limited_dryout_coupling_test.py (ec_consumption_scale=0.1, proper 50% SoH
cutoff, matched solver/var_pts between baseline and dry-out) showed R_dry pinned at EXACTLY
1.0000 for the entire run -- no genuine dry-out at all. This matches ec_dryout_wrapper.py's
documented structural finding: si_gr_expansion's default "SEI partial molar volume"
(9.585e-05 m3/mol, both phases) is larger than EC's own molar volume
(M_EC/RHO_EC = 6.67e-05 m3/mol), so simulated pore volume shrinks at least as fast as
electrolyte volume for ANY reaction rate -- `ec_consumption_scale` can't fix this (it only
scales dn_EC, not the ratio between SEI's own solid volume and EC's molar volume).

Mirrors ruihe_dryout_validation.py's own fix for the identical structural issue under
OKane2022: reduce "SEI partial molar volume" so simulated SEI solid is denser/thinner per
mole reacted, letting pore volume shrink SLOWER relative to electrolyte volume.
SEI_MOLAR_VOLUME_DIV=5.0 here matches ruihe_dryout_validation.py's own chosen divisor as a
first attempt (not re-derived from si_gr_expansion's own kinetics, which differ from
OKane2022's) -- everything else (ec_consumption_scale=0.1, r_eres=0%, soh_target=0.5,
UPDATE_EVERY_N_CYCLES=15, matched solver/var_pts) is identical to the sibling script, so
the two can be compared directly to isolate the molar-volume effect alone.
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

# Anchor recipe -- identical to the sibling script.
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

# NEW vs. the sibling script: reduce "SEI partial molar volume" so simulated SEI solid is
# denser/thinner per mole reacted -- pore volume then shrinks SLOWER relative to
# electrolyte volume, the fix ruihe_dryout_validation.py used for the identical structural
# issue under OKane2022. First attempt at /5.0 (that script's own divisor) DID get R_dry
# moving (min 0.9863) but badly distorted the tuned recipe's own knee/LAM balance (cycles to
# 50% SoH doubled 150->300, LAM_neg ballooned 20%->55%, sharp knee replaced by a smooth
# featureless decline) -- si_gr_expansion's SEI/porosity/LAM coupling is evidently far more
# sensitive to this parameter than OKane2022's (LAM/mechanics-off) setup was. /1.5 here is a
# much gentler attempt, trading off a smaller (possibly negligible) R_dry effect for
# hopefully preserving more of the recipe's own tuned character.
SEI_MOLAR_VOLUME_DIV = 1.5
SEI_MOLAR_VOLUME_DEFAULT = 9.585e-05  # si_gr_expansion.py default, both phases

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
            "Primary: SEI partial molar volume [m3.mol-1]": SEI_MOLAR_VOLUME_DEFAULT / SEI_MOLAR_VOLUME_DIV,
            "Secondary: SEI partial molar volume [m3.mol-1]": SEI_MOLAR_VOLUME_DEFAULT / SEI_MOLAR_VOLUME_DIV,
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


print(f"\n=== Baseline: no solvent consumption (SEI molar volume /{SEI_MOLAR_VOLUME_DIV:g}) ===")
base_param = build_param()
_, traj_base = run_baseline_no_dryout(base_param)

print(f"\n=== Dry-out coupled: r_eres={R_ERES:.0%}, ec_consumption_scale={EC_CONSUMPTION_SCALE:g}, "
      f"SEI molar volume /{SEI_MOLAR_VOLUME_DIV:g} ===")
dryout_param = build_param()
sol_dryout, history, ledger, traj_dryout = run_ec_dryout_degradation(
    dryout_param, formation_exp, ageing_cycle,
    update_every_n_cycles=UPDATE_EVERY_N_CYCLES, max_total_cycles=MAX_TOTAL_CYCLES,
    r_eres=R_ERES, options=MODEL_OPTIONS_BASE, composite=True,
    sei_ec_coupling="ec_reaction", ec_consumption_scale=EC_CONSUMPTION_SCALE,
    soh_target=SOH_TARGET, solver=solver, var_pts=VAR_PTS,
)

nominal_cap = base_param["Nominal cell capacity [A.h]"]
cap0 = traj_base["cap"][0]

# ---------------------------------------------------------------------------
# Figure 1: SoH / LLI / LAM(neg) trajectories, baseline vs. dry-out
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(17, 5))

ax[0].plot(traj_base["thr"], 100 * traj_base["cap"] / cap0, "^-", ms=3,
           color="tab:gray", label="no solvent consumption")
ax[0].plot(traj_dryout["thr"], 100 * traj_dryout["cap"] / cap0, "o-", ms=3,
           color="tab:blue", label=f"dry-out, r_eres={R_ERES:.0%}")
ax[0].axhline(50, color="gray", ls="--", lw=1)
ax[0].set_xlabel("Throughput capacity [A.h]")
ax[0].set_ylabel("SoH (discharge capacity) [%]")
ax[0].set_title("(a) State of Health")
ax[0].legend(fontsize=8)
ax[0].grid(alpha=0.3)

ax[1].plot(traj_base["thr"], 100 * traj_base["LLI"] / nominal_cap, "-", lw=2,
           color="tab:gray", label="LLI, no solvent consumption")
ax[1].plot(traj_base["thr"], traj_base["LAM_neg"], "--", lw=2,
           color="tab:gray", label="LAM, no solvent consumption")
ax[1].plot(traj_dryout["thr"], 100 * traj_dryout["LLI"] / nominal_cap, "-", lw=2,
           color="tab:blue", label=f"LLI, dry-out r_eres={R_ERES:.0%}")
ax[1].plot(traj_dryout["thr"], traj_dryout["LAM_neg"], "--", lw=2,
           color="tab:blue", label=f"LAM, dry-out r_eres={R_ERES:.0%}")
ax[1].set_xlabel("Throughput capacity [A.h]")
ax[1].set_ylabel("Final cumulative [%]  (solid=LLI, dashed=LAM)")
ax[1].set_title("(b) LLI + LAM (negative electrode)")
ax[1].legend(fontsize=7)
ax[1].grid(alpha=0.3)

ax[2].plot(traj_base["thr"], traj_base["cell_amplitude"] * 1e6, "^-", ms=3,
           color="tab:gray", label="no solvent consumption")
ax[2].plot(traj_dryout["thr"], traj_dryout["cell_amplitude"] * 1e6, "o-", ms=3,
           color="tab:blue", label=f"dry-out, r_eres={R_ERES:.0%}")
ax[2].set_xlabel("Throughput capacity [A.h]")
ax[2].set_ylabel("Cell-level within-cycle amplitude [um]")
ax[2].set_title("(c) Reversible expansion amplitude")
ax[2].legend(fontsize=8)
ax[2].grid(alpha=0.3)

fig.suptitle(
    f"'ec reaction limited' + dry-out -- SEI partial molar volume /{SEI_MOLAR_VOLUME_DIV:g}, "
    f"ec_consumption_scale={EC_CONSUMPTION_SCALE:g}, r_eres={R_ERES:.0%}, cut off at {SOH_TARGET:.0%} SoH"
)
plt.tight_layout()
os.makedirs(os.path.join(SCRIPT_DIR, "pics"), exist_ok=True)
out1 = os.path.join(SCRIPT_DIR, "pics", "ec_reaction_limited_dryout_molarvol_test_trajectories.png")
plt.savefig(out1, dpi=150)
print(f"\nSaved: {out1}")

# ---------------------------------------------------------------------------
# Figure 2: every dry-out ledger diagnostic vs. cycle number
# ---------------------------------------------------------------------------
cyc_x = UPDATE_EVERY_N_CYCLES * np.arange(1, len(history) + 1)
fig2, ax2 = plt.subplots(3, 2, figsize=(14, 12), sharex=True)

ax2[0, 0].plot(cyc_x, [h["R_dry"] for h in history], "o-", ms=4, color="tab:blue")
ax2[0, 0].set_ylabel("R_dry\n(V_eJR / V_pore)")
ax2[0, 0].set_title("Electrolyte dry-out ratio")
ax2[0, 0].grid(alpha=0.3)

ax2[0, 1].plot(cyc_x, [h["R_Li"] for h in history], "o-", ms=4, color="tab:blue")
ax2[0, 1].set_ylabel("R_Li\n(Li+ mixing ratio)")
ax2[0, 1].set_title("Lithium-ion concentration change ratio")
ax2[0, 1].grid(alpha=0.3)

ax2[1, 0].plot(cyc_x, [h["c_EC_new"] for h in history], "o-", ms=4, color="tab:blue")
ax2[1, 0].set_ylabel("Bulk EC concentration\n[mol/m3]")
ax2[1, 0].set_title("EC concentration")
ax2[1, 0].grid(alpha=0.3)

ax2[1, 1].plot(cyc_x, [h["V_eres"] for h in history], "o-", ms=4, color="tab:blue")
ax2[1, 1].set_ylabel("Reservoir electrolyte\nvolume V_eres [m3]")
ax2[1, 1].set_title(f"Reservoir depletion (r_eres={R_ERES:.0%} -> always 0 here)")
ax2[1, 1].grid(alpha=0.3)

ax2[2, 0].plot(cyc_x, [h["dn_EC"] for h in history], "o-", ms=4, color="tab:blue")
ax2[2, 0].set_ylabel("EC consumed this batch\ndn_EC [mol]")
ax2[2, 0].set_xlabel("Ageing cycle number")
ax2[2, 0].set_title("Per-batch EC consumption")
ax2[2, 0].grid(alpha=0.3)

ax2[2, 1].plot(cyc_x, [h["n_parallel_new"] for h in history], "o-", ms=4, color="tab:blue")
ax2[2, 1].set_ylabel("Effective parallel electrodes\n(cumulative area scale)")
ax2[2, 1].set_xlabel("Ageing cycle number")
ax2[2, 1].set_title("Cumulative electrode-area shrink")
ax2[2, 1].grid(alpha=0.3)

fig2.suptitle(
    f"Dry-out ledger diagnostics -- SEI partial molar volume /{SEI_MOLAR_VOLUME_DIV:g}, "
    f"ec_consumption_scale={EC_CONSUMPTION_SCALE:g}, r_eres={R_ERES:.0%}"
)
plt.tight_layout()
out2 = os.path.join(SCRIPT_DIR, "pics", "ec_reaction_limited_dryout_molarvol_test_diagnostics.png")
plt.savefig(out2, dpi=150)
print(f"Saved: {out2}")

print("\n--- summary ---")
print(f"Baseline (no solvent consumption): final capacity retention = "
      f"{100 * traj_base['cap'][-1] / cap0:.2f}%  final LLI={100*traj_base['LLI'][-1]/nominal_cap:.2f}%  "
      f"final LAM_neg={traj_base['LAM_neg'][-1]:.2f}%")
print(f"Dry-out (r_eres={R_ERES:.0%}): final capacity retention = "
      f"{100 * traj_dryout['cap'][-1] / cap0:.2f}%  final LLI={100*traj_dryout['LLI'][-1]/nominal_cap:.2f}%  "
      f"final LAM_neg={traj_dryout['LAM_neg'][-1]:.2f}%")
total_dnEC = sum(h["dn_EC"] for h in history)
print(f"Total EC consumed: {total_dnEC:.3e} mol  "
      f"final R_dry={history[-1]['R_dry']:.6f}  min R_dry={min(h['R_dry'] for h in history):.6f}  "
      f"final c_EC={history[-1]['c_EC_new']:.1f} mol/m3 (BOL={ledger.c_EC0:.1f} mol/m3)  "
      f"final effective n_parallel={history[-1]['n_parallel_new']:.6f} (BOL={ledger.n_parallel_0:.6f})")
if min(h["R_dry"] for h in history) > 0.9999:
    print(
        f"\nNOTE: R_dry still pinned at ~1.0000 even after dividing SEI partial molar "
        f"volume by {SEI_MOLAR_VOLUME_DIV:g} -- the ratio still favours pore-volume shrinkage "
        "over electrolyte-volume shrinkage at this divisor; a larger divisor (or boosting "
        "D_ec/c_ec_0 alongside it, mirroring ruihe_dryout_validation.py's paired diffusivity "
        "boost) would be the next thing to try."
    )
else:
    print(
        f"\nGenuine dry-out engaged: min R_dry={min(h['R_dry'] for h in history):.6f} < 1.0 "
        f"with SEI partial molar volume /{SEI_MOLAR_VOLUME_DIV:g}."
    )
