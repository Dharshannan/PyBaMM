"""
ec_reaction_limited_dryout_coupling_test.py -- Phase 1a/1b's first actual coupling of
ec_dryout_wrapper.py's electrolyte dry-out/solvent-consumption loop to the "ec reaction
limited" short (~160-cycle) recipe validated in ../ec_reaction_limited_swap_test.py and
../sweeps/ec_reaction_limited_ksei_sweep_test.py (k_sei x0.0017, f0=0.7, width=0.01,
composite Si/Gr, pore buffering "physical").

Per instruction, this is a FIRST estimate, deliberately unrecalibrated: run the wrapper at
si_gr_expansion.py's own BASE dry-out-relevant parameters ("SEI partial molar volume",
"EC diffusivity", etc. all untouched) and r_eres=0% (no extra/reservoir electrolyte -- the
paper's most severe, buffer-free case), coupled via sei_ec_coupling="ec_reaction" (the mode
implemented in ec_dryout_wrapper.py's implementation_plan.md section 7, matching this
recipe's "ec reaction limited" SEI option). This is NOT expected to be the final calibrated
answer -- ec_dryout_wrapper.py's own module history documents that si_gr_expansion's
default "SEI partial molar volume" (9.585e-05 m3/mol) is structurally larger than EC's own
molar volume (M_EC/RHO_EC = 6.67e-05 m3/mol), which pinned R_dry at exactly 1.0000
(no genuine dry-out) the one other time this composite recipe was tried (under
"solvent-diffusion limited" coupling) -- the same structural ratio applies here regardless
of SEI submodel, since the ledger's mass-balance mechanics (ECDryoutLedger.update) are
shared across both coupling modes. This script's job is to confirm/deny that empirically
for THIS recipe before any recalibration is attempted (Phase 1b), not to assume it.

Bug found and fixed while building this script: ec_dryout_wrapper.py's
sei_capacity_loss_Ah() only summed bulk "...SEI [A.h]" capacity loss, not the separate
"...SEI on cracks [A.h]" pathway -- since every recipe in this project runs "SEI on
cracks": "true", this under-counted EC consumption (and therefore dry-out severity)
whenever cracking is active. Fixed directly in ec_dryout_wrapper.py (wrapper-only, additive,
degrades to +0 wherever the variable isn't built -- see that function's docstring),
confirmed not to affect ruihe_dryout_validation.py's already-validated OKane2022 path
(no particle mechanics there, so the cracks term is always 0).

UPDATE (recalibration pass, per follow-up instruction): the first pass (ec_consumption_
scale=1.0, i.e. the paper's exact 1:1:1 stoichiometry, run to a fixed MAX_TOTAL_CYCLES=240
rather than a proper SoH cutoff) showed R_dry genuinely move (0.9535 min, NOT pinned at
1.0000 like the earlier "solvent_diffusion"-coupling attempt) but consumed essentially the
ENTIRE BOL EC pool (~0.0244 mol consumed against a ~0.0244 mol pool) by cycle 240 -- since
"ec reaction limited"'s j_sei is directly proportional to the fed-back c_EC, this
self-limited the SEI/LLI channel almost to a halt (final LLI 4.46% vs. the no-dry-out
baseline's 23.33% AT THE SAME CYCLE COUNT), which is why "final capacity retention" looked
HIGHER under dry-out (57.28%) than the baseline (34.00%) -- the opposite of the usual
"dry-out worsens degradation" narrative, and not a representative first estimate. Two fixes,
both added to ec_dryout_wrapper.py (wrapper-only, backward-compatible new keyword args,
default values preserve every previously-validated behaviour):
  1. `ec_consumption_scale` (new ECDryoutLedger/run_ec_dryout_degradation kwarg, default
     1.0): a multiplier on dn_EC, stretching how many cycles the same BOL EC pool lasts
     before full depletion without touching anything else. Set to 0.1 here as a first
     recalibration attempt (order-of-magnitude reduction, matching this project's usual
     first-pass convention).
  2. `soh_target` (new run_ec_dryout_degradation kwarg, default None): stops the batch loop
     as soon as running SoH drops to/below this value, exactly like this project's other
     degradation_test_matrix scripts' current_soh()-based termination -- fixes the bug this
     script originally had (and the baseline loop below also had, independently) of running
     to a fixed cycle count regardless of SoH, which is why the first pass's "final capacity
     retention" numbers (57%/34%) were never actually AT 50% SoH.

UPDATE (robustness pass, per follow-up instruction): with the two fixes above, the dry-out
run's SoH curve showed one spurious single-cycle dropout to ~0% right near the 50% SoH
cutoff, absent from the baseline. Root cause found while investigating: `run_ec_dryout_
degradation` never accepted a `solver`/`var_pts` override at all -- every one of its
internal `pybamm.Simulation(...)` calls (formation AND every ageing batch) ran on PyBaMM's
own DEFAULTS, while this script's baseline loop explicitly used the tightened
`IDAKLUSolver` + custom `VAR_PTS` below -- i.e. the two runs were never using matched
numerics, an apples-to-oranges comparison, not just a "tolerance too loose" issue. Fixed by
adding `solver=None, var_pts=None` passthrough parameters to `run_ec_dryout_degradation`
(wrapper-only, default `None` preserves every previously-validated call site's behaviour,
including `ruihe_dryout_validation.py`'s), and this script now passes its own `solver`/
`VAR_PTS` into BOTH the baseline and the dry-out call -- this (matched numerics), not
tolerance level, is the fix; tolerance was tried at `1e-08` then reverted back to the
original `1e-06` per follow-up instruction (unnecessary once the two runs share the same
solver/var_pts). `UPDATE_EVERY_N_CYCLES` reduced (`20` -> `15`, shorter batches near
end-of-life), per instruction.

Produces two figures:
  - ec_reaction_limited_dryout_coupling_test_trajectories.png: SoH / LLI / LAM(neg) vs.
    throughput, baseline (no solvent consumption) vs. dry-out-coupled (r_eres=0%).
  - ec_reaction_limited_dryout_coupling_test_diagnostics.png: every ledger diagnostic
    (R_dry, R_Li, bulk c_EC, reservoir volume V_eres, EC consumed per batch dn_EC,
    cumulative parallel-electrode/area scale) vs. ageing cycle number.
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

# Anchor recipe -- IDENTICAL to ec_reaction_limited_swap_test.py / ksei_sweep's x0.0017
# best-match point. Dry-out-specific parameters (SEI partial molar volume, EC diffusivity,
# EC initial concentration) are left at si_gr_expansion.py's own defaults -- this is the
# "start with the base, no extra electrolyte" first estimate, not a recalibration.
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

UPDATE_EVERY_N_CYCLES = 15  # was 20 -- shorter batches near end-of-life, see module docstring
MAX_TOTAL_CYCLES = 400  # safety cap; SOH_TARGET below is the real stopping criterion
SOH_TARGET = 0.5
R_ERES = 0.0  # "no extra electrolyte" -- the paper's most severe, buffer-free case

# Recalibration: the first pass (ec_consumption_scale=1.0, the paper's exact 1:1:1
# stoichiometry) consumed ~0.0244 mol EC against a ~0.0244 mol BOL pool by cycle 240 --
# i.e. essentially the ENTIRE EC pool was exhausted, which (since "ec reaction limited"'s
# j_sei is directly proportional to the fed-back c_EC) self-limited the SEI/LLI channel
# almost completely rather than showing a representative dry-out trajectory. Reducing this
# scale stretches how many cycles the same BOL pool lasts before full depletion, without
# changing anything else (BOL geometry, R_dry/R_Li mechanics, or the recipe itself).
EC_CONSUMPTION_SCALE = 0.1

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
            # Dry-out-relevant parameters (SEI partial molar volume, EC diffusivity,
            # EC initial concentration) deliberately NOT touched -- base estimate.
        },
        check_already_exists=False,
    )
    return param


def run_baseline_no_dryout(base_param):
    """Same recipe, but c_EC and electrode area are never updated -- the
    "cell without solvent consumption" comparison, mirroring
    ruihe_dryout_validation.py's run_baseline_no_dryout."""
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


print(f"\n=== Baseline: no solvent consumption ===")
base_param = build_param()
_, traj_base = run_baseline_no_dryout(base_param)

print(f"\n=== Dry-out coupled: r_eres={R_ERES:.0%}, ec_consumption_scale={EC_CONSUMPTION_SCALE:g} ===")
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
# Figure 1: SoH / (LLI+LAM combined) / reversible expansion, baseline vs. dry-out
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(17, 5))

ax[0].plot(traj_base["thr"], 100 * traj_base["cap"] / cap0, "^-", ms=3,
           color="tab:gray", label="no solvent consumption")
ax[0].plot(traj_dryout["thr"], 100 * traj_dryout["cap"] / cap0, "o-", ms=3,
           color="tab:red", label=f"dry-out, r_eres={R_ERES:.0%}")
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
           color="tab:red", label=f"LLI, dry-out r_eres={R_ERES:.0%}")
ax[1].plot(traj_dryout["thr"], traj_dryout["LAM_neg"], "--", lw=2,
           color="tab:red", label=f"LAM, dry-out r_eres={R_ERES:.0%}")
ax[1].set_xlabel("Throughput capacity [A.h]")
ax[1].set_ylabel("Final cumulative [%]  (solid=LLI, dashed=LAM)")
ax[1].set_title("(b) LLI + LAM (negative electrode)")
ax[1].legend(fontsize=7)
ax[1].grid(alpha=0.3)

ax[2].plot(traj_base["thr"], traj_base["cell_amplitude"] * 1e6, "^-", ms=3,
           color="tab:gray", label="no solvent consumption")
ax[2].plot(traj_dryout["thr"], traj_dryout["cell_amplitude"] * 1e6, "o-", ms=3,
           color="tab:red", label=f"dry-out, r_eres={R_ERES:.0%}")
ax[2].set_xlabel("Throughput capacity [A.h]")
ax[2].set_ylabel("Cell-level within-cycle amplitude [um]")
ax[2].set_title("(c) Reversible expansion amplitude")
ax[2].legend(fontsize=8)
ax[2].grid(alpha=0.3)

fig.suptitle(
    "'ec reaction limited' + dry-out wrapper coupling -- "
    f"ec_consumption_scale={EC_CONSUMPTION_SCALE:g}, r_eres={R_ERES:.0%} (no extra electrolyte), "
    f"cut off at {SOH_TARGET:.0%} SoH"
)
plt.tight_layout()
os.makedirs(os.path.join(SCRIPT_DIR, "pics"), exist_ok=True)
out1 = os.path.join(SCRIPT_DIR, "pics", "ec_reaction_limited_dryout_coupling_test_trajectories.png")
plt.savefig(out1, dpi=150)
print(f"\nSaved: {out1}")

# ---------------------------------------------------------------------------
# Figure 2: every dry-out ledger diagnostic vs. cycle number
# ---------------------------------------------------------------------------
cyc_x = UPDATE_EVERY_N_CYCLES * np.arange(1, len(history) + 1)
fig2, ax2 = plt.subplots(3, 2, figsize=(14, 12), sharex=True)

ax2[0, 0].plot(cyc_x, [h["R_dry"] for h in history], "o-", ms=4, color="tab:red")
ax2[0, 0].set_ylabel("R_dry\n(V_eJR / V_pore)")
ax2[0, 0].set_title("Electrolyte dry-out ratio")
ax2[0, 0].grid(alpha=0.3)

ax2[0, 1].plot(cyc_x, [h["R_Li"] for h in history], "o-", ms=4, color="tab:red")
ax2[0, 1].set_ylabel("R_Li\n(Li+ mixing ratio)")
ax2[0, 1].set_title("Lithium-ion concentration change ratio")
ax2[0, 1].grid(alpha=0.3)

ax2[1, 0].plot(cyc_x, [h["c_EC_new"] for h in history], "o-", ms=4, color="tab:red")
ax2[1, 0].set_ylabel("Bulk EC concentration\n[mol/m3]")
ax2[1, 0].set_title("EC concentration")
ax2[1, 0].grid(alpha=0.3)

ax2[1, 1].plot(cyc_x, [h["V_eres"] for h in history], "o-", ms=4, color="tab:red")
ax2[1, 1].set_ylabel("Reservoir electrolyte\nvolume V_eres [m3]")
ax2[1, 1].set_title(f"Reservoir depletion (r_eres={R_ERES:.0%} -> always 0 here)")
ax2[1, 1].grid(alpha=0.3)

ax2[2, 0].plot(cyc_x, [h["dn_EC"] for h in history], "o-", ms=4, color="tab:red")
ax2[2, 0].set_ylabel("EC consumed this batch\ndn_EC [mol]")
ax2[2, 0].set_xlabel("Ageing cycle number")
ax2[2, 0].set_title("Per-batch EC consumption")
ax2[2, 0].grid(alpha=0.3)

ax2[2, 1].plot(cyc_x, [h["n_parallel_new"] for h in history], "o-", ms=4, color="tab:red")
ax2[2, 1].set_ylabel("Effective parallel electrodes\n(cumulative area scale)")
ax2[2, 1].set_xlabel("Ageing cycle number")
ax2[2, 1].set_title("Cumulative electrode-area shrink")
ax2[2, 1].grid(alpha=0.3)

fig2.suptitle(
    "Dry-out ledger diagnostics -- 'ec reaction limited', "
    f"ec_consumption_scale={EC_CONSUMPTION_SCALE:g}, r_eres={R_ERES:.0%}"
)
plt.tight_layout()
out2 = os.path.join(SCRIPT_DIR, "pics", "ec_reaction_limited_dryout_coupling_test_diagnostics.png")
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
        "\nNOTE: R_dry stayed pinned at ~1.0000 throughout -- no genuine dry-out engaged. "
        "This matches ec_dryout_wrapper.py's documented structural finding: si_gr_expansion's "
        "default 'SEI partial molar volume' (9.585e-05 m3/mol) is larger than EC's own molar "
        "volume (M_EC/RHO_EC ~= 6.67e-05 m3/mol), so simulated pore volume shrinks at least as "
        "fast as electrolyte volume regardless of reaction rate -- recalibration (e.g. reducing "
        "SEI partial molar volume, mirroring ruihe_dryout_validation.py's OKane2022 calibration) "
        "is the expected next step (Phase 1b), not a wrapper bug."
    )
