"""
ruihe_dryout_validation.py -- electrolyte dry-out demonstration for the
ec_dryout_wrapper.py implementation of Li et al. 2022 (J. Electrochem.
Soc. 169 060516), using the OKane2022 parameter set (a single,
non-composite negative electrode) with everything except SEI growth +
SEI-driven porosity change turned off (no LAM, no mechanics/cracking, no
pore buffering). Compares a baseline (no solvent-consumption wrapper) run
against the dry-out wrapper at three initial extra-electrolyte fractions
(0%, 6%, 9%), matching the paper's own r_eres sweep.

An earlier version of this script used si_gr_expansion (composite Si/Gr
negative electrode, "pore buffering": "true") instead -- that run already
validated the wrapper works correctly WITH pore buffering enabled (see
implementation_plan.md section 4.4), but never showed genuine dry-out:
R_dry stayed pinned at exactly 1.0000 throughout, because si_gr_expansion's
SEI/porosity parameterisation shrinks pore volume faster than electrolyte
volume for any reaction rate (a structural, rate-independent ratio, not a
wrapper bug). This version switches to OKane2022 with "pore buffering":
"false" and calibrated SEI kinetics specifically to demonstrate the
dry-out mechanism ITSELF (R_dry actually dropping below 1).

Produces:
  - ruihe_dryout_validation_fig3.png: capacity retention / LLI (SEI) vs
    throughput, for all four conditions -- cf. paper Fig. 3(a-b). (LAM is
    omitted -- disabled in MODEL_OPTIONS, always zero here.)
  - ruihe_dryout_validation_ratios.png: R_dry / R_Li / bulk EC
    concentration vs cycle number, for the three dry-out conditions -- cf.
    paper Fig. A-3(d-f) / A-6.

--- Why the SEI kinetics needed calibrating (not just "increase the rate") ---

The default OKane2022 "SEI solvent diffusivity" (2.5e-22 m2/s, same order
as si_gr_expansion's) makes SEI growth (and therefore EC consumption) so
slow that no measurable dry-out occurs within any practical cycle count.
The natural fix -- just boosting the diffusivity -- was tried first and
found NOT to work on its own: even at a 3000x boost, over a full 200
cycles, R_dry stayed pinned at exactly 1.0000 throughout (verified). The
reason is structural, not a rate problem: the SEI current density j_SEI is
proportional to c_EC (Eq. 4), so BOTH the EC-consumption rate (which
shrinks V_eJR, per the wrapper's external M_EC/RHO_EC constants) and
PyBaMM's own porosity-decrease rate (which shrinks V_pore, via "SEI
partial molar volume") are driven by the exact SAME j_SEI and hence decay
together, self-limitingly, as c_EC drops -- boosting the OVERALL rate
just makes both decay faster in lockstep, without changing their RATIO.
That ratio is set structurally by OKane2022's "SEI partial molar volume"
(9.585e-05 m3/mol) being large relative to the wrapper's EC molar volume
(M_EC/RHO_EC = 6.67e-05 m3/mol) -- i.e. this parameterisation implies MORE
solid volume (hence pore-volume shrinkage) forms per mole of EC reacted
than the actual EC molecule's own volume, so V_pore always shrinks at
least as fast as V_eJR, regardless of overall rate.

Fix: reduce "SEI partial molar volume [m3.mol-1]" (making PyBaMM's own
simulated SEI solid denser/thinner per mole reacted) alongside boosting
the diffusivity (so the effect is visible within a practical cycle
count). Calibrated by direct testing (see implementation_plan.md):
  - Diffusivity boost alone (3000x): R_dry stayed at 1.0000 for 200 cycles.
  - Diffusivity x1000, molar volume /10: R_dry hit 0.9075 within just 10
    cycles and 20% capacity fade within 40 -- too aggressive for a 200-
    cycle demo.
  - Diffusivity x100, molar volume /5 (CHOSEN): R_dry settles into a
    sustained ~0.97-0.996 range across the full 200 cycles, 86.3%
    capacity retention at cycle 200 (r_eres=0% case) -- a gradual, clearly
    visible but not catastrophic dry-out trajectory.
This is a transparent, documented calibration choice for DEMONSTRATING
the mechanism, not a claim about OKane2022's real physical SEI density.
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import pybamm

from ec_dryout_wrapper import extract_batch_series, run_ec_dryout_degradation

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

pybamm.set_logging_level("NOTICE")

MODEL_OPTIONS = {
    "SEI": "solvent-diffusion limited",
    "SEI porosity change": "true",
    "pore buffering": "false",
}

# Calibrated multipliers -- see module docstring for the calibration history.
SEI_DIFFUSIVITY_MULT = 100.0
SEI_MOLAR_VOLUME_DIV = 5.0

UPDATE_EVERY_N_CYCLES = 20
MAX_TOTAL_CYCLES = 200

FORMATION_EXP = pybamm.Experiment(
    [
        "Discharge at 0.1C until 2.5 V",
        "Charge at C/3 until 4.2 V",
        "Hold at 4.2 V until C/100",
    ]
)
AGEING_CYCLE = (
    "Discharge at 1C until 2.5 V",
    "Charge at 0.3C until 4.2 V",
    "Hold at 4.2 V until C/100",
)


def make_param():
    param = pybamm.ParameterValues("OKane2022")
    param.update(
        {
            "SEI solvent diffusivity [m2.s-1]": (
                param["SEI solvent diffusivity [m2.s-1]"] * SEI_DIFFUSIVITY_MULT
            ),
            "SEI partial molar volume [m3.mol-1]": (
                param["SEI partial molar volume [m3.mol-1]"] / SEI_MOLAR_VOLUME_DIV
            ),
        },
        check_already_exists=False,
    )
    return param


def run_baseline_no_dryout(base_param):
    """Same calibrated SEI kinetics, but c_EC and electrode area are NEVER
    updated (paper's "cell without solvent consumption" comparison)."""
    param = base_param.copy()
    sim0 = pybamm.Simulation(
        pybamm.lithium_ion.DFN(MODEL_OPTIONS), parameter_values=param,
        experiment=FORMATION_EXP,
    )
    sol = sim0.solve(initial_soc=1.0)
    print("[baseline] formation cycle solved.")

    traj_parts = dict(thr=[], cap=[], LLI=[], LAM_neg=[], LAM_pos=[])
    total_cycles = 0
    while total_cycles < MAX_TOTAL_CYCLES:
        batch_exp = pybamm.Experiment(
            [AGEING_CYCLE for _ in range(UPDATE_EVERY_N_CYCLES)]
        )
        model = pybamm.lithium_ion.DFN(MODEL_OPTIONS)
        sim = pybamm.Simulation(model, parameter_values=param, experiment=batch_exp)
        sol = sim.solve(starting_solution=sol.last_state)
        total_cycles += UPDATE_EVERY_N_CYCLES

        batch_series = extract_batch_series(sol, composite=False)
        for key in traj_parts:
            traj_parts[key].append(batch_series[key])

        print(f"[baseline] {total_cycles} cycles")

    trajectory = {
        key: (np.concatenate(parts) if parts else np.array([]))
        for key, parts in traj_parts.items()
    }
    return sol, trajectory


def main():
    base_param = make_param()

    print("\n=== baseline: no solvent consumption (OKane2022, calibrated SEI kinetics) ===")
    _, traj_base = run_baseline_no_dryout(base_param)

    conditions = []
    for r_eres in (0.0, 0.06, 0.09):
        print(f"\n=== dry-out wrapper: r_eres={r_eres:.0%} ===")
        _, history, ledger, traj = run_ec_dryout_degradation(
            base_param, FORMATION_EXP, AGEING_CYCLE,
            update_every_n_cycles=UPDATE_EVERY_N_CYCLES,
            max_total_cycles=MAX_TOTAL_CYCLES, r_eres=r_eres,
            options=MODEL_OPTIONS, composite=False,
        )
        conditions.append((r_eres, traj, history))

    colors = {0.0: "tab:red", 0.06: "tab:orange", 0.09: "tab:green"}
    nominal_cap = base_param["Nominal cell capacity [A.h]"]

    # -------------------------------------------------------------------
    # Figure 1: capacity retention / LLI -- cf. paper Fig. 3(a-c)
    # (LAM omitted -- disabled in MODEL_OPTIONS, always zero here)
    # -------------------------------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))

    cap0 = traj_base["cap"][0]
    ax[0].plot(traj_base["thr"], 100 * traj_base["cap"] / cap0, "^-", ms=4,
               color="tab:gray", label="no solvent consumption")
    for r_eres, traj, _ in conditions:
        ax[0].plot(traj["thr"], 100 * traj["cap"] / cap0, "o-", ms=4,
                   color=colors[r_eres], label=f"{r_eres:.0%} extra electrolyte")
    ax[0].set_xlabel("Throughput capacity [A.h]")
    ax[0].set_ylabel("Discharge capacity retention [%]")
    ax[0].set_title("(a) Capacity retention")
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3)

    ax[1].plot(traj_base["thr"], 100 * traj_base["LLI"] / nominal_cap, "^-", ms=4,
               color="tab:gray", label="no solvent consumption")
    for r_eres, traj, _ in conditions:
        ax[1].plot(traj["thr"], 100 * traj["LLI"] / nominal_cap, "o-", ms=4,
                   color=colors[r_eres], label=f"{r_eres:.0%}")
    ax[1].set_xlabel("Throughput capacity [A.h]")
    ax[1].set_ylabel("Capacity loss to negative SEI [%]")
    ax[1].set_title("(b) LLI (SEI)")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)

    fig.suptitle(
        "OKane2022 dry-out demonstration (single electrode, pore buffering OFF, "
        f"calibrated SEI kinetics, {MAX_TOTAL_CYCLES}-cycle protocol) -- cf. paper Fig. 3(a-b)"
    )
    plt.tight_layout()
    out1 = os.path.join(SCRIPT_DIR, "ruihe_dryout_validation_fig3.png")
    plt.savefig(out1, dpi=150)
    print(f"\nSaved: {out1}")

    # -------------------------------------------------------------------
    # Figure 2: R_dry / R_Li / bulk EC concentration vs cycle -- cf. paper Fig. A-3(d-f)/A-6
    # -------------------------------------------------------------------
    fig2, ax2 = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
    for r_eres, traj, history in conditions:
        cyc_x = UPDATE_EVERY_N_CYCLES * np.arange(1, len(history) + 1)
        ax2[0].plot(cyc_x, [h["R_dry"] for h in history], "o-", ms=4,
                    color=colors[r_eres], label=f"{r_eres:.0%}")
        ax2[1].plot(cyc_x, [h["R_Li"] for h in history], "o-", ms=4, color=colors[r_eres])
        ax2[2].plot(cyc_x, [h["c_EC_new"] for h in history], "o-", ms=4, color=colors[r_eres])
    ax2[0].set_ylabel("R_dry\n(electrolyte dry-out ratio)")
    ax2[0].set_title("(d) Electrode-area / dry-out ratio")
    ax2[0].legend(fontsize=8)
    ax2[0].grid(alpha=0.3)
    ax2[1].set_ylabel("R_Li\n(Li+ mixing ratio)")
    ax2[1].set_title("(e) Lithium-ion concentration change ratio")
    ax2[1].grid(alpha=0.3)
    ax2[2].set_ylabel("Bulk EC\nconcentration [mol/m3]")
    ax2[2].set_xlabel("Ageing cycle number")
    ax2[2].set_title("(f) EC concentration")
    ax2[2].grid(alpha=0.3)
    fig2.suptitle(
        "Electrolyte dry-out description (OKane2022, calibrated SEI kinetics, "
        f"{MAX_TOTAL_CYCLES}-cycle protocol) -- cf. paper Fig. A-3(d-f) / A-6"
    )
    plt.tight_layout()
    out2 = os.path.join(SCRIPT_DIR, "ruihe_dryout_validation_ratios.png")
    plt.savefig(out2, dpi=150)
    print(f"Saved: {out2}")

    print("\n--- summary ---")
    print(f"baseline (no solvent consumption): final capacity retention = "
          f"{100 * traj_base['cap'][-1] / cap0:.2f}%")
    for r_eres, traj, history in conditions:
        total_dnEC = sum(h["dn_EC"] for h in history)
        print(f"r_eres={r_eres:.0%}: final capacity retention = "
              f"{100 * traj['cap'][-1] / cap0:.2f}%  "
              f"final R_dry={history[-1]['R_dry']:.4f}  min R_dry={min(h['R_dry'] for h in history):.4f}  "
              f"total EC consumed={total_dnEC:.3e} mol  "
              f"final c_EC={history[-1]['c_EC_new']:.1f} mol/m3")


if __name__ == "__main__":
    main()
