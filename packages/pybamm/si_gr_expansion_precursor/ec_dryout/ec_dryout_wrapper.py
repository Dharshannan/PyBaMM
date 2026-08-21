"""
ec_dryout_wrapper.py -- standalone electrolyte dry-out / solvent-consumption
wrapper around PyBaMM, implementing Li et al. 2022 (J. Electrochem. Soc.
169 060516) "Modelling Solvent Consumption from SEI Layer Growth in
Lithium-Ion Batteries" as a batch-restart post-processing loop. See
si_gr_expansion_precursor/ec_dryout/implementation_plan.md for the full
derivation, API research and validation history behind this script.

No PyBaMM core changes. Every UPDATE_EVERY_N_CYCLES cycles, the just-
completed batch's solution is used to compute three ratios (R_dry, R_EC,
R_Li) per the paper's Eqs. 26/30/32, applied as:
  - R_dry -> shrinks "Number of electrodes connected in parallel to make
    a cell" (a pure area multiplier, safe for any current-collector mode --
    see plan doc section 3.5) for the NEXT batch's ParameterValues.
  - R_EC  -> sets "{Primary,Secondary}: Bulk solvent concentration
    [mol.m-3]" directly (it is a plain scalar Parameter in PyBaMM's
    "solvent-diffusion limited" SEI submodel, not a PDE state).
  - R_Li  -> rescales the transplanted electrolyte Li+ concentration
    state ("Porosity times concentration [mol.m-3]") in place on
    `solution.last_state.y`, exact and validated (spike_li_rescale.py;
    the state has zero reference offset, so a pure multiplicative
    rescale is correct).
Everything else (particle concentrations, SEI amount, active-material
volume fractions, cracking, pore-buffering's own algebraic outputs, ...)
is carried over unmodified via `sim.solve(starting_solution=...)`, which
this project already uses elsewhere and which was confirmed (this
implementation) to tolerate a DIFFERENT ParameterValues object between
batches -- it is not restricted to identical-parameter continuation as
originally assumed.

Compatible with "pore buffering": "true" with NO special-casing needed:
the pore-buffering submodel (reaction_driven_porosity.py) introduces no
new state variables of its own -- it only recomputes porosity/thickness
outputs from the SAME underlying SEI-thickness state this wrapper already
transplants, and neither electrode area nor bulk EC concentration enter
its partition equations. See implementation_plan.md section 4.4.
"""
import numpy as np
import pybamm

# ---------------------------------------------------------------------------
# EC (ethylene carbonate) physical constants -- NOT PyBaMM Parameters (the
# repo has no "solvent density"/"solvent molar mass" parameter; these are
# fixed physical constants of the specific solvent species, kept as plain
# module constants, matching how the paper's own Table A-II treats them as
# external to the DFN parameter table).
# ---------------------------------------------------------------------------
M_EC = 88.062e-3     # EC molar mass [kg/mol]
RHO_EC = 1321.0      # EC density [kg/m3] (~1.32 g/cm3)

MODEL_OPTIONS = {
    "particle phases": ("2", "1"),
    "open-circuit potential": (("single", "current sigmoid"), "single"),
    "SEI": "solvent-diffusion limited",   # required: matches the paper's
                                           # j_SEI = -F*D_EC*c_EC/L_SEI exactly
    "SEI porosity change": "true",        # gives V_pore(t) via porosity --
                                           # no need to re-derive from dn_EC
    "particle mechanics": ("swelling and cracking", "swelling only"),
    "SEI on cracks": "true",
    "loss of active material": "stress-driven",
    "pore buffering": "true",
    "pore buffering transition": "physical",
}

AREA_PARAM = "Number of electrodes connected in parallel to make a cell"
POROSITY_TIMES_CONC_NAME = "Porosity times concentration [mol.m-3]"


# ---------------------------------------------------------------------------
# R_Li state rescale -- validated mechanism, see spike_li_rescale.py
# ---------------------------------------------------------------------------
def _find_y_slice(solution_last_state, var_name):
    lm = solution_last_state.all_models[-1]
    for var, slices in lm.y_slices.items():
        if getattr(var, "name", None) == var_name:
            return slices[0]
    raise KeyError(
        f"{var_name!r} not found in y_slices -- PyBaMM's internal state "
        f"representation may have changed; re-run spike_li_rescale.py"
    )


def apply_li_rescale(last_state, R_Li):
    """In-place rescale of the transplanted electrolyte Li+ concentration
    state by R_Li (a no-op if R_Li == 1, i.e. whenever no reservoir liquid
    mixed in this interval). Mutates and returns last_state."""
    if R_Li == 1.0:
        return last_state
    sl = _find_y_slice(last_state, POROSITY_TIMES_CONC_NAME)
    last_state.y[sl] = last_state.y[sl] * R_Li
    return last_state


# ---------------------------------------------------------------------------
# External ledger: everything the paper treats as independent of PyBaMM's
# own state (JR/reservoir electrolyte volumes, bulk EC concentration,
# cumulative electrode-area shrink factor). Plain Python floats, updated
# once per batch from that batch's solved output.
# ---------------------------------------------------------------------------
class ECDryoutLedger:
    def __init__(self, param, r_eres, area_param_name=AREA_PARAM, composite=True):
        self.param = param
        self.area_param_name = area_param_name
        self.composite = composite
        self.area_scale = 1.0            # cumulative product of R_dry
        self.n_parallel_0 = param[area_param_name]
        c_ec_key = (
            "Primary: Bulk solvent concentration [mol.m-3]" if composite
            else "Bulk solvent concentration [mol.m-3]"
        )
        self.c_EC = param[c_ec_key]
        self.c_EC0 = self.c_EC   # BOL EC concentration (reservoir liquid's
                                  # concentration, per the paper's assumption
                                  # that reservoir electrolyte matches the
                                  # initial JR concentration -- Eq. 31)
        self.r_eres = r_eres
        self.V_pore = None
        self.V_eJR = None
        self.V_eres = None

    def initialise_from_formation_solution(self, sol0):
        """Call once, right after the formation solve, before any ageing
        batches. Reads BOL geometry/porosity from the SOLVED state
        (robust to porosity being a FunctionParameter, unlike trying to
        hand-compute eps0 from raw Parameters)."""
        p = self.param
        L_n = p["Negative electrode thickness [m]"]
        L_s = p["Separator thickness [m]"]
        L_p = p["Positive electrode thickness [m]"]
        eps_n0 = sol0["X-averaged negative electrode porosity"].entries[0]
        eps_s0 = p["Separator porosity"]
        eps_p0 = p["Positive electrode porosity"]
        A_cell0 = (
            p["Electrode width [m]"] * p["Electrode height [m]"] * self.n_parallel_0
        )
        self.V_pore = A_cell0 * (L_n * eps_n0 + L_s * eps_s0 + L_p * eps_p0)
        self.V_eJR = self.V_pore          # fully wetted at BOL
        self.V_eres = self.r_eres * self.V_eJR

    def update(self, delta_Q_SEI_Ah, eps_n_avg_new, c_Li_avg_end, F):
        """Advance the ledger by one batch. Returns a dict of the new
        absolute values ready to feed into the next batch's ParameterValues
        / initial-condition rescale.

        delta_Q_SEI_Ah : this batch's (end - start) cumulative
            "Loss of capacity to negative {primary,secondary} SEI [A.h]"
            (summed over phases -- both phases share the SEI reaction).
        eps_n_avg_new : "X-averaged negative electrode porosity" at the
            end of this batch.
        c_Li_avg_end : "X-averaged electrolyte concentration [mol.m-3]"
            at the end of this batch (pre-mixing reference for R_Li).
        F : Faraday constant [C/mol].
        """
        p = self.param
        L_n = p["Negative electrode thickness [m]"]
        L_s = p["Separator thickness [m]"]
        L_p = p["Positive electrode thickness [m]"]
        eps_s0 = p["Separator porosity"]
        eps_p0 = p["Positive electrode porosity"]
        A_cell = (
            p["Electrode width [m]"] * p["Electrode height [m]"]
            * self.n_parallel_0 * self.area_scale
        )

        # Step 1: EC consumed this interval (1:1:1 Li:EC:e- stoichiometry,
        # Eq. 1: 2Li+ + 2EC + 2e- -> SEI + gas)
        dn_EC = delta_Q_SEI_Ah * 3600.0 / F

        # Step 2: pore & electrolyte volumes
        V_pore_new = A_cell * (L_n * eps_n_avg_new + L_s * eps_s0 + L_p * eps_p0)
        dV_EC = dn_EC * M_EC / RHO_EC
        V_eJR_predryout = self.V_eJR - dV_EC
        shortfall = max(V_pore_new - V_eJR_predryout, 0.0)

        # Step 3: reservoir refill, if any
        dV_add = min(shortfall, self.V_eres)
        V_eres_new = self.V_eres - dV_add
        V_eJR_new = V_eJR_predryout + dV_add

        # Step 4: dry-out ratio -> cumulative area shrink
        R_dry = float(np.clip(V_eJR_new / V_pore_new, 0.0, 1.0))
        self.area_scale *= R_dry

        # Step 5: bulk EC concentration update
        n_EC_old = self.c_EC * self.V_eJR
        n_EC_new = n_EC_old - dn_EC + self.c_EC0 * dV_add
        c_EC_new = max(n_EC_new / V_eJR_new, 0.0)

        # Step 6: Li+ ratio (only non-trivial when dV_add > 0)
        if dV_add > 0:
            c_Li0 = p["Initial concentration in electrolyte [mol.m-3]"]
            R_Li = (
                c_Li_avg_end * self.V_eJR + c_Li0 * dV_add
            ) / (c_Li_avg_end * V_eJR_new)
        else:
            R_Li = 1.0

        # commit
        self.V_pore, self.V_eJR, self.V_eres, self.c_EC = (
            V_pore_new, V_eJR_new, V_eres_new, c_EC_new,
        )

        return dict(
            R_dry=R_dry, R_Li=R_Li, c_EC_new=c_EC_new,
            n_parallel_new=self.n_parallel_0 * self.area_scale,
            V_pore=V_pore_new, V_eJR=V_eJR_new, V_eres=V_eres_new,
            dn_EC=dn_EC, dV_add=dV_add,
        )


# ---------------------------------------------------------------------------
# batch-restart driver
# ---------------------------------------------------------------------------
def sei_capacity_loss_Ah(sol, composite=True):
    """Cumulative "Loss of capacity to negative SEI [A.h]" -- summed over
    both phases if composite=True (both phases share the SEI reaction,
    Eq. 1, and therefore both consume EC), or the single-phase variant
    otherwise (e.g. OKane2022, a plain single-electrode chemistry)."""
    if composite:
        return (
            sol["Loss of capacity to negative primary SEI [A.h]"].entries[-1]
            + sol["Loss of capacity to negative secondary SEI [A.h]"].entries[-1]
        )
    return sol["Loss of capacity to negative SEI [A.h]"].entries[-1]


def _cycle_ageing_leg_capacity(cyc):
    """(cap, discharge_end_thr) for the ageing-rate (>1 A mean current)
    discharge leg of this cycle, or (0.0, None) if not found."""
    for step in getattr(cyc, "steps", []):
        try:
            I = step["Current [A]"].entries
            q = step["Discharge capacity [A.h]"].entries
            thr_step = step["Throughput capacity [A.h]"].entries
        except (KeyError, TypeError, AttributeError):
            continue
        if I.size >= 2 and np.mean(I) > 1.0 and q.size >= 2:
            return float(q[-1] - q[0]), float(thr_step[-1])
    return 0.0, None


def extract_batch_series(sol, composite=True):
    """Per-cycle discharge capacity + degradation-mode values, for the
    cycles contained in THIS SINGLE BATCH's `sol` alone.

    IMPORTANT: because run_ec_dryout_degradation continues each batch via
    `solution.last_state` (a single-timestep reduction, needed for the
    R_Li rescale -- see apply_li_rescale), each batch's own `sol` does
    NOT retain earlier batches' cycle history the way this project's other
    scripts' `starting_solution=full_solution` pattern does (verified:
    passing `.last_state` gives a `sol.cycles` covering only that batch's
    own cycles, with `sol["Throughput capacity [A.h]"]` starting from
    that batch's own t=0, not the run's). The underlying STATE is still
    correctly carried forward (confirmed: end-of-batch cumulative values
    like SEI capacity loss are physically correct), so this only matters
    for reconstructing a FULL-RUN trajectory -- which is exactly why
    run_ec_dryout_degradation calls this function once per batch and
    accumulates the results into its returned `trajectory` dict, rather
    than callers trying to read a full history off the final returned
    `sol` (which would silently only reflect the last batch)."""
    Qt_full = sol["Throughput capacity [A.h]"].entries
    if composite:
        LLI_full = (
            sol["Loss of capacity to negative primary SEI [A.h]"].entries
            + sol["Loss of capacity to negative secondary SEI [A.h]"].entries
        )
    else:
        LLI_full = sol["Loss of capacity to negative SEI [A.h]"].entries
    # LAM variables only exist if "loss of active material" is enabled --
    # minimal test configs (e.g. the OKane2022 dry-out-only setup) turn
    # this off entirely, so degrade gracefully to all-zero rather than
    # raising a KeyError.
    try:
        LAM_neg_full = sol["Loss of active material in negative electrode [%]"].entries
        LAM_pos_full = sol["Loss of active material in positive electrode [%]"].entries
    except KeyError:
        LAM_neg_full = np.zeros_like(Qt_full)
        LAM_pos_full = np.zeros_like(Qt_full)

    thr_list, cap_list = [], []
    for cyc in sol.cycles:
        cap, thr = _cycle_ageing_leg_capacity(cyc)
        if cap > 0 and thr is not None:
            thr_list.append(thr)
            cap_list.append(cap)
    thr_arr = np.array(thr_list)
    cap_arr = np.array(cap_list)
    if thr_arr.size:
        LLI_arr = np.interp(thr_arr, Qt_full, LLI_full)
        LAM_neg_arr = np.interp(thr_arr, Qt_full, LAM_neg_full)
        LAM_pos_arr = np.interp(thr_arr, Qt_full, LAM_pos_full)
    else:
        LLI_arr = LAM_neg_arr = LAM_pos_arr = np.array([])
    return dict(
        thr=thr_arr, cap=cap_arr, LLI=LLI_arr,
        LAM_neg=LAM_neg_arr, LAM_pos=LAM_pos_arr,
    )


def run_ec_dryout_degradation(
    base_param,
    formation_experiment,
    ageing_cycle_steps,
    update_every_n_cycles,
    max_total_cycles,
    r_eres,
    options=None,
    apply_li_correction=True,
    composite=True,
):
    """Top-level batch-restart loop.

    base_param : pybamm.ParameterValues, NOT yet containing the dry-out
        wrapper's per-batch overrides (those get patched in fresh each
        batch as the loop runs).
    formation_experiment : pybamm.Experiment for the initial (BOL) solve.
    ageing_cycle_steps : the per-cycle experiment step tuple, e.g.
        ("Discharge at 1C until 2.5 V", "Charge at 0.3C until 4.2 V",
         "Hold at 4.2 V until C/100").
    update_every_n_cycles : batch size -- how many ageing cycles run
        between each EC-dryout ledger update (paper's own protocol uses
        ~16; this project's other scripts default to 20-50).
    r_eres : ratio of initial extra (reservoir) electrolyte to initial
        JR electrolyte volume (paper's Eq. 33; sweep 0.0, 0.06, 0.09 to
        reproduce the paper's three cases).
    options : model options dict; defaults to MODEL_OPTIONS (pore
        buffering ON, "physical" transition) if not given.
    composite : True for a composite (Primary:/Secondary:-prefixed)
        negative electrode chemistry (e.g. si_gr_expansion); False for a
        plain single-phase chemistry (e.g. OKane2022) -- controls which
        parameter/variable names are used for bulk solvent concentration
        and SEI capacity loss.

    Returns (final_solution, history, ledger, trajectory).
    `history` is a list of the dict returned by ECDryoutLedger.update()
    for each batch, in order. `trajectory` is a dict of numpy arrays
    (keys: thr, cap, LLI, LAM_neg, LAM_pos) giving the FULL-RUN, cycle-by-
    cycle degradation trajectory, accumulated batch-by-batch as the loop
    runs -- use this for full-run plotting, NOT `final_solution.cycles`
    (which only covers the last batch; see extract_batch_series's
    docstring for why).
    """
    options = dict(options or MODEL_OPTIONS)
    param = base_param.copy()
    F = float(pybamm.constants.F.evaluate())

    sim0 = pybamm.Simulation(
        pybamm.lithium_ion.DFN(options), parameter_values=param,
        experiment=formation_experiment,
    )
    sol = sim0.solve(initial_soc=1.0)
    print("[ec_dryout] formation cycle solved.")

    ledger = ECDryoutLedger(param, r_eres, composite=composite)
    ledger.initialise_from_formation_solution(sol)
    print(
        f"[ec_dryout] BOL: V_pore={ledger.V_pore:.4e} m3  V_eJR={ledger.V_eJR:.4e} m3  "
        f"V_eres={ledger.V_eres:.4e} m3 (r_eres={r_eres:.2%})"
    )

    Q_SEI_prev = sei_capacity_loss_Ah(sol, composite=composite)
    history = []
    traj_parts = dict(thr=[], cap=[], LLI=[], LAM_neg=[], LAM_pos=[])
    total_cycles = 0
    while total_cycles < max_total_cycles:
        batch_exp = pybamm.Experiment(
            [ageing_cycle_steps for _ in range(update_every_n_cycles)]
        )
        model = pybamm.lithium_ion.DFN(options)
        sim = pybamm.Simulation(model, parameter_values=param, experiment=batch_exp)

        last_state = sol.last_state
        if apply_li_correction and history:
            apply_li_rescale(last_state, history[-1]["R_Li"])

        sol = sim.solve(starting_solution=last_state)
        total_cycles += update_every_n_cycles

        batch_series = extract_batch_series(sol, composite=composite)
        for key in traj_parts:
            traj_parts[key].append(batch_series[key])

        Q_SEI_new = sei_capacity_loss_Ah(sol, composite=composite)
        delta_Q_SEI = Q_SEI_new - Q_SEI_prev
        Q_SEI_prev = Q_SEI_new

        eps_n_avg_new = sol["X-averaged negative electrode porosity"].entries[-1]
        c_Li_avg_end = sol["X-averaged electrolyte concentration [mol.m-3]"].entries[-1]
        step = ledger.update(delta_Q_SEI, eps_n_avg_new, c_Li_avg_end, F)
        history.append(step)

        if composite:
            ec_updates = {
                "Primary: Bulk solvent concentration [mol.m-3]": step["c_EC_new"],
                "Secondary: Bulk solvent concentration [mol.m-3]": step["c_EC_new"],
            }
        else:
            ec_updates = {"Bulk solvent concentration [mol.m-3]": step["c_EC_new"]}
        param.update(
            {**ec_updates, AREA_PARAM: step["n_parallel_new"]},
            check_already_exists=False,
        )

        print(
            f"[ec_dryout] batch {len(history)}: {total_cycles} cycles  "
            f"R_dry={step['R_dry']:.4f}  R_Li={step['R_Li']:.4f}  "
            f"c_EC={step['c_EC_new']:.1f} mol/m3  V_eres={step['V_eres']:.3e} m3  "
            f"dn_EC={step['dn_EC']:.3e} mol"
        )

    trajectory = {
        key: (np.concatenate(parts) if parts else np.array([]))
        for key, parts in traj_parts.items()
    }
    return sol, history, ledger, trajectory


if __name__ == "__main__":
    param = pybamm.ParameterValues("si_gr_expansion")
    formation_exp = pybamm.Experiment(
        [
            "Discharge at 0.1C until 2.5 V",
            "Charge at C/3 until 4.2 V",
            "Hold at 4.2 V until C/100",
        ]
    )
    ageing_cycle = (
        "Discharge at 1C until 2.5 V",
        "Charge at 0.3C until 4.2 V",
        "Hold at 4.2 V until C/100",
    )
    sol, history, ledger, trajectory = run_ec_dryout_degradation(
        param, formation_exp, ageing_cycle,
        update_every_n_cycles=20, max_total_cycles=60,
        r_eres=0.06,   # 6% extra electrolyte, matching one of the paper's cases
    )
    print(f"\nFinal: {len(history)} batches, {sum(h['dn_EC'] for h in history):.3e} "
          f"mol EC consumed total, final c_EC={history[-1]['c_EC_new']:.1f} mol/m3, "
          f"cumulative area scale={ledger.area_scale:.4f}")
    print(f"Full-run trajectory: {trajectory['thr'].size} cycles, "
          f"capacity {trajectory['cap'][0]:.4f} -> {trajectory['cap'][-1]:.4f} A.h")
