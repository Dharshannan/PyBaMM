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

SEI submodel coupling mode (sei_ec_coupling, see implementation_plan.md
section 9 for the full rationale): this wrapper's R_EC ratio needs a
scalar Parameter representing "how much EC/solvent is available", to
write back into ParameterValues each batch -- which PyBaMM parameter that
is depends on which SEI submodel the caller's `options["SEI"]` uses:
  - "solvent_diffusion" (default, preserves the originally-validated
    behaviour): options["SEI"] == "solvent-diffusion limited",
    j_sei = -F*D_sol*c_sol/L_sei -- a PURE diffusion law with no reaction-
    kinetic term at all. Self-limiting/negative-feedback by construction
    (as c_sol falls, j_sei falls proportionally, with nothing to offset
    it) -- this is a deliberate, known property of this submodel (Marquis
    thesis eq. 5.91), not a wrapper bug, but it means a cell run under
    this submodel alone will NOT show an accelerating-into-a-knee
    capacity trajectory from SEI growth alone; the knee behaviour this
    project's tuned recipes rely on comes from the SEPARATE stress-driven
    LAM/cracking pathway, not from SEI growth's own dynamics under this
    submodel. Writes "{pref}Bulk solvent concentration [mol.m-3]".
  - "ec_reaction": options["SEI"] == "ec reaction limited" (or "...
    (asymmetric)") -- Yang et al. 2017 (sei_growth.py's "ec reaction
    limited" branch), a self-consistent linear solve for j and the local
    EC surface concentration given a Butler-Volmer-style REACTION-rate
    term (k_exp = k_sei * exp(-alpha_SEI*F_RT*eta_SEI), genuinely
    overpotential/resistance-driven, same qualitative character as
    "reaction limited"'s own j0_sei*exp(...) law) in series with EC
    diffusion to the reaction site. Unlike "solvent-diffusion limited",
    this retains real reaction kinetics that CAN accelerate as cell
    resistance grows (i.e. can still produce knee-like behaviour), while
    remaining genuinely coupled to EC availability via a real parameter
    (c_ec_0) this wrapper can update -- unlike plain "reaction limited"
    (j_sei = -j0_sei*exp(...) only), which has NO concentration term of
    any kind and therefore nothing for a dry-out mechanism to attach to
    without inventing new core-model physics (explicitly out of scope --
    "implement changes into the wrapper, not the core scripts"). Writes
    "{pref}EC initial concentration in electrolyte [mol.m-3]" (the c_ec_0
    parameter -- see lithium_ion_parameters.py:446-448) instead.
Both modes share IDENTICAL mass-balance mechanics (EC consumption
stoichiometry, pore-volume tracking, reservoir refill, R_dry/R_Li) --
only the ONE parameter key that receives the updated EC concentration
each batch differs. si_gr_expansion.py already defines
"{Primary,Secondary}: EC initial concentration in electrolyte [mol.m-3]"
/ "EC diffusivity [m2.s-1]" / "SEI kinetic rate constant [m.s-1]"
(:824-826, :841-843) -- "ec_reaction" mode needs no new parameters, only
this wrapper's own key-selection logic and the caller's own recipe
recalibration (separate concern, tracked in
degradation_test_matrix_plan.md, not this wrapper).
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

SEI_EC_COUPLING_MODES = ("solvent_diffusion", "ec_reaction")


def _ec_concentration_param_name(pref, sei_ec_coupling):
    """The one PyBaMM parameter key that represents "how much EC/solvent is
    available" for the given SEI submodel family -- see the module
    docstring's "SEI submodel coupling mode" section for the full
    rationale. `pref` is "Primary: "/"Secondary: " (composite) or ""
    (single-phase), matching this file's existing pref convention."""
    if sei_ec_coupling == "ec_reaction":
        return f"{pref}EC initial concentration in electrolyte [mol.m-3]"
    if sei_ec_coupling == "solvent_diffusion":
        return f"{pref}Bulk solvent concentration [mol.m-3]"
    raise ValueError(
        f"sei_ec_coupling={sei_ec_coupling!r} not recognised, "
        f"must be one of {SEI_EC_COUPLING_MODES}"
    )


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
    def __init__(
        self, param, r_eres, area_param_name=AREA_PARAM, composite=True,
        sei_ec_coupling="solvent_diffusion", ec_consumption_scale=1.0,
        ec_molar_volume_scale=1.0,
    ):
        self.param = param
        self.area_param_name = area_param_name
        self.composite = composite
        self.sei_ec_coupling = sei_ec_coupling
        self.ec_consumption_scale = ec_consumption_scale
        self.ec_molar_volume_scale = ec_molar_volume_scale
        self.area_scale = 1.0            # cumulative product of R_dry
        self.n_parallel_0 = param[area_param_name]
        pref = "Primary: " if composite else ""
        c_ec_key = _ec_concentration_param_name(pref, sei_ec_coupling)
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
        # Eq. 1: 2Li+ + 2EC + 2e- -> SEI + gas), scaled by ec_consumption_scale
        # (default 1.0 = paper's exact stoichiometry). A calibration lever, not
        # part of the paper: si_gr_expansion's default EC pool (V_eJR * c_EC0)
        # can be small enough that the paper's exact stoichiometry exhausts it
        # within a fraction of the recipe's own knee/lifetime (observed:
        # degradation_test_matrix/test_dry_out's first coupling test consumed
        # ~0.0244 mol against a ~0.0244 mol BOL pool by cycle 240) -- for
        # "ec reaction limited" specifically, where j_sei is directly
        # proportional to the fed-back c_EC, premature exhaustion self-limits
        # the SEI/LLI channel entirely rather than producing a representative
        # dry-out trajectory. Reducing this scale stretches how many cycles
        # the same BOL pool lasts, without changing anything else.
        dn_EC = delta_Q_SEI_Ah * 3600.0 / F * self.ec_consumption_scale

        # Step 2: pore & electrolyte volumes
        V_pore_new = A_cell * (L_n * eps_n_avg_new + L_s * eps_s0 + L_p * eps_p0)
        # ec_molar_volume_scale (default 1.0 = EC's real physical molar volume) is a
        # SEPARATE lever from "SEI partial molar volume": that PyBaMM parameter also drives
        # the actual simulated porosity-decline rate (hence the recipe's own LAM/knee
        # physics), so shrinking it to get R_dry moving distorts the tuned recipe itself
        # (confirmed: dryout_molarvol_sweep_test.py found the knee already eroded by
        # SEI_MOLAR_VOLUME_DIV~2, well before R_dry moved meaningfully). M_EC/RHO_EC here
        # are plain module constants used ONLY in this ledger's own EC-volume bookkeeping --
        # scaling them changes how much electrolyte volume is deemed consumed per mole of
        # EC reacted WITHOUT touching PyBaMM's own porosity/LAM/knee simulation at all.
        dV_EC = dn_EC * M_EC / RHO_EC * self.ec_molar_volume_scale
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
def _cracks_sei_capacity_loss_Ah(sol, key):
    """"Loss of capacity to ... SEI on cracks [A.h]" at end of `sol`, or 0.0
    if that variable wasn't built for the given model options (e.g.
    "particle mechanics": "none", as in ruihe_dryout_validation.py's
    OKane2022 setup) -- degrades gracefully rather than raising."""
    try:
        return sol[key].entries[-1]
    except KeyError:
        return 0.0


def sei_capacity_loss_Ah(sol, composite=True):
    """Cumulative EC-consuming SEI capacity loss: "Loss of capacity to
    negative SEI [A.h]" (bulk SEI) PLUS "...SEI on cracks [A.h]" (crack-
    surface SEI), summed over both phases if composite=True (both phases
    share the SEI reaction, Eq. 1, and therefore both consume EC), or the
    single-phase variant otherwise (e.g. OKane2022, a plain single-
    electrode chemistry).

    The "SEI on cracks" term matters: crack-surface SEI growth is driven
    by the same j_sei reaction as bulk SEI and therefore consumes EC too,
    but every degradation_test_matrix recipe in this project runs with
    "SEI on cracks": "true" -- omitting this term would systematically
    under-count EC consumption (hence under-count dry-out severity)
    whenever cracking is active. Found when building
    degradation_test_matrix/test_dry_out/'s first coupling test; degrades
    to +0 (via _cracks_sei_capacity_loss_Ah) wherever that variable isn't
    built, so ruihe_dryout_validation.py's already-validated OKane2022
    behaviour (no particle mechanics -> no SEI-on-cracks) is unaffected."""
    if composite:
        return (
            sol["Loss of capacity to negative primary SEI [A.h]"].entries[-1]
            + sol["Loss of capacity to negative secondary SEI [A.h]"].entries[-1]
            + _cracks_sei_capacity_loss_Ah(
                sol, "Loss of capacity to negative primary SEI on cracks [A.h]")
            + _cracks_sei_capacity_loss_Ah(
                sol, "Loss of capacity to negative secondary SEI on cracks [A.h]")
        )
    return (
        sol["Loss of capacity to negative SEI [A.h]"].entries[-1]
        + _cracks_sei_capacity_loss_Ah(sol, "Loss of capacity to negative SEI on cracks [A.h]")
    )


def _cycle_ageing_leg_capacity(cyc):
    """(cap, discharge_end_thr, charge_end_thr) for the ageing-rate (>1 A
    mean current) discharge leg of this cycle, or (0.0, None, None) if not
    found. charge_end_thr (added for reversible-expansion-amplitude
    tracking, degradation_test_matrix/test_dry_out) is the throughput at
    the end of this same cycle's charge leg (mean current < -1e-3 A),
    matching this project's other scripts' cycle_ageing_leg convention."""
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
    # "Cell thickness change [m]" (reversible expansion amplitude) only exists if particle
    # mechanics/pore buffering produce a thickness-change output -- degrade gracefully to
    # all-zero rather than raising, same convention as the LAM guard above.
    try:
        tc_cell_full = sol["Cell thickness change [m]"].entries
    except KeyError:
        tc_cell_full = np.zeros_like(Qt_full)

    thr_list, cap_list, discharge_end_list, charge_end_list = [], [], [], []
    for cyc in sol.cycles:
        cap, discharge_end_thr, charge_end_thr = _cycle_ageing_leg_capacity(cyc)
        if cap > 0 and discharge_end_thr is not None:
            thr_list.append(discharge_end_thr)
            cap_list.append(cap)
            discharge_end_list.append(discharge_end_thr)
            charge_end_list.append(charge_end_thr)
    thr_arr = np.array(thr_list)
    cap_arr = np.array(cap_list)
    if thr_arr.size:
        LLI_arr = np.interp(thr_arr, Qt_full, LLI_full)
        LAM_neg_arr = np.interp(thr_arr, Qt_full, LAM_neg_full)
        LAM_pos_arr = np.interp(thr_arr, Qt_full, LAM_pos_full)
        discharge_end_arr = np.array(discharge_end_list)
        # charge_end_thr can be None if a cycle's charge leg wasn't found (e.g. the very
        # last, incomplete cycle of a batch) -- fall back to the discharge-end value
        # (amplitude 0 for that one point) rather than letting None reach np.interp.
        charge_end_arr = np.array(
            [c if c is not None else d for c, d in zip(charge_end_list, discharge_end_list)]
        )
        tc_charge_end = np.interp(charge_end_arr, Qt_full, tc_cell_full)
        tc_discharge_end = np.interp(discharge_end_arr, Qt_full, tc_cell_full)
        cell_amplitude_arr = tc_charge_end - tc_discharge_end
    else:
        LLI_arr = LAM_neg_arr = LAM_pos_arr = cell_amplitude_arr = np.array([])
    return dict(
        thr=thr_arr, cap=cap_arr, LLI=LLI_arr,
        LAM_neg=LAM_neg_arr, LAM_pos=LAM_pos_arr, cell_amplitude=cell_amplitude_arr,
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
    sei_ec_coupling="solvent_diffusion",
    ec_consumption_scale=1.0,
    ec_molar_volume_scale=1.0,
    soh_target=None,
    solver=None,
    var_pts=None,
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
        buffering ON, "physical" transition) if not given. `options["SEI"]`
        must be consistent with `sei_ec_coupling` (see below) -- this
        function does not cross-check the two, since `options` is caller-
        supplied and may already differ from MODEL_OPTIONS in other ways.
    composite : True for a composite (Primary:/Secondary:-prefixed)
        negative electrode chemistry (e.g. si_gr_expansion); False for a
        plain single-phase chemistry (e.g. OKane2022) -- controls which
        parameter/variable names are used for bulk solvent concentration
        and SEI capacity loss.
    sei_ec_coupling : "solvent_diffusion" (default, matches
        `options["SEI"] == "solvent-diffusion limited"`) or "ec_reaction"
        (matches `options["SEI"]` starting with "ec reaction limited") --
        selects which PyBaMM parameter the ledger's EC concentration state
        is read from/written to each batch. See the module docstring's
        "SEI submodel coupling mode" section and
        `_ec_concentration_param_name` for the full rationale: in short,
        "solvent_diffusion" is a pure, self-limiting diffusion law with no
        reaction-kinetic term (won't accelerate into a knee from SEI
        growth alone); "ec_reaction" keeps genuine Butler-Volmer-style
        reaction kinetics (can accelerate/produce a knee) while still
        exposing a real EC-concentration parameter (`c_ec_0`) this wrapper
        can drive dry-out through -- unlike plain "reaction limited",
        which has no such parameter at all.
    ec_consumption_scale : multiplier (default 1.0, the paper's exact 1:1:1
        stoichiometry) applied to dn_EC each batch -- see ECDryoutLedger.update's
        comment. A calibration lever for cases (e.g. "ec_reaction" mode with
        si_gr_expansion's default BOL EC pool) where the exact stoichiometry
        exhausts the entire EC pool well before the recipe's own knee/lifetime.
    ec_molar_volume_scale : multiplier (default 1.0 = EC's real physical molar
        volume, M_EC/RHO_EC) applied to dV_EC (electrolyte volume consumed per
        mole EC reacted) each batch. Distinct from recalibrating "SEI partial
        molar volume" (a real PyBaMM parameter that ALSO drives the recipe's
        own simulated porosity-decline/LAM/knee physics): degradation_test_
        matrix/test_dry_out/dryout_molarvol_sweep_test.py found that lever
        can't get R_dry moving without already eroding the tuned recipe's knee
        shape (both happen in the same divisor range). M_EC/RHO_EC are plain
        module constants used ONLY in this ledger's own EC-volume bookkeeping,
        so scaling them changes how much electrolyte volume is deemed consumed
        per mole reacted WITHOUT touching PyBaMM's own simulated physics at
        all -- a genuinely independent lever on R_dry.
    soh_target : if given (e.g. 0.5 for 50% SoH), the loop stops as soon as
        the running discharge-capacity ratio (this batch's last ageing-leg
        capacity / the very first recorded ageing-leg capacity) drops to or
        below this value, in addition to the max_total_cycles cap -- mirrors
        this project's other degradation_test_matrix scripts' SoH-based
        termination. Default None preserves the original behaviour (run
        exactly to max_total_cycles).
    solver, var_pts : forwarded to every pybamm.Simulation this function
        builds (formation AND every ageing batch). Both default to None,
        i.e. PyBaMM's own defaults -- preserves every previously-validated
        call site's behaviour (ruihe_dryout_validation.py never passed
        these). IMPORTANT for any caller that builds its own "baseline, no
        solvent consumption" comparison run with an explicit solver/var_pts
        (as degradation_test_matrix/test_dry_out's coupling test does): pass
        the SAME values here, or the two runs use different numerics
        (different spatial discretisation and/or solver tolerances) and
        aren't a controlled comparison -- this was found and fixed after
        that script's dry-out run showed a spurious single-cycle SoH
        dropout its differently-discretised baseline did not.

    Returns (final_solution, history, ledger, trajectory).
    `history` is a list of the dict returned by ECDryoutLedger.update()
    for each batch, in order. `trajectory` is a dict of numpy arrays
    (keys: thr, cap, LLI, LAM_neg, LAM_pos, cell_amplitude -- the last is
    reversible expansion amplitude, tc_cell at charge-end minus
    discharge-end per ageing cycle, zeros if "Cell thickness change [m]"
    isn't built for the given model options) giving the FULL-RUN, cycle-by-
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
        experiment=formation_experiment, solver=solver, var_pts=var_pts,
    )
    sol = sim0.solve(initial_soc=1.0)
    print("[ec_dryout] formation cycle solved.")

    ledger = ECDryoutLedger(
        param, r_eres, composite=composite, sei_ec_coupling=sei_ec_coupling,
        ec_consumption_scale=ec_consumption_scale,
        ec_molar_volume_scale=ec_molar_volume_scale,
    )
    ledger.initialise_from_formation_solution(sol)
    print(
        f"[ec_dryout] BOL: V_pore={ledger.V_pore:.4e} m3  V_eJR={ledger.V_eJR:.4e} m3  "
        f"V_eres={ledger.V_eres:.4e} m3 (r_eres={r_eres:.2%})"
    )

    Q_SEI_prev = sei_capacity_loss_Ah(sol, composite=composite)
    history = []
    traj_parts = dict(thr=[], cap=[], LLI=[], LAM_neg=[], LAM_pos=[], cell_amplitude=[])
    total_cycles = 0
    cap0 = None
    while total_cycles < max_total_cycles:
        batch_exp = pybamm.Experiment(
            [ageing_cycle_steps for _ in range(update_every_n_cycles)]
        )
        model = pybamm.lithium_ion.DFN(options)
        sim = pybamm.Simulation(
            model, parameter_values=param, experiment=batch_exp,
            solver=solver, var_pts=var_pts,
        )

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
                _ec_concentration_param_name("Primary: ", sei_ec_coupling): step["c_EC_new"],
                _ec_concentration_param_name("Secondary: ", sei_ec_coupling): step["c_EC_new"],
            }
        else:
            ec_updates = {
                _ec_concentration_param_name("", sei_ec_coupling): step["c_EC_new"],
            }
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

        if soh_target is not None and batch_series["cap"].size:
            if cap0 is None:
                cap0 = traj_parts["cap"][0][0]
            soh_now = batch_series["cap"][-1] / cap0
            print(f"[ec_dryout] batch {len(history)}: SoH={soh_now:.3f}")
            if soh_now <= soh_target:
                print(
                    f"[ec_dryout] reached soh_target ({soh_target:.0%}) at "
                    f"{total_cycles} cycles -- stopping."
                )
                break

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
