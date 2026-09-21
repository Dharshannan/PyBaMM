#
# Class for varying active material volume fraction
#
import pybamm

from .base_active_material import BaseModel


class LossActiveMaterial(BaseModel):
    """Submodel for varying active material volume fraction from :footcite:t:`Ai2019`
    and :footcite:t:`Reniers2019`.

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    domain : str
        The domain of the model either 'Negative' or 'Positive'
    options : dict
        Additional options to pass to the model
    x_average : bool
        Whether to use x-averaged variables (SPM, SPMe, etc) or full variables (DFN)

    """

    def __init__(self, param, domain, options, x_average, phase):
        super().__init__(param, domain, options=options, phase=phase)
        pybamm.citations.register("Reniers2019")
        self.x_average = x_average

    def get_fundamental_variables(self):
        domain, Domain = self.domain_Domain
        phase = self.phase_name

        if self.x_average is True:
            eps_solid_xav = pybamm.Variable(
                f"X-averaged {domain} electrode {phase}active material volume fraction",
                domain="current collector",
            )
            eps_solid = pybamm.PrimaryBroadcast(eps_solid_xav, f"{domain} electrode")
        else:
            eps_solid = pybamm.Variable(
                f"{Domain} electrode {phase}active material volume fraction",
                domain=f"{domain} electrode",
                auxiliary_domains={"secondary": "current collector"},
            )
        variables = self._get_standard_active_material_variables(eps_solid)
        lli_due_to_lam = pybamm.Variable(
            f"Loss of lithium due to loss of {phase}active material "
            f"in {domain} electrode [mol]"
        )

        variables.update(
            {
                f"Loss of lithium due to loss of {phase}active material "
                f"in {domain} electrode [mol]": lli_due_to_lam
            }
        )

        # Porosity-isolation gate as a genuine relaxation STATE, not an
        # algebraic snapshot of porosity -- see the isolation-gate-lag
        # investigation this session: the algebraic gate has no timescale
        # of its own, so tuning its shape/onset/magnitude could never
        # decouple "how wide the post-knee transition is" from "when it
        # starts" (both inherited entirely from porosity's own collapse
        # rate). This state relaxes toward the same algebraic target
        # (isolation_gate_target, computed in get_coupled_variables) with
        # its own time constant tau_LAM_iso, giving the transition width
        # an independent lever for the first time.
        lam_option = getattr(getattr(self.options, domain), self.phase)[
            "loss of active material"
        ]
        if "porosity" in lam_option:
            if self.x_average is True:
                isolation_state_xav = pybamm.Variable(
                    f"X-averaged {domain} electrode {phase}"
                    "porosity-isolation gate state",
                    domain="current collector",
                )
                isolation_state = pybamm.PrimaryBroadcast(
                    isolation_state_xav, f"{domain} electrode"
                )
            else:
                isolation_state = pybamm.Variable(
                    f"{Domain} electrode {phase}porosity-isolation gate state",
                    domain=f"{domain} electrode",
                    auxiliary_domains={"secondary": "current collector"},
                )
            variables.update(
                {
                    f"{Domain} electrode {phase}"
                    "porosity-isolation gate state": isolation_state,
                    f"X-averaged {domain} electrode {phase}"
                    "porosity-isolation gate state": pybamm.x_average(isolation_state),
                }
            )

        return variables

    def get_coupled_variables(self, variables):
        domain, Domain = self.domain_Domain
        phase_name = self.phase_name

        deps_solid_dt = 0
        lam_option = getattr(getattr(self.options, domain), self.phase)[
            "loss of active material"
        ]
        if "stress" in lam_option:
            # obtain the rate of loss of active materials (LAM) by stress
            # This is loss of active material model by mechanical effects
            stress_critical = self.phase_param.stress_critical
            m_LAM = self.phase_param.m_LAM
            eps_solid_stress_init = self.phase_param.epsilon_s
            if self.x_average is True:
                stress_t_surf = variables[
                    f"X-averaged {domain} {phase_name}particle surface tangential stress [Pa]"
                ]
                stress_r_surf = variables[
                    f"X-averaged {domain} {phase_name}particle surface radial stress [Pa]"
                ]
                T = variables[f"X-averaged {domain} electrode temperature [K]"]
                eps_solid_stress = variables[
                    f"X-averaged {domain} electrode {phase_name}active material volume fraction"
                ]
                eps_solid_stress_init = pybamm.x_average(eps_solid_stress_init)
            else:
                stress_t_surf = variables[
                    f"{Domain} {phase_name}particle surface tangential stress [Pa]"
                ]
                stress_r_surf = variables[
                    f"{Domain} {phase_name}particle surface radial stress [Pa]"
                ]
                T = variables[f"{Domain} electrode temperature [K]"]
                eps_solid_stress = variables[
                    f"{Domain} electrode {phase_name}active material volume fraction"
                ]
            # compute hydrostatic stress
            stress_h_surf = (stress_r_surf + 2 * stress_t_surf) / 3
            # Self-limiting: applied to the STRESS itself, before the
            # **m_LAM power law, not to the term's final output -- as
            # eps_solid -> 0 the particle-mechanics stress fields (radius/
            # concentration-gradient terms) can hit a genuine numerical
            # singularity (observed directly this session: stress_h_surf
            # reaching ~1e35 Pa, many orders of magnitude past anything
            # physical, once Si's active material is nearly consumed).
            # Damping the OUTPUT of an already-diverged power-law result
            # barely helps (a huge number times a modest factor is still
            # huge); damping the stress INPUT itself prevents the power
            # law from ever amplifying a diverging value in the first
            # place, and decays this term toward 0 as the phase is
            # exhausted instead -- the same physical argument as
            # remaining_frac elsewhere in this file (no active material
            # left behind should mean no more stress-driven loss left to
            # apply).
            # Floored at 0, not just capped at 1: without the floor, a
            # solver step that overshoots eps_solid_stress past zero flips
            # this ratio negative, which flips the WHOLE term's sign --
            # suddenly ADDING active material back instead of removing it
            # (the same latent bug found and fixed for remaining_frac
            # below, during the v7c isolation-reformulation testing).
            remaining_frac_stress = pybamm.minimum(
                pybamm.maximum(eps_solid_stress / eps_solid_stress_init, 0), 1
            )
            stress_h_surf = stress_h_surf * remaining_frac_stress
            # separate compressive and tensile stresses
            stress_h_surf_compressive = stress_h_surf * (stress_h_surf < 0)
            stress_h_surf_tensile = stress_h_surf * (stress_h_surf > 0)

            if "asymmetric stress" in lam_option:
                pybamm.citations.register("Pannala2024")
                # semi-empirical model for stress-driven LAM that includes both
                # compressive and tensile stresses
                beta_LAM_compressive = self.phase_param.beta_LAM(
                    T, direction="compressive"
                )
                beta_LAM_tensile = self.phase_param.beta_LAM(T, direction="tensile")
                j_stress_LAM = (
                    -beta_LAM_compressive
                    * (abs(stress_h_surf_compressive) / stress_critical) ** m_LAM
                    - beta_LAM_tensile
                    * (abs(stress_h_surf_tensile) / stress_critical) ** m_LAM
                )
            else:
                beta_LAM = self.phase_param.beta_LAM(T)
                # assuming that only tensile stress contributes and that the minimum
                # (tensile) hydrostatic stress is zero for full cycles
                stress_h_surf_min = stress_h_surf * 0
                j_stress_LAM = (
                    -beta_LAM
                    * ((stress_h_surf_tensile - stress_h_surf_min) / stress_critical)
                    ** m_LAM
                )

            # Diagnostic breakdown (this session's LAM-mechanism-attribution
            # investigation): expose this term's own contribution so it can
            # be read directly off the solution, rather than reconstructed
            # by hand from raw current-density/stress variables (error-
            # prone) or inferred indirectly via ablation runs (expensive).
            if self.x_average is True:
                variables.update(
                    {
                        f"X-averaged {domain} electrode {phase_name}"
                        "LAM rate from stress [s-1]": j_stress_LAM
                    }
                )
            else:
                # DFN: this submodel operates on the full (non-averaged)
                # field, so the "X-averaged ..." companion isn't created
                # automatically -- compute it explicitly so diagnostics can
                # read a single representative scalar per time, consistent
                # with every other "X-averaged ..." variable used this
                # session.
                variables.update(
                    {
                        f"{Domain} electrode {phase_name}"
                        "LAM rate from stress [s-1]": j_stress_LAM,
                        f"X-averaged {domain} electrode {phase_name}"
                        "LAM rate from stress [s-1]": pybamm.x_average(j_stress_LAM),
                    }
                )
            deps_solid_dt += j_stress_LAM

        if "reaction" in lam_option:
            beta_LAM_sei = self.phase_param.beta_LAM_sei
            if self.x_average is True:
                a_j_sei = variables[
                    f"X-averaged {domain} electrode {phase_name}SEI "
                    "volumetric interfacial current density [A.m-3]"
                ]
            else:
                a_j_sei = variables[
                    f"{Domain} electrode {phase_name}SEI volumetric "
                    "interfacial current density [A.m-3]"
                ]

            j_stress_reaction = beta_LAM_sei * a_j_sei / self.param.F
            # Diagnostic breakdown: the beta_LAM_sei-only contribution,
            # before any porosity-isolation addition below.
            if self.x_average is True:
                variables.update(
                    {
                        f"X-averaged {domain} electrode {phase_name}"
                        "LAM rate from reaction (beta_LAM_sei) [s-1]": j_stress_reaction
                    }
                )
            else:
                variables.update(
                    {
                        f"{Domain} electrode {phase_name}"
                        "LAM rate from reaction (beta_LAM_sei) [s-1]": j_stress_reaction,
                        f"X-averaged {domain} electrode {phase_name}"
                        "LAM rate from reaction (beta_LAM_sei) [s-1]": pybamm.x_average(
                            j_stress_reaction
                        ),
                    }
                )

            deps_solid_dt += j_stress_reaction

        if "porosity" in lam_option:
            # Porosity-gated LAM ("... and porosity isolation"): an
            # independent term (NOT nested under "reaction" -- a phase can
            # have "porosity isolation" without "reaction", e.g. graphite's
            # "stress-driven and porosity isolation" used for item 24's
            # redirect-gating, which deliberately does NOT want
            # beta_LAM_sei's reaction-driven term) that ramps from ~0 while
            # the electrode still has porosity headroom (relative to its
            # as-set-up value) to its full rate as porosity approaches its
            # floor. Physically: electrolyte-transport limitation as pore
            # volume collapses traps/isolates active material, so this
            # stays small pre-knee (deliberately, so it doesn't perturb the
            # already-tuned pre-knee fit) and grows sharply once the
            # electrode nears its porosity floor (post-knee).
            if self.x_average is True:
                a_j_sei = variables[
                    f"X-averaged {domain} electrode {phase_name}SEI "
                    "volumetric interfacial current density [A.m-3]"
                ]
            else:
                a_j_sei = variables[
                    f"{Domain} electrode {phase_name}SEI volumetric "
                    "interfacial current density [A.m-3]"
                ]
            beta_LAM_iso = self.phase_param.beta_LAM_iso
            eps_floor = self.domain_param.epsilon_min
            eps_solid_init = self.phase_param.epsilon_s
            if self.x_average is True:
                eps_n0 = pybamm.x_average(self.domain_param.epsilon_init)
                eps_n = variables[f"X-averaged {domain} electrode porosity"]
                eps_solid_for_iso = variables[
                    f"X-averaged {domain} electrode {phase_name}"
                    "active material volume fraction"
                ]
                eps_solid_init = pybamm.x_average(eps_solid_init)
            else:
                eps_n0 = self.domain_param.epsilon_init
                eps_n = variables[f"{Domain} electrode porosity"]
                eps_solid_for_iso = variables[
                    f"{Domain} electrode {phase_name}active material volume fraction"
                ]
            # headroom_frac: 1 at the as-set-up porosity, 0 at the floor
            headroom_frac = (eps_n - eps_floor) / (eps_n0 - eps_floor)
            headroom_frac = pybamm.minimum(pybamm.maximum(headroom_frac, 0), 1)
            eta_LAM_iso = self.phase_param.eta_LAM_iso
            # Algebraic TARGET the lagged gate state relaxes toward (set_rhs
            # uses this), not the gate actually used in the LAM term below
            # -- see get_fundamental_variables' isolation gate STATE doc-
            # comment for why this is now a relaxation ODE rather than an
            # instantaneous snapshot of porosity.
            isolation_gate_target = (1 - headroom_frac) ** eta_LAM_iso
            if self.x_average is True:
                isolation_gate = variables[
                    f"X-averaged {domain} electrode {phase_name}"
                    "porosity-isolation gate state"
                ]
                variables.update(
                    {
                        f"X-averaged {domain} electrode {phase_name}"
                        "porosity-isolation gate target": pybamm.x_average(
                            isolation_gate_target
                        )
                    }
                )
            else:
                isolation_gate = variables[
                    f"{Domain} electrode {phase_name}porosity-isolation gate state"
                ]
                variables.update(
                    {
                        f"{Domain} electrode {phase_name}"
                        "porosity-isolation gate target": isolation_gate_target,
                        f"X-averaged {domain} electrode {phase_name}"
                        "porosity-isolation gate target": pybamm.x_average(
                            isolation_gate_target
                        ),
                    }
                )
            # Self-limiting: scale by the REMAINING active-material
            # fraction (relative to its own as-set-up value) so this term
            # decays toward 0 as the phase is consumed, instead of driving
            # eps_solid through/past zero -- the same "current-focusing
            # singularity" (a = 3*eps_solid/R -> 0 forces j = i/a ->
            # infinity) documented for an earlier, unrelated isolation-
            # mechanism attempt this session. Floored at 0, not just capped
            # at 1: without the floor, a solver step that overshoots
            # eps_solid_for_iso past zero flips this ratio negative, which
            # flips the WHOLE isolation term's sign -- suddenly ADDING
            # active material back instead of removing it. a_j_sei's
            # natural decay (v7b) kept the term's magnitude small enough
            # near end-of-life that this was never triggered in practice;
            # L_sei's stronger, non-decaying drive (v7c) produces a large
            # enough instantaneous rate to overshoot in a single solver
            # step, exposing this latent bug (observed directly: LAM_Si
            # going negative, i.e. eps_solid increasing past its own
            # initial value, in the v7c reformulation's first test).
            remaining_frac = pybamm.minimum(
                pybamm.maximum(eps_solid_for_iso / eps_solid_init, 0), 1
            )
            # Driving variable: SEI THICKNESS (a stock), not a_j_sei (a
            # flow) or the main intercalation current -- v7c reformulation,
            # see TUNING STATUS item 20. This model's own porosity submodel
            # (reaction_driven_porosity.py) already makes porosity decline
            # as a multiple of SEI/plating THICKNESS growth, so L_sei is
            # the same mechanistic quantity that causes the pore collapse
            # this gate responds to -- staying internally consistent with
            # the same "SEI deposition closes the pore and traps active
            # material" story as before, just measuring the ACCUMULATED
            # deposit rather than its instantaneous deposition rate. This
            # matters because a_j_sei (tried first) is diffusion-limited
            # and crashes to ~0 once electrolyte access is cut off --
            # exactly the regime this gate is meant to activate in --
            # forcing a race between the (now tau-lagged) gate catching up
            # and a driving current that's simultaneously dying, capping
            # how far widening the gate alone (via tau) can ever reach
            # (confirmed empirically: pure-tau widening plateaus around
            # end_efc~156, well short of the real ~197 anchor). L_sei,
            # being a monotonically non-decreasing accumulated thickness,
            # does not share this problem: once new deposition stops,
            # thickness plateaus rather than collapsing, so it stays
            # available as a persistent, non-vanishing signal for however
            # long the (tau-paced) gate needs to catch up -- decoupling
            # isolation's available DURATION from a_j_sei's own finite
            # lifetime.
            if self.options["SEI reaction redirect to LAM"] == "true":
                # Redirect-to-LAM (item 23): the current-conserving
                # alternative to v7c's L_sei-based term above. Reuses
                # a_j_sei directly, matching sei_growth.py's matching gate
                # on dcdt_sei (which reduces SEI film growth -- hence new
                # LLI -- by the SAME isolation_gate fraction). Unlike v7c's
                # L_sei-based term, this does NOT add a second, independent
                # lithium-consumption channel on top of unchanged SEI
                # growth -- it is the SAME reaction current that
                # sei_growth.py has stopped counting toward SEI/LLI, now
                # counted here toward LAM instead. a_j_sei's own decay
                # (diffusion-limited by L_sei's growth) is itself slowed by
                # the redirect -- less current still going to SEI growth
                # means L_sei grows slower, keeping diffusion resistance
                # lower for longer, keeping a_j_sei alive longer than in
                # the reaction-only (v7b) case.
                redirect_LAM_yield = self.phase_param.redirect_LAM_yield
                j_iso = (
                    redirect_LAM_yield
                    * isolation_gate
                    * a_j_sei
                    / self.param.F
                    * remaining_frac
                )
            else:
                if self.x_average is True:
                    L_sei = variables[
                        f"X-averaged {domain} {phase_name}SEI thickness [m]"
                    ]
                else:
                    L_sei = variables[
                        f"{Domain} {phase_name}SEI thickness [m]"
                    ]
                # Leading minus sign: a_j_sei (v7b's driver) was always
                # negative by convention, so it supplied the "this DEPLETES
                # active material" sign for free; L_sei (a thickness) is
                # always >=0, so that sign has to be made explicit here --
                # missing this was a real, confirmed bug in the first v7c
                # implementation (directly observed: eps_solid_Si GROWING
                # past its own BoL value, i.e. the term was silently adding
                # material back instead of removing it).
                j_iso = -beta_LAM_iso * L_sei * isolation_gate * remaining_frac
            # Diagnostic breakdown: the porosity-isolation term's OWN
            # increment, separate from beta_LAM_sei's contribution above.
            if self.x_average is True:
                variables.update(
                    {
                        f"X-averaged {domain} electrode {phase_name}"
                        "LAM rate from porosity isolation [s-1]": j_iso
                    }
                )
            else:
                variables.update(
                    {
                        f"{Domain} electrode {phase_name}"
                        "LAM rate from porosity isolation [s-1]": j_iso,
                        f"X-averaged {domain} electrode {phase_name}"
                        "LAM rate from porosity isolation [s-1]": pybamm.x_average(j_iso),
                    }
                )
            deps_solid_dt += j_iso

        if "current" in lam_option:
            # obtain the rate of loss of active materials (LAM) driven by current
            if self.x_average is True:
                T = variables[f"X-averaged {domain} electrode temperature [K]"]
            else:
                T = variables[f"{Domain} electrode temperature [K]"]

            j_current_LAM = self.domain_param.LAM_rate_current(
                self.param.current_density_with_time, T
            )
            deps_solid_dt += j_current_LAM

        variables.update(
            self._get_standard_active_material_change_variables(deps_solid_dt)
        )
        return variables

    def set_rhs(self, variables):
        domain, Domain = self.domain_Domain
        phase_name = self.phase_name

        if self.x_average is True:
            eps_solid = variables[
                f"X-averaged {domain} electrode {phase_name}active material volume fraction"
            ]
            deps_solid_dt = variables[
                f"X-averaged {domain} electrode {phase_name}active material "
                "volume fraction change [s-1]"
            ]
        else:
            eps_solid = variables[
                f"{Domain} electrode {phase_name}active material volume fraction"
            ]
            deps_solid_dt = variables[
                f"{Domain} electrode {phase_name}active material volume fraction change [s-1]"
            ]

        # Loss of lithium due to loss of active material
        # See eq 37 in "Sulzer, Valentin, et al. "Accelerated battery lifetime
        # simulations using adaptive inter-cycle extrapolation algorithm."
        # Journal of The Electrochemical Society 168.12 (2021): 120531.
        lli_due_to_lam = variables[
            f"Loss of lithium due to loss of {phase_name}active material "
            f"in {domain} electrode [mol]"
        ]
        # Multiply by mol.m-3 * m3 to get mol
        c_s_rav = variables[
            f"R-averaged {domain} {phase_name}particle concentration [mol.m-3]"
        ]
        V = self.domain_param.L * self.param.A_cc

        self.rhs = {
            # minus sign because eps_solid is decreasing and LLI measures positive
            lli_due_to_lam: -V * pybamm.x_average(c_s_rav * deps_solid_dt),
            eps_solid: deps_solid_dt,
        }

        lam_option = getattr(getattr(self.options, domain), self.phase)[
            "loss of active material"
        ]
        if "porosity" in lam_option:
            tau_LAM_iso = self.phase_param.tau_LAM_iso
            if self.x_average is True:
                isolation_state = variables[
                    f"X-averaged {domain} electrode {phase_name}"
                    "porosity-isolation gate state"
                ]
                isolation_gate_target = variables[
                    f"X-averaged {domain} electrode {phase_name}"
                    "porosity-isolation gate target"
                ]
            else:
                isolation_state = variables[
                    f"{Domain} electrode {phase_name}porosity-isolation gate state"
                ]
                isolation_gate_target = variables[
                    f"{Domain} electrode {phase_name}porosity-isolation gate target"
                ]
            # First-order relaxation toward the algebraic target -- see
            # get_fundamental_variables' doc-comment for why.
            self.rhs[isolation_state] = (
                isolation_gate_target - isolation_state
            ) / tau_LAM_iso

    def set_initial_conditions(self, variables):
        domain, Domain = self.domain_Domain
        phase_name = self.phase_name

        eps_solid_init = self.phase_param.epsilon_s

        if self.x_average is True:
            eps_solid_xav = variables[
                f"X-averaged {domain} electrode {phase_name}active material volume fraction"
            ]
            self.initial_conditions = {eps_solid_xav: pybamm.x_average(eps_solid_init)}
        else:
            eps_solid = variables[
                f"{Domain} electrode {phase_name}active material volume fraction"
            ]
            self.initial_conditions = {eps_solid: eps_solid_init}

        lli_due_to_lam = variables[
            f"Loss of lithium due to loss of {phase_name}active material "
            f"in {domain} electrode [mol]"
        ]
        self.initial_conditions[lli_due_to_lam] = pybamm.Scalar(0)

        lam_option = getattr(getattr(self.options, domain), self.phase)[
            "loss of active material"
        ]
        if "porosity" in lam_option:
            # Starts at 0, matching the algebraic target's own BoL value
            # (full porosity headroom -> isolation_gate_target = 0).
            if self.x_average is True:
                isolation_state = variables[
                    f"X-averaged {domain} electrode {phase_name}"
                    "porosity-isolation gate state"
                ]
            else:
                isolation_state = variables[
                    f"{Domain} electrode {phase_name}porosity-isolation gate state"
                ]
            self.initial_conditions[isolation_state] = pybamm.Scalar(0)
