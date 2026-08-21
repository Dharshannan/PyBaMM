#
# Class for reaction driven porosity changes as a multiple of SEI/plating thicknesses
#
import pybamm

from .base_porosity import BaseModel


class ReactionDriven(BaseModel):
    """Reaction-driven porosity changes as a multiple of SEI/plating thicknesses

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    options : dict
        Options dictionary passed from the full model
    x_average : bool
        Whether to use x-averaged variables (SPM, SPMe, etc) or full variables (DFN)
    """

    def __init__(self, param, options, x_average):
        super().__init__(param, options)
        self.x_average = x_average

    def _transmitted_fraction(self, eps_struct, domain_param):
        """Dispatches on options["pore buffering transition"] between the
        semi-empirical tanh blend (default) and the physically-derived
        compliance-ratio form. See _transmitted_fraction_tanh and
        _transmitted_fraction_physical for the two implementations."""
        if self.options["pore buffering transition"] == "physical":
            return self._transmitted_fraction_physical(eps_struct, domain_param)
        return self._transmitted_fraction_tanh(eps_struct, domain_param)

    def _transmitted_fraction_tanh(self, eps_struct, domain_param):
        """Smooth (C-infinity) tanh blend between the "maximal buffering"
        plateau (f0, at high structural porosity) and "pores closed" (f=1,
        at/below closure porosity), replacing the original three-regime
        piecewise-linear f(eps_struct) (pore-buffering volume partition,
        "Making the Pore-Buffering Model Physical" notes, Eq. 18).

        The piecewise-linear version is continuous but has kinks
        (discontinuous derivative) at eps_min/eps_min+width; once those
        were tuned to a narrow band to time the buffering-saturation
        transition against the capacity knee (see CHANGES.md item 7), the
        short linear ramp between the kinks made the transition look and
        behave like a near step-function. This tanh blend removes the
        kinks entirely (smooth everywhere) and needs no min/max clamping
        (tanh is already bounded).

        Uses eps_min_transfer (closure porosity, a genuine physical
        reference point) and eps_transfer_width (NOT a second porosity
        threshold -- see CHANGES.md item 9 for why the original
        "eps_max_transfer" framing was dropped: once tuned for knee
        alignment, eps_init ends up far above it, so the cell never
        actually occupies Eq. 18's "shared" regime near BoL -- the
        parameter's only real role had become "sets the transition width",
        so it's named that directly). eps_min_transfer + eps_transfer_width
        is approximately the 10% point of the transition, eps_min_transfer
        itself approximately the 90% point (f(eps_min) approx= f0 +
        0.9*(1-f0)) -- not exact endpoints, a true sigmoid only approaches
        its asymptotes.

        Depends only on the *structural* (irreversible) porosity, so it is
        monotonic over life and introduces no algebraic loop with the
        reversible buffered volume computed downstream."""
        eps_min = domain_param.eps_min_transfer
        width = domain_param.eps_transfer_width
        f0 = domain_param.f_transmit_min
        eps_mid = eps_min + width / 2
        # 2*arctanh(0.8) ~= 2.1972: steepness such that eps_min and
        # eps_min+width sit at the 90%/10% points of the blend.
        steepness = 2.1972 / width
        blend = 0.5 * (1 - pybamm.tanh(steepness * (eps_struct - eps_mid)))
        return f0 + (1 - f0) * blend

    def _transmitted_fraction_physical(self, eps_struct, domain_param):
        """Compliance-ratio ("physical") transmitted fraction, derived from
        the pore-network/stack stiffness balance in the pore-buffering
        notes (Eq. 24): f = 1/(1 + K*C_pore(eps_struct)), where C_pore(eps)
        is the pore network's own compliance -- how much it can still
        compress under load.

        Modelled as a smooth, saturating function of the headroom above
        closure:
            C_pore(eps) = C_max * (1 - exp(-headroom(eps) / eps_transfer_width))
            headroom(eps) = softplus(eps - eps_min_transfer, 0, 100)
        -- headroom (hence C_pore) smoothly floors at 0 as
        eps -> eps_min_transfer (percolation closure: the remaining pore
        network rigidifies, transmitting everything, f -> 1), and C_pore
        saturates at C_max as eps grows large (abundant pore space: even
        then, compliance can't exceed a finite value set by the stack's own
        relative stiffness -- this is what gives the f0 plateau, unlike an
        unbounded C_pore(eps) which would drive f -> 0 instead of f0).

        K*C_max is not a free parameter: it's fixed by requiring the
        eps -> infinity limit to equal f0 exactly (f0 = k_BoL, the BoL
        calibration), giving K*C_max = (1-f0)/f0. eps_transfer_width (a
        pore-compliance decay length -- see CHANGES.md item 9 for why this
        replaced a second "eps_max_transfer" porosity threshold) is the
        only remaining shape parameter, same "width sets transition
        sharpness" role as in the tanh version, so eps_min_transfer/
        eps_transfer_width/f_transmit_min don't need separate tuning per
        transition option.

        Unlike the tanh blend this is NOT symmetric in eps_struct (a
        genuine feature of the underlying compliance-ratio physics, not an
        artifact) -- see si_gr_expansion_precursor/test_pore_buffering/
        f_form_comparison.png for a direct comparison."""
        eps_min = domain_param.eps_min_transfer
        width = domain_param.eps_transfer_width
        f0 = domain_param.f_transmit_min
        k_eps = 100.0  # softplus sharpness, matches the porosity-floor convention
        headroom = pybamm.softplus(eps_struct - eps_min, 0, k_eps)
        c_pore_normalised = 1 - pybamm.exp(-headroom / width)
        k_cmax = (1 - f0) / f0
        return 1 / (1 + k_cmax * c_pore_normalised)

    def get_coupled_variables(self, variables):
        eps_dict = {}
        for domain in self.options.whole_cell_domains:
            delta_eps_k = 0
            if domain != "separator":  # separator porosity does not change
                dom = domain.split()[0]
                Domain = dom.capitalize()
                SEI_option = getattr(self.options, dom)["SEI"]
                phases_option = getattr(self.options, dom)["particle phases"]
                phases = self.options.phases[dom]
                for phase in phases:
                    if phases_option == "1" and phase == "primary":
                        # `domain` has one phase
                        phase_name = ""
                        pref = ""
                    else:
                        # `domain` has more than one phase
                        phase_name = phase + " "
                        pref = phase.capitalize() + ": "
                    a_k = variables[
                        f"{Domain} electrode {phase_name}"
                        "surface area to volume ratio [m-1]"
                    ]
                    if SEI_option == "none":
                        L_sei_0 = pybamm.Scalar(0)
                    else:
                        L_sei_0 = pybamm.Parameter(f"{pref}Initial SEI thickness [m]")
                    L_sei_k = variables[f"{Domain} {phase_name}SEI thickness [m]"]
                    L_pl_k = variables[
                        f"{Domain} {phase_name}lithium plating thickness [m]"
                    ]
                    L_dead_k = variables[
                        f"{Domain} {phase_name}dead lithium thickness [m]"
                    ]
                    L_sei_cr_k = variables[
                        f"{Domain} {phase_name}SEI on cracks thickness [m]"
                    ]
                    roughness_k = variables[
                        f"{Domain} {phase_name}electrode roughness ratio"
                    ]

                    L_tot = (
                        (L_sei_k - L_sei_0)
                        + L_pl_k
                        + L_dead_k
                        + L_sei_cr_k * (roughness_k - 1)
                    )

                    # This assumes a thin film so curvature effects are neglected.
                    # They could be included (e.g. for a sphere it is
                    # a_n * (L_tot + L_tot ** 2 / R_n + L_tot ** # 3 / (3 * R_n ** 2)))
                    # but it is not clear if it is relevant or not.
                    delta_eps_k += -a_k * L_tot

            domain_param = self.param.domain_params[domain.split()[0]]
            # Structural (irreversible-only) porosity: SEI/plating/dead-Li/
            # crack terms only, independent of pore buffering. This is the
            # state variable that actually evolves with degradation; the
            # pore-buffering partition below reads it but never feeds back
            # into it, so there is no algebraic loop.
            eps_struct = domain_param.epsilon_init + delta_eps_k
            if domain != "separator":
                variables[f"{Domain} electrode structural porosity"] = eps_struct
            eps_k = eps_struct

            # Pore buffering (volume partition): partitions active-material
            # swelling between pore-volume buffering and electrode thickness
            # change, instead of routing 100% of it to thickness. Opt-in via
            # options["pore buffering"]; scoped to the negative electrode only
            # (the Si/Gr precursor physics this implements is anode-specific --
            # see si_gr_expansion_precursor/pore_buffering_implementation_plan.md).
            # Falls through to the unbuffered eps_struct above when off, or
            # when the electrode has no particle-mechanics thickness-change
            # output to partition (e.g. "particle mechanics": "none").
            thickness_key = f"{Domain} electrode thickness change [m]"
            if (
                domain == "negative electrode"
                and self.options["pore buffering"] == "true"
                and thickness_key in variables
            ):
                thickness_change_unbuffered = variables[thickness_key]
                dv_solid = thickness_change_unbuffered / (
                    self.param.n_electrodes_parallel * domain_param.L
                )

                # eps_struct is a full spatial (x-resolved) field here (it
                # inherits x/y/z-dependence from the "electrode porosity"
                # FunctionParameter, needed by the electrolyte transport
                # PDEs), but dv_solid is already x-averaged -- "electrode
                # thickness change" is inherently a lumped quantity. Use the
                # x-averaged structural porosity for the partition math so
                # dv_buffered/dv_thickness stay scalar (x-averaged), matching
                # dv_solid; eps_struct itself (spatial) is still used below
                # for the porosity output eps_k, so the buffered volume is
                # simply subtracted uniformly across the x-profile.
                eps_struct_avg = pybamm.x_average(eps_struct)
                f = self._transmitted_fraction(eps_struct_avg, domain_param)

                # Hard geometric ceiling: cannot buffer more void than exists
                # above closure. Binds only on net swelling (dv_solid > 0);
                # on delithiation the buffer reopens reversibly (Eq. 19).
                # headroom is floored at 0: eps_struct is the *raw*
                # (unfloored) irreversible porosity, which can go negative
                # under heavy degradation late in life (the numerical
                # epsilon_min softplus floor further down normally hides
                # this from every other consumer). Without the floor here,
                # a negative headroom would make pybamm.minimum pick the
                # negative headroom over a correctly-clamped-to-zero
                # (1 - f) * dv_solid once f saturates at 1, injecting a
                # spurious large negative dv_buffered instead of 0.
                headroom = pybamm.maximum(
                    eps_struct_avg - domain_param.eps_min_transfer, 0
                )
                dv_buffered = pybamm.minimum((1 - f) * dv_solid, headroom)
                dv_thickness = dv_solid - dv_buffered

                thickness_change_buffered = (
                    self.param.n_electrodes_parallel * dv_thickness * domain_param.L
                )
                variables[thickness_key] = thickness_change_buffered
                variables[f"{Domain} electrode buffered volume change"] = dv_buffered
                variables[f"{Domain} electrode solid volume change"] = dv_solid
                variables[f"{Domain} electrode transfer ratio k"] = (
                    pybamm.x_average(dv_thickness)
                    / (pybamm.x_average(dv_solid) + 1e-30)
                )

                # "Cell thickness change [m]" was already built by particle
                # mechanics (base_mechanics.py's _aggregate_cell_thickness_
                # change), which runs earlier in the build order, from the
                # UNBUFFERED thickness_change_unbuffered captured above --
                # pybamm expressions are immutable, so overwriting
                # variables[thickness_key] here does not retroactively
                # propagate into that already-built "Cell thickness change"
                # expression. Correct it with a delta (subtract the stale
                # unbuffered contribution, add the buffered one) rather than
                # re-deriving the full neg+pos+thermal formula, to avoid
                # duplicating/drifting from that formula.
                cell_key = "Cell thickness change [m]"
                if cell_key in variables:
                    variables[cell_key] = (
                        variables[cell_key]
                        - thickness_change_unbuffered
                        + thickness_change_buffered
                    )

                # eps_min_transfer only gates the partition above (via f and
                # headroom) -- it is NOT a floor on the reported/transport
                # porosity itself. Once buffering saturates (dv_buffered ->
                # 0), the actual porosity keeps decreasing with eps_struct,
                # bounded only by the existing numerical softplus floor
                # (domain_param.epsilon_min, e.g. 0.01) applied below -- the
                # same floor already tuned for this recipe's resistance-
                # driven knee. (Using eps_min_transfer as a second, higher
                # floor here would silently override that tuning once
                # buffering saturates, since eps_min_transfer is typically
                # set well above the numerical floor.)
                eps_k = eps_struct - dv_buffered

            # Porosity floor: softplus-smoothed lower bound, set via
            # "{Domain} electrode porosity floor" (default 0.08). Higher
            # values are more stable for the solver but suppress the
            # resistance increase (and knee-sharpening) that a genuinely
            # low porosity would cause via the Bruggeman effective-
            # transport relations; set to 0.0 to let porosity approach
            # zero as SEI/plating consume pore volume. This is a purely
            # numerical backstop, independent of the pore-buffering closure
            # floor (eps_min_transfer) above, which is a physical parameter.
            if domain != "separator":
                eps_min = domain_param.epsilon_min
                k_eps = 100.0  # softplus sharpness
                eps_k = pybamm.softplus(eps_k, eps_min, k_eps)

            eps_dict[domain] = eps_k

        variables.update(self._get_standard_porosity_variables(eps_dict))

        return variables

    def add_events_from(self, variables):
        eps_p = variables["Positive electrode porosity"]
        self.events.append(
            pybamm.Event(
                "Zero positive electrode porosity cut-off",
                pybamm.min(eps_p),
                pybamm.EventType.TERMINATION,
            )
        )
        self.events.append(
            pybamm.Event(
                "Max positive electrode porosity cut-off",
                1 - pybamm.max(eps_p),
                pybamm.EventType.TERMINATION,
            )
        )
        if "negative electrode" in self.options.whole_cell_domains:
            eps_n = variables["Negative electrode porosity"]
            self.events.append(
                pybamm.Event(
                    "Zero negative electrode porosity cut-off",
                    pybamm.min(eps_n),
                    pybamm.EventType.TERMINATION,
                )
            )
            self.events.append(
                pybamm.Event(
                    "Max negative electrode porosity cut-off",
                    1 - pybamm.max(eps_n),
                    pybamm.EventType.TERMINATION,
                )
            )
