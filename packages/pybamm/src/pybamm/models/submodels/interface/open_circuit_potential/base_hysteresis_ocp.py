#
# Base class for open-circuit potential with hysteresis
#
import pybamm

from . import BaseOpenCircuitPotential


class BaseHysteresisOpenCircuitPotential(BaseOpenCircuitPotential):
    """
    Base class for open-circuit potentials with single-state hysteresis

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    domain : str
        The domain to implement the model, either: 'Negative' or 'Positive'.
    reaction : str
        The name of the reaction being implemented
    options: dict
        A dictionary of options to be passed to the model. See
        :class:`pybamm.BaseBatteryModel`
    phase : str, optional
        Phase of the particle (default is "primary")
    x_average : bool
        Whether the particle concentration is averaged over the x-direction. Default is False.
    """

    def __init__(
        self, param, domain, reaction, options, phase="primary", x_average=False
    ):
        super().__init__(
            param, domain, reaction, options=options, phase=phase, x_average=x_average
        )

    def _get_hysteresis_state_variables(self):
        domain, Domain = self.domain_Domain
        domain_options = getattr(self.options, domain)
        phase_name = self.phase_name

        if self.x_average is False:
            h = pybamm.Variable(
                f"{Domain} electrode {phase_name}hysteresis state",
                domains={
                    "primary": f"{domain} electrode",
                    "secondary": "current collector",
                },
            )
            h_x_av = pybamm.x_average(h)
        else:
            h_x_av = pybamm.Variable(
                f"X-averaged {domain} electrode {phase_name}hysteresis state",
                domains={
                    "primary": "current collector",
                },
            )
            h = pybamm.PrimaryBroadcast(h_x_av, [f"{domain} electrode"])
        variables = {
            f"{Domain} electrode {phase_name}hysteresis state": h,
            f"X-averaged {domain} electrode {phase_name}hysteresis state": h_x_av,
        }
        if domain_options["particle size"] == "distribution":
            h_dist = pybamm.PrimaryBroadcast(h, [f"{domain} {phase_name}particle size"])
            h_dist_x_av = pybamm.x_average(h_dist)
            variables.update(
                {
                    f"{Domain} electrode {phase_name}hysteresis state distribution": h_dist,
                    f"X-averaged {domain} electrode {phase_name}hysteresis state distribution": h_dist_x_av,
                }
            )
        return variables

    # ------------------------------------------------------------------
    # "silicon crystallisation": crystalline fraction phi (c-Li15Si4)
    # ------------------------------------------------------------------
    def _crystallisation_on(self):
        domain, _ = self.domain_Domain
        phase_options = getattr(getattr(self.options, domain), self.phase)
        return (
            self.reaction == "lithium-ion main"
            and phase_options["silicon crystallisation"] == "true"
        )

    def _get_crystal_variables(self):
        domain, Domain = self.domain_Domain
        phase_name = self.phase_name
        phi = pybamm.Variable(
            f"{Domain} electrode {phase_name}crystalline fraction",
            domains={"primary": f"{domain} electrode", "secondary": "current collector"},
        )
        return {
            f"{Domain} electrode {phase_name}crystalline fraction": phi,
            f"X-averaged {domain} electrode {phase_name}crystalline fraction": pybamm.x_average(phi),
        }

    def _set_crystal_rhs(self, variables):
        """d phi/dt = k_c (phi + phi_seed)(1 - phi) S(sto > sto_c) S(lithiating) - k_a max(phi - f_plateau, 0),
        f_plateau = (sto - x_t)/(1 - x_t) clipped to [0, 1]: the crystal grows while the phase is
        deeply lithiated and is consumed as its plateau sites delithiate (two-phase c-Li15Si4 ->
        amorphous), as in Jiang et al. (2020)."""
        domain, Domain = self.domain_Domain
        phase_name = self.phase_name
        p = self.phase_param
        phi = variables[f"{Domain} electrode {phase_name}crystalline fraction"]
        sto_surf, _, _, _ = self._get_stoichiometry_and_temperature(variables)
        if sto_surf.domain != phi.domain:  # e.g. x-averaged particles: broadcast onto the electrode
            sto_surf = pybamm.PrimaryBroadcast(sto_surf, [f"{domain} electrode"])
        deep = pybamm.sigmoid(p.cryst_onset_sto, sto_surf, 400)
        x_t = p.cryst_plateau_sto
        f_plateau = pybamm.maximum(pybamm.minimum((sto_surf - x_t) / (1 - x_t), 1), 0)
        # only while the phase is lithiating (charge and CV hold), as Jiang's nucleation step f3: not
        # during the brief Li hand-over into Si at the start of a discharge under the current-sigmoid switch
        current = variables["Total current density [A.m-2]"]
        lithiation_current = -current if domain == "negative" else current
        lithiating = pybamm.sigmoid(0, lithiation_current, 100)
        grow = p.cryst_rate * (phi + p.cryst_seed) * (1 - phi) * deep * lithiating
        consume = p.cryst_consumption_rate * pybamm.maximum(phi - f_plateau, 0)
        self.rhs[phi] = grow - consume

    def _set_crystal_initial_conditions(self, variables):
        domain, Domain = self.domain_Domain
        phi = variables[f"{Domain} electrode {self.phase_name}crystalline fraction"]
        self.initial_conditions[phi] = pybamm.Scalar(0)

    def _get_coupled_variables(self, variables):
        domain, Domain = self.domain_Domain
        domain_options = getattr(self.options, domain)
        phase_name = self.phase_name

        if self.reaction == "lithium-ion main":
            # Get equilibrium, delithiation and lithiation OCPs
            sto_surf, sto_bulk, T, T_bulk = self._get_stoichiometry_and_temperature(
                variables
            )
            U_eq = self.phase_param.U(sto_surf, T)
            U_lith = self.phase_param.U(sto_surf, T, "lithiation")
            U_lith_bulk = self.phase_param.U(sto_bulk, T_bulk, "lithiation")
            U_delith = self.phase_param.U(sto_surf, T, "delithiation")
            U_delith_bulk = self.phase_param.U(sto_bulk, T_bulk, "delithiation")

            # Aging deformation (item 25 / throughput-driver follow-up):
            # applied to the DELITHIATION branch ONLY, before the
            # lithiation/delithiation blend below -- not to the blended
            # ocp_surf/ocp_bulk (as an earlier version of this did). The
            # real per-RPT calibration data behind BOTH drivers (LAM
            # fraction and throughput) comes entirely from RPT DISCHARGE
            # (delithiation) eSOH fits; deforming lithiation too would
            # silently extend an unmeasured correction onto the CHARGE
            # curve as well. Since charge steps are voltage-limited
            # (CC-CV to a fixed cutoff), a lowered lithiation OCP would
            # make the cell reach that cutoff EARLIER (at a lower degree
            # of lithiation) as aging progresses -- throttling how much
            # lithium the model allows into the phase each aging cycle, a
            # spurious, compounding side effect that could itself look
            # like (or interact with) excess overpotential rather than
            # being a validated thermodynamic correction. Doing the
            # deformation here (pre-blend) instead of post-blend also
            # means the published "OCP hysteresis [V]" (H, below) reflects
            # the now-asymmetric gap correctly.
            if self._crystallisation_on():
                # delithiation branch blended with the crystalline fraction (before the aging
                # deformation, which then acts on the blended curve)
                phi = variables[f"{Domain} electrode {phase_name}crystalline fraction"]
                phi_av = variables[f"X-averaged {domain} electrode {phase_name}crystalline fraction"]
                phi_surf = phi
                if sto_surf.domain != phi.domain:
                    phi_surf = phi_av
                U_delith = self.phase_param.U_delith_cryst(sto_surf, phi_surf, T)
                U_delith_bulk = self.phase_param.U_delith_cryst(sto_bulk, phi_av, T_bulk)
            U_delith, U_delith_bulk = self._apply_ocp_aging_deformation(
                variables, U_delith, U_delith_bulk
            )
            # Optional: the lithiation branch too, with its own end
            # scale/shift from the charge eSOH fits (same LAM-fraction ramp).
            if self.options["OCP aging deformation branches"] == "both":
                U_lith, U_lith_bulk = self._apply_ocp_aging_deformation(
                    variables, U_lith, U_lith_bulk, branch="lithiation"
                )

            H = U_lith - U_delith

            # Under particle-size distribution the OCPs above live on the
            # "<domain> particle size" domain. Publish the distribution-named
            # variables at that shape, and size-average before publishing under
            # the electrode-level names so they line up with other
            # electrode-domain variables (e.g. in the thermal submodel).
            if domain_options["particle size"] == "distribution":
                variables.update(
                    {
                        f"{Domain} electrode {phase_name}OCP hysteresis distribution [V]": H,
                        f"X-averaged {domain} electrode {phase_name}OCP hysteresis distribution [V]": pybamm.x_average(
                            H
                        ),
                        f"{Domain} electrode {phase_name}equilibrium open-circuit potential distribution [V]": U_eq,
                        f"X-averaged {domain} electrode {phase_name}equilibrium open-circuit potential distribution [V]": pybamm.x_average(
                            U_eq
                        ),
                    }
                )
                U_eq = pybamm.size_average(U_eq)
                H = pybamm.size_average(H)

            variables.update(
                {
                    f"{Domain} electrode {phase_name}OCP hysteresis [V]": H,
                    f"X-averaged {domain} electrode {phase_name}OCP hysteresis [V]": pybamm.x_average(
                        H
                    ),
                    f"{Domain} electrode {phase_name}equilibrium open-circuit potential [V]": U_eq,
                    f"X-averaged {domain} electrode {phase_name}equilibrium open-circuit potential [V]": pybamm.x_average(
                        U_eq
                    ),
                }
            )

            # Get correct hysteresis state variable
            if domain_options["particle size"] == "distribution":
                h = variables[
                    f"{Domain} electrode {phase_name}hysteresis state distribution"
                ]
                h_x_av = variables[
                    f"X-averaged {domain} electrode {phase_name}hysteresis state distribution"
                ]
            else:
                h = variables[f"{Domain} electrode {phase_name}hysteresis state"]
                h_x_av = variables[
                    f"X-averaged {domain} electrode {phase_name}hysteresis state"
                ]

            # check MPM
            if (
                domain_options["particle size"] == "distribution"
                and "current collector" in sto_surf.domains["secondary"]
            ):
                ocp_surf = (1 + h_x_av) / 2 * U_delith + (1 - h_x_av) / 2 * U_lith
            else:
                ocp_surf = (1 + h) / 2 * U_delith + (1 - h) / 2 * U_lith

            # Size average
            h_s_av = pybamm.size_average(h_x_av)
            ocp_bulk = (1 + h_s_av) / 2 * U_delith_bulk + (1 - h_s_av) / 2 * U_lith_bulk
            dUdT = self.phase_param.dUdT(sto_surf)

        variables.update(self._get_standard_ocp_variables(ocp_surf, ocp_bulk, dUdT))
        return variables

    def set_initial_conditions(self, variables):
        domain, Domain = self.domain_Domain
        phase_name = self.phase_name

        h_init = self.phase_param.h_init
        if self.x_average is False:
            h = variables[f"{Domain} electrode {phase_name}hysteresis state"]
        else:
            h = variables[f"X-averaged {domain} electrode {phase_name}hysteresis state"]
            h_init = pybamm.x_average(h_init)

        self.initial_conditions[h] = h_init
