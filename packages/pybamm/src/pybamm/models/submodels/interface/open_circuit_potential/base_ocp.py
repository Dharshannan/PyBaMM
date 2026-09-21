#
# Base class for open-circuit potential
#
import pybamm
from pybamm.models.submodels.interface.base_interface import BaseInterface


class BaseOpenCircuitPotential(BaseInterface):
    """
    Base class for open-circuit potentials

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
        super().__init__(param, domain, reaction, options=options, phase=phase)
        self.x_average = x_average

    def _alias_ocp_as_equilibrium(self, variables):
        """Publish the "equilibrium open-circuit potential [V]" variables as
        aliases of the main OCP variables. For OCP models without hysteresis
        (single OCP, MSMR) the equilibrium OCP equals the OCP, and the shapes
        are already handled by :meth:`_get_standard_ocp_variables`.
        """
        domain, Domain = self.domain_Domain
        domain_options = getattr(self.options, domain)
        phase_name = self.phase_name

        variables.update(
            {
                f"{Domain} electrode {phase_name}equilibrium open-circuit potential [V]": variables[
                    f"{Domain} electrode {phase_name}open-circuit potential [V]"
                ],
                f"X-averaged {domain} electrode {phase_name}equilibrium open-circuit potential [V]": variables[
                    f"X-averaged {domain} electrode {phase_name}open-circuit potential [V]"
                ],
            }
        )
        if domain_options["particle size"] == "distribution":
            variables.update(
                {
                    f"{Domain} electrode {phase_name}equilibrium open-circuit potential distribution [V]": variables[
                        f"{Domain} electrode {phase_name}open-circuit potential distribution [V]"
                    ],
                    f"X-averaged {domain} electrode {phase_name}equilibrium open-circuit potential distribution [V]": variables[
                        f"X-averaged {domain} electrode {phase_name}open-circuit potential distribution [V]"
                    ],
                }
            )

    def _apply_ocp_aging_deformation(self, variables, ocp_surf, ocp_bulk):
        """Item 25: aging-dependent OCP deformation (U_new = scale*U_base +
        shift), ramped by this phase's own LAM FRACTION (1 - eps_solid /
        eps_solid_BOL, clipped to [0, 1]) -- see lithium_ion_parameters.py's
        ocp_aging_deform_scale/shift doc-comment for the CELL064/manuscript
        motivation. NOT tied to the porosity-isolation gate STATE: that gate
        is deliberately slow-relaxing (tau_LAM_iso can be ~1e8 s) so that a
        small, SUSTAINED gate value integrates into large cumulative LAM
        growth over the whole simulated life -- but the OCP deformation is
        an INSTANTANEOUS multiplier, not a cumulative one, so tying it to the
        same near-zero-most-of-the-time gate value made it numerically
        negligible even at RPT5 despite LAM_Si correctly reaching ~80% there
        (confirmed empirically: switching from the gate to the LAM fraction
        below was needed after the gate-based version moved voltage RMSE by
        <0.001 V). LAM fraction is exactly the "how far along is this phase's
        degradation" progress variable already validated against real
        LAM_Si data, so it's the natural weight for a mechanism ("burn-out")
        that time-integrated LAM growth also represents.
        Only active when this phase's LAM option includes "porosity" (scopes
        this to phases the redirect/isolation mechanism already applies to,
        so ocp_aging_deform_scale/shift only need defining for those) and the
        global option is "true"; returns (ocp_surf, ocp_bulk) unchanged
        otherwise, so fully backward-compatible. Shared by every OCP submodel
        variant (single, hysteresis, ...) rather than duplicated in each.
        """
        if self.reaction != "lithium-ion main":
            return ocp_surf, ocp_bulk
        domain, Domain = self.domain_Domain
        lam_option = getattr(getattr(self.options, domain), self.phase)[
            "loss of active material"
        ]
        if (
            self.options["open-circuit potential aging deformation"] != "true"
            or "porosity" not in lam_option
        ):
            return ocp_surf, ocp_bulk

        phase_name = self.phase_name
        end_scale = self.phase_param.ocp_aging_deform_scale
        end_shift = self.phase_param.ocp_aging_deform_shift
        eps_solid_init = self.phase_param.epsilon_s

        eps_solid_surf = variables[
            f"{Domain} electrode {phase_name}active material volume fraction"
        ]
        lam_frac_surf = pybamm.minimum(
            pybamm.maximum(1 - eps_solid_surf / eps_solid_init, 0), 1
        )
        ocp_surf = ocp_surf + lam_frac_surf * ((end_scale - 1) * ocp_surf + end_shift)

        eps_solid_bulk = variables[
            f"X-averaged {domain} electrode {phase_name}active material volume fraction"
        ]
        lam_frac_bulk = pybamm.minimum(
            pybamm.maximum(1 - eps_solid_bulk / eps_solid_init, 0), 1
        )
        ocp_bulk = ocp_bulk + lam_frac_bulk * ((end_scale - 1) * ocp_bulk + end_shift)

        return ocp_surf, ocp_bulk

    def _get_standard_ocp_variables(self, ocp_surf, ocp_bulk, dUdT):
        domain, Domain = self.domain_Domain
        reaction_name = self.reaction_name

        # Update size variables then size average.
        if ocp_surf.domain in [
            ["negative particle size"],
            ["positive particle size"],
            ["negative primary particle size"],
            ["positive primary particle size"],
            ["negative secondary particle size"],
            ["positive secondary particle size"],
        ]:
            variables = self._get_standard_size_distribution_ocp_variables(
                ocp_surf, dUdT
            )
            ocp_surf = pybamm.size_average(ocp_surf)
            dUdT = pybamm.size_average(dUdT)
        else:
            variables = {}

        # Average, and broadcast if necessary
        dUdT_av = pybamm.x_average(dUdT)
        ocp_surf_av = pybamm.x_average(ocp_surf)
        if self.options.electrode_types[domain] == "planar":
            # Half-cell domain, ocp_surf should not be broadcast
            pass
        elif ocp_surf.domain == []:
            ocp_surf = pybamm.FullBroadcast(
                ocp_surf, f"{domain} electrode", "current collector"
            )
        elif ocp_surf.domain == ["current collector"]:
            ocp_surf = pybamm.PrimaryBroadcast(ocp_surf, f"{domain} electrode")

        # Particle overpotential is the difference between the average(U(c_surf)) and
        # U(c_bulk), i.e. the overpotential due to concentration gradients in the
        # particle
        eta_particle = ocp_surf_av - ocp_bulk
        variables.update(
            {
                f"{Domain} electrode {reaction_name}"
                "open-circuit potential [V]": ocp_surf,
                f"X-averaged {domain} electrode {reaction_name}"
                "open-circuit potential [V]": ocp_surf_av,
                f"{Domain} electrode {reaction_name}"
                "bulk open-circuit potential [V]": ocp_bulk,
                f"{Domain} {reaction_name}particle concentration "
                "overpotential [V]": eta_particle,
            }
        )
        if self.reaction in ["lithium-ion main", "lead-acid main"]:
            variables.update(
                {
                    f"{Domain} electrode {reaction_name}entropic change [V.K-1]": dUdT,
                    f"X-averaged {domain} electrode {reaction_name}entropic change [V.K-1]": dUdT_av,
                }
            )

        return variables

    def _get_standard_size_distribution_ocp_variables(self, ocp, dUdT):
        domain, Domain = self.domain_Domain
        reaction_name = self.reaction_name

        # X-average or broadcast to electrode if necessary
        if ocp.domains["secondary"] != [f"{domain} electrode"]:
            ocp_av = ocp
            ocp = pybamm.SecondaryBroadcast(ocp, f"{domain} electrode")
        else:
            ocp_av = pybamm.x_average(ocp)

        if dUdT.domains["secondary"] != [f"{domain} electrode"]:
            dUdT_av = dUdT
            dUdT = pybamm.SecondaryBroadcast(dUdT, f"{domain} electrode")
        else:
            dUdT_av = pybamm.x_average(dUdT)

        variables = {
            f"{Domain} electrode {reaction_name}"
            "open-circuit potential distribution [V]": ocp,
            f"X-averaged {domain} electrode {reaction_name}"
            "open-circuit potential distribution [V]": ocp_av,
        }
        if self.reaction_name == "":
            variables.update(
                {
                    f"{Domain} electrode entropic change "
                    "(size-dependent) [V.K-1]": dUdT,
                    f"X-averaged {domain} electrode entropic change "
                    "(size-dependent) [V.K-1]": dUdT_av,
                }
            )

        return variables

    def _get_stoichiometry_and_temperature(self, variables):
        domain, Domain = self.domain_Domain
        domain_options = getattr(self.options, domain)
        phase_name = self.phase_name

        sto_bulk = variables[f"{Domain} electrode {phase_name}stoichiometry"]
        T = variables[f"{Domain} electrode temperature [K]"]
        T_bulk = pybamm.xyzs_average(T)

        # For "particle-size distribution" models, take distribution version
        # of sto_surf that depends on particle size.
        if domain_options["particle size"] == "distribution":
            sto_surf = variables[
                f"{Domain} {phase_name}particle surface stoichiometry distribution"
            ]
            # If variable was broadcast, take only the orphan
            if isinstance(sto_surf, pybamm.Broadcast) and isinstance(
                T, pybamm.Broadcast
            ):
                sto_surf = sto_surf.orphans[0]
                T = T.orphans[0]
            T = pybamm.PrimaryBroadcast(T, [f"{domain} {phase_name}particle size"])
        else:
            sto_surf = variables[f"{Domain} {phase_name}particle surface stoichiometry"]
            # If variable was broadcast, take only the orphan
            if isinstance(sto_surf, pybamm.Broadcast) and isinstance(
                T, pybamm.Broadcast
            ):
                sto_surf = sto_surf.orphans[0]
                T = T.orphans[0]

        return sto_surf, sto_bulk, T, T_bulk
