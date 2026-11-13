#
# Class for cracking
#
import numpy as np

import pybamm

from .base_mechanics import BaseMechanics


class CrackPropagation(BaseMechanics):
    """
    Cracking behaviour in electrode particles. See :footcite:t:`Ai2019` for mechanical
    model (thickness change) and :footcite:t:`Deshpande2012` for cracking model.

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    domain : str
        The domain of the model either 'Negative' or 'Positive'
    x_average : bool
        Whether to use x-averaged variables (SPM, SPMe, etc) or full variables (DFN)
    options: dict
        A dictionary of options to be passed to the model.
        See :class:`pybamm.BaseBatteryModel`
    phase : str, optional
        Phase of the particle (default is "primary")

    """

    def __init__(self, param, domain, x_average, options, phase="primary"):
        super().__init__(param, domain, options, phase)
        self.x_average = x_average

        pybamm.citations.register("Ai2019")
        pybamm.citations.register("Deshpande2012")

    def get_fundamental_variables(self):
        domain, Domain = self.domain_Domain
        phase_name = self.phase_param.phase_name
        if self.x_average:
            if self.size_distribution:
                l_cr_av_dist = pybamm.Variable(
                    f"X-averaged {domain} {phase_name}particle crack length distribution [m]",
                    domains={
                        "primary": f"{domain} particle size",
                        "secondary": "current collector",
                    },
                    scale=self.phase_param.l_cr_0,
                )
                l_cr_dist = pybamm.SecondaryBroadcast(
                    l_cr_av_dist, f"{domain} electrode"
                )
                l_cr_av = pybamm.size_average(l_cr_av_dist)
            else:
                l_cr_av = pybamm.Variable(
                    f"X-averaged {domain} {phase_name}particle crack length [m]",
                    domain="current collector",
                    scale=self.phase_param.l_cr_0,
                )
            l_cr = pybamm.PrimaryBroadcast(l_cr_av, f"{domain} electrode")
        else:
            if self.size_distribution:
                l_cr_dist = pybamm.Variable(
                    f"{Domain} {phase_name}particle crack length distribution [m]",
                    domains={
                        "primary": f"{domain} particle size",
                        "secondary": f"{domain} electrode",
                        "tertiary": "current collector",
                    },
                    scale=self.phase_param.l_cr_0,
                )
                l_cr = pybamm.size_average(l_cr_dist)
            else:
                l_cr = pybamm.Variable(
                    f"{Domain} {phase_name}particle crack length [m]",
                    domain=f"{domain} electrode",
                    auxiliary_domains={"secondary": "current collector"},
                    scale=self.phase_param.l_cr_0,
                )

        variables = self._get_standard_variables(l_cr)
        if self.size_distribution:
            variables.update(self._get_standard_size_distribution_variables(l_cr_dist))
        if self._amplitude_gated:
            variables.update(self._get_strain_envelope_variables())

        return variables

    @property
    def _amplitude_gated(self):
        phase_options = getattr(getattr(self.options, self.domain), self.phase)
        return phase_options["crack closure"] in ("strain amplitude", "strain amplitude power")

    def _get_strain_envelope_variables(self):
        """"crack closure" = "strain amplitude": the upper / lower envelopes of the phase's
        x-averaged volumetric strain eps_V = ln(V/V0), their difference (the recent strain
        amplitude A) and the opening fraction omega(A) the crack-face SEI and growth take."""
        domain, Domain = self.domain_Domain
        phase_name = self.phase_param.phase_name
        hi = pybamm.Variable(
            f"X-averaged {domain} {phase_name}particle volumetric strain upper envelope",
            domain="current collector",
        )
        lo = pybamm.Variable(
            f"X-averaged {domain} {phase_name}particle volumetric strain lower envelope",
            domain="current collector",
        )
        amplitude = hi - lo
        phase_options = getattr(getattr(self.options, domain), self.phase)
        if phase_options["crack closure"] == "strain amplitude power":
            # crack-face SEI activity (A / A_ref)^n; the growth is not gated by it
            opening = self.phase_param.crack_activity_amplitude_power(amplitude)
        else:
            opening = self.phase_param.crack_opening_amplitude(amplitude)
        return {
            f"X-averaged {domain} {phase_name}particle volumetric strain upper envelope": hi,
            f"X-averaged {domain} {phase_name}particle volumetric strain lower envelope": lo,
            f"X-averaged {domain} {phase_name}particle volumetric strain amplitude": amplitude,
            f"X-averaged {domain} electrode {phase_name}crack opening fraction": opening,
            f"{Domain} electrode {phase_name}crack opening fraction": pybamm.PrimaryBroadcast(
                opening, f"{domain} electrode"
            ),
        }

    def _volume_ratio_function(self, variables, x_avg=False, bol=False):
        """t(sto) = V/V0 of the phase: the fitted 1 + 3*sto**n_eff (n_eff moving with the
        LAM fraction) when "volume change aging deformation" is on with a porosity LAM
        option, else 1 + t_change(sto). x_avg: n_eff from the x-averaged volume fraction;
        bol: n_eff at no LAM (state-free, for initial conditions)."""
        domain, Domain = self.domain_Domain
        phase_name = self.phase_param.phase_name
        phase_param = self.phase_param
        phase_options = getattr(getattr(self.options, domain), self.phase)
        lam_option = phase_options["loss of active material"]
        if phase_options["volume change aging deformation"] == "true" and (
            "porosity" in lam_option
        ) and bol:
            n_bol = phase_param.volume_change_deform_exponent_bol

            def t_of(s):
                return 1 + 3 * s**n_bol
        elif phase_options["volume change aging deformation"] == "true" and (
            "porosity" in lam_option
        ):
            eps_s = variables[
                f"{Domain} electrode {phase_name}active material volume fraction"
            ]
            if x_avg:
                eps_s = pybamm.x_average(eps_s)
            exp_bol = phase_param.volume_change_deform_exponent_bol
            exp_end = phase_param.volume_change_deform_exponent_end
            lam_frac = pybamm.minimum(
                pybamm.maximum(1 - eps_s / pybamm.x_average(phase_param.epsilon_s)
                               if x_avg else 1 - eps_s / phase_param.epsilon_s, 0), 1
            )
            n_eff = exp_bol + lam_frac * (exp_end - exp_bol)

            def t_of(s):
                return 1 + 3 * s**n_eff
        else:
            # The literature t_change(sto) functions give the relative volume
            # CHANGE dV/V0 (e.g. ~0-0.1 for graphite), whereas the fitted law
            # above is the volume RATIO V/V0; the log strain needs the ratio.

            def t_of(s):
                return 1 + phase_param.t_change(s)

        return t_of

    def get_coupled_variables(self, variables):
        domain, Domain = self.domain_Domain
        phase_name = self.phase_param.phase_name
        variables.update(self._get_standard_surface_variables(variables))
        variables.update(self._get_mechanical_results(variables))
        if self.size_distribution:
            variables.update(self._get_mechanical_size_distribution_results(variables))
        T = variables[f"{Domain} electrode temperature [K]"]
        k_cr = self.phase_param.k_cr(T)
        m_cr = self.phase_param.m_cr
        b_cr = self.phase_param.b_cr
        if self.size_distribution:
            stress_t_surf = variables[
                f"{Domain} {phase_name}particle surface tangential stress distribution [Pa]"
            ]
        else:
            stress_t_surf = variables[
                f"{Domain} {phase_name}particle surface tangential stress [Pa]"
            ]
        if self.size_distribution:
            l_cr = variables[
                f"{Domain} {phase_name}particle crack length distribution [m]"
            ]
        else:
            l_cr = variables[f"{Domain} {phase_name}particle crack length [m]"]
        # # compressive stress will not lead to crack propagation
        # dK_SIF = stress_t_surf * b_cr * pybamm.sqrt(np.pi * l_cr) * (stress_t_surf >= 0)
        # dl_cr = k_cr * (dK_SIF**m_cr) / 3600  # divide by 3600 to replace t0_cr
        # TODO: Anti-Paris law update below (uncomment below, and use above for full Paris-Law):
        R_typ = self.phase_param.R_typ
        stress_relief = pybamm.maximum(1 - l_cr / R_typ, 0)
        stress_eff = stress_t_surf * stress_relief
        dK_SIF = stress_eff * b_cr * pybamm.sqrt(np.pi * l_cr) * (stress_eff >= 0)
        dl_cr = k_cr * (dK_SIF ** m_cr) / 3600
        # The strain-fatigue term (if on) is added in set_rhs, not here. It
        # needs the particle rhs, which (with stress-driven diffusion) needs
        # this submodel's stresses. A KeyError here would defer this submodel
        # in build_coupled_variables AFTER its thickness variables were
        # written to the shared dict, and the retry would then run after the
        # pore-buffering porosity submodel and overwrite its buffered
        # thickness change with the unbuffered one (found 2026-10-04: the
        # CELL009 expansion hump vanished with route B on).
        variables.update(
            {
                f"{Domain} {phase_name}particle cracking rate [m.s-1]": dl_cr,
                f"X-averaged {domain} {phase_name}particle cracking rate [m.s-1]": pybamm.x_average(
                    dl_cr
                ),
            }
        )
        return variables

    def _strain_fatigue_cracking_rate(self, variables, l_cr):
        """Volume-change (strain) fatigue crack growth, added on top of the
        Paris law when "particle cracking growth" includes "strain fatigue":

            dl/dt|_V = k_V(T) * D(d eps_V/dt) * l * (1 - l/R),
            eps_V = ln t(sto_rav),
            d eps_V/dt = (dt/dsto / t) * d sto_rav/dt,

        with D = |.| ("Paris + strain fatigue") or max(-., 0) (contraction
        only). t is the phase's volume RATIO V/V0: the fitted 1 + 3*sto**n_eff
        when "volume change aging deformation" is on (as in base_mechanics),
        else 1 + t_change(sto) (literature t_change is dV/V0). d sto_rav/dt is the r-average of the
        particle rhs over c_max; dt/dsto is a central difference (avoids the
        sto**(n-1) singularity of the power law at sto = 0). Called from
        set_rhs, where every coupled variable exists; it also overwrites the
        total cracking-rate variables with Paris + strain. Returns None
        when the option is "Paris" (default)."""
        domain, Domain = self.domain_Domain
        phase_name = self.phase_param.phase_name
        phase_options = getattr(getattr(self.options, domain), self.phase)
        growth = phase_options["particle cracking growth"]
        if growth == "Paris":
            return None
        if self.size_distribution or self.x_average:
            raise NotImplementedError(
                "strain-fatigue cracking is only implemented for x-resolved "
                "models without particle-size distributions"
            )
        phase_param = self.phase_param
        T = variables[f"{Domain} electrode temperature [K]"]
        R_typ = phase_param.R_typ
        sto_rav = variables[f"R-averaged {domain} {phase_name}particle concentration"]
        rhs = variables[f"{Domain} {phase_name}particle rhs [mol.m-3.s-1]"]
        dsto_dt = pybamm.r_average(rhs) / phase_param.c_max
        t_of = self._volume_ratio_function(variables)

        # sto is clamped to [s_floor, 1]: the solver can overshoot past 0 at
        # the end of a deep discharge (s**n -> NaN), and at exactly 0 the
        # Jacobian of s**n_eff w.r.t. the (LAM-dependent) exponent,
        # s**n * ln(s), is 0 * -inf = NaN, which stalls IDA.
        delta = 1e-3
        s_floor = 1e-6
        s_c = pybamm.minimum(pybamm.maximum(sto_rav, s_floor), 1)
        s_hi = pybamm.minimum(s_c + delta, 1)
        s_lo = pybamm.maximum(s_c - delta, s_floor)
        dt_ds = (t_of(s_hi) - t_of(s_lo)) / (s_hi - s_lo)
        strain_rate = dt_ds / t_of(s_c) * dsto_dt
        if growth == "Paris + strain fatigue (contraction)":
            driver = pybamm.maximum(-strain_rate, 0)
        else:
            driver = pybamm.AbsoluteValue(strain_rate)
        dl_cr_strain = (
            phase_param.k_V(T) * driver * l_cr * pybamm.maximum(1 - l_cr / R_typ, 0)
        )
        variables.update(
            {
                f"{Domain} {phase_name}particle volumetric strain rate [s-1]": strain_rate,
                f"X-averaged {domain} {phase_name}particle volumetric strain rate [s-1]": pybamm.x_average(
                    strain_rate
                ),
                f"X-averaged {domain} {phase_name}particle strain-fatigue driver [s-1]": pybamm.x_average(
                    driver
                ),
                f"{Domain} {phase_name}particle strain-fatigue cracking rate [m.s-1]": dl_cr_strain,
                f"X-averaged {domain} {phase_name}particle strain-fatigue cracking rate [m.s-1]": pybamm.x_average(
                    dl_cr_strain
                ),
            }
        )
        return dl_cr_strain

    def set_rhs(self, variables):
        domain, Domain = self.domain_Domain
        phase_name = self.phase_param.phase_name
        if self.x_average is True:
            if self.size_distribution:
                l_cr = variables[
                    f"X-averaged {domain} {phase_name}particle crack length distribution [m]"
                ]
            else:
                l_cr = variables[
                    f"X-averaged {domain} {phase_name}particle crack length [m]"
                ]
            dl_cr = variables[
                f"X-averaged {domain} {phase_name}particle cracking rate [m.s-1]"
            ]
        else:
            if self.size_distribution:
                l_cr = variables[
                    f"{Domain} {phase_name}particle crack length distribution [m]"
                ]
            else:
                l_cr = variables[f"{Domain} {phase_name}particle crack length [m]"]
            dl_cr = variables[f"{Domain} {phase_name}particle cracking rate [m.s-1]"]
        dl_cr_strain = self._strain_fatigue_cracking_rate(variables, l_cr)
        phase_options = getattr(getattr(self.options, domain), self.phase)
        gated = phase_options["crack closure"] in ("SEI and growth", "strain amplitude")
        if dl_cr_strain is not None or gated:
            if dl_cr_strain is not None:
                dl_cr = dl_cr + dl_cr_strain
            if gated:  # "crack closure": a crack held shut does not propagate
                if self._amplitude_gated:
                    opening = variables[
                        f"X-averaged {domain} electrode {phase_name}crack opening fraction"
                    ]
                else:
                    sto_av = variables[
                        f"R-averaged {domain} {phase_name}particle stoichiometry"
                    ]
                    if self.x_average is True:
                        sto_av = pybamm.x_average(sto_av)
                    opening = self.phase_param.crack_opening(sto_av)
                # growth residual g: closed cracks still propagate at g x the open rate
                g = self.phase_param.crack_closure_growth_residual
                dl_cr = dl_cr * (g + (1 - g) * opening)
            variables.update(
                {
                    f"{Domain} {phase_name}particle cracking rate [m.s-1]": dl_cr,
                    f"X-averaged {domain} {phase_name}particle cracking rate [m.s-1]": pybamm.x_average(
                        dl_cr
                    ),
                }
            )
        self.rhs = {l_cr: dl_cr}
        if self._amplitude_gated:
            self.rhs.update(self._strain_envelope_rhs(variables))

    # rise rate of the strain envelopes onto a new peak / trough [s-1] (~100 s)
    ENVELOPE_RISE_RATE = 1e-2

    def _x_averaged_volume_strain(self, variables, sto=None, bol=False):
        domain = self.domain
        phase_name = self.phase_param.phase_name
        if sto is None:
            sto = pybamm.x_average(
                variables[f"R-averaged {domain} {phase_name}particle stoichiometry"]
            )
        t_of = self._volume_ratio_function(variables, x_avg=True, bol=bol)
        return pybamm.log(t_of(pybamm.minimum(pybamm.maximum(sto, 1e-6), 1)))

    def _strain_envelope_rhs(self, variables):
        """d eps_hi/dt = k max(eps - eps_hi, 0) - (eps_hi - eps)/tau, and the mirror image for
        eps_lo: each envelope jumps onto a new peak (trough) within ~1/k and relaxes towards
        the current strain over the memory time tau."""
        domain = self.domain
        phase_name = self.phase_param.phase_name
        hi = variables[
            f"X-averaged {domain} {phase_name}particle volumetric strain upper envelope"
        ]
        lo = variables[
            f"X-averaged {domain} {phase_name}particle volumetric strain lower envelope"
        ]
        eps = self._x_averaged_volume_strain(variables)
        k, tau = self.ENVELOPE_RISE_RATE, self.phase_param.crack_closure_memory_time
        return {
            hi: k * pybamm.maximum(eps - hi, 0) - (hi - eps) / tau,
            lo: -k * pybamm.maximum(lo - eps, 0) + (eps - lo) / tau,
        }

    def set_initial_conditions(self, variables):
        domain, Domain = self.domain_Domain
        phase_name = self.phase_param.phase_name
        l_cr_0 = self.phase_param.l_cr_0
        if self.x_average is True:
            if self.size_distribution:
                l_cr = variables[
                    f"X-averaged {domain} {phase_name}particle crack length distribution [m]"
                ]
                l_cr_0 = pybamm.PrimaryBroadcast(l_cr_0, f"{domain} particle size")
            else:
                l_cr = variables[
                    f"X-averaged {domain} {phase_name}particle crack length [m]"
                ]
        else:
            if self.size_distribution:
                l_cr = variables[
                    f"{Domain} {phase_name}particle crack length distribution [m]"
                ]
                l_cr_0 = pybamm.PrimaryBroadcast(l_cr_0, f"{domain} electrode")
                l_cr_0 = pybamm.PrimaryBroadcast(l_cr_0, f"{domain} particle size")
            else:
                l_cr = variables[f"{Domain} {phase_name}particle crack length [m]"]
                l_cr_0 = pybamm.PrimaryBroadcast(l_cr_0, f"{domain} electrode")
        self.initial_conditions = {l_cr: l_cr_0}
        if self._amplitude_gated:  # both envelopes start at the initial strain: A = 0
            eps0 = self._x_averaged_volume_strain(
                variables, sto=self.phase_param.sto_init_av, bol=True
            )
            for name in ("upper", "lower"):
                v = variables[
                    f"X-averaged {domain} {phase_name}particle volumetric strain {name} envelope"
                ]
                self.initial_conditions[v] = eps0

    def add_events_from(self, variables):
        domain, Domain = self.domain_Domain
        phase_name = self.phase_param.phase_name
        if self.x_average is True:
            l_cr = variables[
                f"X-averaged {domain} {phase_name}particle crack length [m]"
            ]
        else:
            l_cr = variables[f"{Domain} {phase_name}particle crack length [m]"]
        self.events.append(
            pybamm.Event(
                f"{domain} {phase_name} particle crack length larger than particle radius",
                1 - pybamm.max(l_cr) / self.phase_param.R_typ,
            )
        )
