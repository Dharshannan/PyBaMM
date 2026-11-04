"""Model options and the shared parameter set for all four conditions.

Base set: ``si_gr_expansion`` (in the modified PyBaMM) with CELL064's
beginning-of-life composition, OCPs and voltage window on top (from its RPT1
C/20 fit), in the cells' design geometry (DESIGN below). The degradation
parameters follow (V7). Only temperature, F0 and the pore-buffering width
differ between conditions (conditions.py). The Arrhenius factors are written
out exactly as fitted.
Note the two sign conventions: crack-growth activation energies use
exp[Ea/R (1/T - 1/T_ref)] (positive = slower when hot, as in Ai2020); all
others use exp[Ea/R (1/T_ref - 1/T)] (positive = faster when hot).
"""
import os

import numpy as np
import pybamm
from pybamm.input.parameters.lithium_ion.si_gr_expansion import (
    electrolyte_conductivity_Nyman2008_arrhenius,
    electrolyte_diffusivity_Nyman2008_arrhenius,
    get_parameter_values as si_gr_expansion_values,
    silicon_LGM50_diffusivity_Bonkile2024,
    silicon_LGM50_electrolyte_exchange_current_density_Chen2020,
)

from .conditions import MODEL_CAPACITY_AH, RATE_BASIS_AH, UPPER_CUTOFF_V

T_REF = 298.15

MODEL_OPTIONS = {
    "particle phases": ("2", "1"),  # negative: graphite (primary) + Si (secondary)
    "open-circuit potential": (("single", "current sigmoid"), "single"),  # Si OCP hysteresis
    "SEI": (("ec reaction limited", "solvent-diffusion limited"), "none"),
    "SEI porosity change": "true",
    "SEI film resistance": "distributed",
    "SEI on cracks": "true",
    "particle mechanics": ("swelling and cracking", "swelling only"),
    "particle cracking growth": (("Paris", "Paris + strain fatigue (contraction)"), "Paris"),
    "loss of active material": (
        ("stress-driven and porosity isolation", "stress and reaction-driven and porosity isolation"),
        "stress-driven",
    ),
    "stress-driven LAM damping": (("linear", "power"), "power"),
    "porosity isolation gate": "dose",
    "SEI reaction redirect to LAM": "false",  # no effect with the dose gate; stated for clarity
    "isolation lithium trapping": (("false", "true"), "false"),
    "pore buffering": "true",
    "pore buffering transition": "physical",
    "open-circuit potential aging deformation": (("false", "true"), "false"),
    "OCP aging deformation driver": "LAM fraction",
    "volume change aging deformation": (("false", "true"), "false"),
    "thermal": "isothermal",
}

VAR_PTS = {"x_n": 5, "x_s": 5, "x_p": 5, "r_n_prim": 20, "r_n_sec": 20, "r_n": 20, "r_p": 20}

# ---------------------------------------------------------------------------
# Design geometry [Wan2025 Table 2; Wan2026 Table S4]: the stack (5 double-sided
# positive, 6 negative sheets) as one lumped electrode pair with the facing area
# of the 10 positive faces (69 x 101 mm). Per-face thicknesses estimated from
# the double-sided loadings at the stated 30 % porosity (the papers give none).
# Every active phase keeps the solid volume of the BOL capacity fit (made in an
# area-fitted geometry, BOL_FIT below), so capacities, particle areas and
# particle fluxes are those of that fit: the volume fractions scale by the
# ratio of the electrode volumes.
# ---------------------------------------------------------------------------
DESIGN = dict(area=10 * 0.069 * 0.101, height=0.065, L_n=52e-6, L_p=76e-6, pos_porosity=0.30)
BOL_FIT = dict(width=0.52843063, L_n=85.2e-6, L_p=8.95119348e-05)  # geometry of the BOL capacity fit
SOLID_SCALE_N = BOL_FIT["width"] * DESIGN["height"] * BOL_FIT["L_n"] / (DESIGN["area"] * DESIGN["L_n"])
SOLID_SCALE_P = BOL_FIT["width"] * DESIGN["height"] * BOL_FIT["L_p"] / (DESIGN["area"] * DESIGN["L_p"])

# Effective (porous) SEI: the SEI partial molar volume, EC diffusivity and initial thicknesses of the base set
# times SEI_VOLUME_FACTOR, so the pores of the design geometry (porosity 0.30) fill to their floor with the
# lithium consumption of the BOL-fit geometry; the SEI resistivity and the solvent diffusivities are scaled
# with it in V7 (same resistance and same growth per unit lithium).
SEI_VOLUME_FACTOR = 2.602

# ---------------------------------------------------------------------------
# Shared degradation parameters. Values at T_REF; activation energies [J/mol]
# give the temperature dependence. Refitted 2026-11-04 in the design geometry
# with the 2.5 Ah basis on all four cells (stress-LAM rates are absolute rates,
# so they scale with the solid volume fractions).
# ---------------------------------------------------------------------------
V7 = dict(
    # Si cracking, route A (Paris law on the diffusion-induced surface stress)
    si_crack_rate=2.709e-17, si_crack_Ea=90000.0,  # positive = slower when hot
    # Si cracking, route B (contraction-only strain fatigue)
    si_strain_kV=2.3e-3 * 1.15, si_strain_Ea=150000.0,  # positive = slower when hot
    # Si stress- and reaction-driven LAM
    si_lam_beta=1.274e-08, si_lam_Ea=40000.0, si_lam_m=2.5, si_critical_stress=6.0e6,
    si_lam_damping_exponent=0.25, si_reaction_lam_factor=1e-7,
    # Si SEI (solvent-diffusion limited)
    si_sei_solvent_diffusivity=3.253e-21, si_sei_Ea=34450.0, si_sei_k=3e-16,
    # Porosity-gated isolation LAM, dose form: k_iso = Y0 exp[E_Y/R(1/T_ref - 1/T)] / tau
    si_iso_exponent=20.0, si_iso_Y0=0.3054, si_iso_E_Y=22000.0, si_iso_tau=3.5e8,
    gr_iso_exponent=10.0, gr_iso_Y0=0.0227, gr_iso_E_Y=61210.0, gr_iso_tau=2e8,
    si_iso_li_trapping=1.0,
    # Si aging deformation (fitted Si OCP deformation and volume-change law)
    si_ocp_deform_scale=0.8663, si_ocp_deform_shift=-0.00908,
    si_volume_exp_bol=1.561216606871758, si_volume_exp_end=2.18807468386046,
    # Graphite
    gr_lam_beta=2.299e-07, gr_lam_m=2.0, gr_critical_stress=3.0e7, gr_crack_rate=3.9e-21,
    gr_sei_Ea=32000.0, gr_sei_k=1.546e-17, gr_sei_solvent_diffusivity=6.506e-22,  # D_sol used by the solvent-diffusion law only
    # Positive electrode stress-driven LAM (saturating "power" damping)
    pos_lam_beta=1.922e-05, pos_lam_Ea=-55000.0, pos_lam_damping_exponent=50.0,
    # SEI, porosity, pore buffering
    sei_resistivity=7.686e4, sei_exponent_cap=70.0, neg_porosity=0.30, neg_porosity_floor=0.0163,
    pore_closure_porosity=0.0793,
    # Other
    si_particle_radius=1.52e-6 * 0.2, neg_density=1745.4, c_e0=1000.0, si_diffusivity_mult=3.0,
)


def _arrhenius_constant(Ea, T):
    """Isothermal factor exp[Ea/R (1/T_ref - 1/T)], evaluated once at the
    operating temperature (exactly 1 at 25 degC)."""
    return float(np.exp(Ea / 8.314462618 * (1 / T_REF - 1 / T)))


def _ocp_functions(data_dir):
    ocp_dir = os.path.join(data_dir, "ocp")

    def interp(csv_name):
        name, (x, y) = pybamm.parameters.process_1D_data(csv_name, path=ocp_dir)
        return lambda sto: pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")

    graphite = interp("graphite_ocp_bol_cell064_native.csv")
    si_delith = interp("silicon_ocp_bol_deformed_cell064.csv")
    si_lith = interp("silicon_ocp_bol_lithiation_deformed_cell064_cell076transplant.csv")
    nmc = interp("nmc622_ocp_bol_cell064_native.csv")
    return graphite, si_lith, si_delith, (lambda sto: (si_lith(sto) + si_delith(sto)) / 2), nmc


def base_parameter_values(data_dir):
    """si_gr_expansion + CELL064's beginning-of-life cell data, in the design geometry."""
    graphite, si_lith, si_delith, si_avg, nmc = _ocp_functions(data_dir)
    values = si_gr_expansion_values()
    pos_fraction = values["Positive electrode active material volume fraction"] * SOLID_SCALE_P
    values.update({
        "Electrode width [m]": DESIGN["area"] / DESIGN["height"],
        "Electrode height [m]": DESIGN["height"],
        "Negative electrode thickness [m]": DESIGN["L_n"],
        "Positive electrode thickness [m]": DESIGN["L_p"],
        "Positive electrode porosity": DESIGN["pos_porosity"],
        "Positive electrode active material volume fraction": pos_fraction,
        "Cell cooling surface area [m2]": 0.00170857,
        "Cell volume [m3]": 7.7867e-06,
        "Nominal cell capacity [A.h]": MODEL_CAPACITY_AH,
        "Current function [A]": MODEL_CAPACITY_AH,
        "Negative electrode porosity": 0.35,
        "Primary: Negative electrode active material volume fraction": 0.592132 * SOLID_SCALE_N,
        "Primary: Initial concentration in negative electrode [mol.m-3]": 25676.9209,
        "Primary: Negative electrode OCP [V]": graphite,
        "Secondary: Negative electrode active material volume fraction": 0.0578685 * SOLID_SCALE_N,
        "Secondary: Initial concentration in negative electrode [mol.m-3]": 276020.0224,
        "Secondary: Negative electrode lithiation OCP [V]": si_lith,
        "Secondary: Negative electrode delithiation OCP [V]": si_delith,
        "Secondary: Negative electrode OCP [V]": si_avg,
        "Positive electrode OCP [V]": nmc,
        "Maximum concentration in positive electrode [mol.m-3]": 50060.0,
        "Upper voltage cut-off [V]": UPPER_CUTOFF_V,
        "Open-circuit voltage at 100% SOC [V]": UPPER_CUTOFF_V,
        "Initial concentration in negative electrode [mol.m-3]": 25676.9209,
        "Initial concentration in positive electrode [mol.m-3]": 262.5541,
    })
    for name in ("Primary: SEI partial molar volume [m3.mol-1]", "Secondary: SEI partial molar volume [m3.mol-1]",
                 "Primary: EC diffusivity [m2.s-1]", "Secondary: EC diffusivity [m2.s-1]",
                 "Primary: Initial SEI thickness [m]", "Secondary: Initial SEI thickness [m]",
                 "Primary: Initial SEI on cracks thickness [m]", "Secondary: Initial SEI on cracks thickness [m]",
                 "Initial SEI thickness [m]"):
        values[name] = values[name] * SEI_VOLUME_FACTOR
    return values


def parameter_values(cond, data_dir, rate_basis_ah=RATE_BASIS_AH, p=V7):
    """pybamm.ParameterValues for one condition."""
    T = cond.temperature_K

    def si_crack_rate(T_dim):  # route A, positive Ea = slower when hot
        return p["si_crack_rate"] * np.exp(p["si_crack_Ea"] / pybamm.constants.R * (1 / T_dim - 1 / T_REF))

    def si_strain_kV(T_dim):  # route B, same convention
        return p["si_strain_kV"] * pybamm.exp(p["si_strain_Ea"] / 8.314 * (1 / T_dim - 1 / T_REF))

    def si_lam_beta(T_dim):
        return p["si_lam_beta"] * np.exp(p["si_lam_Ea"] / 8.314 * (1 / T_REF - 1 / T_dim))

    def si_diffusivity(sto, T_dim):
        return p["si_diffusivity_mult"] * silicon_LGM50_diffusivity_Bonkile2024(sto, T_dim)

    def k_iso(Y0, E_Y, tau):
        return Y0 * _arrhenius_constant(E_Y, T) / tau

    param = pybamm.ParameterValues(base_parameter_values(data_dir))
    param.update({
        "Nominal cell capacity [A.h]": rate_basis_ah,  # sets the A of "C/x" steps
        "Ambient temperature [K]": T,
        "Initial temperature [K]": T,
        "Initial concentration in electrolyte [mol.m-3]": p["c_e0"],
        "Negative electrode density [kg.m-3]": p["neg_density"],
        "Negative electrode porosity": p["neg_porosity"],
        "Negative electrode porosity floor": p["neg_porosity_floor"],
        "Electrolyte diffusivity [m2.s-1]": electrolyte_diffusivity_Nyman2008_arrhenius,
        "Electrolyte conductivity [S.m-1]": electrolyte_conductivity_Nyman2008_arrhenius,
        # --- Si (secondary phase) ---
        "Secondary: Negative particle radius [m]": p["si_particle_radius"],
        "Secondary: Negative particle diffusivity [m2.s-1]": si_diffusivity,
        "Secondary: Negative electrode exchange-current density [A.m-2]":
            silicon_LGM50_electrolyte_exchange_current_density_Chen2020,
        "Secondary: Negative electrode cracking rate": si_crack_rate,
        "Secondary: Negative electrode strain-fatigue cracking constant": si_strain_kV,
        "Secondary: Negative electrode critical stress [Pa]": p["si_critical_stress"],
        "Secondary: Negative electrode LAM constant proportional term [s-1]": si_lam_beta,
        "Secondary: Negative electrode LAM constant exponential term": p["si_lam_m"],
        "Secondary: Negative electrode stress-driven LAM damping exponent": p["si_lam_damping_exponent"],
        "Secondary: Negative electrode reaction-driven LAM factor [m3.mol-1]": p["si_reaction_lam_factor"],
        "Secondary: SEI solvent diffusivity [m2.s-1]": p["si_sei_solvent_diffusivity"],
        "Secondary: SEI growth activation energy [J.mol-1]": p["si_sei_Ea"],
        "Secondary: SEI kinetic rate constant [m.s-1]": p["si_sei_k"],
        "Secondary: SEI resistivity [Ohm.m]": p["sei_resistivity"],
        "Secondary: SEI reaction exponent cap": p["sei_exponent_cap"],
        "Secondary: Negative electrode porosity-isolation LAM exponent": p["si_iso_exponent"],
        "Secondary: Negative electrode porosity-isolation dose LAM factor [m3.mol-1.s-1]":
            k_iso(p["si_iso_Y0"], p["si_iso_E_Y"], p["si_iso_tau"]),
        "Secondary: Negative electrode isolation lithium trapping fraction": p["si_iso_li_trapping"],
        "Secondary: Negative electrode OCP aging-deformation end scale": p["si_ocp_deform_scale"],
        "Secondary: Negative electrode OCP aging-deformation end shift [V]": p["si_ocp_deform_shift"],
        "Secondary: Negative electrode volume change aging-deformation BOL exponent": p["si_volume_exp_bol"],
        "Secondary: Negative electrode volume change aging-deformation end exponent": p["si_volume_exp_end"],
        # --- graphite (primary phase) ---
        "Primary: Negative electrode LAM constant proportional term [s-1]": p["gr_lam_beta"],
        "Primary: Negative electrode LAM constant exponential term": p["gr_lam_m"],
        "Primary: Negative electrode critical stress [Pa]": p["gr_critical_stress"],
        "Primary: Negative electrode cracking rate": p["gr_crack_rate"],
        "Primary: SEI growth activation energy [J.mol-1]": p["gr_sei_Ea"],
        "Primary: SEI kinetic rate constant [m.s-1]": p["gr_sei_k"],
        "Primary: SEI solvent diffusivity [m2.s-1]": p["gr_sei_solvent_diffusivity"],
        "Primary: SEI resistivity [Ohm.m]": p["sei_resistivity"],
        "Primary: SEI reaction exponent cap": p["sei_exponent_cap"],
        "Primary: Negative electrode porosity-isolation LAM exponent": p["gr_iso_exponent"],
        "Primary: Negative electrode porosity-isolation dose LAM factor [m3.mol-1.s-1]":
            k_iso(p["gr_iso_Y0"], p["gr_iso_E_Y"], p["gr_iso_tau"]),
        # --- positive electrode ---
        "Positive electrode LAM constant proportional term [s-1]":
            p["pos_lam_beta"] * _arrhenius_constant(p["pos_lam_Ea"], T),
        "Positive electrode stress-driven LAM damping exponent": p["pos_lam_damping_exponent"],
        # --- pore buffering (width and F0 per condition) ---
        "Negative electrode pore buffering closure porosity": p["pore_closure_porosity"],
        "Negative electrode pore buffering transition width": cond.width,
        "Negative electrode transmitted fraction plateau": cond.F0,
    }, check_already_exists=False)
    return param
