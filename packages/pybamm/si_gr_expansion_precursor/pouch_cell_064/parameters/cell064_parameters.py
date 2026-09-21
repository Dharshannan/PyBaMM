"""CELL064 (GMFEB23S NMC622 / Si-Gr composite pouch cell) parameter script.

Not a registered PyBaMM parameter set (no entry point) -- a plain module with
a ``get_parameter_values()`` function, per the project's own convention for
one-off cell parameterisations. Usage::

    import pybamm
    from cell064_parameters import get_parameter_values

    param = pybamm.ParameterValues(get_parameter_values())

Built as ``si_gr_expansion`` (this project's own composite Si/graphite
degradation recipe -- see
``src/pybamm/input/parameters/lithium_ion/si_gr_expansion.py``) with CELL064's
own measured/fitted BoL values layered on top via ``dict.update()``, rather
than a from-scratch parameter set. This keeps every degradation submodel
default (SEI-on-cracks, plating, stress-driven LAM, pore buffering, the
project's stability-fix parameters) exactly as already validated -- only the
BoL-specific entries below (geometry, composition, OCPs, capacity, voltage
window) are cell064-specific.

Source and provenance
----------------------
CELL064's own BoL fit lives in a separate project
(``Si_Gr_Expansion_Precursor/Pouch_Data/model_cell064/``), specifically the
"FINAL v6" parameter set
(``pybamm_parameters/CELL064_BoL_parameter_set_hysteresis_cell076transplant.py``,
built by ``stage2_build_pybamm_parameters.py`` + ``stage8b_...`` from the
RPT001 discharge fit, EFC~50.6 -- not pristine BoL, but the earliest
QC-valid, internally-consistent checkpoint CELL064 has at both C/20 and C/3
rates; RPT000's true-BoL C/20 discharge failed QC). Every value below is
copied from that file's ``get_parameter_values()``, with provenance tags
preserved in comments ([DATA] = measured for this cell, [SOLVED GEOMETRY] =
solved to match measured capacity/area, [DATA x LITERATURE] = measured
stoichiometry x assumed c_s_max).

What is NOT overridden (and why): particle radii, Young's moduli, critical
stresses, LAM proportional/exponential terms, cracking rates, and SEI/
plating kinetics all carry the "[LITERATURE - TO FIT]" tag in the source
file -- CELL064's own Stage 3 template copied these UNCHANGED from the same
Chen2020/Ai2020/Ai2022/OKane2022 literature family ``si_gr_expansion``
already uses, and cross-checking confirms the numeric values are identical
to ``si_gr_expansion``'s own defaults (e.g. both use 5.86e-6 m graphite
particle radius, 720 MPa silicon critical stress, 50 GPa silicon Young's
modulus). These are exactly the parameters this project's degradation test
matrix (``degradation_test_matrix/``) exists to re-tune against CELL064's
own real capacity-fade data -- inheriting ``si_gr_expansion``'s literature
placeholders here is the correct starting point, not an oversight.
``si_gr_expansion``'s own retuned positive-electrode partial molar volume
(1.58e-6, vs the naive 1.25e-5 literature default CELL064's template still
carries -- see si_gr_expansion.py's inline comment) and its silicon
volume-change function (measured V/V0 curve vs. the analytic Ai2020 CELL064
still uses) are likewise kept as improvements, not reverted.

Not yet retuned for CELL064: pore-buffering closure porosity / transition
width / plateau (still ``si_gr_expansion``'s own values, tuned against a
different cell's synthetic recipe -- refit against CELL064's real reversible
expansion in ``expansion_pattern_test/``, per ``../README.md``'s pipeline
order).

Model options: unlike ``si_gr_expansion``'s own test scripts, CELL064's own
BoL fit activates ``"contact resistance": "true"`` (with base value 0 Ohm,
fit separately in stage 3). This is NOT enabled by default here -- it is
an optional resistance-growth lever, not part of the existing LLI/LAM
degradation coupling this pipeline reuses; add
``"contact resistance": "true"`` to a simulation's ``MODEL_OPTIONS`` only if
fitting voltage/resistance growth explicitly, and give
``"Contact resistance [Ohm]"`` a nonzero value if so (0.0 here is inert).
"""

import os

import pybamm
from pybamm.input.parameters.lithium_ion.si_gr_expansion import (
    get_parameter_values as _si_gr_expansion_values,
)

DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def _process_1d(csv_name):
    return pybamm.parameters.process_1D_data(csv_name, path=DATA_PATH)


_graphite_ocp_data = _process_1d("graphite_ocp_bol_cell064_native.csv")
_silicon_ocp_delithiation_data = _process_1d("silicon_ocp_bol_deformed_cell064.csv")
_silicon_ocp_lithiation_data = _process_1d(
    "silicon_ocp_bol_lithiation_deformed_cell064_cell076transplant.csv"
)
_nmc622_ocp_data = _process_1d("nmc622_ocp_bol_cell064_native.csv")


def graphite_ocp_cell064(sto):
    """Negative electrode (graphite) OCP: CELL064's own native delithiation-
    branch measurement. [DATA]"""
    name, (x, y) = _graphite_ocp_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


def silicon_ocp_delithiation_cell064(sto):
    """Silicon delithiation-branch OCP: CELL064's own native measurement,
    deformed by this cell's own fitted U_si_new = scale*U_si_base + shift
    (scale=0.96750, shift=-0.01720 V). [DATA]"""
    name, (x, y) = _silicon_ocp_delithiation_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


def silicon_ocp_lithiation_cell064(sto):
    """Silicon lithiation-branch OCP: transplanted from CELL076's own real,
    charge-calibrated lithiation curve (same material family), flat-
    extrapolated onto CELL064's padded stoichiometry domain -- CELL064's own
    narrower charge/discharge data could not resolve this branch directly.
    [DATA, cross-cell transplant]"""
    name, (x, y) = _silicon_ocp_lithiation_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


def silicon_ocp_average_cell064(sto):
    return (silicon_ocp_lithiation_cell064(sto) + silicon_ocp_delithiation_cell064(sto)) / 2


def nmc622_ocp_cell064(sto):
    """Positive electrode (NMC622) OCP: CELL064's own native lithiation-
    branch measurement. [DATA]"""
    name, (x, y) = _nmc622_ocp_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


def get_parameter_values():
    """``si_gr_expansion`` with CELL064's own BoL(-ish, RPT001, EFC~50.6)
    geometry, composition, OCPs, capacity and voltage window layered on top.

    Key numbers for this cell (RPT001 C/20 discharge, eSOH reference):
        Cn_Si = 1.2618 Ah, Cn_Gr = 1.3329 Ah, Cn = 2.5947 Ah, Cp = 2.7432 Ah
        Si share of anode capacity = 0.4863 (high-silicon design)
        Electrode area (measured) = 69.69 cm^2; (solved, geometry convention
        matching CELL076/si_gr_expansion) = 343.48 cm^2
    """
    values = _si_gr_expansion_values()
    values.update(
        {
            # --- cell geometry: measured/solved for THIS cell -------------
            "Positive electrode thickness [m]": 8.95119348e-05,  # [SOLVED GEOMETRY]
            "Electrode width [m]": 0.52843063,  # [DATA] = measured area / assumed height
            "Cell cooling surface area [m2]": 0.00170857,  # [SCALED from literature reference cell]
            "Cell volume [m3]": 7.7867e-06,  # [SCALED from literature reference cell]
            "Nominal cell capacity [A.h]": 2.5947,  # [DATA] Cn_Gr + Cn_Si
            "Current function [A]": 2.5947,  # [DATA] 1C default; override per experiment
            # Negative current collector/separator/positive current collector
            # thickness, electrode height, contact resistance base value and
            # both mechanical [LITERATURE - TO FIT] blocks are identical to
            # si_gr_expansion's own defaults already -- no override needed.

            # --- negative electrode: measured composition + OCPs ----------
            "Negative electrode porosity": 0.350000,  # [SOLVED GEOMETRY]
            "Primary: Negative electrode active material volume fraction": 0.592132,  # [SOLVED GEOMETRY]
            "Primary: Initial concentration in negative electrode [mol.m-3]": 25676.9209,  # [DATA x LITERATURE]
            "Primary: Negative electrode OCP [V]": graphite_ocp_cell064,  # [DATA]
            "Secondary: Negative electrode active material volume fraction": 0.0578685,  # [SOLVED GEOMETRY]
            "Secondary: Initial concentration in negative electrode [mol.m-3]": 276020.0224,  # [DATA x LITERATURE]
            "Secondary: Negative electrode lithiation OCP [V]": silicon_ocp_lithiation_cell064,
            "Secondary: Negative electrode delithiation OCP [V]": silicon_ocp_delithiation_cell064,
            "Secondary: Negative electrode OCP [V]": silicon_ocp_average_cell064,
            # Secondary: Maximum concentration (278000 mol/m3) matches
            # si_gr_expansion's own default exactly -- no override needed.

            # --- positive electrode: measured OCP + capacity ---------------
            "Positive electrode OCP [V]": nmc622_ocp_cell064,  # [DATA]
            "Maximum concentration in positive electrode [mol.m-3]": 50060.0,  # [LITERATURE - TO FIT]
            # Positive electrode porosity/active material volume fraction/
            # particle radius are already identical to si_gr_expansion's
            # defaults -- no override needed.

            # --- voltage window: measured cycler cutoffs -------------------
            "Upper voltage cut-off [V]": 4.195,  # [DATA] protocol CV cutoff
            "Open-circuit voltage at 100% SOC [V]": 4.195,  # [DATA]
            # Lower cutoff / OCV at 0% SOC (2.5 V) already match
            # si_gr_expansion's own defaults -- no override needed.

            # --- composite-average fallback initial concentrations --------
            # Not authoritative when "Primary:"/"Secondary:" keys are set
            # (composite negative electrode with 2 particle phases), but set
            # for consistency with the source parameter set.
            "Initial concentration in negative electrode [mol.m-3]": 25676.9209,
            "Initial concentration in positive electrode [mol.m-3]": 262.5541,  # [DATA x LITERATURE]
        }
    )
    return values
