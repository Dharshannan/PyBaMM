"""
dma_ocp_fit_final.py -- consolidated composite-DMA fit core: electrode OCP functions, the
pseudo-OCV reconstruction model, degradation-mode formulas, and the per-RPT nonlinear fit.
This is the single surviving fit module for this pipeline (formerly "Method C"/v2, ported
from the user's own working notebook nuDMA_GrSi_fresh 1.ipynb) -- the additive-cost
"Method B" (dma_ocp_fit.py) and the notebook-faithful delta-split-only "Method A"-style
comparisons that motivated this file's own design choices have been retired now that this
version's cost function, search strategy (see fit_rpt_v2's "global" search_style), and
constraint scheme have been validated as the best-performing of everything tried. The
electrode OCP functions and blend_anode_v_of_z below are inlined from the retired
dma_ocp_fit.py (unchanged) rather than imported, so this module has no other dependency.

WHAT MAKES THIS FIT'S DESIGN (still true, condensed from the original "Method B"
comparison that produced it -- kept as design record, not a live comparison):
  1. Parameterization: (nu_ne, nu_pe, sigma_ne, sigma_pe, phi_si, delta_v) -- a slope/
     intercept form (nu_e = x_hi-x_lo, sigma_e = x_lo) rather than PyProBE-style raw
     (x_lo, x_hi) limits, NOT re-bounded to [0,1] (the notebook's own bounds allow mild
     excursions: nu up to 1.05, sigma negative), matching the original notebook's own
     approach on this project's simulated data rather than a hybrid.
  2. delta_v is a flat voltage offset added to the PE curve -- algebraically -I_rpt*R for
     a constant per-RPT current, just without an explicit resistance parameter.
  3. A "delta-split" cost option (fit VOLTAGE on one side of the sharpest dV/dQ feature
     and dV/dQ, Huber-weighted, on the other) -- CURRENTLY INACTIVE in the validated
     config (dma_final_plot.py sets W_V_SPLIT=W_DVDQ_SPLIT=0), superseded by the
     whole-domain OCV+DV+dV/dQ cost (w_ocv_b/w_dv_b/w_full_dvdq) that turned out to
     track electrode attribution more reliably. Kept as a toggle, not deleted, since it's
     still a legitimate alternative cost worth revisiting.
  4. A soft V_min/V_max penalty, a feasibility pre-scan before differential_evolution, and
     TWO selectable search strategies (fit_rpt_v2's search_style): "restart" (adaptive
     shrinking-spread multi-restart, this file's own original scheme -- fine for RPT0-13,
     which converge to the same optimum regardless of search style or budget) and
     "global" (single/ensemble full-bounds Latin Hypercube search with wider mutation,
     RMSE-selected across independent seeds -- REQUIRED for RPT14, whose cost landscape
     has multiple competing basins that "restart" can never escape once seeded near the
     previous RPT's fit; see dma_final_plot.py's SEARCH_STYLE_LAST for the full story).

Degradation-mode formulas (compute_reservoirs / compute_degradation_modes_v2): ported
directly from the notebook's cell 13 / build_reservoir_table; PyProBE-cross-checked.
"""
import os

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import differential_evolution
from scipy.signal import savgol_filter

BIG_COST = 1e6

# ---------------------------------------------------------------------------
# Electrode OCP functions + blend-anode model -- inlined from the retired dma_ocp_fit.py
# ("Method B"), unchanged. See that file's own git history for the full revision story
# (v1-v5 bug fixes to the clipping epsilons and the blend mixing rule); kept here only as
# validated, working code, not re-derived.
# ---------------------------------------------------------------------------
_GRAPHITE_CSV = os.path.join(
    os.path.dirname(os.path.abspath(__import__("pybamm").__file__)),
    "input", "parameters", "lithium_ion", "data", "graphite_ocp_Enertech_Ai2020.csv",
)

# Clipping epsilons: _Z_EPS_TIGHT (silicon/NMC) is pure numerical safety -- both raw
# polynomials are smooth and finite at z=0/1 on their own, no physically-motivated margin
# needed. _Z_EPS_GRAPHITE_LO is different: the graphite CSV's first row is an explicit
# "extra point to avoid extrapolation" sentinel (z=0 -> 3.5V, not real data), with genuine
# data only starting at z=0.0005 -- this margin stays clear of that sentinel-influenced
# region (its high side, z->1, has ordinary finite data and needs no margin).
_Z_EPS_TIGHT = 1e-6
_Z_EPS_GRAPHITE_LO = 0.003


def _clip_tight(z):
    return np.clip(z, _Z_EPS_TIGHT, 1.0 - _Z_EPS_TIGHT)


def _clip_graphite(z):
    return np.clip(z, _Z_EPS_GRAPHITE_LO, 1.0 - _Z_EPS_TIGHT)


def silicon_ocp_delithiation_Mark2016(sto):
    """Verbrugge/Baker/Xiao silicon DELITHIATION-branch OCP -- the branch every RPT this
    pipeline fits actually needs (full-cell discharge -> anode delithiates), per PyBaMM's
    own "current sigmoid" hysteresis submodel (h->+1 during discharge collapses its
    (1+h)/2*U_delith + (1-h)/2*U_lith blend to U_delith exactly)."""
    p = [-51.02, 161.3, -205.7, 140.2, -58.76, 16.87, -3.792, 0.9937]
    return np.polyval(p, _clip_tight(sto))


def _load_graphite_ocp_spline():
    data = np.genfromtxt(_GRAPHITE_CSV, delimiter=",", comments="#")
    x, y = data[:, 0], data[:, 1]
    order = np.argsort(x)
    return CubicSpline(x[order], y[order])


_graphite_spline = _load_graphite_ocp_spline()


def graphite_ocp_Enertech_Ai2020(sto):
    """Ai2020/Enertech graphite OCP, cubic-spline interpolated from the same CSV
    si_gr_expansion.py uses."""
    return _graphite_spline(_clip_graphite(sto))


def nmc_LGM50_ocp_Chen2020(sto):
    """Chen2020 NMC OCP -- the raw analytic fit from si_gr_expansion.py, without that
    file's differentiability-only singular term (irrelevant for a numerical scipy fit)."""
    sto = _clip_tight(sto)
    return (
        -0.8090 * sto + 4.4875
        - 0.0428 * np.tanh(18.5138 * (sto - 0.5542))
        - 17.7326 * np.tanh(15.7890 * (sto - 0.3117))
        + 17.5842 * np.tanh(15.9308 * (sto - 0.3120))
    )


# Blend-anode OCP: Eq. 1 (Rehm et al.) is z_blend(V) = xi*z_Si(V) + (1-xi)*z_Gr(V) --
# capacity as a function of voltage (Si/Gr particles sit in parallel at the same local
# voltage; it's their capacities at a given voltage that combine, not their voltages at a
# given capacity). Implementation inverts each material's own OCP(z) onto a shared voltage
# grid once at import time (both are monotonically decreasing in z over their usable
# range, so the inversion is well-posed), then combines at any xi and inverts back.
_Z_GRID_SI = np.linspace(_Z_EPS_TIGHT, 1 - _Z_EPS_TIGHT, 4000)
_Z_GRID_GR = np.linspace(_Z_EPS_GRAPHITE_LO, 1 - _Z_EPS_TIGHT, 4000)
_V_SI_GRID = silicon_ocp_delithiation_Mark2016(_Z_GRID_SI)   # descending in V
_V_GR_GRID = graphite_ocp_Enertech_Ai2020(_Z_GRID_GR)        # descending in V

_V_GRID = np.linspace(0.001, 1.3, 4000)  # shared voltage grid, covers both materials' range
_Z_SI_OF_V = np.interp(_V_GRID, _V_SI_GRID[::-1], _Z_GRID_SI[::-1])
_Z_GR_OF_V = np.interp(_V_GRID, _V_GR_GRID[::-1], _Z_GRID_GR[::-1])


def blend_anode_v_of_z(z_target, xi):
    """Invert z_blend(V) = xi*z_Si(V) + (1-xi)*z_Gr(V) to get V(z_blend) -- the anode's
    OWN pOCV as a function of its blended stoichiometry z_target (0-1 fraction of the
    BLEND's total capacity, not either material's own)."""
    z_blend_grid = xi * _Z_SI_OF_V + (1 - xi) * _Z_GR_OF_V
    return np.interp(np.clip(z_target, 0.0, 1.0), z_blend_grid[::-1], _V_GRID[::-1])


# ---------------------------------------------------------------------------
# Cross-RPT physical constraints (paper Table 2's "max per-RPT decrease/increase" bounds)
# -- inlined from the retired dma_ocp_fit.py, unchanged.
# ---------------------------------------------------------------------------
CONSTRAINT_BOUNDS = {
    # key: (max_decrease, max_increase) -- either may be None ("not set" in the paper)
    "LAM_an": (0.01, None),
    "LAM_ca": (0.01, None),
    "LAM_Gr": (0.005, 0.05),
    "LAM_Si": (0.01, 0.05),
    "LI": (0.03, None),
}
CONSTRAINT_PENALTY_SCALE = 4.0  # in the same normalised-MSE units as ocv_term/dv_term;
                                  # large enough to dominate the cost once a bound is
                                  # breached, without needing a hard constraint
PARAM_NAMES = ("nu_ne", "nu_pe", "sigma_ne", "sigma_pe", "phi_si", "delta_v")
PARAM_NAMES_SOC_R = ("nu_ne", "nu_pe", "sigma_ne", "sigma_pe", "phi_si",
                      "delta_v0", "delta_v1", "delta_v2")

# Notebook's own bounds (theta = [nu_NE, nu_PE, sigma_NE, sigma_PE, phi_Si, delta_v]),
# unchanged -- see module docstring point 1 for why these are intentionally not
# re-bounded to [0, 1] the way dma_ocp_fit.DEFAULT_BOUNDS is.
DEFAULT_BOUNDS = [
    (0.7, 1.05),   # nu_ne
    (0.7, 1.05),   # nu_pe
    (-0.6, 0.2),   # sigma_ne
    (-0.2, 0.6),   # sigma_pe
    (0.0, 0.3),    # phi_si
    (-0.08, 0.08), # delta_v
]

LOW, HIGH = 0.02, 0.98  # notebook's own "low"/"high" ROI margin

# BOUNDS FIX: DEFAULT_BOUNDS above (the notebook's own values) are calibrated to their
# specific real Samsung cell's electrode geometry -- e.g. phi_si upper bound 0.3 assumes
# a cell with at most 30% Si capacity share, and nu_pe in [0.7, 1.05] assumes a cathode
# window spanning ~70-100% of its own stoichiometry. Neither holds for this project's own
# recipe (dma_baseline_run.py's known truth: xi=0.45 already exceeds 0.3; the cathode
# window spans x_pe_lo-x_pe_hi=0.7618-0.2660=0.496, well under 0.7) -- confirmed by a
# smoke-test fit pinning at exactly nu_pe=0.700 and phi_si=0.300 (both bounds) with a
# 63 mV RMSE. Adapted from dma_ocp_fit.DEFAULT_BOUNDS's own (x_lo, x_hi, xi) ranges via
# nu_e=x_hi-x_lo, sigma_ne=x_ne_lo, sigma_pe=1-x_pe_lo. delta_v is widened well past the
# notebook's +-0.08 V since dma_ocp_fit.py's own R fit reached up to ~550 mOhm at 0.25 A
# late-life (~140 mV) -- this project's own resistance growth is real and must fit inside
# the bound, not be a real-cell-specific assumption.
OUR_RECIPE_BOUNDS = [
    (0.05, 1.0),   # nu_ne   = x_ne_hi - x_ne_lo
    (0.05, 1.0),   # nu_pe   = x_pe_lo - x_pe_hi
    (0.0, 1.0),    # sigma_ne = x_ne_lo
    (0.0, 1.0),    # sigma_pe = 1 - x_pe_lo
    (0.02, 0.98),  # phi_si  = xi (Si capacity share)
    (-0.25, 0.05), # delta_v (V) -- flat offset, algebraically -I_rpt*R
]

# SOC-DEPENDENT RESISTANCE variant (8 params instead of 6) -- see reconstruct_pocv_v2's
# delta_v_params docstring. delta_v0's own range matches OUR_RECIPE_BOUNDS above;
# delta_v1 (linear-in-centred-SOC coefficient) and delta_v2 (quadratic) are each allowed
# +-0.3 V, which alone could contribute up to +-0.15 V at the SOC extremes -- ample
# headroom for the ~43 mV peak-to-peak SOC-dependent residual found with a flat delta_v,
# without being so wide it can freely fake unrelated shape.
OUR_RECIPE_BOUNDS_SOC_R = OUR_RECIPE_BOUNDS[:5] + [
    (-0.25, 0.05), # delta_v0 (V)
    (-0.3, 0.3),   # delta_v1 (V, linear-in-(soc-0.5) coefficient)
    (-0.3, 0.3),   # delta_v2 (V, quadratic-in-(soc-0.5) coefficient)
]


# ---------------------------------------------------------------------------
# Reconstruction (notebook's compute_model_components, adapted to reuse this project's
# own validated OCP/blend functions instead of the notebook's real-cell interpolators)
# ---------------------------------------------------------------------------
def reconstruct_pocv_v2(soc, nu_ne, nu_pe, sigma_ne, sigma_pe, phi_si, *delta_v_params):
    """V_model(soc) = [OCP_pe(soc_pe) + delta_v(soc)] - OCP_ne_blend(soc_ne, phi_si), with
    soc_ne = nu_ne*soc + sigma_ne, soc_pe = 1-(nu_pe*soc+sigma_pe) -- notebook's own
    coordinate transform (cell 7), unchanged. soc=1 is fully charged (4.2V), soc=0 is
    fully discharged (2.5V), matching this project's own z_full convention (confirmed
    consistent with the notebook's own "discharge" direction handling, cell 8).

    delta_v_params : 1 value -> FLAT offset delta_v0 (the notebook's own convention,
    algebraically -I_rpt*R since current is constant per RPT). 3 values -> SOC-DEPENDENT
    quadratic delta_v(soc) = delta_v0 + delta_v1*(soc-0.5) + delta_v2*(soc-0.5)**2.
    Auto-detected from how many trailing params are passed (via *theta), so existing
    6-param theta/bounds keep working unchanged -- see OUR_RECIPE_BOUNDS_SOC_R for the
    8-param version's bounds.

    Added because dvdq_rmse_constrained's own residual (measured-reconstructed) with a
    flat delta_v shows a strong, structured SOC-dependent pattern (checked directly on
    RPT14: +16 to +23 mV for soc<0.32, swinging to -20 mV around soc=0.4-0.5, recovering
    toward 0 by soc=0.86-0.95) -- a ~43 mV peak-to-peak swing a single scalar cannot
    represent at all, left as pure residual error. Since nothing else in the model can
    absorb an SOC-varying discrepancy except the cathode window, this is a strong
    candidate for what forces nu_pe/sigma_pe to distort in compensation (showing up as
    inflated LAM_pos)."""
    soc = np.asarray(soc, dtype=float)
    soc_ne = np.clip(nu_ne * soc + sigma_ne, 0.0, 1.0)
    soc_pe = np.clip(1.0 - (nu_pe * soc + sigma_pe), 0.0, 1.0)
    if len(delta_v_params) == 1:
        delta_v = delta_v_params[0]
    elif len(delta_v_params) == 3:
        d0, d1, d2 = delta_v_params
        s = soc - 0.5
        delta_v = d0 + d1 * s + d2 * s ** 2
    else:
        raise ValueError(f"expected 1 (flat) or 3 (SOC-dependent) delta_v params, got {len(delta_v_params)}")
    U_pe = nmc_LGM50_ocp_Chen2020(soc_pe) + delta_v
    U_ne = blend_anode_v_of_z(soc_ne, phi_si)
    return U_pe - U_ne


def is_theta_valid(theta, soc_exp):
    """Notebook's own feasibility check (cell 10): both electrodes' stoichiometry must
    stay within [0, 1] across the measured SOC range -- unchanged."""
    nu_ne, nu_pe, sigma_ne, sigma_pe = theta[0], theta[1], theta[2], theta[3]
    z_min, z_max = np.min(soc_exp), np.max(soc_exp)
    soc_ne_min = nu_ne * z_min + sigma_ne
    soc_ne_max = nu_ne * z_max + sigma_ne
    soc_pe_min = 1 - (nu_pe * z_max + sigma_pe)
    soc_pe_max = 1 - (nu_pe * z_min + sigma_pe)
    return (
        0 <= soc_ne_min <= 1 and 0 <= soc_ne_max <= 1 and
        0 <= soc_pe_min <= 1 and 0 <= soc_pe_max <= 1
    )


# ---------------------------------------------------------------------------
# Degradation modes -- ported from the notebook's compute_reservoirs/build_reservoir_
# table; see module docstring for the algebraic-equivalence check against dma_ocp_fit.py.
# ---------------------------------------------------------------------------
def compute_reservoirs(theta, q_full):
    nu_ne, nu_pe, sigma_ne, sigma_pe = theta[0], theta[1], theta[2], theta[3]
    q_ne = q_full / nu_ne
    q_pe = q_full / nu_pe
    q_li = q_full * (sigma_ne / nu_ne + (1 - sigma_pe) / nu_pe)
    return q_ne, q_pe, q_li


def compute_degradation_modes_v2(theta_ini, theta_now, q_full_ini, q_full_now,
                                  li_formula="pyprobe"):
    """Same dict shape/keys as dma_ocp_fit.compute_degradation_modes, for drop-in use
    with dma_common's plotting helpers.

    li_formula : "pyprobe" (default) uses q_li from compute_reservoirs (PyProBE's own
    total-lithium-at-the-SOC=0-reference definition, Qli = Q_pe*x_pe_lo + Q_ne*x_ne_lo --
    see that function's own history for why this replaced the paper's Eq. 8 for Method
    C). "rehm" uses the ORIGINAL Rehm et al. 2026 Eq. 8 instead, re-derived here in v2's
    own (nu, sigma) parameterisation for a direct side-by-side comparison:
        Eq. 8's own (alpha, beta) come from x_e(x_full)=(x_full-beta_e)/alpha_e (Eqs.
        4/5), giving alpha_ne=1/nu_ne, beta_ne=-sigma_ne/nu_ne, alpha_pe=-1/nu_pe,
        beta_pe=(1-sigma_pe)/nu_pe (note alpha_pe is NEGATIVE -- the cathode's opposite
        lithiation direction falls out of the algebra, same as Method B's own finding).
        Substituting into Eq. 8's own (alpha_pe+beta_pe-beta_ne)*Q_full term simplifies
        to exactly: A = -sigma_pe*Q_pe + sigma_ne*Q_ne, LI = 1 - A_now/A_ini. Contrast
        with PyProBE's Qli = (1-sigma_pe)*Q_pe + sigma_ne*Q_ne = Q_pe + A -- i.e. the two
        formulas differ by exactly Q_pe (the cathode's own total capacity), not a
        constant factor, since Q_pe itself also changes between RPTs via LAM_ca."""
    phi_i = theta_ini[4]
    phi = theta_now[4]
    q_ne_i, q_pe_i, q_li_i = compute_reservoirs(theta_ini, q_full_ini)
    q_ne, q_pe, q_li = compute_reservoirs(theta_now, q_full_now)

    lam_an = 1 - q_ne / q_ne_i
    lam_ca = 1 - q_pe / q_pe_i
    if li_formula == "pyprobe":
        li = 1 - q_li / q_li_i
    elif li_formula == "rehm":
        sigma_pe_i, sigma_ne_i = theta_ini[3], theta_ini[2]
        sigma_pe, sigma_ne = theta_now[3], theta_now[2]
        a_i = -sigma_pe_i * q_pe_i + sigma_ne_i * q_ne_i
        a_now = -sigma_pe * q_pe + sigma_ne * q_ne
        li = 1 - a_now / a_i
    else:
        raise ValueError(f"li_formula must be 'pyprobe' or 'rehm', got {li_formula!r}")
    lam_gr = 1 - ((1 - phi) * q_ne) / ((1 - phi_i) * q_ne_i)
    lam_si = 1 - (phi * q_ne) / (phi_i * q_ne_i)
    return dict(LAM_an=lam_an, LAM_ca=lam_ca, LI=li, LAM_Gr=lam_gr, LAM_Si=lam_si)


# ---------------------------------------------------------------------------
# Cross-RPT physical constraint -- OPTIONAL (toggle via fit_rpt_v2's constraint_args),
# so Method C can be run WITH or WITHOUT it and compared directly. Not part of the
# notebook's own optimise_rpts (which only warm-starts from the previous RPT's theta,
# with no penalty on how far the resulting degradation modes move) -- added here because
# without it, Method C's per-RPT-good-but-cross-RPT-inconsistent fits let LLI/LAM_pos
# swing non-monotonically late-life (confirmed: LLI climbed to ~1.2 A.h by RPT 12 then
# reversed down to 0.69 at RPT 14 -- physically impossible for a cumulative-loss metric).
# Reuses dma_ocp_fit.py's own CONSTRAINT_BOUNDS/CONSTRAINT_PENALTY_SCALE verbatim (same
# DM dict keys/semantics, since compute_degradation_modes_v2 returns the identical shape)
# rather than redefining them, so both methods are held to the same physical bounds.
# ---------------------------------------------------------------------------
def _constraint_penalty_v2(theta_now, theta_prev, theta_ini, q_full_prev, q_full_now,
                            q_full_ini, bounds=CONSTRAINT_BOUNDS, phi_si_bounds=None,
                            li_formula="pyprobe"):
    """v2 analogue of dma_ocp_fit._constraint_penalty -- same soft per-mode
    [-max_decrease, +max_increase] penalty between consecutive RPTs, computed from
    compute_degradation_modes_v2 instead of dma_ocp_fit.compute_degradation_modes.

    phi_si_bounds : optional (max_decrease, max_increase) applied to phi_si's own
    PER-RPT-STEP movement (theta index 4), on top of the LAM_Gr/LAM_Si bounds above.
    BUG FOUND testing this alone: it bounds the STEP size between consecutive RPTs, but
    a small per-step drift that never breaches the bound can still accumulate, and worse,
    even a TINY overshoot above phi_si's initial (RPT0) value is enough to flip the
    Si/Gr ordering (see _phi_si_vs_ini_penalty below, applied separately in
    dvdq_rmse_constrained with its own much larger scale) -- confirmed empirically: with
    only this per-step bound (0.01), phi_si sat at 0.451-0.457 for RPTs 1-8 (each
    individual step tiny and within bound) while the true reference is 0.450, giving
    LAM_Gr > LAM_Si
    (backwards) for 8 of the first 9 RPTs despite "smooth" phi_si movement.

    phi_si_vs_ini_tol : the ACTUAL correct constraint, derived exactly (not empirically):
    LAM_Si - LAM_Gr = [(1-phi)/(1-phi_i) - phi/phi_i] * (Q_ne/Q_ne_i). Since Q_ne/Q_ne_i
    is always positive, the SIGN of (LAM_Si-LAM_Gr) depends ONLY on whether
    phi_si(now) is below or above phi_si(RPT0/theta_ini) -- entirely independent of
    nu_ne. So the per-step phi_si_bounds above cannot guarantee correct ordering (it
    only bounds the rate of drift, not the absolute position relative to the anchor);
    this bounds phi_si(now) directly against phi_si(theta_ini), which is what the
    ordering actually depends on. Since this recipe's Si degrades faster than Gr at
    EVERY RPT (confirmed against Method A's own ground truth), phi_si should never
    exceed its RPT0 value by more than this tolerance."""
    dm_prev = compute_degradation_modes_v2(theta_ini, theta_prev, q_full_ini, q_full_prev,
                                            li_formula=li_formula)
    dm_now = compute_degradation_modes_v2(theta_ini, theta_now, q_full_ini, q_full_now,
                                           li_formula=li_formula)
    penalty = 0.0
    for key in dm_prev:
        max_decrease, max_increase = bounds[key]
        delta = dm_now[key] - dm_prev[key]
        if max_decrease is not None and delta < -max_decrease:
            penalty += (delta + max_decrease) ** 2
        elif max_increase is not None and delta > max_increase:
            penalty += (delta - max_increase) ** 2
    if phi_si_bounds is not None:
        max_decrease, max_increase = phi_si_bounds
        delta_phi = theta_now[4] - theta_prev[4]
        if max_decrease is not None and delta_phi < -max_decrease:
            penalty += (delta_phi + max_decrease) ** 2
        elif max_increase is not None and delta_phi > max_increase:
            penalty += (delta_phi - max_increase) ** 2
    return penalty


PHI_SI_VS_INI_PENALTY_SCALE = 5000.0  # matching soft_lambda's order of magnitude (the
# vbound_pen convention already in this module) rather than CONSTRAINT_PENALTY_SCALE
# (4.0, tuned for the LAM-based bounds above) -- this constraint is an exact algebraic
# identity (see _constraint_penalty_v2's phi_si_vs_ini_tol docstring), not a heuristic,
# so it should behave close to a hard constraint: any violation is simply wrong, and a
# weak penalty (tried first, at CONSTRAINT_PENALTY_SCALE) let several RPTs' phi_si sit
# 1-7 thousandths above the anchor anyway, which was enough to flip the ordering.


def _phi_si_vs_ini_penalty(theta_now, theta_ini, tol):
    """Standalone penalty for phi_si(now) exceeding phi_si(ini) by more than `tol` --
    kept separate from _constraint_penalty_v2 (and its DM-dependent, more lightly
    scaled bounds) so it can use its own much larger PHI_SI_VS_INI_PENALTY_SCALE."""
    overshoot = theta_now[4] - theta_ini[4] - tol
    return overshoot ** 2 if overshoot > 0 else 0.0


LAM_AN_STEP_PENALTY_SCALE = 2000.0  # same reasoning as PHI_SI_VS_INI_PENALTY_SCALE above:
# V2_CONSTRAINT_BOUNDS["LAM_an"]'s max-increase cap, run through _constraint_penalty_v2
# at the shared CONSTRAINT_PENALTY_SCALE=4.0, was confirmed NOT binding -- RPT14's fitted
# step still overshot a 0.06 cap by 0.0164 (checked directly via compute_degradation_
# modes_v2 on the fitted thetas), because (0.0164**2)*4.0 = 0.0011 is negligible next to
# the voltage-fit terms' own scale. Rather than re-tuning CONSTRAINT_PENALTY_SCALE (shared
# by every other LAM_Gr/LAM_Si/LAM_ca/LI bound, all already tuned against it), give only
# this one an independent, much larger scale -- exactly the fix already applied once
# before for phi_si_vs_ini_tol when the same shared-scale-too-weak failure mode showed up.


def _lam_an_step_penalty(theta_now, theta_prev, theta_ini, q_full_prev, q_full_now,
                          q_full_ini, max_increase, li_formula="pyprobe"):
    """Standalone, heavily-scaled penalty for LAM_an's own RPT-to-RPT step exceeding
    `max_increase` -- see LAM_AN_STEP_PENALTY_SCALE above for why this needs to be
    separate from _constraint_penalty_v2's shared-scale LAM_an bound."""
    dm_prev = compute_degradation_modes_v2(theta_ini, theta_prev, q_full_ini, q_full_prev,
                                            li_formula=li_formula)
    dm_now = compute_degradation_modes_v2(theta_ini, theta_now, q_full_ini, q_full_now,
                                           li_formula=li_formula)
    overshoot = (dm_now["LAM_an"] - dm_prev["LAM_an"]) - max_increase
    return overshoot ** 2 if overshoot > 0 else 0.0


# ---------------------------------------------------------------------------
# Cost-function machinery -- ported from the notebook's cell 11, near-verbatim (only
# renamed for module-local clarity and to take `direction`/`low`/`high` as explicit
# arguments instead of notebook globals).
# ---------------------------------------------------------------------------
def _safe_grad(x, eps=1e-12):
    g = np.gradient(np.asarray(x, dtype=float))
    g[np.abs(g) < eps] = eps
    return g


def _safe_savgol(y, window=11, poly=2, deriv=0, delta=1.0):
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 3:
        return y.copy() if deriv == 0 else np.zeros_like(y)
    w = int(window) | 1
    w = min(w, n if n % 2 == 1 else n - 1)
    min_valid_w = poly + 2
    if min_valid_w % 2 == 0:
        min_valid_w += 1
    w = max(w, min_valid_w)
    if w > n:
        w = n if n % 2 == 1 else n - 1
    if w <= poly or w < 3:
        return y.copy() if deriv == 0 else np.gradient(y)
    return savgol_filter(y, window_length=w, polyorder=poly, deriv=deriv, delta=delta)


def _clamp(x, lo, hi):
    return max(lo, min(hi, x))


def _rmse(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return float(np.sqrt(np.mean((a - b) ** 2)))


def _huber(a, b, delta=0.02):
    r = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    ar = np.abs(r)
    return float(np.sqrt(np.mean(np.where(ar <= delta, 0.5 * r * r, delta * (ar - 0.5 * delta)))))


# NOTE: no longer called anywhere in this module (see dvdq_rmse_constrained's own
# comment for why the voltage-based re-interpolation it enabled was a bug) -- kept
# defined for reference/fidelity to the notebook it was ported from.
def _interp_sorted_unique(x_src, y_src, x_tgt):
    x_src, y_src, x_tgt = (np.asarray(a, dtype=float) for a in (x_src, y_src, x_tgt))
    order = np.argsort(x_src)
    x_sorted, y_sorted = x_src[order], y_src[order]
    x_unique, inv = np.unique(x_sorted, return_inverse=True)
    y_accum = np.zeros_like(x_unique, dtype=float)
    counts = np.zeros_like(x_unique, dtype=float)
    np.add.at(y_accum, inv, y_sorted)
    np.add.at(counts, inv, 1.0)
    y_unique = y_accum / np.maximum(counts, 1.0)
    if len(x_unique) < 2:
        return np.full_like(x_tgt, np.nan, dtype=float)
    return np.interp(x_tgt, x_unique, y_unique)


def _compute_delta(soc, V, direction="discharge", lo=LOW, hi=HIGH, smooth=True,
                    sg_win=11, sg_poly=2):
    """SOC at the sharpest dV/dQ feature (argmin of the SIGNED derivative) within
    [lo, hi] -- this is where the Si/Gr voltage-plateau crossover notch sits. Splits the
    domain for _window_masks below."""
    soc, V = np.asarray(soc, dtype=float), np.asarray(V, dtype=float)
    if smooth:
        V = _safe_savgol(V, window=sg_win, poly=sg_poly)
    dVdQ = _safe_grad(V) / _safe_grad(soc)
    s = -1.0 if str(direction).lower() == "discharge" else 1.0
    metric = s * dVdQ
    mask = (soc >= lo) & (soc <= hi) & np.isfinite(metric)
    if not np.any(mask):
        raise ValueError("No valid points to compute delta.")
    soc_win, metric_win = soc[mask], metric[mask]
    idx = np.nanargmin(metric_win)
    return _clamp(float(soc_win[idx]), lo + 0.02, hi - 0.02)


def _window_masks(soc, delta, direction="discharge", lo=LOW, hi=HIGH):
    """discharge: V-fit window = [lo, delta], dV/dQ-fit window = [delta, hi] (opposite
    for charge) -- i.e. fit raw voltage on the smooth side of the crossover notch, and
    the (Huber-weighted, smoothed) derivative on the side containing it."""
    if str(direction).lower() == "discharge":
        a_lo, a_hi, b_lo, b_hi = lo, delta, delta, hi
    else:
        a_lo, a_hi, b_lo, b_hi = delta, hi, lo, delta
    mask_v = (soc >= min(a_lo, a_hi)) & (soc <= max(a_lo, a_hi))
    mask_dvdq = (soc >= min(b_lo, b_hi)) & (soc <= max(b_lo, b_hi))
    return mask_v, mask_dvdq


def dvdq_rmse_constrained(theta, soc, v_meas, direction="discharge", delta=None,
                           use_span_check=True, vspan_tol=0.03, soft_lambda=5000.0,
                           smooth=True, sg_win=11, sg_poly=2, min_pts=6,
                           theta_valid_fn=is_theta_valid, low=LOW, high=HIGH,
                           constraint_args=None, constraint_bounds=CONSTRAINT_BOUNDS,
                           w_full_dvdq=0.0, phi_si_bounds=None, phi_si_vs_ini_tol=None,
                           w_ocv_b=0.0, w_dv_b=0.0, ocv_scale_b=1.0, dv_scale_b=1.0,
                           w_v_split=0.2, w_dvdq_split=0.8, dvdq_huber_delta=0.15,
                           li_formula="pyprobe", lam_an_step_cap=None,
                           full_dvdq_smooth_window=15, full_dvdq_smooth_poly=3):
    """Notebook's own dvdq_rmse_constrained (cell 11), ported: cost = w_v_split*RMSE(V) +
    w_dvdq_split*Huber(dV/dQ) + soft V-endpoint penalty, computed over the delta-split
    windows above (w_v_split/w_dvdq_split default to the notebook's own 0.2/0.8; set
    either to 0 to drop that delta-split term entirely -- e.g. for a clean whole-domain
    OCV+DV+dV/dQ cost, zero both of these and use w_ocv_b/w_dv_b/w_full_dvdq instead,
    matching the paper's own Eq. 3 three-term structure with all three actually turned on
    rather than its own lambda_IC=0 default). Returns BIG_COST for any infeasible/
    degenerate theta (out-of-bounds stoichiometry, too few points in a window, non-finite
    reconstruction).

    constraint_args : optional (theta_prev, theta_ini, q_full_prev, q_full_now,
    q_full_ini) tuple -- see _constraint_penalty_v2. None (the default, and the
    notebook's own behaviour) disables the cross-RPT constraint entirely; pass it to
    enable dma_ocp_fit's soft per-mode bound between consecutive RPTs, so Method C can be
    run with or without it and compared directly.

    w_full_dvdq : optional EXTRA dV/dQ term over the WHOLE [low, high] ROI (Huber loss,
    same analytic-savgol-derivative machinery as the delta-split term). OFF by default
    (0.0, matching the notebook exactly) -- added because the notebook's own delta-split
    design only ever scores dV/dQ agreement on ONE side of the crossover point (mask_dvdq
    -- [delta, hi] for discharge); the OTHER side (mask_v) is scored on raw voltage only,
    with no derivative term at all, which left dma_method_c_dvdq_fits.png visibly poor
    across roughly half the domain. This term does not replace the delta-split terms
    (which still do the heavy lifting where the sharp crossover notch sits) -- it just
    ensures the WHOLE curve's derivative shape is actually optimized, not just one side.

    w_ocv_b, w_dv_b, ocv_scale_b, dv_scale_b : "bring in Method B's cost" option. Adds
    Method B's OWN additive, whole-[low,high]-domain, MSE-normalised OCV term
    (mean(((V_model-V_meas)/ocv_scale_b)**2)) and/or raw-finite-difference dV/dSOC term
    (same normalisation, using np.gradient -- NOT the savgol-smoothed dV/dQ above, to be
    a faithful port of dma_ocp_fit._dV_dSOC exactly) as EXTRA weighted terms alongside
    Method C's own delta-split cost. OFF by default (0.0). Motivation: Method B's own
    per-RPT trace keeps xi anchored in a narrow band (0.431-0.450) and gets the
    LAM_Si>LAM_Gr ordering right at every single RPT, which this project's own analysis
    traced to that whole-domain, voltage-heavy cost giving dense, uniform signal that
    pins phi_si/xi tightly -- vs. Method C's own derivative-heavy, delta-split cost,
    which left phi_si more free to wander (see phi_si_bounds' docstring for the direct
    fix). This lets B's anchoring effect be tested as a cost-function ingredient in its
    own right, instead of only via the (more ad hoc) direct phi_si_bounds constraint.
    ocv_scale_b/dv_scale_b should be precomputed ONCE per RPT from v_meas (see
    fit_rpt_v2), matching dma_ocp_fit.fit_rpt's own ocv_scale/dv_scale convention."""
    if not theta_valid_fn(theta, soc):
        return BIG_COST
    try:
        v_model = reconstruct_pocv_v2(soc, *theta)
    except Exception:
        return BIG_COST
    if v_model is None or np.any(~np.isfinite(v_model)):
        return BIG_COST

    v_e = _safe_savgol(v_meas, window=sg_win, poly=sg_poly) if smooth else v_meas
    w = max(7, int(sg_win) | 1)  # shared Savitzky-Golay derivative window, used below

    vbound_pen = 0.0
    if use_span_check:
        err_min = abs(np.nanmin(v_model) - np.nanmin(v_meas))
        err_max = abs(np.nanmax(v_model) - np.nanmax(v_meas))
        viol_min = max(0.0, err_min - vspan_tol)
        viol_max = max(0.0, err_max - vspan_tol)
        vbound_pen = soft_lambda * (viol_min ** 2 + viol_max ** 2)

    rmse_volt, rmse_dvdq = 0.0, 0.0
    if w_v_split or w_dvdq_split:
        if delta is None:
            try:
                delta = _compute_delta(soc, v_e, direction, low, high, smooth=False)
            except Exception:
                return BIG_COST

        mask_v, mask_dvdq = _window_masks(soc, delta, direction, low, high)
        if (np.sum(mask_v) < min_pts) or (np.sum(mask_dvdq) < min_pts):
            pad = 0.03
            delta = min(delta, high - pad) if str(direction).lower() == "charge" else max(delta, low + pad)
            mask_v, mask_dvdq = _window_masks(soc, delta, direction, low, high)
            if (np.sum(mask_v) < min_pts) or (np.sum(mask_dvdq) < min_pts):
                mid = 0.5
                if str(direction).lower() == "charge":
                    mask_v = (soc >= mid) & (soc <= high)
                    mask_dvdq = (soc >= low) & (soc < mid)
                else:
                    mask_v = (soc >= low) & (soc <= mid)
                    mask_dvdq = (soc > mid) & (soc <= high)
                if (np.sum(mask_v) < min_pts) or (np.sum(mask_dvdq) < min_pts):
                    return BIG_COST

        try:
            rmse_volt = _rmse(v_e[mask_v], v_model[mask_v])
        except Exception:
            return BIG_COST

        try:
            v_exp_b, v_mod_b, soc_b = v_e[mask_dvdq], v_model[mask_dvdq], soc[mask_dvdq]
            if len(v_exp_b) < min_pts:
                return BIG_COST
            dV_model_b = _safe_savgol(v_mod_b, window=w, poly=sg_poly, deriv=1)
            dV_exp_b = _safe_savgol(v_exp_b, window=w, poly=sg_poly, deriv=1)
            dQ_b = _safe_grad(soc_b)
            dvdq_model_b = dV_model_b / dQ_b
            dvdq_exp_b = dV_exp_b / dQ_b
            # BUG FIX: dvdq_model_b and dvdq_exp_b are already on the IDENTICAL soc_b
            # grid (same index = same SOC position) -- comparing them directly is a
            # proper SOC-aligned comparison. The notebook's own code instead
            # re-interpolated dvdq_model_b onto v_mod_b (the model's own voltage) and
            # evaluated that at v_exp_b (the MEASURED voltage) -- i.e. it compared the
            # two curves at matched VOLTAGE, not matched SOC. That only agrees with a
            # SOC-matched comparison when the underlying voltage fit has ~zero residual;
            # once there's a genuine voltage discrepancy (confirmed present and growing
            # at later RPTs -- see the residual-vs-SOC diagnostic), evaluating the
            # model's dV/dQ "at the voltage the measurement reaches" pulls the comparison
            # to a DIFFERENT SOC location than intended, which is exactly the "shape is
            # there but squeezed, peaks don't align" symptom seen from RPT11 onward.
            # (This project's own full_dvdq_term, just below, already compares directly
            # by SOC index and never had this bug.)
            if np.any(~np.isfinite(dvdq_model_b)):
                return BIG_COST
            rmse_dvdq = _huber(dvdq_exp_b, dvdq_model_b, delta=0.02)
        except Exception:
            return BIG_COST

    full_dvdq_term = 0.0
    if w_full_dvdq:
        try:
            roi_mask = (soc >= low) & (soc <= high)
            soc_full = soc[roi_mask]
            v_e_full, v_mod_full = v_e[roi_mask], v_model[roi_mask]
            dV_e_full = _safe_savgol(v_e_full, window=w, poly=sg_poly, deriv=1)
            dV_mod_full = _safe_savgol(v_mod_full, window=w, poly=sg_poly, deriv=1)
            dQ_full = _safe_grad(soc_full)
            dvdq_e_full = dV_e_full / dQ_full
            dvdq_mod_full = dV_mod_full / dQ_full
            # HUBER DELTA FIX: the notebook's own delta=0.02 is calibrated to THEIR real
            # cell's dV/dQ scale. This recipe's own dV/dQ ranges roughly -1.5 to -0.15
            # (std~0.25) -- 0.02 is only ~8% of that std, so almost every meaningful
            # residual (peak-misalignment error included) already falls in Huber's
            # LINEAR/downweighted regime instead of being penalized quadratically like a
            # plain MSE would.
            #
            # RPT14 INVESTIGATION (tried derivative-method + Huber-delta + plain-MSE swaps,
            # settled on a DIFFERENT fix -- see fit_rpt_v2's "global" branch): confirmed
            # directly that v2's DERIVATIVE-based cost terms (this one and dv_term) score
            # the MISALIGNED RPT14 solution as better than the correctly-notched one EVEN
            # WITH plain RMSE and Method B's own smooth-then-gradient derivative -- the two
            # objectives (raw voltage fit vs. v2's own weighted derivative cost) genuinely
            # disagree at RPT14, not just because of a loss-shape/derivative-method quirk.
            # Reverted this term to its original form (no benefit once the real fix is in
            # fit_rpt_v2's own ensemble selection criterion, which now picks by raw voltage
            # RMSE instead of total cost -- see there for why).
            full_dvdq_term = _huber(dvdq_e_full, dvdq_mod_full, delta=dvdq_huber_delta)
        except Exception:
            return BIG_COST

    b_style_term = 0.0
    if w_ocv_b or w_dv_b:
        try:
            roi_mask = (soc >= low) & (soc <= high)
            soc_full = soc[roi_mask]
            v_e_full, v_mod_full = v_meas[roi_mask], v_model[roi_mask]
            if w_ocv_b:
                b_style_term += w_ocv_b * np.mean(((v_mod_full - v_e_full) / ocv_scale_b) ** 2)
            if w_dv_b:
                dv_e_full = np.gradient(v_e_full, soc_full)
                dv_mod_full = np.gradient(v_mod_full, soc_full)
                b_style_term += w_dv_b * np.mean(((dv_mod_full - dv_e_full) / dv_scale_b) ** 2)
        except Exception:
            return BIG_COST

    cost = (w_v_split * rmse_volt + w_dvdq_split * rmse_dvdq + w_full_dvdq * full_dvdq_term
            + b_style_term + vbound_pen)
    if constraint_args is not None:
        theta_prev, theta_ini, q_full_prev, q_full_now, q_full_ini = constraint_args
        cost += CONSTRAINT_PENALTY_SCALE * _constraint_penalty_v2(
            theta, theta_prev, theta_ini, q_full_prev, q_full_now, q_full_ini,
            bounds=constraint_bounds, phi_si_bounds=phi_si_bounds, li_formula=li_formula)
        if phi_si_vs_ini_tol is not None:
            cost += PHI_SI_VS_INI_PENALTY_SCALE * _phi_si_vs_ini_penalty(
                theta, theta_ini, phi_si_vs_ini_tol)
        if lam_an_step_cap is not None:
            cost += LAM_AN_STEP_PENALTY_SCALE * _lam_an_step_penalty(
                theta, theta_prev, theta_ini, q_full_prev, q_full_now, q_full_ini,
                lam_an_step_cap, li_formula=li_formula)
    return cost


def rmse_of_fit_v2(soc, v_meas, theta):
    resid = reconstruct_pocv_v2(soc, *theta) - v_meas
    return float(np.sqrt(np.mean(resid ** 2)))


# ---------------------------------------------------------------------------
# Fitting driver -- notebook's find_feasible_theta / make_init_population /
# optimise_one_rpt, condensed into one function with the same shrinking-restart
# structure. maxiter/max_restarts scaled down from the notebook's own (800/5, tuned for
# noisy real experimental data) since this project's simulated RPT curves are much
# cleaner -- override if a closer replication is wanted.
# ---------------------------------------------------------------------------
def _random_theta(bounds, rng):
    return np.array([rng.uniform(lo, hi) for lo, hi in bounds], dtype=float)


def find_feasible_theta(soc, v_meas, bounds, direction="discharge", trials=300, rng_seed=0):
    rng = np.random.default_rng(rng_seed)
    best_theta, best_cost = None, np.inf
    for _ in range(trials):
        theta = _random_theta(bounds, rng)
        cost = dvdq_rmse_constrained(theta, soc, v_meas, direction=direction,
                                      use_span_check=False, smooth=True)
        if np.isfinite(cost) and cost < BIG_COST and cost < best_cost:
            best_theta, best_cost = theta, cost
    return best_theta, best_cost


def make_init_population(bounds, center, popsize=60, spread=0.10, rng_seed=1):
    rng = np.random.default_rng(rng_seed)
    lo = np.array([b[0] for b in bounds], float)
    hi = np.array([b[1] for b in bounds], float)
    center = np.asarray(center, float)
    scale = (hi - lo) * spread
    pop = np.clip(center + rng.normal(0.0, scale, size=(popsize, len(bounds))), lo, hi)
    pop[0] = np.clip(center, lo, hi)
    return pop


def fit_rpt_v2(soc, v_meas, bounds=DEFAULT_BOUNDS, direction="discharge",
               init_center=None, seed=0, target_rmse=0.005, max_restarts=2,
               maxiter=300, popsize=15, feasibility_trials=300,
               constraint_args=None, constraint_bounds=CONSTRAINT_BOUNDS,
               w_full_dvdq=0.0, phi_si_bounds=None, phi_si_vs_ini_tol=None,
               w_ocv_b=0.0, w_dv_b=0.0, low=LOW, high=HIGH,
               w_v_split=0.2, w_dvdq_split=0.8, dvdq_huber_delta=0.15,
               li_formula="pyprobe", lam_an_step_cap=None,
               full_dvdq_smooth_window=15, full_dvdq_smooth_poly=3,
               search_style="restart"):
    """Global fit of (nu_ne, nu_pe, sigma_ne, sigma_pe, phi_si, delta_v) to one RPT's
    (soc, v_meas) curve via the notebook's delta-split dV/dQ-weighted cost, with an
    adaptive shrinking-spread multi-restart differential_evolution loop (stops early
    once full-curve RMSE <= target_rmse). Returns (best_theta, rmse).

    search_style : "restart" (default, this function's own original notebook-derived
    scheme) vs "global" (dma_ocp_fit.fit_rpt's own scheme -- ONE differential_evolution
    call over a full-bounds Latin Hypercube population with init_center injected as just
    one member, mutation=(0.4,1.5), recombination=0.8, no restarts/narrowing).

    FOUND (RPT14 investigation): "restart"'s population is always a GAUSSIAN centered on
    the current best with spread SHRINKING 0.7x every restart, starting at only 8% of the
    bounds range -- more restarts/iterations only refine WITHIN whatever basin restart 0
    landed in, never re-widen to escape it. Directly tested "global" on v2's own cost
    function at RPT14 (identical constraint_args/weights to the "restart" search that
    found bu6arl7hh's theta): found a solution with voltage RMSE 6.77mV (vs "restart"'s
    11.03mV -- and even better than Method B's own 10.84mV) with the dV/dQ notch visibly
    correctly positioned, essentially zero constraint-penalty violation, and a smooth,
    physically continuous LAM trajectory from RPT13 (LAM_Si +2.86pp, LAM_Gr -0.69pp, LAM_an
    +0.91pp -- no erratic jumps). This is the actual root cause the whole session's
    constraint/cost-function tuning was fighting: not a wrong penalty or a masked
    residual, but a fundamentally too-local search strategy that can never reach a
    genuinely different (better) basin once seeded near the previous RPT's fit.

    constraint_args, constraint_bounds : see dvdq_rmse_constrained -- constraint_args
    defaults to None, matching the notebook's own optimise_rpts (no cross-RPT
    constraint); pass it (fits[i-1], fits[0], q_full_prev, q_full_now, q_full_ini]) to
    enable the SAME soft cross-RPT bound Method B uses, so runs with and without it can
    be compared directly (this is the toggle -- see dma_method_c_plot.py's
    USE_CROSS_RPT_CONSTRAINT).

    w_full_dvdq : see dvdq_rmse_constrained -- OFF (0.0) by default.

    phi_si_bounds : see _constraint_penalty_v2 -- direct (max_decrease, max_increase) on
    phi_si's own per-RPT movement, only applied when constraint_args is not None.

    w_ocv_b, w_dv_b : "bring in Method B's cost" option, see dvdq_rmse_constrained.
    ocv_scale_b/dv_scale_b are computed HERE, once per RPT from v_meas (matching
    dma_ocp_fit.fit_rpt's own ocv_scale=ptp(v_meas), dv_scale=std(dV/dSOC) convention),
    not re-derived on every cost evaluation."""
    ocv_scale_b = float(np.ptp(v_meas)) or 1.0
    roi_mask_b = (soc >= low) & (soc <= high)
    dv_scale_b = float(np.std(np.gradient(v_meas[roi_mask_b], soc[roi_mask_b]))) or 1.0

    delta_fixed = _compute_delta(soc, v_meas, direction, LOW, HIGH, smooth=True)

    if init_center is None:
        theta_seed, _ = find_feasible_theta(soc, v_meas, bounds, direction,
                                             trials=feasibility_trials, rng_seed=seed)
        center = theta_seed if theta_seed is not None else np.array(
            [(lo + hi) / 2 for lo, hi in bounds])
    else:
        center = np.asarray(init_center, dtype=float).copy()

    def objective(theta):
        return dvdq_rmse_constrained(theta, soc, v_meas, direction=direction,
                                      delta=delta_fixed, use_span_check=True,
                                      vspan_tol=0.03, soft_lambda=5000.0, smooth=True,
                                      low=low, high=high,
                                      constraint_args=constraint_args,
                                      constraint_bounds=constraint_bounds,
                                      w_full_dvdq=w_full_dvdq,
                                      phi_si_bounds=phi_si_bounds,
                                      phi_si_vs_ini_tol=phi_si_vs_ini_tol,
                                      w_ocv_b=w_ocv_b, w_dv_b=w_dv_b,
                                      ocv_scale_b=ocv_scale_b, dv_scale_b=dv_scale_b,
                                      w_v_split=w_v_split, w_dvdq_split=w_dvdq_split,
                                      dvdq_huber_delta=dvdq_huber_delta,
                                      li_formula=li_formula,
                                      lam_an_step_cap=lam_an_step_cap,
                                      full_dvdq_smooth_window=full_dvdq_smooth_window,
                                      full_dvdq_smooth_poly=full_dvdq_smooth_poly)

    if search_style == "global":
        # Method B's own scheme (dma_ocp_fit.fit_rpt): ONE DE call, full-bounds Latin
        # Hypercube population (center injected as just one member, not the population's
        # basis), wider mutation/more recombination -- see this function's own docstring
        # for why this matters (never narrows, so it can reach a genuinely different basin
        # a Gaussian-around-center population never samples).
        # ENSEMBLE FIX (found AFTER a single "global" call turned out fragile): a single
        # call's result is EXTREMELY sensitive to microscopic (~3e-4) differences in the
        # seed point -- confirmed directly at RPT14: the exact cached RPT13 theta as
        # center reproducibly lands in the SAME basin "restart" search always finds
        # (11.02mV), while that same theta rounded to 3 decimal places (a ~3e-4 nudge)
        # lands in a dramatically better one (6.78mV, correctly-positioned notch). This
        # is a knife-edge basin boundary, not a robust preference of one search style --
        # so a single deterministic call cannot be trusted. Run (1+max_restarts)
        # INDEPENDENT global searches (different RNG seed each).
        #
        # SELECTION BY RMSE, NOT COST (found AFTER the ensemble above still reliably
        # returned the misaligned basin): ran all 5 seeds individually at RPT14 -- one
        # found the correctly-notched basin (RMSE 6.79mV) and four found the misaligned
        # one (RMSE 11.02mV) every time (fully deterministic per seed, not chance), but
        # the misaligned basin's total COST is LOWER (0.048 vs 0.053) than the correct
        # one's -- confirmed this holds with plain-MSE dV/dQ too, not just Huber: v2's own
        # derivative-based cost terms (dv_term, full_dvdq_term) genuinely score the WRONG
        # basin as better, because they weight derivative-SHAPE agreement more than raw
        # voltage-fit accuracy, and this particular misaligned solution happens to have a
        # marginally smoother/closer-average derivative match despite being visibly wrong.
        # Raw voltage RMSE has correctly identified the better basin in every comparison
        # done this session (matching the plotted notch alignment every time), so select
        # ensemble candidates by RMSE directly instead of by the objective's own cost --
        # skipping any candidate whose cost signals outright infeasibility (BIG_COST).
        lb = np.array([b[0] for b in bounds], dtype=float)
        ub = np.array([b[1] for b in bounds], dtype=float)
        n_params = len(bounds)
        best_x, best_rmse, best_fun = None, np.inf, np.inf
        for k in range(1 + max_restarts):
            rng = np.random.default_rng(seed * 100 + k)
            pop = lb + rng.random((popsize * n_params, n_params)) * (ub - lb)
            pop[0] = np.clip(center, lb + 1e-9, ub - 1e-9)
            result = differential_evolution(
                objective, bounds=bounds, init=pop, seed=seed * 100 + k, maxiter=maxiter,
                popsize=popsize, tol=1e-10, mutation=(0.4, 1.5), recombination=0.8,
                polish=True, updating="deferred",
            )
            if float(result.fun) >= BIG_COST:
                continue
            candidate_rmse = rmse_of_fit_v2(soc, v_meas, result.x)
            if candidate_rmse < best_rmse:
                best_x, best_rmse, best_fun = result.x.copy(), candidate_rmse, float(result.fun)
        if best_x is None:
            # every candidate was infeasible -- fall back to whatever minimises cost
            for k in range(1 + max_restarts):
                rng = np.random.default_rng(seed * 100 + k)
                pop = lb + rng.random((popsize * n_params, n_params)) * (ub - lb)
                pop[0] = np.clip(center, lb + 1e-9, ub - 1e-9)
                result = differential_evolution(
                    objective, bounds=bounds, init=pop, seed=seed * 100 + k, maxiter=maxiter,
                    popsize=popsize, tol=1e-10, mutation=(0.4, 1.5), recombination=0.8,
                    polish=True, updating="deferred",
                )
                if float(result.fun) < best_fun:
                    best_x, best_fun = result.x.copy(), float(result.fun)
        rmse = rmse_of_fit_v2(soc, v_meas, best_x)
        return best_x, rmse

    best_x, best_fun = None, np.inf
    spread = 0.08
    for k in range(1 + max_restarts):
        init_pop = make_init_population(bounds, center, popsize=max(80, popsize * 5),
                                         spread=spread, rng_seed=seed * 100 + k)
        result = differential_evolution(
            objective, bounds=bounds, init=init_pop, maxiter=maxiter, popsize=popsize,
            tol=8e-6, strategy="rand1bin", polish=True, updating="deferred", seed=seed,
        )
        if float(result.fun) < best_fun:
            best_x, best_fun = result.x.copy(), float(result.fun)
        rmse_full = rmse_of_fit_v2(soc, v_meas, best_x)
        if rmse_full <= target_rmse:
            break
        center = best_x.copy()
        spread *= 0.7

    rmse = rmse_of_fit_v2(soc, v_meas, best_x)
    return best_x, rmse
