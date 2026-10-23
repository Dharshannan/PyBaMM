"""Composite-anode degradation-mode analysis (DMA) of the C/20 RPT discharges,
fitted to voltage and dV/dQ only: no expansion term anywhere in the cost.

The electrode model is PyProBE's (ImperialCollegeLondon/PyProBE):
- the cathode is an `OCP`;
- the Si/graphite anode is a `CompositeOCP` (both components at one common
  potential, the anode stoichiometry being the capacity-weighted sum).

Everything else follows the DMA we were given (the base fit), so the only
difference is the cost:
- the same discharge half-cell OCPs: NMC622 lithiating, and graphite and Si
  (Verbrugge) delithiating;
- the same per-RPT deformation of the Si OCP, U' = scale U + shift;
- the same seven voltage parameters (Cn_Si, Cn_Gr, Cp, x_start, y_start,
  scale, shift);
- the same objective window (the last 5 % of the discharge dropped);
- the same dV/dQ weight (0.5) and LLI definition (x_start Cn + y_start Cp).

The base cost also had the cell expansion and dE/dQ in it (and two expansion
parameters); this one has

    J = mean[((V_fit - V)/0.01 V)^2] + 0.5 mean[((dVdQ_fit - dVdQ)/0.1 V/Ah)^2].
"""
import os

import numpy as np
import pandas as pd
import polars as pl
from pyprobe.analysis import degradation_mode_analysis as dma
from pyprobe.analysis import smoothing
from pyprobe.result import Result
from scipy import optimize

W_DVDQ = 0.5  # dV/dQ weight, as in the base fit
V_SCALE, D_SCALE = 0.01, 0.1  # residual scales [V], [V/Ah]: both terms O(1) at a good fit
N_POINTS = 300  # fit points per discharge, as in the base fit
PARAMS = ("Cn_Si_Ah", "Cn_Gr_Ah", "Cp_Ah", "x_start", "y_start", "si_ocp_scale", "si_ocp_shift_V")
# Bounds enclose every base value of the four cells with margin.
BOUNDS = ((0.05, 1.8), (0.9, 1.9), (2.2, 3.1), (0.6, 1.0), (0.0, 0.06), (0.75, 1.1), (-0.10, 0.03))
FOLDERS = {"CELL064": "CELL064_25C_full", "CELL017": "CELL017_45C", "CELL009": "CELL009_25C_15-95",
           "CELL026": "CELL026_25C_20-80"}


# ---------------------------------------------------------------------------
# Half-cell OCPs
# ---------------------------------------------------------------------------
def _monotone(sto, p, decreasing):
    """Sort by stoichiometry, force the potential monotone (tabulation noise), and
    average the stoichiometry over repeated potentials (the flat graphite stages),
    so that x(U) can be interpolated."""
    o = np.argsort(sto)
    sto, p = sto[o], p[o]
    p = np.minimum.accumulate(p) if decreasing else np.maximum.accumulate(p)
    d = pd.DataFrame({"sto": sto, "p": p}).groupby("p", sort=True)["sto"].mean()
    return d.to_numpy(), d.index.to_numpy()


class HalfCells:
    """The three discharge half-cell OCPs (data/ocp_halfcell)."""

    def __init__(self, data_dir):
        d = os.path.join(data_dir, "ocp_halfcell")
        nmc = pd.read_csv(os.path.join(d, "nmc622_lithiation.csv")).sort_values("sto")
        # the file's "sto" is 1 - y (y: cathode lithium fraction)
        self.cathode = dma.OCP.from_data(1 - nmc["sto"].to_numpy(), nmc["p"].to_numpy())
        gr = pd.read_csv(os.path.join(d, "graphite_delithiation.csv"))
        si = pd.read_csv(os.path.join(d, "silicon_delithiation_verbrugge.csv"))
        self.x_gr, self.u_gr = _monotone(gr["sto"].to_numpy(), gr["p"].to_numpy(), decreasing=True)
        self.x_si, self.u_si = _monotone(si["sto"].to_numpy(), si["p"].to_numpy(), decreasing=True)

    def anode(self, frac_si, scale, shift):
        """PyProBE composite anode with the deformed Si OCP; frac_si = Cn_Si / Cn."""
        comp = dma.CompositeOCP.from_data(self.x_si, scale * self.u_si + shift, self.x_gr, self.u_gr)
        comp.comp_fraction = frac_si
        return comp

    def voltage(self, theta, q):
        """Full-cell voltage at discharged capacity q [Ah] for parameters theta (PARAMS order)."""
        cn_si, cn_gr, cp, x0, y0, scale, shift = theta
        cn = cn_si + cn_gr
        x = x0 - q / cn  # composite anode stoichiometry
        y = y0 + q / cp  # cathode lithium fraction
        return self.cathode.eval(y) - self.anode(cn_si / cn, scale, shift).eval(x)


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
def rpt_files(data_dir, cell):
    d = os.path.join(data_dir, FOLDERS[cell], "low_rate_c20")
    return {int(f.split("_RPT")[1][:3]): os.path.join(d, f) for f in sorted(os.listdir(d))
            if f.endswith("_lowrate_c20_discharge.csv")}


SG_WINDOW_AH = 0.04  # Savitzky-Golay window for dV/dQ [Ah] (~1.6 % of the BoL capacity)


def load_discharge(path):
    """(q, V, dV/dQ, q_full, base_V) on N_POINTS uniform capacities inside the base
    objective window. V is the measured voltage; dV/dQ is from PyProBE's
    Savitzky-Golay filter (cubic, SG_WINDOW_AH wide) on a uniform 1500-point grid
    (a smoothing spline with an automatic penalty bends away from the data where
    the voltage falls steeply at the window end). base_V is the base fit's voltage."""
    d = pd.read_csv(path, usecols=["Q_discharge_Ah", "voltage_fit_used_V", "fit_objective_included",
                                   "model_voltage_V"])
    d = d.sort_values("Q_discharge_Ah").drop_duplicates("Q_discharge_Ah")
    q_full = float(d["Q_discharge_Ah"].max())
    w = d[d["fit_objective_included"].astype(bool)]
    qa, va = w["Q_discharge_Ah"].to_numpy(), w["voltage_fit_used_V"].to_numpy()
    q_dense = np.linspace(qa.min(), qa.max(), 1500)
    dq = q_dense[1] - q_dense[0]
    n = max(5, int(SG_WINDOW_AH / dq) // 2 * 2 + 1)
    res = Result(lf=pl.LazyFrame({"Voltage [V]": np.interp(q_dense, qa, va)}), info={})
    g = smoothing.savgol_smoothing(res, "Voltage [V]", window_length=n, polyorder=3, derivative=1).data
    q = np.linspace(qa.min(), qa.max(), N_POINTS)
    v = np.interp(q, qa, va)
    dvdq = np.interp(q, q_dense, g["Voltage [V]"].to_numpy() / dq)
    base_v = np.interp(q, w["Q_discharge_Ah"].to_numpy(), w["model_voltage_V"].to_numpy())
    return q, v, dvdq, q_full, base_v


# ---------------------------------------------------------------------------
# Fit
# ---------------------------------------------------------------------------
def residuals(hc, theta, q, v, dvdq):
    vm = hc.voltage(theta, q)
    return vm - v, np.gradient(vm, q) - dvdq


def cost(theta, hc, q, v, dvdq):
    rv, rd = residuals(hc, theta, q, v, dvdq)
    j = np.mean((rv / V_SCALE) ** 2) + W_DVDQ * np.mean((rd / D_SCALE) ** 2)
    return j if np.isfinite(j) else 1e12


def fit_rpt(data_dir, cell, rpt, base_theta=None, seeds=(0, 1, 2), maxiter=400, popsize=20):
    """Differential evolution (+ L-BFGS-B polish) from several seeds; the base
    solution, when given, is one member of the first seeds' initial population.
    Returns a dict: parameters, LLI, RMSEs and cost of the best run."""
    hc = HalfCells(data_dir)
    q, v, dvdq, q_full, _ = load_discharge(rpt_files(data_dir, cell)[rpt])
    best = None
    for k, seed in enumerate(seeds):
        x0 = None if base_theta is None or k == len(seeds) - 1 else np.clip(base_theta, *np.array(BOUNDS).T)
        r = optimize.differential_evolution(cost, BOUNDS, args=(hc, q, v, dvdq), seed=seed, x0=x0,
                                            popsize=popsize, maxiter=maxiter, tol=1e-8, init="sobol",
                                            polish=True, updating="deferred", workers=1)
        if best is None or r.fun < best.fun:
            best = r
    theta = best.x
    rv, rd = residuals(hc, theta, q, v, dvdq)
    out = dict(zip(PARAMS, theta))
    out.update(cell=cell, rpt=rpt, cost=float(best.fun), rmse_voltage_V=float(np.sqrt(np.mean(rv ** 2))),
               rmse_dVdQ_V_per_Ah=float(np.sqrt(np.mean(rd ** 2))), measured_capacity_Ah=q_full,
               Cn_total_Ah=theta[0] + theta[1], LLI_Ah=theta[3] * (theta[0] + theta[1]) + theta[4] * theta[2])
    return out


def base_check(data_dir, cell, rpt, base_row):
    """RMSE between this forward model, evaluated with the base parameters, and the
    base fit's own voltage: validates the OCP orientations and the composite model."""
    hc = HalfCells(data_dir)
    q, v, dvdq, _, base_v = load_discharge(rpt_files(data_dir, cell)[rpt])
    theta = [base_row[p] for p in PARAMS]
    return float(np.sqrt(np.mean((hc.voltage(theta, q) - base_v) ** 2)))


def degradation_modes(table, ref_rpt=1):
    """LAM of Si, graphite and the positive electrode, and LLI [%], relative to ref_rpt
    (the formula of PyProBE's calculate_dma_parameters, applied to each capacity)."""
    t = table.sort_values("rpt").reset_index(drop=True)
    ref = t[t["rpt"] == ref_rpt].iloc[0]
    out = t[["rpt"]].copy()
    for col, name in (("Cn_Si_Ah", "LAM_Si"), ("Cn_Gr_Ah", "LAM_Gr"), ("Cp_Ah", "LAM_PE"), ("LLI_Ah", "LLI")):
        out[name] = 100 * (1 - t[col] / ref[col])
    return out
