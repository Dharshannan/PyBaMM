"""Composite DMA of the C/20 RPT discharges with pyprobe-dma (ImperialCollegeLondon/pyprobe-dma,
built on the PyProBE "version 3" branch). Benchmark only: not used in the paper.

Runs in its own environment (benchmark/environment-v3.yml): pyprobe-dma pins a PyProBE branch
that would replace the PyProBE 2.6.1 the paper's DMA (../composite_dma.py) uses.

What pyprobe-dma does differently from ../composite_dma.py:
- The voltage window is part of the model: 0 % and 100 % SoC are where the full-cell OCV equals
  the cell's voltage limits, so the free variables are only the negative electrode's bottom
  stoichiometry, its capacity and the cell capacity (plus the graphite fraction and a resistance);
  the positive electrode's window follows from the limits.
- A resistance model (here constant) carries the C/20 overpotential, and the measured curve is
  placed on the SoC axis through its first and last points (TwoPointSOCEstimator).
- No Si OCP deformation: with `deform=True` it is fitted around pyprobe-dma by an outer search
  over (scale, shift), each candidate being a full pyprobe-dma fit with the deformed Si OCP.
"""
import os
import sys

import numpy as np
import pandas as pd
import polars as pl
import pyprobe_dma as dma
from pyprobe.result import Table
from scipy import optimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import tables  # noqa: E402  (pandas only)

VOLTAGE_LIMITS = (2.5, 4.2)  # cell voltage window (V_min, V_max): the OCV at 0 % and 100 % SoC
W_DVDQ, V_SCALE, D_SCALE = 0.5, 0.01, 0.25  # as ../composite_dma.py; dV/dSoC scale ~ 0.1 V/Ah x 2.5 Ah
FIT_POINTS = 300


def _ocp(data_dir, name):
    """Half-cell OCP, downsampled to one point every 2 mV as in the pyprobe-dma example: the
    graphite and NMC tables are quantised (repeated potentials), which otherwise makes the
    analytic dV/dSoC of the piecewise-linear OCP spiky."""
    from pyprobe.analysis.smoothing import downsample

    d = pd.read_csv(os.path.join(data_dir, "ocp_halfcell", name)).sort_values("sto").drop_duplicates("sto")
    t = downsample(input_data=Table(lf=pl.DataFrame({"z": d["sto"].to_numpy(), "OCV": d["p"].to_numpy()}),
                                    metadata={}), target_column="OCV", sampling_interval=0.002)
    s, p = t.get("z"), t.get("OCV")
    o = np.argsort(s)
    s, p = s[o], p[o]
    keep = np.r_[True, np.diff(s) > 0] & np.r_[True, np.diff(p) != 0]
    return s[keep], p[keep]


ANODE_TAIL_V = 1.5  # anode OCP at zero stoichiometry (see _tail)
CATHODE_TAIL_V = 4.45  # cathode OCP at zero lithium: the NMC table stops at 4.293 V


def _tail(sto, ocp):
    """Close the delithiation curve at zero stoichiometry with a steep tail to ANODE_TAIL_V.
    The half-cell tables stop at ~1.0 V (graphite) and 0.98 V (Si), but with the cell's lower
    voltage limit enforced, a cell that has lost much of its lithium must reach 2.5 V with the
    anode nearly empty, i.e. above 1 V: without the tail no window fits late in life."""
    if sto[0] > 0:
        return np.r_[0.0, sto], np.r_[ANODE_TAIL_V, ocp]
    return sto, np.r_[ANODE_TAIL_V, ocp[1:]]


class HalfCells:
    """Discharge half-cell OCPs (as the paper's DMA): NMC622 lithiating (file sto = 1 - y),
    graphite and Si (Verbrugge) delithiating."""

    def __init__(self, data_dir, tails=True):
        s, p = _ocp(data_dir, "nmc622_lithiation.csv")
        y, u = (1 - s)[::-1], p[::-1]  # y = 1 - sto, ascending
        if tails:
            y, u = (np.r_[0.0, y], np.r_[CATHODE_TAIL_V, u]) if y[0] > 0 else (y, np.r_[CATHODE_TAIL_V, u[1:]])
        self.pe = dma.OCP.from_arrays(stoichiometry=y, ocp=u)
        self.x_gr, self.u_gr = _ocp(data_dir, "graphite_delithiation.csv")
        self.x_si, self.u_si = _ocp(data_dir, "silicon_delithiation_verbrugge.csv")
        if tails:
            self.x_gr, self.u_gr = _tail(self.x_gr, self.u_gr)
            self.x_si, self.u_si = _tail(self.x_si, self.u_si)

    def ocv_model(self, scale=1.0, shift=0.0):
        ne = dma.CompositeOCP.from_arrays(ocp_comp1=self.u_gr, sto_comp1=self.x_gr,
                                          ocp_comp2=scale * self.u_si + shift, sto_comp2=self.x_si)
        return dma.OCVModel(ocp_pe=self.pe, ocp_ne=ne, voltage_limits=VOLTAGE_LIMITS,
                            soc_estimator=dma.TwoPointSOCEstimator(),
                            resistance_model=dma.ResistanceModel.from_template(template="constant"))


TRIM = True  # fit the discharge without its last 5 % of capacity, as the paper's DMA does


def curves(path, trim=None):
    """Voltage and current curves on the net-capacity axis (decreasing on discharge). trim (default
    TRIM): only the paper DMA's objective window, i.e. the discharge without its last 5 % of
    capacity; the model's SoC window stays pinned to the voltage limits either way."""
    trim = TRIM if trim is None else trim
    d = pd.read_csv(path, usecols=["Q_discharge_Ah", "voltage_fit_used_V", "current_A", "fit_objective_included"])
    d = d.sort_values("Q_discharge_Ah").drop_duplicates("Q_discharge_Ah")
    if trim:
        d = d[d["fit_objective_included"].astype(bool)]
    q = d["Q_discharge_Ah"].to_numpy()
    q_grid = np.linspace(q.min(), q.max(), 1500)
    t = Table(lf=pl.DataFrame({"Net Capacity / Ah": -q_grid[::-1],
                               "Voltage / V": np.interp(q_grid, q, d["voltage_fit_used_V"].to_numpy())[::-1],
                               "Current / A": np.interp(q_grid, q, d["current_A"].to_numpy())[::-1]}))
    return t.to_curve(y="Voltage / V", x="Net Capacity / Ah"), t.to_curve(y="Current / A", x="Net Capacity / Ah")


def cost_function():
    rmse_v = dma.cost_functions.RMSE(dma.cost_functions.voltage)
    rmse_d = dma.cost_functions.RMSE(dma.cost_functions.dv_dq)
    return (rmse_v / V_SCALE) ** 2 + W_DVDQ * (rmse_d / D_SCALE) ** 2


def constraints():
    return [dma.Constraint("ne_capacity", lo=0.5, hi=4.0),  # aged cells exceed the default 1.2 x cell capacity
            # A bound on pe_capacity (a derived quantity, sent as a nonlinear constraint) would match
            # ../composite_dma.py, but pyprobe-dma 0.1.0 then calls the resistance model without its
            # parameter and fails; the best of several seeds (fit_best) avoids the wrong basin instead.
            dma.Constraint("r_params", lo=[0.0], hi=[0.3])]


def fit(hc, vc, ic, scale=1.0, shift=0.0, seed=0, maxiter=300):
    return dma.run_ocv_curve_fit(voltage_curve=vc, current_curve=ic, ocv_model=hc.ocv_model(scale, shift),
                                 cost_function=cost_function(), constraints=constraints(),
                                 optimizer_options={"seed": seed, "maxiter": maxiter, "popsize": 20, "tol": 1e-8})


SEEDS = (0, 1, 2, 3)


def fit_best(hc, vc, ic, scale=1.0, shift=0.0):
    """Lowest-cost fit over SEEDS (as ../composite_dma.py, which keeps the best of three seeds)."""
    return min((fit(hc, vc, ic, scale, shift, seed=k) for k in SEEDS), key=lambda r: r.optimization_result.fun)


def summarise(res, cell, rpt, scale, shift):
    """One row: capacities (Si and graphite from the fitted graphite fraction), lithium inventory,
    resistance, the Si OCP deformation used, and the cost."""
    lim = res.stoichiometry_limits_result.data.to_dicts()[0]
    frac = lim["NE Component 1 Fraction"]
    frac_gr = float(frac[0] if isinstance(frac, (list, tuple)) else frac)
    cn = lim["NE Capacity / Ah"]
    return dict(cell=cell, rpt=rpt, Cn_Si_Ah=(1 - frac_gr) * cn, Cn_Gr_Ah=frac_gr * cn, Cn_total_Ah=cn,
                Cp_Ah=lim["PE Capacity / Ah"], LLI_Ah=lim["Li Inventory / Ah"],
                cell_capacity_Ah=lim["Cell Capacity / Ah"], R_ohm=lim["Resistance Model Parameter 0"],
                soc_min=lim["Min SOC"], soc_max=lim["Max SOC"], si_ocp_scale=scale, si_ocp_shift_V=shift,
                cost=float(res.optimization_result.fun))


def fit_curve(res):
    """(discharged capacity [Ah], fitted terminal voltage [V], fitted OCV [V]) of a fit."""
    f = res.fit_result.data.sort("Net Capacity / Ah", descending=True)
    q = -f["Net Capacity / Ah"].to_numpy()
    return q - np.nanmin(q), f["Fitted Voltage / V"].to_numpy(), f["Fitted OCV / V"].to_numpy()


def refit_curves(data_dir, out_dir, jobs=12):
    """Rebuild the fitted curves of each variant from its table: one pyprobe-dma fit per RPT at the
    stored Si OCP scale and shift (deterministic, so the same fit as the search found)."""
    from joblib import Parallel, delayed

    sys.path.insert(0, os.path.join(HERE, ".."))
    import composite_dma as C

    def one(cell, rpt, scale, shift, tails):
        hc = HalfCells(data_dir, tails=tails)
        vc, ic = curves(C.rpt_files(data_dir, cell)[rpt])
        return fit_curve(fit(hc, vc, ic, scale, shift))

    for label, tails in (("new_nodef_notails", False), ("new_nodef", True), ("new_def", True)):
        t = pd.read_csv(os.path.join(out_dir, f"{label}.csv"))
        out = Parallel(n_jobs=jobs)(delayed(one)(r["cell"], int(r["rpt"]), r["si_ocp_scale"], r["si_ocp_shift_V"], tails)
                                    for _, r in t.iterrows())
        np.savez_compressed(os.path.join(out_dir, f"{label}_curves.npz"),
                            **{f"{r['cell']}_{int(r['rpt'])}_{k}": arr for (_, r), o in zip(t.iterrows(), out)
                               for k, arr in zip(("q", "v", "ocv"), o)})
        print(label, "curves rebuilt", flush=True)


SCALES = (0.80, 0.86, 0.92, 0.98, 1.04)
SHIFTS = (-0.06, -0.035, -0.01, 0.015, 0.04)
DEFORM_BOUNDS = ((0.75, 1.1), (-0.10, 0.03))  # as ../composite_dma.py


def fit_rpt(data_dir, cell, rpt, deform, tails=True):
    """pyprobe-dma fit of one RPT; with deform=True the Si OCP scale and shift are found by a 5 x 5
    grid and a bounded Nelder-Mead refinement, each evaluation a full pyprobe-dma fit."""
    sys.path.insert(0, os.path.join(HERE, ".."))
    import composite_dma as C

    hc = HalfCells(data_dir, tails=tails)
    vc, ic = curves(C.rpt_files(data_dir, cell)[rpt])
    if not deform:
        res = fit_best(hc, vc, ic)
        return summarise(res, cell, rpt, 1.0, 0.0), fit_curve(res)
    memo = {}

    def outer(p):
        p = tuple(np.round(np.clip(p, *np.array(DEFORM_BOUNDS).T), 4))
        if p not in memo:
            memo[p] = fit(hc, vc, ic, *p)
        return memo[p].optimization_result.fun

    grid = [(a, b) for a in SCALES for b in SHIFTS]
    start = min(grid, key=outer)
    r = optimize.minimize(outer, start, method="Nelder-Mead",
                          options={"maxfev": 25, "xatol": 2e-3, "fatol": 1e-4,
                                   "initial_simplex": [start, (start[0] + 0.03, start[1]), (start[0], start[1] + 0.015)]})
    best = min(memo, key=lambda k: memo[k].optimization_result.fun)
    res = fit_best(hc, vc, ic, *best)  # the search used one seed; the chosen deformation gets them all
    return summarise(res, cell, rpt, *best), fit_curve(res)


def run_all(data_dir, out_dir, jobs=12):
    """Every variant for every C/20 RPT of the four cells: <out>/new_nodef_notails.csv (half-cell
    OCPs as measured), new_nodef.csv and new_def.csv (with the tails), and the fitted curves
    new_<variant>_curves.npz."""
    from joblib import Parallel, delayed

    sys.path.insert(0, os.path.join(HERE, ".."))
    import composite_dma as C

    os.makedirs(out_dir, exist_ok=True)
    todo = [(cell, rpt) for cell in C.FOLDERS for rpt in sorted(C.rpt_files(data_dir, cell))]
    for label, deform, tails in (("new_nodef_notails", False, False), ("new_nodef", False, True),
                                 ("new_def", True, True)):
        out = Parallel(n_jobs=jobs, verbose=2)(delayed(fit_rpt)(data_dir, c, r, deform, tails) for c, r in todo)
        pd.DataFrame([o[0] | {"variant": label} for o in out]).to_csv(os.path.join(out_dir, f"{label}.csv"),
                                                                       index=False)
        np.savez_compressed(os.path.join(out_dir, f"{label}_curves.npz"),
                            **{f"{c}_{r}_{k}": arr for (c, r), o in zip(todo, out)
                               for k, arr in zip(("q", "v", "ocv"), o[1])})
        print(label, "written", flush=True)


if __name__ == "__main__":
    _data, _out = os.path.join(HERE, "..", "..", "data"), os.path.join(HERE, "..", "results", "benchmark")
    if "--curves-only" in sys.argv:
        refit_curves(_data, _out)
    else:
        run_all(_data, _out)
