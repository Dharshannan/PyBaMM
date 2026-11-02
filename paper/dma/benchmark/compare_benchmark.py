"""Benchmark comparison (PyProBE 2.6.1 environment): the paper's composite DMA (old PyProBE) and
pyprobe-dma, each with and without the Si OCP deformation, on the measured C/20 discharges.

    python compare_benchmark.py

Reads ../results/benchmark/{old_def,old_nodef,new_nodef_notails,new_nodef,new_def}.csv and the
new_* fitted curves; writes the figures (JES style) and bench_table.tex to the same folder.
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from scipy.signal import savgol_filter  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import composite_dma as C  # noqa: E402
import run_dma as R  # noqa: E402

DATA = os.path.join(HERE, "..", "..", "data")
OUT = os.path.join(HERE, "..", "results", "benchmark")
VARIANTS = {  # label, style
    "old_def": ("PyProBE 2.6.1, Si OCP deformed", dict(color="tab:blue", ls="-", marker="o", ms=3.6)),
    "old_nodef": ("PyProBE 2.6.1, Si OCP fixed", dict(color="tab:blue", ls="--", marker="o", ms=3.6, mfc="white")),
    "new_def": ("pyprobe-dma (trimmed), Si OCP deformed", dict(color="tab:orange", ls="-", marker="s", ms=3.4)),
    "new_nodef": ("pyprobe-dma (trimmed), Si OCP fixed", dict(color="tab:orange", ls="--", marker="s", ms=3.4, mfc="white")),
}
EXTRA = {  # pyprobe-dma on the full discharge (table only)
    "new_def_full": ("pyprobe-dma (full curve), Si OCP deformed", dict(color="tab:red", ls="-")),
    "new_nodef_full": ("pyprobe-dma (full curve), Si OCP fixed", dict(color="tab:red", ls="--")),
}
LIMITS = (2.5, 4.2)


def load():
    t = {k: pd.read_csv(os.path.join(OUT, f"{k}.csv"))
         for k in (*VARIANTS, *EXTRA, "new_nodef_notails", "new_nodef_notails_full")}
    curves = {k: dict(np.load(os.path.join(OUT, f"{k}_curves.npz"))) for k in ("new_def", "new_nodef", *EXTRA)}
    return t, curves


def model_voltage(variant, row, q, hc, curves):
    """Fitted terminal voltage at discharged capacity q (old: OCV, as fitted; new: OCV + i R)."""
    if variant.startswith("old"):
        return hc.voltage([row[p] for p in C.PARAMS], q)
    c = curves[variant]
    key = f"{row['cell']}_{int(row['rpt'])}"
    qq, vv = c[f"{key}_q"], c[f"{key}_v"]
    ok = np.isfinite(qq) & np.isfinite(vv)
    o = np.argsort(qq[ok])
    return np.interp(q, qq[ok][o], vv[ok][o], left=np.nan, right=np.nan)


def scores(t, curves):
    """Voltage and dV/dQ RMSE of every variant on the common window (the paper DMA's 300 points)."""
    hc = C.HalfCells(DATA)
    rows = []
    for k in (*VARIANTS, *EXTRA):
        for _, r in t[k].iterrows():
            q, v, dvdq, _, _ = C.load_discharge(C.rpt_files(DATA, r["cell"])[int(r["rpt"])])
            vm = model_voltage(k, r, q, hc, curves)
            rows.append(dict(variant=k, cell=r["cell"], rpt=int(r["rpt"]),
                             rmse_V_mV=1e3 * np.sqrt(np.nanmean((vm - v) ** 2)),
                             rmse_dVdQ=np.sqrt(np.nanmean((np.gradient(vm, q) - dvdq) ** 2))))
    return pd.DataFrame(rows)


def modes(t):
    out = []
    base = R.base_tables(DATA)
    for k in (*VARIANTS, *EXTRA):
        for cell in C.FOLDERS:
            b = base[cell]
            efc0 = float(b.loc[b["rpt"] == 1, "EFC"].iloc[0])
            m = C.degradation_modes(t[k][t[k]["cell"] == cell])
            m["efc"] = m["rpt"].map(dict(zip(b["rpt"], b["EFC"] - efc0)))
            m["cell"], m["variant"] = cell, k
            out.append(m)
    return pd.concat(out, ignore_index=True)


def plot_modes(m, variants=None, name="bench_modes"):
    variants = variants or VARIANTS
    st = R._style()
    cells = list(C.FOLDERS)
    fig = plt.figure(figsize=(st.DOUBLE_IN, 7.6), layout="constrained")
    rows = fig.subfigures(len(cells), 1)
    n = 0
    for i, cell in enumerate(cells):
        ax = rows[i].subplots(1, 4)
        rows[i].suptitle(R.TAG[cell], x=0.01, ha="left", fontsize=8, fontweight="bold")
        for j, (col, lab) in enumerate(R.MODES):
            for k, (_, sty) in variants.items():
                d = m[(m["cell"] == cell) & (m["variant"] == k) & (m["rpt"] >= 1)].sort_values("efc")  # from RPT1
                ax[j].plot(d["efc"], d[col], lw=1.0, **sty)
            ax[j].axhline(0, color="k", lw=0.5)
            ax[j].set_ylabel(lab)
            if i == len(cells) - 1:
                ax[j].set_xlabel("EFC since RPT1 [-]")
            st.panel_label(ax[j], n)
            n += 1
    fig.legend(handles=[Line2D([], [], label=lab, **sty) for lab, sty in variants.values()],
               loc="outside upper center", ncol=2, frameon=False)
    return st.save(fig, OUT, name)[0]


def plot_rmse(s, variants=None, name="bench_rmse"):
    variants = variants or VARIANTS
    st = R._style()
    fig, ax = plt.subplots(1, 2, figsize=(st.DOUBLE_IN, 2.5), layout="constrained")
    order = [(c, r) for c in C.FOLDERS for r in sorted(s.loc[s["cell"] == c, "rpt"].unique())]
    x = np.arange(len(order))
    for k, (lab, sty) in variants.items():
        d = s[s["variant"] == k].set_index(["cell", "rpt"]).loc[order]
        ax[0].plot(x, d["rmse_V_mV"], lw=0.8, **sty)
        ax[1].plot(x, d["rmse_dVdQ"], lw=0.8, **sty)
    for a, yl in zip(ax, ("Voltage RMSE [mV]", "dV/dQ RMSE [V/A h]")):
        a.set_yscale("log")
        a.set_ylabel(yl)
        a.set_xticks(x, [f"{c[4:]}-{r}" for c, r in order], rotation=90, fontsize=7)
        a.set_xlabel("Cell–RPT")
    for i, a in enumerate(ax):
        st.panel_label(a, i)
    fig.legend(handles=[Line2D([], [], label=lab, **sty) for lab, sty in variants.values()],
               loc="outside upper center", ncol=2, frameon=False)
    return st.save(fig, OUT, name)[0]


def plot_examples(t, curves, picks=(("CELL064", 1), ("CELL064", 5), ("CELL009", 8), ("CELL026", 8)), variants=None,
                  name="bench_examples", whole=False):
    """Voltage (top) and dV/dQ (bottom) of selected RPTs: measured and the fits in `variants`; whole: over the
    whole discharge, the dropped last 5 % shaded."""
    variants = variants or VARIANTS
    st = R._style()
    hc = C.HalfCells(DATA)
    fig, ax = plt.subplots(2, len(picks), figsize=(st.DOUBLE_IN, 4.3), layout="constrained")
    for j, (cell, rpt) in enumerate(picks):
        if whole:
            _whole_discharge_fits(ax[0, j], ax[1, j], t, curves, hc, cell, rpt, variants)
        else:
            q, v, dvdq, _, _ = C.load_discharge(C.rpt_files(DATA, cell)[rpt])
            ax[0, j].plot(q, v, color="k", lw=1.6)
            ax[1, j].plot(q, dvdq, color="k", lw=1.6)
            for k, (_, sty) in variants.items():
                r = t[k][(t[k]["cell"] == cell) & (t[k]["rpt"] == rpt)].iloc[0]
                vm = model_voltage(k, r, q, hc, curves)
                s = {kk: vv for kk, vv in sty.items() if kk not in ("marker", "ms", "mfc")}
                ax[0, j].plot(q, vm, lw=1.0, **s)
                ax[1, j].plot(q, np.gradient(vm, q), lw=1.0, **s)
            ax[1, j].set_ylim(np.percentile(dvdq, 2) - 0.3, max(0.05, dvdq.max() + 0.1))
        ax[1, j].set_xlabel("Discharge capacity [A h]")
        st.panel_label(ax[0, j], 2 * j, f"{cell} RPT{rpt}")
        st.panel_label(ax[1, j], 2 * j + 1)
    ax[0, 0].set_ylabel("Voltage [V]")
    ax[1, 0].set_ylabel("dV/dQ [V/A h]")
    fig.legend(handles=[Line2D([], [], color="k", lw=1.6, label="Measured")]
               + [Line2D([], [], label=lab, **{kk: vv for kk, vv in sty.items() if kk not in ("marker", "ms", "mfc")})
                  for lab, sty in variants.values()], loc="outside upper center", ncol=3, frameon=False)
    return st.save(fig, OUT, name)[0]


def plot_cell_fits(t, curves, cell, variants=None, name="bench_fits", whole=False):
    """Every RPT of one cell: voltage (upper panel) and dV/dQ (lower panel), measured and the fits in
    `variants`, with each fit's voltage RMSE on the common window; four RPTs per block of two rows. whole: over
    the whole discharge, the dropped last 5 % shaded."""
    variants = variants or VARIANTS
    st = R._style()
    hc = C.HalfCells(DATA)
    rpts = sorted(t["old_def"].loc[t["old_def"]["cell"] == cell, "rpt"].astype(int))
    ncol = 4
    nblk = int(np.ceil(len(rpts) / ncol))
    fig, ax = plt.subplots(2 * nblk, ncol, figsize=(st.DOUBLE_IN, 2.7 * nblk + 0.5), layout="constrained",
                           squeeze=False)
    k = 0
    for n, rpt in enumerate(rpts):
        av, ad = ax[2 * (n // ncol), n % ncol], ax[2 * (n // ncol) + 1, n % ncol]
        if whole:
            rm = [f"{x:.1f}" for x in _whole_discharge_fits(av, ad, t, curves, hc, cell, rpt, variants).values()]
        else:
            q, v, dvdq, _, _ = C.load_discharge(C.rpt_files(DATA, cell)[rpt])
            av.plot(q, v, color="k", lw=1.6)
            ad.plot(q, dvdq, color="k", lw=1.6)
            rm = []
            for var, (_, sty) in variants.items():
                r = t[var][(t[var]["cell"] == cell) & (t[var]["rpt"] == rpt)].iloc[0]
                vm = model_voltage(var, r, q, hc, curves)
                line = {kk: vv for kk, vv in sty.items() if kk not in ("marker", "ms", "mfc")}
                av.plot(q, vm, lw=1.0, **line)
                ad.plot(q, np.gradient(vm, q), lw=1.0, **line)
                rm.append(f"{1e3 * np.sqrt(np.nanmean((vm - v) ** 2)):.1f}")
            ad.set_ylim(np.percentile(dvdq, 2) - 0.3, max(0.05, dvdq.max() + 0.1))
        av.text(0.03, 0.04, "RMSE [mV]\n" + " / ".join(rm), transform=av.transAxes, ha="left", va="bottom",
                fontsize=7)
        st.panel_label(av, k, f"RPT{rpt}")
        st.panel_label(ad, k + 1)
        k += 2
        ad.set_xlabel("Discharge capacity [A h]")
    for a in ax.flat:
        if not a.lines:
            a.set_visible(False)
    for b_ in range(nblk):
        ax[2 * b_, 0].set_ylabel("Voltage [V]")
        ax[2 * b_ + 1, 0].set_ylabel("dV/dQ [V/A h]")
    fig.legend(handles=[Line2D([], [], color="k", lw=1.6, label="Measured")]
               + [Line2D([], [], label=lab, **{kk: vv for kk, vv in sty.items() if kk not in ("marker", "ms", "mfc")})
                  for lab, sty in variants.values()], loc="outside upper center", ncol=3, frameon=False)
    return st.save(fig, OUT, f"{name}_{cell}")[0]


def plot_window_full(t, curves, cell="CELL064", rpt=5):
    """pyprobe-dma's fitted OCV on its SoC axis, trimmed against full curve (Si OCP deformed): both run from
    V_min to V_max; the trimmed data end above 0 % SoC, the full-curve data at it."""
    st = R._style()
    fig, ax = plt.subplots(1, 2, figsize=(st.DOUBLE_IN, 2.7), layout="constrained", sharey=True)
    d = pd.read_csv(C.rpt_files(DATA, cell)[rpt], usecols=["Q_discharge_Ah", "voltage_fit_used_V"])
    for j, (k, (lab, sty)) in enumerate(FULL_STYLE.items()):
        a = ax[j]
        rn = t[k][(t[k]["cell"] == cell) & (t[k]["rpt"] == rpt)].iloc[0]
        c, key, cap = curves[k], f"{cell}_{rpt}", rn["cell_capacity_Ah"]
        a.plot(rn["soc_max"] - d["Q_discharge_Ah"] / cap, d["voltage_fit_used_V"], color="k", lw=1.0,
               label="Measured (on the fit's SoC axis)")
        a.plot(rn["soc_max"] - c[f"{key}_q"] / cap, c[f"{key}_ocv"], color=sty["color"], label="Fitted OCV")
        a.plot([0, 1], LIMITS, "o", color=sty["color"], mfc="white", ms=4.5, label="SoC 0 / 100 %")
        for lim in LIMITS:
            a.axhline(lim, color="0.5", ls=":", lw=0.8)
        a.set_xlabel("SoC [-]")
        st.panel_label(a, j, lab.replace("pyprobe-dma, ", ""))
    ax[0].set_ylabel("Voltage [V]")
    ax[0].legend(loc="lower right", fontsize=7)
    return st.save(fig, OUT, "full_window")[0]


def plot_deformation_full(t):
    """Fitted Si OCP scale and shift: full curve against trimmed (pyprobe-dma, Si OCP deformed)."""
    st = R._style()
    fig, ax = plt.subplots(1, 2, figsize=(st.DOUBLE_IN, 2.5), layout="constrained")
    a, b = t["new_def"].set_index(["cell", "rpt"]), t["new_def_full"].set_index(["cell", "rpt"])
    for cell, mk in zip(C.FOLDERS, ("o", "s", "^", "D")):
        idx = [i for i in a.index if i[0] == cell]
        ax[0].plot(a.loc[idx, "si_ocp_scale"], b.loc[idx, "si_ocp_scale"], mk, ms=4, label=cell, ls="none")
        ax[1].plot(1e3 * a.loc[idx, "si_ocp_shift_V"], 1e3 * b.loc[idx, "si_ocp_shift_V"], mk, ms=4, ls="none")
    for x, lab, lim in ((ax[0], "scale $a$ [-]", (0.74, 1.11)), (ax[1], "shift $b$ [mV]", (-105, 35))):
        x.plot(lim, lim, color="0.5", lw=0.8)
        x.set_xlim(lim)
        x.set_ylim(lim)
        x.set_xlabel(f"Trimmed {lab}")
        x.set_ylabel(f"Full curve {lab}")
    for i, x in enumerate(ax):
        st.panel_label(x, i)
    ax[0].legend(loc="lower right")  # the upper left holds a point on the scale bound
    return st.save(fig, OUT, "full_deformation")[0]


def plot_window(t, curves, cell="CELL064", rpt=5):
    """The OCV window of each fit against SoC: pyprobe-dma's OCV runs from V_min to V_max by
    construction; the old fits' OCV at the ends of the measured discharge is left free."""
    st = R._style()
    hc = C.HalfCells(DATA)
    fig, ax = plt.subplots(1, 2, figsize=(st.DOUBLE_IN, 2.7), layout="constrained", sharey=True)
    for j, k in enumerate(("def", "nodef")):
        a = ax[j]
        ro = t[f"old_{k}"][(t[f"old_{k}"]["cell"] == cell) & (t[f"old_{k}"]["rpt"] == rpt)].iloc[0]
        q_full = ro["measured_capacity_Ah"]
        qq = np.linspace(0, q_full, 300)
        c = curves[f"new_{k}"]
        key = f"{cell}_{rpt}"
        rn = t[f"new_{k}"][(t[f"new_{k}"]["cell"] == cell) & (t[f"new_{k}"]["rpt"] == rpt)].iloc[0]
        a.plot(rn["soc_max"] - qq / rn["cell_capacity_Ah"], hc.voltage([ro[p] for p in C.PARAMS], qq),
               color="tab:blue", label="PyProBE 2.6.1: OCV over the measured discharge")
        qn, ocv = c[f"{key}_q"], c[f"{key}_ocv"]
        cap = rn["cell_capacity_Ah"]
        soc_meas = rn["soc_max"] - qn / cap
        a.plot(soc_meas, ocv, color="tab:orange", label="pyprobe-dma: fitted OCV")
        a.plot([0, 1], LIMITS, "o", color="tab:orange", mfc="white", ms=4.5, label="pyprobe-dma: SoC 0 / 100 %")
        d = pd.read_csv(C.rpt_files(DATA, cell)[rpt], usecols=["Q_discharge_Ah", "voltage_fit_used_V"])
        a.plot(rn["soc_max"] - d["Q_discharge_Ah"] / cap, d["voltage_fit_used_V"], color="k", lw=1.0,
               label="Measured (on pyprobe-dma's SoC axis)")
        for lim in LIMITS:
            a.axhline(lim, color="0.5", ls=":", lw=0.8)
        a.set_xlabel("SoC [-]")
        st.panel_label(a, j, f"Si OCP {'deformed' if k == 'def' else 'fixed'}")
    ax[0].set_ylabel("Voltage [V]")
    ax[0].legend(loc="lower right", fontsize=7)
    return st.save(fig, OUT, "bench_window")[0]


TRIM_STYLE = {  # trimmed against full-curve pyprobe-dma (Si OCP deformed), with the PyProBE 2.6.1 reference
    "old_def": ("PyProBE 2.6.1 (trimmed)", dict(color="tab:blue", ls="-", marker="o", ms=3.4)),
    "new_def": ("pyprobe-dma, trimmed", dict(color="tab:orange", ls="-", marker="s", ms=3.4)),
    "new_def_full": ("pyprobe-dma, full curve", dict(color="tab:red", ls="--", marker="^", ms=3.6, mfc="white")),
}
R_BOUNDS = (0.0, 0.3)  # pyprobe-dma's constant-resistance box [ohm]
FULL_STYLE = {k: TRIM_STYLE[k] for k in ("new_def", "new_def_full")}  # the appendix: pyprobe-dma only


def _whole_discharge_fits(av, ad, t, curves, hc, cell, rpt, variants):
    """Measured voltage and dV/dQ over the whole discharge (uniform grid, the fits' Savitzky-Golay window) and
    the fits in `variants` on it; the dropped last 5 % shaded. Returns {variant: voltage RMSE [mV] on the
    common fit window}."""
    q, v, dvdq, _, _ = C.load_discharge(C.rpt_files(DATA, cell)[rpt])
    d = pd.read_csv(C.rpt_files(DATA, cell)[rpt], usecols=["Q_discharge_Ah", "voltage_fit_used_V"])
    d = d.sort_values("Q_discharge_Ah").drop_duplicates("Q_discharge_Ah")
    qa = np.linspace(0, d["Q_discharge_Ah"].max(), 1500)
    va = np.interp(qa, d["Q_discharge_Ah"], d["voltage_fit_used_V"])
    dq = qa[1] - qa[0]
    n = int(C.SG_WINDOW_AH / dq) // 2 * 2 + 1
    av.plot(qa, va, color="k", lw=1.6)
    ad.plot(qa, savgol_filter(va, n, 3, deriv=1, delta=dq), color="k", lw=1.6)
    rmse = {}
    for k, (_, sty) in variants.items():
        r = t[k][(t[k]["cell"] == cell) & (t[k]["rpt"] == rpt)].iloc[0]
        line = {kk: vv for kk, vv in sty.items() if kk not in ("marker", "ms", "mfc")}
        vm = model_voltage(k, r, qa, hc, curves)
        ok = np.isfinite(vm)  # the trimmed pyprobe-dma curve stops at the end of its fit window
        av.plot(qa, vm, lw=1.0, **line)
        ad.plot(qa[ok], savgol_filter(vm[ok], n, 3, deriv=1, delta=dq), lw=1.0, **line)
        rmse[k] = 1e3 * np.sqrt(np.nanmean((model_voltage(k, r, q, hc, curves) - v) ** 2))
    for a in (av, ad):
        a.axvspan(q.max(), qa.max(), color="0.88", lw=0)
    ad.set_ylim(np.percentile(dvdq, 2) - 0.3, max(0.05, dvdq.max() + 0.1))
    return rmse


def plot_trim(t, curves, s, m, pick=("CELL064", 5)):
    """Why the fit window is trimmed (Si OCP deformed): (a) voltage RMSE of every RPT on the common window and
    (b) fitted resistance, trimmed against full curve; (c, d) one late-life RPT over the whole discharge, the
    dropped last 5 % shaded; (e-h) Si LAM of each cell, with the PyProBE 2.6.1 fit as reference."""
    st = R._style()
    hc = C.HalfCells(DATA)
    fig = plt.figure(figsize=(st.DOUBLE_IN, 6.6), layout="constrained")
    top, mid, bot = fig.subfigures(3, 1, height_ratios=[1, 1, 0.9])
    ax_r = top.subplots(1, 2)
    order = [(c, r) for c in C.FOLDERS for r in sorted(s.loc[s["cell"] == c, "rpt"].unique())]
    x = np.arange(len(order))
    for k in ("new_def", "new_def_full"):
        lab, sty = TRIM_STYLE[k]
        d = s[s["variant"] == k].set_index(["cell", "rpt"]).loc[order]
        ax_r[0].plot(x, d["rmse_V_mV"], lw=0.8, **sty)
        rr = t[k].set_index(["cell", "rpt"]).reindex(order)["R_ohm"]
        ax_r[1].plot(x, 1e3 * rr.to_numpy(), lw=0.8, **sty)
    for b in R_BOUNDS:
        ax_r[1].axhline(1e3 * b, color="0.5", ls=":", lw=0.8)
    ax_r[0].set_yscale("log")
    ax_r[0].set_ylabel("Voltage RMSE [mV]")
    ax_r[1].set_ylabel(r"Fitted $R$ [m$\Omega$]")
    for a in ax_r:
        a.set_xticks(x, [f"{c[4:]}-{r}" for c, r in order], rotation=90, fontsize=7)
        a.set_xlabel("Cell–RPT")
    ax_e = mid.subplots(1, 2)
    cell, rpt = pick
    _whole_discharge_fits(ax_e[0], ax_e[1], t, curves, hc, cell, rpt, TRIM_STYLE)
    ax_e[0].text(0.985, 0.55, "dropped\nlast 5 %", transform=ax_e[0].transAxes, ha="right", fontsize=7)
    for a in ax_e:
        a.set_xlabel("Discharge capacity [A h]")
    ax_e[0].set_ylabel("Voltage [V]")
    ax_e[1].set_ylabel("dV/dQ [V/A h]")
    ax_l = bot.subplots(1, len(C.FOLDERS))
    for j, c in enumerate(C.FOLDERS):
        for k, (_, sty) in TRIM_STYLE.items():
            dd = m[(m["cell"] == c) & (m["variant"] == k) & (m["rpt"] >= 1)].sort_values("efc")
            ax_l[j].plot(dd["efc"], dd["LAM_Si"], lw=1.0, **sty)
        ax_l[j].set_xlabel("EFC since RPT1 [-]")
        ax_l[j].text(0.04, 0.96, R.TAG[c].replace(", ", "\n", 1), transform=ax_l[j].transAxes, ha="left", va="top",
                     fontsize=7)
    ax_l[0].set_ylabel("LAM, Si [%]")
    n = 0
    for a in [*ax_r, *ax_e, *ax_l]:
        st.panel_label(a, n, f"{cell} RPT{rpt}" if a is ax_e[0] else None)
        n += 1
    fig.legend(handles=[Line2D([], [], color="k", lw=1.6, label="Measured")]
               + [Line2D([], [], label=lab, **sty) for lab, sty in TRIM_STYLE.values()],
               loc="outside upper center", ncol=4, frameon=False)
    return st.save(fig, OUT, "bench_trim")[0]


def plot_deformation(t):
    st = R._style()
    fig, ax = plt.subplots(1, 2, figsize=(st.DOUBLE_IN, 2.5), layout="constrained")
    a, b = t["old_def"].set_index(["cell", "rpt"]), t["new_def"].set_index(["cell", "rpt"])
    for cell, mk in zip(C.FOLDERS, ("o", "s", "^", "D")):
        idx = [i for i in a.index if i[0] == cell]
        ax[0].plot(a.loc[idx, "si_ocp_scale"], b.loc[idx, "si_ocp_scale"], mk, ms=4, label=cell, ls="none")
        ax[1].plot(1e3 * a.loc[idx, "si_ocp_shift_V"], 1e3 * b.loc[idx, "si_ocp_shift_V"], mk, ms=4, ls="none")
    for x, lab, lim in ((ax[0], "scale $a$ [-]", (0.74, 1.11)), (ax[1], "shift $b$ [mV]", (-105, 35))):
        x.plot(lim, lim, color="0.5", lw=0.8)
        x.set_xlim(lim)
        x.set_ylim(lim)
        x.set_xlabel(f"PyProBE 2.6.1 {lab}")
        x.set_ylabel(f"pyprobe-dma {lab}")
    for i, x in enumerate(ax):
        st.panel_label(x, i)
    ax[0].legend(loc="upper left")
    return st.save(fig, OUT, "bench_deformation")[0]


def latex_table(t, s, m):
    """Per variant: fits that converged, median RMSE, and the differences from the paper DMA's
    setting (PyProBE 2.6.1 with deformation) in LLI and Si LAM."""
    ref = m[m["variant"] == "old_def"].set_index(["cell", "rpt"])
    lines = []
    for key, lab in (("new_nodef_notails", "pyprobe-dma (trimmed), OCPs as measured"),
                     ("new_nodef_notails_full", "pyprobe-dma (full curve), OCPs as measured")):
        nt = t[key]
        lines.append(lab + f" & {np.isfinite(nt['cost']).sum()} / {len(nt)}" + r" & -- & -- & -- & -- \\")
    lines.append(r"\midrule")
    for k, (lab, _) in {**VARIANTS, **EXTRA}.items():
        ok = np.isfinite(t[k]["cost"]).sum() if "cost" in t[k] else len(t[k])
        sk = s[s["variant"] == k]
        mk = m[m["variant"] == k].set_index(["cell", "rpt"])
        dl = (mk["LLI"] - ref["LLI"]).abs()
        ds = (mk["LAM_Si"] - ref["LAM_Si"]).abs()
        lines.append(f"{lab} & {ok} / {len(t[k])} & {sk['rmse_V_mV'].median():.1f} & {sk['rmse_dVdQ'].median():.3f} & "
                     f"{dl[dl.index.get_level_values(1) >= 1].mean():.1f} & {ds[ds.index.get_level_values(1) >= 1].mean():.1f} "
                     + r"\\")
    body = "\n".join(lines)
    with open(os.path.join(OUT, "bench_table.tex"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(r"\begin{tabular}{lrrrrr}" + "\n" + r"\toprule" + "\n"
                 + r"Variant & Fits & V RMSE [mV] & dV/dQ RMSE & $|\Delta$LLI$|$ [pp] & $|\Delta$LAM Si$|$ [pp] \\"
                 + "\n" + r"\midrule" + "\n" + body + "\n" + r"\bottomrule" + "\n" + r"\end{tabular}" + "\n")


def main():
    t, curves = load()
    s = scores(t, curves)
    s.to_csv(os.path.join(OUT, "bench_scores.csv"), index=False)
    m = modes(t)
    m.to_csv(os.path.join(OUT, "bench_modes.csv"), index=False)
    figs = [plot_modes(m), plot_rmse(s), plot_examples(t, curves), plot_window(t, curves), plot_deformation(t),
            plot_trim(t, curves, s, m)]
    # trimmed against full curve, Si OCP deformed: every figure again
    figs += [plot_rmse(s, FULL_STYLE, "full_rmse"), plot_modes(m, FULL_STYLE, "full_modes"),
             plot_examples(t, curves, variants=FULL_STYLE, name="full_examples", whole=True),
             plot_window_full(t, curves), plot_deformation_full(t)]
    figs += [plot_cell_fits(t, curves, cell, FULL_STYLE, "full_fits", whole=True) for cell in C.FOLDERS]
    for f in figs + [plot_cell_fits(t, curves, cell) for cell in C.FOLDERS]:
        print(f)
    latex_table(t, s, m)
    pd.set_option("display.width", 200)
    print(s.groupby("variant")[["rmse_V_mV", "rmse_dVdQ"]].median().round(3))
    w = m[m["rpt"] >= 1].pivot_table(index=["cell", "rpt"], columns="variant", values=["LAM_Si", "LLI", "LAM_Gr", "LAM_PE"])
    print(w.round(1).to_string())


if __name__ == "__main__":
    main()
