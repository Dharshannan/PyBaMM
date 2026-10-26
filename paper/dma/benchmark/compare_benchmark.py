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

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import composite_dma as C  # noqa: E402
import run_dma as R  # noqa: E402

DATA = os.path.join(HERE, "..", "..", "data")
OUT = os.path.join(HERE, "..", "results", "benchmark")
VARIANTS = {  # label, style
    "old_def": ("PyProBE 2.6.1, Si OCP deformed", dict(color="tab:blue", ls="-", marker="o", ms=3.6)),
    "old_nodef": ("PyProBE 2.6.1, Si OCP fixed", dict(color="tab:blue", ls="--", marker="o", ms=3.6, mfc="white")),
    "new_def": ("pyprobe-dma, Si OCP deformed", dict(color="tab:orange", ls="-", marker="s", ms=3.4)),
    "new_nodef": ("pyprobe-dma, Si OCP fixed", dict(color="tab:orange", ls="--", marker="s", ms=3.4, mfc="white")),
}
LIMITS = (2.5, 4.2)


def load():
    t = {k: pd.read_csv(os.path.join(OUT, f"{k}.csv")) for k in (*VARIANTS, "new_nodef_notails")}
    curves = {k: dict(np.load(os.path.join(OUT, f"{k}_curves.npz"))) for k in ("new_def", "new_nodef")}
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
    for k in VARIANTS:
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
    for k in VARIANTS:
        for cell in C.FOLDERS:
            b = base[cell]
            efc0 = float(b.loc[b["rpt"] == 1, "EFC"].iloc[0])
            m = C.degradation_modes(t[k][t[k]["cell"] == cell])
            m["efc"] = m["rpt"].map(dict(zip(b["rpt"], b["EFC"] - efc0)))
            m["cell"], m["variant"] = cell, k
            out.append(m)
    return pd.concat(out, ignore_index=True)


def plot_modes(m):
    st = R._style()
    cells = list(C.FOLDERS)
    fig = plt.figure(figsize=(st.DOUBLE_IN, 7.6), layout="constrained")
    rows = fig.subfigures(len(cells), 1)
    n = 0
    for i, cell in enumerate(cells):
        ax = rows[i].subplots(1, 4)
        rows[i].suptitle(R.TAG[cell], x=0.01, ha="left", fontsize=8, fontweight="bold")
        for j, (col, lab) in enumerate(R.MODES):
            for k, (_, sty) in VARIANTS.items():
                d = m[(m["cell"] == cell) & (m["variant"] == k)].sort_values("efc")
                ax[j].plot(d["efc"], d[col], lw=1.0, **sty)
            ax[j].axhline(0, color="k", lw=0.5)
            ax[j].set_ylabel(lab)
            if i == len(cells) - 1:
                ax[j].set_xlabel("EFC since RPT1 [-]")
            st.panel_label(ax[j], n)
            n += 1
    fig.legend(handles=[Line2D([], [], label=lab, **sty) for lab, sty in VARIANTS.values()],
               loc="outside upper center", ncol=4, frameon=False)
    return st.save(fig, OUT, "bench_modes")[0]


def plot_rmse(s):
    st = R._style()
    fig, ax = plt.subplots(1, 2, figsize=(st.DOUBLE_IN, 2.5), layout="constrained")
    order = [(c, r) for c in C.FOLDERS for r in sorted(s.loc[s["cell"] == c, "rpt"].unique())]
    x = np.arange(len(order))
    for k, (lab, sty) in VARIANTS.items():
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
    fig.legend(handles=[Line2D([], [], label=lab, **sty) for lab, sty in VARIANTS.values()],
               loc="outside upper center", ncol=4, frameon=False)
    return st.save(fig, OUT, "bench_rmse")[0]


def plot_examples(t, curves, picks=(("CELL064", 1), ("CELL064", 5), ("CELL009", 8), ("CELL026", 8))):
    """Voltage (top) and dV/dQ (bottom) of selected RPTs: measured and the four fits."""
    st = R._style()
    hc = C.HalfCells(DATA)
    fig, ax = plt.subplots(2, len(picks), figsize=(st.DOUBLE_IN, 4.3), layout="constrained")
    for j, (cell, rpt) in enumerate(picks):
        q, v, dvdq, _, _ = C.load_discharge(C.rpt_files(DATA, cell)[rpt])
        ax[0, j].plot(q, v, color="k", lw=1.6)
        ax[1, j].plot(q, dvdq, color="k", lw=1.6)
        for k, (_, sty) in VARIANTS.items():
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
                  for lab, sty in VARIANTS.values()], loc="outside upper center", ncol=3, frameon=False)
    return st.save(fig, OUT, "bench_examples")[0]


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
    nt = t["new_nodef_notails"]
    lines.append(r"pyprobe-dma, half-cell OCPs as measured & " + f"{np.isfinite(nt['cost']).sum()} / {len(nt)}"
                 + r" & -- & -- & -- & -- \\")
    for k, (lab, _) in VARIANTS.items():
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
    for f in (plot_modes(m), plot_rmse(s), plot_examples(t, curves), plot_window(t, curves), plot_deformation(t)):
        print(f)
    latex_table(t, s, m)
    pd.set_option("display.width", 200)
    print(s.groupby("variant")[["rmse_V_mV", "rmse_dVdQ"]].median().round(3))
    w = m[m["rpt"] >= 1].pivot_table(index=["cell", "rpt"], columns="variant", values=["LAM_Si", "LLI", "LAM_Gr", "LAM_PE"])
    print(w.round(1).to_string())


if __name__ == "__main__":
    main()
