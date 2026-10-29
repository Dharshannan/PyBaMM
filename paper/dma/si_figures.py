"""Supplementary figures of the paper's DMA (ours: pyprobe-dma, trimmed C/20 discharge, Si OCP deformed;
results/dma_pyprobe_dma.csv, fitted curves in results/benchmark/new_def_curves.npz), one per cell: LAM and LLI
against EFC, and the voltage and dV/dQ of each C/20 RPT from RPT1 on, measured against the fit, with the RMSEs.

    python si_figures.py [--data ../data] [--out ../paper_pics/si]

Writes si_dma_<cell>.png|pdf and prints the RMSE ranges. The LAM and LLI are the tables the paper reads
(<data>/<folder>/<CELL>_LAM_derived.csv and _LLI.csv, written by tables.py --source ours); the RMSEs are on the
fitted part of each discharge (composite_dma.load_discharge), as in compare_ours.scores.
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import composite_dma as C  # noqa: E402
import run_dma as R  # noqa: E402
from compare_ours import CURVES, ours_voltage  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIT = dict(color="tab:orange", lw=1.1, ls="--")
MODE_COLOURS = (("LAM_NE_silicon_pct", "Si", "darkorange"), ("LAM_NE_graphite_pct", "graphite", "saddlebrown"),
                ("LAM_PE_pct", "positive", "seagreen"))


def modes(data_dir, cell):
    """LAM (Si, graphite, positive) and LLI [%] against EFC since RPT1, from the paper's tables."""
    d = os.path.join(data_dir, C.FOLDERS[cell])
    lam = pd.read_csv(os.path.join(d, f"{cell}_LAM_derived.csv"))
    lli = pd.read_csv(os.path.join(d, f"{cell}_LLI.csv"))
    efc0 = float(lam.loc[lam["rpt"] == 1, "efc"].iloc[0])
    lam["efc"] -= efc0
    lli["efc"] -= efc0
    lli["LLI_pct"] = 100 * lli["charge_LLI_loss"] / lli.loc[lli["rpt"] == 1, "charge_LLI"].iloc[0]
    return lam[lam["efc"] >= 0], lli[lli["efc"] >= 0]


def figure(data_dir, cell, curves):
    """(figure, [(rpt, voltage RMSE [mV], dV/dQ RMSE [V/(A h)])])."""
    st = R._style()
    files = {r: p for r, p in C.rpt_files(data_dir, cell).items() if r >= 1}
    rpts = sorted(files)
    ncol = 4
    nblk = int(np.ceil(len(rpts) / ncol))
    fig = plt.figure(figsize=(st.DOUBLE_IN, 1.9 + 2.55 * nblk + 0.3), layout="constrained")
    top, bottom = fig.subfigures(2, 1, height_ratios=[1.9, 2.55 * nblk])
    am = top.subplots(1, 2)
    lam, lli = modes(data_dir, cell)
    for col, lab, c in MODE_COLOURS:
        am[0].plot(lam["efc"], lam[col], "s-", ms=3.5, lw=1.0, color=c, label=f"LAM {lab}")
    am[1].plot(lli["efc"], lli["LLI_pct"], "D-", ms=3.2, lw=1.0, color="crimson")
    am[0].set_ylabel("LAM [%]")
    am[1].set_ylabel("LLI [%]")
    am[0].legend(loc="upper left")
    for a in am:
        a.axhline(0, color="k", lw=0.5)
        a.set_xlabel("EFC since RPT1 [-]")
    st.panel_label(am[0], 0)
    st.panel_label(am[1], 1)
    ax = bottom.subplots(2 * nblk, ncol, squeeze=False)
    k, scores = 2, []
    for n, r in enumerate(rpts):
        av, ad = ax[2 * (n // ncol), n % ncol], ax[2 * (n // ncol) + 1, n % ncol]
        q, v, dvdq, _, _ = C.load_discharge(files[r])
        vm = ours_voltage(curves, cell, r, q)
        dm = np.gradient(vm, q)
        rv = 1e3 * np.sqrt(np.nanmean((vm - v) ** 2))
        rd = np.sqrt(np.nanmean((dm - dvdq) ** 2))
        scores.append((r, rv, rd))
        av.plot(q, v, color="k", lw=1.5)
        av.plot(q, vm, **FIT)
        ad.plot(q, dvdq, color="k", lw=1.5)
        ad.plot(q, dm, **FIT)
        ad.set_ylim(np.percentile(dvdq, 2) - 0.25, max(0.05, dvdq.max() + 0.1))
        av.text(0.97, 0.95, f"RMSE {rv:.1f} mV", transform=av.transAxes, ha="right", va="top", fontsize=7)
        ad.text(0.5, 0.05, f"RMSE {rd:.3f} V A$^{{-1}}$ h$^{{-1}}$", transform=ad.transAxes, ha="center",
                va="bottom", fontsize=7)
        st.panel_label(av, k, f"RPT{r}")
        st.panel_label(ad, k + 1)
        k += 2
        ad.set_xlabel("Discharge capacity [A h]")
    for a in ax.flat:
        if not a.lines:
            a.set_visible(False)
    for b in range(nblk):
        ax[2 * b, 0].set_ylabel("Voltage [V]")
        ax[2 * b + 1, 0].set_ylabel("dV/dQ [V A$^{-1}$ h$^{-1}$]")
    bottom.legend(handles=[Line2D([], [], color="k", lw=1.5, label="Measured, C/20 discharge"),
                           Line2D([], [], label="DMA fit", **FIT)], loc="outside upper center", ncol=2, frameon=False)
    return fig, scores


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(HERE, "..", "data"))
    ap.add_argument("--out", default=os.path.join(HERE, "..", "paper_pics", "si"))
    a = ap.parse_args()
    curves = np.load(CURVES)
    st = R._style()
    rows = []
    for cell in C.FOLDERS:
        fig, sc = figure(a.data, cell, curves)
        print(st.save(fig, a.out, f"si_dma_{cell}")[0])
        plt.close(fig)
        rows += [dict(cell=cell, rpt=r, rmse_V_mV=rv, rmse_dVdQ=rd) for r, rv, rd in sc]
    df = pd.DataFrame(rows)
    print(df.round(3).to_string(index=False))
    print(f"voltage RMSE {df.rmse_V_mV.min():.1f}-{df.rmse_V_mV.max():.1f} mV; "
          f"dV/dQ RMSE {df.rmse_dVdQ.min():.3f}-{df.rmse_dVdQ.max():.3f} V/(A h)")


if __name__ == "__main__":
    main()
