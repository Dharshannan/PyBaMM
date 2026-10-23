"""Refit the composite DMA of every C/20 RPT of the four cells with voltage + dV/dQ
only, and compare with the DMA we were given (base: voltage + dV/dQ + expansion).

    python run_dma.py [--data ../data] [--out results] [--jobs 10]

Writes <out>/dma_new.csv (fitted parameters), <out>/dma_modes.csv (LAM and LLI,
base and new) and the comparison figures <out>/dma_modes_compare.png and
<out>/dma_fits_<cell>.png.
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import composite_dma as C  # noqa: E402

MODES = (("LAM_Si", "LAM, Si [%]"), ("LAM_Gr", "LAM, graphite [%]"), ("LAM_PE", "LAM, positive [%]"),
         ("LLI", "LLI [%]"))
TAG = {"CELL064": "CELL064, 25 °C, full", "CELL017": "CELL017, 45 °C, full", "CELL009": "CELL009, 25 °C, 15–95%",
       "CELL026": "CELL026, 25 °C, 20–80%"}


def base_tables(data_dir):
    return {cell: pd.read_csv(os.path.join(data_dir, f, f"{cell}_dma_base.csv")) for cell, f in C.FOLDERS.items()}


def fit_all(data_dir, jobs):
    base = base_tables(data_dir)
    todo = [(cell, int(r.rpt), np.array([r[p] for p in C.PARAMS])) for cell, b in base.items()
            for _, r in b.iterrows() if int(r.rpt) in C.rpt_files(data_dir, cell)]
    rows = Parallel(n_jobs=jobs, verbose=5)(delayed(C.fit_rpt)(data_dir, cell, rpt, th) for cell, rpt, th in todo)
    return pd.DataFrame(rows)


def base_rmse(data_dir, cell, row):
    """Base fit scored on the same data, points and smoothing as the new fit."""
    hc = C.HalfCells(data_dir)
    q, v, dvdq, _, _ = C.load_discharge(C.rpt_files(data_dir, cell)[int(row["rpt"])])
    rv, rd = C.residuals(hc, [row[p] for p in C.PARAMS], q, v, dvdq)
    return float(np.sqrt(np.mean(rv ** 2))), float(np.sqrt(np.mean(rd ** 2)))


def modes_table(data_dir, new):
    out = []
    for cell, b in base_tables(data_dir).items():
        b = b[b["rpt"].isin(C.rpt_files(data_dir, cell))].copy()
        efc0 = float(b.loc[b["rpt"] == 1, "EFC"].iloc[0])
        n = new[new["cell"] == cell]
        for label, t in (("base", b), ("new", n)):
            m = C.degradation_modes(t)
            m["efc"] = m["rpt"].map(dict(zip(b["rpt"], b["EFC"] - efc0)))
            m["cell"], m["dma"] = cell, label
            out.append(m)
    return pd.concat(out, ignore_index=True)


def plot_modes(modes, path):
    cells = list(C.FOLDERS)
    fig, ax = plt.subplots(len(cells), 4, figsize=(13, 11), sharex="row", layout="constrained")
    for i, cell in enumerate(cells):
        for j, (col, lab) in enumerate(MODES):
            a = ax[i, j]
            for dma_, sty in (("base", dict(marker="s", ls="--", color="0.45", mfc="white", label="Given DMA (V + dV/dQ + expansion)")),
                              ("new", dict(marker="o", ls="-", color="tab:blue", label="New DMA (V + dV/dQ)"))):
                m = modes[(modes["cell"] == cell) & (modes["dma"] == dma_)].sort_values("efc")
                a.plot(m["efc"], m[col], ms=5, lw=1.2, **sty)
            a.axhline(0, color="k", lw=0.5)
            a.set_ylabel(lab)
            if i == len(cells) - 1:
                a.set_xlabel("EFC since RPT1 [-]")
            if j == 0:
                a.set_title(TAG[cell], loc="left", fontsize=10)
            a.grid(alpha=0.3)
    ax[0, 0].legend(fontsize=8, loc="upper left")
    fig.suptitle("Composite DMA of the C/20 RPTs: given (voltage + dV/dQ + expansion) vs new (voltage + dV/dQ only)")
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_fits(data_dir, cell, base, new, path):
    hc = C.HalfCells(data_dir)
    files = C.rpt_files(data_dir, cell)
    rpts = [r for r in sorted(files) if r in set(new["rpt"])]
    fig, ax = plt.subplots(2, len(rpts), figsize=(2.9 * len(rpts), 5.6), layout="constrained", squeeze=False)
    for k, r in enumerate(rpts):
        q, v, dvdq, _, _ = C.load_discharge(files[r])
        tb = [base.loc[base["rpt"] == r, p].iloc[0] for p in C.PARAMS]
        tn = [new.loc[(new["cell"] == cell) & (new["rpt"] == r), p].iloc[0] for p in C.PARAMS]
        vb, vn = hc.voltage(tb, q), hc.voltage(tn, q)
        ax[0, k].plot(q, v, color="k", lw=2, label="Measured (C/20)")
        ax[0, k].plot(q, vb, "--", color="0.45", lw=1.3, label="Given DMA")
        ax[0, k].plot(q, vn, color="tab:blue", lw=1.3, label="New DMA")
        ax[1, k].plot(q, dvdq, color="k", lw=2)
        ax[1, k].plot(q, np.gradient(vb, q), "--", color="0.45", lw=1.3)
        ax[1, k].plot(q, np.gradient(vn, q), color="tab:blue", lw=1.3)
        ax[1, k].set_ylim(np.percentile(dvdq, 1) - 0.3, 0.1)
        rb = 1e3 * np.sqrt(np.mean((vb - v) ** 2))
        rn = 1e3 * np.sqrt(np.mean((vn - v) ** 2))
        ax[0, k].set_title(f"RPT{r}: RMSE given {rb:.1f} / new {rn:.1f} mV", fontsize=8.5)
        ax[1, k].set_xlabel("Discharge capacity [A h]")
        for a in ax[:, k]:
            a.grid(alpha=0.3)
    ax[0, 0].set_ylabel("Voltage [V]")
    ax[1, 0].set_ylabel("dV/dQ [V/A h]")
    ax[0, 0].legend(fontsize=8)
    fig.suptitle(f"{TAG[cell]}: C/20 discharge and dV/dQ, given vs new DMA fit")
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--data", default=os.path.join(here, "..", "data"))
    ap.add_argument("--out", default=os.path.join(here, "results"))
    ap.add_argument("--jobs", type=int, default=10)
    ap.add_argument("--plots-only", action="store_true", help="reuse <out>/dma_new.csv")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, "dma_new.csv")
    new = pd.read_csv(path) if a.plots_only else fit_all(a.data, a.jobs)
    new.to_csv(path, index=False)

    base = base_tables(a.data)
    score = []
    for _, r in new.iterrows():
        b = base[r["cell"]]
        rb = b[b["rpt"] == r["rpt"]].iloc[0]
        bv, bd = base_rmse(a.data, r["cell"], rb)
        score.append(dict(cell=r["cell"], rpt=int(r["rpt"]), rmse_V_mV_given=1e3 * bv, rmse_V_mV_new=1e3 * r["rmse_voltage_V"],
                          rmse_dVdQ_given=bd, rmse_dVdQ_new=r["rmse_dVdQ_V_per_Ah"],
                          Cn_Si_given=rb["Cn_Si_Ah"], Cn_Si_new=r["Cn_Si_Ah"], Cn_Gr_given=rb["Cn_Gr_Ah"],
                          Cn_Gr_new=r["Cn_Gr_Ah"], Cp_given=rb["Cp_Ah"], Cp_new=r["Cp_Ah"],
                          LLI_given=rb["LLI_Ah"], LLI_new=r["LLI_Ah"]))
    score = pd.DataFrame(score)
    score.to_csv(os.path.join(a.out, "dma_fit_scores.csv"), index=False)
    modes = modes_table(a.data, new)
    modes.to_csv(os.path.join(a.out, "dma_modes.csv"), index=False)
    plot_modes(modes, os.path.join(a.out, "dma_modes_compare.png"))
    for cell in C.FOLDERS:
        plot_fits(a.data, cell, base[cell], new, os.path.join(a.out, f"dma_fits_{cell}.png"))
    pd.set_option("display.width", 200)
    print(score.round(3).to_string(index=False))
    wide = modes.pivot_table(index=["cell", "rpt"], columns="dma", values=["LAM_Si", "LAM_Gr", "LAM_PE", "LLI"])
    print(wide.round(1).to_string())


if __name__ == "__main__":
    main()
