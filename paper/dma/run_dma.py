"""Refit the composite DMA of every C/20 RPT of the four cells with voltage + dV/dQ
only, and compare with the DMA we were given (base: voltage + dV/dQ + expansion).

    python run_dma.py [--data ../data] [--out results] [--jobs 10]

Writes <out>/dma_new.csv (fitted parameters), <out>/dma_modes.csv (LAM and LLI,
given and new), <out>/dma_fit_scores.csv, the comparison figures (JES style, PNG
and PDF) <out>/dma_modes_compare.* and <out>/dma_fits_<cell>.*, and the LaTeX
table <out>/dma_table.tex. Then `python tables.py` writes the paper's LAM/LLI
tables from dma_new.csv.
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


def _style():
    """The paper's JES figure style (sigr/style.py), loaded by path: importing the
    sigr package would pull in PyBaMM, which this environment does not have."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "sigr_style", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "sigr", "style.py"))
    st = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(st)
    st.use()
    return st


GIVEN = dict(marker="s", ls="--", color="0.45", mfc="white", mew=1.0, ms=4, lw=1.0)
NEW = dict(marker="o", ls="-", color="tab:blue", ms=3.8, lw=1.1)


def plot_modes(modes, out_dir):
    """LAM (Si, graphite, positive) and LLI against EFC: one row per cell, given vs new DMA."""
    from matplotlib.lines import Line2D

    st = _style()
    cells = list(C.FOLDERS)
    fig = plt.figure(figsize=(st.DOUBLE_IN, 7.6), layout="constrained")
    rows = fig.subfigures(len(cells), 1)
    k = 0
    for i, cell in enumerate(cells):
        ax = rows[i].subplots(1, 4)
        rows[i].suptitle(TAG[cell], x=0.01, ha="left", fontsize=8, fontweight="bold")
        for j, (col, lab) in enumerate(MODES):
            a = ax[j]
            for dma_, sty in (("base", GIVEN), ("new", NEW)):
                m = modes[(modes["cell"] == cell) & (modes["dma"] == dma_)].sort_values("efc")
                a.plot(m["efc"], m[col], **sty)
            a.axhline(0, color="k", lw=0.5)
            a.set_ylabel(lab)
            if i == len(cells) - 1:
                a.set_xlabel("EFC since RPT1 [-]")
            st.panel_label(a, k)
            k += 1
    fig.legend(handles=[Line2D([], [], label="Given DMA (voltage + dV/dQ + expansion)", **GIVEN),
                        Line2D([], [], label="New DMA, PyProBE (voltage + dV/dQ)", **NEW)],
               loc="outside upper center", ncol=2, frameon=False)
    return [st.save(fig, out_dir, "dma_modes_compare")[0]]


def plot_fits(data_dir, cell, base, new, out_dir, collab=None):
    """Each C/20 RPT: voltage (upper panel) and dV/dQ (lower panel), measured vs given
    and new DMA fits (and the collaborators' discharge fits, when given); four RPTs per
    block of two rows."""
    from matplotlib.lines import Line2D

    st = _style()
    hc = C.HalfCells(data_dir)
    files = C.rpt_files(data_dir, cell)
    rpts = [r for r in sorted(files) if r in set(new.loc[new["cell"] == cell, "rpt"])]
    ncol = 4
    nblk = int(np.ceil(len(rpts) / ncol))
    fig, ax = plt.subplots(2 * nblk, ncol, figsize=(st.DOUBLE_IN, 2.6 * nblk + 0.35), layout="constrained",
                           squeeze=False)
    k = 0
    for n, r in enumerate(rpts):
        av, ad = ax[2 * (n // ncol), n % ncol], ax[2 * (n // ncol) + 1, n % ncol]
        q, v, dvdq, _, _ = C.load_discharge(files[r])
        tb = [base.loc[base["rpt"] == r, p].iloc[0] for p in C.PARAMS]
        tn = [new.loc[(new["cell"] == cell) & (new["rpt"] == r), p].iloc[0] for p in C.PARAMS]
        fits = [(hc.voltage(tb, q), dict(ls="--", color="0.45")), (hc.voltage(tn, q), dict(color="tab:blue"))]
        if collab is not None:
            tc = [collab.loc[(collab["cell"] == cell) & (collab["rpt"] == r), p].iloc[0] for p in C.PARAMS]
            fits.append((hc.voltage(tc, q), dict(ls=":", color="tab:red")))
        av.plot(q, v, color="k", lw=1.6)
        ad.plot(q, dvdq, color="k", lw=1.6)
        for vm, sty in fits:
            av.plot(q, vm, lw=1.0, **sty)
            ad.plot(q, np.gradient(vm, q), lw=1.0, **sty)
        ad.set_ylim(np.percentile(dvdq, 2) - 0.25, max(0.05, dvdq.max() + 0.1))
        rmse = " / ".join(f"{1e3 * np.sqrt(np.mean((vm - v) ** 2)):.1f}" for vm, _ in fits)
        av.text(0.03, 0.04, f"RMSE {rmse} mV", transform=av.transAxes, ha="left", va="bottom", fontsize=7)
        st.panel_label(av, k, f"RPT{r}")
        st.panel_label(ad, k + 1)
        k += 2
        ad.set_xlabel("Discharge capacity [A h]")
    for a in list(ax.flat):
        if not a.lines:
            a.set_visible(False)
    for b in range(nblk):
        ax[2 * b, 0].set_ylabel("Voltage [V]")
        ax[2 * b + 1, 0].set_ylabel("dV/dQ [V/A h]")
    handles = [Line2D([], [], color="k", lw=1.6, label="Measured, C/20 discharge"),
               Line2D([], [], color="0.45", ls="--", label="Given DMA"),
               Line2D([], [], color="tab:blue", label="New DMA (PyProBE)")]
    if collab is not None:
        handles.append(Line2D([], [], color="tab:red", ls=":", label="Collaborators' DMA (no expansion)"))
    fig.legend(handles=handles, loc="outside upper center", ncol=len(handles), frameon=False)
    return [st.save(fig, out_dir, f"dma_fits_{cell}")[0]]


def latex_table(score, modes, path):
    """Per-RPT table for the note: RMSEs and the degradation modes, given / new."""
    m = modes.pivot_table(index=["cell", "rpt"], columns="dma", values=["LAM_Si", "LAM_Gr", "LAM_PE", "LLI"])
    lines = []
    for cell in C.FOLDERS:
        sc = score[score["cell"] == cell].sort_values("rpt")
        for i, r in enumerate(sc.itertuples()):
            g = m.loc[(cell, r.rpt)]
            name = TAG[cell].replace("%", r"\%").replace("°C", r"$^\circ$C").replace("–", "--") if i == 0 else ""
            lines.append(f"{name} & {r.rpt} & {r.rmse_V_mV_given:.1f} / {r.rmse_V_mV_new:.1f} & "
                         f"{g[('LAM_Si', 'base')]:.0f} / {g[('LAM_Si', 'new')]:.0f} & "
                         f"{g[('LAM_Gr', 'base')]:.0f} / {g[('LAM_Gr', 'new')]:.0f} & "
                         f"{g[('LAM_PE', 'base')]:.0f} / {g[('LAM_PE', 'new')]:.0f} & "
                         f"{g[('LLI', 'base')]:.0f} / {g[('LLI', 'new')]:.0f} " + r"\\")
        lines.append(r"\midrule")
    lines[-1] = r"\bottomrule"
    head = (r"\begin{tabular}{lrrrrrr}" + "\n" + r"\toprule" + "\n"
            + r"Cell & RPT & V RMSE [mV] & LAM Si [\%] & LAM Gr [\%] & LAM PE [\%] & LLI [\%] \\" + "\n"
            + r"\midrule" + "\n")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(head + "\n".join(lines) + "\n" + r"\end{tabular}" + "\n")
    return path


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
    plot_modes(modes, a.out)
    for cell in C.FOLDERS:
        plot_fits(a.data, cell, base[cell], new, a.out)
    latex_table(score, modes, os.path.join(a.out, "dma_table.tex"))
    pd.set_option("display.width", 200)
    print(score.round(3).to_string(index=False))
    wide = modes.pivot_table(index=["cell", "rpt"], columns="dma", values=["LAM_Si", "LAM_Gr", "LAM_PE", "LLI"])
    print(wide.round(1).to_string())


if __name__ == "__main__":
    main()
