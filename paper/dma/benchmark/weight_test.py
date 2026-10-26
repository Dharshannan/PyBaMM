"""dV/dQ weight test on one cell: refit RPT1 and a late RPT at several dV/dQ weights.

    python weight_test.py --pkg old   (PyProBE 2.6.1 environment)
    python weight_test.py --pkg new   (pyprobe-dma environment)
    python weight_test.py --report    (PyProBE 2.6.1 environment: table and figure)

Writes ../results/benchmark/wt_{old,new}.csv, wt_new_curves.npz, wt_report.csv and wt_fits.png/pdf.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
DATA = os.path.join(HERE, "..", "..", "data")
OUT = os.path.join(HERE, "..", "results", "benchmark")
CELL, RPTS = "CELL017", (1, 7)
WEIGHTS = (0.5, 2.0, 5.0, 20.0)


def _old(rpt, w, deform):
    import composite_dma as C

    C.W_DVDQ = w
    row = C.fit_rpt(DATA, CELL, rpt, deform=deform)
    return row | {"w": w, "deform": deform}


def _new(rpt, w, deform):
    import pyprobe_dma_v3 as B

    B.W_DVDQ = w
    row, curve = B.fit_rpt(DATA, CELL, rpt, deform)
    return row | {"w": w, "deform": deform}, curve


def run(pkg, jobs):
    from joblib import Parallel, delayed

    todo = [(r, w, d) for r in RPTS for w in WEIGHTS for d in (False, True)]
    if pkg == "old":
        rows = Parallel(n_jobs=jobs, verbose=2)(delayed(_old)(*t) for t in todo)
        pd.DataFrame(rows).to_csv(os.path.join(OUT, "wt_old.csv"), index=False)
    else:
        out = Parallel(n_jobs=jobs, verbose=2)(delayed(_new)(*t) for t in todo)
        pd.DataFrame([o[0] for o in out]).to_csv(os.path.join(OUT, "wt_new.csv"), index=False)
        np.savez_compressed(os.path.join(OUT, "wt_new_curves.npz"),
                            **{f"{r}_{w}_{int(d)}_{k}": a for (r, w, d), o in zip(todo, out)
                               for k, a in zip(("q", "v", "ocv"), o[1])})
    print(pkg, "done", flush=True)


def report():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    import composite_dma as C
    import run_dma as R

    st = R._style()
    hc = C.HalfCells(DATA)
    old = pd.read_csv(os.path.join(OUT, "wt_old.csv"))
    new = pd.read_csv(os.path.join(OUT, "wt_new.csv"))
    cv = dict(np.load(os.path.join(OUT, "wt_new_curves.npz")))
    data = {r: C.load_discharge(C.rpt_files(DATA, CELL)[r]) for r in RPTS}

    def model(pkg, row, q):
        if pkg == "old":
            return hc.voltage([row[p] for p in C.PARAMS], q)
        key = f"{int(row['rpt'])}_{row['w']}_{int(row['deform'])}"
        qq, vv = cv[f"{key}_q"], cv[f"{key}_v"]
        ok = np.isfinite(qq) & np.isfinite(vv)
        o = np.argsort(qq[ok])
        return np.interp(q, qq[ok][o], vv[ok][o], left=np.nan, right=np.nan)

    rows = []
    for pkg, t in (("old", old), ("new", new)):
        for (w, d), g in t.groupby(["w", "deform"]):
            ref = g[g["rpt"] == RPTS[0]].iloc[0]
            r = g[g["rpt"] == RPTS[1]].iloc[0]
            q, v, dvdq, _, _ = data[RPTS[1]]
            vm = model(pkg, r, q)
            rows.append(dict(pkg=pkg, deform=bool(d), w=w,
                             rmse_V_mV=1e3 * np.sqrt(np.nanmean((vm - v) ** 2)),
                             rmse_dVdQ=np.sqrt(np.nanmean((np.gradient(vm, q) - dvdq) ** 2)),
                             LAM_Si=100 * (1 - r["Cn_Si_Ah"] / ref["Cn_Si_Ah"]),
                             LAM_Gr=100 * (1 - r["Cn_Gr_Ah"] / ref["Cn_Gr_Ah"]),
                             LAM_PE=100 * (1 - r["Cp_Ah"] / ref["Cp_Ah"]),
                             LLI=100 * (1 - r["LLI_Ah"] / ref["LLI_Ah"]),
                             si_ocp_scale=r["si_ocp_scale"], si_ocp_shift_V=r["si_ocp_shift_V"]))
    rep = pd.DataFrame(rows)
    rep.to_csv(os.path.join(OUT, "wt_report.csv"), index=False)

    # figure: RPT7 dV/dQ, one panel per variant, one line per weight
    variants = [("old", False, "PyProBE 2.6.1, Si OCP fixed"), ("old", True, "PyProBE 2.6.1, Si OCP deformed"),
                ("new", False, "pyprobe-dma, Si OCP fixed"), ("new", True, "pyprobe-dma, Si OCP deformed")]
    cols = plt.cm.viridis(np.linspace(0, 0.85, len(WEIGHTS)))
    fig, ax = plt.subplots(2, 2, figsize=(st.DOUBLE_IN, 4.6), sharex=True, sharey=True, layout="constrained")
    q, v, dvdq, _, _ = data[RPTS[1]]
    for i, (pkg, d, lab) in enumerate(variants):
        a = ax.flat[i]
        a.plot(q, dvdq, color="k", lw=1.6)
        t = (old if pkg == "old" else new)
        for w, c in zip(WEIGHTS, cols):
            r = t[(t["rpt"] == RPTS[1]) & (t["w"] == w) & (t["deform"] == d)].iloc[0]
            a.plot(q, np.gradient(model(pkg, r, q), q), color=c, lw=1.0)
        a.set_ylim(np.percentile(dvdq, 2) - 0.3, max(0.05, dvdq.max() + 0.1))
        st.panel_label(a, i, lab)
        if i >= 2:
            a.set_xlabel("Discharge capacity [A h]")
        if i % 2 == 0:
            a.set_ylabel("dV/dQ [V/A h]")
    fig.legend(handles=[Line2D([], [], color="k", lw=1.6, label="Measured")]
               + [Line2D([], [], color=c, label=f"dV/dQ weight {w:g}") for w, c in zip(WEIGHTS, cols)],
               loc="outside upper center", ncol=5, frameon=False)
    st.save(fig, OUT, "wt_fits")
    pd.set_option("display.width", 200)
    print(rep.round(3).to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pkg", choices=("old", "new"))
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--jobs", type=int, default=16)
    a = ap.parse_args()
    report() if a.report else run(a.pkg, a.jobs)
