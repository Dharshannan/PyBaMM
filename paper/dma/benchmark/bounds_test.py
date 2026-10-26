"""CELL017 RPT1/RPT7 refitted (PyProBE 2.6.1, Si OCP deformed) with the collaborators' bounds:
wider Si OCP deformation and/or their cathode-capacity box. Compared with their own fit.

    python bounds_test.py
"""
import os
import sys

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
DATA = os.path.join(HERE, "..", "..", "data")
OUT = os.path.join(HERE, "..", "results", "benchmark")
CELL = "CELL017"
# PARAMS order: Cn_Si, Cn_Gr, Cp, x_start, y_start, si_ocp_scale, si_ocp_shift_V
VARIANTS = {
    "ours": None,
    "wide_deform": dict(scale=(0.5, 1.1), shift=(-0.10, 0.08)),
    "wide_deform+Cp": dict(scale=(0.5, 1.1), shift=(-0.10, 0.08), cp=(2.6, 2.8)),
}


def one(name, rpt):
    import composite_dma as C

    v = VARIANTS[name]
    if v:
        b = list(C.BOUNDS)
        b[5], b[6] = v["scale"], v["shift"]
        if "cp" in v:
            b[2] = v["cp"]
        C.BOUNDS = tuple(b)
    return C.fit_rpt(DATA, CELL, rpt, deform=True) | {"variant": name}


def main():
    import composite_dma as C

    rows = Parallel(n_jobs=6)(delayed(one)(n, r) for n in VARIANTS for r in (1, 7))
    t = pd.DataFrame(rows)
    col = pd.read_csv(os.path.join(DATA, "dma_collaborators_discharge.csv"))
    col = col[col["cell"] == CELL].assign(variant="collaborators")
    t = pd.concat([t, col[col["rpt"].isin((1, 7))]], ignore_index=True)
    hc = C.HalfCells(DATA)
    q, v, dvdq, _, _ = C.load_discharge(C.rpt_files(DATA, CELL)[7])
    out = []
    for name, g in t.groupby("variant", sort=False):
        ref, r = g[g["rpt"] == 1].iloc[0], g[g["rpt"] == 7].iloc[0]
        vm = hc.voltage([r[p] for p in C.PARAMS], q)
        out.append(dict(variant=name, rmse_V_mV=1e3 * np.sqrt(np.mean((vm - v) ** 2)),
                        rmse_dVdQ=np.sqrt(np.mean((np.gradient(vm, q) - dvdq) ** 2)),
                        LAM_Si=100 * (1 - r["Cn_Si_Ah"] / ref["Cn_Si_Ah"]),
                        LAM_Gr=100 * (1 - r["Cn_Gr_Ah"] / ref["Cn_Gr_Ah"]),
                        LAM_PE=100 * (1 - r["Cp_Ah"] / ref["Cp_Ah"]), LLI=100 * (1 - r["LLI_Ah"] / ref["LLI_Ah"]),
                        Cp=r["Cp_Ah"], scale=r["si_ocp_scale"], shift=r["si_ocp_shift_V"]))
    rep = pd.DataFrame(out)
    t.to_csv(os.path.join(OUT, "bounds_fits.csv"), index=False)
    rep.to_csv(os.path.join(OUT, "bounds_report.csv"), index=False)
    pd.set_option("display.width", 200)
    print(rep.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
