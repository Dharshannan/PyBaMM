"""pyprobe-dma CELL017 RPT7 refitted with the collaborators' cathode box (2.6-2.8 Ah), as the PyProBE 2.6.1
fits adopted for it. The cathode is made the free electrode so that its capacity is a direct box bound (a
bound on a derived quantity fails in pyprobe-dma 0.1.0). With the Si OCP deformed, the deformation is held
at a = 0.80, b = +80 mV, the wider-bounds optimum without the box: a free deformation search under the box
ends at a = 0.92, b = -20 mV with 40 % graphite LAM and misses the dV/dQ dip at ~0.2 Ah. Patches
new_def.csv, new_nodef.csv and their curves.

    python refit_rpt7_cpbox_v3.py
"""
import os
import sys

import numpy as np
import pandas as pd
import pyprobe_dma as dma

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
import pyprobe_dma_v3 as B  # noqa: E402

DATA = os.path.join(HERE, "..", "..", "data")
OUT = os.path.join(HERE, "..", "results", "benchmark")
CP_BOX = (2.6, 2.8)
DEFORM = (0.80, 0.08)


def ocv_model(self, scale=1.0, shift=0.0):
    ne = dma.CompositeOCP.from_arrays(ocp_comp1=self.u_gr, sto_comp1=self.x_gr,
                                      ocp_comp2=scale * self.u_si + shift, sto_comp2=self.x_si)
    return dma.OCVModel(ocp_pe=self.pe, ocp_ne=ne, voltage_limits=B.VOLTAGE_LIMITS,
                        soc_estimator=dma.TwoPointSOCEstimator(), free_electrode="pe",
                        resistance_model=dma.ResistanceModel.from_template(template="constant"))


def one(deform):
    """Best-of-seeds fit with the cathode box, set in this (worker) process."""
    sys.path.insert(0, os.path.join(HERE, ".."))
    import composite_dma as C

    B.HalfCells.ocv_model = ocv_model
    B.constraints = lambda: [dma.Constraint("pe_capacity", lo=CP_BOX[0], hi=CP_BOX[1]),
                             dma.Constraint("r_params", lo=[0.0], hi=[0.3])]
    scale, shift = DEFORM if deform else (1.0, 0.0)
    hc = B.HalfCells(DATA)
    vc, ic = B.curves(C.rpt_files(DATA, "CELL017")[7])
    res = B.fit_best(hc, vc, ic, scale, shift)
    return B.summarise(res, "CELL017", 7, scale, shift), B.fit_curve(res)


def patch(label, row, curve):
    t = pd.read_csv(os.path.join(OUT, f"{label}.csv"))
    i = t.index[(t["cell"] == "CELL017") & (t["rpt"] == 7)][0]
    for k, v in row.items():
        if k in t.columns:
            t.loc[i, k] = v
    t.to_csv(os.path.join(OUT, f"{label}.csv"), index=False)
    c = dict(np.load(os.path.join(OUT, f"{label}_curves.npz")))
    for k, a in zip(("q", "v", "ocv"), curve):
        c[f"CELL017_7_{k}"] = a
    np.savez_compressed(os.path.join(OUT, f"{label}_curves.npz"), **c)


if __name__ == "__main__":
    from joblib import Parallel, delayed

    out = Parallel(n_jobs=2)(delayed(one)(d) for d in (False, True))
    for label, (row, curve) in zip(("new_nodef", "new_def"), out):
        print(label, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
        if "--dry" not in sys.argv:
            patch(label, row, curve)
