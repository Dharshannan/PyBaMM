"""Undo the CELL017 RPT7 refits of the Si-fixed PyProBE 2.6.1 fit and the deformed pyprobe-dma fit:
the former is refitted exactly as run_old.py did (default bounds, started from the original deformed fit),
the latter is restored from the weight test at the default weight (wt_new.csv, same procedure)."""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
import composite_dma as C  # noqa: E402
from run_old import ocv_ends  # noqa: E402

DATA = os.path.join(HERE, "..", "..", "data")
OUT = os.path.join(HERE, "..", "results", "benchmark")
# original deformed fit of CELL017 RPT7 (before the wider-bounds refit), the start run_old.py used
th0 = np.array([0.4041, 1.7771, 2.3802, 0.6426, 0.0, 0.8473, 0.03])
row = C.fit_rpt(DATA, "CELL017", 7, th0, deform=False)
t = pd.read_csv(os.path.join(OUT, "old_nodef.csv"))
i = t.index[(t["cell"] == "CELL017") & (t["rpt"] == 7)][0]
for k, v in row.items():
    if k in t.columns and k not in ("cell", "rpt"):
        t.loc[i, k] = v
t.loc[i, "ocv_start_V"], t.loc[i, "ocv_end_V"] = ocv_ends(C.HalfCells(DATA), t.loc[i], t.loc[i, "measured_capacity_Ah"])
t.to_csv(os.path.join(OUT, "old_nodef.csv"), index=False)
print("old_nodef RPT7:", {p: round(float(row[p]), 4) for p in C.PARAMS}, "rmse", round(1e3 * row["rmse_voltage_V"], 2))

w = pd.read_csv(os.path.join(OUT, "wt_new.csv"))
r = w[(w["rpt"] == 7) & (w["w"] == 0.5) & (w["deform"])].iloc[0]
t = pd.read_csv(os.path.join(OUT, "new_def.csv"))
i = t.index[(t["cell"] == "CELL017") & (t["rpt"] == 7)][0]
for k in t.columns:
    if k in r.index and k not in ("cell", "rpt", "variant"):
        t.loc[i, k] = r[k]
t.to_csv(os.path.join(OUT, "new_def.csv"), index=False)
cw = dict(np.load(os.path.join(OUT, "wt_new_curves.npz")))
c = dict(np.load(os.path.join(OUT, "new_def_curves.npz")))
for k in ("q", "v", "ocv"):
    c[f"CELL017_7_{k}"] = cw[f"7_0.5_1_{k}"]
np.savez_compressed(os.path.join(OUT, "new_def_curves.npz"), **c)
print("new_def RPT7 restored: scale %.3f shift %.3f Cn_Si %.3f Cn_Gr %.3f" % (r.si_ocp_scale, r.si_ocp_shift_V, r.Cn_Si_Ah, r.Cn_Gr_Ah))
