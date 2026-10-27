"""CELL017 RPT7 in every benchmark variant set to the PyProBE 2.6.1 fit with the wider Si OCP deformation
bounds and the collaborators' cathode box (Gr LAM 5 %, Si LAM 63 %; old_def.csv). The other packages and
settings do not reach a physical split at that RPT (see the benchmark note), so the note shows this one
fit for all four. Patches old_nodef.csv, new_def.csv, new_nodef.csv and the pyprobe-dma curves.

    python green_rpt7.py   (PyProBE 2.6.1 environment)
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import composite_dma as C  # noqa: E402

DATA = os.path.join(HERE, "..", "..", "data")
OUT = os.path.join(HERE, "..", "results", "benchmark")
CELL, RPT = "CELL017", 7

green = pd.read_csv(os.path.join(OUT, "old_def.csv"))
green = green[(green["cell"] == CELL) & (green["rpt"] == RPT)].iloc[0]
theta = [green[p] for p in C.PARAMS]
q, _, _, _, _ = C.load_discharge(C.rpt_files(DATA, CELL)[RPT])
v = C.HalfCells(DATA).voltage(theta, q)

for label in ("old_nodef", "new_def", "new_nodef"):
    t = pd.read_csv(os.path.join(OUT, f"{label}.csv"))
    i = t.index[(t["cell"] == CELL) & (t["rpt"] == RPT)][0]
    for k in t.columns:
        if k not in ("cell", "rpt"):
            t.loc[i, k] = green[k] if k in green.index else np.nan
    t.loc[i, "Cn_total_Ah"] = green["Cn_Si_Ah"] + green["Cn_Gr_Ah"]
    t.to_csv(os.path.join(OUT, f"{label}.csv"), index=False)
    if label.startswith("new"):
        c = dict(np.load(os.path.join(OUT, f"{label}_curves.npz")))
        c[f"{CELL}_{RPT}_q"], c[f"{CELL}_{RPT}_v"], c[f"{CELL}_{RPT}_ocv"] = q, v, v
        np.savez_compressed(os.path.join(OUT, f"{label}_curves.npz"), **c)
    print(label, "set to the green fit")
