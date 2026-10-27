"""pyprobe-dma CELL017 RPT7 (Si OCP deformed) refitted with the wider deformation bounds used for the
PyProBE 2.6.1 fit (scale 0.5-1.1, shift -0.10-0.08 V); patches new_def.csv and its curves.
(No cathode-capacity box: pyprobe-dma 0.1.0 fails on a pe_capacity constraint with a resistance model.)"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
import pyprobe_dma_v3 as B  # noqa: E402

OUT = os.path.join(HERE, "..", "results", "benchmark")
B.DEFORM_BOUNDS = ((0.5, 1.1), (-0.10, 0.08))
B.SCALES = (0.6, 0.7, 0.8, 0.9, 1.0)
B.SHIFTS = (-0.06, -0.02, 0.02, 0.05, 0.08)
row, curve = B.fit_rpt(os.path.join(HERE, "..", "..", "data"), "CELL017", 7, True)
print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()})
t = pd.read_csv(os.path.join(OUT, "new_def.csv"))
i = t.index[(t["cell"] == "CELL017") & (t["rpt"] == 7)][0]
for k, v in row.items():
    if k in t.columns:
        t.loc[i, k] = v
t.to_csv(os.path.join(OUT, "new_def.csv"), index=False)
c = dict(np.load(os.path.join(OUT, "new_def_curves.npz")))
for k, a in zip(("q", "v", "ocv"), curve):
    c[f"CELL017_7_{k}"] = a
np.savez_compressed(os.path.join(OUT, "new_def_curves.npz"), **c)
print("patched")
