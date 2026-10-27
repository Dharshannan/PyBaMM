"""PyProBE 2.6.1 CELL017 RPT7 with the Si OCP fixed, refitted with the collaborators' cathode box
(2.6-2.8 Ah), as the deformed fit; patches old_nodef.csv."""
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
import composite_dma as C  # noqa: E402
from run_old import ocv_ends  # noqa: E402

DATA = os.path.join(HERE, "..", "..", "data")
OUT = os.path.join(HERE, "..", "results", "benchmark")
b = list(C.BOUNDS)
b[2] = (2.6, 2.8)
C.BOUNDS = tuple(b)
row = C.fit_rpt(DATA, "CELL017", 7, deform=False)
print({k: (round(float(v), 4) if isinstance(v, float) else v) for k, v in row.items()})
t = pd.read_csv(os.path.join(OUT, "old_nodef.csv"))
i = t.index[(t["cell"] == "CELL017") & (t["rpt"] == 7)][0]
for k, v in row.items():
    if k in t.columns and k not in ("cell", "rpt"):
        t.loc[i, k] = v
t.loc[i, "ocv_start_V"], t.loc[i, "ocv_end_V"] = ocv_ends(C.HalfCells(DATA), t.loc[i], t.loc[i, "measured_capacity_Ah"])
t.to_csv(os.path.join(OUT, "old_nodef.csv"), index=False)
print("patched")
