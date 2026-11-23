"""Measured capacity knee of each paper cell: the breakpoint of a continuous two-segment (hinge) least-squares fit
y = b + s1*EFC + ds*max(0, EFC - knee) to the measured RPT SoH from RPT1 on (C/20 RPTs, plus CELL064's high-rate
RPTs). The knee is searched between the second and the second-to-last RPT, so that each segment rests on at least
two RPTs. The model knee of the paper uses the same fit on the model's C/3 SoH.

Prints the knee (EFC since RPT1) and the raw test EFC to store as Condition.knee_efc_raw in sigr/conditions.py.

Run from paper/:  python tools/measured_knee.py [data_dir]
"""
import sys

import numpy as np

sys.path.insert(0, ".")
from sigr import data  # noqa: E402
from sigr.conditions import CONDITIONS  # noqa: E402


def hinge_knee(x, y, lo, hi, n=20001):
    best = (np.inf, np.nan)
    for xb in np.linspace(lo, hi, n):
        a = np.column_stack([np.ones_like(x), x, np.clip(x - xb, 0, None)])
        c, *_ = np.linalg.lstsq(a, y, rcond=None)
        r = float(np.sum((a @ c - y) ** 2))
        if r < best[0]:
            best = (r, xb)
    return best[1]


def measured_knee(cond, data_dir):
    cap = data.load(cond, data_dir, "capacity_fade")  # EFC since RPT1, earlier rows dropped
    x, y = cap["efc"].to_numpy(float), cap["capacity_retention_pct"].to_numpy(float)
    return hinge_knee(x, y, x[1], x[-2])


if __name__ == "__main__":
    dd = sys.argv[1] if len(sys.argv) > 1 else "data"
    for cond in CONDITIONS.values():
        k = measured_knee(cond, dd)
        print(f"{cond.cell}: knee {k:.1f} EFC since RPT1; knee_efc_raw = {k + data.rpt1_efc(cond, dd):.6f}")
