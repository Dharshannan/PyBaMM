"""Benchmark, PyProBE 2.6.1 side (environment ../environment.yml): the paper's composite DMA
(../composite_dma.py) with the Si OCP deformation (../results/dma_new.csv, already fitted) and
without it (fitted here), plus the OCV each fit implies at the two ends of the measured discharge.

    python run_old.py [--jobs 10]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import composite_dma as C  # noqa: E402

DATA = os.path.join(HERE, "..", "..", "data")
OUT = os.path.join(HERE, "..", "results", "benchmark")


def ocv_ends(hc, row, q_full):
    """Model OCV at the start (Q = 0) and the end (Q = full discharge) of the measured discharge."""
    v = hc.voltage([row[p] for p in C.PARAMS], np.array([0.0, q_full]))
    return float(v[0]), float(v[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=10)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    with_def = pd.read_csv(os.path.join(HERE, "..", "results", "dma_new.csv"))
    todo = [(r["cell"], int(r["rpt"]), np.array([r[p] for p in C.PARAMS])) for _, r in with_def.iterrows()]
    rows = Parallel(n_jobs=a.jobs, verbose=2)(
        delayed(C.fit_rpt)(DATA, cell, rpt, th, deform=False) for cell, rpt, th in todo)
    no_def = pd.DataFrame(rows)
    hc = C.HalfCells(DATA)
    for label, t in (("old_def", with_def), ("old_nodef", no_def)):
        ends = [ocv_ends(hc, r, r["measured_capacity_Ah"]) for _, r in t.iterrows()]
        t = t.assign(ocv_start_V=[e[0] for e in ends], ocv_end_V=[e[1] for e in ends], variant=label)
        t.to_csv(os.path.join(OUT, f"{label}.csv"), index=False)
        print(label, "written")


if __name__ == "__main__":
    main()
