"""pyprobe-dma (trimmed, Si OCP deformed) fits of CELL064's high-rate neighbour discharges: the C/3 cycling
discharges just before and after each RPT (the RPT0 and RPT3 capacities in the data bundle are their means,
as no C/20 discharge was recorded there). Same setting as new_def.csv.

    python highrate_v3.py <high_rate_neighbor dir> <CELL064_capacity_fade_highrate_all.csv>
    python -c "import highrate_v3; highrate_v3.table(<data dir>)"   (the plotted LAM/LLI table)

Writes ../results/benchmark/highrate_new_def.csv (one row per discharge, with side and raw EFC).
"""
import os
import re
import sys

import pandas as pd
from joblib import Parallel, delayed

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
DATA = os.path.join(HERE, "..", "..", "data")
OUT = os.path.join(HERE, "..", "results", "benchmark")
CELL = "CELL064"


def one(path, rpt, side):
    import pyprobe_dma_v3 as B

    row, _ = B.fit_rpt(DATA, CELL, rpt, True, path=path)
    return row | {"side": side, "file": os.path.basename(path)}


def table(data_dir, rpts=(3,)):
    """<CELL>_DMA_highrate.csv in data_dir's CELL064 folder: LAM and LLI of the substitute RPTs (the mean of
    the fits before and after), relative to the mean of the RPT1 high-rate fits, so that the rate's bias
    cancels; EFC as the bundle's substitute capacity point."""
    t = pd.read_csv(os.path.join(OUT, "highrate_new_def.csv")).groupby("rpt")[
        ["Cn_Si_Ah", "Cn_Gr_Ah", "Cp_Ah", "LLI_Ah"]].mean()
    folder = os.path.join(data_dir, "CELL064_25C_full")
    efc = pd.read_csv(os.path.join(folder, f"{CELL}_capacity_fade.csv")).set_index("rpt")["efc"]
    ref, f = t.loc[1], t.loc[list(rpts)]
    out = pd.DataFrame({"rpt": f.index, "efc": efc.loc[f.index].to_numpy(),
                        "LAM_NE_graphite_pct": 100 * (1 - f["Cn_Gr_Ah"] / ref["Cn_Gr_Ah"]).to_numpy(),
                        "LAM_NE_silicon_pct": 100 * (1 - f["Cn_Si_Ah"] / ref["Cn_Si_Ah"]).to_numpy(),
                        "LAM_PE_pct": 100 * (1 - f["Cp_Ah"] / ref["Cp_Ah"]).to_numpy(),
                        "LLI_pct": 100 * (1 - f["LLI_Ah"] / ref["LLI_Ah"]).to_numpy()})
    out.to_csv(os.path.join(folder, f"{CELL}_DMA_highrate.csv"), index=False)
    return out


if __name__ == "__main__":
    src, cap = sys.argv[1], sys.argv[2]
    todo = []
    for f in sorted(os.listdir(src)):
        m = re.match(rf"{CELL}_RPT(\d{{3}})_(before|after)_CYC\d+_highrate_discharge\.csv$", f)
        if m:
            todo.append((os.path.join(src, f), int(m.group(1)), m.group(2)))
    rows = Parallel(n_jobs=len(todo), verbose=2)(delayed(one)(*t) for t in todo)
    t = pd.DataFrame(rows).merge(pd.read_csv(cap)[["rpt", "side", "efc"]], on=["rpt", "side"], how="left")
    t.to_csv(os.path.join(OUT, "highrate_new_def.csv"), index=False)
    print(t[["rpt", "side", "efc", "Cn_Si_Ah", "Cn_Gr_Ah", "Cp_Ah", "LLI_Ah", "si_ocp_scale", "si_ocp_shift_V",
             "cost"]].round(4).to_string(index=False))
