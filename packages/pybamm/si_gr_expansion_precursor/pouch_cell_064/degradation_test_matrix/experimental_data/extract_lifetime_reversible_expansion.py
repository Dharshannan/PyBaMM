# -*- coding: utf-8 -*-
"""Per-cycle reversible expansion amplitude from an IC_PYBAMM_FIG1 data-package
lifetime pickle (cells/CELLxxx/lifetime/*_cleaned.pkl.gz).

Amplitude per ageing cycle = max - min of `expansion_cleaned_from_bol_um`
over that cycle's discharge half (full charge -> end of discharge), i.e. the
reversible breathing swing. Only ageing-rate discharges are kept: low-rate
(C/20 RPT) discharges and short/partial segments are dropped, as are
segments spanning an acquisition-session boundary (the package stitches
steps only within a session).

Validated against CELL064's own pipeline output (CELL064_reversible_expansion.csv,
segment_stitch_amp_um) via --validate.
"""
import argparse
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


def per_cycle_amplitude(pkl_path, min_abs_current=0.3, min_segment_ah=0.5):
    df = pd.read_pickle(pkl_path)
    I = df["current_A"].to_numpy()
    # Discharge = negative current in this package's sign convention (C/20 RPT
    # discharges sit at ~-0.125 A).
    dis = I < -0.02
    # Segment boundaries: discharge on/off transitions or session changes.
    sess = df["vdf_session"].cat.codes.to_numpy()
    edge = np.flatnonzero(np.diff(dis.astype(np.int8)) != 0) + 1
    starts = np.r_[0, edge]
    ends = np.r_[edge, len(df)]
    q = df["capacity_Ah_cumulative"].to_numpy()
    efc = df["EFC"].to_numpy()
    exp_um = df["expansion_cleaned_from_bol_um"].to_numpy()
    rows = []
    for s, e in zip(starts, ends):
        if not dis[s] or e - s < 20:
            continue
        if sess[s] != sess[e - 1]:
            continue
        seg_I = np.median(np.abs(I[s:e]))
        seg_q = q[e - 1] - q[s]
        if seg_I < min_abs_current or seg_q < min_segment_ah:
            continue  # low-rate RPT or partial discharge
        ex = exp_um[s:e]
        ex = ex[np.isfinite(ex)]
        if ex.size < 20:
            continue
        rows.append({"efc": float(efc[s]), "discharge_Ah": float(seg_q),
                     "median_abs_current_A": float(seg_I),
                     "reversible_expansion_um": float(ex.max() - ex.min())})
    out = pd.DataFrame(rows).sort_values("efc").reset_index(drop=True)
    out.insert(0, "cyc_seq", np.arange(len(out)))
    return out


def validate_cell064(pkg_root):
    pkl = os.path.join(pkg_root, "cells", "CELL064", "lifetime",
                       "CELL064_lifetime_current_voltage_expansion_cleaned.pkl.gz")
    mine = per_cycle_amplitude(pkl)
    ref = pd.read_csv(os.path.join(HERE, "CELL064_reversible_expansion.csv"))
    m = np.interp(ref["efc"], mine["efc"], mine["reversible_expansion_um"])
    rel = (m - ref["reversible_expansion_um"]) / ref["reversible_expansion_um"]
    print(f"CELL064 validation: {len(mine)} cycles extracted vs {len(ref)} reference; "
          f"median rel. diff {100*np.median(rel):+.1f}%, IQR "
          f"[{100*np.percentile(rel,25):+.1f}, {100*np.percentile(rel,75):+.1f}]%; "
          f"corr {np.corrcoef(m, ref['reversible_expansion_um'])[0,1]:.3f}")
    return mine, ref


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pkg", required=True, help="IC_PYBAMM_FIG1 package root")
    ap.add_argument("--cell", default="CELL017")
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--out")
    a = ap.parse_args()
    if a.validate:
        validate_cell064(a.pkg)
    else:
        pkl = os.path.join(a.pkg, "cells", a.cell, "lifetime",
                           f"{a.cell}_lifetime_current_voltage_expansion_cleaned.pkl.gz")
        out = per_cycle_amplitude(pkl)
        path = a.out or os.path.join(HERE, f"{a.cell}_reversible_expansion_raw.csv")
        out.to_csv(path, index=False)
        print(f"{a.cell}: {len(out)} cycles -> {path}")
