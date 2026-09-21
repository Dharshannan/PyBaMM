"""build_final_rpt45_comparison.py -- final 3-way RPT4/RPT5 V(Q) and dV/dQ
comparison: real vs. true old baseline (2026-09-14 config, stability-fixed)
vs. new baseline (radius x0.2/diffusivity x3/BOL 0.17/K_SEI retuned).
Reads the raw curve CSVs already produced by build_true_old_baseline.py and
build_new_baseline_curves.py -- no simulation needed.
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REAL_RPT_DISCHARGE_DIR = (
    "C:/Shannan_PhD_Stuff_Local/Si_Gr_Expansion_Precursor/"
    "Pouch_Data/cell064_data/data_timeseries/low_rate_c20"
)


def load_real_rpt(rpt_num):
    path = os.path.join(REAL_RPT_DISCHARGE_DIR, f"CELL064_RPT{rpt_num:03d}_lowrate_c20_discharge.csv")
    df = pd.read_csv(path, usecols=["Q_discharge_Ah", "voltage_terminal_V"])
    q = df["Q_discharge_Ah"].to_numpy()
    v = df["voltage_terminal_V"].to_numpy()
    return q - q[0], v


def smooth_dvdq(q, v, n_bins=120):
    order = np.argsort(q)
    q, v = np.asarray(q)[order], np.asarray(v)[order]
    q_edges = np.linspace(q.min(), q.max(), n_bins + 1)
    q_mid = 0.5 * (q_edges[:-1] + q_edges[1:])
    v_binned = np.full(n_bins, np.nan)
    for i in range(n_bins):
        mask = (q >= q_edges[i]) & (q < q_edges[i + 1])
        if mask.any():
            v_binned[i] = v[mask].mean()
    valid = ~np.isnan(v_binned)
    q_mid, v_binned = q_mid[valid], v_binned[valid]
    dvdq = np.gradient(v_binned, q_mid)
    return q_mid, dvdq


def main():
    old = pd.read_csv(os.path.join(HERE, "true_old_baseline_rpt45_curves.csv"))
    new = pd.read_csv(os.path.join(HERE, "new_baseline_rpt45_curves.csv"))

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, rpt_num in zip(axes, [4, 5]):
        q_real, v_real = load_real_rpt(rpt_num)
        ax.plot(q_real, v_real, "-", color="black", lw=2, label=f"Real RPT{rpt_num}")
        o = old[old["rpt"] == rpt_num]
        ax.plot(o["q"], o["v"], "--", color="tab:red", lw=1.6, label="Old baseline (true, 2026-09-14)")
        n = new[new["rpt"] == rpt_num]
        ax.plot(n["q"], n["v"], "--", color="tab:blue", lw=1.6, label="New baseline")
        ax.set_xlabel("Discharge capacity Q [A.h]"); ax.set_ylabel("Terminal voltage [V]")
        ax.set_title(f"RPT{rpt_num} (EFC~{171.4 if rpt_num == 4 else 197.4:.0f}) discharge curve")
        ax.legend(fontsize=9); ax.grid(alpha=0.3)
    fig.suptitle("CELL064: RPT4/RPT5 discharge V(Q) -- real vs. true old vs. new baseline")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "final_rpt45_voltage_3way.png"), dpi=150)
    print("Saved: final_rpt45_voltage_3way.png")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, rpt_num in zip(axes, [4, 5]):
        q_real, v_real = load_real_rpt(rpt_num)
        qb, dvdq = smooth_dvdq(q_real, v_real)
        ax.plot(qb, dvdq, "-", color="black", lw=2, label=f"Real RPT{rpt_num}")
        o = old[old["rpt"] == rpt_num]
        qb, dvdq = smooth_dvdq(o["q"].to_numpy(), o["v"].to_numpy())
        ax.plot(qb, dvdq, "--", color="tab:red", lw=1.6, label="Old baseline (true, 2026-09-14)")
        n = new[new["rpt"] == rpt_num]
        qb, dvdq = smooth_dvdq(n["q"].to_numpy(), n["v"].to_numpy())
        ax.plot(qb, dvdq, "--", color="tab:blue", lw=1.6, label="New baseline")
        ax.set_xlabel("Discharge capacity Q [A.h]"); ax.set_ylabel("dV/dQ [V/A.h]")
        ax.set_title(f"RPT{rpt_num} (EFC~{171.4 if rpt_num == 4 else 197.4:.0f}) dV/dQ")
        ax.set_ylim(-3, 0.5)
        ax.legend(fontsize=9); ax.grid(alpha=0.3)
    fig.suptitle("CELL064: RPT4/RPT5 dV/dQ -- real vs. true old vs. new baseline")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "final_rpt45_dvdq_3way.png"), dpi=150)
    print("Saved: final_rpt45_dvdq_3way.png")


if __name__ == "__main__":
    main()
