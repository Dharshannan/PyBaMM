"""plot_soc_4rpt_dvdq_configs.py -- same 3-figure x 4-subplot layout as
plot_soc_4rpt_configs.py, but dV/d(SOC) instead of V: both real and model
curves normalised to their own depth-of-discharge fraction [0,1] first
(so a difference in total capacity doesn't distort the derivative's
x-axis, mirroring the voltage plot's own SOC-normalisation), THEN
differentiated -- comparable across curves of different total capacity.
Reuses the same rpt45_curves_<tag>.csv files, no new simulation needed.
"""
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
REAL_RPT_DISCHARGE_DIR = ("C:/Shannan_PhD_Stuff_Local/Si_Gr_Expansion_Precursor/"
                          "Pouch_Data/cell064_data/data_timeseries/low_rate_c20")
REAL_EFC = {1: 0.0, 2: 49.4, 4: 171.4, 5: 197.4}

CONFIGS = [
    ("true_old_baseline_4rpt", "Old baseline (pre LAM/LLI-split retune)",
     "old_baseline_soc_4rpt_dvdq.png"),
    ("new_baseline_no_ocp", "New baseline, Si OCP deformation OFF",
     "new_baseline_no_ocp_soc_4rpt_dvdq.png"),
    ("new_baseline_ocp", "New baseline, Si OCP deformation ON (throughput driver)",
     "new_baseline_ocp_soc_4rpt_dvdq.png"),
]


def load_real_rpt(rpt_num):
    path = os.path.join(REAL_RPT_DISCHARGE_DIR,
                        f"CELL064_RPT{rpt_num:03d}_lowrate_c20_discharge.csv")
    df = pd.read_csv(path, usecols=["Q_discharge_Ah", "voltage_terminal_V"])
    q = df["Q_discharge_Ah"].to_numpy()
    v = df["voltage_terminal_V"].to_numpy()
    return q - q[0], v


def smooth_dvdsoc(soc, v, n_bins=120):
    order = np.argsort(soc)
    soc, v = np.asarray(soc)[order], np.asarray(v)[order]
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    mid = 0.5 * (edges[:-1] + edges[1:])
    v_binned = np.full(n_bins, np.nan)
    for i in range(n_bins):
        mask = (soc >= edges[i]) & (soc < edges[i + 1])
        if mask.any():
            v_binned[i] = v[mask].mean()
    valid = ~np.isnan(v_binned)
    mid, v_binned = mid[valid], v_binned[valid]
    if len(mid) < 3:
        return mid, np.full(len(mid), np.nan)
    dvdsoc = np.gradient(v_binned, mid)
    return mid, dvdsoc


def main():
    for tag, title, out_name in CONFIGS:
        csv_path = os.path.join(HERE, f"rpt45_curves_{tag}.csv")
        if not os.path.exists(csv_path):
            print(f"SKIP {tag}: {csv_path} not found")
            continue
        model = pd.read_csv(csv_path)

        fig, axes = plt.subplots(2, 2, figsize=(11, 9))
        for ax, rpt_num in zip(axes.flat, (1, 2, 4, 5)):
            q_real, v_real = load_real_rpt(rpt_num)
            soc_real = (q_real - q_real.min()) / (q_real.max() - q_real.min())
            soc_r, dvdsoc_r = smooth_dvdsoc(soc_real, v_real)
            ax.plot(soc_r, dvdsoc_r, "-", color="black", lw=2,
                    label=f"Real RPT{rpt_num} (EFC~{REAL_EFC[rpt_num]:.0f})")

            d = model[model["rpt"] == rpt_num].sort_values("q")
            if len(d):
                soc_model = (d["q"] - d["q"].min()) / (d["q"].max() - d["q"].min())
                soc_m, dvdsoc_m = smooth_dvdsoc(soc_model.to_numpy(), d["v"].to_numpy())
                ax.plot(soc_m, dvdsoc_m, "--", color="tab:blue", lw=1.8,
                        label=f"Model RPT{rpt_num}")
            else:
                ax.text(0.5, 0.5, "no matched model RPT", transform=ax.transAxes,
                        ha="center", va="center", color="tab:red")

            ax.set_xlabel("Depth of discharge [-]")
            ax.set_ylabel("dV/d(SOC) [V]")
            ax.set_title(f"RPT{rpt_num}")
            ax.set_ylim(-3, 0.5)
            ax.legend(fontsize=9)
            ax.grid(alpha=0.3)

        fig.suptitle(f"CELL064: SOC-normalised C/20 RPT dV/d(SOC) -- {title}")
        fig.tight_layout()
        out_path = os.path.join(HERE, out_name)
        fig.savefig(out_path, dpi=150)
        print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
