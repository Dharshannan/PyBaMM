"""plot_soc_4rpt_configs.py -- three separate figures (old baseline / new
baseline without Si OCP deformation / new baseline with it), each with 4
subplots (RPT1, RPT2, RPT4, RPT5), real vs. model, SOC-normalised (depth
of discharge), all at C/20. Reads the rpt45_curves_<tag>.csv files
task_soc_4rpt_plots.py's worker run produces.
"""
import os

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
REAL_RPT_DISCHARGE_DIR = ("C:/Shannan_PhD_Stuff_Local/Si_Gr_Expansion_Precursor/"
                          "Pouch_Data/cell064_data/data_timeseries/low_rate_c20")
REAL_EFC = {1: 0.0, 2: 49.4, 4: 171.4, 5: 197.4}

CONFIGS = [
    ("true_old_baseline_4rpt", "Old baseline (pre LAM/LLI-split retune)",
     "old_baseline_soc_4rpt.png"),
    ("new_baseline_no_ocp", "New baseline, Si OCP deformation OFF",
     "new_baseline_no_ocp_soc_4rpt.png"),
    ("new_baseline_ocp", "New baseline, Si OCP deformation ON (throughput driver)",
     "new_baseline_ocp_soc_4rpt.png"),
]


def load_real_rpt(rpt_num):
    path = os.path.join(REAL_RPT_DISCHARGE_DIR,
                        f"CELL064_RPT{rpt_num:03d}_lowrate_c20_discharge.csv")
    df = pd.read_csv(path, usecols=["Q_discharge_Ah", "voltage_terminal_V"])
    q = df["Q_discharge_Ah"].to_numpy()
    v = df["voltage_terminal_V"].to_numpy()
    return q - q[0], v


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
            ax.plot(soc_real, v_real, "-", color="black", lw=2,
                    label=f"Real RPT{rpt_num} (EFC~{REAL_EFC[rpt_num]:.0f})")

            d = model[model["rpt"] == rpt_num].sort_values("q")
            if len(d):
                soc_model = (d["q"] - d["q"].min()) / (d["q"].max() - d["q"].min())
                ax.plot(soc_model, d["v"], "--", color="tab:blue", lw=1.8,
                        label=f"Model RPT{rpt_num}")
            else:
                ax.text(0.5, 0.5, "no matched model RPT", transform=ax.transAxes,
                        ha="center", va="center", color="tab:red")

            ax.set_xlabel("Depth of discharge [-]")
            ax.set_ylabel("Terminal voltage [V]")
            ax.set_title(f"RPT{rpt_num}")
            ax.legend(fontsize=9)
            ax.grid(alpha=0.3)

        fig.suptitle(f"CELL064: SOC-normalised C/20 RPT discharge -- {title}")
        fig.tight_layout()
        out_path = os.path.join(HERE, out_name)
        fig.savefig(out_path, dpi=150)
        print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
