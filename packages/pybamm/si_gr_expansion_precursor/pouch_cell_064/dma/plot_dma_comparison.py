"""plot_dma_comparison.py -- CELL064 Si stoichiometric-window collapse
comparison, now split into three separate figures per feedback (a single
combined figure got too crowded once the new baseline was added):

  1. dma_si_real_window.png -- real per-RPT Si/Gr stoichiometric window
     (unchanged reference plot, from real_si_gr_sto_window_by_rpt.csv).
  2. dma_si_c20rpt_comparison.png -- the apples-to-apples comparison: real
     Si (C/20) vs. old baseline's own C/20 RPT legs vs. new baseline's own
     C/20 RPT legs, all normalised to their first point.
  3. dma_si_c3_comparison.png -- old vs. new baseline's C/3 ageing-leg
     collapse (no real curve here -- the real DMA data has no C/3-rate
     equivalent).

Reads:
  - real_si_gr_sto_window_by_rpt.csv (real data, always present)
  - model_si_sto_window_c3_and_c20rpt.csv (old baseline, from
    extract_model_si_sto_window_c20rpt.py)
  - model_si_sto_window_new_baseline.csv (new baseline, from
    extract_new_baseline_sto_collapse.py)
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REAL_CSV = os.path.join(HERE, "real_si_gr_sto_window_by_rpt.csv")
OLD_CSV = os.path.join(HERE, "model_si_sto_window_c3_and_c20rpt.csv")
NEW_CSV = os.path.join(HERE, "model_si_sto_window_new_baseline.csv")


def main():
    real = pd.read_csv(REAL_CSV)
    real_norm = real["x_si_delta"] / real["x_si_delta"].iloc[real["rpt"].tolist().index(1)]

    old = pd.read_csv(OLD_CSV) if os.path.exists(OLD_CSV) else None
    new = pd.read_csv(NEW_CSV) if os.path.exists(NEW_CSV) else None

    # --- Figure 1: real-only reference (unchanged) ---
    fig, ax = plt.subplots(figsize=(7, 5.5))
    ax.plot(real["efc"], real["x_si_delta"], "o-", color="tab:blue", label="Real Si (DMA, x_si_delta, C/20)")
    ax.plot(real["efc"], real["x_gr_delta"], "s--", color="tab:green", label="Real Gr (DMA, x_gr_delta, C/20)")
    ax.set_xlabel("EFC")
    ax.set_ylabel("Stoichiometric window, delta_sto (x_max - x_min)")
    ax.set_title("CELL064 real per-RPT Si/Gr stoichiometric window\n(from discharge_esoh_fit model_x_si/model_x_gr traces)")
    ax.set_ylim(0, 1.05)
    ax.legend(loc="lower left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "dma_si_real_window.png"), dpi=150)
    print("Saved: dma_si_real_window.png")

    # --- Figure 2: C/20 RPT comparison (real, old, new) ---
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(real["efc"], real_norm, "o-", color="black", lw=2, label="Real Si, C/20 (norm. to RPT1)")
    if old is not None:
        rpt_old = old[old["kind"] == "rpt_c20"].sort_values("efc")
        if len(rpt_old):
            norm = rpt_old["si_delta_sto"] / rpt_old["si_delta_sto"].iloc[0]
            ax.plot(rpt_old["efc"], norm, "^-", color="tab:red", alpha=0.9,
                     label="Old baseline, C/20 RPT legs")
    if new is not None:
        rpt_new = new[new["kind"] == "rpt_c20"].sort_values("efc")
        if len(rpt_new):
            norm = rpt_new["si_delta_sto"] / rpt_new["si_delta_sto"].iloc[0]
            ax.plot(rpt_new["efc"], norm, "v-", color="tab:blue", alpha=0.9,
                     label="New baseline, C/20 RPT legs")
    ax.set_xlabel("EFC")
    ax.set_ylabel("Normalised stoichiometric window (-)")
    ax.set_title("Si stoichiometric-window collapse: C/20 RPT\n(real vs. old vs. new baseline -- the apples-to-apples comparison)")
    ax.legend(loc="lower left", fontsize=9)
    ax.grid(alpha=0.3)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "dma_si_c20rpt_comparison.png"), dpi=150)
    print("Saved: dma_si_c20rpt_comparison.png")

    # --- Figure 3: C/3 ageing comparison (old, new -- no real equivalent) ---
    fig, ax = plt.subplots(figsize=(8, 6))
    if old is not None:
        c3_old = old[old["kind"] == "ageing_c3"].sort_values("efc")
        if len(c3_old):
            norm = c3_old["si_delta_sto"] / c3_old["si_delta_sto"].iloc[0]
            ax.plot(c3_old["efc"], norm, "-", color="tab:red", alpha=0.8, lw=1.2,
                     label="Old baseline, C/3 ageing legs")
    if new is not None:
        c3_new = new[new["kind"] == "ageing_c3"].sort_values("efc")
        if len(c3_new):
            norm = c3_new["si_delta_sto"] / c3_new["si_delta_sto"].iloc[0]
            ax.plot(c3_new["efc"], norm, "-", color="tab:blue", alpha=0.8, lw=1.2,
                     label="New baseline, C/3 ageing legs")
    ax.set_xlabel("EFC")
    ax.set_ylabel("Normalised stoichiometric window (-)")
    ax.set_title("Si stoichiometric-window collapse: C/3 ageing legs\n(old vs. new baseline -- no real-data equivalent at this rate)")
    ax.legend(loc="lower left", fontsize=9)
    ax.grid(alpha=0.3)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "dma_si_c3_comparison.png"), dpi=150)
    print("Saved: dma_si_c3_comparison.png")


if __name__ == "__main__":
    main()
