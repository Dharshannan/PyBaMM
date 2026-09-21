"""
cell064_rpt_discharge_overlay.py -- overlay the model's own C/20 RPT
discharge curves (voltage vs. discharge capacity) against the REAL
experimental C/20 RPT discharge curves, to check how well the model
predicts both capacity fade AND voltage-profile shape as the cell degrades
-- not just the single-number SoH-at-RPT check score_rpt_gap does.

Real data: Si_Gr_Expansion_Precursor/Pouch_Data/cell064_data/data_timeseries/
low_rate_c20/CELL064_RPT00{1,2,4,5}_lowrate_c20_discharge.csv (RPT0/RPT3 have
no real C/20 discharge -- QC failed there, see experimental_data/NOTES.md).
Columns used: Q_discharge_Ah, voltage_terminal_V.

Runs the model once (via cell064_degradation_fit's own run_degradation(),
using whatever C64_* env vars this process is launched with -- pass the
same ones as the final adopted run) to get its own RPT voltage curves, each
tagged with its own EFC. Each real RPT's EFC comes from
load_experimental_capacity_fade() (already EFC_START_OFFSET-shifted); real
and model RPTs are then paired by NEAREST EFC (not by RPT index -- the
model's own RPT cadence drifts from the real RPT schedule post-knee, per
TUNING STATUS item 19's EFC-compression finding), so the pairing is
explicitly approximate and each subplot/legend entry states both EFCs.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import cell064_degradation_fit as c64

REAL_DIR = ("C:/Shannan_PhD_Stuff_Local/Si_Gr_Expansion_Precursor/Pouch_Data/"
            "cell064_data/data_timeseries/low_rate_c20")
REAL_RPTS = [1, 2, 4, 5]  # RPT0/RPT3 have no real C/20 discharge (QC failed)


def load_real_rpt(rpt_num):
    path = os.path.join(REAL_DIR, f"CELL064_RPT{rpt_num:03d}_lowrate_c20_discharge.csv")
    df = pd.read_csv(path, usecols=["Q_discharge_Ah", "voltage_terminal_V"])
    q = df["Q_discharge_Ah"].to_numpy()
    v = df["voltage_terminal_V"].to_numpy()
    return q - q[0], v


def main():
    print("Running model (this reuses cell064_degradation_fit's own config/env)...", flush=True)
    sol = c64.run_degradation()
    results = c64.extract_results(sol)

    exp_cap = c64.load_experimental_capacity_fade()
    real_efc_by_rpt = dict(zip(exp_cap["rpt"], exp_cap["efc"]))

    model_curves = results["rpt_voltage_curves"]
    if not model_curves:
        raise SystemExit("No model RPT voltage curves available -- run didn't reach an RPT.")
    model_efcs = c64.efc_from_throughput(np.array([c["thr"] for c in model_curves]))

    fig, ax = plt.subplots(figsize=(9, 7))
    cmap = plt.cm.viridis
    real_efcs = [real_efc_by_rpt[r] for r in REAL_RPTS if r in real_efc_by_rpt]
    norm = plt.Normalize(min(real_efcs), max(real_efcs))

    print(f"{'real RPT':>8} {'real EFC':>9} {'matched model EFC':>18} {'|delta|':>8}", flush=True)
    for rpt_num in REAL_RPTS:
        if rpt_num not in real_efc_by_rpt:
            continue
        real_efc = real_efc_by_rpt[rpt_num]
        q_real, v_real = load_real_rpt(rpt_num)
        j = int(np.argmin(np.abs(model_efcs - real_efc)))
        model_efc = float(model_efcs[j])
        print(f"{rpt_num:>8} {real_efc:>9.1f} {model_efc:>18.1f} {abs(model_efc - real_efc):>8.1f}",
              flush=True)
        color = cmap(norm(real_efc))
        ax.plot(q_real, v_real, "-", color=color, lw=2.5, alpha=0.85,
                label=f"Real RPT{rpt_num} (EFC {real_efc:.0f})")
        ax.plot(model_curves[j]["q"], model_curves[j]["v"], "--", color=color, lw=1.5,
                label=f"Model RPT (EFC {model_efc:.0f}, matched to RPT{rpt_num})")

    ax.set_xlabel("Discharge capacity [A.h] (C/20)")
    ax.set_ylabel("Terminal voltage [V]")
    ax.set_title("CELL064: C/20 RPT discharge curves, model vs. experimental\n"
                  "(solid = real, dashed = model, matched by nearest EFC -- see legend for both EFCs)")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3)
    plt.tight_layout()

    out_tag = os.environ.get("C64_OUT_TAG", "")
    out_name = f"cell064_rpt_discharge_overlay{('_' + out_tag) if out_tag else ''}.png"
    outpath = os.path.join(c64.SCRIPT_DIR, out_name)
    fig.savefig(outpath, dpi=150)
    print(f"Saved: {outpath}", flush=True)


if __name__ == "__main__":
    main()
