"""build_true_old_baseline.py -- reconstructs the EXACT canonical, pre-
session CELL064 baseline (the run that produced cell064_degradation_fit_
result.png / cell064_expansion_fit_result.png / cell064_overall_summary_
result.png, dated 2026-09-14, and motivated the whole Si-LAM investigation)
by running cell064_degradation_fit.py with ZERO env var overrides.

Confirmed via README.md's documented "item 19 final: run07 recipe" (C64_K_
SEI_MULT=1.0e-3, C64_NEG_POROSITY_FLOOR=0.023, C64_EXPONENT_MAX_SEI=100,
C64_LAM_PROP_MULT=0.6, C64_SI_BETA_LAM_SEI=7e-6, BoL porosity 0.25,
F0/width=0.7/0.05) that these are EXACTLY the file's own current hardcoded
defaults (verified by reading each os.environ.get(...) call directly) --
none of this session's newer additions (redirect-to-LAM, porosity-isolation
beta, diffusivity/radius/exchange-current multipliers) existed at that
point, and all default OFF/neutral when unset, so leaving every C64_* env
var completely unset reconstructs that exact historical config.

Produces: the standard plot_all() suite (tagged "true_old_baseline", NOT
overwriting the canonical untagged files), plus raw RPT4/RPT5 V(Q) curves
and dV/dQ for direct comparison against the new baseline and real data.
"""
import os
import sys

for k in list(os.environ):
    if k.startswith("C64_"):
        del os.environ[k]
os.environ["C64_OUT_TAG"] = "true_old_baseline"
os.environ["C64_NO_PLOT"] = "0"
# Stability fix (2026-09-20 debugging): the exact item-19/run07 floor
# (0.023) now crashes (IDA_BAD_K) shortly post-knee -- something in the
# shared submodel code has shifted the numerically-safe margin since then
# (this session's own later work needed 0.035 for a similar reason).
# 0.028 is the smallest bump found that survives the full run without
# changing any of the other five documented run07 parameters. MAX_CYCLES
# raised to 350 to bring RPT5 back into range (documented in README as
# item 19 phase 4's own reason for the same raise).
os.environ["C64_NEG_POROSITY_FLOOR"] = "0.028"
os.environ["C64_MAX_CYCLES"] = "350"

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DEG_DIR = os.path.join(os.path.dirname(HERE), "degradation_test_matrix")
sys.path.insert(0, DEG_DIR)
import cell064_degradation_fit as c64  # noqa: E402

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
    q, v = q[order], v[order]
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


def nearest_curve(sim_efc_rpt, curves, target_efc):
    if not curves or sim_efc_rpt.size == 0:
        return None, None
    idx = int(np.argmin(np.abs(sim_efc_rpt - target_efc)))
    c = curves[idx]
    return c["q"], c["v"]


def main():
    print("Config check (should all be defaults):", flush=True)
    print(f"  K_SEI_MULT={c64.K_SEI_MULT}, NEG_POROSITY_FLOOR={c64.NEG_POROSITY_FLOOR}, "
          f"EXPONENT_MAX_SEI={c64.EXPONENT_MAX_SEI}, LAM_PROP_MULT={c64.LAM_PROP_MULT}, "
          f"SI_BETA_LAM_SEI={c64.SI_BETA_LAM_SEI}, NEG_POROSITY_BOL={c64.NEG_POROSITY_BOL}, "
          f"F0={c64.F0_BASELINE}, WIDTH={c64.WIDTH_BASELINE}", flush=True)

    sol = c64.run_degradation()
    results = c64.extract_results(sol)
    c64.plot_all(results, sol)

    sim_efc_rpt = c64.efc_from_throughput(results["rpt_thr"]) if results["rpt_thr"].size else np.array([])
    curves = results["rpt_voltage_curves"]
    real_efc = {4: 171.4, 5: 197.4}

    # Save raw curves to CSV for reuse.
    rows = []
    for rpt_num, efc in real_efc.items():
        q, v = nearest_curve(sim_efc_rpt, curves, efc)
        if q is None:
            continue
        for qi, vi in zip(q, v):
            rows.append({"rpt": rpt_num, "q": qi, "v": vi})
    pd.DataFrame(rows).to_csv(os.path.join(HERE, "true_old_baseline_rpt45_curves.csv"), index=False)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, rpt_num in zip(axes, [4, 5]):
        q_real, v_real = load_real_rpt(rpt_num)
        ax.plot(q_real, v_real, "-", color="black", lw=2, label=f"Real RPT{rpt_num}")
        q_old, v_old = nearest_curve(sim_efc_rpt, curves, real_efc[rpt_num])
        if q_old is not None:
            ax.plot(q_old, v_old, "--", color="tab:red", lw=1.6, label="Old baseline (true, 2026-09-14 config)")
        ax.set_xlabel("Discharge capacity Q [A.h]"); ax.set_ylabel("Terminal voltage [V]")
        ax.set_title(f"RPT{rpt_num} (EFC~{real_efc[rpt_num]:.0f}) discharge curve")
        ax.legend(fontsize=9); ax.grid(alpha=0.3)
    fig.suptitle("CELL064: RPT4/RPT5 V(Q) -- real vs. TRUE old baseline (2026-09-14 config)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "true_old_baseline_rpt45_voltage.png"), dpi=150)
    print("Saved: true_old_baseline_rpt45_voltage.png", flush=True)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, rpt_num in zip(axes, [4, 5]):
        q_real, v_real = load_real_rpt(rpt_num)
        qb, dvdq = smooth_dvdq(q_real, v_real)
        ax.plot(qb, dvdq, "-", color="black", lw=2, label=f"Real RPT{rpt_num}")
        q_old, v_old = nearest_curve(sim_efc_rpt, curves, real_efc[rpt_num])
        if q_old is not None:
            qb, dvdq = smooth_dvdq(np.asarray(q_old), np.asarray(v_old))
            ax.plot(qb, dvdq, "--", color="tab:red", lw=1.6, label="Old baseline (true, 2026-09-14 config)")
        ax.set_xlabel("Discharge capacity Q [A.h]"); ax.set_ylabel("dV/dQ [V/A.h]")
        ax.set_title(f"RPT{rpt_num} (EFC~{real_efc[rpt_num]:.0f}) dV/dQ")
        ax.set_ylim(-3, 0.5)
        ax.legend(fontsize=9); ax.grid(alpha=0.3)
    fig.suptitle("CELL064: RPT4/RPT5 dV/dQ -- real vs. TRUE old baseline (2026-09-14 config)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "true_old_baseline_rpt45_dvdq.png"), dpi=150)
    print("Saved: true_old_baseline_rpt45_dvdq.png", flush=True)


if __name__ == "__main__":
    main()
