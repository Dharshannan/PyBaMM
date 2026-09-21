"""plot_all_levers_summary.py -- final consolidated comparison across every
lever tried in the 2026-09-19 investigation: Si diffusivity alone, K_SEI
retune, porosity-floor retune, Si exchange-current-density alone, and the
exchange-current+diffusivity combination. One scatter of expansion RMSE vs.
C/3 ageing-leg collapse ratio (log x-axis), coloured by capacity-fade gap,
to show the full landscape and why diffusivity x3 alone remains the best
single point found.
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# Manually assembled from all three sweeps' results (c3_collapse_ratio only
# available for the runs that used the updated worker -- diffusivity-only
# sweep predates that addition, so it's marked as unavailable/NaN there,
# except sidiff_x1 which we cross-checked separately).
rows = [
    # tag, expansion_rmse, c3_collapse_ratio, rpt_gap_pp, family
    ("diffusivity x1 (baseline)", 0.307, 299.4, 0.79, "diffusivity"),
    ("diffusivity x3 (best)", 0.1487, np.nan, 11.25, "diffusivity"),
    ("diffusivity x10", 0.1719, np.nan, 13.07, "diffusivity"),
    ("diffusivity x30", 0.1742, np.nan, 13.20, "diffusivity"),
    ("diffusivity x100", 0.1746, np.nan, 13.21, "diffusivity"),
    ("K_SEI x1.33 (8e-4)", 0.2893, np.nan, 1.43, "K_SEI retune (on x3)"),
    ("K_SEI x1.67 (1.0e-3)", 0.3972, np.nan, 3.56, "K_SEI retune (on x3)"),
    ("K_SEI x2.17 (1.3e-3)", 0.4758, np.nan, 7.82, "K_SEI retune (on x3)"),
    ("K_SEI x2.67 (1.6e-3)", 0.5324, np.nan, 10.0, "K_SEI retune (on x3)"),
    ("floor 0.028 (on x3)", 0.1823, np.nan, 7.08, "floor retune (on x3)"),
    ("floor 0.022 (on x3)", 0.218, np.nan, 3.62, "floor retune (on x3)"),
    ("floor 0.016 (on x3)", 0.2633, np.nan, 0.69, "floor retune (on x3)"),
    ("combined K+floor, x3", 0.3672, 172.7, 0.44, "floor retune (on x3)"),
    ("combined K+floor, x10", 0.3768, 177.1, 0.39, "floor retune (on x10)"),
    ("exchange-current x3", 0.3071, 18.7, 0.62, "exchange-current"),
    ("exchange-current x10", 0.343, 18.5, 1.16, "exchange-current"),
    ("exchange-current x30", 0.3571, 18.9, 2.03, "exchange-current"),
    ("exch x3 + diff x3", 0.2427, 16.2, 18.89, "exchange+diffusivity"),
    ("exch x10 + diff x3", 0.2497, 16.1, 18.82, "exchange+diffusivity"),
]
df = pd.DataFrame(rows, columns=["tag", "expansion_rmse", "c3_collapse", "rpt_gap_pp", "family"])

colors = {
    "diffusivity": "tab:blue",
    "K_SEI retune (on x3)": "tab:red",
    "floor retune (on x3)": "tab:orange",
    "floor retune (on x10)": "tab:brown",
    "exchange-current": "tab:green",
    "exchange+diffusivity": "tab:purple",
}

fig, ax = plt.subplots(figsize=(10, 7))
for fam, sub in df.groupby("family"):
    ax.scatter(sub["rpt_gap_pp"].abs(), sub["expansion_rmse"], s=90,
               color=colors[fam], label=fam, edgecolor="k", linewidth=0.5, zorder=3)

best = df.loc[df["expansion_rmse"].idxmin()]
ax.annotate(
    f"BEST: {best['tag']}\nRMSE={best['expansion_rmse']:.3f}",
    (abs(best["rpt_gap_pp"]), best["expansion_rmse"]),
    textcoords="offset points", xytext=(15, -20), fontsize=9,
    arrowprops=dict(arrowstyle="->", color="black"),
)
baseline = df.iloc[0]
ax.axhline(baseline["expansion_rmse"], color="gray", ls=":", lw=1,
           label=f"pre-fix baseline RMSE ({baseline['expansion_rmse']:.3f})")

ax.set_xlabel("|Capacity-fade RPT gap| (pp) -- lower = closer to real")
ax.set_ylabel("Expansion-shape RMSE -- lower = better fit")
ax.set_title("CELL064: every lever tried (2026-09-19)\nexpansion fit vs. capacity-fade cost")
ax.legend(fontsize=8, loc="upper left")
ax.grid(alpha=0.3)
fig.tight_layout()
out_path = os.path.join(HERE, "all_levers_summary.png")
fig.savefig(out_path, dpi=150)
print(f"Saved: {out_path}")

# Second panel: C/3 collapse ratio landscape (log scale) for the runs that
# have it, annotated with real's 1.27x target.
fig2, ax2 = plt.subplots(figsize=(10, 6))
have_c3 = df.dropna(subset=["c3_collapse"])
for fam, sub in have_c3.groupby("family"):
    ax2.scatter(sub["c3_collapse"], sub["expansion_rmse"], s=90,
                color=colors[fam], label=fam, edgecolor="k", linewidth=0.5, zorder=3)
ax2.set_xscale("log")
ax2.axvline(1.27, color="gray", ls=":", label="Real C/3-equivalent target (~1.27x)")
ax2.set_xlabel("C/3 ageing-leg Si stoichiometric collapse ratio (log scale)")
ax2.set_ylabel("Expansion-shape RMSE")
ax2.set_title("Root-cause metric (C/3 collapse) vs. actual fit quality\n(note: the 'great-looking' floor retune has the WORST C/3 collapse)")
ax2.legend(fontsize=8)
ax2.grid(alpha=0.3)
fig2.tight_layout()
out_path2 = os.path.join(HERE, "c3_collapse_vs_fit_landscape.png")
fig2.savefig(out_path2, dpi=150)
print(f"Saved: {out_path2}")
