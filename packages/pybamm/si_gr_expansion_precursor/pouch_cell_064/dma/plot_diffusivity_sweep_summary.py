"""plot_diffusivity_sweep_summary.py -- summarises si_diffusivity_sweep_results.csv
(mirrored from sweeps/baseline/) into one figure: expansion-shape RMSE,
Si sto-window collapse ratio, and capacity-fade RPT gap, all vs. the Si
diffusivity multiplier, plus the item-30-on-top-of-best-diffusivity results
as a separate annotated panel.
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "si_diffusivity_sweep_results.csv")


def main():
    df = pd.read_csv(CSV)
    pure_diff = df[df["tag"].str.match(r"^sidiff_x\d+$")].copy()
    pure_diff["mult"] = pure_diff["tag"].str.extract(r"sidiff_x(\d+)").astype(int)
    pure_diff = pure_diff.sort_values("mult")

    residual_rows = df[df["tag"].isin(
        ["sidiff_x3", "sidiff_x3_residual05", "sidiff_x3_residual09", "residual09_only", "sidiff_x1"]
    )].copy()

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    ax = axes[0]
    ax.plot(pure_diff["mult"], pure_diff["expansion_rmse"], "o-", color="tab:blue")
    ax.set_xscale("log")
    ax.set_xlabel("Si diffusivity multiplier")
    ax.set_ylabel("Expansion-shape RMSE (lower=better)")
    ax.set_title("Expansion fit vs. Si diffusivity")
    ax.grid(alpha=0.3)
    for _, r in pure_diff.iterrows():
        ax.annotate(f"{r['expansion_rmse']:.3f}", (r["mult"], r["expansion_rmse"]),
                    textcoords="offset points", xytext=(0, 8), fontsize=8, ha="center")

    ax = axes[1]
    ax2 = ax.twinx()
    l1, = ax.plot(pure_diff["mult"], pure_diff["collapse_ratio"], "s-", color="tab:red", label="C/20 RPT si_delta_sto collapse ratio")
    l2, = ax2.plot(pure_diff["mult"], pure_diff["rpt_gap_pp"], "^--", color="tab:purple", label="Capacity-fade RPT gap (pp)")
    ax.axhline(1.27, color="gray", ls=":", label="Real collapse ratio (1.27x)")
    ax.set_xscale("log")
    ax.set_xlabel("Si diffusivity multiplier")
    ax.set_ylabel("Si sto-window collapse ratio", color="tab:red")
    ax2.set_ylabel("Capacity-fade RPT gap [pp]", color="tab:purple")
    ax.set_title("Root-cause metric & capacity-fade side-effect")
    ax.legend(handles=[l1, l2, plt.Line2D([], [], color="gray", ls=":", label="Real collapse ratio (1.27x)")],
              fontsize=7, loc="upper right")
    ax.grid(alpha=0.3)

    ax = axes[2]
    order = ["sidiff_x1", "sidiff_x3", "sidiff_x3_residual05", "sidiff_x3_residual09", "residual09_only"]
    labels = ["baseline\n(x1, no residual)", "diffusivity x3\n(no residual)",
              "x3 + residual\n(frac=0.5)", "x3 + residual\n(frac=0.9)",
              "residual only\n(x1, frac=0.9)"]
    vals = [df[df["tag"] == t]["expansion_rmse"].iloc[0] for t in order]
    colors = ["gray", "tab:green", "tab:orange", "tab:red", "tab:blue"]
    ax.bar(labels, vals, color=colors)
    ax.set_ylabel("Expansion-shape RMSE (lower=better)")
    ax.set_title("Item 30 (isolated-Si residual) on top of\nthe best diffusivity fix -- does it help?")
    ax.tick_params(axis="x", labelrotation=20, labelsize=8)
    for i, v in enumerate(vals):
        ax.annotate(f"{v:.3f}", (i, v), textcoords="offset points", xytext=(0, 5), ha="center", fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    fig.suptitle("CELL064: Si diffusivity fix + item-30 residual sweep (2026-09-19)", fontsize=12)
    fig.tight_layout()
    out_path = os.path.join(HERE, "si_diffusivity_sweep_summary.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
