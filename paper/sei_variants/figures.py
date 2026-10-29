"""Figures and scores of the single-law SEI variants against the paper runs (results/: kinetic term on graphite,
diffusion limit on silicon).

    python paper/sei_variants/figures.py

writes, for each variant with runs, figures/<variant>/compare_<cell>.png|pdf (SoH, graphite SEI, Si SEI, LLI),
figures/<variant>/<value>/<cell>_overall_summary.png|pdf, and scores.csv for all of them (mean |SoH gap| at the
C/20 RPTs, mean C/20 voltage RMSE, SEI and SoH at end of test).
"""
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
sys.path.insert(0, PAPER)
import sigr  # noqa: E402
from sigr import analysis as an, data, plots, style  # noqa: E402
from sweep import KEYS, VARIANTS  # noqa: E402

DATA = os.path.join(PAPER, "outputs", "dma_data", "dma_ours")  # our pyprobe-dma DMA, as in the paper
RUNS = os.path.join(HERE, "runs")
FIGS = os.path.join(HERE, "figures")
PAPER_RUNS = os.path.join(PAPER, "results")
TITLES = {"diffusion_both": r"Diffusion limit on both phases: $D_\mathrm{sol,Gr}$ in m$^2$ s$^{-1}$",
          "kinetic_both": r"Kinetic (EC) law on both phases: $k_\mathrm{SEI,Si}$ in m s$^{-1}$"}
SYMBOL = {"diffusion_both": "D", "kinetic_both": "k"}


def si_sei(res):
    """Si SEI (particle surface and crack faces) [A h]."""
    return res["Q_si_sei"] + res["Q_si_sei_cracks"]


def gr_sei(res):
    """Graphite SEI (particle surface and crack faces) [A h]: side-reaction capacity less the Si SEI."""
    return res["Q_side"] - si_sei(res)


def value_label(variant, d):
    """'D = 1.25e-23' etc.: the swept value nearest the folder name (which keeps only two digits)."""
    x = float(os.path.basename(d)[1:])
    v = min(VARIANTS[variant]["values"], key=lambda y: abs(y - x) / y)
    return f"{SYMBOL[variant]} = {v:.3g}"


def sweep_dirs(variant, key):
    dirs = sorted(glob.glob(os.path.join(RUNS, variant, "*")), key=lambda p: float(os.path.basename(p)[1:]))
    return [d for d in dirs if os.path.exists(os.path.join(d, f"{key}_results.npz"))]


def compare(variant, key):
    cond = sigr.CONDITIONS[key]
    style.use()
    fig, ax = plt.subplots(1, 4, figsize=(style.DOUBLE_IN, 2.7), layout="constrained")
    dirs = sweep_dirs(variant, key)
    runs = [("Paper run", PAPER_RUNS, "black", "--")]
    runs += [(value_label(variant, d), d, plt.cm.viridis(i / max(len(dirs) - 1, 1)), "-")
             for i, d in enumerate(dirs)]
    for label, d, c, ls in runs:
        res, _ = an.load(key, d)
        ax[0].plot(an.efc(res, res["rpt_thr"]), an.rpt_soh(res), ls, lw=1.0, marker="^", ms=2.5, color=c, label=label)
        m = an.life_mask(res)
        e = an.efc(res, res["thr"])[m]
        ax[1].plot(e, 100 * gr_sei(res)[m] / res["efc_basis_ah"], ls, lw=1.0, color=c)
        ax[2].plot(e, 100 * si_sei(res)[m] / res["efc_basis_ah"], ls, lw=1.0, color=c)
        ax[3].plot(e, an.lli_pct(res, True)[m], ls, lw=1.0, color=c)
    cap = data.load(cond, DATA, "capacity_fade")
    c20 = cap[cap["source"] == "C/20 RPT"]
    ax[0].plot(c20["efc"], c20["capacity_retention_pct"], "o", ms=4, color="red", label="Measured, C/20 RPT", zorder=6)
    hr = cap[cap["source"] != "C/20 RPT"]
    if len(hr):
        ax[0].plot(hr["efc"], hr["capacity_retention_pct"], "o", ms=4, mfc="white", mec="red", mew=0.8,
                   label="Measured, high-rate RPT", zorder=6)
    lli = data.load(cond, DATA, "LLI")
    ax[3].plot(lli["efc"], 100 * lli["charge_LLI_loss"] / lli["charge_LLI"].iloc[0], "D", ms=3.2, color="red",
               label="Measured (DMA)")
    for a, lab in zip(ax, ("SoH, C/20 RPT [%]", "Graphite SEI [% of 2.59 A h]", "Si SEI [% of 2.59 A h]", "LLI [%]")):
        a.set_ylabel(lab)
        a.set_xlabel(plots.EFC_LABEL)
    ax[3].legend(loc="upper left", fontsize=7)
    for i, a in enumerate(ax):
        style.panel_label(a, i)
    fig.legend(*ax[0].get_legend_handles_labels(), loc="outside lower center", ncol=6, fontsize=7,
               title=TITLES[variant], title_fontsize=7)
    return style.save(fig, os.path.join(FIGS, variant), f"compare_{key}")[0]


CELL_TITLES = {"25deg_full_soc": "25 °C, full window", "45deg_full_soc": "45 °C, full window",
               "25deg_15-95_soc": "25 °C, 15–95 %", "25deg_20-80_soc": "25 °C, 20–80 %"}


def si_figure(variant):
    """Supplementary figure: one row per cell (SoH, graphite SEI, Si SEI, LLI), every swept value against the
    paper run and the data."""
    style.use()
    fig, ax = plt.subplots(len(KEYS), 4, figsize=(style.DOUBLE_IN, 8.4), layout="constrained")
    values = sorted({os.path.basename(d) for k in KEYS for d in sweep_dirs(variant, k)}, key=lambda n: float(n[1:]))
    colour = {n: plt.cm.viridis(i / max(len(values) - 1, 1)) for i, n in enumerate(values)}
    for r, key in enumerate(KEYS):
        cond = sigr.CONDITIONS[key]
        runs = [("Paper run", PAPER_RUNS, "black", "--")]
        runs += [(value_label(variant, d), d, colour[os.path.basename(d)], "-")
                 for d in sweep_dirs(variant, key)]
        a = ax[r]
        for label, d, c, ls in runs:
            res, _ = an.load(key, d)
            a[0].plot(an.efc(res, res["rpt_thr"]), an.rpt_soh(res), ls, lw=0.8, marker="^", ms=2.5, color=c, label=label)
            m = an.life_mask(res)
            e = an.efc(res, res["thr"])[m]
            a[1].plot(e, 100 * gr_sei(res)[m] / res["efc_basis_ah"], ls, lw=0.8, color=c)
            a[2].plot(e, 100 * si_sei(res)[m] / res["efc_basis_ah"], ls, lw=0.8, color=c)
            a[3].plot(e, an.lli_pct(res, True)[m], ls, lw=0.8, color=c)
        cap = data.load(cond, DATA, "capacity_fade")
        c20 = cap[cap["source"] == "C/20 RPT"]
        a[0].plot(c20["efc"], c20["capacity_retention_pct"], "o", ms=3, color="red", label="Measured, C/20 RPT",
                  zorder=6)
        hr = cap[cap["source"] != "C/20 RPT"]
        if len(hr):
            a[0].plot(hr["efc"], hr["capacity_retention_pct"], "o", ms=3, mfc="white", mec="red", mew=0.8,
                      label="Measured, high-rate RPT", zorder=6)
        lli = data.load(cond, DATA, "LLI")
        a[3].plot(lli["efc"], 100 * lli["charge_LLI_loss"] / lli["charge_LLI"].iloc[0], "D", ms=2.5, color="red",
                  label="Measured (DMA)", zorder=6)
        for j, lab in enumerate(("SoH, C/20 RPT [%]", "Gr SEI [%]", "Si SEI [%]", "LLI [%]")):
            a[j].set_ylabel(lab)
            a[j].set_xlabel(plots.EFC_LABEL if r == len(KEYS) - 1 else "")
            style.panel_label(a[j], 4 * r + j, CELL_TITLES[key] if j == 3 else None)
    handles, labels = ax[0][0].get_legend_handles_labels()
    h3, l3 = ax[0][3].get_legend_handles_labels()
    fig.legend(handles + h3[-1:], labels + l3[-1:], loc="outside lower center", ncol=6, fontsize=7,
               title=TITLES[variant] + "; SEI in % of 2.59 A h (surface and crack faces)", title_fontsize=7)
    return style.save(fig, FIGS, f"si_{variant}")[0]


LOW_D = 1.25e-23  # diffusion_both run with a negligible graphite SEI (Supplementary Figs. S7-S11)
VARIANT = "Variant (diffusion limit on Gr)"  # legend name, so these figures are not mistaken for the paper model


def _relabel(ax, **kw):
    """Redraw the legend of ax with "Model" renamed to VARIANT."""
    leg = ax.get_legend()
    ax.legend(handles=leg.legend_handles, labels=[t.get_text().replace("Model", VARIANT) for t in leg.get_texts()],
              **kw)


def si_low_d():
    """Supplementary figures of the negligible-graphite-SEI run: a summary grid (one row per cell: SoH, LAM/LLI
    with our DMA, C/20 voltage) and, per cell, the C/20 voltage and dV/dSOC at each RPT."""
    d = os.path.join(RUNS, "diffusion_both", f"D{LOW_D:.1e}")
    out = os.path.join(FIGS, "low_d")
    style.use()
    fig, ax = plt.subplots(len(KEYS), 3, figsize=(style.DOUBLE_IN, 8.2), layout="constrained")
    paths = []
    for r, key in enumerate(KEYS):
        cond = sigr.CONDITIONS[key]
        res, curves = an.load(key, d)
        plots.soh(ax[r][0], res, cond, DATA)
        plots.lam_lli(ax[r][1], res, cond, DATA)
        plots.rpt_voltage(ax[r][2], res, curves, cond, DATA)
        _relabel(ax[r][2], loc="lower left")  # lower left: clear of the colour bar at this panel height
        if r:
            for a in ax[r][:2]:
                a.get_legend().remove()
        else:
            _relabel(ax[r][0], loc="upper right")
            _relabel(ax[r][1], loc="upper left")
        for j in range(3):
            style.panel_label(ax[r][j], 3 * r + j, CELL_TITLES[key] if j == 2 else None)
        for kind in ("voltage", "dvdsoc"):
            f = plots.rpt_dod(res, curves, cond, DATA, kind, model_label=VARIANT)
            paths.append(style.save(f, out, f"{key}_rpt_{kind}")[0])
            plt.close(f)
    paths.insert(0, style.save(fig, out, "summary")[0])
    plt.close(fig)
    return paths


def score(variant, label, d, key):
    cond = sigr.CONDITIONS[key]
    res, curves = an.load(key, d)
    s = an.summary(res, curves, cond, DATA)
    gaps = an.rpt_gaps(res, cond, DATA)
    return dict(variant=variant, run=label, condition=key, mean_abs_soh_gap_pp=round(s["mean_abs_soh_gap_pp"], 2),
                mean_voltage_rmse_mV=round(1e3 * s["mean_voltage_rmse_V"], 1),
                gr_sei_end_Ah=round(float(gr_sei(res)[-1]), 4), si_sei_end_Ah=round(float(si_sei(res)[-1]), 4),
                end_rpt_soh_pct=round(float(an.rpt_soh(res)[-1]), 1),
                soh_gaps_pp=" ".join(f"{g:+.1f}" for g in gaps["gap_pp"]))


if __name__ == "__main__":
    rows = [score("paper", "paper run (kinetic on Gr, diffusion limit on Si)", PAPER_RUNS, k) for k in KEYS]
    for variant in VARIANTS:
        for key in KEYS:
            dirs = sweep_dirs(variant, key)
            if not dirs:
                continue
            cond = sigr.CONDITIONS[key]
            for d in dirs:
                name = os.path.basename(d)
                res, curves = an.load(key, d)
                plots.save_all(res, curves, cond, DATA, os.path.join(FIGS, variant, name), only=("overall_summary",))
                rows.append(score(variant, value_label(variant, d), d, key))
            print(compare(variant, key))
        if any(sweep_dirs(variant, k) for k in KEYS):
            print(si_figure(variant))
    if os.path.isdir(os.path.join(RUNS, "diffusion_both", f"D{LOW_D:.1e}")):
        print("\n".join(si_low_d()))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(HERE, "scores.csv"), index=False)
    print(df.to_string(index=False))
