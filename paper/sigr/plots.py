"""Figures for one condition (same panels as the study's diagnostic plots).
Each function returns the matplotlib Figure; save_all writes them as
<key>_<figure>.png."""
import os

import matplotlib.pyplot as plt
import numpy as np

from . import analysis as an
from . import data

C_MODEL, C_RPT, C_GR, C_SI, C_PE, C_LLI, C_EXP = (
    "tab:blue", "tab:orange", "saddlebrown", "darkorange", "seagreen", "crimson", "tab:purple")


def _finish(ax, xlabel="EFC since RPT1 [-]"):
    ax.set_xlabel(xlabel)
    ax.grid(alpha=0.3)


def soh(ax, res, cond, data_dir):
    ax.plot(an.efc(res, res["age_thr"]), an.ageing_soh(res), "o-", ms=3, color=C_MODEL, label="Model (C/3 ageing)")
    ax.scatter(an.efc(res, res["rpt_thr"]), an.rpt_soh(res), marker="^", s=40, color=C_RPT,
               label="Model RPT (C/20)", zorder=5)
    cap = data.load(cond, data_dir, "capacity_fade")
    c20, sub = cap[cap["source"] == "C/20 RPT"], cap[cap["source"] != "C/20 RPT"]
    ax.scatter(c20["efc"], c20["capacity_retention_pct"], s=70, color="black", label="Experimental (C/20 RPT)", zorder=6)
    if len(sub):
        ax.scatter(sub["efc"], sub["capacity_retention_pct"], s=70, facecolors="none", edgecolors="black",
                   linewidths=1.4, label="Experimental (high-rate substitute)", zorder=6)
    k = data.knee_efc(cond, data_dir)
    ax.axvline(k, color="gray", ls=":", label=f"{cond.cell} knee (EFC ~{k:.0f})")
    ax.set_ylabel("SoH [%]")
    ax.set_title("SoH: model vs. experimental")
    ax.legend(fontsize=7)
    _finish(ax)


def lam_lli(ax, res, cond, data_dir):
    m = an.life_mask(res)
    e = an.efc(res, res["thr"])[m]
    for name, col, lab in (("LAM_gr", C_GR, "Gr"), ("LAM_si", C_SI, "Si"), ("LAM_pos", C_PE, "positive")):
        ax.plot(e, res[name][m], color=col, label=f"Model LAM {lab}")
    lam = data.load(cond, data_dir, "LAM_derived")
    for col_name, col, lab in (("LAM_NE_graphite_pct", C_GR, "Gr"), ("LAM_NE_silicon_pct", C_SI, "Si"),
                               ("LAM_PE_pct", C_PE, "PE")):
        ax.scatter(lam["efc"], lam[col_name], color=col, marker="s", label=f"Exp. LAM {lab}")
    ax.plot(e, an.lli_pct(res, False)[m], color=C_LLI, ls="--", label="Model LLI (side reactions)")
    ax.plot(e, an.lli_pct(res, True)[m], color=C_LLI, label="Model LLI (side reactions + LAM-trapped Li)")
    lli = data.load(cond, data_dir, "LLI")
    ax.scatter(lli["efc"], 100 * lli["charge_LLI_loss"] / lli["charge_LLI"].iloc[0], color=C_LLI, marker="D",
               label="Exp. LLI")
    ax.set_ylabel("LAM / LLI [%]")
    ax.set_title("LAM by phase + LLI: model vs. experimental")
    ax.legend(fontsize=6)
    _finish(ax)


def rpt_voltage(ax, res, curves, cond, data_dir):
    efcs = an.efc(res, [c["thr"] for c in curves])
    norm = plt.Normalize(efcs.min(), efcs.max())
    for c, e in zip(curves, efcs):
        ax.plot(c["q"], c["v"], color=plt.cm.viridis(norm(e)), label=f"model RPT @ EFC~{e:.0f}")
    for r, e_r, j, e_m in an.match_rpt_curves(res, curves, cond, data_dir, max_gap=60.0):
        q, v = data.load_rpt_discharge(cond, data_dir, r)
        note = f" (Δ{abs(e_m - e_r):.0f} EFC)" if abs(e_m - e_r) > 20 else ""
        ax.plot(q, v, color=plt.cm.viridis(norm(e_r)), ls="--", lw=1.3, label=f"real RPT{r} @ EFC~{e_r:.0f}{note}")
    ax.set_ylabel("Terminal voltage [V]")
    ax.set_title("C/20 RPT discharge (solid = model, dashed = real)")
    ax.legend(fontsize=6.5)
    _finish(ax, "Discharge capacity [A.h] (C/20)")


def expansion(ax, res, cond, data_dir, n_ref=5):
    ok = ~np.isnan(res["age_expansion_um"])
    e, x = an.efc(res, res["age_thr"])[ok], res["age_expansion_um"][ok]
    real = data.load(cond, data_dir, "reversible_expansion")
    rx = real["reversible_expansion_um"].to_numpy()
    ax.plot(e, x / np.median(x[:n_ref]), lw=1, color=C_EXP, label="Model (C/3 ageing)")
    ax.scatter(real["efc"], rx / np.median(rx[:n_ref]), s=8, alpha=0.4, color="black", label="Experimental")
    ax.set_ylabel("Reversible expansion / early-life value [-]")
    ax.set_title("Reversible expansion (normalised)")
    ax.legend(fontsize=8)
    _finish(ax)


def k_ratio(ax, res, cond, data_dir):
    ax.scatter(an.efc(res, res["rpt_thr"]), res["rpt_k"], marker="^", s=110, facecolors="none",
               edgecolors=C_RPT, linewidths=1.8, label="Model RPT (C/20)", zorder=7)
    k = data.load(cond, data_dir, "k_expansion_scale")
    c20, sub = k[k["source"] == "C/20 RPT"], k[k["source"] != "C/20 RPT"]
    ax.scatter(c20["efc"], c20["k"], s=70, color="black", label="Experimental (C/20 RPT)", zorder=6)
    if len(sub):
        ax.scatter(sub["efc"], sub["k"], s=70, facecolors="none", edgecolors="black", linewidths=1.4,
                   label="Experimental (high-rate substitute)", zorder=6)
    ax.set_ylabel("k = δ_rev,cell / δ_rev,particle [-]")
    ax.set_title("Expansion scale k")
    ax.legend(fontsize=8)
    _finish(ax)


def overall_summary(res, curves, cond, data_dir):
    fig, ax = plt.subplots(2, 3, figsize=(18, 9))
    soh(ax[0, 0], res, cond, data_dir)
    lam_lli(ax[0, 1], res, cond, data_dir)
    rpt_voltage(ax[0, 2], res, curves, cond, data_dir)
    expansion(ax[1, 0], res, cond, data_dir)
    k_ratio(ax[1, 1], res, cond, data_dir)
    ax[1, 2].axis("off")
    fig.suptitle(f"{cond.label}: overall summary (degradation + expansion)")
    fig.tight_layout()
    return fig


def _smooth_dvdsoc(soc, v, n_bins=120):
    """Bin V on a uniform depth-of-discharge grid, then differentiate."""
    o = np.argsort(soc)
    soc, v = np.asarray(soc)[o], np.asarray(v)[o]
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    mid = 0.5 * (edges[:-1] + edges[1:])
    vb = np.array([v[(soc >= lo) & (soc < hi)].mean() if ((soc >= lo) & (soc < hi)).any() else np.nan
                   for lo, hi in zip(edges[:-1], edges[1:])])
    ok = ~np.isnan(vb)
    return (mid[ok], np.gradient(vb[ok], mid[ok])) if ok.sum() >= 3 else (mid[ok], np.full(ok.sum(), np.nan))


def rpt_dod(res, curves, cond, data_dir, kind="voltage", trim=0.01):
    """Each real C/20 RPT vs its nearest model RPT, normalised to depth of
    discharge: kind = "voltage" (V) or "dvdsoc" (dV/dSOC)."""
    pairs = an.match_rpt_curves(res, curves, cond, data_dir)
    ncol = 2 if len(pairs) <= 4 else 4
    nrow = int(np.ceil(len(pairs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.5 * ncol, 4.5 * nrow), squeeze=False)
    for ax, (r, e_r, j, e_m) in zip(axes.flat, pairs):
        q_r, v_r = data.load_rpt_discharge(cond, data_dir, r)
        s_r, s_m, v_m = q_r / q_r.max(), curves[j]["q"] / curves[j]["q"].max(), curves[j]["v"]
        lab_r, lab_m = f"Real RPT{r} (EFC~{e_r:.0f})", f"Model RPT (EFC~{e_m:.0f})"
        if kind == "voltage":
            ax.plot(s_r, v_r, color="black", lw=2, label=lab_r)
            ax.plot(s_m, v_m, "--", color=C_MODEL, lw=1.8, label=lab_m)
            ax.set_ylabel("Terminal voltage [V]")
        else:  # the model's first ~1% DoD is the step-start transient; dropped
            keep = s_m >= trim
            ax.plot(*_smooth_dvdsoc(s_r, v_r), color="black", lw=2, label=lab_r)
            ax.plot(*_smooth_dvdsoc(s_m[keep], v_m[keep]), "--", color=C_MODEL, lw=1.8, label=lab_m)
            ax.set_ylabel("dV/d(SOC) [V]")
            ax.set_ylim(-3, 0.5)
        ax.set_title(f"RPT{r}")
        ax.legend(fontsize=9)
        _finish(ax, "Depth of discharge [-]")
    for ax in list(axes.flat)[len(pairs):]:
        ax.set_visible(False)
    fig.suptitle(f"{cond.label}: SOC-normalised C/20 RPT {'discharge' if kind == 'voltage' else 'dV/d(SOC)'}")
    fig.tight_layout()
    return fig


def crack_diagnostics(res, cond):
    e = an.efc(res, res["thr"])
    to_pct = 100 / 2.5947
    fig, ax = plt.subplots(4, 2, figsize=(11, 14))
    ax[0, 0].plot(e, res["crack_length"] / res["crack_length"][0])
    ax[0, 0].set_ylabel("Si crack length / initial [-]")
    ax[0, 1].plot(e, res["roughness"], color="tab:red")
    ax[0, 1].set_ylabel("Si roughness ratio [-]")
    ax[1, 0].plot(an.efc(res, res["cycle_thr"]), res["cycle_peak_stress_MPa"], ".", ms=3)
    ax[1, 0].set_ylabel("peak |Si surface tangential stress| per cycle [MPa]")
    ax[1, 1].plot(e, res["porosity"])
    ax[1, 1].set_ylabel("X-averaged negative electrode porosity [-]")
    ax[2, 0].plot(e, res["Q_si_sei"] * to_pct, label="Si SEI (bulk)")
    ax[2, 0].plot(e, res["Q_si_sei_cracks"] * to_pct, label="Si SEI on cracks")
    ax[2, 0].set_ylabel("LLI [% of nominal]")
    ax[2, 0].legend()
    ax[2, 1].plot(e, res["strain_E"])
    ax[2, 1].set_ylabel("cumulative Si contraction strain E [-]")
    ax[3, 0].plot(e, 1e9 * res["dl_paris"], label="route A (Paris, diffusion stress)")
    ax[3, 0].plot(e, 1e9 * res["dl_strain"], label="route B (strain fatigue)")
    ax[3, 0].set_ylabel("cumulative crack growth [nm]")
    ax[3, 0].legend()
    ax[3, 1].set_visible(False)
    for a in ax.flat:
        _finish(a)
    fig.suptitle(f"{cond.label}: Si crack / porosity diagnostics")
    fig.tight_layout()
    return fig


def save_all(res, curves, cond, data_dir, out_dir, dpi=150):
    os.makedirs(out_dir, exist_ok=True)
    figs = {
        "overall_summary": overall_summary(res, curves, cond, data_dir),
        "rpt_voltage_dod": rpt_dod(res, curves, cond, data_dir, "voltage"),
        "rpt_dvdsoc_dod": rpt_dod(res, curves, cond, data_dir, "dvdsoc"),
        "crack_diagnostics": crack_diagnostics(res, cond),
    }
    paths = {}
    for name, fig in figs.items():
        paths[name] = os.path.join(out_dir, f"{cond.key}_{name}.png")
        fig.savefig(paths[name], dpi=dpi)
        plt.close(fig)
    return paths
