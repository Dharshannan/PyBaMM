"""Figures for one condition, drawn at their printed JES size (see style.py):
no titles, lettered panels. Each function returns the matplotlib Figure;
save_all writes them as <key>_<figure>.png and .pdf."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from . import analysis as an
from . import data
from . import style

C_MODEL, C_RPT, C_GR, C_SI, C_PE, C_LLI = (
    "tab:blue", "tab:orange", "saddlebrown", "darkorange", "seagreen", "crimson")
EFC_LABEL = "EFC since RPT1 [-]"


def _finish(ax, xlabel=EFC_LABEL):
    ax.set_xlabel(xlabel)


def soh(ax, res, cond, data_dir):
    ax.plot(an.efc(res, res["age_thr"]), an.ageing_soh(res), "-", lw=1.0, color=C_MODEL, label="Model, C/3 ageing")
    ax.plot(an.efc(res, res["rpt_thr"]), an.rpt_soh(res), "^", ms=4, color=C_RPT, label="Model, C/20 RPT", zorder=5)
    cap = data.load(cond, data_dir, "capacity_fade")
    c20, sub = cap[cap["source"] == "C/20 RPT"], cap[cap["source"] != "C/20 RPT"]
    ax.plot(c20["efc"], c20["capacity_retention_pct"], "o", ms=4.5, color="black", label="Measured, C/20 RPT",
            zorder=6)
    if len(sub):
        ax.plot(sub["efc"], sub["capacity_retention_pct"], "o", ms=4.5, mfc="white", mec="black", mew=1.0,
                label="Measured, high-rate RPT", zorder=6)
    ax.axvline(data.knee_efc(cond, data_dir), color="gray", ls=":", lw=1.0, label="Measured knee")
    ax.set_ylabel("SoH [%]")
    ax.legend(loc="lower left")
    _finish(ax)


def lam_lli(ax, res, cond, data_dir):
    m = an.life_mask(res)
    e = an.efc(res, res["thr"])[m]
    for name, col in (("LAM_gr", C_GR), ("LAM_si", C_SI), ("LAM_pos", C_PE)):
        ax.plot(e, res[name][m], color=col)
    lam = data.load(cond, data_dir, "LAM_derived")
    for col_name, col in (("LAM_NE_graphite_pct", C_GR), ("LAM_NE_silicon_pct", C_SI), ("LAM_PE_pct", C_PE)):
        ax.plot(lam["efc"], lam[col_name], "s", ms=3.5, color=col)
    ax.plot(e, an.lli_pct(res, False)[m], color=C_LLI, ls="--")
    ax.plot(e, an.lli_pct(res, True)[m], color=C_LLI)
    lli = data.load(cond, data_dir, "LLI")
    ax.plot(lli["efc"], 100 * lli["charge_LLI_loss"] / lli["charge_LLI"].iloc[0], "D", ms=3.2, color=C_LLI)
    handles = [Line2D([], [], color=c, marker=mk, ms=3.5, label=lab) for c, mk, lab in
               ((C_SI, "s", "LAM Si"), (C_GR, "s", "LAM graphite"), (C_PE, "s", "LAM positive"), (C_LLI, "D", "LLI"))]
    handles += [Line2D([], [], color="0.3", label="Model"),
                Line2D([], [], color="0.3", ls="none", marker="s", ms=3.5, label="Measured (DMA)"),
                Line2D([], [], color=C_LLI, ls="--", label="LLI, side reactions only")]
    ax.legend(handles=handles, loc="upper left")
    ax.set_ylabel("LAM, LLI [%]")
    _finish(ax)


def rpt_voltage(ax, res, curves, cond, data_dir):
    efcs = an.efc(res, [c["thr"] for c in curves])
    norm = plt.Normalize(efcs.min(), efcs.max())
    cmap = plt.cm.viridis
    for c, e in zip(curves, efcs):
        ax.plot(c["q"], c["v"], color=cmap(norm(e)), lw=1.0)
    for r, e_r, j, e_m in an.match_rpt_curves(res, curves, cond, data_dir, max_gap=60.0):
        q, v = data.load_rpt_discharge(cond, data_dir, r)
        ax.plot(q, v, color=cmap(norm(e_r)), ls="--", lw=1.0)
    ax.legend(handles=[Line2D([], [], color="0.3", label="Model"), Line2D([], [], color="0.3", ls="--", label="Measured")],
              loc="lower left", bbox_to_anchor=(0.0, 0.27))
    cax = ax.inset_axes([0.06, 0.11, 0.42, 0.045])
    cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cb.set_label("EFC since RPT1", fontsize=7, labelpad=1)
    cb.ax.xaxis.set_label_position("top")
    top = 10 * np.floor(efcs.max() / 10)
    cb.set_ticks([0, 10 * np.round(top / 20), top])
    cb.ax.tick_params(labelsize=7, length=1.5, pad=1)
    cb.outline.set_linewidth(0.4)
    cax.grid(False)
    ax.set_ylabel("Terminal voltage [V]")
    _finish(ax, "Discharge capacity, C/20 [A h]")


def overall_summary(res, curves, cond, data_dir):
    fig, ax = plt.subplots(1, 3, figsize=(style.DOUBLE_IN, 2.55), layout="constrained")
    soh(ax[0], res, cond, data_dir)
    lam_lli(ax[1], res, cond, data_dir)
    rpt_voltage(ax[2], res, curves, cond, data_dir)
    for i, a in enumerate(ax):
        style.panel_label(a, i)
    return fig


_smooth_dvdsoc = an.smooth_dvdsoc  # one definition, shared with the dV/dSOC score


def rpt_dod(res, curves, cond, data_dir, kind="voltage", trim=0.01):
    """Each real C/20 RPT vs its nearest model RPT, normalised to depth of
    discharge: kind = "voltage" (V) or "dvdsoc" (dV/dSOC). Up to four RPTs fit
    one column (2 x 2); more take two columns (2 x 3 or 2 x 4)."""
    pairs = an.match_rpt_curves(res, curves, cond, data_dir)
    ncol = 2 if len(pairs) <= 4 else (3 if len(pairs) <= 6 else 4)
    nrow = int(np.ceil(len(pairs) / ncol))
    width = style.SINGLE_IN if ncol == 2 else style.DOUBLE_IN
    fig, axes = plt.subplots(nrow, ncol, figsize=(width, 1.45 * nrow + 0.35), sharex=True, sharey=True,
                             squeeze=False, layout="constrained")
    for i, (ax, (r, e_r, j, e_m)) in enumerate(zip(axes.flat, pairs)):
        q_r, v_r = data.load_rpt_discharge(cond, data_dir, r)
        s_r, s_m, v_m = q_r / q_r.max(), curves[j]["q"] / curves[j]["q"].max(), curves[j]["v"]
        if kind == "voltage":
            ax.plot(s_r, v_r, color="black", lw=1.3, label="Measured")
            ax.plot(s_m, v_m, "--", color=C_MODEL, lw=1.2, label="Model")
        else:  # the model's first ~1% DoD is the step-start transient; dropped
            keep = s_m >= trim
            ax.plot(*_smooth_dvdsoc(s_r, v_r), color="black", lw=1.3, label="Measured")
            ax.plot(*_smooth_dvdsoc(s_m[keep], v_m[keep]), "--", color=C_MODEL, lw=1.2, label="Model")
            ax.set_ylim(-3, 0.5)
        style.panel_label(ax, i, f"RPT{r}, EFC {e_r:.0f}")
        if i == 0:
            ax.legend(loc="lower left")
    for ax in list(axes.flat)[len(pairs):]:
        ax.set_visible(False)
    for ax in axes[:, 0]:
        ax.set_ylabel("Terminal voltage [V]" if kind == "voltage" else "dV/dSoC [V]")
    for c in range(ncol):  # x label on the lowest visible panel of each column
        col = [a for a in axes[:, c] if a.get_visible()]
        if col:
            col[-1].set_xlabel("Depth of discharge [-]")
            col[-1].xaxis.set_tick_params(labelbottom=True)
    return fig


def crack_diagnostics(res, cond):
    e = an.efc(res, res["thr"])
    to_pct = 100 / 2.5947
    fig, axs = plt.subplots(4, 2, figsize=(style.DOUBLE_IN, 6.2), sharex=True, layout="constrained")
    ax = axs.flat
    ax[0].plot(e, res["crack_length"] / res["crack_length"][0])
    ax[0].set_ylabel("Si crack length / initial [-]")
    ax[1].plot(e, res["roughness"], color="tab:red")
    ax[1].set_ylabel("Si roughness ratio [-]")
    ax[2].plot(an.efc(res, res["cycle_thr"]), res["cycle_peak_stress_MPa"], ".", ms=2)
    ax[2].set_ylabel("Peak Si surface\ntangential stress [MPa]")
    # Per-cycle minimum (charged end of each ageing cycle). The first ageing
    # cycle after an RPT starts from the RPT's full charge rather than the
    # window top, so for partial windows its minimum is lower: left out.
    after_rpt = np.searchsorted(res["age_thr"], res["rpt_thr"], side="right")
    keep = np.ones(res["age_thr"].size, bool)
    keep[after_rpt[after_rpt < keep.size]] = not cond.partial
    ax[3].plot(np.r_[0.0, an.efc(res, res["age_thr"])[keep]],
               np.r_[res["initial_porosity"], res["age_min_porosity"][keep]], ".", ms=2)
    ax[3].set_ylabel("Negative electrode porosity,\nper-cycle minimum [-]")
    ax[4].plot(e, res["Q_si_sei"] * to_pct, label="Si SEI, particle surface")
    ax[4].plot(e, res["Q_si_sei_cracks"] * to_pct, label="Si SEI, crack faces")
    ax[4].set_ylabel("LLI [% of nominal]")
    ax[4].legend(loc="upper left")
    ax[5].plot(e, res["strain_E"])
    ax[5].set_ylabel("Cumulative Si\ncontraction strain [-]")
    ax[6].plot(e, 1e9 * res["dl_paris"], label="Route A (Paris, diffusion stress)")
    ax[6].plot(e, 1e9 * res["dl_strain"], label="Route B (volume-strain fatigue)")
    ax[6].set_ylabel("Cumulative crack growth [nm]")
    ax[6].legend(loc="upper left")
    ax[7].set_visible(False)
    ax[5].xaxis.set_tick_params(labelbottom=True)
    for i, a in enumerate(ax[:7]):
        style.panel_label(a, i)
    for a in (ax[5], ax[6]):
        _finish(a)
    return fig


def save_all(res, curves, cond, data_dir, out_dir):
    """Write the four per-condition figures (PNG at 300 dpi and PDF) to out_dir."""
    style.use()
    figs = {
        "overall_summary": overall_summary(res, curves, cond, data_dir),
        "rpt_voltage_dod": rpt_dod(res, curves, cond, data_dir, "voltage"),
        "rpt_dvdsoc_dod": rpt_dod(res, curves, cond, data_dir, "dvdsoc"),
        "crack_diagnostics": crack_diagnostics(res, cond),
    }
    paths = {}
    for name, fig in figs.items():
        paths[name] = style.save(fig, out_dir, f"{cond.key}_{name}")[0]
        plt.close(fig)
    return paths
