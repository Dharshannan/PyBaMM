"""The paper's figures, all at their printed JES size (see style.py):

    import sigr.figures as F
    F.make_all("data", "results", "outputs_ablations", "paper_pics")

Per condition (plots.py): overall summary, C/20 RPT voltage and dV/dSoC
against depth of discharge, Si crack/porosity diagnostics. Shared: the
measured-data overview, the two mechanism schematics, the ablations
(ablations.py) and the partial-window comparison of the Discussion.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle

from . import ablations
from . import analysis as an
from . import data, plots, style
from .conditions import CONDITIONS

C_GR, C_SI, C_SEI, C_EL, C_ISO = "#8c8c8c", "#f39c34", "#7a4a1e", "#d6ecfa", "#c0392b"
CELL_COLOR = {"25deg_full_soc": "tab:blue", "45deg_full_soc": "tab:red", "25deg_15-95_soc": "tab:green",
              "25deg_20-80_soc": "tab:purple"}
EPS0, EPS_MIN = 0.17, 0.035  # initial and floor porosity of the negative electrode


# ---------------------------------------------------------------------------
# Measured data overview
# ---------------------------------------------------------------------------
def data_overview(data_dir):
    """Measured SoH, negative- and positive-electrode LAM and LLI of the four cells."""
    fig, axs = plt.subplots(2, 2, figsize=(style.DOUBLE_IN, 4.3), sharex=True, layout="constrained")
    ax = axs.ravel()
    for k, c in CONDITIONS.items():
        col = CELL_COLOR[k]
        cap = data.load(c, data_dir, "capacity_fade").sort_values("efc")
        c20, sub = cap[cap.source == "C/20 RPT"], cap[cap.source != "C/20 RPT"]
        ax[0].plot(cap.efc, cap.capacity_retention_pct, "-", color=col, lw=0.9)
        ax[0].plot(c20.efc, c20.capacity_retention_pct, "o", ms=3.5, color=col)
        if len(sub):
            ax[0].plot(sub.efc, sub.capacity_retention_pct, "o", ms=3.5, mfc="white", mec=col, mew=1.0)
        ax[0].axvline(data.knee_efc(c, data_dir), color=col, ls=":", lw=0.9)
        lam = data.load(c, data_dir, "LAM_derived")
        ax[1].plot(lam.efc, lam.LAM_NE_silicon_pct, "s-", ms=3, color=col)
        ax[1].plot(lam.efc, lam.LAM_NE_graphite_pct, "s--", ms=3, color=col, alpha=0.55)
        ax[2].plot(lam.efc, lam.LAM_PE_pct, "^-", ms=3, color=col)
        lli = data.load(c, data_dir, "LLI")
        ax[3].plot(lli.efc, 100 * lli.charge_LLI_loss / lli.charge_LLI.iloc[0], "D-", ms=3, color=col)
    ax[2].axhline(0, color="k", lw=0.5)
    for a, yl in zip(ax, ("SoH [%]", "LAM, negative electrode [%]", "LAM, positive electrode [%]", "LLI [%]")):
        a.set_ylabel(yl)
    for i, a in enumerate(ax):
        style.panel_label(a, i)
    for a in ax[2:]:
        a.set_xlabel(plots.EFC_LABEL)
    cells = [Line2D([], [], color=CELL_COLOR[k], marker="o", ms=3.5, label=ablations.cond_tag(c))
             for k, c in CONDITIONS.items()]
    marks = [Line2D([], [], color="0.3", ls="none", marker="o", ms=3.5, label="C/20 RPT"),
             Line2D([], [], color="0.3", ls="none", marker="o", ms=3.5, mfc="white", label="High-rate RPT (CELL064)"),
             Line2D([], [], color="0.3", ls=":", label="Knee")]
    fig.legend(handles=cells + marks, loc="outside upper center", ncol=4, frameon=False)
    ax[1].legend(handles=[Line2D([], [], color="0.3", marker="s", ms=3, label="Si"),
                          Line2D([], [], color="0.3", marker="s", ms=3, ls="--", alpha=0.55, label="Graphite")],
                 loc="upper left")
    return fig


# ---------------------------------------------------------------------------
# Schematic helpers
# ---------------------------------------------------------------------------
def _box(ax, x, y, w, h, lines, fc="#f4f6f8", ec="#555"):
    """Rounded box in axes coordinates with centred lines [(text, fontsize), ...]."""
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.02",
                                fc=fc, ec=ec, lw=0.7, transform=ax.transAxes, clip_on=False))
    n = len(lines)
    for i, (txt, size) in enumerate(lines):
        ax.text(x + w / 2, y + h * (1 - (i + 0.5) / n), txt, ha="center", va="center", fontsize=size,
                transform=ax.transAxes)


def _arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=8, lw=0.9,
                                 color="#333", transform=ax.transAxes))


def _electrode(ax, sei, isolate=()):
    """Cartoon cross-section: graphite and Si particles in electrolyte, with SEI
    of thickness `sei` (0-1) filling the pores; particles in `isolate` are cut off."""
    ax.add_patch(Rectangle((0, 0), 1, 1, fc=C_EL, ec="none"))
    rng = np.random.default_rng(3)
    pts = [(0.12 + 0.19 * i + 0.04 * (j % 2), 0.15 + 0.21 * j) for j in range(4) for i in range(5)]
    for k, (x, y) in enumerate(pts):
        si = k in (2, 6, 8, 12, 17)
        r = (0.055 if si else 0.075) + 0.004 * rng.standard_normal()
        ax.add_patch(Circle((x, y), r + 0.045 * sei, fc=C_SEI, ec="none", alpha=0.85))
        ax.add_patch(Circle((x, y), r, fc=C_SI if si else C_GR, ec="k", lw=0.4))
        if k in isolate:
            ax.add_patch(Circle((x, y), r, fc="none", ec=C_ISO, lw=1.4))
            for s in (1, -1):
                ax.plot([x - r * 0.6, x + r * 0.6], [y - s * r * 0.6, y + s * r * 0.6], color=C_ISO, lw=1.0)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)


# ---------------------------------------------------------------------------
# Mechanism I schematic
# ---------------------------------------------------------------------------
def schematic_mech1(results_dir, key="25deg_15-95_soc"):
    fig = plt.figure(figsize=(style.DOUBLE_IN, 4.9), layout="constrained")
    sf = fig.subfigures(2, 2, width_ratios=[1.1, 1.0], height_ratios=[1.0, 1.05])

    a0 = sf[0, 0].subplots()
    a0.axis("off")
    sub = a0.inset_axes([0.0, 0.2, 0.48, 0.68])
    sub2 = a0.inset_axes([0.52, 0.2, 0.48, 0.68])
    _electrode(sub, 0.0)
    _electrode(sub2, 1.0, isolate=(6, 7, 12))
    sub.set_title("Early life", fontsize=7.5, pad=2)
    sub2.set_title("After SEI fills the pores", fontsize=7.5, pad=2)
    handles = [Circle((0, 0), 1, fc=C_GR, ec="k", lw=0.4), Circle((0, 0), 1, fc=C_SI, ec="k", lw=0.4),
               Circle((0, 0), 1, fc=C_SEI), Rectangle((0, 0), 1, 1, fc=C_EL),
               Circle((0, 0), 1, fc="none", ec=C_ISO, lw=1.4)]
    a0.legend(handles, ["Graphite", "Si (SiO$_x$)", "SEI", "Electrolyte (pores)", "Isolated"],
              loc="lower center", bbox_to_anchor=(0.5, -0.06), ncol=3, frameon=False, handlelength=1.0)
    style.panel_label(a0, 0)

    a1 = sf[0, 1].subplots()
    h = np.linspace(0, 1, 400)
    for eta, col, lab in ((20, C_SI, r"Si, $\eta$ = 20"), (10, C_GR, r"Graphite, $\eta$ = 10")):
        a1.plot(h, (1 - h) ** eta, color=col, lw=1.5, label=lab)
    a1.set_xlabel(r"Porosity headroom $h$ [-]")
    a1.set_ylabel(r"Gate $g^\ast = (1-h)^\eta$ [-]")
    a1.invert_xaxis()
    a1.legend(loc="upper left")
    top = a1.secondary_xaxis("top", functions=(lambda hh: EPS_MIN + hh * (EPS0 - EPS_MIN),
                                                lambda e: (e - EPS_MIN) / (EPS0 - EPS_MIN)))
    top.set_xlabel(r"Porosity $\varepsilon$ [-]")
    top.tick_params(labelsize=7.5)
    style.panel_label(a1, 1)

    a2 = sf[1, 0].subplots()
    res, _ = an.load(key, results_dir)
    e_age = an.efc(res, res["age_thr"])
    hh = np.clip((res["age_min_porosity"] - EPS_MIN) / (EPS0 - EPS_MIN), 0, 1)
    gate = (1 - hh) ** 20
    dose = np.cumsum(gate)
    m = an.life_mask(res)
    a2.plot(e_age, hh, color="#2c7fb8", lw=1.4, label=r"Headroom $h$")
    a2.plot(e_age, gate, color=C_SI, lw=1.2, ls="--", label=r"Gate $g^\ast$")
    a2.plot(e_age, dose / dose.max(), color="k", lw=1.2, ls="-.", label=r"Dose $D$ (normalised)")
    a2.plot(an.efc(res, res["thr"])[m], res["LAM_si"][m] / 100, color=C_ISO, lw=1.5, label="Si LAM (fraction)")
    a2.set_xlabel(plots.EFC_LABEL)
    a2.set_ylabel("Normalised value [-]")
    a2.legend(loc="lower left", bbox_to_anchor=(0.0, 0.07))
    style.panel_label(a2, 2)

    a3 = sf[1, 1].subplots()
    a3.axis("off")
    _box(a3, 0.0, 0.81, 1.0, 0.17, [("SEI growth (particle surfaces and crack faces)", 7.5),
                                    (r"reduces the porosity $\varepsilon$", 7.5)])
    _arrow(a3, 0.5, 0.81, 0.5, 0.745)
    _box(a3, 0.0, 0.575, 1.0, 0.17, [(r"Gate $g^\ast = (1-h)^\eta$", 7.5),
                                     ("(shut until the pores near their floor)", 7.5)])
    _arrow(a3, 0.5, 0.575, 0.5, 0.51)
    _box(a3, 0.0, 0.34, 1.0, 0.17, [(r"Dose $\mathrm{d}D/\mathrm{d}t = g^\ast$", 7.5),
                                    (r"(after closure, $D \propto \Delta t$)", 7.5)])
    _arrow(a3, 0.5, 0.34, 0.5, 0.275)
    _box(a3, 0.0, 0.0, 1.0, 0.275, [
        (r"Isolation LAM: $\dot\varepsilon_s = k_\mathrm{iso}(T)\,D\,a j_\mathrm{SEI}/F\,(\varepsilon_s/\varepsilon_{s,0})$", 7.5),
        (r"$\Rightarrow$ LAM $\propto \Delta t^2$ (knee), rounded off by Si depletion", 7.5),
        (r"Trapped Li: $\theta\,(c_\mathrm{max}-\bar c)$ per isolated volume", 7.5)], fc="#fdecea", ec=C_ISO)
    style.panel_label(a3, 3)
    return fig


# ---------------------------------------------------------------------------
# Mechanism II schematic
# ---------------------------------------------------------------------------
def _gradient_particle(ax, x, y, r):
    """Route A: a delithiating Si particle with a Li-poor shell and a Li-rich
    core; the shell wants to shrink onto the core, so its surface is in
    tangential (hoop) tension and a surface crack grows."""
    for k, frac in enumerate(np.linspace(1, 0.15, 12)):
        f = k / 11  # 0 at the surface (Li-poor), 1 at the core (Li-rich)
        ax.add_patch(Circle((x, y), r * frac, fc=(1.0, 0.85 - 0.45 * f, 0.55 - 0.5 * f), ec="none"))
    ax.add_patch(Circle((x, y), r, fc="none", ec="k", lw=0.6))
    t = np.deg2rad(60)
    ax.plot([x + r * np.cos(t), x + 0.62 * r * np.cos(t)], [y + r * np.sin(t), y + 0.62 * r * np.sin(t)],
            color="k", lw=1.6)
    for ang in (-20, 140, 230):  # tangential tension arrows along the surface
        t0, t1 = np.deg2rad(ang), np.deg2rad(ang + 32)
        ax.annotate("", xy=(x + 1.16 * r * np.cos(t1), y + 1.16 * r * np.sin(t1)),
                    xytext=(x + 1.16 * r * np.cos(t0), y + 1.16 * r * np.sin(t0)),
                    arrowprops=dict(arrowstyle="<->", color=C_ISO, lw=0.8, mutation_scale=6))
    ax.text(x, y, "Li-rich\ncore", ha="center", va="center", fontsize=7, color="#5a2d00")
    ax.text(x + 1.2 * r, y + 0.95 * r, r"$\sigma_t > 0$", fontsize=7.5, color=C_ISO)


def _strain_cluster(ax, x, y, r):
    """Route B: a Si particle held by binder bridges and graphite neighbours.
    Lithiation swells it (dashed outline); on delithiation it contracts while
    the binder and neighbours hold its surface, so cracks open on contraction."""
    for d, ang in ((1.95, 25), (1.95, 100), (2.0, 175), (1.95, 250), (2.0, 320)):  # graphite neighbours
        t = np.deg2rad(ang)
        cx, cy = x + d * r * np.cos(t), y + d * r * np.sin(t)
        ax.add_patch(Circle((cx, cy), 0.78 * r, fc=C_GR, ec="k", lw=0.4))
        ax.plot([x + r * np.cos(t), cx - 0.78 * r * np.cos(t)], [y + r * np.sin(t), cy - 0.78 * r * np.sin(t)],
                color="#2e8b57", lw=2.4, solid_capstyle="round", alpha=0.9)  # binder bridge
    ax.add_patch(Circle((x, y), 1.32 * r, fc="none", ec=C_SI, lw=1.1, ls="--"))  # swollen (lithiated)
    ax.add_patch(Circle((x, y), r, fc=C_SI, ec="k", lw=0.6))  # contracted (delithiated)
    for ang in (60, 200, 300):  # cracks opening on contraction
        t = np.deg2rad(ang)
        ax.plot([x + r * np.cos(t), x + 0.7 * r * np.cos(t)], [y + r * np.sin(t), y + 0.7 * r * np.sin(t)],
                color="k", lw=1.4)
    for ang in (0, 120, 240):  # contraction arrows
        t = np.deg2rad(ang + 30)
        ax.annotate("", xy=(x + 1.03 * r * np.cos(t), y + 1.03 * r * np.sin(t)),
                    xytext=(x + 1.3 * r * np.cos(t), y + 1.3 * r * np.sin(t)),
                    arrowprops=dict(arrowstyle="-|>", color="#2c7fb8", lw=0.9, mutation_scale=7))


def schematic_mech2(results_dir):
    fig = plt.figure(figsize=(style.DOUBLE_IN, 5.1), layout="constrained")
    sf = fig.subfigures(2, 2, width_ratios=[1.1, 1.0], height_ratios=[1.0, 1.05])

    a0 = sf[0, 0].subplots()
    a0.set_xlim(0, 3.2)
    a0.set_ylim(-0.78, 1.38)
    a0.set_aspect("equal")
    a0.axis("off")
    _gradient_particle(a0, 0.62, 0.66, 0.36)
    _strain_cluster(a0, 2.25, 0.5, 0.26)
    a0.text(0.62, -0.76, "A: Paris law\nhoop tension from the\nconcentration gradient",
            ha="center", va="bottom", fontsize=7)
    a0.text(2.25, -0.76, "B: volume-strain fatigue\nSi swells and contracts against\nbinder and neighbours",
            ha="center", va="bottom", fontsize=7)
    handles = [Circle((0, 0), 1, fc=C_SI, ec="k", lw=0.4), Circle((0, 0), 1, fc=C_GR, ec="k", lw=0.4),
               Line2D([], [], color="#2e8b57", lw=2.4), Line2D([], [], color=C_SI, lw=1.1, ls="--")]
    a0.legend(handles, ["Si", "Graphite", "Binder", "Swollen Si"], loc="upper center", frameon=False, ncol=4,
              bbox_to_anchor=(0.5, 1.04), handlelength=1.1, columnspacing=0.7)
    style.panel_label(a0, 0)

    a1 = sf[0, 1].subplots()
    s = np.linspace(1e-4, 1, 400)
    for n, ls in ((1.56, "-"), (2.19, "--")):
        v = 1 + 3 * s ** n
        a1.plot(s, np.log(v), color=C_SI, lw=1.4, ls=ls)
        a1.plot(s, 3 * n * s ** (n - 1) / v, color="#2c7fb8", lw=1.1, ls=ls)
    a1.set_xlabel(r"Si stoichiometry $\bar s$ [-]")
    a1.set_ylabel(r"$\varepsilon_V$, d$\varepsilon_V$/d$\bar s$ [-]")
    a1.legend(handles=[Line2D([], [], color=C_SI, label=r"$\varepsilon_V = \ln(1+3\bar s^{\,n})$"),
                       Line2D([], [], color="#2c7fb8", label=r"d$\varepsilon_V$/d$\bar s$"),
                       Line2D([], [], color="0.3", label=r"$n$ = 1.56 (beginning of life)"),
                       Line2D([], [], color="0.3", ls="--", label=r"$n$ = 2.19 (aged)")], loc="upper left", ncol=2)
    a1.set_ylim(0, 2.45)
    style.panel_label(a1, 1)

    a2 = sf[1, 0].subplots()
    for key in ("25deg_full_soc", "25deg_15-95_soc", "25deg_20-80_soc"):
        res, _ = an.load(key, results_dir)
        m = an.life_mask(res)
        e = an.efc(res, res["thr"])[m]
        a2.plot(e, 1e9 * res["dl_paris"][m], color=CELL_COLOR[key], lw=1.4)
        a2.plot(e, 1e9 * res["dl_strain"][m], color=CELL_COLOR[key], lw=1.4, ls="--")
    a2.set_xlabel(plots.EFC_LABEL)
    a2.set_ylabel("Cumulative crack growth [nm]")
    a2.legend(handles=[Line2D([], [], color=CELL_COLOR[k], label=ablations.cond_tag(CONDITIONS[k]))
                       for k in ("25deg_full_soc", "25deg_15-95_soc", "25deg_20-80_soc")]
              + [Line2D([], [], color="0.3", label="Route A (Paris)"),
                 Line2D([], [], color="0.3", ls="--", label="Route B (strain)")], loc="upper left")
    style.panel_label(a2, 2)

    a3 = sf[1, 1].subplots()
    a3.axis("off")
    _box(a3, 0.0, 0.83, 1.0, 0.16, [("Route A: Paris law (tension only)", 7.5),
                                    (r"$\mathrm{d}l/\mathrm{d}t = k_{cr}(T)\,[\sigma_t(1-l/R)\,b_{cr}\sqrt{\pi l}\,]^{m_{cr}}$", 7.5)])
    a3.text(0.5, 0.807, "+", ha="center", va="center", fontsize=8, transform=a3.transAxes)
    _box(a3, 0.0, 0.555, 1.0, 0.225, [("Route B: volume-strain fatigue (contraction only)", 7.5),
                                      (r"$\mathrm{d}l/\mathrm{d}t = k_V(T)\,\max(-\dot\varepsilon_V, 0)\,l\,(1-l/R)$", 7.5),
                                      (r"Miner-linear: $\Delta\ln l \approx k_V\,\Delta\varepsilon_V$ per contraction", 7.5)])
    _arrow(a3, 0.5, 0.555, 0.5, 0.505)
    _box(a3, 0.06, 0.395, 0.88, 0.11, [(r"Crack length $l$ $\Rightarrow$ roughness $\rho = 1 + 2\,l\,w_{cr}\rho_{cr}$", 7.5)])
    _arrow(a3, 0.5, 0.395, 0.5, 0.345)
    _box(a3, 0.06, 0.235, 0.88, 0.11, [("SEI on the new crack faces: LLI and pore filling", 7.5)])
    _arrow(a3, 0.5, 0.235, 0.5, 0.185)
    _box(a3, 0.06, 0.075, 0.88, 0.11, [(r"Pores close $\Rightarrow$ isolation (mechanism I)", 7.5)], fc="#fdecea", ec=C_ISO)
    a3.text(0.5, 0.0, r"Both routes slower when hot: $E_{cr}$ = +90, $E_V$ = +150 kJ mol$^{-1}$",
            ha="center", va="bottom", fontsize=7, transform=a3.transAxes, color="#444")
    style.panel_label(a3, 3)
    return fig


# ---------------------------------------------------------------------------
# Discussion: partial windows late in life
# ---------------------------------------------------------------------------
def partial_windows(data_dir, results_dir):
    """Why the 20-80 % cell collapses through route A late in life while the
    15-95 % cell does not (fixed-Ah vs voltage-limited ageing discharge)."""
    cells = (("25deg_20-80_soc", "fixed 1.5 A h discharge"), ("25deg_15-95_soc", "voltage-limited, 3.14 V"))
    fig = plt.figure(figsize=(style.DOUBLE_IN, 4.6), layout="constrained")
    gs = fig.add_gridspec(2, 6)
    ax = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]), fig.add_subplot(gs[0, 4:6]),
          fig.add_subplot(gs[1, 0:3]), fig.add_subplot(gs[1, 3:6])]
    for key, how in cells:
        cond, col = CONDITIONS[key], CELL_COLOR[key]
        r, _ = an.load(key, results_dir)
        ax[0].plot(an.efc(r, r["rpt_thr"]), an.rpt_soh(r), "-^", color=col, ms=3.5, lw=1.1)
        cap = data.load(cond, data_dir, "capacity_fade")
        c20 = cap[cap.source == "C/20 RPT"]
        ax[0].plot(c20.efc, c20.capacity_retention_pct, "o", mfc="white", mec=col, mew=1.1, ms=4.5, ls="none")
        ax[1].plot(an.efc(r, r["age_thr"]), r["age_cap"], color=col, lw=1.2)
        ax[2].plot(an.efc(r, r["cycle_thr"]), r["cycle_peak_stress_MPa"], ".", ms=1.6, color=col)
        m = an.life_mask(r)
        e = an.efc(r, r["thr"])[m]
        ax[3].plot(e, 1e9 * r["dl_paris"][m], color=col, lw=1.4)
        ax[3].plot(e, 1e9 * r["dl_strain"][m], color=col, lw=1.4, ls="--")
        ax[4].plot(e, r["LAM_si"][m], color=col, lw=1.4)
    ax[1].axhline(1.5, color="k", ls=":", lw=0.8)
    ax[1].text(5, 1.52, "1.5 A h set-point", fontsize=7)
    ax[2].set_ylim(0, 190)
    for a, yl in zip(ax, ("SoH, C/20 RPT [%]", "C/3 ageing discharge [A h]", "Peak Si surface\ntangential stress [MPa]",
                          "Cumulative crack growth [nm]", "Si LAM [%]")):
        a.set_ylabel(yl)
        a.set_xlabel(plots.EFC_LABEL)
    for i, a in enumerate(ax):
        style.panel_label(a, i)
    fig.legend(handles=[Line2D([], [], color=CELL_COLOR[k], lw=1.4, label=f"{ablations.cond_tag(CONDITIONS[k])} ({how})")
                        for k, how in cells], loc="outside upper center", ncol=2, frameon=False)
    ax[0].legend(handles=[Line2D([], [], color="0.3", marker="^", ms=3.5, label="Model"),
                          Line2D([], [], color="0.3", ls="none", marker="o", ms=4.5, mfc="white", label="Measured")],
                 loc="lower left")
    ax[3].legend(handles=[Line2D([], [], color="0.3", label="Route A (Paris)"),
                          Line2D([], [], color="0.3", ls="--", label="Route B (strain)")], loc="upper left")
    return fig


# ---------------------------------------------------------------------------
def make_all(data_dir, results_dir, ablations_dir, out_dir, log=print):
    """Write every paper figure (PNG at 300 dpi and PDF, printed size) to out_dir.
    The ablation figure is skipped when ablations_dir has no results."""
    style.use()
    for key, cond in CONDITIONS.items():
        res, curves = an.load(key, results_dir)
        for p in plots.save_all(res, curves, cond, data_dir, out_dir).values():
            log("written", p)
    figs = {"data_overview": lambda: data_overview(data_dir),
            "schematic_mech1": lambda: schematic_mech1(results_dir),
            "schematic_mech2": lambda: schematic_mech2(results_dir),
            "discussion_partial_windows": lambda: partial_windows(data_dir, results_dir)}
    if ablations_dir and os.path.isdir(ablations_dir):
        figs["ablations"] = lambda: ablations.figure(data_dir, ablations_dir, results_dir)
    for name, make in figs.items():
        fig = make()
        log("written", style.save(fig, out_dir, name)[0])
        plt.close(fig)
