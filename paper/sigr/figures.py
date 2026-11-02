"""The paper's figures, all at their printed JES size (see style.py):

    import sigr.figures as F
    F.make_all("data", "results", "outputs_ablations", "paper_pics")

Per condition (plots.py): overall summary, C/20 RPT voltage and dV/dSoC
against depth of discharge, Si crack/porosity diagnostics. Shared: the
measured-data overview, the two mechanism schematics, the ablations
(ablations.py) and the partial-window comparison of the Discussion. The
figures with DMA data go to one subfolder per DMA (make_all's dma).
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
    for eta, col, lab in ((20, C_SI, r"Si, $\omega$ = 20"), (10, C_GR, r"Graphite, $\omega$ = 10")):
        a1.plot(h, (1 - h) ** eta, color=col, lw=1.5, label=lab)
    a1.set_xlabel(r"Porosity headroom $h$ [-]")
    a1.set_ylabel(r"Gate $g^\ast = (1-h)^\omega$ [-]")
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
    a2.plot(e_age, dose / dose.max(), color="k", lw=1.2, ls="-.", label=r"Dose $\Gamma$ (normalised)")
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
    _box(a3, 0.0, 0.575, 1.0, 0.17, [(r"Gate $g^\ast = (1-h)^\omega$", 7.5),
                                     ("(shut until the pores near their floor)", 7.5)])
    _arrow(a3, 0.5, 0.575, 0.5, 0.51)
    _box(a3, 0.0, 0.34, 1.0, 0.17, [(r"Dose $\mathrm{d}\Gamma/\mathrm{d}t = g^\ast$", 7.5),
                                    (r"(after closure, $\Gamma \propto \Delta t$)", 7.5)])
    _arrow(a3, 0.5, 0.34, 0.5, 0.275)
    _box(a3, 0.0, 0.0, 1.0, 0.275, [
        (r"Isolation LAM: $\dot\varepsilon_s = k_\mathrm{iso}(T)\,\Gamma\,a j_\mathrm{SEI}/F\,(\varepsilon_s/\varepsilon_{s,0})$", 7.5),
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
    a1.legend(handles=[Line2D([], [], color=C_SI, label=r"$\varepsilon_V = \ln(1+3\bar s^{\,\lambda})$"),
                       Line2D([], [], color="#2c7fb8", label=r"d$\varepsilon_V$/d$\bar s$"),
                       Line2D([], [], color="0.3", label=r"$\lambda$ = 1.56 (beginning of life)"),
                       Line2D([], [], color="0.3", ls="--", label=r"$\lambda$ = 2.19 (aged)")], loc="upper left", ncol=2)
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
                                    (r"$\mathrm{d}\ell \,/\,\mathrm{d}t = k_{cr}(T)\,[\sigma_t(1-\ell \,/\,R_\mathrm{Si})\,b_{cr}\sqrt{\pi \ell}\,]^{m_{cr}}$", 7.5)])
    a3.text(0.5, 0.8035, "+", ha="center", va="center_baseline", fontsize=8, transform=a3.transAxes)
    _box(a3, 0.0, 0.555, 1.0, 0.225, [("Route B: volume-strain fatigue (contraction only)", 7.5),
                                      (r"$\mathrm{d}\ell \,/\,\mathrm{d}t = k_V(T)\,\max(-\dot\varepsilon_V, 0)\,\ell \,(1-\ell \,/\,R_\mathrm{Si})$", 7.5),
                                      (r"Miner-linear: $\Delta\ln \ell  \approx k_V\,\Delta\varepsilon_V$ per contraction", 7.5)])
    _arrow(a3, 0.5, 0.555, 0.5, 0.505)
    _box(a3, 0.06, 0.395, 0.88, 0.11, [(r"Crack length $\ell $ $\Rightarrow$ roughness $\Lambda = 1 + 2\,\ell \,w_{cr}N_{cr}$", 7.5)])
    _arrow(a3, 0.5, 0.395, 0.5, 0.345)
    _box(a3, 0.06, 0.235, 0.88, 0.11, [("SEI on the new crack faces: LLI and pore filling", 7.5)])
    _arrow(a3, 0.5, 0.235, 0.5, 0.185)
    _box(a3, 0.06, 0.075, 0.88, 0.11, [(r"Pores close $\Rightarrow$ isolation (mechanism I)", 7.5)], fc="#fdecea", ec=C_ISO)
    a3.text(0.5, 0.0, r"Both routes slower when hot: $E_{a,cr}$ = +90, $E_{a,V}$ = +150 kJ mol$^{-1}$",
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
# Results: mechanism chain, Paris-led against strain-led cracking
# ---------------------------------------------------------------------------
def crack_comparison(data_dir, results_dir, keys=("25deg_full_soc", "25deg_15-95_soc")):
    """Si cracking and pore closure of two cells overlaid (one colour per cell): crack length, roughness, peak
    surface tangential stress of each cycle, per-cycle minimum porosity, Si SEI LLI (surface and crack faces),
    cumulative contraction strain, crack growth by route, and Si LAM against the DMA."""
    fig, axs = plt.subplots(4, 2, figsize=(style.DOUBLE_IN, 6.6), sharex=True, layout="constrained")
    ax = axs.flat
    to_pct = 100 / 2.5947
    for key in keys:
        cond, col = CONDITIONS[key], CELL_COLOR[key]
        res, _ = an.load(key, results_dir)
        e = an.efc(res, res["thr"])
        m = an.life_mask(res)
        ax[0].plot(e, res["crack_length"] / res["crack_length"][0], color=col)
        ax[1].plot(e, res["roughness"], color=col)
        cyc = np.isin(res["cycle_thr"], res["age_thr"])  # ageing cycles only (formation and RPTs left out)
        ax[2].plot(an.efc(res, res["cycle_thr"])[cyc], res["cycle_peak_stress_MPa"][cyc], ".", ms=1.5, color=col)
        # per-cycle minimum porosity; the first cycle after an RPT starts from a full charge (left out for partial)
        after_rpt = np.searchsorted(res["age_thr"], res["rpt_thr"], side="right")
        keep = np.ones(res["age_thr"].size, bool)
        keep[after_rpt[after_rpt < keep.size]] = not cond.partial
        ax[3].plot(np.r_[0.0, an.efc(res, res["age_thr"])[keep]],
                   np.r_[res["initial_porosity"], res["age_min_porosity"][keep]], color=col)
        ax[4].plot(e, res["Q_si_sei_cracks"] * to_pct, color=col)
        ax[4].plot(e, res["Q_si_sei"] * to_pct, color=col, ls="--")
        ax[5].plot(e, res["strain_E"], color=col)
        ax[6].plot(e, 1e9 * res["dl_paris"], color=col)
        ax[6].plot(e, 1e9 * res["dl_strain"], color=col, ls="--")
        ax[7].plot(e[m], res["LAM_si"][m], color=col)
        lam = data.load(cond, data_dir, "LAM_derived")
        ax[7].plot(lam["efc"], lam["LAM_NE_silicon_pct"], "s", ms=3.5, color=col)
    ax[3].axhline(EPS_MIN, color="gray", ls=":", lw=0.8)
    for a, lab in zip(ax, ("Si crack length / initial [-]", "Si roughness ratio [-]",
                           "Peak Si surface tangential\nstress, ageing cycles [MPa]",
                           "Negative electrode porosity,\nper-cycle minimum [-]", "LLI to Si SEI [%]",
                           "Cumulative Si\ncontraction strain [-]", "Cumulative crack growth [nm]", "Si LAM [%]")):
        a.set_ylabel(lab)
    for a in axs[-1]:
        a.set_xlabel(plots.EFC_LABEL)
    ax[4].legend(handles=[Line2D([], [], color="0.3", label="Crack faces"),
                          Line2D([], [], color="0.3", ls="--", label="Particle surface")], loc="upper left")
    ax[6].legend(handles=[Line2D([], [], color="0.3", label="Route A (Paris, diffusion stress)"),
                          Line2D([], [], color="0.3", ls="--", label="Route B (volume-strain fatigue)")],
                 loc="upper left")
    ax[7].legend(handles=[Line2D([], [], color="0.3", label="Model"),
                          Line2D([], [], color="0.3", ls="none", marker="s", ms=3.5, label="Measured (DMA)")],
                 loc="upper left")
    for i, a in enumerate(ax):
        style.panel_label(a, i)
    fig.legend(handles=[Line2D([], [], color=CELL_COLOR[k], lw=1.4, label=ablations.cond_tag(CONDITIONS[k]))
                        for k in keys], loc="outside upper center", ncol=len(keys), frameon=False)
    return fig


# ---------------------------------------------------------------------------
# Results: C/20 voltage and dV/dSoC at three RPTs of each cell
# ---------------------------------------------------------------------------
RPT_ROLES = (("RPT1", "black"), ("last before the knee", "tab:blue"), ("last RPT", "tab:red"))
DVDSOC_FROM = 0.03  # dV/dSoC drawn from 3 % depth of discharge (step-start transient, edge of the binning)


def rpt_selection(data_dir, results_dir, trim=0.01):
    """One row per cell, against depth of discharge at three measured RPTs (RPT1, the last before the measured
    knee, the last): the C/20 discharge voltage of all three (first column; measured solid, model dashed), and the
    dV/dSoC of each in its own panel (measured black, model dashed in the RPT's colour). Returns
    (figure, {cell: [(rpt, EFC)]})."""
    keys = list(CONDITIONS)
    fig, axs = plt.subplots(len(keys), 4, figsize=(style.DOUBLE_IN, 7.8), sharex=True, layout="constrained",
                            gridspec_kw=dict(width_ratios=[1.25, 1, 1, 1]))
    chosen = {}
    for r_, key in enumerate(keys):
        cond = CONDITIONS[key]
        res, curves = an.load(key, results_dir)
        pairs = an.match_rpt_curves(res, curves, cond, data_dir)
        knee = data.knee_efc(cond, data_dir)
        before = [p for p in pairs if p[1] < knee]
        sel = [pairs[0], before[-1], pairs[-1]]
        chosen[cond.cell] = [(r, e) for r, e, _, _ in sel]
        av = axs[r_, 0]
        for c, ((r, e_r, j, _), (_, col)) in enumerate(zip(sel, RPT_ROLES)):
            q_r, v_r = data.load_rpt_discharge(cond, data_dir, r)
            s_r, s_m, v_m = q_r / q_r.max(), curves[j]["q"] / curves[j]["q"].max(), curves[j]["v"]
            keep = s_m >= trim  # the model's first ~1 % is the step-start transient
            av.plot(s_r, v_r, color=col, lw=1.2, label=f"RPT{r}, {e_r:.0f} EFC")
            av.plot(s_m[keep], v_m[keep], "--", color=col, lw=1.0)
            ad = axs[r_, c + 1]
            for (x, d), kw in ((an.smooth_dvdsoc(s_r, v_r), dict(color="black", lw=1.2)),
                               (an.smooth_dvdsoc(s_m[keep], v_m[keep]), dict(color=col, lw=1.1, ls="--"))):
                ad.plot(x[x >= DVDSOC_FROM], d[x >= DVDSOC_FROM], **kw)
            ad.set_ylim(-2.6, 0)
            if c:
                ad.tick_params(labelleft=False)
            style.panel_label(ad, 4 * r_ + c + 1, f"RPT{r}")
        av.set_ylabel("Terminal voltage [V]")
        axs[r_, 1].set_ylabel("dV/dSoC [V]")
        av.legend(loc="lower left", title=ablations.cond_tag(cond), title_fontsize=7, fontsize=7)
        style.panel_label(av, 4 * r_)
    for a in axs[-1]:
        a.set_xlabel("Depth of discharge [-]")
    fig.legend(handles=[Line2D([], [], color="0.3", lw=1.2, label="Measured (dV/dSoC: black)"),
                        Line2D([], [], color="0.3", ls="--", lw=1.1, label="Model")]
               + [Line2D([], [], color=c, lw=1.2, label=lab) for lab, c in RPT_ROLES],
               loc="outside upper center", ncol=5, frameon=False)
    return fig, chosen


# ---------------------------------------------------------------------------
# Graphical abstract (ECS: a single image, uploaded separately as GA*.png)
# ---------------------------------------------------------------------------
def graphical_abstract(data_dir, results_dir, key_v="25deg_full_soc", trim=0.01):
    """Top: the mechanism chain with its key equations. Bottom: the four cells with one parameter set (C/20 SoH and
    Si LAM, model against measured), and the C/20 voltage and dV/dSoC of one cell (key_v) at the three RPTs of
    Fig. 9 (RPT_ROLES)."""
    fig = plt.figure(figsize=(style.DOUBLE_IN, 5.0))
    ax = fig.add_axes((0.015, 0.555, 0.97, 0.36))
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    route_a, route_b, sei, iso, knee = "#dce9f7", "#e0f2e0", "#f3ead9", "#f8dcdc", "#eeeeee"
    _box(ax, 0.005, 0.56, 0.255, 0.43, [("A: Paris cracking", 8), ("(full window, diffusion stress)", 7),
                                        (r"$\dot\ell_A \propto \left(\sigma_t\sqrt{\pi\ell}\right)^{m_{cr}}$", 8)],
         fc=route_a)
    _box(ax, 0.005, 0.01, 0.255, 0.43, [("B: volume-strain fatigue", 8), ("(partial windows, contraction)", 7),
                                        (r"$\dot\ell_B = k_V\,\max(-\dot\varepsilon_V,0)\,\ell\,(1-\ell/R)$", 7.5)],
         fc=route_b)
    _box(ax, 0.305, 0.17, 0.235, 0.66, [("SEI on fresh crack faces", 8), ("fills the pores", 8),
                                        (r"$\Lambda = 1 + 2\,\ell\,w_{cr}N_{cr}$", 8),
                                        (r"$\varepsilon = \varepsilon_0 - a\left[\Delta L_{\mathrm{SEI}}"
                                         r" + (\Lambda-1)\,L_{\mathrm{SEI,cr}}\right]$", 7.5)], fc=sei)
    _box(ax, 0.585, 0.17, 0.19, 0.66, [("Pores close:", 8), ("the gate opens", 8),
                                       (r"$h = \frac{\varepsilon-\varepsilon_{\min}}{\varepsilon_0-\varepsilon_{\min}}$", 8.5),
                                       (r"$g^\ast = (1-h)^{\omega}$", 8)], fc=iso)
    _box(ax, 0.82, 0.17, 0.175, 0.66, [("Si isolated with", 8), ("its lithium:", 8), ("capacity knee", 8),
                                       (r"$\dot\Gamma = g^\ast$", 8),
                                       (r"$\dot\varepsilon_{s,\mathrm{iso}} = k_{\mathrm{iso}}\,\Gamma\,\frac{a\,j_{\mathrm{SEI}}}{F}$", 8)],
         fc=knee)
    _arrow(ax, 0.268, 0.775, 0.3, 0.6)
    _arrow(ax, 0.268, 0.225, 0.3, 0.4)
    _arrow(ax, 0.548, 0.5, 0.578, 0.5)
    _arrow(ax, 0.783, 0.5, 0.813, 0.5)
    fig.text(0.5, 0.985, "Two silicon cracking routes, one isolation gate", ha="center", va="top", fontsize=9,
             fontweight="bold")
    gs = fig.add_gridspec(2, 2, left=0.085, right=0.985, top=0.44, bottom=0.14, wspace=0.3, hspace=0.12)
    ax_s = fig.add_subplot(gs[0, 0])
    ax_l = fig.add_subplot(gs[1, 0], sharex=ax_s)
    ax_v = fig.add_subplot(gs[0, 1])
    ax_d = fig.add_subplot(gs[1, 1], sharex=ax_v)
    for key, cond in CONDITIONS.items():
        col = CELL_COLOR[key]
        res, _ = an.load(key, results_dir)
        ax_s.plot(an.efc(res, res["rpt_thr"]), an.rpt_soh(res), "-", color=col, lw=1.2)
        cap = data.load(cond, data_dir, "capacity_fade")
        c20 = cap[cap["source"] == "C/20 RPT"]
        ax_s.plot(c20["efc"], c20["capacity_retention_pct"], "o", ms=3.5, mfc="white", mec=col, mew=1.0)
        m = an.life_mask(res)
        ax_l.plot(an.efc(res, res["thr"])[m], res["LAM_si"][m], "-", color=col, lw=1.2)
        lam = data.load(cond, data_dir, "LAM_derived")
        ax_l.plot(lam["efc"], lam["LAM_NE_silicon_pct"], "s", ms=3.2, mfc="white", mec=col, mew=1.0)
    ax_s.set_ylabel("SoH, C/20 [%]")
    ax_l.set_ylabel("Si LAM [%]")
    ax_l.set_xlabel(plots.EFC_LABEL)
    ax_s.tick_params(labelbottom=False)
    ax_l.legend(handles=[Line2D([], [], color="0.3", lw=1.2, label="Model"),
                         Line2D([], [], color="0.3", ls="none", marker="s", ms=3.2, mfc="white",
                                label="Measured (DMA)")], loc="upper left", fontsize=7)
    cond = CONDITIONS[key_v]
    res, curves = an.load(key_v, results_dir)
    pairs = an.match_rpt_curves(res, curves, cond, data_dir)
    knee = data.knee_efc(cond, data_dir)
    sel = [pairs[0], [p for p in pairs if p[1] < knee][-1], pairs[-1]]  # as Fig. 9: RPT1, last before the knee, last
    for (r, e_r, j, _), (_, col) in zip(sel, RPT_ROLES):
        q_r, v_r = data.load_rpt_discharge(cond, data_dir, r)
        s_r, s_m, v_m = q_r / q_r.max(), curves[j]["q"] / curves[j]["q"].max(), curves[j]["v"]
        keep = s_m >= trim  # the model's first ~1 % is the step-start transient
        ax_v.plot(s_r, v_r, color=col, lw=1.2, label=f"RPT{r}, {e_r:.0f} EFC")
        ax_v.plot(s_m[keep], v_m[keep], "--", color=col, lw=1.0)
        for (x, d), ls in ((an.smooth_dvdsoc(s_r, v_r), "-"), (an.smooth_dvdsoc(s_m[keep], v_m[keep]), "--")):
            ax_d.plot(x[x >= DVDSOC_FROM], d[x >= DVDSOC_FROM], color=col, lw=1.2 if ls == "-" else 1.0, ls=ls)
    ax_d.set_ylim(-2.6, 0)
    ax_v.set_ylabel("Voltage, C/20 [V]")
    ax_d.set_ylabel("dV/dSoC [V]")
    ax_d.set_xlabel("Depth of discharge [-]")
    ax_v.tick_params(labelbottom=False)
    ax_v.legend(loc="lower left", fontsize=7, handlelength=1.4)
    ax_d.legend(handles=[Line2D([], [], color="0.3", lw=1.2, label="Measured"),
                         Line2D([], [], color="0.3", ls="--", lw=1.0, label="Model")], loc="lower left", fontsize=7)
    fig.text(0.285, 0.475, "One parameter set, four cells", ha="center", va="bottom", fontsize=9, fontweight="bold")
    fig.text(0.79, 0.475, ablations.cond_tag(cond).split(", ", 1)[1] + ": voltage and dV/dSoC", ha="center",
             va="bottom", fontsize=9, fontweight="bold")
    fig.legend(handles=[Line2D([], [], color=CELL_COLOR[k], lw=1.2, label=ablations.cond_tag(c).split(", ", 1)[1])
                        for k, c in CONDITIONS.items()],
               loc="lower center", bbox_to_anchor=((ax_l.get_position().x0 + ax_l.get_position().x1) / 2, 0.0),
               ncol=2, fontsize=7.5, frameon=False, handlelength=1.4, columnspacing=1.5)
    return fig


# ---------------------------------------------------------------------------
# Parameterisation: beginning-of-life fit
# ---------------------------------------------------------------------------
BOL_DVDQ_WINDOW = {"discharge": (0.02, 0.95), "charge": (0.05, 0.98)}  # normalised capacity scored


def bol_fit(data_dir, results_dir, key="25deg_full_soc", rpt=1, trim=0.01):
    """Model against measured C/20 discharge (left) and constant-current C/20 charge (right) at the reference
    RPT of CELL064: voltage (top) and dV/dQ (bottom) against capacity. The model curves are the formation of
    the paper run (simulate.bol_curves), cached in results_dir. Returns (figure, {kind: (RMSE [mV], model
    capacity [A h], measured capacity [A h], dV/dQ RMSE [V/(A h)])}). The RMSEs, shown in the panels, are
    taken on a common grid of the capacity normalised by each curve's own; dV/dQ away from the steep end of
    each step (BOL_DVDQ_WINDOW)."""
    from . import simulate  # pybamm import only when this figure is made
    cond = CONDITIONS[key]
    cache = os.path.join(results_dir, f"{key}_bol_curves.npz")
    if not os.path.exists(cache):
        c = simulate.bol_curves(cond, data_dir)
        np.savez(cache, **{f"{k}_{x}": v for k, (q, vv) in c.items() for x, v in (("q", q), ("v", vv))})
    z = np.load(cache)
    meas = {"discharge": data.load_rpt_discharge(cond, data_dir, rpt), "charge": data.load_rpt_charge(cond, data_dir, rpt)}
    fig, ax = plt.subplots(2, 2, figsize=(style.DOUBLE_IN, 3.9), sharex="col", layout="constrained")
    stats = {}
    for c, kind in enumerate(("discharge", "charge")):
        q_r, v_r = meas[kind]
        q_m, v_m = z[f"{kind}_q"], z[f"{kind}_v"]
        keep = q_m >= trim * q_m.max()  # the model's first ~1 % is the step-start transient
        q_m, v_m = q_m[keep], v_m[keep]
        ax[0, c].plot(q_r, v_r, color="black", lw=1.3, label="Measured")
        ax[0, c].plot(q_m, v_m, "--", color=plots.C_MODEL, lw=1.2, label="Model")
        dvdq = []
        for q, v, kw in ((q_r, v_r, dict(color="black", lw=1.3)), (q_m, v_m, dict(color=plots.C_MODEL, lw=1.2, ls="--"))):
            x, d = an.smooth_dvdsoc(q / q.max(), v)  # binned on a uniform grid, then differentiated
            ax[1, c].plot(x * q.max(), d / q.max(), **kw)
            dvdq.append((x, d / q.max()))
        ax[1, c].set_xlabel(f"{kind.capitalize()} capacity, C/20 [A h]")
        x = np.linspace(0.02, 0.98, 200)
        rmse_v = np.sqrt(np.mean((np.interp(x, q_m / q_m.max(), v_m) - np.interp(x, q_r / q_r.max(), v_r)) ** 2))
        (x_r, d_r), (x_m, d_m) = dvdq
        lo, hi = BOL_DVDQ_WINDOW[kind]
        sel = (x_r >= max(lo, x_m.min())) & (x_r <= min(hi, x_m.max()))
        rmse_d = np.sqrt(np.mean((np.interp(x_r[sel], x_m, d_m) - d_r[sel]) ** 2))
        stats[kind] = (1e3 * rmse_v, float(q_m.max()), float(q_r.max()), rmse_d)
        where = dict(discharge=(0.97, 0.95, "right", "top"), charge=(0.03, 0.95, "left", "top"))[kind]
        ax[0, c].text(*where[:2], f"RMSE {1e3 * rmse_v:.0f} mV", transform=ax[0, c].transAxes, ha=where[2], va=where[3])
        ax[1, c].text(0.03 if kind == "discharge" else 0.97, 0.95, f"RMSE {rmse_d:.3f} V A$^{{-1}}$ h$^{{-1}}$",
                      transform=ax[1, c].transAxes, ha="left" if kind == "discharge" else "right", va="top")
    ax[1, 0].set_ylim(-1.0, 0.1)
    ax[1, 1].set_ylim(0.0, 1.1)
    ax[0, 0].set_ylabel("Terminal voltage [V]")
    ax[1, 0].set_ylabel("dV/dQ [V A$^{-1}$ h$^{-1}$]")
    ax[0, 0].legend(loc="lower left")
    for i, a in enumerate(ax.flat):
        style.panel_label(a, i)
    return fig, stats


# ---------------------------------------------------------------------------
DMA_FIGURES = ("overall_summary",)  # per condition; with data_overview and ablations, the figures with DMA data


def make_all(data_dir, results_dir, ablations_dir, out_dir, log=print, dma=None):
    """Write every paper figure (PNG at 300 dpi and PDF, printed size). Figures without DMA data go to
    out_dir; those with it (overall summaries, data overview, ablations) to out_dir/<folder> for each
    {folder: data directory with that DMA's LAM/LLI tables} in dma (default {"dma_ours": data_dir}).
    The ablation figure is skipped when ablations_dir has no results."""
    style.use()
    dma = dma or {"dma_ours": data_dir}
    have_ablations = ablations_dir and os.path.isdir(ablations_dir)
    for key, cond in CONDITIONS.items():
        res, curves = an.load(key, results_dir)
        for p in plots.save_all(res, curves, cond, data_dir, out_dir, exclude=DMA_FIGURES).values():
            log("written", p)
        for folder, d in dma.items():
            for p in plots.save_all(res, curves, cond, d, os.path.join(out_dir, folder), only=DMA_FIGURES).values():
                log("written", p)
    figs = {"schematic_mech1": lambda: schematic_mech1(results_dir),
            "schematic_mech2": lambda: schematic_mech2(results_dir),
            "discussion_partial_windows": lambda: partial_windows(data_dir, results_dir),
            "bol_fit": lambda: bol_fit(data_dir, results_dir)[0],
            "crack_comparison": lambda: crack_comparison(data_dir, results_dir)}
    for folder, d in dma.items():
        figs[os.path.join(folder, "data_overview")] = lambda d=d: data_overview(d)
        if have_ablations:
            figs[os.path.join(folder, "ablations")] = lambda d=d: ablations.figure(d, ablations_dir, results_dir)
    for name, make in figs.items():
        fig = make()
        log("written", style.save(fig, os.path.join(out_dir, os.path.dirname(name)), os.path.basename(name))[0])
        plt.close(fig)
