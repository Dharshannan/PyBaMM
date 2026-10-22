"""Fixed-parameter ablations: switch one mechanism off, keep every other
parameter at its fitted (v7) value, and rerun the affected conditions.

    import sigr.ablations as ab
    ab.run_ablations(data_dir, "outputs_ablations")       # all of them, in parallel
    res, curves = sigr.analysis.load("25deg_15-95_soc", "outputs_ablations/no_strain_cracking")
"""
import os

from . import _run_jobs
from .conditions import EFC_BASIS_AH, RATE_BASIS_AH

ABLATIONS = {
    "no_isolation": dict(
        # Rate constant set to zero rather than removing "porosity isolation" from the
        # LAM option: the Si volume-law exponents and OCP aging deformation are scoped
        # to that option in the core model, so removing it would change them too.
        description="Mechanism I off: isolation rate k_iso = 0 for Si and Gr (so no lithium trapping either)",
        keys=["25deg_full_soc", "45deg_full_soc", "25deg_15-95_soc", "25deg_20-80_soc"],
        params={"si_iso_Y0": 0.0, "gr_iso_Y0": 0.0},
    ),
    "no_strain_cracking": dict(
        description="Mechanism II off: Paris-law cracking only (k_V = 0)",
        keys=["25deg_full_soc", "25deg_15-95_soc", "25deg_20-80_soc"],
        params={"si_strain_kV": 0.0},
    ),
    "strain_cracking_no_temperature": dict(
        description="Mechanism II without its temperature dependence (E_V = 0)",
        keys=["45deg_full_soc"],
        params={"si_strain_Ea": 0.0},
    ),
    "no_li_trapping": dict(
        description="No extra lithium trapped in isolated Si (theta = 0)",
        keys=["25deg_full_soc", "25deg_15-95_soc"],
        params={"si_iso_li_trapping": 0.0},
    ),
}


def jobs(data_dir, out_dir, names=None, rate_basis_ah=RATE_BASIS_AH, efc_basis_ah=EFC_BASIS_AH):
    out = []
    for name in names or ABLATIONS:
        a = ABLATIONS[name]
        for key in a["keys"]:
            out.append(dict(key=key, label=f"{name}/{key}", data_dir=data_dir, out_dir=os.path.join(out_dir, name),
                            rate_basis_ah=rate_basis_ah, efc_basis_ah=efc_basis_ah,
                            options=a.get("options"), params=a.get("params")))
    return out


def run_ablations(data_dir, out_dir, names=None, max_workers=None, report_every_s=300,
                  rate_basis_ah=RATE_BASIS_AH, efc_basis_ah=EFC_BASIS_AH):
    """Run the ablations (all by default): results in <out_dir>/<ablation>/."""
    _run_jobs(jobs(data_dir, out_dir, names, rate_basis_ah, efc_basis_ah), max_workers, report_every_s)


# ---------------------------------------------------------------------------
# Comparison with the full model
# ---------------------------------------------------------------------------
ABLATION_STYLE = {"no_isolation": ("tab:red", "--"), "no_strain_cracking": ("tab:purple", "-."),
                  "strain_cracking_no_temperature": ("tab:green", ":"), "no_li_trapping": ("tab:cyan", "--")}


def _available(out_dir, name, key):
    return os.path.exists(os.path.join(out_dir, name, f"{key}_results.npz"))


def compare(data_dir, out_dir, baseline_dir):
    """One row per (ablation, condition) that has results, next to the full model:
    mean |SoH gap| to the experiment, final RPT SoH, and final Si LAM / LLI."""
    import numpy as np
    import pandas as pd

    from . import analysis as an
    from .conditions import CONDITIONS

    def scores(res, curves, cond):
        s = an.summary(res, curves, cond, data_dir)
        m = an.life_mask(res)
        return dict(soh_gap_pp=s["mean_abs_soh_gap_pp"], final_rpt_soh=float(an.rpt_soh(res)[-1]),
                    final_efc=float(an.efc(res, res["rpt_thr"])[-1]), final_lam_si=float(res["LAM_si"][m][-1]),
                    final_lli=float(an.lli_pct(res, True)[m][-1]),
                    knee_efc=an.knee_summary(res, cond, data_dir)["model_knee_efc"])

    rows = []
    for name, a in ABLATIONS.items():
        for key in a["keys"]:
            if not _available(out_dir, name, key):
                continue
            cond = CONDITIONS[key]
            full = scores(*an.load(key, baseline_dir), cond)
            abl = scores(*an.load(key, os.path.join(out_dir, name)), cond)
            rows.append(dict(ablation=name, condition=key, **{f"{k}_full": v for k, v in full.items()},
                             **{f"{k}_ablation": v for k, v in abl.items()}))
    return pd.DataFrame(rows)


ABLATION_LABEL = {"no_isolation": "No isolation (mechanism I off)",
                  "no_strain_cracking": "No strain cracking (route B off)",
                  "strain_cracking_no_temperature": "Route B without temperature dependence",
                  "no_li_trapping": r"No lithium trapping ($\theta$ = 0)"}


def figure(data_dir, out_dir, baseline_dir, keys=None):
    """SoH (C/20 RPTs, top row) and Si LAM (bottom row) against EFC, one column
    per condition: measurement, full model, and each ablation that has results.
    Printed size, lettered panels, one shared legend (see style.py)."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    from . import analysis as an
    from . import data, style
    from .conditions import CONDITIONS

    style.use()
    keys = list(keys or CONDITIONS)
    fig, axes = plt.subplots(2, len(keys), figsize=(style.DOUBLE_IN, 3.9), squeeze=False, layout="constrained")
    shown = {}
    for c, key in enumerate(keys):
        cond = CONDITIONS[key]
        ax_s, ax_l = axes[0, c], axes[1, c]
        cap = data.load(cond, data_dir, "capacity_fade")
        c20, sub = cap[cap["source"] == "C/20 RPT"], cap[cap["source"] != "C/20 RPT"]
        ax_s.plot(c20["efc"], c20["capacity_retention_pct"], "o", ms=4, color="black", zorder=6)
        if len(sub):
            ax_s.plot(sub["efc"], sub["capacity_retention_pct"], "o", ms=4, mfc="white", mec="black", mew=1.0,
                      zorder=6)
            shown["sub"] = True
        lam = data.load(cond, data_dir, "LAM_derived")
        ax_l.plot(lam["efc"], lam["LAM_NE_silicon_pct"], "s", ms=3.5, color="black", zorder=6)
        runs = [("full", baseline_dir, ("tab:blue", "-"))]
        runs += [(name, os.path.join(out_dir, name), ABLATION_STYLE.get(name, ("gray", "--")))
                 for name, a in ABLATIONS.items() if key in a["keys"] and _available(out_dir, name, key)]
        for name, d, (col, ls) in runs:
            res, _ = an.load(key, d)
            ax_s.plot(an.efc(res, res["rpt_thr"]), an.rpt_soh(res), ls=ls, marker="o", ms=2.2, color=col)
            m = an.life_mask(res)
            ax_l.plot(an.efc(res, res["thr"])[m], res["LAM_si"][m], ls=ls, color=col)
            shown[name] = (col, ls)
        style.panel_label(ax_s, c)
        ax_s.text(0.04, 0.04, cond_tag(cond).replace(", ", "\n", 1), transform=ax_s.transAxes, ha="left",
                  va="bottom", fontsize=7.5)
        style.panel_label(ax_l, len(keys) + c)
        ax_l.set_xlabel("EFC since RPT1 [-]")
    axes[0, 0].set_ylabel("SoH, C/20 RPT [%]")
    axes[1, 0].set_ylabel("Si LAM [%]")
    handles = [Line2D([], [], color="black", ls="none", marker="o", ms=4, label="Measured")]
    if shown.pop("sub", False):
        handles.append(Line2D([], [], color="black", ls="none", marker="o", ms=4, mfc="white",
                              label="Measured, high-rate RPT"))
    col, ls = shown.pop("full")
    handles.append(Line2D([], [], color=col, ls=ls, label="Full model"))
    handles += [Line2D([], [], color=c_, ls=l_, label=ABLATION_LABEL.get(n, n)) for n, (c_, l_) in shown.items()]
    fig.legend(handles=handles, loc="outside upper center", ncol=4, frameon=False)
    return fig


def cond_tag(cond):
    """Short panel label for a condition, e.g. 'CELL064, 25 °C, full'."""
    window = cond.window.replace(" SoC", "").replace("-", "–")
    return f"{cond.cell}, {cond.temperature_K - 273.15:.0f} °C, {window}"
