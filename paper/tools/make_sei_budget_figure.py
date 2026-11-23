"""SEI volume budget of the paper runs (Supplementary figure, internal review 2026-11-21).

    python paper/tools/make_sei_budget_figure.py [--results paper/results] [--out paper/paper_pics]

Per cell: cumulative SEI volume in the negative electrode (lithium lost to SEI / z x partial molar volume, z = 2),
split into silicon crack faces, silicon particle surface and graphite, against EFC; horizontal lines at the initial
pore volume eps_0 L A and at the volume that brings the porosity to its floor (eps_0 - eps_min) L A; dotted line at
the model knee (two-segment least-squares fit to the model's C/3 capacity). Also prints the numbers quoted in the text.
"""
import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, HERE)
from sigr import analysis as an, simulate as S, style  # noqa: E402
from sigr.ablations import cond_tag  # noqa: E402
from sigr.conditions import CONDITIONS  # noqa: E402
from sigr.figures import C_GR, C_SEI, C_SI  # noqa: E402

F = 96485.33212
Z_SEI = 2
KEYS = ["25deg_full_soc", "45deg_full_soc", "25deg_15-95_soc", "25deg_20-80_soc"]


def two_segment_knee(x, y, n_min=3):
    """Breakpoint of a continuous two-segment least-squares fit."""
    best = (np.inf, np.nan)
    for xb in np.linspace(x[n_min], x[-n_min], 400):
        A = np.column_stack([np.ones_like(x), x, np.clip(x - xb, 0, None)])
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        r = float(np.sum((A @ c - y) ** 2))
        if r < best[0]:
            best = (r, xb)
    return best[1]


def budget(key, results, data_dir):
    cond = CONDITIONS[key]
    res, _ = an.load(key, results)
    pv = S.parameter_values(cond, data_dir, S.RATE_BASIS_AH, p=S.V7)
    vol = pv["Negative electrode thickness [m]"] * pv["Electrode height [m]"] * pv["Electrode width [m]"] * 1e6
    eps0, eps_min = pv["Negative electrode porosity"], pv["Negative electrode porosity floor"]
    cm3_per_ah = 3600 / F / Z_SEI * pv["Secondary: SEI partial molar volume [m3.mol-1]"] * 1e6
    m = an.life_mask(res)
    e = an.efc(res, res["thr"])[m]
    si, cr = res["Q_si_sei"][m], res["Q_si_sei_cracks"][m]
    gr = res["Q_side"][m] - si - cr
    knee = two_segment_knee(an.efc(res, res["age_thr"]), an.ageing_soh(res))
    pmin, e_age = res["age_min_porosity"], an.efc(res, res["age_thr"])
    floor_efc = next((float(e_age[i]) for i in range(pmin.size) if pmin[i] <= eps_min + 0.002), np.nan)
    return dict(cond=cond, efc=e, v_cr=cr * cm3_per_ah, v_si=si * cm3_per_ah, v_gr=gr * cm3_per_ah,
                pore=eps0 * vol, usable=(eps0 - eps_min) * vol, knee=knee, floor_efc=floor_efc,
                cm3_per_ah=cm3_per_ah)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    ap.add_argument("--data", default=os.path.join(HERE, "data"))
    ap.add_argument("--out", default=os.path.join(HERE, "paper_pics"))
    a = ap.parse_args()
    style.use()
    fig, axs = plt.subplots(2, 2, figsize=(style.DOUBLE_IN, 4.2), sharey=True, constrained_layout=True)
    for i, (key, ax) in enumerate(zip(KEYS, axs.ravel())):
        b = budget(key, a.results, a.data)
        e = b["efc"]
        ax.stackplot(e, b["v_cr"], b["v_si"], b["v_gr"], colors=[C_SEI, C_SI, C_GR], alpha=0.85, lw=0,
                     labels=["Si crack faces", "Si particle surface", "graphite"])
        ax.axhline(b["pore"], color="k", lw=0.9, label=r"initial pore volume $\varepsilon_0 L A$")
        ax.axhline(b["usable"], color="k", lw=0.9, ls="--",
                   label=r"pore volume down to the floor $(\varepsilon_0-\varepsilon_\mathrm{min}) L A$")
        ax.axvline(b["knee"], color="0.3", lw=0.9, ls=":", label="model knee")
        ax.set_xlim(0, e[-1])
        ax.set_xlabel("EFC")
        if i % 2 == 0:
            ax.set_ylabel(r"SEI volume (cm$^3$)")
        style.panel_label(ax, i, cond_tag(b["cond"]))
        tot = b["v_cr"][-1] + b["v_si"][-1] + b["v_gr"][-1]
        print(f"{b['cond'].cell}: end SEI {tot:.2f} cm3 = {tot / b['pore']:.2f} x pore volume "
              f"(cracks {b['v_cr'][-1]:.2f}, Si surface {b['v_si'][-1]:.2f}, Gr {b['v_gr'][-1]:.2f}); "
              f"pore {b['pore']:.2f}, to floor {b['usable']:.2f} cm3; {b['cm3_per_ah']:.2f} cm3/Ah; "
              f"model knee {b['knee']:.0f} EFC; porosity at the floor from {b['floor_efc']:.0f} EFC")
    axs[0, 1].legend(loc="upper left", fontsize=7)
    style.save(fig, a.out, "si_sei_volume_budget")


if __name__ == "__main__":
    main()
