"""Supplementary figure: the silicon delithiation and lithiation OCPs as they deform with the silicon LAM
(Eq. 6 of the main text; the lithiation branch as in the 'both deformed' variant of the Supplementary section).
From paper/:

    python charge_branch/ocp_deform_figure.py

U -> U + f_LAM [(S_end - 1) U + dU_end], f_LAM = 1 - eps_Si / eps_Si,0 (base_ocp._apply_ocp_aging_deformation).
Writes charge_branch/figures/si_ocp_deformation.{png,pdf}: one panel, delithiation solid and lithiation dashed,
coloured by the Si LAM from 0 to the largest reached by the adopted runs (CELL064, results/); the undeformed
lithiation branch of the main-text model dotted.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

sys.path.insert(0, os.getcwd())
import sigr  # noqa: E402
from sigr import analysis as an, parameters, style  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
OCP_DIR = os.path.join("data", "ocp")
RESULTS = "results"
DELITH = "silicon_ocp_bol_deformed_cell064.csv"
LITH = "silicon_ocp_bol_lithiation_deformed_cell064_cell076transplant.csv"
# lithiation end of the 'both deformed' variant (refit_local/cmaes/OCPL_RPT4.json): CELL064 charge, RPT4 vs RPT1
LITH_SCALE, LITH_SHIFT = 0.6898, 0.02219
N_LINES = 10


def base(name):
    d = pd.read_csv(os.path.join(OCP_DIR, name))
    keep = (d["sto"] >= 0.0) & (d["sto"] <= 1.0)
    return d["sto"].to_numpy()[keep], d["p"].to_numpy()[keep]


def deform(u, f, scale, shift):
    return u + f * ((scale - 1.0) * u + shift)


def max_lam_si():
    """Largest Si LAM [%] reached by any of the four adopted runs."""
    out = 0.0
    for key in sigr.CONDITIONS:
        if os.path.exists(os.path.join(RESULTS, f"{key}_results.npz")):
            res, _ = an.load(key, RESULTS)
            out = max(out, float(res["LAM_si"][an.life_mask(res)].max()))
    return out


def main():
    style.use()
    d_scale, d_shift = parameters.V7["si_ocp_deform_scale"], parameters.V7["si_ocp_deform_shift"]
    x_d, u_d = base(DELITH)
    x_l, u_l = base(LITH)
    lam_max = max_lam_si()
    cmap, norm = plt.cm.viridis, plt.Normalize(0.0, lam_max)

    fig, ax = plt.subplots(figsize=(style.SINGLE_IN, 2.6), layout="constrained")
    for lam in np.linspace(0.0, lam_max, N_LINES):
        c, f = cmap(norm(lam)), lam / 100.0
        ax.plot(x_d, deform(u_d, f, d_scale, d_shift), color=c, lw=0.8)
        ax.plot(x_l, deform(u_l, f, LITH_SCALE, LITH_SHIFT), color=c, lw=0.8, ls="--")
    ax.plot(x_l, u_l, color="k", lw=1.0, ls=":")
    ax.legend(handles=[Line2D([], [], color="0.3", label="Delithiation"),
                       Line2D([], [], color="0.3", ls="--", label="Lithiation, deformed"),
                       Line2D([], [], color="k", ls=":", label="Lithiation, undeformed (main text)")],
              loc="upper right")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel(r"Si stoichiometry [-]")
    ax.set_ylabel("Si OCP [V vs Li/Li$^+$]")
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label="Si LAM [%]", pad=0.02, aspect=25)
    paths = style.save(fig, OUT, "si_ocp_deformation")
    print(f"Si LAM up to {lam_max:.1f} %")
    print("\n".join(paths))


if __name__ == "__main__":
    main()
