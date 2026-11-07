"""Supplementary figure: the silicon delithiation and lithiation OCPs as they deform with the silicon LAM fraction
(Eq. 6 of the main text; the lithiation branch as in the 'both deformed' variant of the Supplementary section).
From paper/:

    python charge_branch/ocp_deform_figure.py

U -> U + f_LAM [(S_end - 1) U + dU_end], f_LAM = 1 - eps_Si / eps_Si,0 (base_ocp._apply_ocp_aging_deformation).
Writes charge_branch/figures/si_ocp_deformation.{png,pdf}:
  a) delithiation branch (main text), b) lithiation branch (both deformed), coloured by f_LAM;
  c) the two branches at f_LAM = 0 and 1 (the hysteresis gap).
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
from sigr import parameters, style  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
OCP_DIR = os.path.join("data", "ocp")
DELITH = "silicon_ocp_bol_deformed_cell064.csv"
LITH = "silicon_ocp_bol_lithiation_deformed_cell064_cell076transplant.csv"
# lithiation end of the 'both deformed' variant (refit_local/cmaes/OCPL_RPT4.json): CELL064 charge, RPT4 vs RPT1
LITH_SCALE, LITH_SHIFT = 0.6898, 0.02219
F_LAM = np.linspace(0.0, 1.0, 6)


def base(name):
    d = pd.read_csv(os.path.join(OCP_DIR, name))
    keep = (d["sto"] >= 0.0) & (d["sto"] <= 1.0)
    return d["sto"].to_numpy()[keep], d["p"].to_numpy()[keep]


def deform(u, f, scale, shift):
    return u + f * ((scale - 1.0) * u + shift)


def main():
    style.use()
    d_scale, d_shift = parameters.V7["si_ocp_deform_scale"], parameters.V7["si_ocp_deform_shift"]
    x_d, u_d = base(DELITH)
    x_l, u_l = base(LITH)
    cmap, norm = plt.cm.viridis, plt.Normalize(0.0, 1.0)

    fig, ax = plt.subplots(1, 3, figsize=(style.DOUBLE_IN, 2.5), layout="constrained", sharey=True)
    for a, (x, u, s, dv, lab) in zip(ax[:2], ((x_d, u_d, d_scale, d_shift, "delithiation"),
                                               (x_l, u_l, LITH_SCALE, LITH_SHIFT, "lithiation"))):
        for f in F_LAM:
            a.plot(x, deform(u, f, s, dv), color=cmap(norm(f)), lw=1.0)
        a.set_xlabel(r"Si stoichiometry $x_\mathrm{Si}$ [-]")
    ax[0].set_ylabel("Si OCP [V vs Li/Li$^+$]")
    for f, ls in ((0.0, "-"), (1.0, "--")):
        ax[2].plot(x_d, deform(u_d, f, d_scale, d_shift), color="tab:blue", ls=ls, lw=1.0)
        ax[2].plot(x_l, deform(u_l, f, LITH_SCALE, LITH_SHIFT), color="tab:red", ls=ls, lw=1.0)
    ax[2].legend(handles=[Line2D([], [], color="tab:blue", label="Delithiation"),
                          Line2D([], [], color="tab:red", label="Lithiation"),
                          Line2D([], [], color="0.3", label=r"$f_\mathrm{LAM}=0$"),
                          Line2D([], [], color="0.3", ls="--", label=r"$f_\mathrm{LAM}=1$")], loc="upper right")
    ax[2].set_xlabel(r"Si stoichiometry $x_\mathrm{Si}$ [-]")
    for a in ax:
        a.set_xlim(0.0, 1.0)
    ax[0].set_ylim(0.0, None)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(sm, ax=ax[:2], label=r"Si LAM fraction $f_\mathrm{LAM}$ [-]", pad=0.02, aspect=25)
    for i, extra in enumerate(("delithiation", "lithiation", "both branches")):
        style.panel_label(ax[i], i, extra)
    paths = style.save(fig, OUT, "si_ocp_deformation")
    print("\n".join(paths))


if __name__ == "__main__":
    main()
