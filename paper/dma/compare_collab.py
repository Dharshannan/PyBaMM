"""Compare the collaborators' expansion-free DMA (discharge and charge fits, one CSV)
with our new DMA and with the given DMA, and the Si OCP deformation on both branches.

    python compare_collab.py --collab <charge_discharge_esoh_combined_4cells.csv> --package <Fig. 1 data package>

Writes to results/: dma_collab_vs_new.*, dma_collab_vs_given.*, si_ocp_deformation_params.*,
si_ocp_branches.*, and prints a check of the charge fits (our forward model on the
C/20 charge data with their parameters and the lithiation-branch OCPs).
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import composite_dma as C  # noqa: E402
import run_dma as R  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
COLLAB_D = dict(marker="^", ls=":", color="tab:red", ms=4.2, lw=1.1, mfc="white", mew=1.0)
COLLAB_C = dict(marker="v", ls=":", color="tab:green", ms=4.2, lw=1.1)


def collab_table(path):
    """Their fits in our column names, without any expansion column."""
    d = pd.read_csv(path)
    d = d[[c for c in d.columns if "expansion" not in c]]
    d = d.rename(columns={"Cn_Ah": "Cn_total_Ah", "x_100_or_start": "x_start", "y_100_or_start": "y_start",
                          "lithium_inventory_Ah": "LLI_Ah", "si_ocp_scale_Sv": "si_ocp_scale",
                          "si_ocp_shift_Uoff_V": "si_ocp_shift_V"})
    return d


def collab_modes(collab, direction, base):
    out = []
    for cell, b in base.items():
        efc0 = float(b.loc[b["rpt"] == 1, "EFC"].iloc[0])
        t = collab[(collab["cell"] == cell) & (collab["direction"] == direction)]
        m = C.degradation_modes(t)
        m["efc"] = m["rpt"].map(dict(zip(t["rpt"], t["EFC"] - efc0)))
        m["cell"], m["dma"] = cell, f"collab_{direction}"
        out.append(m)
    return pd.concat(out, ignore_index=True)


def plot_pair(modes, other, other_style, other_label, out_dir, name):
    st = R._style()
    cells = list(C.FOLDERS)
    fig = plt.figure(figsize=(st.DOUBLE_IN, 7.6), layout="constrained")
    rows = fig.subfigures(len(cells), 1)
    k = 0
    for i, cell in enumerate(cells):
        ax = rows[i].subplots(1, 4)
        rows[i].suptitle(R.TAG[cell], x=0.01, ha="left", fontsize=8, fontweight="bold")
        for j, (col, lab) in enumerate(R.MODES):
            a = ax[j]
            for dma_, sty in ((other, other_style), ("collab_discharge", COLLAB_D)):
                m = modes[(modes["cell"] == cell) & (modes["dma"] == dma_)].sort_values("efc")
                a.plot(m["efc"], m[col], **sty)
            a.axhline(0, color="k", lw=0.5)
            a.set_ylabel(lab)
            if i == len(cells) - 1:
                a.set_xlabel("EFC since RPT1 [-]")
            st.panel_label(a, k)
            k += 1
    fig.legend(handles=[Line2D([], [], label=other_label, **other_style),
                        Line2D([], [], label="Collaborators' DMA, discharge (voltage only, no expansion)", **COLLAB_D)],
               loc="outside upper center", ncol=2, frameon=False)
    return st.save(fig, out_dir, name)[0]


def plot_deformation_params(new, base, collab, out_dir):
    """Fitted Si OCP scale and shift against EFC: discharge fits (given, new, collaborators')
    and the collaborators' charge fits."""
    st = R._style()
    cells = list(C.FOLDERS)
    fig, ax = plt.subplots(2, 4, figsize=(st.DOUBLE_IN, 4.2), layout="constrained", sharex="col")
    for j, cell in enumerate(cells):
        b = base[cell]
        efc0 = float(b.loc[b["rpt"] == 1, "EFC"].iloc[0])
        efc_of = dict(zip(b["rpt"], b["EFC"] - efc0))
        n = new[new["cell"] == cell].sort_values("rpt")
        series = [(b.sort_values("rpt"), b.sort_values("rpt")["EFC"] - efc0, R.GIVEN),
                  (n, n["rpt"].map(efc_of), R.NEW)]
        for direction, sty in (("discharge", COLLAB_D), ("charge", COLLAB_C)):
            t = collab[(collab["cell"] == cell) & (collab["direction"] == direction)].sort_values("rpt")
            series.append((t, t["EFC"] - efc0, sty))
        for t, e, sty in series:
            ax[0, j].plot(e, t["si_ocp_scale"], **sty)
            ax[1, j].plot(e, 1e3 * t["si_ocp_shift_V"], **sty)
        ax[1, j].set_xlabel("EFC since RPT1 [-]")
        st.panel_label(ax[0, j], j, cell)
        st.panel_label(ax[1, j], 4 + j)
    ax[0, 0].set_ylabel("Si OCP scale $a$ [-]")
    ax[1, 0].set_ylabel("Si OCP shift $b$ [mV]")
    fig.legend(handles=[Line2D([], [], label="Given DMA, discharge", **R.GIVEN),
                        Line2D([], [], label="New DMA, discharge", **R.NEW),
                        Line2D([], [], label="Collaborators, discharge", **COLLAB_D),
                        Line2D([], [], label="Collaborators, charge", **COLLAB_C)],
               loc="outside upper center", ncol=4, frameon=False)
    return st.save(fig, out_dir, "si_ocp_deformation_params")[0]


def _read_ocp(path):
    d = pd.read_excel(path) if path.endswith("xlsx") else pd.read_csv(path)
    d = d.apply(pd.to_numeric, errors="coerce").dropna().sort_values("sto")
    return d["sto"].to_numpy(), d["p"].to_numpy()


def plot_branches(package, data_dir, new, collab, out_dir, model_end=(0.8663, -0.00908), cell="CELL064"):
    """Si OCP, delithiation (discharge) and lithiation (charge) branches: the half-cell data,
    the model's beginning-of-life curves and its end-of-life aging deformation (delithiation
    only), and the fitted deformations of the first and last RPT of `cell`."""
    st = R._style()
    mt = os.path.join(package, "shared", "model_tables")
    xd, ud = _read_ocp(os.path.join(mt, "Silicon_OCP_Delithiation_verbrugge.xlsx"))
    xl, ul = _read_ocp(os.path.join(mt, "Silicon_OCP_Lithiation_updated.csv"))
    xmd, umd = _read_ocp(os.path.join(data_dir, "ocp", "silicon_ocp_bol_deformed_cell064.csv"))
    xml, uml = _read_ocp(os.path.join(data_dir, "ocp", "silicon_ocp_bol_lithiation_deformed_cell064_cell076transplant.csv"))
    x = np.linspace(0.02, 0.98, 300)
    fig, ax = plt.subplots(1, 2, figsize=(st.DOUBLE_IN, 2.9), layout="constrained", sharey=True)
    a = ax[0]
    a.plot(x, np.interp(x, xd, ud), color="k", lw=1.4, label="Half-cell data (Verbrugge)")
    ubol = np.interp(x, xmd, umd)
    a.plot(x, ubol, color="tab:purple", lw=1.2, label="Model, beginning of life")
    a.plot(x, model_end[0] * ubol + model_end[1], color="tab:purple", lw=1.2, ls="--",
           label="Model, end of life (Si LAM = 100 %)")
    n = new[new["cell"] == cell].sort_values("rpt")
    t = collab[(collab["cell"] == cell) & (collab["direction"] == "discharge")].sort_values("rpt")
    for df, sty, lab in ((n, dict(color="tab:blue"), "New DMA"), (t, dict(color="tab:red"), "Collaborators")):
        last = df.iloc[-1]
        a.plot(x, last["si_ocp_scale"] * np.interp(x, xd, ud) + last["si_ocp_shift_V"], lw=1.0, ls=":",
               label=f"{lab}, RPT{int(last['rpt'])} fit", **sty)
    a.set_title("Delithiation (discharge)", fontsize=7.5)
    a = ax[1]
    a.plot(x, np.interp(x, xl, ul), color="k", lw=1.4, label="Half-cell data")
    a.plot(x, np.interp(x, xml, uml), color="tab:purple", lw=1.2,
           label="Model, all of life (no aging deformation)")
    # The collaborators' charge-fit deformations are not drawn: their charge fits use lithiation
    # OCPs we do not have (check_charge_fits cannot reproduce them with the package's OCPs).
    a.set_title("Lithiation (charge)", fontsize=7.5)
    for i, a in enumerate(ax):
        a.set_xlabel("Si stoichiometry [-]")
        a.legend(loc="upper right")
        st.panel_label(a, i, cell)
    ax[0].set_ylabel("Si OCP [V]")
    return st.save(fig, out_dir, "si_ocp_branches")[0]


def check_charge_fits(package, collab, cells=("CELL064", "CELL009")):
    """Our forward model on the C/20 charge data with the collaborators' charge parameters and
    the lithiation-branch OCPs (NMC622 delithiating, graphite and Si lithiating): RMSE [mV]."""
    mt = os.path.join(package, "shared", "model_tables")
    xn, un = _read_ocp(os.path.join(mt, "NMC622_OCP_Delithiation_ds.csv"))
    xg, ug = C._monotone(*_read_ocp(os.path.join(mt, "Graphite_OCP_Lithiation.csv")), decreasing=True)
    xs, us = C._monotone(*_read_ocp(os.path.join(mt, "Silicon_OCP_Lithiation_updated.csv")), decreasing=True)
    from pyprobe.analysis import degradation_mode_analysis as dma

    out = []
    for cell in cells:
        for _, r in collab[(collab["cell"] == cell) & (collab["direction"] == "charge")].iterrows():
            f = os.path.join(package, "cells", cell, "charge_c20", f"{cell}_RPT{int(r['rpt']):03d}_charge_c20_cc.csv")
            if not os.path.exists(f):
                continue
            d = pd.read_csv(f, usecols=["Q_charge_Ah", "voltage_measured_V"]).dropna()
            q, v = d["Q_charge_Ah"].to_numpy(), d["voltage_measured_V"].to_numpy()
            qt = q.max()
            anode = dma.CompositeOCP.from_data(xs, r["si_ocp_scale"] * us + r["si_ocp_shift_V"], xg, ug)
            anode.comp_fraction = r["Cn_Si_Ah"] / (r["Cn_Si_Ah"] + r["Cn_Gr_Ah"])
            x = r["x_start"] - (qt - q) / (r["Cn_Si_Ah"] + r["Cn_Gr_Ah"])  # x_start: anode at 100 % SoC
            y = r["y_start"] + (qt - q) / r["Cp_Ah"]
            vm = np.interp(y, xn, un) - anode.eval(x)
            keep = (q > 0.02 * qt) & (q < 0.98 * qt)
            out.append((cell, int(r["rpt"]), 1e3 * np.sqrt(np.mean((vm - v)[keep] ** 2)), 1e3 * r["rmse_V_V"]))
    return pd.DataFrame(out, columns=["cell", "rpt", "rmse_mV_ours", "rmse_mV_reported"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--collab", required=True)
    ap.add_argument("--package", required=True)
    ap.add_argument("--data", default=os.path.join(HERE, "..", "data"))
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    a = ap.parse_args()
    collab = collab_table(a.collab)
    base = R.base_tables(a.data)
    new = pd.read_csv(os.path.join(a.out, "dma_new.csv"))
    modes = pd.concat([R.modes_table(a.data, new), collab_modes(collab, "discharge", base)], ignore_index=True)
    modes.to_csv(os.path.join(a.out, "dma_modes_with_collab.csv"), index=False)
    print(plot_pair(modes, "new", R.NEW, "New DMA, PyProBE (voltage + dV/dQ)", a.out, "dma_collab_vs_new"))
    print(plot_pair(modes, "base", R.GIVEN, "Given DMA (voltage + dV/dQ + expansion)", a.out, "dma_collab_vs_given"))
    print(plot_deformation_params(new, base, collab, a.out))
    print(plot_branches(a.package, a.data, new, collab, a.out))
    for cell in C.FOLDERS:  # the per-RPT fits, now with the collaborators' discharge fits as well
        print(R.plot_fits(a.data, cell, base[cell], new, a.out, collab[collab["direction"] == "discharge"])[0])
    pd.set_option("display.width", 200)
    print(check_charge_fits(a.package, collab).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
