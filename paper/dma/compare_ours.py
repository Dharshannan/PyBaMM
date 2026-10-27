"""The DMA comparison note's figures and numbers: our DMA (pyprobe-dma, trimmed discharge, Si OCP deformed;
results/dma_pyprobe_dma.csv, fitted curves in results/benchmark/new_def_curves.npz) against the given DMA
(with expansion) and the collaborators' expansion-free discharge DMA.

    python compare_ours.py [--data ../data] [--out results] [--collab-all <charge_discharge_esoh_combined_4cells.csv>]

Writes to results/: dma_modes_compare.* (given vs ours), dma_collab_vs_ours.*, dma_collab_vs_given.*,
dma_ours_compare.* (all three), si_ocp_deformation_params.* (with the collaborators' charge fits when
--collab-all is given), dma_fits_<cell>.* (voltage and dV/dQ of every RPT), dma_table.tex (per RPT, given /
ours) and dma_ours_diff.csv; prints the numbers quoted in the note.
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
import tables  # noqa: E402
from compare_collab import COLLAB_C, COLLAB_D, collab_table  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CURVES = os.path.join(HERE, "results", "benchmark", "new_def_curves.npz")
OURS = dict(marker="o", ls="-", color="tab:orange", ms=3.8, lw=1.1)
STYLE = {"given": R.GIVEN, "collaborators": COLLAB_D, "ours": OURS}
LABEL = {"ours": "Our DMA, pyprobe-dma (voltage + dV/dQ)",
         "collaborators": "Collaborators' DMA, discharge (voltage only, no expansion)",
         "given": "Given DMA (voltage + dV/dQ + expansion)"}
DEFORM_BOUNDS = ((0.75, 1.1), (-0.10, 0.03))  # pyprobe_dma_v3.DEFORM_BOUNDS, all RPTs but CELL017 RPT7


def fits(data_dir):
    return {"ours": tables.load_fits("ours", data_dir), "collaborators": tables.load_fits("collaborators", data_dir),
            "given": tables.load_fits("collaborators-expansion", data_dir)}


def all_modes(data_dir):
    """LAM and LLI of every DMA, relative to RPT1, on the EFC axis of the given DMA."""
    base = R.base_tables(data_dir)
    out = []
    for cell, b in base.items():
        efc0 = float(b.loc[b["rpt"] == 1, "EFC"].iloc[0])
        for name, f in fits(data_dir).items():
            m = C.degradation_modes(f[f["cell"] == cell])
            m["efc"] = m["rpt"].map(dict(zip(b["rpt"], b["EFC"] - efc0)))
            m["cell"], m["dma"] = cell, name
            out.append(m)
    return pd.concat(out, ignore_index=True)


def plot_modes(modes, dmas, out_dir, name):
    """LAM (Si, graphite, positive) and LLI against EFC since RPT1, one row per cell, for the DMAs named
    (RPT1 onwards: RPT0 precedes the reference)."""
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
            for d in dmas:
                m = modes[(modes["cell"] == cell) & (modes["dma"] == d) & (modes["rpt"] >= 1)].sort_values("efc")
                a.plot(m["efc"], m[col], **STYLE[d])
            a.axhline(0, color="k", lw=0.5)
            a.set_ylabel(lab)
            if i == len(cells) - 1:
                a.set_xlabel("EFC since RPT1 [-]")
            st.panel_label(a, k)
            k += 1
    fig.legend(handles=[Line2D([], [], label=LABEL[d], **STYLE[d]) for d in dmas[::-1]],
               loc="outside upper center", ncol=2, frameon=False)
    return st.save(fig, out_dir, name)[0]


def ours_voltage(curves, cell, rpt, q):
    """Our fitted terminal voltage at discharged capacity q (pyprobe-dma curve)."""
    qq, vv = curves[f"{cell}_{rpt}_q"], curves[f"{cell}_{rpt}_v"]
    ok = np.isfinite(qq) & np.isfinite(vv)
    o = np.argsort(qq[ok])
    return np.interp(q, qq[ok][o], vv[ok][o], left=np.nan, right=np.nan)


def model_curves(data_dir, f, curves, cell, rpt, q, hc):
    """(given, ours, collaborators) fitted voltage of one RPT at q."""
    def params(t):
        return [t.loc[(t["cell"] == cell) & (t["rpt"] == rpt), p].iloc[0] for p in C.PARAMS]
    return hc.voltage(params(f["given"]), q), ours_voltage(curves, cell, rpt, q), hc.voltage(params(f["collaborators"]), q)


def scores(data_dir, f, curves):
    """Voltage and dV/dQ RMSE of the three DMAs on the 300 points of the fit window."""
    hc = C.HalfCells(data_dir)
    rows = []
    for cell in C.FOLDERS:
        for rpt, path in C.rpt_files(data_dir, cell).items():
            q, v, dvdq, _, _ = C.load_discharge(path)
            row = dict(cell=cell, rpt=rpt)
            for name, vm in zip(("given", "ours", "collaborators"), model_curves(data_dir, f, curves, cell, rpt, q, hc)):
                row[f"rmse_V_mV_{name}"] = 1e3 * np.sqrt(np.nanmean((vm - v) ** 2))
                row[f"rmse_dVdQ_{name}"] = np.sqrt(np.nanmean((np.gradient(vm, q) - dvdq) ** 2))
            rows.append(row)
    return pd.DataFrame(rows)


def plot_fits(data_dir, cell, f, curves, out_dir):
    """Each C/20 RPT: voltage (upper panel) and dV/dQ (lower panel), measured and the three DMA fits; four RPTs
    per block of two rows."""
    st = R._style()
    hc = C.HalfCells(data_dir)
    files = C.rpt_files(data_dir, cell)
    rpts = sorted(files)
    ncol = 4
    nblk = int(np.ceil(len(rpts) / ncol))
    fig, ax = plt.subplots(2 * nblk, ncol, figsize=(st.DOUBLE_IN, 2.6 * nblk + 0.35), layout="constrained",
                           squeeze=False)
    styles = (dict(ls="--", color="0.45"), dict(color="tab:orange"), dict(ls=":", color="tab:red"))
    k = 0
    for n, r in enumerate(rpts):
        av, ad = ax[2 * (n // ncol), n % ncol], ax[2 * (n // ncol) + 1, n % ncol]
        q, v, dvdq, _, _ = C.load_discharge(files[r])
        av.plot(q, v, color="k", lw=1.6)
        ad.plot(q, dvdq, color="k", lw=1.6)
        vms = model_curves(data_dir, f, curves, cell, r, q, hc)
        for vm, sty in zip(vms, styles):
            av.plot(q, vm, lw=1.0, **sty)
            ad.plot(q, np.gradient(vm, q), lw=1.0, **sty)
        ad.set_ylim(np.percentile(dvdq, 2) - 0.25, max(0.05, dvdq.max() + 0.1))
        rmse = " / ".join(f"{1e3 * np.sqrt(np.nanmean((vm - v) ** 2)):.1f}" for vm in vms)
        av.text(0.03, 0.04, f"RMSE [mV]\n{rmse}", transform=av.transAxes, ha="left", va="bottom", fontsize=7)
        st.panel_label(av, k, f"RPT{r}")
        st.panel_label(ad, k + 1)
        k += 2
        ad.set_xlabel("Discharge capacity [A h]")
    for a in list(ax.flat):
        if not a.lines:
            a.set_visible(False)
    for b in range(nblk):
        ax[2 * b, 0].set_ylabel("Voltage [V]")
        ax[2 * b + 1, 0].set_ylabel("dV/dQ [V/A h]")
    handles = [Line2D([], [], color="k", lw=1.6, label="Measured, C/20 discharge"),
               Line2D([], [], label="Given DMA", **styles[0]), Line2D([], [], label="Our DMA (pyprobe-dma)", **styles[1]),
               Line2D([], [], label="Collaborators' DMA (no expansion)", **styles[2])]
    fig.legend(handles=handles, loc="outside upper center", ncol=4, frameon=False)
    return st.save(fig, out_dir, f"dma_fits_{cell}")[0]


def plot_deformation(data_dir, f, out_dir, collab_all=None):
    """Fitted Si OCP scale and shift against EFC: given, ours, collaborators' discharge (and charge) fits."""
    st = R._style()
    base = R.base_tables(data_dir)
    fig, ax = plt.subplots(2, 4, figsize=(st.DOUBLE_IN, 4.2), layout="constrained", sharex="col")
    for j, cell in enumerate(C.FOLDERS):
        b = base[cell]
        efc0 = float(b.loc[b["rpt"] == 1, "EFC"].iloc[0])
        efc_of = dict(zip(b["rpt"], b["EFC"] - efc0))
        series = []
        for name in ("given", "ours", "collaborators"):
            t = f[name][(f[name]["cell"] == cell) & (f[name]["rpt"] >= 1)].sort_values("rpt")
            series.append((t, t["rpt"].map(efc_of), STYLE[name]))
        if collab_all is not None:
            t = collab_all[(collab_all["cell"] == cell) & (collab_all["direction"] == "charge")
                           & (collab_all["rpt"] >= 1)].sort_values("rpt")
            series.append((t, t["EFC"] - efc0, COLLAB_C))
        for t, e, sty in series:
            ax[0, j].plot(e, t["si_ocp_scale"], **sty)
            ax[1, j].plot(e, 1e3 * t["si_ocp_shift_V"], **sty)
        ax[1, j].set_xlabel("EFC since RPT1 [-]")
        st.panel_label(ax[0, j], j, cell)
        st.panel_label(ax[1, j], 4 + j)
    ax[0, 0].set_ylabel("Si OCP scale $a$ [-]")
    ax[1, 0].set_ylabel("Si OCP shift $b$ [mV]")
    handles = [Line2D([], [], label="Given DMA, discharge", **R.GIVEN), Line2D([], [], label="Our DMA, discharge", **OURS),
               Line2D([], [], label="Collaborators, discharge", **COLLAB_D)]
    if collab_all is not None:
        handles.append(Line2D([], [], label="Collaborators, charge", **COLLAB_C))
    fig.legend(handles=handles, loc="outside upper center", ncol=len(handles), frameon=False)
    return st.save(fig, out_dir, "si_ocp_deformation_params")[0]


def latex_table(score, modes, path):
    """Per-RPT table for the note: voltage RMSE and the degradation modes, given / ours."""
    m = modes.pivot_table(index=["cell", "rpt"], columns="dma", values=["LAM_Si", "LAM_Gr", "LAM_PE", "LLI"])
    lines = []
    for cell in C.FOLDERS:
        sc = score[score["cell"] == cell].sort_values("rpt")
        for i, r in enumerate(sc.itertuples()):
            g = m.loc[(cell, r.rpt)]
            name = R.TAG[cell].replace("%", r"\%").replace("°C", r"$^\circ$C").replace("–", "--") if i == 0 else ""
            cells = [f"{g[(c, 'given')]:.0f} / {g[(c, 'ours')]:.0f}" for c in ("LAM_Si", "LAM_Gr", "LAM_PE", "LLI")]
            lines.append(f"{name} & {r.rpt} & {r.rmse_V_mV_given:.1f} / {r.rmse_V_mV_ours:.1f} & " + " & ".join(cells)
                         + r" \\")
        lines.append(r"\midrule")
    lines[-1] = r"\bottomrule"
    head = (r"\begin{tabular}{lrrrrrr}" + "\n" + r"\toprule" + "\n"
            + r"Cell & RPT & V RMSE [mV] & LAM Si [\%] & LAM Gr [\%] & LAM PE [\%] & LLI [\%] \\" + "\n"
            + r"\midrule" + "\n")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(head + "\n".join(lines) + "\n" + r"\end{tabular}" + "\n")
    return path


def differences(modes, ref):
    """Mean / max absolute difference of the other DMAs from ref, RPT1 onwards."""
    m = modes[modes["rpt"] >= 1].set_index(["cell", "rpt", "dma"])
    r = m.xs(ref, level="dma")
    rows = []
    for name in sorted(set(modes["dma"]) - {ref}):
        d = (m.xs(name, level="dma") - r).abs()
        rows.append({"ref": ref, "against": name, **{f"{c}_{s}": getattr(d[c], s)()
                                                     for c in ("LLI", "LAM_Si", "LAM_Gr", "LAM_PE") for s in ("mean", "max")}})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(HERE, "..", "data"))
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    ap.add_argument("--collab-all", default=None, help="the collaborators' combined charge/discharge CSV")
    a = ap.parse_args()
    f = fits(a.data)
    curves = dict(np.load(CURVES))
    modes = all_modes(a.data)
    for dmas, name in ((("given", "ours"), "dma_modes_compare"), (("ours", "collaborators"), "dma_collab_vs_ours"),
                       (("given", "collaborators"), "dma_collab_vs_given"),
                       (("given", "collaborators", "ours"), "dma_ours_compare")):
        print(plot_modes(modes, dmas, a.out, name))
    print(plot_deformation(a.data, f, a.out, collab_table(a.collab_all) if a.collab_all else None))
    for cell in C.FOLDERS:
        print(plot_fits(a.data, cell, f, curves, a.out))
    score = scores(a.data, f, curves)
    print(latex_table(score, modes[modes["dma"].isin(("given", "ours"))], os.path.join(a.out, "dma_table.tex")))
    d = pd.concat([differences(modes, "ours"), differences(modes, "collaborators")], ignore_index=True)
    d.to_csv(os.path.join(a.out, "dma_ours_diff.csv"), index=False)
    pd.set_option("display.width", 220)
    print(d.round(1).to_string(index=False))
    s1 = score[score["rpt"] >= 0]
    print("voltage RMSE [mV] min/median/max:", {n: [round(s1[f"rmse_V_mV_{n}"].agg(x), 1) for x in ("min", "median", "max")]
                                               for n in ("given", "ours", "collaborators")})
    print("dV/dQ RMSE ours < given in", int((s1["rmse_dVdQ_ours"] < s1["rmse_dVdQ_given"]).sum()), "of", len(s1))
    o, g = f["ours"].set_index(["cell", "rpt"]), f["given"].set_index(["cell", "rpt"])
    cn = (o["Cn_Si_Ah"] + o["Cn_Gr_Ah"]) / (g["Cn_Si_Ah"] + g["Cn_Gr_Ah"]) - 1
    print("total anode capacity, ours vs given: median %.1f %%, range %.1f to %.1f %%" % (100 * cn.median(), 100 * cn.min(), 100 * cn.max()))
    (a0, a1), (b0, b1) = DEFORM_BOUNDS
    at = o[(np.isclose(o["si_ocp_scale"], a0, atol=2e-3) | np.isclose(o["si_ocp_scale"], a1, atol=2e-3)
            | np.isclose(o["si_ocp_shift_V"], b0, atol=1e-3) | np.isclose(o["si_ocp_shift_V"], b1, atol=1e-3))]
    print("ours: deformation on a bound in", len(at), "fits:", list(at.index))
    last = modes[modes["rpt"] == modes.groupby("cell")["rpt"].transform("max")]
    print(last.pivot(index="cell", columns="dma", values=["LAM_Si", "LAM_Gr", "LAM_PE", "LLI"]).round(1).to_string())
    print(modes[(modes["rpt"] >= 1)].groupby("dma")[["LAM_Gr", "LAM_PE"]].agg(["min", "max"]).round(1))


if __name__ == "__main__":
    main()
