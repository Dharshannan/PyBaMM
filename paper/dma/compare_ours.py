"""Our DMA for the paper (pyprobe-dma, trimmed discharge, Si OCP deformed; results/dma_pyprobe_dma.csv)
against the collaborators' expansion-free discharge DMA, the given DMA (with expansion) and our PyProBE 2.6.1
refit (results/dma_new.csv), for the DMA comparison note.

    python compare_ours.py [--data ../data] [--out results]

Writes results/dma_ours_compare.* (modes: ours, collaborators, given), results/dma_ours_diff.csv (mean / max
absolute difference from ours, RPT1 onwards) and prints the last-RPT values.
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import composite_dma as C  # noqa: E402
import run_dma as R  # noqa: E402
import tables  # noqa: E402
from compare_collab import COLLAB_D  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OURS = dict(marker="o", ls="-", color="tab:orange", ms=3.8, lw=1.1)
DMAS = {"ours": "Our DMA, pyprobe-dma (voltage + dV/dQ)",
        "collaborators": "Collaborators' DMA, discharge (voltage only, no expansion)",
        "given": "Given DMA (voltage + dV/dQ + expansion)",
        "pyprobe": "Our PyProBE 2.6.1 refit (voltage + dV/dQ)"}


def all_modes(data_dir):
    """LAM and LLI of every DMA, relative to RPT1, on the EFC axis of the given DMA."""
    base = R.base_tables(data_dir)
    fits = {"ours": tables.load_fits("ours", data_dir), "collaborators": tables.load_fits("collaborators", data_dir),
            "given": tables.load_fits("collaborators-expansion", data_dir), "pyprobe": tables.load_fits("pyprobe", data_dir)}
    out = []
    for cell, b in base.items():
        efc0 = float(b.loc[b["rpt"] == 1, "EFC"].iloc[0])
        for name, f in fits.items():
            m = C.degradation_modes(f[f["cell"] == cell])
            m["efc"] = m["rpt"].map(dict(zip(b["rpt"], b["EFC"] - efc0)))
            m["cell"], m["dma"] = cell, name
            out.append(m)
    return pd.concat(out, ignore_index=True)


def plot(modes, out_dir, name="dma_ours_compare"):
    st = R._style()
    cells = list(C.FOLDERS)
    styles = (("given", R.GIVEN), ("collaborators", COLLAB_D), ("ours", OURS))
    fig = plt.figure(figsize=(st.DOUBLE_IN, 7.6), layout="constrained")
    rows = fig.subfigures(len(cells), 1)
    k = 0
    for i, cell in enumerate(cells):
        ax = rows[i].subplots(1, 4)
        rows[i].suptitle(R.TAG[cell], x=0.01, ha="left", fontsize=8, fontweight="bold")
        for j, (col, lab) in enumerate(R.MODES):
            a = ax[j]
            for dma_, sty in styles:
                m = modes[(modes["cell"] == cell) & (modes["dma"] == dma_)].sort_values("efc")
                a.plot(m["efc"], m[col], **sty)
            a.axhline(0, color="k", lw=0.5)
            a.set_ylabel(lab)
            if i == len(cells) - 1:
                a.set_xlabel("EFC since RPT1 [-]")
            st.panel_label(a, k)
            k += 1
    fig.legend(handles=[Line2D([], [], label=DMAS[d], **s) for d, s in styles[::-1]],
               loc="outside upper center", ncol=2, frameon=False)
    return st.save(fig, out_dir, name)[0]


def differences(modes):
    """Mean / max absolute difference of each DMA from ours, RPT1 onwards."""
    m = modes[modes["rpt"] >= 1].set_index(["cell", "rpt", "dma"])
    ours = m.xs("ours", level="dma")
    rows = []
    for name in ("collaborators", "given", "pyprobe"):
        d = (m.xs(name, level="dma") - ours).abs()
        rows.append({"against": name, **{f"{c}_{s}": getattr(d[c], s)() for c in ("LLI", "LAM_Si", "LAM_Gr", "LAM_PE")
                                          for s in ("mean", "max")}})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(HERE, "..", "data"))
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    a = ap.parse_args()
    modes = all_modes(a.data)
    print(plot(modes, a.out))
    d = differences(modes)
    d.to_csv(os.path.join(a.out, "dma_ours_diff.csv"), index=False)
    pd.set_option("display.width", 200)
    print(d.round(1).to_string(index=False))
    last = modes[modes["rpt"] == modes.groupby("cell")["rpt"].transform("max")]
    print(last.pivot(index="cell", columns="dma", values=["LAM_Si", "LAM_Gr", "LAM_PE", "LLI"]).round(1).to_string())


if __name__ == "__main__":
    main()
