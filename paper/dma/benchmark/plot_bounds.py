import os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from matplotlib.lines import Line2D
B = r"C:\Shannan_PhD_Stuff_Local\PyBaMM_Pressure\paper\dma"
sys.path.insert(0, B)
import composite_dma as C, run_dma as R
st = R._style(); DATA = os.path.join(B, "..", "data"); OUT = os.path.join(B, "results", "benchmark")
t = pd.read_csv(os.path.join(OUT, "bounds_fits.csv")); rep = pd.read_csv(os.path.join(OUT, "bounds_report.csv")).set_index("variant")
hc = C.HalfCells(DATA); q, v, d, _, _ = C.load_discharge(C.rpt_files(DATA, "CELL017")[7])
sty = {"ours": ("tab:blue", "-", "Ours (shift $\leq$ 30 mV)"), "wide_deform": ("tab:orange", "-", "Ours, shift $\leq$ 80 mV, scale $\geq$ 0.5"),
       "wide_deform+Cp": ("tab:green", "-", "... and $C_p$ 2.6–2.8 A h"), "collaborators": ("tab:red", ":", "Collaborators' fit")}
fig, ax = plt.subplots(1, 2, figsize=(st.DOUBLE_IN, 2.9), layout="constrained")
ax[0].plot(q, v, "k", lw=1.6); ax[1].plot(q, d, "k", lw=1.6)
for k, (c, ls, lab) in sty.items():
    r = t[(t.variant == k) & (t.rpt == 7)].iloc[0]; vm = hc.voltage([r[p] for p in C.PARAMS], q)
    ax[0].plot(q, vm, color=c, ls=ls, lw=1.0); ax[1].plot(q, np.gradient(vm, q), color=c, ls=ls, lw=1.0)
ax[1].set_ylim(np.percentile(d, 2) - 0.3, 0.05)
ax[0].set_ylabel("Voltage [V]"); ax[1].set_ylabel("dV/dQ [V/A h]")
for i, a in enumerate(ax):
    a.set_xlabel("Discharge capacity [A h]"); st.panel_label(a, i, "CELL017 RPT7")
fig.legend(handles=[Line2D([], [], color="k", lw=1.6, label="Measured")] + [
    Line2D([], [], color=c, ls=ls, label=f"{lab}: Gr LAM {rep.loc[k,'LAM_Gr']:.0f} %, Si LAM {rep.loc[k,'LAM_Si']:.0f} %")
    for k, (c, ls, lab) in sty.items()], loc="outside upper center", ncol=2, frameon=False)
print(st.save(fig, OUT, "bounds_rpt7")[0])
