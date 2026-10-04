# -*- coding: utf-8 -*-
"""Overpotential breakdown at every model C/20 RPT discharge, to explain why
late RPT curves sit lower in one recipe than another.

Usage: python overpotential_breakdown.py study|fork TAG
  study = cell064_degradation_fit_crack_dominated.py (25degC_baseline1 defaults)
  fork  = high_temp_45C/cell064_degradation_fit_45C.py (+ whatever C64_* env)
Prints, at Q = 0.05 Ah and 0.6 Ah into each RPT discharge: terminal V, OCV,
and each loss term (V), plus remaining active-material fractions. Also saves
a stacked bar chart of the loss terms at Q = 0.6 Ah vs RPT EFC.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
which, TAG = sys.argv[1], sys.argv[2]
LOSS = {
    "neg reaction": "X-averaged battery negative reaction overpotential [V]",
    "pos reaction": "X-averaged battery positive reaction overpotential [V]",
    "neg particle": "Battery negative particle concentration overpotential [V]",
    "pos particle": "Battery positive particle concentration overpotential [V]",
    "electrolyte ohmic": "X-averaged battery electrolyte ohmic losses [V]",
    "electrolyte conc.": "X-averaged battery concentration overpotential [V]",
    "solid ohmic": "X-averaged battery solid phase ohmic losses [V]",
    "SEI film": "X-averaged SEI film overpotential [V]",
}
EXTRA = list(LOSS.values()) + [
    "Battery open-circuit voltage [V]",
    "X-averaged negative electrode porosity",
    "X-averaged negative electrode primary active material volume fraction",
    "X-averaged negative electrode secondary active material volume fraction",
]
os.environ["C64_DIAG_EXTRA_VARS"] = "|".join(EXTRA)
os.environ.setdefault("C64_NO_PLOT", "1")
os.environ["C64_OUT_TAG"] = TAG
if which == "fork":
    sys.path.insert(0, os.path.join(os.path.dirname(HERE), "high_temp_45C"))
    os.environ.setdefault("C64_T_AMBIENT_K", "298.15")
    import cell064_degradation_fit_45C as m  # noqa: E402
else:
    sys.path.insert(0, HERE)
    import cell064_degradation_fit_crack_dominated as m  # noqa: E402

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sol = m.run_degradation()
rows = []
for cyc in sol.cycles:
    st = cyc.steps[0]
    I = st["Current [A]"].entries
    if not (0.05 < np.median(np.abs(I)) < 0.3):
        continue  # not a C/20 RPT discharge
    q = st["Discharge capacity [A.h]"].entries
    q = q - q[0]
    efc = float(m.efc_from_throughput(st["Throughput capacity [A.h]"].entries[0]))
    row = {"efc": efc, "Q_end": float(q[-1])}
    for qq in (0.05, 0.6):
        if q[-1] < qq:
            continue
        i = int(np.searchsorted(q, qq))
        row[f"V@{qq}"] = float(st["Terminal voltage [V]"].entries[i])
        row[f"OCV@{qq}"] = float(st["Battery open-circuit voltage [V]"].entries[i])
        for k, v in LOSS.items():
            row[f"{k}@{qq}"] = float(np.ravel(st[v].entries)[i] if np.ndim(st[v].entries) == 1
                                     else st[v].entries[..., i].mean())
    for lab, v in (("eps_Gr", EXTRA[-2]), ("eps_Si", EXTRA[-1]), ("porosity", EXTRA[-3])):
        row[lab] = float(np.ravel(st[v].entries)[0])
    rows.append(row)

print(f"\n=== Overpotential breakdown ({which}, {TAG}) ===", flush=True)
for r in rows:
    print(f"RPT @ EFC {r['efc']:.0f}: Q_end={r['Q_end']:.3f} Ah, eps_Gr={r['eps_Gr']:.4f}, "
          f"eps_Si={r['eps_Si']:.4f}, porosity={r['porosity']:.4f}", flush=True)
    for qq in (0.05, 0.6):
        if f"V@{qq}" not in r:
            continue
        losses = "  ".join(f"{k}={1000*r[f'{k}@{qq}']:+.1f}" for k in LOSS)
        print(f"   Q={qq}: V={r[f'V@{qq}']:.4f} OCV={r[f'OCV@{qq}']:.4f} "
              f"(OCV-V={1000*(r[f'OCV@{qq}']-r[f'V@{qq}']):.1f} mV) | {losses} [mV]", flush=True)

ok = [r for r in rows if "V@0.6" in r]
fig, ax = plt.subplots(figsize=(9, 5))
bottom = np.zeros(len(ok))
x = np.arange(len(ok))
for k in LOSS:
    vals = np.array([abs(1000 * r[f"{k}@0.6"]) for r in ok])
    ax.bar(x, vals, bottom=bottom, label=k)
    bottom += vals
ax.set_xticks(x, [f"{r['efc']:.0f}" for r in ok])
ax.set_xlabel("RPT EFC [-]")
ax.set_ylabel("|loss| at Q = 0.6 Ah into C/20 discharge [mV]")
ax.set_title(f"Overpotential breakdown ({which}: {TAG})")
ax.legend(fontsize=8)
fig.tight_layout()
out = os.path.join(HERE, f"overpotential_breakdown_{TAG}.png")
fig.savefig(out, dpi=130)
print("Saved:", out, flush=True)
