# -*- coding: utf-8 -*-
"""Lithium budget at every model C/20 RPT (diagnostic only, no physics change).

Usage: python lithium_budget.py study|fork TAG
Prints, relative to BoL and as % of nominal capacity: change in lithium in the
negative and positive particles and in the electrolyte, and the counted sinks
(negative SEI, negative SEI on cracks, LAM-trapped lithium per electrode).
Closure check: particles + electrolyte + sinks should sum to ~0.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
which, TAG = sys.argv[1], sys.argv[2]
V = {
    "Li neg particles": "Total lithium in negative electrode [mol]",
    "Li pos particles": "Total lithium in positive electrode [mol]",
    "Li electrolyte": "Total lithium in electrolyte [mol]",
    "SEI": "Loss of lithium to negative SEI [mol]",
    "SEI on cracks": "Loss of lithium to negative SEI on cracks [mol]",
    "LAM-trapped Gr": "Loss of lithium due to loss of primary active material in negative electrode [mol]",
    "LAM-trapped Si": "Loss of lithium due to loss of secondary active material in negative electrode [mol]",
    "LAM-trapped pos": "Loss of lithium due to loss of active material in positive electrode [mol]",
}
os.environ["C64_DIAG_EXTRA_VARS"] = "|".join(V.values())
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

sol = m.run_degradation()
F = 96485.33212
Q = m.NOMINAL_CAP_AH
pct = lambda mol: 100 * mol * F / 3600 / Q
d = {k: sol[v].entries for k, v in V.items()}
efc = m.efc_from_throughput(sol["Throughput capacity [A.h]"].entries)
print(f"\n=== Lithium budget ({which}, {TAG}); % of nominal, change since BoL ===", flush=True)
hdr = ["EFC"] + list(V) + ["closure"]
print("  ".join(f"{h:>15s}" for h in hdr), flush=True)
for e in (0, 50, 100, 130, 150, 171.4, 197.4):
    if e > efc[-1]:
        continue
    i = int(np.argmin(np.abs(efc - e)))
    row = [pct(d[k][i] - d[k][0]) if k.startswith("Li ") else pct(d[k][i]) for k in V]
    closure = row[0] + row[1] + row[2] + sum(row[3:])  # stocks change + counted sinks
    print("  ".join(f"{x:15.2f}" for x in [efc[i]] + row + [closure]), flush=True)
