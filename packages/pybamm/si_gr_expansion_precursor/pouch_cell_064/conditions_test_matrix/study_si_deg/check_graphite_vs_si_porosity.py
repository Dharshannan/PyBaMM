# -*- coding: utf-8 -*-
"""study_si_deg follow-up (2026-09-28, iter13's finding): SI_CRACK_MULT had
ZERO measurable effect on this recipe's knee timing or any other metric --
surprising, since cracks were found to dominate (65.6%) porosity closure at
the 45 degC recipe. Directly checks whether GRAPHITE's own (unmodified, "ec
reaction limited") SEI growth is actually the dominant contributor to
electrode porosity closure in THIS 25 degC recipe, which would undercut the
whole "cracking-dominated Si SEI" premise this study is testing.

Compares graphite's own porosity-closing SEI thickness (Primary, no
cracks contribution needed since graphite's SI_CRACK_MULT-equivalent,
GR_CRACK_RATE_MULT, is untouched at its small default 0.1) against Si's
COMBINED (bulk + cracks) contribution, using iter08's adopted config
(now the script's own defaults).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("C64_EFC_LIMIT", "220")
os.environ.setdefault("C64_NO_PLOT", "1")
os.environ.setdefault("C64_OUT_TAG", "graphite_vs_si_diag")

import numpy as np  # noqa: E402
import pybamm  # noqa: E402
import cell064_degradation_fit_crack_dominated as m  # noqa: E402

_DIAG_VARS = list(m._CORE_OUTPUT_VARIABLES) + [
    "X-averaged negative primary SEI thickness [m]",
    "X-averaged negative secondary SEI thickness [m]",
    "X-averaged negative secondary SEI on cracks thickness [m]",
    "X-averaged negative secondary electrode roughness ratio",
    "X-averaged negative primary electrode roughness ratio",
]
m.solver = pybamm.IDAKLUSolver(
    root_tol=m.SOLVER_TOL, atol=m.SOLVER_TOL, rtol=m.SOLVER_TOL,
    output_variables=_DIAG_VARS,
)

print("Config check (iter08 adopted defaults):", flush=True)
print(f"  yield={m.SI_REDIRECT_LAM_YIELD}, exponent={m.SI_LAM_ISO_EXPONENT}, "
      f"SOLVENT_MULT={m.SOLVENT_MULT}, SI_CRACK_MULT env={os.environ.get('C64_SI_CRACK_MULT', '1.0 (default)')}",
      flush=True)

sol = m.run_degradation()
results = m.extract_results(sol)
rpt_thr = results["rpt_thr"]
rpt_efc = m.efc_from_throughput(rpt_thr) if rpt_thr.size else np.array([])

print("\n--- Graphite (Primary) vs Si (Secondary, bulk+cracks) porosity-closing SEI thickness ---", flush=True)
print(f"{'EFC':>8} {'L_gr_bulk [nm]':>15} {'gr_rough':>9} {'gr_contrib[nm]':>15} "
      f"{'L_si_bulk[nm]':>14} {'L_si_cr_raw[nm]':>16} {'si_rough':>9} "
      f"{'si_contrib[nm]':>15} {'gr_frac':>8} {'si_frac':>8}", flush=True)
try:
    L_gr_full = sol["X-averaged negative primary SEI thickness [m]"].entries
    gr_rough_full = sol["X-averaged negative primary electrode roughness ratio"].entries
    L_si_bulk_full = sol["X-averaged negative secondary SEI thickness [m]"].entries
    L_si_cr_full = sol["X-averaged negative secondary SEI on cracks thickness [m]"].entries
    si_rough_full = sol["X-averaged negative secondary electrode roughness ratio"].entries
    Qt_full = sol["Throughput capacity [A.h]"].entries

    for efc, thr in zip(rpt_efc, rpt_thr):
        idx = int(np.argmin(np.abs(Qt_full - thr)))
        L_gr = float(L_gr_full[idx])
        gr_rough = float(gr_rough_full[idx])
        gr_contrib = L_gr * gr_rough  # graphite's own crack contribution folded in if roughness>1
        L_si_bulk = float(L_si_bulk_full[idx])
        L_si_cr_raw = float(L_si_cr_full[idx])
        si_rough = float(si_rough_full[idx])
        si_cracks_contrib = L_si_cr_raw * (si_rough - 1.0)
        si_contrib = L_si_bulk + si_cracks_contrib
        total = gr_contrib + si_contrib
        gr_frac = gr_contrib / total if total > 0 else float("nan")
        si_frac = si_contrib / total if total > 0 else float("nan")
        print(f"{efc:8.1f} {L_gr * 1e9:15.4f} {gr_rough:9.4f} {gr_contrib * 1e9:15.4f} "
              f"{L_si_bulk * 1e9:14.4f} {L_si_cr_raw * 1e9:16.4f} {si_rough:9.4f} "
              f"{si_contrib * 1e9:15.4f} {gr_frac:8.3f} {si_frac:8.3f}", flush=True)

    L_gr_end = float(L_gr_full[-1])
    gr_rough_end = float(gr_rough_full[-1])
    gr_contrib_end = L_gr_end * gr_rough_end
    L_si_bulk_end = float(L_si_bulk_full[-1])
    L_si_cr_end = float(L_si_cr_full[-1])
    si_rough_end = float(si_rough_full[-1])
    si_contrib_end = L_si_bulk_end + L_si_cr_end * (si_rough_end - 1.0)
    total_end = gr_contrib_end + si_contrib_end
    print(f"\nEnd of run: graphite={100 * gr_contrib_end / total_end:.1f}%, "
          f"Si(bulk+cracks)={100 * si_contrib_end / total_end:.1f}% of total "
          f"porosity-closing SEI thickness", flush=True)
except Exception as exc:  # noqa: BLE001
    print(f"Diagnostic readout failed: {type(exc).__name__}: {exc}", flush=True)
