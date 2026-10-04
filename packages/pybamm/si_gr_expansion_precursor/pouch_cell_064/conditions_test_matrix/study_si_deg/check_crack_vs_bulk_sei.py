# -*- coding: utf-8 -*-
"""study_si_deg item 1: quantify the split between Si's "normal" (bulk) SEI
growth and its "SEI on cracks" growth in the CURRENT high_temp_45C fork
(latest recipe: solvent-diffusion-limited SEI for Si, mult=5/positive E_sei,
stress-LAM bump, SI_CRACK_MULT=1.0 genuine cracking), to inform how much of
pre-knee porosity closure should be attributed to cracking vs bulk SEI when
designing the new 25 degC "cracking-dominated" model (item 2 of this study).

Imports the high_temp_45C fork UNMODIFIED (no monkey-patching of its
mechanism) -- this is a pure diagnostic on the EXISTING fork, not a new
mechanism test. Reads Si's bulk SEI thickness, SEI-on-cracks thickness, and
crack roughness ratio directly from the solution at each RPT checkpoint, and
reports both contributions to porosity closure using the SAME formula
reaction_driven_porosity.py itself uses (bulk: L_sei; cracks: L_sei_cr *
(roughness - 1)).
"""
import os
import sys

sys.path.insert(
    0,
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "high_temp_45C",
    ),
)
os.environ.setdefault("C64_EFC_LIMIT", "450")
os.environ.setdefault("C64_MAX_CYCLES", "550")
os.environ.setdefault("C64_NO_PLOT", "1")
os.environ.setdefault("C64_OUT_TAG", "crack_vs_bulk_diag")
# 2026-09-28: this fork now has its own output_variables memory
# optimization (ported from the shared script after a plain run got
# OOM-killed) -- this flag adds the 3 extra variables this diagnostic
# needs to that restricted list, applied correctly to BOTH the main and
# retry solvers (unlike a manual solver override in this script alone,
# which would miss the retry path).
os.environ.setdefault("C64_DIAG_CRACK_SPLIT", "1")

import numpy as np  # noqa: E402
import cell064_degradation_fit_45C as m  # noqa: E402

print("Config check (high_temp_45C latest recipe, unmodified):", flush=True)
print(f"  SEI={m.MODEL_OPTIONS_BASE['SEI']}, SOLVENT_MULT={m.SOLVENT_MULT}, "
      f"SI_ESEI={m.SI_ESEI}, SI_CRIT_STRESS={m.SI_CRIT_STRESS:.2e}, "
      f"SI_LAM_EAC={m.SI_LAM_EAC}, SI_CRACK_MULT env={os.environ.get('C64_SI_CRACK_MULT')}, "
      f"WIDTH={m.WIDTH_BASELINE}, F0={m.F0_BASELINE}", flush=True)

sol = m.run_degradation()
results = m.extract_results(sol)

rpt_thr = results["rpt_thr"]
rpt_efc = m.efc_from_throughput(rpt_thr) if rpt_thr.size else np.array([])

print("\n--- Bulk vs cracks SEI contribution to porosity closure, at each RPT ---", flush=True)
print(f"{'EFC':>8} {'L_bulk [nm]':>12} {'L_cracks_raw [nm]':>18} {'roughness':>10} "
      f"{'cracks_contrib [nm]':>20} {'bulk_frac':>10} {'cracks_frac':>12}", flush=True)

try:
    L_bulk_full = sol["X-averaged negative secondary SEI thickness [m]"].entries
    L_cracks_full = sol["X-averaged negative secondary SEI on cracks thickness [m]"].entries
    roughness_full = sol["X-averaged negative secondary electrode roughness ratio"].entries
    Qt_full = sol["Throughput capacity [A.h]"].entries

    for efc, thr in zip(rpt_efc, rpt_thr):
        idx = int(np.argmin(np.abs(Qt_full - thr)))
        L_bulk = float(L_bulk_full[idx])
        L_cracks_raw = float(L_cracks_full[idx])
        roughness = float(roughness_full[idx])
        cracks_contrib = L_cracks_raw * (roughness - 1.0)
        total = L_bulk + cracks_contrib
        bulk_frac = L_bulk / total if total > 0 else float("nan")
        cracks_frac = cracks_contrib / total if total > 0 else float("nan")
        print(f"{efc:8.1f} {L_bulk * 1e9:12.4f} {L_cracks_raw * 1e9:18.4f} {roughness:10.4f} "
              f"{cracks_contrib * 1e9:20.4f} {bulk_frac:10.3f} {cracks_frac:12.3f}", flush=True)

    # Final summary using the last available point.
    L_bulk_end = float(L_bulk_full[-1])
    L_cracks_raw_end = float(L_cracks_full[-1])
    roughness_end = float(roughness_full[-1])
    cracks_contrib_end = L_cracks_raw_end * (roughness_end - 1.0)
    total_end = L_bulk_end + cracks_contrib_end
    print(f"\nEnd of run: bulk={100 * L_bulk_end / total_end:.1f}%, "
          f"cracks={100 * cracks_contrib_end / total_end:.1f}% of total porosity-closing "
          f"SEI thickness (roughness={roughness_end:.3f})", flush=True)
except Exception as exc:  # noqa: BLE001
    print(f"Diagnostic readout failed: {type(exc).__name__}: {exc}", flush=True)
