# Si OCP ageing deformation: C/20 discharge and charge (Supplementary Figs. S20-S31, Table S6)

Three variants of the adopted shared set, everything else unchanged:

| variant | delithiation branch (Eq. 6) | lithiation branch | runs (local) |
|---|---|---|---|
| main text | S_end 0.8663, dU_end -9.08 mV | undeformed | refit_local/runs/CHG_ADOPT2 (= ADOPT2) |
| both | as main text | S 0.6898, dU +22.19 mV (option "OCP aging deformation branches": "both") | refit_local/runs/CHG_OCPL_RPT4 |
| undeformed | S_end 1, dU_end 0 | undeformed | refit_local/runs/CHG_NODEF |

The lithiation end values are the collaborators' charge-curve fit of CELL064 at RPT4 relative to RPT1 (method of
Wan2026). Candidates: refit_local/cmaes/CHG_*.json; all three keep the RPT charge (<key>_rpt_charge_curves.csv).
The paper keeps the main-text variant.

- Measured charges: `python tools/import_c20_charges.py` (FIG1 package -> data/<folder>/low_rate_c20/, the
  expansion columns dropped).
- Figures and scores (from paper/):
  `python -W ignore charge_branch/figures.py refit_local/runs/CHG_ADOPT2 refit_local/runs/CHG_OCPL_RPT4 refit_local/runs/CHG_NODEF`
  -> figures/si_ocp_* (PDFs copied to manuscript_overleaf/figures/paper_pics/si/), scores.csv, scores_summary.csv.

Result (2026-11-05), means over the four cells (main / both / undeformed):
- SoH gap 2.82 / 2.83 / 2.81 pp; per cell within 0.1 pp.
- Discharge V RMSE 75.8 / 76.4 / 81.6 mV (CELL064 85 / 85 / 97), discharge dV/dSoC 0.200 / 0.200 / 0.216.
- Charge V RMSE (all RPTs) 32.1 / 30.3 / 32.2 mV, dV/dSoC 0.230 / 0.221 / 0.235; the lithiation gain is late in
  life (CELL009 RPT7 36 -> 22, RPT8 53 -> 37 mV).
- Largest mode change at the last RPT: 1.7 pp (CELL017 LLI, both) and 2.0 pp (CELL026 Si LAM, undeformed).
So the deformation shapes the C/20 curves (each branch the curve it sets) but not the degradation.

Pairing: full-window runs keep no formation charge (its own cycle), so the EFC-0 charge comes from
simulate.bol_curves (cached as <key>_bol_charge.npz in the run folder); RPTs are paired within 20 EFC.
