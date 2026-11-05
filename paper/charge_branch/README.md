# C/20 charge: Si lithiation OCP fixed vs deformed (Supplementary Figs. S20-S25, Table S6)

The paper runs deform only the Si delithiation OCP (option "OCP aging deformation branches": "delithiation").
This check also deforms the lithiation branch ("both") with S_end,lith 0.6898 and dU_end,lith +22.19 mV (the
collaborators' charge-curve fit of CELL064 at RPT4 relative to RPT1, method of Wan2026), on the same LAM ramp;
everything else is the adopted shared set. The paper keeps the lithiation branch undeformed.

- Measured charges: `python tools/import_c20_charges.py` (FIG1 package -> data/<folder>/low_rate_c20/, the
  expansion columns dropped).
- Runs (local): refit_local/runs/CHG_ADOPT2 (fixed, = ADOPT2) and CHG_OCPL_RPT4 (deformed), from
  refit_local/cmaes/CHG_*.json; both keep the RPT charge (<key>_rpt_charge_curves.csv).
- Figures and scores: `python -W ignore charge_branch/figures.py refit_local/runs/CHG_ADOPT2 refit_local/runs/CHG_OCPL_RPT4`
  -> figures/ (PDFs copied to manuscript_overleaf/figures/paper_pics/si/), scores.csv, scores_summary.csv.

Result (2026-11-05): mean charge V RMSE 32.1 -> 30.3 mV, dV/dSoC 0.230 -> 0.221, almost all late in life
(CELL009 RPT7 36 -> 22, RPT8 53 -> 37 mV); discharge and SoH unchanged (SoH gaps within 0.1 pp, largest mode
change at the last RPT 1.7 pp, CELL017 LLI).

Pairing: full-window runs keep no formation charge (its own cycle), so the EFC-0 charge comes from
simulate.bol_curves (cached as <key>_bol_charge.npz in the run folder); RPTs are paired within 20 EFC.
