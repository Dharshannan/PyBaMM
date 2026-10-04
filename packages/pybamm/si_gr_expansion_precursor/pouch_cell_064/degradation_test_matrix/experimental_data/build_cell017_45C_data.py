# -*- coding: utf-8 -*-
"""Build the CELL017 (45 degC, 15 psi / 103 kPa) experimental set in the same
formats as the CELL064 25 degC files, from the IC_PYBAMM_FIG1 data package
(cells/CELL017 + parameters/all_cells_parameter_summary_combined.csv).

Outputs (CELL017_45C/ next to this file), same columns as the CELL064 files:
  CELL017_capacity_fade.csv       rpt, efc, discharge_capacity_Ah, source, capacity_retention_pct
  CELL017_LLI.csv                 rpt, efc, charge_LLI, joule_LLI, charge_LLI_loss, joule_LLI_loss
  CELL017_LAM_derived.csv         rpt, efc, discharge_Cn_Si/Gr/Cp, LAM_NE_graphite/silicon_pct, LAM_PE_pct
  CELL017_k_expansion_scale.csv   rpt, efc, k, source
  CELL017_reversible_expansion.csv cyc_seq, efc, reversible_expansion_um
  low_rate_c20/CELL017_RPTnnn_lowrate_c20_discharge.csv  (copied as-is)

Conventions vs CELL064 (differences are forced by what the package contains):
  * RPT0 (EFC 6.6) is a 25 degC pre-characterisation (test name "P25C"), so
    RPT1 (first 45 degC RPT, EFC 8.72) is the 100% / zero reference, like
    CELL064's RPT1. RPT0 is kept in the files but the model loaders drop it.
  * Capacity retention is relative to RPT1's C/20 capacity (same rule as the
    CELL064 builder); the user's earlier quick overlay used capacity / 2.5 Ah
    instead -- the two differ by ~0.2%.
  * LLI: CELL064's charge_LLI came from a charge-curve fit not shipped in the
    package. Here both LLI columns are the package's common `LLI_Ah`
    (x_start*Cn_total + y_start*Cp from the discharge fit), which tracks
    CELL064's charge_LLI within ~1-5%.
  * Reversible expansion from extract_lifetime_reversible_expansion.py
    (validated on CELL064: median -2.9% vs its own pipeline, corr 0.92).
  * No high-rate-neighbor substitutes exist for CELL017 (its summary is
    empty) -- every RPT has a valid C/20 fit.
"""
import argparse
import os
import shutil

import pandas as pd

from extract_lifetime_reversible_expansion import per_cycle_amplitude

HERE = os.path.dirname(os.path.abspath(__file__))
CELL = "CELL017"


def build(pkg):
    out_dir = os.path.join(HERE, "CELL017_45C")
    os.makedirs(os.path.join(out_dir, "low_rate_c20"), exist_ok=True)
    u = pd.read_csv(os.path.join(pkg, "parameters", "all_cells_parameter_summary_combined.csv"))
    u = u[u["cell"] == 17].sort_values("EFC").reset_index(drop=True)
    ref = u[u["rpt"] == 1].iloc[0]

    cap = pd.DataFrame({"rpt": u["rpt"], "efc": u["EFC"],
                        "discharge_capacity_Ah": u["measured_capacity_Ah"],
                        "source": "C/20 RPT"})
    cap["capacity_retention_pct"] = 100 * cap["discharge_capacity_Ah"] / ref["measured_capacity_Ah"]
    cap.to_csv(os.path.join(out_dir, f"{CELL}_capacity_fade.csv"), index=False)

    lli = pd.DataFrame({"rpt": u["rpt"], "efc": u["EFC"],
                        "charge_LLI": u["LLI_Ah"], "joule_LLI": u["LLI_Ah"]})
    lli["charge_LLI_loss"] = ref["LLI_Ah"] - lli["charge_LLI"]
    lli["joule_LLI_loss"] = lli["charge_LLI_loss"]
    lli.to_csv(os.path.join(out_dir, f"{CELL}_LLI.csv"), index=False)

    lam = pd.DataFrame({"rpt": u["rpt"], "efc": u["EFC"], "discharge_Cn_Si": u["Cn_Si_Ah"],
                        "discharge_Cn_Gr": u["Cn_Gr_Ah"], "discharge_Cp": u["Cp_Ah"]})
    lam["LAM_NE_graphite_pct"] = 100 * (1 - lam["discharge_Cn_Gr"] / ref["Cn_Gr_Ah"])
    lam["LAM_NE_silicon_pct"] = 100 * (1 - lam["discharge_Cn_Si"] / ref["Cn_Si_Ah"])
    lam["LAM_PE_pct"] = 100 * (1 - lam["discharge_Cp"] / ref["Cp_Ah"])
    lam.to_csv(os.path.join(out_dir, f"{CELL}_LAM_derived.csv"), index=False)

    k = pd.DataFrame({"rpt": u["rpt"], "efc": u["EFC"], "k": u["expansion_scale_k"],
                      "source": "C/20 RPT"})
    k.to_csv(os.path.join(out_dir, f"{CELL}_k_expansion_scale.csv"), index=False)

    pkl = os.path.join(pkg, "cells", CELL, "lifetime",
                       f"{CELL}_lifetime_current_voltage_expansion_cleaned.pkl.gz")
    rev = per_cycle_amplitude(pkl)
    rev[["cyc_seq", "efc", "reversible_expansion_um"]].to_csv(
        os.path.join(out_dir, f"{CELL}_reversible_expansion.csv"), index=False)

    src = os.path.join(pkg, "cells", CELL, "low_rate_c20")
    for f in sorted(os.listdir(src)):
        shutil.copy2(os.path.join(src, f), os.path.join(out_dir, "low_rate_c20", f))

    print(cap.round(3).to_string(index=False))
    print(lli.round(3).to_string(index=False))
    print(lam.round(2).to_string(index=False))
    print(k.round(3).to_string(index=False))
    print(f"reversible expansion: {len(rev)} cycles, EFC {rev['efc'].min():.1f}-{rev['efc'].max():.1f}")
    print("->", out_dir)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pkg", required=True, help="IC_PYBAMM_FIG1 package root")
    build(ap.parse_args().pkg)
