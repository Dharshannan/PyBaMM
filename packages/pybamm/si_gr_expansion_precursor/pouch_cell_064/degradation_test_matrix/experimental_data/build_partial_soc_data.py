# -*- coding: utf-8 -*-
"""Build the partial-SoC-window cycling sets, CELL009 (15-95% SoC) and CELL026
(20-80% SoC), both 25 degC / 15 psi (103 kPa), in the same formats as
CELL017_45C/ (see build_cell017_45C_data.py). Source: the IC_PYBAMM_FIG1 data
package (cells/CELLxxx, parameters/all_cells_parameter_summary_combined.csv,
manifest/cells.csv).

Outputs, one folder per cell next to this file (CELL009_25C_15-95/, CELL026_25C_20-80/):
  CELLxxx_capacity_fade.csv        rpt, efc, discharge_capacity_Ah, source, capacity_retention_pct
  CELLxxx_LLI.csv                  rpt, efc, charge_LLI, joule_LLI, charge_LLI_loss, joule_LLI_loss
  CELLxxx_LAM_derived.csv          rpt, efc, discharge_Cn_Si/Gr/Cp, LAM_NE_graphite/silicon_pct, LAM_PE_pct
  CELLxxx_k_expansion_scale.csv    rpt, efc, k, source
  CELLxxx_reversible_expansion.csv cyc_seq, efc, reversible_expansion_um (per partial-window ageing cycle)
  CELLxxx_cycling_protocol.json    the ageing-cycle control, detected from the lifetime data
  low_rate_c20/CELLxxx_RPTnnn_lowrate_c20_discharge.csv  (copied as-is)

Conventions (the same as CELL017 unless noted):
  * RPT1 is the 100% / zero reference (RPT0 is an earlier pre-characterisation
    that the model loaders drop).
  * Some RPTs have no eSoH fit: CELL009 RPT2 and RPT4; CELL026 RPT2, RPT3 and
    RPT5. They were marked invalid only for missing raw expansion; their C/20
    discharge ran normally. Their capacity (the QC table's v_q_span_Ah) is kept
    in capacity_fade with source "C/20 RPT", but there is no LLI/LAM/k row and
    no discharge-curve file for them.
  * EFC is the package's, throughput / (2 x 2.5 Ah), the same basis as the
    CELL064 and CELL017 files.
  * Ageing control, detected here and written to *_cycling_protocol.json:
      CELL009: C/3 CC-CV charge to 4.150 V (CV to ~14 mA), C/3 CC discharge to 3.140 V
               (voltage-limited at both ends; Ah per cycle falls with ageing).
      CELL026: C/3 CC-CV charge to 3.969 V (CV to ~14 mA), C/3 discharge of a fixed
               1.500 Ah (coulomb counting, 60% of 2.5 Ah) with a 2.5 V safety floor.
"""
import argparse
import json
import os
import shutil

import numpy as np
import pandas as pd

from extract_lifetime_reversible_expansion import per_cycle_amplitude

HERE = os.path.dirname(os.path.abspath(__file__))
CELLS = {9: "CELL009_25C_15-95", 26: "CELL026_25C_20-80"}


def detect_protocol(pkl_path):
    """Summarise the partial-window ageing steps (|I| ~ C/3) from the lifetime
    data: CV voltage and end current of the charges, end voltage and Ah of the
    discharges, early vs late in life."""
    df = pd.read_pickle(pkl_path)
    I, V, q = df["current_A"].to_numpy(), df["voltage_V"].to_numpy(), df["capacity_Ah_cumulative"].to_numpy()
    efc = df["EFC"].to_numpy()
    state = np.where(I > 0.01, 1, np.where(I < -0.01, -1, 0))
    edge = np.flatnonzero(np.diff(state) != 0) + 1
    rows = []
    for s, e in zip(np.r_[0, edge], np.r_[edge, len(df)]):
        if e - s < 3 or state[s] == 0:
            continue
        med = np.median(np.abs(I[s:e]))
        dq = abs(q[e - 1] - q[s])
        if 0.7 < med < 1.0 and dq > 0.3:
            rows.append(dict(kind="chg" if state[s] > 0 else "dis", efc=efc[s], I=med, dQ=dq,
                             V_max=V[s:e].max(), V_min=V[s:e].min(), I_end=abs(I[e - 1])))
    st = pd.DataFrame(rows)
    ch, di = st[st.kind == "chg"], st[st.kind == "dis"]
    early, late = di[di.efc < 100], di[di.efc > 250]
    return {
        "ageing_current_A": round(float(st["I"].median()), 4),
        "charge_cv_voltage_V": round(float(ch["V_max"].median()), 4),
        "charge_cv_end_current_A": round(float(ch["I_end"].median()), 4),
        "discharge_end_voltage_V_early_median": round(float(early["V_min"].median()), 4),
        "discharge_end_voltage_V_late_median": round(float(late["V_min"].median()), 4),
        "discharge_Ah_early_median": round(float(early["dQ"].median()), 4),
        "discharge_Ah_late_median": round(float(late["dQ"].median()), 4),
        "n_ageing_discharges": int(len(di)),
    }


def build(pkg, cell_num):
    cell = f"CELL{cell_num:03d}"
    out_dir = os.path.join(HERE, CELLS[cell_num])
    os.makedirs(os.path.join(out_dir, "low_rate_c20"), exist_ok=True)
    u = pd.read_csv(os.path.join(pkg, "parameters", "all_cells_parameter_summary_combined.csv"))
    u = u[u["cell"] == cell_num].sort_values("EFC").reset_index(drop=True)
    ref = u[u["rpt"] == 1].iloc[0]
    qc = pd.read_csv(os.path.join(pkg, "cells", cell, "qc", f"{cell}_valid_rpt_and_cleaning_table.csv"))

    cap = pd.DataFrame({"rpt": qc["rpt"], "efc": qc["efc"], "discharge_capacity_Ah": qc["v_q_span_Ah"],
                        "source": "C/20 RPT"})
    cap["capacity_retention_pct"] = 100 * cap["discharge_capacity_Ah"] / ref["measured_capacity_Ah"]
    cap["esoh_fit"] = cap["rpt"].isin(u["rpt"])
    cap.to_csv(os.path.join(out_dir, f"{cell}_capacity_fade.csv"), index=False)

    lli = pd.DataFrame({"rpt": u["rpt"], "efc": u["EFC"], "charge_LLI": u["LLI_Ah"], "joule_LLI": u["LLI_Ah"]})
    lli["charge_LLI_loss"] = ref["LLI_Ah"] - lli["charge_LLI"]
    lli["joule_LLI_loss"] = lli["charge_LLI_loss"]
    lli.to_csv(os.path.join(out_dir, f"{cell}_LLI.csv"), index=False)

    lam = pd.DataFrame({"rpt": u["rpt"], "efc": u["EFC"], "discharge_Cn_Si": u["Cn_Si_Ah"],
                        "discharge_Cn_Gr": u["Cn_Gr_Ah"], "discharge_Cp": u["Cp_Ah"]})
    lam["LAM_NE_graphite_pct"] = 100 * (1 - lam["discharge_Cn_Gr"] / ref["Cn_Gr_Ah"])
    lam["LAM_NE_silicon_pct"] = 100 * (1 - lam["discharge_Cn_Si"] / ref["Cn_Si_Ah"])
    lam["LAM_PE_pct"] = 100 * (1 - lam["discharge_Cp"] / ref["Cp_Ah"])
    lam.to_csv(os.path.join(out_dir, f"{cell}_LAM_derived.csv"), index=False)

    k = pd.DataFrame({"rpt": u["rpt"], "efc": u["EFC"], "k": u["expansion_scale_k"], "source": "C/20 RPT"})
    k.to_csv(os.path.join(out_dir, f"{cell}_k_expansion_scale.csv"), index=False)

    pkl = os.path.join(pkg, "cells", cell, "lifetime", f"{cell}_lifetime_current_voltage_expansion_cleaned.pkl.gz")
    rev = per_cycle_amplitude(pkl)
    rev[["cyc_seq", "efc", "reversible_expansion_um"]].to_csv(
        os.path.join(out_dir, f"{cell}_reversible_expansion.csv"), index=False)

    man = pd.read_csv(os.path.join(pkg, "manifest", "cells.csv"))
    man = man[man["cell_label"] == cell].iloc[0]
    proto = detect_protocol(pkl)
    proto.update(cell=cell, soc_window_pct=[float(man["soc_window_low_pct"]), float(man["soc_window_high_pct"])],
                 temperature_C=float(man["cycling_temperature_C"]), pressure_kPa=float(man["stack_pressure_kPa"]),
                 capacity_knee_efc=float(man["capacity_knee_efc"]),
                 expansion_transition_efc=float(man["expansion_transition_efc"]),
                 rpt1_efc=float(ref["EFC"]))
    with open(os.path.join(out_dir, f"{cell}_cycling_protocol.json"), "w") as fh:
        json.dump(proto, fh, indent=2)

    src = os.path.join(pkg, "cells", cell, "low_rate_c20")
    for f in sorted(os.listdir(src)):
        shutil.copy2(os.path.join(src, f), os.path.join(out_dir, "low_rate_c20", f))

    print(f"===== {cell} -> {out_dir}")
    print(cap.round(3).to_string(index=False))
    print(lli.round(3).to_string(index=False))
    print(lam.round(2).to_string(index=False))
    print(k.round(3).to_string(index=False))
    print(f"reversible expansion: {len(rev)} cycles, EFC {rev['efc'].min():.1f}-{rev['efc'].max():.1f}")
    print(json.dumps(proto, indent=2))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pkg", required=True, help="IC_PYBAMM_FIG1 package root")
    for c in CELLS:
        build(ap.parse_args().pkg, c)
