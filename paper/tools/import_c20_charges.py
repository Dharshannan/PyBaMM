"""Copy the measured C/20 RPT charges (constant-current part) of the four cells into data/ (from paper/):

    python tools/import_c20_charges.py [<FIG1 package cells folder>]

Source: the collaborators' FIG1 data package, cells/<CELL>/charge_c20/<CELL>_RPTnnn_charge_c20_cc.csv. Its RPT
numbering is ours; its EFC column (twice ours) and the expansion and dV/dQ columns are dropped. Writes
data/<folder>/low_rate_c20/<CELL>_RPTnnn_lowrate_c20_charge.csv (cell, rpt, test_name, time_h, Q_charge_Ah,
current_A, voltage_terminal_V), the format sigr.data.load_rpt_charge reads.
"""
import glob
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.getcwd())
import sigr  # noqa: E402

PKG = (r"C:\Users\ds3420\Downloads\IC_PYBAMM_FIG1_DATA_PACKAGE_v1_20260927"
       r"\IC_PYBAMM_FIG1_DATA_PACKAGE_v1_20260927\cells")
COLS = ["cell", "rpt", "test_name", "time_h", "Q_charge_Ah", "current_A", "voltage_terminal_V"]


def main(pkg=PKG, data_dir="data"):
    for cond in sigr.CONDITIONS.values():
        dest = os.path.join(data_dir, cond.data_folder, "low_rate_c20")
        n = 0
        for f in sorted(glob.glob(os.path.join(pkg, cond.cell, "charge_c20", f"{cond.cell}_RPT*_charge_c20_cc.csv"))):
            r = int(re.search(r"_RPT(\d{3})_", f).group(1))
            d = pd.read_csv(f).rename(columns={"voltage_measured_V": "voltage_terminal_V"})
            d["rpt"] = r
            d[COLS].to_csv(os.path.join(dest, f"{cond.cell}_RPT{r:03d}_lowrate_c20_charge.csv"), index=False)
            n += 1
        print(f"{cond.cell}: {n} C/20 charges -> {dest}")


if __name__ == "__main__":
    main(*sys.argv[1:2])
