"""Write the measured LAM and LLI tables the paper code reads (sigr.data) from an
expansion-free composite DMA.

    python tables.py [--source collaborators|new] [--data ../data]

source "collaborators" (the paper's choice): the collaborators' discharge fits
(voltage-based, no expansion in the objective), kept in the data bundle as
<data>/dma_collaborators_discharge.csv. source "new": our PyProBE refit,
results/dma_new.csv (the cross-check).

<folder>/<CELL>_LAM_derived.csv and <folder>/<CELL>_LLI.csv keep their columns;
LAM and LLI are relative to RPT1, and the EFC of each RPT is taken from the
cell's capacity_fade.csv. Needs pandas only (also called by the bundle builder).
"""
import argparse
import os

import pandas as pd

FOLDERS = {"CELL064": "CELL064_25C_full", "CELL017": "CELL017_45C", "CELL009": "CELL009_25C_15-95",
           "CELL026": "CELL026_25C_20-80"}
HERE = os.path.dirname(os.path.abspath(__file__))
COLLAB_FILE = "dma_collaborators_discharge.csv"


def collaborators_discharge(src_csv):
    """The discharge rows of the collaborators' combined charge/discharge eSoH table, in
    our column names, without any expansion column (for the data bundle)."""
    d = pd.read_csv(src_csv)
    d = d[d["direction"] == "discharge"]
    d = d[[c for c in d.columns if "expansion" not in c]]
    return d.rename(columns={"Cn_Ah": "Cn_total_Ah", "x_100_or_start": "x_start", "y_100_or_start": "y_start",
                             "lithium_inventory_Ah": "LLI_Ah", "si_ocp_scale_Sv": "si_ocp_scale",
                             "si_ocp_shift_Uoff_V": "si_ocp_shift_V"}).reset_index(drop=True)


def load_fits(source, data_dir):
    if source == "collaborators":
        return pd.read_csv(os.path.join(data_dir, COLLAB_FILE))
    return pd.read_csv(os.path.join(HERE, "results", "dma_new.csv"))


def write_tables(fits, data_dir, label):
    for cell, folder in FOLDERS.items():
        d = os.path.join(data_dir, folder)
        efc = pd.read_csv(os.path.join(d, f"{cell}_capacity_fade.csv")).set_index("rpt")["efc"]
        f = fits[fits["cell"] == cell].sort_values("rpt").reset_index(drop=True)
        ref = f[f["rpt"] == 1].iloc[0]
        lam = pd.DataFrame({"rpt": f["rpt"], "efc": f["rpt"].map(efc), "discharge_Cn_Si": f["Cn_Si_Ah"],
                            "discharge_Cn_Gr": f["Cn_Gr_Ah"], "discharge_Cp": f["Cp_Ah"]})
        lam["LAM_NE_graphite_pct"] = 100 * (1 - f["Cn_Gr_Ah"] / ref["Cn_Gr_Ah"])
        lam["LAM_NE_silicon_pct"] = 100 * (1 - f["Cn_Si_Ah"] / ref["Cn_Si_Ah"])
        lam["LAM_PE_pct"] = 100 * (1 - f["Cp_Ah"] / ref["Cp_Ah"])
        lam.to_csv(os.path.join(d, f"{cell}_LAM_derived.csv"), index=False)
        lli = pd.DataFrame({"rpt": f["rpt"], "efc": f["rpt"].map(efc), "charge_LLI": f["LLI_Ah"],
                            "joule_LLI": f["LLI_Ah"]})
        lli["charge_LLI_loss"] = ref["LLI_Ah"] - lli["charge_LLI"]
        lli["joule_LLI_loss"] = lli["charge_LLI_loss"]
        lli.to_csv(os.path.join(d, f"{cell}_LLI.csv"), index=False)
        print(f"{cell}: LAM/LLI tables from the {label} DMA, RPTs {list(f['rpt'])}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", choices=("collaborators", "new"), default="collaborators")
    ap.add_argument("--data", default=os.path.join(HERE, "..", "data"))
    a = ap.parse_args()
    write_tables(load_fits(a.source, a.data), a.data, a.source)


if __name__ == "__main__":
    main()
