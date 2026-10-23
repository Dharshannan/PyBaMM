"""Write the measured LAM and LLI tables the paper code reads (sigr.data) from the
new composite DMA (results/dma_new.csv: voltage + dV/dQ fits, no expansion).

    python tables.py [--fits results/dma_new.csv] [--data ../data]

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


def write_tables(fits_csv, data_dir):
    fits = pd.read_csv(fits_csv)
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
        print(f"{cell}: LAM/LLI tables from the new DMA, RPTs {list(f['rpt'])}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fits", default=os.path.join(HERE, "results", "dma_new.csv"))
    ap.add_argument("--data", default=os.path.join(HERE, "..", "data"))
    a = ap.parse_args()
    write_tables(a.fits, a.data)


if __name__ == "__main__":
    main()
