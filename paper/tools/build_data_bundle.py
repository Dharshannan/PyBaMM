"""Assemble the data directory the paper code reads (the Zenodo data record)
from the study's local files. For the author; readers download the bundle.

    python paper/tools/build_data_bundle.py [--dest paper/data] [--pouch-data <Pouch_Data dir>]
                                            [--package <collaborators' Fig. 1 data package>]

From the collaborators' data package it also takes, for the composite DMA
(paper/dma): the three half-cell OCPs used on discharge (as CSV), and the
given (base) DMA fit table of each cell, without any expansion column.
"""
import argparse
import glob
import os
import shutil

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
STUDY = os.path.join(REPO, "packages", "pybamm", "si_gr_expansion_precursor", "pouch_cell_064")
EXP = os.path.join(STUDY, "degradation_test_matrix", "experimental_data")
OCP = os.path.join(STUDY, "parameters", "data")
POUCH_DATA = "C:/Shannan_PhD_Stuff_Local/Si_Gr_Expansion_Precursor/Pouch_Data"
TABLES = ("capacity_fade", "LLI", "LAM_derived")  # no expansion data in the paper release
PACKAGE = "C:/Users/ds3420/Downloads/IC_PYBAMM_FIG1_DATA_PACKAGE_v1_20260927/IC_PYBAMM_FIG1_DATA_PACKAGE_v1_20260927"
# Half-cell OCPs on discharge: NMC622 lithiating ("sto" in the file is 1 - y), graphite and Si delithiating.
HALF_CELL_OCPS = {"NMC622_OCP_Lithiation.xlsx": "nmc622_lithiation.csv",
                  "Graphite_OCP_Delithiation.xlsx": "graphite_delithiation.csv",
                  "Silicon_OCP_Delithiation_verbrugge.xlsx": "silicon_delithiation_verbrugge.csv"}
# Base DMA columns kept (the voltage-and-capacity part of the given fit; nothing on expansion).
BASE_DMA_COLUMNS = ["cell_label", "rpt", "EFC", "measured_capacity_Ah", "Cn_Si_Ah", "Cn_Gr_Ah", "Cn_total_Ah", "Cp_Ah",
                    "x_start", "y_start", "LLI_Ah", "si_ocp_scale", "si_ocp_shift_V", "rmse_voltage_V",
                    "rmse_dVdQ_V_per_Ah"]
FOLDERS = {"CELL064": "CELL064_25C_full", "CELL017": "CELL017_45C", "CELL009": "CELL009_25C_15-95",
           "CELL026": "CELL026_25C_20-80"}


def copy(src, dst_dir):
    os.makedirs(dst_dir, exist_ok=True)
    shutil.copy2(src, dst_dir)


def dma_inputs(package, dest):
    """Half-cell OCPs (CSV) and the base DMA table (low-rate C/20 fits only) of each cell."""
    import pandas as pd

    out = os.path.join(dest, "ocp_halfcell")
    os.makedirs(out, exist_ok=True)
    for src, name in HALF_CELL_OCPS.items():
        pd.read_excel(os.path.join(package, "shared", "model_tables", src))[["sto", "p"]].to_csv(
            os.path.join(out, name), index=False)
    summ = pd.read_csv(os.path.join(package, "parameters", "all_cells_parameter_summary_combined.csv"))
    summ = summ[summ["source_type"] == "low_rate_c20_rpt_fit"]
    for cell, folder in FOLDERS.items():
        summ.loc[summ["cell_label"] == cell, BASE_DMA_COLUMNS].sort_values("rpt").to_csv(
            os.path.join(dest, folder, f"{cell}_dma_base.csv"), index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dest", default=os.path.join(REPO, "paper", "data"))
    ap.add_argument("--pouch-data", default=POUCH_DATA, help="source of CELL064's C/20 RPT discharge curves")
    ap.add_argument("--package", default=PACKAGE, help="collaborators' data package (half-cell OCPs, base DMA)")
    a = ap.parse_args()

    for f in sorted(glob.glob(os.path.join(OCP, "*.csv"))):
        copy(f, os.path.join(a.dest, "ocp"))

    d064 = os.path.join(a.dest, "CELL064_25C_full")
    for t in TABLES:
        copy(os.path.join(EXP, f"CELL064_{t}.csv"), d064)
    for f in sorted(glob.glob(os.path.join(a.pouch_data, "cell064_data", "data_timeseries", "low_rate_c20",
                                           "CELL064_RPT*_lowrate_c20_discharge.csv"))):
        copy(f, os.path.join(d064, "low_rate_c20"))

    for folder in ("CELL017_45C", "CELL009_25C_15-95", "CELL026_25C_20-80"):
        cell = folder.split("_")[0]
        for t in TABLES:
            copy(os.path.join(EXP, folder, f"{cell}_{t}.csv"), os.path.join(a.dest, folder))
        for f in sorted(glob.glob(os.path.join(EXP, folder, "low_rate_c20", "*_discharge.csv"))):
            copy(f, os.path.join(a.dest, folder, "low_rate_c20"))

    dma_inputs(a.package, a.dest)
    n = sum(len(files) for _, _, files in os.walk(a.dest))
    print(f"{n} files in {a.dest}")


if __name__ == "__main__":
    main()
