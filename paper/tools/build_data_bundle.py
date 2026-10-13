"""Assemble the data directory the paper code reads (the Zenodo data record)
from the study's local files. For the author; readers download the bundle.

    python paper/tools/build_data_bundle.py [--dest paper/data] [--pouch-data <Pouch_Data dir>]
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
TABLES = ("capacity_fade", "LLI", "LAM_derived", "reversible_expansion", "k_expansion_scale")


def copy(src, dst_dir):
    os.makedirs(dst_dir, exist_ok=True)
    shutil.copy2(src, dst_dir)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dest", default=os.path.join(REPO, "paper", "data"))
    ap.add_argument("--pouch-data", default=POUCH_DATA, help="source of CELL064's C/20 RPT discharge curves")
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
        for f in [os.path.join(EXP, folder, f"{cell}_{t}.csv") for t in TABLES] + \
                 glob.glob(os.path.join(EXP, folder, f"{cell}_cycling_protocol.json")):
            copy(f, os.path.join(a.dest, folder))
        for f in sorted(glob.glob(os.path.join(EXP, folder, "low_rate_c20", "*_discharge.csv"))):
            copy(f, os.path.join(a.dest, folder, "low_rate_c20"))

    n = sum(len(files) for _, _, files in os.walk(a.dest))
    print(f"{n} files in {a.dest}")


if __name__ == "__main__":
    main()
