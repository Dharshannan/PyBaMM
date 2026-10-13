"""Experimental data for one condition, on the model's EFC axis (EFC since
the cell's RPT1; RPT1 is the 100% / zero reference).

Layout under the data directory (the Zenodo data record)::

    <data_folder>/<CELL>_capacity_fade.csv         rpt, efc, discharge_capacity_Ah, source, capacity_retention_pct
    <data_folder>/<CELL>_LLI.csv                   rpt, efc, charge_LLI, joule_LLI, charge_LLI_loss, joule_LLI_loss
    <data_folder>/<CELL>_LAM_derived.csv           rpt, efc, ..., LAM_NE_graphite_pct, LAM_NE_silicon_pct, LAM_PE_pct
    <data_folder>/<CELL>_reversible_expansion.csv  cyc_seq, efc, reversible_expansion_um
    <data_folder>/<CELL>_k_expansion_scale.csv     rpt, efc, k, source
    <data_folder>/low_rate_c20/<CELL>_RPT<nnn>_lowrate_c20_discharge.csv   Q_discharge_Ah, voltage_terminal_V
    ocp/*.csv                                      beginning-of-life OCPs (parameters.py)

EFC = throughput / (2 x 2.5 Ah) throughout.
"""
import os
import re

import pandas as pd


def _folder(cond, data_dir):
    return os.path.join(data_dir, cond.data_folder)


def rpt1_efc(cond, data_dir):
    cap = pd.read_csv(os.path.join(_folder(cond, data_dir), f"{cond.cell}_capacity_fade.csv"))
    return float(cap.loc[cap["rpt"] == 1, "efc"].iloc[0])


def load(cond, data_dir, name):
    """One data table, EFC shifted to 'since RPT1' and earlier rows dropped."""
    df = pd.read_csv(os.path.join(_folder(cond, data_dir), f"{cond.cell}_{name}.csv"))
    df["efc"] = df["efc"] - rpt1_efc(cond, data_dir)
    return df[df["efc"] >= 0].reset_index(drop=True)


def load_all(cond, data_dir):
    return {name: load(cond, data_dir, name)
            for name in ("capacity_fade", "LLI", "LAM_derived", "reversible_expansion", "k_expansion_scale")}


def rpt_targets(cond, data_dir):
    """EFCs (since RPT1, excluding RPT1 itself) at which the model runs its
    C/20 RPTs: every real RPT, C/20 or high-rate substitute."""
    cap = load(cond, data_dir, "capacity_fade")
    return sorted(float(e) for e in cap["efc"] if e > 0.5)


def knee_efc(cond, data_dir):
    return cond.knee_efc_raw - rpt1_efc(cond, data_dir)


def rpt_discharge_numbers(cond, data_dir):
    """Real RPTs with a C/20 discharge curve, within the shifted data."""
    d = os.path.join(_folder(cond, data_dir), "low_rate_c20")
    nums = sorted(int(m.group(1)) for f in os.listdir(d)
                  if (m := re.match(rf"{cond.cell}_RPT(\d{{3}})_lowrate_c20_discharge\.csv$", f)))
    c20 = load(cond, data_dir, "capacity_fade")
    c20 = set(c20.loc[c20["source"] == "C/20 RPT", "rpt"].astype(int))
    return [n for n in nums if n in c20]


def load_rpt_discharge(cond, data_dir, rpt):
    """(Q from 0 [A.h], V [V]) of a real C/20 RPT discharge."""
    path = os.path.join(_folder(cond, data_dir), "low_rate_c20",
                        f"{cond.cell}_RPT{rpt:03d}_lowrate_c20_discharge.csv")
    df = pd.read_csv(path, usecols=["Q_discharge_Ah", "voltage_terminal_V"])
    q = df["Q_discharge_Ah"].to_numpy()
    return q - q[0], df["voltage_terminal_V"].to_numpy()
