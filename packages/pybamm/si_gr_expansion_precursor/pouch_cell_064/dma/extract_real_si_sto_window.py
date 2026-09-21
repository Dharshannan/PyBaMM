"""extract_real_si_sto_window.py -- pulls CELL064's REAL per-RPT Si (and Gr)
stoichiometric-window trend directly from the project's own authoritative
eSOH/DMA fit outputs, to check whether the PyBaMM model's severe post-knee
Si stoichiometry-range collapse (~14.3x, delta_sto 0.9405 -> 0.0659) matches
reality.

Data source (NOT a fresh fit -- reuses Pouch_Data's own existing, already
-validated `discharge_esoh_fit.py` results): the raw per-point RPT discharge
timeseries CSVs already carry a `model_x_si` / `model_x_gr` column, i.e. the
eSOH-fit's own reconstructed Si/Gr stoichiometry at every point of the
measured discharge curve. Si's per-RPT stoichiometric window is simply
max(model_x_si) - min(model_x_si) over one full discharge.

RPT1/2/4/5: low_rate_c20 (C/20 RPT discharge, the primary/most reliable
    source per the project's own eSOH summary CSV).
RPT0/3: no C/20 RPT exists for these (see the project's own
    CELL064_capacity_fade.csv, source="high-rate neighbor substitute") --
    use the high_rate_neighbor discharge CSVs instead (averaging
    before/after variants when both exist).

RPT->EFC anchors are taken verbatim from the project's own
`degradation_test_matrix/experimental_data/CELL064_capacity_fade.csv` so
this lines up exactly with everything else already fitted for CELL064.
"""
import glob
import os

import numpy as np
import pandas as pd

POUCH_DATA = (
    r"C:\Shannan_PhD_Stuff_Local\Si_Gr_Expansion_Precursor\Pouch_Data\cell064_data"
)
LOW_RATE_DIR = os.path.join(POUCH_DATA, "data_timeseries", "low_rate_c20")
HIGH_RATE_DIR = os.path.join(POUCH_DATA, "data_timeseries", "high_rate_neighbor")

HERE = os.path.dirname(os.path.abspath(__file__))
EXP_DATA = os.path.join(
    os.path.dirname(HERE), "degradation_test_matrix", "experimental_data"
)

# RPT -> EFC, verbatim from the project's own real-data anchor table.
EFC_BY_RPT = {}
_cap_fade = pd.read_csv(os.path.join(EXP_DATA, "CELL064_capacity_fade.csv"))
for _, row in _cap_fade.iterrows():
    EFC_BY_RPT[int(row["rpt"])] = float(row["efc"])
SOURCE_BY_RPT = {int(r["rpt"]): r["source"] for _, r in _cap_fade.iterrows()}


def window_from_trace(df, col):
    x = df[col].to_numpy(dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan, np.nan, np.nan
    return float(x.max()), float(x.min()), float(x.max() - x.min())


def extract_low_rate(rpt):
    path = os.path.join(
        LOW_RATE_DIR, f"CELL064_RPT{rpt:03d}_lowrate_c20_discharge.csv"
    )
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    si_max, si_min, si_delta = window_from_trace(df, "model_x_si")
    gr_max, gr_min, gr_delta = window_from_trace(df, "model_x_gr")
    return {
        "rpt": rpt,
        "source": "low_rate_c20 (C/20 RPT)",
        "n_files": 1,
        "x_si_max": si_max,
        "x_si_min": si_min,
        "x_si_delta": si_delta,
        "x_gr_max": gr_max,
        "x_gr_min": gr_min,
        "x_gr_delta": gr_delta,
    }


def extract_high_rate(rpt):
    pattern = os.path.join(HIGH_RATE_DIR, f"CELL064_RPT{rpt:03d}_*highrate_discharge.csv")
    files = sorted(glob.glob(pattern))
    if not files:
        return None
    si_deltas, gr_deltas, si_maxs, si_mins, gr_maxs, gr_mins = [], [], [], [], [], []
    for f in files:
        df = pd.read_csv(f)
        si_max, si_min, si_delta = window_from_trace(df, "model_x_si")
        gr_max, gr_min, gr_delta = window_from_trace(df, "model_x_gr")
        si_maxs.append(si_max)
        si_mins.append(si_min)
        si_deltas.append(si_delta)
        gr_maxs.append(gr_max)
        gr_mins.append(gr_min)
        gr_deltas.append(gr_delta)
    return {
        "rpt": rpt,
        "source": f"high_rate_neighbor (avg of {len(files)})",
        "n_files": len(files),
        "x_si_max": float(np.mean(si_maxs)),
        "x_si_min": float(np.mean(si_mins)),
        "x_si_delta": float(np.mean(si_deltas)),
        "x_gr_max": float(np.mean(gr_maxs)),
        "x_gr_min": float(np.mean(gr_mins)),
        "x_gr_delta": float(np.mean(gr_deltas)),
    }


def main():
    rows = []
    for rpt in range(0, 6):
        rec = extract_low_rate(rpt)
        if rec is None:
            rec = extract_high_rate(rpt)
        if rec is None:
            print(f"RPT{rpt}: no discharge data found, skipping")
            continue
        rec["efc"] = EFC_BY_RPT.get(rpt, np.nan)
        rec["efc_source"] = SOURCE_BY_RPT.get(rpt, "")
        rows.append(rec)

    out = pd.DataFrame(rows)
    out = out[
        [
            "rpt",
            "efc",
            "efc_source",
            "source",
            "n_files",
            "x_si_max",
            "x_si_min",
            "x_si_delta",
            "x_gr_max",
            "x_gr_min",
            "x_gr_delta",
        ]
    ]
    # normalise to RPT1 (earliest available, pre-knee) as the baseline, since
    # that's what the model-side collapse ratio (14.3x) was measured against
    ref = out.loc[out["rpt"] == 1, "x_si_delta"]
    ref_val = float(ref.iloc[0]) if len(ref) else np.nan
    out["x_si_delta_norm_to_rpt1"] = out["x_si_delta"] / ref_val
    ref_gr = out.loc[out["rpt"] == 1, "x_gr_delta"]
    ref_gr_val = float(ref_gr.iloc[0]) if len(ref_gr) else np.nan
    out["x_gr_delta_norm_to_rpt1"] = out["x_gr_delta"] / ref_gr_val

    out_path = os.path.join(HERE, "real_si_gr_sto_window_by_rpt.csv")
    out.to_csv(out_path, index=False)
    print(out.to_string(index=False))
    print(f"\nSaved: {out_path}")

    collapse_ratio = ref_val / float(out["x_si_delta"].min())
    print(
        f"\nReal Si delta_sto: RPT1={ref_val:.4f} -> min observed="
        f"{out['x_si_delta'].min():.4f} (collapse ratio {collapse_ratio:.2f}x "
        f"across real life)"
    )


if __name__ == "__main__":
    main()
