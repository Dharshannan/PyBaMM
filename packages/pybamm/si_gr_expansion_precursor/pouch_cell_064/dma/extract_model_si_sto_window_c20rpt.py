"""extract_model_si_sto_window_c20rpt.py -- same CELL064 degradation run as
extract_model_si_sto_window.py, but this time extracts Si's stoichiometric
window from the model's own periodic C/20 RPT legs (built into
cell064_degradation_fit.py's ageing_batch_exp every RPT_INTERVAL cycles,
classified the same way the fit script's own rpt_leg() does: 0.05 < mean_I
<= 0.5 A) instead of the C/3 ageing legs.

This is the true apples-to-apples analogue of the REAL x_si_delta numbers
in real_si_gr_sto_window_by_rpt.csv, which also come from a slow C/20 full
discharge -- removing the C/3-vs-C/20 rate-mismatch caveat noted when
comparing against extract_model_si_sto_window.py's ageing-leg numbers.

Also records the ageing-leg (C/3) delta_sto alongside each RPT for direct
side-by-side comparison in one file.
"""
import os
import sys

os.environ["C64_SI_REDIRECT_TO_LAM"] = "1"
os.environ["C64_SI_BETA_LAM_ISO"] = "1.0"
os.environ["C64_SI_TAU_LAM_ISO"] = "2e8"
os.environ["C64_SI_K_SEI_MULT"] = "6e-4"
os.environ["C64_GR_K_SEI_MULT"] = "1.25e-5"
os.environ["C64_NEG_POROSITY_BOL"] = "0.130"
os.environ["C64_SI_OCP_AGING_DEFORM"] = "1"
os.environ["C64_SI_REDIRECT_LAM_YIELD"] = "0.02"
os.environ["C64_NEG_POROSITY_FLOOR"] = "0.035"
os.environ["C64_SI_VOLUME_CHANGE_AGING_DEFORM"] = "1"
os.environ["C64_SI_VOLUME_CHANGE_EXP_BOL"] = "1.561216606871758"
os.environ["C64_SI_VOLUME_CHANGE_EXP_END"] = "2.18807468386046"
os.environ["C64_MAX_CYCLES"] = "700"
os.environ["C64_SOH_FLOOR"] = "30"
os.environ["C64_EFC_LIMIT"] = "260"
os.environ["C64_NO_PLOT"] = "1"

DEG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + \
    os.sep + "degradation_test_matrix"
sys.path.insert(0, DEG_DIR)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import cell064_degradation_fit as c64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

SI_STO_VAR = "X-averaged negative secondary particle surface stoichiometry"
GR_STO_VAR = "X-averaged negative primary particle surface stoichiometry"
EPS_S_VAR = (
    "X-averaged negative electrode secondary active material volume fraction"
)


def _sto_window_from_step(step, var_name):
    try:
        x = step[var_name].entries
    except (KeyError, TypeError, AttributeError):
        return None
    if x.size < 2:
        return None
    return float(x.min()), float(x.max()), float(x.max() - x.min())


def probe_cycle(cyc):
    """Classify each step of this cycle the same way cell064_degradation_fit's
    own cycle_ageing_leg()/rpt_leg() do, and return whichever kind of
    discharge leg (ageing C/3, or RPT C/20) is present, plus Si/Gr sto
    windows for that leg."""
    result = {"kind": None, "efc": None}
    for step in getattr(cyc, "steps", []):
        try:
            I = step["Current [A]"].entries
            thr = step["Throughput capacity [A.h]"].entries
        except (KeyError, TypeError, AttributeError):
            continue
        if I.size < 2:
            continue
        mean_I = np.mean(I)
        if mean_I > 0.5:
            kind = "ageing_c3"
        elif 0.05 < mean_I <= 0.5:
            kind = "rpt_c20"
        else:
            continue
        si_win = _sto_window_from_step(step, SI_STO_VAR)
        gr_win = _sto_window_from_step(step, GR_STO_VAR)
        try:
            eps_s = step[EPS_S_VAR].entries
            eps_s_min = float(eps_s.min())
        except (KeyError, TypeError, AttributeError):
            eps_s_min = None
        if si_win is None:
            continue
        return {
            "kind": kind,
            "efc": float(thr[-1]) / (2.0 * c64.NOMINAL_CAP_AH),
            "si_sto_min": si_win[0],
            "si_sto_max": si_win[1],
            "si_delta_sto": si_win[2],
            "gr_delta_sto": gr_win[2] if gr_win else None,
            "eps_s_min": eps_s_min,
        }
    return result


def main():
    sol = c64.run_degradation()
    print(f"[MODEL] {len(sol.cycles)} total cycles", flush=True)
    rows = []
    for i, cyc in enumerate(sol.cycles):
        rec = probe_cycle(cyc)
        if rec.get("kind") is None:
            continue
        rec["cycle_index"] = i
        rows.append(rec)
        if rec["kind"] == "rpt_c20":
            print(
                f"[MODEL RPT] cycle {i}: efc={rec['efc']:.1f}, "
                f"si_delta_sto={rec['si_delta_sto']:.4f}, "
                f"eps_s_min={rec['eps_s_min']:.5f}",
                flush=True,
            )

    df = pd.DataFrame(rows)
    out_path = os.path.join(HERE, "model_si_sto_window_c3_and_c20rpt.csv")
    df.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path}")

    rpt_rows = df[df["kind"] == "rpt_c20"]
    if len(rpt_rows):
        print("\n=== Model's own C/20 RPT legs (apples-to-apples with real data) ===")
        print(rpt_rows[["cycle_index", "efc", "si_delta_sto", "gr_delta_sto", "eps_s_min"]].to_string(index=False))
        ref = rpt_rows["si_delta_sto"].iloc[0]
        worst = rpt_rows["si_delta_sto"].min()
        print(
            f"\n[MODEL RPT] si_delta_sto: first RPT={ref:.4f} -> min="
            f"{worst:.4f} (collapse ratio {ref / worst:.2f}x)"
        )
    else:
        print("\nNo RPT-classified (C/20) legs found in solution -- check "
              "RPT_INTERVAL/BATCH_SIZE alignment.")


if __name__ == "__main__":
    main()
