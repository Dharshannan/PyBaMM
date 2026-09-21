"""extract_new_baseline_sto_collapse.py -- same per-cycle Si sto-window
extraction as extract_model_si_sto_window_c20rpt.py, but for the NEW
adopted baseline (pbuf2_f0p7_w0p03: radius x0.2, diffusivity x3, BOL=0.17,
SI_K_SEI=3e-4/GR_K_SEI=6e-6, F0=0.7/width=0.03), so
dma_si_sto_window_comparison.png can show old vs. new vs. real together.
EFC_LIMIT=205 per memory/cell064-efc-limit-205.md.
"""
import os
import sys

os.environ["C64_SI_REDIRECT_TO_LAM"] = "1"
os.environ["C64_SI_OCP_AGING_DEFORM"] = "1"
os.environ["C64_SI_VOLUME_CHANGE_AGING_DEFORM"] = "1"
os.environ["C64_SI_REDIRECT_LAM_YIELD"] = "0.02"
os.environ["C64_NEG_POROSITY_FLOOR"] = "0.035"
os.environ["C64_EXPONENT_MAX_SEI"] = "70"
os.environ["C64_SI_BETA_LAM_SEI"] = "1e-7"
os.environ["C64_SI_BETA_LAM_ISO"] = "1.0"
os.environ["C64_SI_LAM_ISO_EXPONENT"] = "3.0"
os.environ["C64_SI_TAU_LAM_ISO"] = "2e8"
os.environ["C64_GR_REDIRECT_TO_LAM"] = "0"
os.environ["C64_GR_REDIRECT_LAM_YIELD"] = "0.0"
os.environ["C64_MAX_CYCLES"] = "700"
os.environ["C64_SOH_FLOOR"] = "30"
os.environ["C64_EFC_LIMIT"] = "205"
os.environ["C64_NO_PLOT"] = "1"
# new baseline
os.environ["C64_SI_PARTICLE_RADIUS_MULT"] = "0.2"
os.environ["C64_SI_DIFFUSIVITY_MULT"] = "3"
os.environ["C64_NEG_POROSITY_BOL"] = "0.17"
os.environ["C64_SI_K_SEI_MULT"] = "3e-4"
os.environ["C64_GR_K_SEI_MULT"] = "6e-6"
os.environ["C64_F0"] = "0.7"
os.environ["C64_WIDTH"] = "0.03"

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
    out_path = os.path.join(HERE, "model_si_sto_window_new_baseline.csv")
    df.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path}")

    rpt_rows = df[df["kind"] == "rpt_c20"]
    if len(rpt_rows):
        ref = rpt_rows["si_delta_sto"].iloc[0]
        worst = rpt_rows["si_delta_sto"].min()
        print(f"[MODEL RPT] si_delta_sto: first={ref:.4f} -> min={worst:.4f} "
              f"(collapse ratio {ref / worst:.2f}x)")
    c3_rows = df[df["kind"] == "ageing_c3"]
    if len(c3_rows):
        ref = c3_rows["si_delta_sto"].iloc[0]
        worst = c3_rows["si_delta_sto"].min()
        print(f"[MODEL C3] si_delta_sto: first={ref:.4f} -> min={worst:.4f} "
              f"(collapse ratio {ref / worst:.2f}x)")


if __name__ == "__main__":
    main()
