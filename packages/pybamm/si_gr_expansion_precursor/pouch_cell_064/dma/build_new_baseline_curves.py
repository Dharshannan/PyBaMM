"""build_new_baseline_curves.py -- generates the new-baseline's raw RPT4/
RPT5 V(Q) curves and LAM_Si/LAM_Gr/LLI trajectories (correct config, with
the porosity-gated isolation mechanism properly ON -- the earlier
build_writeup_comparison.py run had a bug where NEW_ENV silently inherited
OLD_ENV's isolation-disabled state). Combines with true_old_baseline_
rpt45_curves.csv (from build_true_old_baseline.py) for the final 3-way
(real / true old / new) comparison plots.
"""
import os
import sys

for k in list(os.environ):
    if k.startswith("C64_"):
        del os.environ[k]
os.environ["C64_SI_REDIRECT_TO_LAM"] = "1"
os.environ["C64_SI_BETA_LAM_ISO"] = "1.0"
os.environ["C64_SI_TAU_LAM_ISO"] = "2e8"
os.environ["C64_SI_REDIRECT_LAM_YIELD"] = "0.02"
os.environ["C64_SI_LAM_ISO_EXPONENT"] = "3.0"
os.environ["C64_SI_OCP_AGING_DEFORM"] = "1"
os.environ["C64_SI_VOLUME_CHANGE_AGING_DEFORM"] = "1"
os.environ["C64_NEG_POROSITY_FLOOR"] = "0.035"
os.environ["C64_EXPONENT_MAX_SEI"] = "70"
os.environ["C64_SI_BETA_LAM_SEI"] = "1e-7"
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

HERE = os.path.dirname(os.path.abspath(__file__))
DEG_DIR = os.path.join(os.path.dirname(HERE), "degradation_test_matrix")
sys.path.insert(0, DEG_DIR)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import cell064_degradation_fit as c64  # noqa: E402


def nearest_curve(sim_efc_rpt, curves, target_efc):
    if not curves or sim_efc_rpt.size == 0:
        return None, None
    idx = int(np.argmin(np.abs(sim_efc_rpt - target_efc)))
    c = curves[idx]
    return c["q"], c["v"]


def main():
    print("Config check (new baseline):", flush=True)
    print(f"  radius_mult={c64.SI_PARTICLE_RADIUS_MULT}, diff_mult={c64.SI_DIFFUSIVITY_MULT}, "
          f"BOL={c64.NEG_POROSITY_BOL}, SI_K={c64.SI_K_SEI_MULT}, GR_K={c64.GR_K_SEI_MULT}, "
          f"F0={c64.F0_BASELINE}, WIDTH={c64.WIDTH_BASELINE}, "
          f"SI_LAM_OPTION={c64.SI_LAM_OPTION}, redirect={c64.SEI_REDIRECT_TO_LAM_OPTION}", flush=True)

    sol = c64.run_degradation()
    results = c64.extract_results(sol)
    sim_efc_rpt = c64.efc_from_throughput(results["rpt_thr"]) if results["rpt_thr"].size else np.array([])
    sim_efc_full = c64.efc_from_throughput(results["Qt"])
    curves = results["rpt_voltage_curves"]
    real_efc = {4: 171.4, 5: 197.4}

    rows = []
    for rpt_num, efc in real_efc.items():
        q, v = nearest_curve(sim_efc_rpt, curves, efc)
        if q is None:
            continue
        for qi, vi in zip(q, v):
            rows.append({"rpt": rpt_num, "q": qi, "v": vi})
    pd.DataFrame(rows).to_csv(os.path.join(HERE, "new_baseline_rpt45_curves.csv"), index=False)
    print(f"Saved: new_baseline_rpt45_curves.csv ({len(curves)} RPT curves total)", flush=True)

    # LAM/LLI trajectories
    lam_df = pd.DataFrame({
        "efc": sim_efc_full,
        "LAM_si": results["LAM_si"],
        "LAM_gr": results["LAM_gr"],
        "lli_ah": results["Q_side"] + results["Q_lli_lam_trap"],
    })
    lam_df.to_csv(os.path.join(HERE, "new_baseline_lam_lli_trajectory.csv"), index=False)
    print("Saved: new_baseline_lam_lli_trajectory.csv", flush=True)
    print(f"DONE, final SoH={100 * results['rpt_cap_arr'][-1] / c64.NOMINAL_CAP_AH:.1f}%", flush=True)


if __name__ == "__main__":
    main()
