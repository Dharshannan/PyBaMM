"""worker_diffusivity_probe.py -- like worker_corrected_lli.py, but for the
2026-09-19 Si-diffusivity/LAM-expansion-residual theory test: runs the
current C64_* env config, reports the standard expansion-shape RMSE, AND
probes the model's own C/20 RPT legs for Si's stoichiometric window
(si_delta_sto) the same way dma/extract_model_si_sto_window_c20rpt.py does,
so a single run reports both "did the fit get better" and "did the root-
cause mechanism (sto-range collapse) actually shrink". Always produces the
full plot suite via plot_all() unless C64_NO_PLOT=1.

Also reports score_rpt_gap (model vs. real C/20 RPT SoH, mean |gap| in
percentage points) per the user's own flagged concern (2026-09-19): raising
Si diffusivity relieves diffusion polarisation, which could shift the
model's ALREADY-TUNED post-knee C/3 capacity-fade behaviour even though the
slow-rate C/20 RPT points themselves are expected to move much less.

Also tracks Si's C/3 AGEING-leg stoichiometric window (not just the C/20
RPT legs) -- added 2026-09-19 per user request, since the C/3 collapse is
what actually drives the fitted "reversible expansion" ptp metric, and the
baseline diagnostic run showed the C/3 collapse (17.7x) is notably WORSE
than the C/20 RPT collapse (9.3x) for the same config -- worth tracking
both across the sweep, not just the RPT one.

Prints one machine-parseable summary line at the end:
    [SWEEP_SUMMARY] tag=<tag> expansion_rmse=<f> si_delta_first=<f>
    si_delta_min=<f> collapse_ratio=<f> c3_delta_first=<f> c3_delta_min=<f>
    c3_collapse_ratio=<f> final_soh_pct=<f> rpt_gap_pp=<f>
"""
import os
import sys

import numpy as np

DEG_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, DEG_DIR)
import cell064_degradation_fit as c64  # noqa: E402

SI_STO_VAR = "X-averaged negative secondary particle surface stoichiometry"


def probe_si_sto_windows(sol):
    """Classify each cycle's discharge leg the same way
    cell064_degradation_fit's own rpt_leg()/cycle_ageing_leg() do, and
    return (rpt_deltas, c3_deltas) -- Si's own sto window for each kind."""
    rpt_deltas, c3_deltas = [], []
    for cyc in sol.cycles:
        for step in getattr(cyc, "steps", []):
            try:
                I = step["Current [A]"].entries
                sto = step[SI_STO_VAR].entries
            except (KeyError, TypeError, AttributeError):
                continue
            if I.size < 2 or sto.size < 2:
                continue
            mean_I = np.mean(I)
            if 0.05 < mean_I <= 0.5:
                rpt_deltas.append(float(sto.max() - sto.min()))
                break
            elif mean_I > 0.5:
                c3_deltas.append(float(sto.max() - sto.min()))
                break
    return np.array(rpt_deltas), np.array(c3_deltas)


def _collapse_stats(deltas):
    if deltas.size >= 2:
        first, worst = float(deltas[0]), float(deltas.min())
        ratio = first / worst if worst > 1e-9 else float("inf")
    else:
        first, worst, ratio = float("nan"), float("nan"), float("nan")
    return first, worst, ratio


def main():
    tag = os.environ.get("C64_OUT_TAG", "untagged")
    sol = c64.run_degradation()
    results = c64.extract_results(sol)
    exp_exp = c64.load_experimental_reversible_expansion()
    exp_cap = c64.load_experimental_capacity_fade()
    sim_efc_age = c64.efc_from_throughput(results["age_thr"])

    c64.plot_all(results, sol)

    rmse = c64.score_expansion_shape(results, exp_exp, sim_efc_age)
    rpt_gap = c64.score_rpt_gap(results, exp_cap)

    rpt_deltas, c3_deltas = probe_si_sto_windows(sol)
    si_first, si_min, collapse_ratio = _collapse_stats(rpt_deltas)
    c3_first, c3_min, c3_collapse_ratio = _collapse_stats(c3_deltas)

    final_soh = (
        100.0 * float(results["rpt_cap_arr"][-1]) / c64.NOMINAL_CAP_AH
        if results["rpt_cap_arr"].size
        else float("nan")
    )

    print(
        f"[SWEEP_SUMMARY] tag={tag} expansion_rmse={rmse if rmse is not None else float('nan'):.4f} "
        f"si_delta_first={si_first:.4f} si_delta_min={si_min:.4f} "
        f"collapse_ratio={collapse_ratio:.3f} c3_delta_first={c3_first:.4f} "
        f"c3_delta_min={c3_min:.4f} c3_collapse_ratio={c3_collapse_ratio:.3f} "
        f"final_soh_pct={final_soh:.2f} "
        f"rpt_gap_pp={rpt_gap if rpt_gap is not None else float('nan'):.2f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
