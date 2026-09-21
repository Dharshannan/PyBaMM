"""worker_full_retune_probe.py -- 2026-09-19 phase 2: once radius/diffusivity
(fixed at the caller's chosen best combo, e.g. radius x0.2 + diffusivity x3)
have addressed the stoichiometric-collapse root cause, this worker scores
the FOUR things that still need re-tuning against real CELL064 data now
that the underlying electrochemistry has changed:
    1. knee timing        (score_knee, |EFC error|)
    2. capacity-fade gap  (score_rpt_gap, mean |pp| at real RPTs)
    3. LLI                (NEW: mean |Ah error| vs CELL064_LLI.csv's
                            charge_LLI_loss, at matched real RPT EFCs)
    4. LAM_Gr / LAM_Si split (NEW: mean |pp error| for each vs
                            CELL064_LAM_derived.csv's LAM_NE_graphite_pct /
                            LAM_NE_silicon_pct, at matched real RPT EFCs)

Also still reports expansion RMSE and the C/3/C/20 sto-collapse ratios for
continuity with the earlier sweeps, but this phase's PRIMARY targets are
the four above.

Prints one machine-parseable summary line:
    [RETUNE2_SUMMARY] tag=<tag> expansion_rmse=<f> c3_collapse_ratio=<f>
    rpt_collapse_ratio=<f> knee_err_efc=<f> rpt_gap_pp=<f> lli_err_ah=<f>
    lam_gr_err_pp=<f> lam_si_err_pp=<f> final_soh_pct=<f>
"""
import os
import sys

import numpy as np

DEG_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, DEG_DIR)
import cell064_degradation_fit as c64  # noqa: E402

SI_STO_VAR = "X-averaged negative secondary particle surface stoichiometry"
F_CONST = 96485.33212  # Faraday constant [C/mol]


def probe_si_sto_windows(sol):
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


def _collapse_ratio(deltas):
    if deltas.size < 2:
        return float("nan")
    worst = deltas.min()
    return deltas[0] / worst if worst > 1e-9 else float("inf")


def score_lli(sol, results, exp_cap):
    """Mean |Ah error| between model's corrected LLI (Q_side + LAM-trapped
    lithium in both electrodes, same convention as worker_corrected_lli.py)
    and real CELL064_LLI.csv's charge_LLI_loss, at each real C/20 RPT's
    matched EFC."""
    cutoff_idx = len(results["Qt"])

    def full_res(name):
        return sol[name].entries[:cutoff_idx]

    lli_lam_neg_primary_Ah = full_res(
        "Loss of lithium due to loss of primary active material in negative electrode [mol]") * F_CONST / 3600.0
    lli_lam_neg_secondary_Ah = full_res(
        "Loss of lithium due to loss of secondary active material in negative electrode [mol]") * F_CONST / 3600.0
    lli_lam_pos_Ah = full_res(
        "Loss of lithium due to loss of active material in positive electrode [mol]") * F_CONST / 3600.0
    lli_corrected_Ah = (
        results["Q_side"] + lli_lam_neg_primary_Ah + lli_lam_neg_secondary_Ah + lli_lam_pos_Ah
    )
    Qt_efc = c64.efc_from_throughput(results["Qt"])

    try:
        lli_real = c64.load_experimental_lli()
    except Exception as e:  # pragma: no cover
        print(f"LLI score: unavailable ({e})", flush=True)
        return None
    in_range = (lli_real["efc"] >= Qt_efc[0]) & (lli_real["efc"] <= Qt_efc[-1])
    if not in_range.any():
        print("LLI score: unavailable (no EFC overlap)", flush=True)
        return None
    model_interp = np.interp(lli_real.loc[in_range, "efc"], Qt_efc, lli_corrected_Ah)
    errs = model_interp - lli_real.loc[in_range, "charge_LLI_loss"].to_numpy()
    mean_abs_err = float(np.mean(np.abs(errs)))
    print(f"LLI: mean |error| = {mean_abs_err:.4f} Ah across {int(in_range.sum())} real RPT(s)", flush=True)
    return mean_abs_err


def score_lam_split(results, exp_lam):
    """Mean |pp error| for LAM_Gr and LAM_Si separately, model vs. real
    CELL064_LAM_derived.csv, at matched real RPT EFCs."""
    Qt_efc = c64.efc_from_throughput(results["Qt"])
    in_range = (exp_lam["efc"] >= Qt_efc[0]) & (exp_lam["efc"] <= Qt_efc[-1])
    if not in_range.any():
        print("LAM split score: unavailable (no EFC overlap)", flush=True)
        return None, None
    real = exp_lam.loc[in_range]
    model_gr = np.interp(real["efc"], Qt_efc, results["LAM_gr"])
    model_si = np.interp(real["efc"], Qt_efc, results["LAM_si"])
    err_gr = float(np.mean(np.abs(model_gr - real["LAM_NE_graphite_pct"].to_numpy())))
    err_si = float(np.mean(np.abs(model_si - real["LAM_NE_silicon_pct"].to_numpy())))
    print(f"LAM_Gr: mean |error| = {err_gr:.2f}pp, LAM_Si: mean |error| = {err_si:.2f}pp "
          f"across {int(in_range.sum())} real RPT(s)", flush=True)
    return err_gr, err_si


def main():
    tag = os.environ.get("C64_OUT_TAG", "untagged")
    sol = c64.run_degradation()
    results = c64.extract_results(sol)
    exp_exp = c64.load_experimental_reversible_expansion()
    exp_cap = c64.load_experimental_capacity_fade()
    exp_lam = c64.load_experimental_lam()
    sim_efc_age = c64.efc_from_throughput(results["age_thr"])

    c64.plot_all(results, sol)

    rmse = c64.score_expansion_shape(results, exp_exp, sim_efc_age)
    rpt_gap = c64.score_rpt_gap(results, exp_cap)
    knee_err = c64.score_knee(results, sim_efc_age)
    lli_err = score_lli(sol, results, exp_cap)
    lam_gr_err, lam_si_err = score_lam_split(results, exp_lam)

    rpt_deltas, c3_deltas = probe_si_sto_windows(sol)
    rpt_collapse = _collapse_ratio(rpt_deltas)
    c3_collapse = _collapse_ratio(c3_deltas)

    final_soh = (
        100.0 * float(results["rpt_cap_arr"][-1]) / c64.NOMINAL_CAP_AH
        if results["rpt_cap_arr"].size
        else float("nan")
    )

    def fmt(x):
        return float("nan") if x is None else x

    print(
        f"[RETUNE2_SUMMARY] tag={tag} expansion_rmse={fmt(rmse):.4f} "
        f"c3_collapse_ratio={c3_collapse:.3f} rpt_collapse_ratio={rpt_collapse:.3f} "
        f"knee_err_efc={fmt(knee_err):.2f} rpt_gap_pp={fmt(rpt_gap):.2f} "
        f"lli_err_ah={fmt(lli_err):.4f} lam_gr_err_pp={fmt(lam_gr_err):.2f} "
        f"lam_si_err_pp={fmt(lam_si_err):.2f} final_soh_pct={final_soh:.2f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
