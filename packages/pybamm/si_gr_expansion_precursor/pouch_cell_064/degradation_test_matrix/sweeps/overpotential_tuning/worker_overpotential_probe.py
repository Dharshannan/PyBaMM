"""worker_overpotential_probe.py -- 2026-09-21: post-knee (and slightly
pre-knee) C/20 RPT discharge overpotential is visibly overpredicted by the
model vs. the real cell (voltage sags more than it should). This worker is
worker_full_retune_probe.py's usual non-regression scorer set (knee timing,
capacity-fade RPT gap, LLI, LAM_Gr/LAM_Si split) PLUS TWO voltage-shape
metrics, kept side by side deliberately (per feedback: a fix that improves
one by quietly making the other worse is over-fitting, not a real fix):
    - voltage_rmse: c64.score_voltage_shape() -- raw discharge-capacity
      match. Conflates true overpotential/resistance shape with capacity-
      fade-RATE mismatch (a model that over-predicts capacity at an aged
      RPT can "improve" this metric without fixing anything real).
    - voltage_rmse_soc: c64.score_voltage_shape_soc() -- SOC-normalised
      (matched depth-of-discharge FRACTION, not raw Ah) match. Isolates
      the resistance/overpotential shape from the capacity mismatch.
    - rpt_gap_pp (below, unchanged) is the existing capacity-fade-fit
      tripwire -- keep watching it alongside both voltage metrics so a
      candidate that "fixes" voltage_rmse by over-predicting capacity
      again gets caught here, not missed.

Prints one machine-parseable summary line:
    [OVERPOT_SUMMARY] tag=<tag> voltage_rmse=<f> rpt1_rmse=<f> rpt2_rmse=<f>
    rpt4_rmse=<f> rpt5_rmse=<f> voltage_rmse_soc=<f> rpt1_rmse_soc=<f>
    rpt2_rmse_soc=<f> rpt4_rmse_soc=<f> rpt5_rmse_soc=<f>
    expansion_rmse=<f> knee_err_efc=<f> rpt_gap_pp=<f> lli_err_ah=<f>
    lam_gr_err_pp=<f> lam_si_err_pp=<f> final_soh_pct=<f>

rpt<N>_rmse[_soc] are parsed back out of the scorers' own per-RPT print
lines (they don't return them individually) so the sweep can see whether a
fix is uniform across RPTs or concentrated post-knee (RPT4/5) as expected.
"""
import io
import os
import re
import sys
from contextlib import redirect_stdout

import numpy as np

DEG_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, DEG_DIR)
import cell064_degradation_fit as c64  # noqa: E402

F_CONST = 96485.33212  # Faraday constant [C/mol]

RPT_RMSE_RE = re.compile(r"^\s*RPT(\d+) \(EFC [\d.]+, model EFC [\d.]+\): RMSE=(\S+) V")


def score_voltage_shape_verbose(results, exp_cap):
    """Wraps c64.score_voltage_shape(), also echoing its per-RPT prints to
    our own stdout (so the log still shows them) while capturing the
    per-RPT RMSEs it doesn't return directly."""
    buf = io.StringIO()
    with redirect_stdout(buf):
        mean_rmse = c64.score_voltage_shape(results, exp_cap)
    text = buf.getvalue()
    print(text, end="", flush=True)
    per_rpt = {}
    for line in text.splitlines():
        m = RPT_RMSE_RE.match(line)
        if m:
            per_rpt[int(m.group(1))] = float(m.group(2))
    return mean_rmse, per_rpt


def score_voltage_shape_soc_verbose(results, exp_cap):
    """Same wrapping as score_voltage_shape_verbose(), for the SOC-
    normalised metric (c64.score_voltage_shape_soc())."""
    buf = io.StringIO()
    with redirect_stdout(buf):
        mean_rmse = c64.score_voltage_shape_soc(results, exp_cap)
    text = buf.getvalue()
    print(text, end="", flush=True)
    per_rpt = {}
    for line in text.splitlines():
        m = RPT_RMSE_RE.match(line)
        if m:
            per_rpt[int(m.group(1))] = float(m.group(2))
    return mean_rmse, per_rpt


def score_lli(sol, results, exp_cap):
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


def dump_rpt45_curves(results, exp_cap, tag):
    """Save the model's own RPT4/RPT5 V(Q) curves (nearest-EFC matched to
    real, same convention as score_voltage_shape) to CSV, for building a
    dedicated zoomed-in comparison plot afterward -- the full summary
    plot's voltage panel is at full 2.5-4.25V axis scale, where a
    real-but-modest (tens of mV) improvement is hard to see by eye."""
    model_curves = results["rpt_voltage_curves"]
    if not model_curves:
        return
    model_efcs = c64.efc_from_throughput(np.array([c["thr"] for c in model_curves]))
    c20 = exp_cap[exp_cap["source"] == "C/20 RPT"]
    real_efc_by_rpt = dict(zip(c20["rpt"], c20["efc"]))
    rows = []
    for rpt_num in (1, 2, 4, 5):
        if rpt_num not in real_efc_by_rpt:
            continue
        real_efc = real_efc_by_rpt[rpt_num]
        j = int(np.argmin(np.abs(model_efcs - real_efc)))
        c = model_curves[j]
        for qi, vi in zip(c["q"], c["v"]):
            rows.append({"rpt": rpt_num, "q": qi, "v": vi})
    out_dir = os.path.dirname(os.path.abspath(__file__))
    out_path = os.path.join(out_dir, f"rpt45_curves_{tag}.csv")
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("rpt,q,v\n")
        for r in rows:
            fh.write(f"{r['rpt']},{r['q']},{r['v']}\n")
    print(f"Saved: {out_path}", flush=True)


def main():
    tag = os.environ.get("C64_OUT_TAG", "untagged")
    sol = c64.run_degradation()
    results = c64.extract_results(sol)
    exp_exp = c64.load_experimental_reversible_expansion()
    exp_cap = c64.load_experimental_capacity_fade()
    exp_lam = c64.load_experimental_lam()
    sim_efc_age = c64.efc_from_throughput(results["age_thr"])

    c64.plot_all(results, sol)
    dump_rpt45_curves(results, exp_cap, tag)

    rmse = c64.score_expansion_shape(results, exp_exp, sim_efc_age)
    rpt_gap = c64.score_rpt_gap(results, exp_cap)
    knee_err = c64.score_knee(results, sim_efc_age)
    lli_err = score_lli(sol, results, exp_cap)
    lam_gr_err, lam_si_err = score_lam_split(results, exp_lam)
    voltage_rmse, per_rpt_rmse = score_voltage_shape_verbose(results, exp_cap)
    voltage_rmse_soc, per_rpt_rmse_soc = score_voltage_shape_soc_verbose(results, exp_cap)

    final_soh = (
        100.0 * float(results["rpt_cap_arr"][-1]) / c64.NOMINAL_CAP_AH
        if results["rpt_cap_arr"].size
        else float("nan")
    )

    def fmt(x):
        return float("nan") if x is None else x

    print(
        f"[OVERPOT_SUMMARY] tag={tag} voltage_rmse={fmt(voltage_rmse):.5f} "
        f"rpt1_rmse={fmt(per_rpt_rmse.get(1)):.5f} rpt2_rmse={fmt(per_rpt_rmse.get(2)):.5f} "
        f"rpt4_rmse={fmt(per_rpt_rmse.get(4)):.5f} rpt5_rmse={fmt(per_rpt_rmse.get(5)):.5f} "
        f"voltage_rmse_soc={fmt(voltage_rmse_soc):.5f} "
        f"rpt1_rmse_soc={fmt(per_rpt_rmse_soc.get(1)):.5f} "
        f"rpt2_rmse_soc={fmt(per_rpt_rmse_soc.get(2)):.5f} "
        f"rpt4_rmse_soc={fmt(per_rpt_rmse_soc.get(4)):.5f} "
        f"rpt5_rmse_soc={fmt(per_rpt_rmse_soc.get(5)):.5f} "
        f"expansion_rmse={fmt(rmse):.4f} knee_err_efc={fmt(knee_err):.2f} "
        f"rpt_gap_pp={fmt(rpt_gap):.2f} lli_err_ah={fmt(lli_err):.4f} "
        f"lam_gr_err_pp={fmt(lam_gr_err):.2f} lam_si_err_pp={fmt(lam_si_err):.2f} "
        f"final_soh_pct={final_soh:.2f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
