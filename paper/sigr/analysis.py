"""Saving/loading reduced results, and the fit scores (model vs data)."""
import os

import numpy as np
import pandas as pd

from . import data
from .conditions import MODEL_CAPACITY_AH


# ---------------------------------------------------------------------------
# I/O: arrays in <key>_results.npz (no pickle), RPT curves in <key>_rpt_curves.csv
# ---------------------------------------------------------------------------
def save(res, curves, key, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    np.savez_compressed(os.path.join(out_dir, f"{key}_results.npz"), **res)
    pd.concat([pd.DataFrame({"rpt": i, "thr": c["thr"], "q": c["q"], "v": c["v"]})
               for i, c in enumerate(curves)]).to_csv(os.path.join(out_dir, f"{key}_rpt_curves.csv"), index=False)


def load(key, out_dir):
    with np.load(os.path.join(out_dir, f"{key}_results.npz"), allow_pickle=False) as z:
        res = {k: z[k] for k in z.files}
    df = pd.read_csv(os.path.join(out_dir, f"{key}_rpt_curves.csv"))
    curves = [dict(thr=float(d["thr"].iloc[0]), q=d["q"].to_numpy(), v=d["v"].to_numpy())
              for _, d in df.groupby("rpt", sort=True)]
    return res, curves


# ---------------------------------------------------------------------------
# Derived model quantities
# ---------------------------------------------------------------------------
def efc(res, thr):
    return np.asarray(thr, dtype=float) / (2 * float(res["efc_basis_ah"]))


def soh_reference(res):
    """Full windows: the first C/3 ageing discharge. Partial windows: the
    conditioned C/20 reference RPT (as the data's RPT1)."""
    i = int(res["soh_ref_index"])
    return res["rpt_cap"][i] if i >= 0 and res["rpt_cap"].size > i else res["age_cap"][0]


def rpt_soh(res):
    return 100 * res["rpt_cap"] / soh_reference(res)


def ageing_soh(res):
    return 100 * res["age_cap"] / res["age_cap"][0]


def lli_pct(res, with_lam_trapped=True):
    q = res["Q_side"] + (res["Q_lli_lam_trap"] if with_lam_trapped else 0)
    return 100 * q / MODEL_CAPACITY_AH


def life_mask(res):
    """Time-series points up to the last ageing cycle (the RPT after it excluded)."""
    return res["thr"] <= res["age_thr"][-1]


def match_rpt_curves(res, curves, cond, data_dir, max_gap=None):
    """[(real RPT number, real EFC, model curve index, model EFC)], each real
    C/20 discharge paired with the nearest unused model RPT."""
    c20 = data.load(cond, data_dir, "capacity_fade")
    c20 = c20[c20["source"] == "C/20 RPT"]
    real_efc = dict(zip(c20["rpt"].astype(int), c20["efc"]))
    model_efc = efc(res, [c["thr"] for c in curves])
    pairs, used = [], set()
    for r in data.rpt_discharge_numbers(cond, data_dir):
        order = np.argsort(np.abs(model_efc - real_efc[r]))
        j = next((int(k) for k in order if int(k) not in used), int(order[0]))
        used.add(j)
        if max_gap is not None and abs(model_efc[j] - real_efc[r]) > max_gap:
            continue
        pairs.append((r, float(real_efc[r]), j, float(model_efc[j])))
    return pairs


# ---------------------------------------------------------------------------
# Scores
# ---------------------------------------------------------------------------
def rpt_gaps(res, cond, data_dir):
    """Model C/20 RPT SoH interpolated at each real C/20 RPT EFC (model - real, pp)."""
    cap = data.load(cond, data_dir, "capacity_fade")
    c20 = cap[cap["source"] == "C/20 RPT"]
    e = efc(res, res["rpt_thr"])
    o = np.argsort(e)
    e, soh = e[o], rpt_soh(res)[o]
    rows = [(r.efc, r.capacity_retention_pct, float(np.interp(r.efc, e, soh)))
            for r in c20.itertuples() if e[0] <= r.efc <= e[-1]]
    return pd.DataFrame(rows, columns=["efc", "real_soh", "model_soh"]).assign(
        gap_pp=lambda d: d.model_soh - d.real_soh)


def knee(res, cond, data_dir, onset_soh=91.0):
    """Model knee = steepest descent of the C/3 ageing SoH; also the SoH-91% crossing."""
    soh, e = ageing_soh(res), efc(res, res["age_thr"])
    d = np.convolve(np.gradient(soh, e), np.ones(3) / 3.0, mode="same")
    i = int(np.argmin(d))
    cross = next((e[j - 1] + (soh[j - 1] - onset_soh) / (soh[j - 1] - soh[j]) * (e[j] - e[j - 1])
                  for j in range(1, soh.size) if soh[j - 1] >= onset_soh > soh[j]), np.nan)
    return dict(knee_efc=float(e[i]), soh_at_knee=float(soh[i]), soh91_efc=float(cross),
                real_knee_efc=data.knee_efc(cond, data_dir))


def voltage_rmse(res, curves, cond, data_dir):
    """V(Q) RMSE of each real C/20 RPT discharge against the nearest model RPT
    (pairs > 20 EFC apart skipped), over the overlapping capacity."""
    rows = []
    for r, e_r, j, e_m in match_rpt_curves(res, curves, cond, data_dir, max_gap=20.0):
        q_r, v_r = data.load_rpt_discharge(cond, data_dir, r)
        q_m, v_m = curves[j]["q"], curves[j]["v"]
        inr = q_r <= q_m.max()
        if inr.sum() >= 5:
            rows.append((r, e_r, e_m, float(np.sqrt(np.mean((np.interp(q_r[inr], q_m, v_m) - v_r[inr]) ** 2)))))
    return pd.DataFrame(rows, columns=["rpt", "real_efc", "model_efc", "rmse_V"])


DVDSOC_TRIM = 0.01  # the model's first ~1% DoD is the step-start transient (dropped, as in the plots)
DVDSOC_WINDOW = (0.02, 0.95)  # DoD range scored: excludes the end-of-discharge cliff


def smooth_dvdsoc(soc, v, n_bins=120):
    """Bin V on a uniform depth-of-discharge grid, then differentiate."""
    o = np.argsort(soc)
    soc, v = np.asarray(soc)[o], np.asarray(v)[o]
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    mid = 0.5 * (edges[:-1] + edges[1:])
    vb = np.array([v[(soc >= lo) & (soc < hi)].mean() if ((soc >= lo) & (soc < hi)).any() else np.nan
                   for lo, hi in zip(edges[:-1], edges[1:])])
    ok = ~np.isnan(vb)
    return (mid[ok], np.gradient(vb[ok], mid[ok])) if ok.sum() >= 3 else (mid[ok], np.full(ok.sum(), np.nan))


def dvdsoc_rmse(res, curves, cond, data_dir):
    """dV/dSOC RMSE [V per unit SOC] of each real C/20 RPT discharge against the
    nearest model RPT (pairs > 20 EFC apart skipped). Both curves are normalised
    to their own depth of discharge, binned and differentiated (smooth_dvdsoc),
    and compared over DVDSOC_WINDOW."""
    rows = []
    for r, e_r, j, e_m in match_rpt_curves(res, curves, cond, data_dir, max_gap=20.0):
        q_r, v_r = data.load_rpt_discharge(cond, data_dir, r)
        s_m, v_m = curves[j]["q"] / curves[j]["q"].max(), curves[j]["v"]
        keep = s_m >= DVDSOC_TRIM
        x_r, d_r = smooth_dvdsoc(q_r / q_r.max(), v_r)
        x_m, d_m = smooth_dvdsoc(s_m[keep], v_m[keep])
        lo, hi = DVDSOC_WINDOW
        sel = (x_r >= max(lo, x_m.min())) & (x_r <= min(hi, x_m.max()))
        if sel.sum() >= 5:
            rows.append((r, e_r, e_m, float(np.sqrt(np.mean((np.interp(x_r[sel], x_m, d_m) - d_r[sel]) ** 2)))))
    return pd.DataFrame(rows, columns=["rpt", "real_efc", "model_efc", "rmse_dvdsoc"])


def lli_lam_table(res, cond, data_dir):
    """Real vs model LLI (particle inventory) and LAM (Gr, Si, PE) at each real RPT."""
    lli = data.load(cond, data_dir, "LLI")
    lam = data.load(cond, data_dir, "LAM_derived")
    e_all = efc(res, res["thr"])
    m = life_mask(res)
    e_life = e_all[m]
    lli_model = 100 * res["Q_lli_particles"] / MODEL_CAPACITY_AH
    lli_real = 100 * lli["charge_LLI_loss"] / lli["charge_LLI"].iloc[0]
    rows = []
    for e, lr, (_, row) in zip(lli["efc"], lli_real, lam.iterrows()):
        if 0 < e <= e_all[-1]:
            rows.append(dict(
                efc=e, LLI_real=lr, LLI_model=np.interp(e, e_all, lli_model),
                LAM_gr_real=row["LAM_NE_graphite_pct"], LAM_gr_model=np.interp(row["efc"], e_life, res["LAM_gr"][m]),
                LAM_si_real=row["LAM_NE_silicon_pct"], LAM_si_model=np.interp(row["efc"], e_life, res["LAM_si"][m]),
                LAM_pe_real=row["LAM_PE_pct"], LAM_pe_model=np.interp(row["efc"], e_life, res["LAM_pos"][m])))
    return pd.DataFrame(rows)


def crack_table(res):
    """Si crack, porosity (minimum of the preceding ageing cycle) and SEI
    state at each model RPT."""
    e_all = efc(res, res["thr"])
    e_age = efc(res, res["age_thr"])
    to_pct = 100 / MODEL_CAPACITY_AH
    rows = []
    for e in efc(res, res["rpt_thr"]):
        i = int(np.argmin(np.abs(e_all - e)))
        before = e_age <= e
        porosity = res["age_min_porosity"][before][-1] if before.any() else float(res["initial_porosity"])
        rows.append(dict(efc=e, l_over_l0=res["crack_length"][i] / res["crack_length"][0],
                         roughness=res["roughness"][i], porosity=porosity,
                         lli_si_sei_pct=res["Q_si_sei"][i] * to_pct,
                         lli_si_sei_cracks_pct=res["Q_si_sei_cracks"][i] * to_pct,
                         strain_E=res["strain_E"][i], dl_paris_nm=1e9 * res["dl_paris"][i],
                         dl_strain_nm=1e9 * res["dl_strain"][i]))
    return pd.DataFrame(rows)


def summary(res, curves, cond, data_dir):
    """Headline scores for one condition (the paper's): mean |SoH gap| at the
    real C/20 RPTs [pp], and the mean C/20 voltage and dV/dSOC RMSEs."""
    gaps = rpt_gaps(res, cond, data_dir)
    vr = voltage_rmse(res, curves, cond, data_dir)
    dr = dvdsoc_rmse(res, curves, cond, data_dir)
    return dict(condition=cond.key, cell=cond.cell,
                mean_abs_soh_gap_pp=float(gaps["gap_pp"].abs().mean()),
                mean_voltage_rmse_V=float(vr["rmse_V"].mean()) if len(vr) else np.nan,
                mean_dvdsoc_rmse=float(dr["rmse_dvdsoc"].mean()) if len(dr) else np.nan)
