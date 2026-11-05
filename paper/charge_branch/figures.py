"""Supplementary figures: the C/20 RPT charge with the silicon lithiation OCP fixed (the paper runs) and
deformed (option "OCP aging deformation branches": "both"), against the measured charges (from paper/):

    python charge_branch/figures.py <runs, lithiation OCP fixed> <runs, lithiation OCP deformed>

Both run folders hold <key>_results.npz, <key>_rpt_curves.csv and <key>_rpt_charge_curves.csv (sigr.simulate
keeps the constant-current part of each RPT charge). Measured charges: data/<folder>/low_rate_c20/
<CELL>_RPTnnn_lowrate_c20_charge.csv (tools/import_c20_charges.py). Writes to charge_branch/figures/:
  si_charge_overlay_full, _partial   two cells each (rows): C/20 discharge (lithiation OCP fixed), charge
                         with the lithiation OCP fixed, charge with it deformed
  si_charge_<key>_voltage, si_charge_<key>_dvdsoc     each measured charge against the model RPT at the same EFC
and charge_branch/scores.csv (per RPT) and scores_summary.csv (per cell, as the Supplementary table).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.getcwd())
import sigr  # noqa: E402
from sigr import analysis as an, data, plots, style  # noqa: E402

DATA = "data"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
KEYS = ["25deg_full_soc", "45deg_full_soc", "25deg_15-95_soc", "25deg_20-80_soc"]
LABELS = ("Model, lithiation OCP fixed", "Model, lithiation OCP deformed")
COLORS = (plots.C_MODEL, "tab:red")
STYLES = ("--", ":")
MAX_GAP = 20.0  # EFC: measured charge and model RPT paired only if closer


def measured_charges(cond):
    """[(rpt, EFC since RPT1, Q from 0, V)] of the measured C/20 CC charges within the shifted data."""
    cap = data.load(cond, DATA, "capacity_fade")
    efc = dict(zip(cap["rpt"].astype(int), cap["efc"]))
    folder = os.path.join(DATA, cond.data_folder, "low_rate_c20")
    out = []
    for f in sorted(os.listdir(folder)):
        if f.startswith(f"{cond.cell}_RPT") and f.endswith("_lowrate_c20_charge.csv"):
            r = int(f[len(cond.cell) + 4:len(cond.cell) + 7])
            if r in efc:
                out.append((r, float(efc[r]), *data.load_rpt_charge(cond, DATA, r)))
    return out


def bol_charge(key, cache_dir):
    """The formation C/20 CC charge (simulate.bol_curves): its own cycle, so not kept by full-window runs.
    No LAM yet, so the same for both runs (the OCP deformation ramps with the LAM fraction)."""
    path = os.path.join(cache_dir, f"{key}_bol_charge.npz")
    if not os.path.exists(path):
        from sigr import simulate
        q, v = simulate.bol_curves(sigr.CONDITIONS[key], DATA)["charge"]
        np.savez(path, q=q, v=v)
    z = np.load(path)
    return z["q"], z["v"]


def load_run(folder, key):
    res, curves = an.load(key, folder)
    chg = an.load_charge(key, folder)
    if chg is None:
        raise FileNotFoundError(f"{folder}: no {key}_rpt_charge_curves.csv")
    if 0 not in chg and an.efc(res, curves[0]["thr"]) < 1.0:
        chg[0] = bol_charge(key, folder)
    return res, curves, chg


def pairs(run, meas):
    """[(measured entry, model RPT index)] nearest unused model RPT with a charge, within MAX_GAP."""
    res, curves, chg = run
    m_efc = an.efc(res, [c["thr"] for c in curves])
    used, out = set(), []
    for m in meas:
        for k in np.argsort(np.abs(m_efc - m[1])):
            if abs(m_efc[k] - m[1]) > MAX_GAP:
                break
            if int(k) in chg and int(k) not in used:
                used.add(int(k))
                out.append((m, int(k)))
                break
    return out


def charge_scores(q_r, v_r, q_m, v_m):
    """V(Q) RMSE [mV] over the overlapping capacity; dV/dSoC RMSE over an.DVDSOC_WINDOW (each charge
    normalised to its own CC capacity, the model's first 1 % dropped), as the discharge scores."""
    inr = q_r <= q_m.max()
    v_rmse = 1e3 * float(np.sqrt(np.mean((np.interp(q_r[inr], q_m, v_m) - v_r[inr]) ** 2)))
    s_m = q_m / q_m.max()
    keep = s_m >= an.DVDSOC_TRIM
    x_r, d_r = an.smooth_dvdsoc(q_r / q_r.max(), v_r)
    x_m, d_m = an.smooth_dvdsoc(s_m[keep], v_m[keep])
    lo, hi = an.DVDSOC_WINDOW
    sel = (x_r >= max(lo, x_m.min())) & (x_r <= min(hi, x_m.max()))
    return v_rmse, float(np.sqrt(np.mean((np.interp(x_r[sel], x_m, d_m) - d_r[sel]) ** 2)))


def charge_overlay(ax, run, meas, norm, cmap):
    """Model charges solid, measured dashed, coloured by EFC as the discharge panel (whose legend and colour
    bar they share)."""
    res, curves, chg = run
    for k in sorted(chg):
        ax.plot(*chg[k], color=cmap(norm(an.efc(res, curves[k]["thr"]))), lw=1.0)
    for (r, e_r, q_r, v_r), _ in pairs(run, meas):
        ax.plot(q_r, v_r, color=cmap(norm(e_r)), ls="--", lw=1.0)
    ax.set_ylabel("Terminal voltage [V]")
    ax.set_xlabel("Charge capacity, C/20 CC [A h]")


def overlay_figure(runs, meas, keys):
    fig, ax = plt.subplots(len(keys), 3, figsize=(style.DOUBLE_IN, 2.55 * len(keys)), layout="constrained")
    for row, key in enumerate(keys):
        cond = sigr.CONDITIONS[key]
        res0, curves0, _ = runs[0][key]
        plots.rpt_voltage(ax[row, 0], res0, curves0, cond, DATA)
        efcs = an.efc(res0, [c["thr"] for c in curves0])
        norm, cmap = plt.Normalize(efcs.min(), efcs.max()), plt.cm.viridis
        for i in (0, 1):
            charge_overlay(ax[row, 1 + i], runs[i][key], meas[key], norm, cmap)
        lo = min(a.get_ylim()[0] for a in ax[row])
        hi = max(a.get_ylim()[1] for a in ax[row])
        for i, (a, what) in enumerate(zip(ax[row], ("discharge", "charge, fixed", "charge, deformed"))):
            a.set_ylim(lo, hi)
            style.panel_label(a, 3 * row + i, f"{cond.cell} {what}")
    return fig


def dod_figure(runs, meas, key, kind):
    """As plots.rpt_dod, for the charge: each measured charge against both models, against SoC."""
    cond = sigr.CONDITIONS[key]
    p = [pairs(run[key], meas[key]) for run in runs]
    rpts = [m[0] for m, _ in p[0]]
    ncol = 2 if len(rpts) <= 4 else (3 if len(rpts) <= 6 else 4)
    nrow = int(np.ceil(len(rpts) / ncol))
    width = style.SINGLE_IN if ncol == 2 else style.DOUBLE_IN
    fig, axes = plt.subplots(nrow, ncol, figsize=(width, 1.45 * nrow + 0.45), sharex=True, sharey=True,
                             squeeze=False, layout="constrained")
    for i, (ax, ((r, e_r, q_r, v_r), _)) in enumerate(zip(axes.flat, p[0])):
        if kind == "voltage":
            ax.plot(q_r / q_r.max(), v_r, color="black", lw=1.3, label="Measured")
        else:
            ax.plot(*an.smooth_dvdsoc(q_r / q_r.max(), v_r), color="black", lw=1.3, label="Measured")
        for n, run in enumerate(runs):
            k = dict((m[0], j) for m, j in p[n]).get(r)
            if k is None:
                continue
            q_m, v_m = run[key][2][k]
            s_m = q_m / q_m.max()
            if kind == "voltage":
                ax.plot(s_m, v_m, STYLES[n], color=COLORS[n], lw=1.2, label=LABELS[n])
            else:
                keep = s_m >= an.DVDSOC_TRIM
                ax.plot(*an.smooth_dvdsoc(s_m[keep], v_m[keep]), STYLES[n], color=COLORS[n], lw=1.2, label=LABELS[n])
        if kind != "voltage":
            ax.set_ylim(-0.3, 2.5)
        style.panel_label(ax, i, f"RPT{r}, EFC {e_r:.0f}")
    for ax in list(axes.flat)[len(rpts):]:
        ax.set_visible(False)
    for ax in axes[:, 0]:
        ax.set_ylabel("Terminal voltage [V]" if kind == "voltage" else "dV/dSoC [V]")
    for c in range(ncol):
        col = [a for a in axes[:, c] if a.get_visible()]
        if col:
            col[-1].set_xlabel("State of charge [-]")
            col[-1].xaxis.set_tick_params(labelbottom=True)
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="outside lower center", ncol=3, frameon=False)
    return fig


def score_tables(runs, meas):
    rows, summ = [], []
    for key in KEYS:
        cond = sigr.CONDITIONS[key]
        s = dict(cell=cond.cell, condition=cond.label)
        for n, (name, run) in enumerate(zip(("fixed", "deformed"), runs)):
            res, curves, chg = run[key]
            sc = []
            for (r, e_r, q_r, v_r), k in pairs(run[key], meas[key]):
                v_rmse, d_rmse = charge_scores(q_r, v_r, *chg[k])
                sc.append((v_rmse, d_rmse))
                rows.append(dict(cell=cond.cell, lithiation_ocp=name, rpt=r, efc=e_r,
                                 model_efc=float(an.efc(res, curves[k]["thr"])), charge_v_rmse_mV=v_rmse,
                                 charge_dvdsoc_rmse=d_rmse))
            s[f"soh_gap_pp_{name}"] = float(an.rpt_gaps(res, cond, DATA)["gap_pp"].abs().mean())
            s[f"discharge_v_rmse_mV_{name}"] = 1e3 * float(an.voltage_rmse(res, curves, cond, DATA)["rmse_V"].mean())
            s[f"charge_v_rmse_mV_{name}"] = float(np.mean([a for a, _ in sc]))
            s[f"charge_dvdsoc_rmse_{name}"] = float(np.mean([b for _, b in sc]))
        summ.append(s)
    return pd.DataFrame(rows), pd.DataFrame(summ)


def main(fixed, deformed):
    style.use()
    runs = [{k: load_run(f, k) for k in KEYS} for f in (fixed, deformed)]
    meas = {k: measured_charges(sigr.CONDITIONS[k]) for k in KEYS}
    paths = style.save(overlay_figure(runs, meas, KEYS[:2]), OUT, "si_charge_overlay_full")
    paths += style.save(overlay_figure(runs, meas, KEYS[2:]), OUT, "si_charge_overlay_partial")
    for key in KEYS:
        for kind in ("voltage", "dvdsoc"):
            paths += style.save(dod_figure(runs, meas, key, kind), OUT, f"si_charge_{key}_{kind}")
    plt.close("all")
    rows, summ = score_tables(runs, meas)
    rows.to_csv(os.path.join(HERE, "scores.csv"), index=False, float_format="%.4g")
    summ.to_csv(os.path.join(HERE, "scores_summary.csv"), index=False, float_format="%.4g")
    pd.set_option("display.width", 200)
    print(summ.round(3).to_string())
    print("\n".join(p for p in paths if p.endswith(".png")))


if __name__ == "__main__":
    main(*[os.path.abspath(a) for a in sys.argv[1:3]])
