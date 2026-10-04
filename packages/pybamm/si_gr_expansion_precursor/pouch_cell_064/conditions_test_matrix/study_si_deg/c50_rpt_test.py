# -*- coding: utf-8 -*-
"""c50_rpt_test: is the late-life 25 degC C/20 voltage offset a lumped-model
porosity-polarisation artefact?

Theory: a lumped (pseudo-2D, x-averaged gate) model concentrates pore
closure uniformly, so it overpredicts electrolyte-transport polarisation
compared with the real cell. If the model's late-life offset is mostly that
polarisation, a slower C/50 discharge from the same state should bring the
model curve up towards the real C/20 curve.

Method: R4c (crack_baseline_v2 + WIDTH 5e-3) at 25 degC, unchanged.
After formation and after every 50-cycle batch, the fork's BATCH_HOOK
branches off one "Discharge at C/50 until 2.5 V" from the same fully
charged state (end of the CV hold) that the next C/20 RPT starts from.
The branch is never fed back, so the degradation trajectory is exactly R4c.

Outputs (this folder, TAG = c50_test_25degC_R4c):
    c50_voltage_dod_<TAG>.png   real C/20 / model C/20 / model C/50 / model OCV vs DoD, per real RPT
    c50_dvdsoc_dod_<TAG>.png    the same as dV/d(SOC) (model first 1% DoD trimmed)
    c50_polarisation_<TAG>.png  mean (OCV - V) at C/20 and C/50 vs EFC
    c50_curves_<TAG>.csv        every model curve (rate, efc, q, v, ocv), for replotting

Usage (from this folder):  python c50_rpt_test.py [--replot]
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TAG = "c50_test_25degC_R4c"
CSV = os.path.join(HERE, f"c50_curves_{TAG}.csv")
OCV_VAR = "Battery open-circuit voltage [V]"

os.environ.update(C64_T_AMBIENT_K="298.15", C64_F0="0.7", C64_WIDTH="5e-3", C64_OUT_TAG=TAG)
with open(os.path.join(HERE, "crack_baseline_v2.env")) as fh:
    for line in fh:
        line = line.split("#", 1)[0].strip()
        if "=" in line:
            key, val = line.split("=", 1)
            os.environ.setdefault(key.strip(), val.strip())
os.environ.setdefault("C64_DIAG_EXTRA_VARS", OCV_VAR)

sys.path.insert(0, os.path.join(os.path.dirname(HERE), "high_temp_45C"))
os.chdir(HERE)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pybamm  # noqa: E402

import cell064_degradation_fit_45C as m  # noqa: E402

m.SCRIPT_DIR = HERE
DVDSOC_TRIM = 0.01
branches = []  # C/50 side discharges collected by the hook


def c50_hook(last_sol, batch_num, make_sim):
    efc = float(m.efc_from_throughput(last_sol["Throughput capacity [A.h]"].entries[-1]))
    exp = pybamm.Experiment([f"Discharge at C/50 until {m.LOWER_CUTOFF_V} V"])
    try:
        side = make_sim(exp).solve(starting_solution=last_sol)
    except Exception as exc:  # noqa: BLE001
        print(f"[C/50 branch {batch_num}] failed at EFC {efc:.1f}: {type(exc).__name__}: {str(exc)[:150]}",
              flush=True)
        return
    step = side.cycles[-1].steps[-1]
    q = step["Discharge capacity [A.h]"].entries
    branches.append(dict(efc=efc, q=q - q[0], v=step["Terminal voltage [V]"].entries.copy(),
                         ocv=step[OCV_VAR].entries.copy()))
    print(f"[C/50 branch {batch_num}] EFC {efc:.1f}: Q = {q[-1] - q[0]:.3f} Ah", flush=True)
    del side


def simulate():
    m.BATCH_HOOK = c50_hook
    sol = m.run_degradation()
    results = m.extract_results(sol)
    os.environ["C64_NO_PLOT"] = "1"
    m.plot_all(results, sol)  # scores only: must reproduce R4c (gap 1.5pp, knee 128.8)

    rows = []
    for cyc in sol.cycles:  # the main trajectory's own C/20 RPT discharges, with OCV
        cap, _, d_end, _ = m.rpt_leg(cyc)
        if cap <= 0 or d_end is None:
            continue
        best = None
        for step in cyc.steps:
            I = step["Current [A]"].entries
            if I.size >= 2 and 0.05 < np.mean(I) <= 0.5:
                q = step["Discharge capacity [A.h]"].entries
                if best is None or q[-1] - q[0] > best[0][-1] - best[0][0]:
                    best = (q, step)
        q, step = best
        efc = float(m.efc_from_throughput(cyc["Throughput capacity [A.h]"].entries[-1]))
        rows.append(pd.DataFrame(dict(rate="C/20", efc=efc, q=q - q[0],
                                      v=step["Terminal voltage [V]"].entries, ocv=step[OCV_VAR].entries)))
    for b in branches:
        rows.append(pd.DataFrame(dict(rate="C/50", **b)))
    pd.concat(rows).to_csv(CSV, index=False)
    print("Saved:", CSV, flush=True)


def smooth_dvdsoc(soc, v, n_bins=120):
    order = np.argsort(soc)
    soc, v = np.asarray(soc)[order], np.asarray(v)[order]
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    mid = 0.5 * (edges[:-1] + edges[1:])
    vb = np.array([v[(soc >= lo) & (soc < hi)].mean() if ((soc >= lo) & (soc < hi)).any() else np.nan
                   for lo, hi in zip(edges[:-1], edges[1:])])
    ok = ~np.isnan(vb)
    return (mid[ok], np.gradient(vb[ok], mid[ok])) if ok.sum() >= 3 else (mid[ok], np.full(ok.sum(), np.nan))


def plot():
    df = pd.read_csv(CSV)
    curves = {rate: [dict(efc=g["efc"].iloc[0], q=g["q"].to_numpy(), v=g["v"].to_numpy(), ocv=g["ocv"].to_numpy())
                     for _, g in d.groupby("efc", sort=True)] for rate, d in df.groupby("rate")}
    c20, c50 = curves["C/20"], curves["C/50"]
    exp_cap = m.load_experimental_capacity_fade()
    real = exp_cap[exp_cap["source"] == "C/20 RPT"]
    real_efc = dict(zip(real["rpt"].astype(int), real["efc"]))

    def nearest(lst, e):
        return lst[int(np.argmin([abs(c["efc"] - e) for c in lst]))]

    # Polarisation (OCV - V), mean over 5-95% DoD, per model RPT
    pol = lambda c: float(np.mean((c["ocv"] - c["v"])[(c["q"] / c["q"].max() > 0.05) & (c["q"] / c["q"].max() < 0.95)]))  # noqa: E731
    print("\nModel RPT: EFC, Qmax C/20 vs C/50 [Ah], mean polarisation C/20 vs C/50 [mV]:", flush=True)
    for c in c20:
        b = nearest(c50, c["efc"])
        print(f"  EFC {c['efc']:6.1f} (C/50 branch EFC {b['efc']:6.1f}): Q {c['q'].max():.3f} / {b['q'].max():.3f}"
              f"  pol {1e3 * pol(c):5.1f} / {1e3 * pol(b):5.1f}", flush=True)

    rpts = [r for r in m.REAL_RPT_NUMS_WITH_DISCHARGE if r in real_efc]
    for kind in ("voltage", "dvdsoc"):
        fig, axes = plt.subplots(2, 2, figsize=(11, 9))
        for ax, r in zip(axes.flat, rpts):
            q_r, v_r = m.load_real_rpt_discharge(r)
            a, b = nearest(c20, real_efc[r]), nearest(c50, real_efc[r])
            sets = [("Real C/20", q_r / q_r.max(), v_r, None, dict(color="black", lw=2)),
                    (f"Model C/20 (EFC~{a['efc']:.0f})", a["q"] / a["q"].max(), a["v"], a["ocv"],
                     dict(color="tab:blue", ls="--", lw=1.8)),
                    (f"Model C/50 (EFC~{b['efc']:.0f})", b["q"] / b["q"].max(), b["v"], b["ocv"],
                     dict(color="tab:green", lw=1.6))]
            for lab, soc, v, ocv, st in sets:
                if kind == "voltage":
                    ax.plot(soc, v, label=lab, **st)
                else:
                    keep = soc >= (DVDSOC_TRIM if lab.startswith("Model") else 0.0)
                    ax.plot(*smooth_dvdsoc(soc[keep], v[keep]), label=lab, **st)
            if kind == "voltage":
                ax.plot(b["q"] / b["q"].max(), b["ocv"], ":", color="grey", lw=1.4, label="Model OCV (C/50 branch)")
                ax.set_ylabel("Terminal voltage [V]")
            else:
                ax.set_ylabel("dV/d(SOC) [V]")
                ax.set_ylim(-3, 0.5)
            ax.set_title(f"RPT{r} (real EFC~{real_efc[r]:.0f})")
            ax.set_xlabel("Depth of discharge [-]")
            ax.legend(fontsize=8)
            ax.grid(alpha=0.3)
        fig.suptitle(f"CELL064 25 degC, R4c: model C/20 vs C/50 RPT from the same state "
                     f"({'discharge' if kind == 'voltage' else 'dV/d(SOC)'})")
        fig.tight_layout()
        out = os.path.join(HERE, f"c50_{kind}_dod_{TAG}.png")
        fig.savefig(out, dpi=150)
        print("Saved:", out, flush=True)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot([c["efc"] for c in c20], [1e3 * pol(c) for c in c20], "o-", color="tab:blue", label="Model C/20")
    ax.plot([c["efc"] for c in c50], [1e3 * pol(c) for c in c50], "s-", color="tab:green", label="Model C/50")
    ax.set_xlabel("EFC [-]")
    ax.set_ylabel("Mean polarisation OCV - V over 5-95% DoD [mV]")
    ax.set_title("R4c 25 degC: RPT polarisation at C/20 vs C/50")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    out = os.path.join(HERE, f"c50_polarisation_{TAG}.png")
    fig.savefig(out, dpi=150)
    print("Saved:", out, flush=True)


if __name__ == "__main__":
    if "--replot" not in sys.argv or not os.path.exists(CSV):
        simulate()
    plot()
