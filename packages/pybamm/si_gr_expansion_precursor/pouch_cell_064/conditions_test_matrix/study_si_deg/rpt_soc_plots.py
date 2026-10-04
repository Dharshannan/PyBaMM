# -*- coding: utf-8 -*-
"""rpt_soc_plots: run the adopted per-temperature configs and plot every real
C/20 RPT against the model's matched RPT, on a depth-of-discharge axis. The
layout follows the DMA FINDINGS_2026-09-20 SOC-normalised figures.

Configs are the current baseline, crack_baseline_v3.env (= R7 / R7a; set
RPT_BASELINE_ENV=crack_baseline_v2.env for the previous one), plus, per
temperature, only:
    25 degC: T = 298.15 K, F0 = 0.70, WIDTH = 5e-3
    45 degC: T = 318.15 K, F0 = 0.80, WIDTH = 3e-4
Every other parameter is shared. Run the baseline with --aligned.

Outputs (in this folder), with TAG = v3_25degC / v3_45degC (or RPT_TAG):
    the fork's standard plots (overall_summary etc.) for TAG
    rpt_voltage_dod_<TAG>.png   V vs depth of discharge, one panel per real RPT
    rpt_dvdsoc_dod_<TAG>.png    dV/d(SOC) vs depth of discharge, same panels
    rpt_curves_<TAG>.csv        the model RPT curves (rpt, efc, q, v), for replotting

Both curves are normalised to their own discharged capacity (DoD 0-1) before
plotting and differentiating, so a capacity mismatch doesn't distort the
shape comparison. SoH gaps are tracked separately by the summary plot.
Model RPTs are matched to real RPTs by nearest EFC, each model RPT used once
(the same rule as the fork's score_voltage_shape).

The model's first 1% of DoD (a step-start numerical transient) is trimmed
from the dV/d(SOC) panels only.

Usage (from this folder):  python rpt_soc_plots.py 25|45 [--replot] [--aligned]
--replot skips the simulation and replots from rpt_curves_<TAG>.csv.
--aligned runs the model RPTs at the real RPT EFCs (aligned_rpt_cycling.py)
instead of every 50 cycles.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# Output folder (plots, CSV, pickle, cwd). Default: this folder; RPT_OUT_DIR
# redirects a study's runs into its own subfolder (e.g. si_volume_cracks/).
OUT = os.path.abspath(os.environ.get("RPT_OUT_DIR", HERE))
os.makedirs(OUT, exist_ok=True)
T_ARG = sys.argv[1] if len(sys.argv) > 1 else "25"
CONFIGS = {
    "25": dict(C64_T_AMBIENT_K="298.15", C64_F0="0.7", C64_WIDTH="5e-3", tag="v3_25degC"),
    "45": dict(C64_T_AMBIENT_K="318.15", C64_F0="0.8", C64_WIDTH="3e-4", tag="v3_45degC"),
    # Partial-SoC-window cells at 25 degC (the 25 degC values: F0 = RPT1 k = 0.70,
    # WIDTH 5e-3), overlaid on their own data via the fork's C64_REAL_CELL_DIR.
    # Ageing control per the cells' *_cycling_protocol.json:
    "25p1595": dict(C64_T_AMBIENT_K="298.15", C64_F0="0.7", C64_WIDTH="5e-3",
                    C64_REAL_CELL_DIR="CELL009_25C_15-95", tag="v3_CELL009_15-95_25degC"),
    "25p2080": dict(C64_T_AMBIENT_K="298.15", C64_F0="0.7", C64_WIDTH="5e-3",
                    C64_REAL_CELL_DIR="CELL026_25C_20-80", tag="v3_CELL026_20-80_25degC"),
}
# Partial-window ageing cycles (model C/3 = Q_nom/3, as for every other run):
#   CELL009 15-95%: voltage-limited at both ends (C/3 CC-CV to 4.15 V, CV to
#     14 mA; C/3 CC discharge to 3.14 V).
#   CELL026 20-80%: CC-CV to 3.969 V (CV to 14 mA) on top; coulomb-counted
#     1.500 Ah C/3 discharge at the bottom, with the 2.5 V safety floor.
PARTIAL_AGEING = {
    "25p1595": ("Discharge at C/3 until 3.14 V", "Charge at C/3 until 4.15 V", "Hold at 4.15 V until 14 mA"),
    "25p2080": (f"Discharge at C/3 for {1.5 / (2.5947 / 3) * 3600:.1f} seconds or until 2.5 V",
                "Charge at C/3 until 3.969 V", "Hold at 3.969 V until 14 mA"),
}
cfg = CONFIGS[T_ARG]
TAG = cfg.pop("tag") + os.environ.get("RPT_TAG_SUFFIX", "")  # suffix for variant runs
TAG = os.environ.get("RPT_TAG", TAG)  # or a full override, e.g. R5a_25degC
for _k, _v in cfg.items():  # per-condition values win over the env, unless overridden as RPT_<key> (e.g. RPT_C64_WIDTH)
    os.environ[_k] = os.environ.get("RPT_" + _k, _v)
with open(os.path.join(HERE, os.environ.get("RPT_BASELINE_ENV", "crack_baseline_v3.env"))) as fh:
    for line in fh:
        line = line.split("#", 1)[0].strip()
        if "=" in line:
            key, val = line.split("=", 1)
            os.environ.setdefault(key.strip(), val.strip())
os.environ["C64_OUT_TAG"] = TAG
LLI_VAR = "Total lithium lost from particles [mol]"
CRACK_VARS = [
    "X-averaged negative secondary particle crack length [m]",
    "X-averaged negative secondary electrode roughness ratio",
    "X-averaged negative secondary particle surface tangential stress [Pa]",
    "X-averaged negative electrode porosity",
    "Loss of lithium to negative secondary SEI [mol]",
    "Loss of lithium to negative secondary SEI on cracks [mol]",
]
# Route-B (strain-fatigue cracking) diagnostics exist only when it is on.
STRAIN_VARS = [
    "X-averaged negative secondary particle strain-fatigue driver [s-1]",
    "X-averaged negative secondary particle strain-fatigue cracking rate [m.s-1]",
    "X-averaged negative secondary particle cracking rate [m.s-1]",
] if os.environ.get("C64_SI_STRAIN_CRACK_KV") else []
os.environ.setdefault("C64_DIAG_EXTRA_VARS", "|".join([LLI_VAR] + CRACK_VARS + STRAIN_VARS))

sys.path.insert(0, os.path.join(os.path.dirname(HERE), "high_temp_45C"))
os.chdir(OUT)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import cell064_degradation_fit_45C as m  # noqa: E402

m.SCRIPT_DIR = OUT
CELL = (m.REAL_CELL_DIR_NAME.split("_")[0] if m.REAL_CELL_DIR_NAME
        else "CELL064" if T_ARG == "25" else "CELL017")
T_LABEL = "45" if T_ARG == "45" else "25"


def real(name):
    """This condition's real data set (EFC since RPT1)."""
    if m.REAL_CELL_DIR_NAME:
        return m.load_real_cell(name)
    if T_ARG == "45":
        return m.load_cell017(name)
    return {"capacity_fade": m.load_experimental_capacity_fade, "LLI": m.load_experimental_lli,
            "LAM_derived": m.load_experimental_lam}[name]()
CSV = os.path.join(OUT, f"rpt_curves_{TAG}.csv")
DVDSOC_TRIM = 0.01  # model DoD below this is dropped from dV/d(SOC) only


def smooth_dvdsoc(soc, v, n_bins=120):
    """Bin V onto a uniform DoD grid, then differentiate (same as the DMA
    findings' plot_soc_4rpt_dvdq_configs.py)."""
    order = np.argsort(soc)
    soc, v = np.asarray(soc)[order], np.asarray(v)[order]
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    mid = 0.5 * (edges[:-1] + edges[1:])
    v_binned = np.array([v[(soc >= lo) & (soc < hi)].mean() if ((soc >= lo) & (soc < hi)).any() else np.nan
                         for lo, hi in zip(edges[:-1], edges[1:])])
    ok = ~np.isnan(v_binned)
    mid, v_binned = mid[ok], v_binned[ok]
    if mid.size < 3:
        return mid, np.full(mid.size, np.nan)
    return mid, np.gradient(v_binned, mid)


def _use_cell017_rpts():
    """What plot_all does at 45 degC: point the real-RPT loaders at CELL017."""
    m.REAL_RPT_DISCHARGE_DIR = os.path.join(m.CELL017_DIR, "low_rate_c20")
    m.REAL_RPT_FILE_PREFIX = "CELL017"
    m.REAL_RPT_NUMS_WITH_DISCHARGE = tuple(int(r) for r in m.load_cell017("capacity_fade")["rpt"])


def print_lli_lam_checks(sol, results):
    """Particle-inventory LLI and Gr/Si LAM vs the real DMA at each RPT
    (same comparison as crack_growth_sweep.py)."""
    efc_full = m.efc_from_throughput(sol["Throughput capacity [A.h]"].entries)
    lli = 100 * sol[LLI_VAR].entries * 96485.33212 / 3600 / m.NOMINAL_CAP_AH
    exp_lli = real("LLI")
    exp_lam = real("LAM_derived")
    efc_lam = m.efc_from_throughput(results["Qt"])
    print("LLI check (real / model particle inventory) and LAM check (Gr, Si real / model):", flush=True)
    lli_real = 100 * exp_lli["charge_LLI_loss"] / exp_lli["charge_LLI"].iloc[0]
    for (e, r), (_, row) in zip(zip(exp_lli["efc"], lli_real), exp_lam.iterrows()):
        if 0 < e <= efc_full[-1]:
            print(f"  EFC={e:.1f}: LLI real={r:.1f}% model={np.interp(e, efc_full, lli):.1f}% | "
                  f"Gr real={row['LAM_NE_graphite_pct']:.1f}% model={np.interp(row['efc'], efc_lam, results['LAM_gr']):.1f}% | "
                  f"Si real={row['LAM_NE_silicon_pct']:.1f}% model={np.interp(row['efc'], efc_lam, results['LAM_si']):.1f}% | "
                  f"PE real={row['LAM_PE_pct']:.1f}% model={np.interp(row['efc'], efc_lam, results['LAM_pos']):.1f}%",
                  flush=True)


def plot_crack_diagnostics(sol, results):
    """Si crack growth, roughness, peak tangential stress per cycle, porosity
    and Si SEI LLI (bulk vs cracks) vs EFC, plus a table at each model RPT."""
    try:
        v = {name: np.asarray(sol[name].entries, dtype=float) for name in CRACK_VARS}
    except Exception as exc:  # noqa: BLE001
        print(f"crack diagnostics unavailable: {exc}", flush=True)
        return
    efc = m.efc_from_throughput(sol["Throughput capacity [A.h]"].entries)
    l_cr, rough, sig, eps, lli_b, lli_c = (v[n] for n in CRACK_VARS)
    to_pct = 96485.33212 / 3600 / m.NOMINAL_CAP_AH * 100
    # peak |tangential stress| per cycle (cycles can be long; take the max over each)
    cyc_efc, cyc_sig = [], []
    for cyc in sol.cycles:
        try:
            sc = np.asarray(cyc[CRACK_VARS[2]].entries, dtype=float)
            cyc_efc.append(float(m.efc_from_throughput(cyc["Throughput capacity [A.h]"].entries[-1])))
            cyc_sig.append(float(np.max(np.abs(sc))) / 1e6)
        except Exception:  # noqa: BLE001
            continue
    print("Crack diagnostics at each model RPT (EFC, l/l0, roughness, porosity, Si SEI LLI bulk / cracks %):", flush=True)
    for e in m.efc_from_throughput(np.asarray(results["rpt_thr"], dtype=float)):
        i = int(np.argmin(np.abs(efc - e)))
        print(f"  EFC {e:6.1f}: l/l0 {l_cr[i] / l_cr[0]:7.3f}  rough {rough[i]:6.3f}  eps {eps[i]:.4f}  "
              f"LLI Si bulk {lli_b[i] * to_pct:5.2f}%  cracks {lli_c[i] * to_pct:5.2f}%", flush=True)
    strain = None
    if STRAIN_VARS:
        try:
            t = sol["Time [s]"].entries
            drv, dl_v, dl_tot = (np.asarray(sol[n].entries, dtype=float) for n in STRAIN_VARS)
            cum = lambda y: np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(t))])  # noqa: E731
            strain = dict(E=cum(drv), dl_v=cum(dl_v), dl_a=cum(dl_tot - dl_v))
            print("Route-B diagnostics at each model RPT (EFC, cumulative Si strain E, "
                  "crack growth from Paris / strain [nm]):", flush=True)
            for e in m.efc_from_throughput(np.asarray(results["rpt_thr"], dtype=float)):
                i = int(np.argmin(np.abs(efc - e)))
                print(f"  EFC {e:6.1f}: E {strain['E'][i]:8.2f}  dl Paris {1e9 * strain['dl_a'][i]:8.2f}  "
                      f"dl strain {1e9 * strain['dl_v'][i]:8.2f}", flush=True)
        except Exception as exc:  # noqa: BLE001
            print(f"route-B diagnostics unavailable: {exc}", flush=True)
    fig, ax = plt.subplots(3 if strain else 2, 2, figsize=(11, 11 if strain else 7.5))
    ax[0, 0].plot(efc, l_cr / l_cr[0]); ax[0, 0].set_ylabel("Si crack length / initial [-]")
    ax2 = ax[0, 0].twinx(); ax2.plot(efc, rough, color="tab:red", lw=1); ax2.set_ylabel("roughness ratio", color="tab:red")
    ax[0, 1].plot(cyc_efc, cyc_sig, ".", ms=3); ax[0, 1].set_ylabel("peak |Si surface tangential stress| per cycle [MPa]")
    ax[1, 0].plot(efc, eps); ax[1, 0].set_ylabel("X-averaged negative porosity [-]")
    ax[1, 1].plot(efc, lli_b * to_pct, label="Si SEI (bulk)"); ax[1, 1].plot(efc, lli_c * to_pct, label="Si SEI on cracks")
    ax[1, 1].set_ylabel("LLI [% of nominal]"); ax[1, 1].legend()
    if strain:
        ax[2, 0].plot(efc, strain["E"]); ax[2, 0].set_ylabel("cumulative Si strain driver E [-]")
        ax[2, 1].plot(efc, 1e9 * strain["dl_a"], label="route A (Paris)")
        ax[2, 1].plot(efc, 1e9 * strain["dl_v"], label="route B (strain fatigue)")
        ax[2, 1].set_ylabel("cumulative crack growth [nm]"); ax[2, 1].legend()
    for a in ax.flat:
        a.set_xlabel("EFC [-]"); a.grid(alpha=0.3)
    fig.suptitle(f"Si crack / porosity diagnostics ({TAG})")
    fig.tight_layout()
    out = os.path.join(OUT, f"crack_diagnostics_{TAG}.png")
    fig.savefig(out, dpi=130)
    print("Saved:", out, flush=True)


def simulate():
    print(f"Config {TAG}:", {k: os.environ[k] for k in sorted(os.environ) if k.startswith("C64_")}, flush=True)
    if "--aligned" in sys.argv:  # RPTs at the real RPT EFCs (aligned_rpt_cycling.py)
        import aligned_rpt_cycling as arc
        sol = arc.run_degradation_aligned(m, arc.real_rpt_efcs(m, T_ARG), PARTIAL_AGEING.get(T_ARG))
    else:
        sol = m.run_degradation()
    results = m.extract_results(sol)
    import pickle
    with open(os.path.join(OUT, f"results_{TAG}.pkl"), "wb") as fh:  # for re-plotting without re-simulating
        pickle.dump(results, fh)
    curves = list(results["rpt_voltage_curves"])  # plot_all can empty results lists
    m.plot_all(results, sol)  # standard summary plots + scores
    print_lli_lam_checks(sol, results)
    try:
        plot_crack_diagnostics(sol, results)
    except Exception as exc:  # noqa: BLE001 -- a diagnostic must never cost the run's outputs
        print(f"crack diagnostics failed: {type(exc).__name__}: {exc}", flush=True)
    model_efc = m.efc_from_throughput(np.array([c["thr"] for c in curves]))
    pd.concat([pd.DataFrame({"rpt": i, "efc": e, "q": c["q"], "v": c["v"]})
               for i, (c, e) in enumerate(zip(curves, model_efc))]).to_csv(CSV, index=False)


def plot_rpts():
    if m.REAL_CELL_DIR_NAME:
        m.use_real_cell_rpts()
    elif T_ARG == "45":
        _use_cell017_rpts()
    df = pd.read_csv(CSV)
    curves = [dict(q=d["q"].to_numpy(), v=d["v"].to_numpy()) for _, d in df.groupby("rpt", sort=True)]
    model_efc = df.groupby("rpt", sort=True)["efc"].first().to_numpy()

    exp_cap = real("capacity_fade")
    c20 = exp_cap[exp_cap["source"] == "C/20 RPT"]
    real_efc = dict(zip(c20["rpt"].astype(int), c20["efc"]))
    rpts = [r for r in m.REAL_RPT_NUMS_WITH_DISCHARGE if r in real_efc]

    pairs, used = [], set()
    for r in rpts:
        order = np.argsort(np.abs(model_efc - real_efc[r]))
        j = next((int(k) for k in order if int(k) not in used), int(order[0]))
        used.add(j)
        q_r, v_r = m.load_real_rpt_discharge(r)
        pairs.append((r, real_efc[r], q_r / q_r.max(), v_r, float(model_efc[j]),
                      curves[j]["q"] / curves[j]["q"].max(), curves[j]["v"]))

    ncol = 2 if len(pairs) <= 4 else 4
    nrow = int(np.ceil(len(pairs) / ncol))
    title = (f"{CELL} ({T_LABEL} degC): SOC-normalised C/20 RPT {{}} -- "
             f"{TAG} ({os.environ.get('RPT_BASELINE_ENV', 'crack_baseline_v3.env').removesuffix('.env')}, width {os.environ['C64_WIDTH']}, F0 {os.environ['C64_F0']})")
    for kind in ("voltage", "dvdsoc"):
        fig, axes = plt.subplots(nrow, ncol, figsize=(5.5 * ncol, 4.5 * nrow), squeeze=False)
        for ax, (r, e_r, soc_r, v_r, e_m, soc_m, v_m) in zip(axes.flat, pairs):
            lab_r, lab_m = f"Real RPT{r} (EFC~{e_r:.0f})", f"Model RPT (EFC~{e_m:.0f})"
            if kind == "voltage":
                ax.plot(soc_r, v_r, "-", color="black", lw=2, label=lab_r)
                ax.plot(soc_m, v_m, "--", color="tab:blue", lw=1.8, label=lab_m)
                ax.set_ylabel("Terminal voltage [V]")
            else:
                # The model's first ~1% of DoD is the step-start numerical
                # transient (visible in V too); drop it before differentiating.
                keep = soc_m >= DVDSOC_TRIM
                ax.plot(*smooth_dvdsoc(soc_r, v_r), "-", color="black", lw=2, label=lab_r)
                ax.plot(*smooth_dvdsoc(soc_m[keep], v_m[keep]), "--", color="tab:blue", lw=1.8, label=lab_m)
                ax.set_ylabel("dV/d(SOC) [V]")
                ax.set_ylim(-3, 0.5)
            ax.set_xlabel("Depth of discharge [-]")
            ax.set_title(f"RPT{r}")
            ax.legend(fontsize=9)
            ax.grid(alpha=0.3)
        for ax in list(axes.flat)[len(pairs):]:
            ax.set_visible(False)
        fig.suptitle(title.format("discharge" if kind == "voltage" else "dV/d(SOC)"))
        fig.tight_layout()
        out = os.path.join(OUT, f"rpt_{'voltage' if kind == 'voltage' else 'dvdsoc'}_dod_{TAG}.png")
        fig.savefig(out, dpi=150)
        print("Saved:", out, flush=True)


if __name__ == "__main__":
    if "--replot" not in sys.argv or not os.path.exists(CSV):
        simulate()
    plot_rpts()
