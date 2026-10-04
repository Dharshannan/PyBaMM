# -*- coding: utf-8 -*-
"""k_compliance_test: does the MEASURED k, applied to the model's own particle
swelling, reproduce the CELL017 (45 degC) reversible-expansion curve?

Option 2 of the 2026-09-30 k discussion: a post-hoc diagnostic with no
feedback. Physics, SoH, LAM and LLI are exactly crack_baseline_v2 (R3c).
Only the cell thickness is rebuilt afterwards.

In the model (reaction_driven_porosity.py) the negative electrode's swelling
is partitioned by the compliance fraction f (= model k):
    dv_buffered  = min((1 - f) dv_solid, headroom)   -> fills pores
    neg_thick    = n L (dv_solid - dv_buffered)       -> reaches thickness
Here the reported thickness is recomputed with the measured k instead,
leaving the rest of the cell (positive electrode, thermal) untouched:
    cell_k = cell - n L (dv_solid - dv_buffered) + k_exp(EFC) n L dv_solid
k_exp(EFC) is CELL017_k_expansion_scale.csv, linearly interpolated between
RPTs and held flat outside them. Variant B instead applies k_exp to the
whole electrode stack (negative particle swelling + positive electrode),
since the positive electrode is ~40% of the model's cell breathing and the
measured rev_particle basis may be the full stack. Expansion is then scored exactly as in the
fork (per-cycle peak-to-peak, normalised to the first-5-cycle median).

If cell_k matches the data, the whole expansion gap is k, and a proper
"prescribed thickness transfer" opt-in (option 3) is worth building. If it
doesn't, the model's particle-level swelling (LAM split / stoichiometry) is
also off. Panel (c) shows that directly as real rev_cell / k_exp against the
model's particle swelling.

Usage (from this folder):  python k_compliance_test.py [--reuse]
--reuse skips the simulation and replots from k_compliance_test_cache.npz.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TAG = "k_compliance_test"
CACHE = os.path.join(HERE, f"{TAG}_cache.npz")

# crack_baseline_v2 (R3c) at 45 degC, CELL017 F0. Explicit env vars win.
with open(os.path.join(HERE, "crack_baseline_v2.env")) as fh:
    for line in fh:
        line = line.split("#", 1)[0].strip()
        if "=" in line:
            key, val = line.split("=", 1)
            os.environ.setdefault(key.strip(), val.strip())
os.environ.setdefault("C64_T_AMBIENT_K", "318.15")
os.environ.setdefault("C64_F0", "0.8")
os.environ.setdefault("C64_OUT_TAG", TAG)

EXTRA = [
    "Negative electrode solid volume change",
    "Negative electrode buffered volume change",
    "Negative electrode thickness change [m]",
    "Positive electrode thickness change [m]",
]
os.environ.setdefault("C64_DIAG_EXTRA_VARS", "|".join(EXTRA))

sys.path.insert(0, os.path.join(os.path.dirname(HERE), "high_temp_45C"))
os.chdir(HERE)

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import cell064_degradation_fit_45C as m  # noqa: E402

m.SCRIPT_DIR = HERE
N_REF = 5  # same early-life normalisation window as the fork


def _scalar_series(cyc, name):
    e = np.asarray(cyc[name].entries, dtype=float)
    return e.mean(axis=0) if e.ndim > 1 else e


def simulate():
    """Run R3c at 45 degC and keep, per C/3 ageing cycle, the time series
    needed to rebuild the thickness with any k."""
    sol = m.run_degradation()
    results = m.extract_results(sol)
    os.environ["C64_NO_PLOT"] = "1"
    m.plot_all(results, sol)  # scores only; confirms this run matches R3c

    pv = m.build_parameter_values()
    nL = float(pv["Number of electrodes connected in parallel to make a cell"]) * float(
        pv["Negative electrode thickness [m]"])

    thr, cell, dvs, dvb, neg, pos, cyc_k = [], [], [], [], [], [], []
    for cyc in sol.cycles:
        cap, _, d_end, _ = m.cycle_ageing_leg(cyc)
        if cap <= 0 or d_end is None:
            continue
        thr.append(float(cyc["Throughput capacity [A.h]"].entries[-1]))
        cell.append(_scalar_series(cyc, "Cell thickness change [m]"))
        dvs.append(_scalar_series(cyc, EXTRA[0]))
        dvb.append(_scalar_series(cyc, EXTRA[1]))
        neg.append(_scalar_series(cyc, EXTRA[2]))
        pos.append(_scalar_series(cyc, EXTRA[3]))
        cyc_k.append(float(np.max(np.abs(_scalar_series(cyc, "Negative electrode transfer ratio k")))))

    # Sanity check: n L (dv_solid - dv_buffered) must equal the model's own
    # negative thickness change, or the reconstruction below is wrong.
    err = max(float(np.max(np.abs(nL * (s - b) - n))) for s, b, n in zip(dvs, dvb, neg))
    scale = max(float(np.max(np.abs(n))) for n in neg)
    print(f"Sanity: max |nL(dv_solid - dv_buffered) - neg thickness| = {err:.3e} m "
          f"(neg thickness scale {scale:.3e} m)", flush=True)

    lengths = np.array([len(c) for c in cell])
    cat = lambda xs: np.concatenate(xs)  # noqa: E731
    np.savez_compressed(
        CACHE, nL=nL, thr=np.array(thr), lengths=lengths, cell=cat(cell), dvs=cat(dvs),
        dvb=cat(dvb), neg=cat(neg), pos=cat(pos), cyc_k=np.array(cyc_k),
        rpt_thr=np.asarray(results["rpt_thr"], dtype=float),
        rpt_k=np.asarray(results["rpt_k_arr"], dtype=float),
    )
    print("Saved cache:", CACHE, flush=True)


def _split(arr, lengths):
    return np.split(arr, np.cumsum(lengths)[:-1])


def _norm(x):
    return x / np.median(x[:N_REF])


def _rmse(model_efc, model_norm, exp_efc, exp_norm):
    ok = (exp_efc >= model_efc[0]) & (exp_efc <= model_efc[-1])
    return float(np.sqrt(np.mean((np.interp(exp_efc[ok], model_efc, model_norm) - exp_norm[ok]) ** 2)))


def analyse():
    c = np.load(CACHE)
    nL, lengths = float(c["nL"]), c["lengths"]
    cell, dvs, dvb, neg, pos = (_split(c[k], lengths) for k in ("cell", "dvs", "dvb", "neg", "pos"))
    efc = m.efc_from_throughput(c["thr"])

    exp_k = m.load_cell017("k_expansion_scale")
    exp_x = m.load_cell017("reversible_expansion")
    k_efc, k_val = exp_k["efc"].to_numpy(), exp_k["k"].to_numpy()
    x_efc, x_val = exp_x["efc"].to_numpy(), exp_x["reversible_expansion_um"].to_numpy()
    k_at = lambda e: np.interp(e, k_efc, k_val)  # noqa: E731  flat outside the RPTs

    ptp = lambda s: float(s.max() - s.min()) * 1e6  # noqa: E731  um
    own = np.array([ptp(cl) for cl in cell])
    presc = np.array([ptp(cl - nL * (s - b) + k_at(e) * nL * s)
                      for cl, s, b, e in zip(cell, dvs, dvb, efc)])
    # Variant B: k_exp scales the whole electrode stack (negative particle
    # swelling + positive electrode), in case the measured rev_particle
    # basis is the full stack rather than the negative electrode alone.
    presc_all = np.array([ptp((cl - n - p) + k_at(e) * (nL * s + p))
                          for cl, n, p, s, e in zip(cell, neg, pos, dvs, efc)])
    particle = np.array([ptp(nL * s) for s in dvs])  # model particle-level (unbuffered) swelling
    neg_own = np.array([ptp(nL * (s - b)) for s, b in zip(dvs, dvb)])
    pos_own = np.array([ptp(p) for p in pos])

    x_norm = _norm(x_val)
    r_own = _rmse(efc, _norm(own), x_efc, x_norm)
    r_presc = _rmse(efc, _norm(presc), x_efc, x_norm)
    r_all = _rmse(efc, _norm(presc_all), x_efc, x_norm)
    # Real particle-level swelling implied by the measurement: rev_cell / k.
    real_particle = x_val / k_at(x_efc)
    r_part = _rmse(efc, _norm(particle), x_efc, _norm(real_particle))

    print("\n--- k compliance test (45 degC, crack_baseline_v2 / R3c) ---", flush=True)
    print(f"Expansion shape RMSE, model's own k      : {r_own:.4f}", flush=True)
    print(f"Expansion shape RMSE, measured k on negative only (A): {r_presc:.4f}", flush=True)
    print(f"Expansion shape RMSE, measured k on neg + pos (B)   : {r_all:.4f}", flush=True)
    print(f"Particle swelling RMSE, model vs real rev_cell/k: {r_part:.4f}", flush=True)
    print(f"Positive electrode share of cell breathing (first {N_REF} cycles): "
          f"{100 * np.median(pos_own[:N_REF]) / np.median(own[:N_REF]):.1f}%", flush=True)
    print("At each real RPT (normalised expansion: real / own k / A / B):", flush=True)
    for e in k_efc:
        if efc[0] <= e <= efc[-1]:
            print(f"  EFC={e:6.1f}: real={np.interp(e, x_efc, x_norm):.3f}  own={np.interp(e, efc, _norm(own)):.3f}"
                  f"  A={np.interp(e, efc, _norm(presc)):.3f}  B={np.interp(e, efc, _norm(presc_all)):.3f}  (k_exp={k_at(e):.3f})", flush=True)

    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    a = ax[0, 0]
    a.scatter(x_efc, x_norm, s=8, alpha=0.4, color="black", label="CELL017 real")
    a.plot(efc, _norm(own), color="tab:purple", label=f"Model, own compliance k (RMSE {r_own:.3f})")
    a.plot(efc, _norm(presc), color="tab:green", label=f"A: measured k on negative (RMSE {r_presc:.3f})")
    a.plot(efc, _norm(presc_all), color="tab:olive", ls="--", label=f"B: measured k on neg + pos (RMSE {r_all:.3f})")
    a.set_ylabel("Reversible expansion, normalised [-]")
    a.set_title("(a) Cell reversible expansion")

    a = ax[0, 1]
    a.scatter(k_efc, k_val, s=60, color="black", zorder=5, label="CELL017 real k (RPTs)")
    ee = np.linspace(efc[0], efc[-1], 400)
    a.plot(ee, k_at(ee), "k--", lw=1, label="k_exp interpolant (used post-hoc)")
    a.plot(efc, c["cyc_k"], color="tab:purple", lw=1, label="Model compliance k (per C/3 cycle peak)")
    a.set_ylabel("k [-]")
    a.set_title("(b) Expansion scale k")

    a = ax[1, 0]
    a.scatter(x_efc, _norm(real_particle), s=8, alpha=0.4, color="black", label="Real rev_cell / k_exp")
    a.plot(efc, _norm(particle), color="tab:red", label=f"Model particle swelling (RMSE {r_part:.3f})")
    a.set_ylabel("Particle-level swelling, normalised [-]")
    a.set_title("(c) Particle-level swelling (what k multiplies)")

    a = ax[1, 1]
    a.plot(efc, own, color="tab:purple", label="Cell, own k")
    a.plot(efc, presc, color="tab:green", label="Cell, measured k (A)")
    a.plot(efc, presc_all, color="tab:olive", ls="--", label="Cell, measured k (B)")
    a.plot(efc, neg_own, color="tab:blue", lw=1, label="Negative electrode, own k")
    a.plot(efc, particle, color="tab:red", lw=1, label="Negative particle swelling (k = 1)")
    a.plot(efc, pos_own, color="tab:orange", lw=1, label="Positive electrode")
    a.set_ylabel("Per-cycle peak-to-peak [um]")
    a.set_title("(d) Model breathing contributions (absolute)")

    for a in ax.flat:
        a.set_xlabel("EFC since RPT1 [-]")
        a.grid(alpha=0.3)
        a.legend(fontsize=8)
    fig.suptitle("k compliance test, 45 degC: measured k applied post-hoc to crack_baseline_v2 (R3c); "
                 "degradation unchanged")
    fig.tight_layout()
    out = os.path.join(HERE, f"{TAG}_45degC.png")
    fig.savefig(out, dpi=130)
    print("Saved:", out, flush=True)


if __name__ == "__main__":
    if "--reuse" not in sys.argv or not os.path.exists(CACHE):
        simulate()
    analyse()
