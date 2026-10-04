# -*- coding: utf-8 -*-
"""study_si_deg: does Si crack GROWTH ever matter in the high_temp_45C recipe?

Finding that motivated this (2026-09-29): in both the study_si_deg baseline
and the high_temp_45C fork, Si crack length stays at its initial 2e-8 m and
roughness at ~2.908 (set by initial crack length x crack density), so every
crack-rate knob (SI_CRACK_MULT, SI_CRACK_EAC) is inert and "SEI on cracks"
only ever acts on the fixed initial crack area.

This runs the high_temp_45C fork UNMODIFIED apart from env overrides
(typically C64_T_AMBIENT_K=298.15 and C64_SI_CRACK_MULT), and reports, at
each RPT, how much the cracks have grown and how much crack SEI contributes
to porosity closure and LLI, versus bulk SEI.

Physical coupling being tested: cracking does not change the SEI growth law
per unit area, but it adds fresh surface. "SEI on cracks" grows on
a_cr = a * (roughness - 1) with its own thickness L_sei_cr, so crack growth
raises the SEI growth RATE per electrode volume in two ways: more area, and
a thinner film on the newly exposed area (diffusion-limited j ~ 1/L). Both
feed pore closure (reaction_driven_porosity.py's L_sei_cr * (roughness - 1)
term) and LLI.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "high_temp_45C"))
os.chdir(HERE)  # plots land in study_si_deg/, not the fork's folder

EXTRA = [
    "X-averaged negative secondary particle crack length [m]",
    "X-averaged negative secondary electrode roughness ratio",
    "X-averaged negative secondary particle surface tangential stress [Pa]",
    "X-averaged negative secondary SEI thickness [m]",
    "X-averaged negative secondary SEI on cracks thickness [m]",
    "Loss of lithium to negative secondary SEI [mol]",
    "Loss of lithium to negative secondary SEI on cracks [mol]",
    "Total lithium lost from particles [mol]",
    "Total capacity lost to side reactions [A.h]",
]
os.environ.setdefault("C64_DIAG_EXTRA_VARS", "|".join(EXTRA))

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import cell064_degradation_fit_45C as m  # noqa: E402

# _savefig() reads SCRIPT_DIR at call time; data paths were already resolved
# at import, so this only redirects the output PNGs into study_si_deg/.
m.SCRIPT_DIR = HERE

TAG = os.environ.get("C64_OUT_TAG", "crack_sweep")
print("Config: T_amb =", os.environ.get("C64_T_AMBIENT_K", "318.15 (fork default)"),
      "| SOLVENT_MULT =", m.SOLVENT_MULT, "| SI_CRACK_MULT =", m.SI_CRACK_RATE_MULT,
      "| SI_CRACK_EAC =", m.SI_CRACK_EAC, flush=True)

sol = m.run_degradation()
results = m.extract_results(sol)
# Grab RPT throughputs before plot_all(), which was observed to leave
# results["rpt_thr"] empty afterwards.
rpt_thr = np.array(results["rpt_thr"], dtype=float).copy()
m.plot_all(results, sol)

Qt = sol["Throughput capacity [A.h]"].entries
efc_full = m.efc_from_throughput(Qt)
l_cr = sol[EXTRA[0]].entries
rough = sol[EXTRA[1]].entries
sig_t = sol[EXTRA[2]].entries
L_bulk = sol[EXTRA[3]].entries
L_cr = sol[EXTRA[4]].entries
lli_bulk = sol[EXTRA[5]].entries
lli_cr = sol[EXTRA[6]].entries
cr_contrib = L_cr * (rough - 1.0)
F = 96485.33212
Q_nom = m.NOMINAL_CAP_AH

rpt_efc = m.efc_from_throughput(rpt_thr) if rpt_thr.size else np.arange(0.0, float(efc_full[-1]), 25.0)
print("\n--- Si crack growth and crack-vs-bulk SEI, at each RPT ---", flush=True)
print(f"{'EFC':>7} {'l_cr [nm]':>10} {'l/l0':>7} {'rough':>7} {'peak sig_t [MPa]':>17} "
      f"{'pore: cracks %':>15} {'LLI bulk %':>11} {'LLI cracks %':>13}", flush=True)
prev = 0
for e in rpt_efc:
    idx = int(np.argmin(np.abs(efc_full - e)))
    tot = L_bulk[idx] + cr_contrib[idx]
    peak = float(np.max(sig_t[prev:idx + 1])) if idx >= prev else float("nan")
    prev = idx
    print(f"{e:7.1f} {l_cr[idx] * 1e9:10.3f} {l_cr[idx] / l_cr[0]:7.3f} {rough[idx]:7.3f} "
          f"{peak / 1e6:17.2f} {100 * cr_contrib[idx] / tot:15.1f} "
          f"{100 * lli_bulk[idx] * F / 3600 / Q_nom:11.2f} {100 * lli_cr[idx] * F / 3600 / Q_nom:13.2f}",
          flush=True)

# LLI like-for-like vs DMA: the plotted "side reactions + LAM-trapped" LLI
# overstates the model's true lithium-inventory loss by ~2-3pp; the DMA's
# charge_LLI is an inventory quantity, so compare against particle Li lost.
lli_part = 100 * sol[EXTRA[7]].entries * F / 3600 / Q_nom
_T = float(os.environ.get("C64_T_AMBIENT_K", 318.15))
exp_lli = None
if abs(_T - 298.15) < 0.01:
    exp_lli = m.load_experimental_lli()
elif abs(_T - 318.15) < 0.01 and os.path.isdir(m.CELL017_DIR):
    exp_lli = m.load_cell017("LLI")
if exp_lli is not None:
    print("\n--- LLI check at each real DMA EFC: real / model particle-inventory LLI ---", flush=True)
    for e, r in zip(exp_lli["efc"], 100 * exp_lli["charge_LLI_loss"] / exp_lli["charge_LLI"].iloc[0]):
        if 0 < e <= efc_full[-1]:
            print(f"  EFC={e:.1f}: real={r:.1f}%, model particles={np.interp(e, efc_full, lli_part):.1f}%", flush=True)

# LAM check: model Gr/Si LAM at each real DMA RPT, plus end-of-run values.
exp_lam = None
if abs(_T - 298.15) < 0.01:
    exp_lam = m.load_experimental_lam()
elif abs(_T - 318.15) < 0.01 and os.path.isdir(m.CELL017_DIR):
    exp_lam = m.load_cell017("LAM_derived")
efc_lam = m.efc_from_throughput(results["Qt"])
if exp_lam is not None:
    print("\n--- LAM check at each real DMA EFC: real / model (Gr, Si) ---", flush=True)
    for _, r in exp_lam.iterrows():
        if 0 < r["efc"] <= efc_lam[-1]:
            print(f"  EFC={r['efc']:.1f}: Gr real={r['LAM_NE_graphite_pct']:.1f}% model="
                  f"{np.interp(r['efc'], efc_lam, results['LAM_gr']):.1f}% | Si real="
                  f"{r['LAM_NE_silicon_pct']:.1f}% model={np.interp(r['efc'], efc_lam, results['LAM_si']):.1f}%",
                  flush=True)
print(f"  end of run (EFC {efc_lam[-1]:.0f}): Gr LAM {results['LAM_gr'][-1]:.1f}%, "
      f"Si LAM {results['LAM_si'][-1]:.1f}%", flush=True)

fig, ax = plt.subplots(2, 2, figsize=(11, 7))
ax[0, 0].plot(efc_full, l_cr / l_cr[0])
ax[0, 0].set_ylabel("Si crack length / initial [-]")
ax[0, 1].plot(efc_full, sig_t / 1e6, lw=0.3)
ax[0, 1].set_ylabel("Si surface tangential stress [MPa]")
ax[1, 0].plot(efc_full, 100 * cr_contrib / (L_bulk + cr_contrib))
ax[1, 0].set_ylabel("cracks' share of pore-closing SEI [%]")
ax[1, 1].plot(efc_full, 100 * lli_bulk * F / 3600 / Q_nom, label="bulk Si SEI")
ax[1, 1].plot(efc_full, 100 * lli_cr * F / 3600 / Q_nom, label="Si SEI on cracks")
ax[1, 1].set_ylabel("LLI [% of nominal]")
ax[1, 1].legend()
for a in ax.flat:
    a.set_xlabel("EFC [-]")
    a.grid(alpha=0.3)
fig.suptitle(f"Si crack growth diagnostics ({TAG})")
fig.tight_layout()
out = os.path.join(HERE, f"crack_diagnostics_{TAG}.png")
fig.savefig(out, dpi=130)
print("Saved:", out, flush=True)
