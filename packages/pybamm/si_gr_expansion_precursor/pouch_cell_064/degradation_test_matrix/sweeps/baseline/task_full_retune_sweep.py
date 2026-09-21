"""task_full_retune_sweep.py -- 2026-09-19 phase 2: with the stoichiometric-
collapse root cause addressed (new baseline: Si particle radius x0.2 +
diffusivity x3 -- exchange-current added negligible further benefit and is
left at x1 to keep this phase simpler), re-tune the SEI/porosity/LAM-
isolation parameters that likely drifted out of fit once the underlying
electrochemistry changed: knee timing, capacity-fade RPT gap, LLI, and the
LAM_Gr/LAM_Si split.

EFC_LIMIT=205 (not 260) per explicit compute-budget instruction -- covers
real RPT4/RPT5's range (171/197 in the shifted convention) without burning
extra cycles past it.

Phase A: NEG_POROSITY_BOL (initial BoL porosity -- sets how much pore
    volume exists before closure even starts) over 0.10/0.13/0.16/0.20.
Phase B: SI_K_SEI_MULT/GR_K_SEI_MULT (paired) around the current baseline,
    on top of the best BOL from phase A.
Phase C: SI_BETA_LAM_ISO (LAM-isolation rate, drives LAM_Si specifically)
    around baseline 1.0, on top of the best of phases A+B.

Sequential (memory discipline), each run tagged + full plot suite kept.
Uses worker_full_retune_probe.py (knee/RPT-gap/LLI/LAM composite scorer).
"""
import csv
import os
import re
import shutil
import subprocess

PHASE_DIR = os.path.dirname(os.path.abspath(__file__))
LOG_DIR = os.path.join(PHASE_DIR, "logs")
DMA_DIR = os.path.join(PHASE_DIR, "..", "..", "..", "dma")
PYEXE = r"C:\Users\ds3420\AppData\Local\anaconda3\envs\PyBaMM_Pressure\python.exe"
WORKER = os.path.join(PHASE_DIR, "worker_full_retune_probe.py")
os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(DMA_DIR, exist_ok=True)

BASE_ENV_KEYS_TO_CLEAR = [
    "C64_K_SEI_MULT", "C64_GR_K_SEI_MULT", "C64_SI_K_SEI_MULT",
    "C64_TIMESCALE_MULT", "C64_LAM_PROP_MULT", "C64_SI_LAM_PROP",
    "C64_SI_LAM_OPTION", "C64_SI_LAM_EMPIRICAL", "C64_EXPONENT_MAX_SEI",
    "C64_SI_BETA_LAM_SEI", "C64_NEG_POROSITY_FLOOR", "C64_SOH_FLOOR",
    "C64_SI_REDIRECT_TO_LAM", "C64_SI_REDIRECT_LAM_YIELD",
    "C64_GR_REDIRECT_TO_LAM", "C64_GR_REDIRECT_LAM_YIELD",
    "C64_F0", "C64_WIDTH", "C64_THERMAL_OPTION",
    "C64_SI_VOLUME_CHANGE_AGING_DEFORM", "C64_SI_DIFFUSIVITY_MULT",
    "C64_SI_EXCHANGE_CURRENT_MULT", "C64_SI_PARTICLE_RADIUS_MULT",
    "C64_SI_LAM_EXPANSION_RESIDUAL", "C64_SI_LAM_EXPANSION_RESIDUAL_FRAC",
    "C64_NEG_POROSITY_BOL", "C64_SI_BETA_LAM_ISO",
    "C64_EFC_LIMIT",
]


def base_env():
    env = dict(os.environ)
    for k in BASE_ENV_KEYS_TO_CLEAR:
        env.pop(k, None)
    env["C64_EXPONENT_MAX_SEI"] = "70"
    env["C64_SI_BETA_LAM_SEI"] = "1e-7"
    env["C64_SI_BETA_LAM_ISO"] = "1.0"
    env["C64_SI_LAM_ISO_EXPONENT"] = "3.0"
    env["C64_SI_TAU_LAM_ISO"] = "2e8"
    env["C64_NEG_POROSITY_BOL"] = "0.130"
    env["C64_SI_K_SEI_MULT"] = "6e-4"
    env["C64_GR_K_SEI_MULT"] = "1.25e-5"
    env["C64_SI_OCP_AGING_DEFORM"] = "1"
    env["C64_SI_VOLUME_CHANGE_AGING_DEFORM"] = "1"
    env["C64_SI_REDIRECT_TO_LAM"] = "1"
    env["C64_SI_REDIRECT_LAM_YIELD"] = "0.02"
    env["C64_GR_REDIRECT_TO_LAM"] = "0"
    env["C64_GR_REDIRECT_LAM_YIELD"] = "0.0"
    env["C64_NEG_POROSITY_FLOOR"] = "0.035"
    env["C64_MAX_CYCLES"] = "700"
    env["C64_SOH_FLOOR"] = "30"
    env["C64_EFC_LIMIT"] = "205"
    env["C64_NO_PLOT"] = "0"
    # new baseline from the radius/diffusivity investigation
    env["C64_SI_PARTICLE_RADIUS_MULT"] = "0.2"
    env["C64_SI_DIFFUSIVITY_MULT"] = "3"
    return env


SUMMARY_RE = re.compile(
    r"\[RETUNE2_SUMMARY\] tag=(\S+) expansion_rmse=(\S+) c3_collapse_ratio=(\S+) "
    r"rpt_collapse_ratio=(\S+) knee_err_efc=(\S+) rpt_gap_pp=(\S+) "
    r"lli_err_ah=(\S+) lam_gr_err_pp=(\S+) lam_si_err_pp=(\S+) final_soh_pct=(\S+)"
)

FIELDNAMES = [
    "tag", "exit_code", "expansion_rmse", "c3_collapse_ratio",
    "rpt_collapse_ratio", "knee_err_efc", "rpt_gap_pp", "lli_err_ah",
    "lam_gr_err_pp", "lam_si_err_pp", "final_soh_pct",
]


def run_one(tag, extra_env):
    env = base_env()
    env.update(extra_env)
    env["C64_OUT_TAG"] = tag
    logpath = os.path.join(LOG_DIR, f"{tag}.log")
    print(f"[RETUNE2] launching {tag} ({extra_env}), log -> {logpath}", flush=True)
    with open(logpath, "w", encoding="utf-8") as fh:
        proc = subprocess.run(
            [PYEXE, WORKER], cwd=PHASE_DIR, env=env, stdout=fh, stderr=subprocess.STDOUT
        )
    print(f"[RETUNE2] {tag} exit_code={proc.returncode}", flush=True)
    with open(logpath, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    m = SUMMARY_RE.search(text)
    if m is None:
        print(f"[RETUNE2] {tag}: no summary line found, tail of log:", flush=True)
        print("\n".join(text.splitlines()[-40:]), flush=True)
        return {k: (proc.returncode if k == "exit_code" else (tag if k == "tag" else None)) for k in FIELDNAMES}
    vals = m.groups()
    row = {"tag": tag, "exit_code": proc.returncode}
    for k, v in zip(FIELDNAMES[2:], vals[1:]):
        row[k] = float(v) if v not in ("nan", "None") else None
    print(f"[RETUNE2] {tag}: {row}", flush=True)
    return row


def composite_score(row):
    """Lower = better. Combines the four phase-2 targets on comparable
    scales (knee EFC error /10, capacity-fade pp, LLI Ah *20 to bring to a
    similar pp-like scale, LAM pp errors) -- a rough single number to rank
    candidates by, not a physically meaningful unit."""
    if row is None:
        return float("inf")
    parts = [
        abs(row.get("knee_err_efc") or 100) / 10.0,
        abs(row.get("rpt_gap_pp") or 50),
        abs(row.get("lli_err_ah") or 2.0) * 20.0,
        abs(row.get("lam_gr_err_pp") or 50),
        abs(row.get("lam_si_err_pp") or 50),
    ]
    if any(p != p for p in parts):  # NaN check
        return float("inf")
    return sum(parts)


def main():
    rows = []

    # Phase A: BoL porosity
    bol_values = ["0.10", "0.13", "0.16", "0.20"]
    phaseA = []
    for bol in bol_values:
        tag = f"retune2_bol{bol.replace('.', 'p')}"
        row = run_one(tag, {"C64_NEG_POROSITY_BOL": bol})
        row["_bol"] = bol
        phaseA.append(row)
        rows.append(row)
    best_bol = min(phaseA, key=composite_score)
    best_bol_val = best_bol["_bol"]
    print(f"[RETUNE2] phase A best: {best_bol}", flush=True)

    # Phase B: SEI kinetics (paired Si/Gr), on top of best BOL
    k_sei_pairs = [
        ("3e-4", "6e-6"), ("6e-4", "1.25e-5"), ("1.2e-3", "2.5e-5"),
    ]
    phaseB = []
    for si_k, gr_k in k_sei_pairs:
        tag = f"retune2_ksi{si_k.replace('.', 'p').replace('-', 'm')}"
        row = run_one(tag, {
            "C64_NEG_POROSITY_BOL": best_bol_val,
            "C64_SI_K_SEI_MULT": si_k, "C64_GR_K_SEI_MULT": gr_k,
        })
        row["_bol"] = best_bol_val
        row["_si_k"] = si_k
        row["_gr_k"] = gr_k
        phaseB.append(row)
        rows.append(row)
    best_kb = min(phaseB, key=composite_score)
    best_si_k = best_kb.get("_si_k", "6e-4")
    best_gr_k = best_kb.get("_gr_k", "1.25e-5")
    print(f"[RETUNE2] phase B best: {best_kb}", flush=True)

    # Phase C: LAM-isolation rate, on top of best BOL + K_SEI
    beta_values = ["0.5", "1.0", "2.0"]
    phaseC = []
    for beta in beta_values:
        tag = f"retune2_betaiso{beta.replace('.', 'p')}"
        row = run_one(tag, {
            "C64_NEG_POROSITY_BOL": best_bol_val,
            "C64_SI_K_SEI_MULT": best_si_k, "C64_GR_K_SEI_MULT": best_gr_k,
            "C64_SI_BETA_LAM_ISO": beta,
        })
        row["_bol"] = best_bol_val
        row["_si_k"] = best_si_k
        row["_gr_k"] = best_gr_k
        row["_beta"] = beta
        phaseC.append(row)
        rows.append(row)
    best_final = min(phaseC, key=composite_score)
    print(f"[RETUNE2] phase C (final) best: {best_final}", flush=True)

    out_csv = os.path.join(PHASE_DIR, "full_retune_sweep_results.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k) for k in FIELDNAMES})
    print(f"\n[RETUNE2] all done, results -> {out_csv}", flush=True)
    shutil.copy(out_csv, os.path.join(DMA_DIR, "full_retune_sweep_results.csv"))
    print(f"[RETUNE2] mirrored -> {DMA_DIR}", flush=True)

    print("\n=== FINAL SUMMARY TABLE ===", flush=True)
    for r in rows:
        print(r, flush=True)
    print(f"\nBest overall config: BOL={best_bol_val}, SI_K={best_si_k}, "
          f"GR_K={best_gr_k}, BETA_ISO={best_final.get('_beta')}", flush=True)


if __name__ == "__main__":
    main()
