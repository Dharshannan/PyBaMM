"""task_si_diffusivity_sweep.py -- 2026-09-19 theory test (user's own
hypothesis): as Si LAM shrinks eps_s, the remaining active Si sees higher
effective current density; if the (never-fitted, no high-C-rate pouch data
available) Si solid diffusivity is too low, this over-estimates solid-
diffusion polarisation and artificially truncates Si's own per-cycle
stoichiometry excursion well before the real accessible capacity is used --
explaining the model's severe post-knee sto-range collapse (confirmed this
session: even at the model's own slow C/20 RPT rate, si_delta_sto collapses
9.3x between EFC~199 and EFC~220, right where real data barely moves).

Runs, SEQUENTIALLY (memory discipline -- repeated OOM history this
session), a baseline (mult=1, no residual) then a log-spaced Si diffusivity
multiplier sweep (3x/10x/30x/100x), each on top of the exact config used in
dma/extract_model_si_sto_window_c20rpt.py (GR_K=1.25e-5/BOL=0.130, SI
redirect+OCP-deform+volchange-deform item29 on). Then, on top of whichever
multiplier gives the best expansion-shape RMSE, additionally re-enables
item 30 (LAM expansion residual fraction, frac=0.9) to test the user's
follow-up request: does restoring isolated-Si's own expansion contribution
help further NOW that the dominant sto-collapse factor is (hopefully)
addressed by the diffusivity fix.

Each run's plots are tagged and kept (C64_NO_PLOT=0); results are collected
into si_diffusivity_sweep_results.csv (this directory) and mirrored to
../../dma/ for the consolidated overnight report.
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
WORKER = os.path.join(PHASE_DIR, "worker_diffusivity_probe.py")
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
    "C64_SI_LAM_EXPANSION_RESIDUAL", "C64_SI_LAM_EXPANSION_RESIDUAL_FRAC",
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
    env["C64_EFC_LIMIT"] = "260"
    env["C64_NO_PLOT"] = "0"
    return env


SUMMARY_RE = re.compile(
    r"\[SWEEP_SUMMARY\] tag=(\S+) expansion_rmse=(\S+) si_delta_first=(\S+) "
    r"si_delta_min=(\S+) collapse_ratio=(\S+) final_soh_pct=(\S+) rpt_gap_pp=(\S+)"
)


def run_one(tag, extra_env):
    env = base_env()
    env.update(extra_env)
    env["C64_OUT_TAG"] = tag
    logpath = os.path.join(LOG_DIR, f"{tag}.log")
    print(f"[SWEEP] launching {tag} ({extra_env}), log -> {logpath}", flush=True)
    with open(logpath, "w", encoding="utf-8") as fh:
        proc = subprocess.run(
            [PYEXE, WORKER], cwd=PHASE_DIR, env=env, stdout=fh, stderr=subprocess.STDOUT
        )
    print(f"[SWEEP] {tag} exit_code={proc.returncode}", flush=True)
    with open(logpath, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    m = SUMMARY_RE.search(text)
    if m is None:
        print(f"[SWEEP] {tag}: no summary line found, tail of log:", flush=True)
        print("\n".join(text.splitlines()[-40:]), flush=True)
        return {
            "tag": tag, "exit_code": proc.returncode, "expansion_rmse": None,
            "si_delta_first": None, "si_delta_min": None, "collapse_ratio": None,
            "final_soh_pct": None, "rpt_gap_pp": None,
        }
    tag_, rmse, si_first, si_min, ratio, soh, gap = m.groups()
    row = {
        "tag": tag,
        "exit_code": proc.returncode,
        "expansion_rmse": float(rmse) if rmse != "nan" else None,
        "si_delta_first": float(si_first),
        "si_delta_min": float(si_min),
        "collapse_ratio": float(ratio),
        "final_soh_pct": float(soh),
        "rpt_gap_pp": float(gap) if gap != "nan" else None,
    }
    print(f"[SWEEP] {tag}: {row}", flush=True)
    return row


def main():
    rows = []

    # Phase 1: baseline + Si diffusivity multiplier sweep (no residual term)
    diffusivity_mults = [1, 3, 10, 30, 100]
    for mult in diffusivity_mults:
        tag = f"sidiff_x{mult}"
        rows.append(run_one(tag, {"C64_SI_DIFFUSIVITY_MULT": str(mult)}))

    # Pick the best-so-far by expansion RMSE (lower is better), but only
    # among runs that didn't also blow up the ALREADY-TUNED capacity-fade
    # fit (per the user's own flagged concern: raising Si diffusivity could
    # shift post-knee C/3 capacity behaviour even if slow-rate RPTs move
    # much less) -- reject anything with |rpt_gap_pp| > 20 as a sanity cap,
    # historical gaps for accepted configs have been single-digit-to-low-
    # teens percentage points.
    valid = [
        r for r in rows
        if r["expansion_rmse"] is not None
        and r["rpt_gap_pp"] is not None
        and abs(r["rpt_gap_pp"]) < 20
    ]
    if valid:
        best = min(valid, key=lambda r: r["expansion_rmse"])
        best_mult = best["tag"].split("_x")[-1]
        print(f"[SWEEP] best diffusivity multiplier so far: {best_mult} "
              f"(expansion_rmse={best['expansion_rmse']:.4f}, "
              f"rpt_gap_pp={best['rpt_gap_pp']:.2f})", flush=True)
    else:
        best_mult = "10"
        print("[SWEEP] no valid (capacity-fade-safe) runs in phase 1, "
              "defaulting best_mult=10", flush=True)

    # Phase 2: on top of the best diffusivity multiplier, re-enable item 30
    # (LAM expansion residual fraction) at two fractions to test the user's
    # follow-up request.
    for frac in [0.5, 0.9]:
        tag = f"sidiff_x{best_mult}_residual{str(frac).replace('.', '')}"
        rows.append(
            run_one(
                tag,
                {
                    "C64_SI_DIFFUSIVITY_MULT": best_mult,
                    "C64_SI_LAM_EXPANSION_RESIDUAL": "1",
                    "C64_SI_LAM_EXPANSION_RESIDUAL_FRAC": str(frac),
                },
            )
        )

    # Also test residual alone (mult=1) for a clean before/after comparison.
    rows.append(
        run_one(
            "residual09_only",
            {
                "C64_SI_DIFFUSIVITY_MULT": "1",
                "C64_SI_LAM_EXPANSION_RESIDUAL": "1",
                "C64_SI_LAM_EXPANSION_RESIDUAL_FRAC": "0.9",
            },
        )
    )

    out_csv = os.path.join(PHASE_DIR, "si_diffusivity_sweep_results.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "tag", "exit_code", "expansion_rmse", "si_delta_first",
                "si_delta_min", "collapse_ratio", "final_soh_pct", "rpt_gap_pp",
            ],
        )
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    print(f"\n[SWEEP] all done, results -> {out_csv}", flush=True)
    shutil.copy(out_csv, os.path.join(DMA_DIR, "si_diffusivity_sweep_results.csv"))
    print(f"[SWEEP] mirrored -> {DMA_DIR}", flush=True)

    print("\n=== FINAL SUMMARY TABLE ===", flush=True)
    for r in rows:
        print(r, flush=True)


if __name__ == "__main__":
    main()
