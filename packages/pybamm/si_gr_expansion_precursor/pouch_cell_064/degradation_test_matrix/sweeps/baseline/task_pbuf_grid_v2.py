"""task_pbuf_grid_v2.py -- 2026-09-19 phase 5: pore-buffering (F0 x width)
re-sweep on top of the NEW best baseline (grid_bol0p17_ksi3em4: radius x0.2,
diffusivity x3, BOL=0.17, SI_K_SEI=3e-4/GR_K_SEI=6e-6, knee EFC=112.6). The
earlier-session F0 x width grid (pbuf_grid_results/, F0 in 0.5-0.7 x width
in 0.03-0.15) was tuned against the OLD, unfixed baseline -- the whole
radius/diffusivity/BOL/K_SEI retune has changed the underlying
electrochemistry (polarisation, knee timing, capacity fade) enough that the
pore-buffering fit likely needs re-doing from scratch, not assumed still
valid.

EFC_LIMIT=205 throughout (see [[cell064-efc-limit-205]]).

Grid: F0 (transmitted fraction plateau) in {0.5, 0.6, 0.7} x width
(transition width) in {0.03, 0.05, 0.10}.
"""
import csv
import os
import shutil
import subprocess
import sys

PHASE_DIR = os.path.dirname(os.path.abspath(__file__))
LOG_DIR = os.path.join(PHASE_DIR, "logs")
DMA_DIR = os.path.join(PHASE_DIR, "..", "..", "..", "dma")
PYEXE = r"C:\Users\ds3420\AppData\Local\anaconda3\envs\PyBaMM_Pressure\python.exe"
WORKER = os.path.join(PHASE_DIR, "worker_full_retune_probe.py")
os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(DMA_DIR, exist_ok=True)

sys.path.insert(0, PHASE_DIR)
from task_full_retune_sweep import (  # noqa: E402
    BASE_ENV_KEYS_TO_CLEAR, FIELDNAMES, SUMMARY_RE,
)


def base_env():
    env = dict(os.environ)
    for k in BASE_ENV_KEYS_TO_CLEAR:
        env.pop(k, None)
    env.pop("C64_F0", None)
    env.pop("C64_WIDTH", None)
    env["C64_EXPONENT_MAX_SEI"] = "70"
    env["C64_SI_BETA_LAM_SEI"] = "1e-7"
    env["C64_SI_BETA_LAM_ISO"] = "1.0"
    env["C64_SI_LAM_ISO_EXPONENT"] = "3.0"
    env["C64_SI_TAU_LAM_ISO"] = "2e8"
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
    # new best baseline (grid_bol0p17_ksi3em4)
    env["C64_SI_PARTICLE_RADIUS_MULT"] = "0.2"
    env["C64_SI_DIFFUSIVITY_MULT"] = "3"
    env["C64_NEG_POROSITY_BOL"] = "0.17"
    env["C64_SI_K_SEI_MULT"] = "3e-4"
    env["C64_GR_K_SEI_MULT"] = "6e-6"
    return env


def run_one(tag, extra_env):
    env = base_env()
    env.update(extra_env)
    env["C64_OUT_TAG"] = tag
    logpath = os.path.join(LOG_DIR, f"{tag}.log")
    print(f"[PBUF2] launching {tag} ({extra_env}), log -> {logpath}", flush=True)
    with open(logpath, "w", encoding="utf-8") as fh:
        proc = subprocess.run(
            [PYEXE, WORKER], cwd=PHASE_DIR, env=env, stdout=fh, stderr=subprocess.STDOUT
        )
    print(f"[PBUF2] {tag} exit_code={proc.returncode}", flush=True)
    with open(logpath, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    m = SUMMARY_RE.search(text)
    if m is None:
        print(f"[PBUF2] {tag}: no summary line found, tail of log:", flush=True)
        print("\n".join(text.splitlines()[-40:]), flush=True)
        return {k: (proc.returncode if k == "exit_code" else (tag if k == "tag" else None)) for k in FIELDNAMES}
    vals = m.groups()
    row = {"tag": tag, "exit_code": proc.returncode}
    for k, v in zip(FIELDNAMES[2:], vals[1:]):
        row[k] = float(v) if v not in ("nan", "None") else None
    print(f"[PBUF2] {tag}: {row}", flush=True)
    return row


def composite_score(row):
    if row is None or row.get("expansion_rmse") is None:
        return float("inf")
    parts = [
        row["expansion_rmse"] * 50,
        abs(row.get("knee_err_efc") or 100),
        abs(row.get("rpt_gap_pp") or 50),
    ]
    if any(p != p for p in parts):
        return float("inf")
    return sum(parts)


def main():
    rows = []
    f0_values = ["0.5", "0.6", "0.7"]
    width_values = ["0.03", "0.05", "0.10"]

    for f0 in f0_values:
        for width in width_values:
            tag = f"pbuf2_f{f0.replace('.', 'p')}_w{width.replace('.', 'p')}"
            row = run_one(tag, {"C64_F0": f0, "C64_WIDTH": width})
            row["_f0"] = f0
            row["_width"] = width
            rows.append(row)

    best = min(rows, key=composite_score)
    print(f"\n[PBUF2] best overall: {best}", flush=True)

    out_csv = os.path.join(PHASE_DIR, "pbuf_grid_v2_results.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k) for k in FIELDNAMES})
    print(f"[PBUF2] all done, results -> {out_csv}", flush=True)
    shutil.copy(out_csv, os.path.join(DMA_DIR, "pbuf_grid_v2_results.csv"))
    print(f"[PBUF2] mirrored -> {DMA_DIR}", flush=True)

    print("\n=== FINAL SUMMARY TABLE ===", flush=True)
    for r in rows:
        print(r, flush=True)


if __name__ == "__main__":
    main()
