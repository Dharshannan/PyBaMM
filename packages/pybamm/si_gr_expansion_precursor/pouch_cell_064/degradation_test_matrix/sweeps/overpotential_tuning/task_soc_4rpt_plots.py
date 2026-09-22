"""task_soc_4rpt_plots.py -- 2026-09-22, generates RPT1/2/4/5 curve CSVs
(via worker_overpotential_probe.py's dump_rpt45_curves, now covering all
four real-discharge RPTs, not just 4/5) for the three configs needed for
the SOC-normalised 4-subplot comparison plots:
    A. true_old_baseline    -- the pre-session canonical config (zero C64_*
                                overrides except the NEG_POROSITY_FLOOR=
                                0.028/MAX_CYCLES=350 stability fix -- see
                                dma/build_true_old_baseline.py), i.e.
                                BEFORE the LAM/LLI-split retune.
    B. new_baseline_no_ocp  -- current best baseline (pbuf2_f0p7_w0p03),
                                SI_OCP_AGING_DEFORM=0.
    C. new_baseline_ocp     -- same, SI_OCP_AGING_DEFORM=1, throughput
                                driver (best candidate from the OCP-
                                deformation investigation).
All at C/20 (no C64_RPT_RATE override).
"""
import os
import subprocess

PHASE_DIR = os.path.dirname(os.path.abspath(__file__))
LOG_DIR = os.path.join(PHASE_DIR, "logs")
PYEXE = r"C:\Users\ds3420\AppData\Local\anaconda3\envs\PyBaMM_Pressure\python.exe"
WORKER = os.path.join(PHASE_DIR, "worker_overpotential_probe.py")
os.makedirs(LOG_DIR, exist_ok=True)


def clear_env():
    env = dict(os.environ)
    for k in list(env):
        if k.startswith("C64_"):
            del env[k]
    return env


def new_baseline_env():
    env = clear_env()
    env["C64_SI_REDIRECT_TO_LAM"] = "1"
    env["C64_SI_BETA_LAM_ISO"] = "1.0"
    env["C64_SI_TAU_LAM_ISO"] = "2e8"
    env["C64_SI_REDIRECT_LAM_YIELD"] = "0.02"
    env["C64_SI_LAM_ISO_EXPONENT"] = "3.0"
    env["C64_SI_VOLUME_CHANGE_AGING_DEFORM"] = "1"
    env["C64_NEG_POROSITY_FLOOR"] = "0.035"
    env["C64_EXPONENT_MAX_SEI"] = "70"
    env["C64_SI_BETA_LAM_SEI"] = "1e-7"
    env["C64_GR_REDIRECT_TO_LAM"] = "0"
    env["C64_GR_REDIRECT_LAM_YIELD"] = "0.0"
    env["C64_MAX_CYCLES"] = "700"
    env["C64_SOH_FLOOR"] = "30"
    env["C64_EFC_LIMIT"] = "205"
    env["C64_NO_PLOT"] = "0"
    env["C64_SI_PARTICLE_RADIUS_MULT"] = "0.2"
    env["C64_SI_DIFFUSIVITY_MULT"] = "3"
    env["C64_NEG_POROSITY_BOL"] = "0.17"
    env["C64_SI_K_SEI_MULT"] = "3e-4"
    env["C64_GR_K_SEI_MULT"] = "6e-6"
    env["C64_F0"] = "0.7"
    env["C64_WIDTH"] = "0.03"
    return env


def old_baseline_env():
    env = clear_env()
    env["C64_NEG_POROSITY_FLOOR"] = "0.028"
    env["C64_MAX_CYCLES"] = "350"
    env["C64_NO_PLOT"] = "0"
    return env


def run_one(tag, env):
    env = dict(env)
    env["C64_OUT_TAG"] = tag
    logpath = os.path.join(LOG_DIR, f"{tag}.log")
    print(f"[SOC4RPT] launching {tag}, log -> {logpath}", flush=True)
    with open(logpath, "w", encoding="utf-8") as fh:
        proc = subprocess.run(
            [PYEXE, WORKER], cwd=PHASE_DIR, env=env, stdout=fh, stderr=subprocess.STDOUT
        )
    print(f"[SOC4RPT] {tag} exit_code={proc.returncode}", flush=True)
    with open(logpath, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    print("\n".join(text.splitlines()[-10:]), flush=True)


def main():
    run_one("true_old_baseline_4rpt", old_baseline_env())

    new_no_ocp = new_baseline_env()
    new_no_ocp["C64_SI_OCP_AGING_DEFORM"] = "0"
    run_one("new_baseline_no_ocp", new_no_ocp)

    new_ocp = new_baseline_env()
    new_ocp["C64_SI_OCP_AGING_DEFORM"] = "1"
    new_ocp["C64_OCP_DEFORM_DRIVER"] = "throughput"
    run_one("new_baseline_ocp", new_ocp)


if __name__ == "__main__":
    main()
