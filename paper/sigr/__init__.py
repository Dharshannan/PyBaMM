"""sigr: reproduce the Si/graphite degradation-model results for the four
test conditions (CELL064 25 degC full SoC, CELL017 45 degC full SoC,
CELL009 25 degC 15-95% SoC, CELL026 25 degC 20-80% SoC) with the modified
PyBaMM in this repository.

    import sigr
    sigr.run_all(data_dir, out_dir)                    # all four, in parallel
    res, curves = sigr.analysis.load("25deg_full_soc", out_dir)
    sigr.plots.save_all(res, curves, sigr.CONDITIONS["25deg_full_soc"], data_dir, out_dir)
"""
import multiprocessing
import os
import time
from concurrent.futures import ProcessPoolExecutor, wait

from . import analysis, data, plots, simulate  # noqa: F401
from .conditions import CONDITIONS, EFC_BASIS_AH, RATE_BASIS_AH  # noqa: F401
from .parameters import MODEL_OPTIONS, V7, parameter_values  # noqa: F401


def run_condition(key, data_dir, out_dir, rate_basis_ah=RATE_BASIS_AH, efc_basis_ah=EFC_BASIS_AH, log=print,
                  options=None, params=None):
    """Simulate one condition (~20-40 min on one core), save <key>_results.npz and
    <key>_rpt_curves.csv to out_dir, and return (results, rpt_curves).
    options / params override MODEL_OPTIONS / V7 (see ablations.py)."""
    cond = CONDITIONS[key]
    res, curves = simulate.run(cond, data_dir, rate_basis_ah, efc_basis_ah, log=log, options=options, params=params)
    analysis.save(res, curves, key, out_dir)
    return res, curves


def _log_path(key, out_dir):
    return os.path.join(out_dir, f"{key}_log.txt")


def _run_logged(job):
    """One job in a worker process; progress goes to <out_dir>/<key>_log.txt."""
    os.makedirs(job["out_dir"], exist_ok=True)
    with open(_log_path(job["key"], job["out_dir"]), "w", encoding="utf-8") as fh:
        def log(*a):
            print(*a, file=fh, flush=True)
        t0 = time.time()
        run_condition(job["key"], job["data_dir"], job["out_dir"], job["rate_basis_ah"], job["efc_basis_ah"],
                      log=log, options=job.get("options"), params=job.get("params"))
        log(f"done in {(time.time() - t0) / 3600:.2f} h")
    return job["key"]


def _last_line(path):
    try:
        with open(path, encoding="utf-8") as fh:
            lines = [ln.strip() for ln in fh if ln.strip()]
        return lines[-1] if lines else "starting"
    except OSError:
        return "queued"


def _run_jobs(jobs, max_workers=None, report_every_s=300):
    """Run jobs in parallel worker processes, printing each one's latest
    progress line every report_every_s seconds."""
    workers = max_workers or max(1, min(len(jobs), os.cpu_count() or 1))
    print(f"Running {len(jobs)} job(s) on {workers} worker process(es); logs: <out_dir>/<condition>_log.txt",
          flush=True)
    t0 = time.time()
    # "spawn": fresh worker processes. Forking a kernel that already runs threads
    # (e.g. JAX on Colab) can deadlock.
    with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        pending = {pool.submit(_run_logged, j) for j in jobs}
        while pending:
            done, pending = wait(pending, timeout=report_every_s)
            for f in done:
                f.result()  # re-raises a failed run's exception
            print(f"--- {(time.time() - t0) / 60:.0f} min", flush=True)
            for j in jobs:
                label = j.get("label", j["key"])
                print(f"  {label:<40} {_last_line(_log_path(j['key'], j['out_dir']))[:100]}", flush=True)


def run_all(data_dir, out_dir, keys=None, rate_basis_ah=RATE_BASIS_AH, efc_basis_ah=EFC_BASIS_AH,
            max_workers=None, report_every_s=300):
    """Simulate several conditions in parallel worker processes (~20-40 min each
    with one core per condition). max_workers defaults to min(#conditions,
    #CPUs); each run needs ~1 GB of RAM."""
    keys = list(keys or CONDITIONS)
    _run_jobs([dict(key=k, data_dir=data_dir, out_dir=out_dir, rate_basis_ah=rate_basis_ah,
                    efc_basis_ah=efc_basis_ah) for k in keys], max_workers, report_every_s)
    return keys
