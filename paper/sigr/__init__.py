"""sigr: reproduce the Si/graphite degradation-model results for the four
test conditions (CELL064 25 degC full SoC, CELL017 45 degC full SoC,
CELL009 25 degC 15-95% SoC, CELL026 25 degC 20-80% SoC) with the modified
PyBaMM in this repository.

    import sigr
    res, curves = sigr.run_condition("25deg_full_soc", data_dir, out_dir)
    sigr.plots.save_all(res, curves, sigr.CONDITIONS["25deg_full_soc"], data_dir, out_dir)
"""
from concurrent.futures import ProcessPoolExecutor

from . import analysis, data, plots, simulate  # noqa: F401
from .conditions import CONDITIONS, EFC_BASIS_AH, RATE_BASIS_AH  # noqa: F401
from .parameters import MODEL_OPTIONS, V6, parameter_values  # noqa: F401


def run_condition(key, data_dir, out_dir, rate_basis_ah=RATE_BASIS_AH, efc_basis_ah=EFC_BASIS_AH, log=print):
    """Simulate one condition (~1-1.5 h), save <key>_results.npz and
    <key>_rpt_curves.csv to out_dir, and return (results, rpt_curves)."""
    cond = CONDITIONS[key]
    sol = simulate.run(cond, data_dir, rate_basis_ah, efc_basis_ah, log=log)
    res, curves = simulate.reduce(sol, cond, efc_basis_ah)
    analysis.save(res, curves, key, out_dir)
    return res, curves


def _run_quiet(args):
    key, data_dir, out_dir, rate_basis_ah, efc_basis_ah = args
    run_condition(key, data_dir, out_dir, rate_basis_ah, efc_basis_ah, log=lambda *a: None)
    return key


def run_all(data_dir, out_dir, keys=None, rate_basis_ah=RATE_BASIS_AH, efc_basis_ah=EFC_BASIS_AH, parallel=True):
    """Simulate several conditions, in parallel processes by default."""
    keys = list(keys or CONDITIONS)
    jobs = [(k, data_dir, out_dir, rate_basis_ah, efc_basis_ah) for k in keys]
    if not parallel:
        return [_run_quiet(j) for j in jobs]
    with ProcessPoolExecutor(max_workers=len(jobs)) as pool:
        return list(pool.map(_run_quiet, jobs))
