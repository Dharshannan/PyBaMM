"""Run one condition: formation, then C/3 ageing in chunks with a C/20 RPT
placed at each real RPT's EFC (aligned RPTs), then reduce the solution to
compact results.

Protocol:
- Formation (from 100% SoC): C/20 discharge to 2.5 V, C/3 CC-CV charge to
  4.195 V (CV to C/50). Its discharge is the model's RPT at EFC 0.
- Partial windows only: a conditioned C/20 reference RPT straight after
  formation (the SoH reference, like the real RPT1). Before every later RPT,
  a 10 min rest and a CC-CV charge to 4.195 V (CC dropped when the window
  top is within 0.1 V, as for CELL009).
- Ageing: the condition's cycle (conditions.ageing_cycle), in chunks of at
  most 50 cycles. Each chunk is sized from the last cycle's EFC so the RPT
  discharge starts within +-0.5 ageing cycles of the real RPT's EFC.
- On a solver failure, the chunk is retried at tolerances x10, x0.1, x100,
  x0.01 before the run stops.
"""
import numpy as np
import pybamm

from . import data
from .conditions import (LOWER_CUTOFF_V, MODEL_CAPACITY_AH, RATE_BASIS_AH, EFC_BASIS_AH, RPT_RATE,
                         UPPER_CUTOFF_V, ageing_cycle)
from .parameters import MODEL_OPTIONS, VAR_PTS, parameter_values

SOLVER_TOL = 1e-6
RETRY_TOL_MULTIPLIERS = [10.0, 0.1, 100.0, 0.01]
BATCH_SIZE = 50
MAX_TOTAL_CYCLES = 1800
SOH_FLOOR_PERCENT = 30.0
F = 96485.33212

CRACK_VARS = [
    "X-averaged negative secondary particle crack length [m]",
    "X-averaged negative secondary electrode roughness ratio",
    "X-averaged negative secondary particle surface tangential stress [Pa]",
    "X-averaged negative electrode porosity",
    "Loss of lithium to negative secondary SEI [mol]",
    "Loss of lithium to negative secondary SEI on cracks [mol]",
]
STRAIN_VARS = [
    "X-averaged negative secondary particle strain-fatigue driver [s-1]",
    "X-averaged negative secondary particle strain-fatigue cracking rate [m.s-1]",
    "X-averaged negative secondary particle cracking rate [m.s-1]",
]
OUTPUT_VARIABLES = [
    "Current [A]",
    "Discharge capacity [A.h]",
    "Throughput capacity [A.h]",
    "Terminal voltage [V]",
    "Cell thickness change [m]",
    "Negative electrode transfer ratio k",
    "Volume-averaged cell temperature [K]",
    "Local ECM resistance [Ohm]",
    "X-averaged negative primary particle surface stoichiometry",
    "X-averaged negative secondary particle surface stoichiometry",
    "Loss of active material in negative electrode [%]",
    "Loss of active material in positive electrode [%]",
    "Loss of active material in primary phase in negative electrode [%]",
    "Loss of active material in secondary phase in negative electrode [%]",
    "Total capacity lost to side reactions [A.h]",
    "Loss of lithium due to loss of primary active material in negative electrode [mol]",
    "Loss of lithium due to loss of secondary active material in negative electrode [mol]",
    "Loss of lithium due to loss of active material in positive electrode [mol]",
    "Total lithium lost from particles [mol]",
] + CRACK_VARS + STRAIN_VARS

MIN_VALID_CAP_AH = 0.10 * MODEL_CAPACITY_AH


def _legs(cyc, lo, hi):
    """Largest discharge step with lo < mean current <= hi: (capacity, end throughput)."""
    cap, end_thr = 0.0, None
    for step in getattr(cyc, "steps", []):
        try:
            current = step["Current [A]"].entries
            q = step["Discharge capacity [A.h]"].entries
            thr = step["Throughput capacity [A.h]"].entries
        except (KeyError, TypeError, AttributeError):
            continue
        if current.size < 2 or q.size < 2:
            continue
        if lo < np.mean(current) <= hi and q[-1] - q[0] > cap:
            cap, end_thr = float(q[-1] - q[0]), float(thr[-1])
    return cap, end_thr


def ageing_discharge(cyc):
    """C/3-type discharge (mean current > 0.5 A) of an ageing cycle."""
    cap, end_thr = _legs(cyc, 0.5, np.inf)
    return (cap, end_thr) if cap >= MIN_VALID_CAP_AH else (0.0, None)


def rpt_discharge(cyc):
    """C/20-type discharge (0.05 < mean current <= 0.5 A) of an RPT cycle."""
    return _legs(cyc, 0.05, 0.5)


def rpt_discharge_curve(cyc):
    best_cap, best = 0.0, (None, None)
    for step in getattr(cyc, "steps", []):
        try:
            current = step["Current [A]"].entries
            q = step["Discharge capacity [A.h]"].entries
            v = step["Terminal voltage [V]"].entries
        except (KeyError, TypeError, AttributeError):
            continue
        if current.size >= 2 and q.size >= 2 and 0.05 < np.mean(current) <= 0.5 and q[-1] - q[0] > best_cap:
            best_cap, best = float(q[-1] - q[0]), (q - q[0], v)
    return best


def run(cond, data_dir, rate_basis_ah=RATE_BASIS_AH, efc_basis_ah=EFC_BASIS_AH, log=print):
    """Simulate one condition; returns the pybamm Solution."""
    model = pybamm.lithium_ion.DFN(dict(MODEL_OPTIONS))
    param = parameter_values(cond, data_dir, rate_basis_ah)
    solver = pybamm.IDAKLUSolver(root_tol=SOLVER_TOL, atol=SOLVER_TOL, rtol=SOLVER_TOL,
                                 output_variables=OUTPUT_VARIABLES)

    def make_sim(experiment, this_solver):
        return pybamm.Simulation(model, parameter_values=param, experiment=experiment,
                                 solver=this_solver, var_pts=VAR_PTS)

    def solve_with_retry(experiment, last_sol, label):
        last_exc = None
        for attempt, mult in enumerate([1.0] + RETRY_TOL_MULTIPLIERS):
            tol = SOLVER_TOL * mult
            this_solver = solver if mult == 1.0 else pybamm.IDAKLUSolver(
                root_tol=tol, atol=tol, rtol=tol, output_variables=OUTPUT_VARIABLES)
            try:
                sol = make_sim(experiment, this_solver).solve(starting_solution=last_sol)
                if attempt:
                    log(f"[{label}] recovered on retry {attempt} (tol={tol:.0e})")
                return sol
            except Exception as exc:  # noqa: BLE001
                last_exc = exc
        log(f"[{label}] solver failure, retries exhausted: {type(last_exc).__name__}: {str(last_exc)[:200]}")
        return None

    def efc_now(sol):
        return float(sol["Throughput capacity [A.h]"].entries[-1]) / (2 * efc_basis_ah)

    def efc_per_cycle(sol):
        for cyc in reversed(sol.cycles):
            cap, _ = ageing_discharge(cyc)
            if cap > 0:
                thr = cyc["Throughput capacity [A.h]"].entries
                return float(thr[-1] - thr[0]) / (2 * efc_basis_ah)
        q = sol.cycles[0]["Discharge capacity [A.h]"].entries
        return float(q.max() - q.min()) / MODEL_CAPACITY_AH

    age = ageing_cycle(cond, rate_basis_ah)
    rpt = (f"Discharge at {RPT_RATE} until {LOWER_CUTOFF_V} V", f"Charge at C/3 until {UPPER_CUTOFF_V} V",
           f"Hold at {UPPER_CUTOFF_V} V until C/50")
    if cond.partial:
        top = float(age[2].split("Hold at ")[1].split(" V")[0])
        cc = () if UPPER_CUTOFF_V - top < 0.1 else (f"Charge at C/3 until {UPPER_CUTOFF_V} V",)
        rpt = ("Rest for 10 minutes",) + cc + (f"Hold at {UPPER_CUTOFF_V} V until C/50",) + rpt
    log(f"{cond.label}: ageing cycle: " + " | ".join(age))
    log(f"{cond.label}: RPT cycle: " + " | ".join(rpt))

    # Three one-step cycles, so the formation RPT is placed at the end of its discharge.
    formation = pybamm.Experiment([f"Discharge at {RPT_RATE} until {LOWER_CUTOFF_V} V",
                                   f"Charge at C/3 until {UPPER_CUTOFF_V} V",
                                   f"Hold at {UPPER_CUTOFF_V} V until C/50"])
    sol = make_sim(formation, solver).solve(initial_soc=1.0)
    if cond.partial:  # conditioned reference RPT (the SoH reference)
        new = solve_with_retry(pybamm.Experiment([rpt[-3:]]), sol, "RPT@ref")
        if new is not None:
            sol = new

    n_age, stop, first_cap = 0, "all real RPT EFCs reached", None
    for target in data.rpt_targets(cond, data_dir):
        while True:
            e, epc = efc_now(sol), efc_per_cycle(sol)
            if target - e <= 0.5 * epc:
                break
            n = int(min(BATCH_SIZE, max(1, np.floor((target - e) / epc + 0.5))))
            if n_age + n > MAX_TOTAL_CYCLES:
                stop = "reached MAX_TOTAL_CYCLES"
                break
            new = solve_with_retry(pybamm.Experiment([age] * n), sol, f"AGE x{n}")
            if new is None:
                stop = f"solver failure before RPT at EFC {target:.0f}"
                break
            sol, n_age = new, n_age + n
            caps = [c for c, _ in (ageing_discharge(cyc) for cyc in sol.cycles) if c > 0]
            first_cap = first_cap or caps[0]
            soh = 100 * caps[-1] / first_cap
            log(f"[AGE] +{n} cycles ({n_age} total), EFC={efc_now(sol):.1f} -> {target:.1f}, SoH={soh:.2f}%")
            if soh <= SOH_FLOOR_PERCENT:
                stop = f"reached {SOH_FLOOR_PERCENT:.0f}% SoH floor"
                break
        if stop != "all real RPT EFCs reached":
            break
        new = solve_with_retry(pybamm.Experiment([rpt]), sol, f"RPT@{target:.1f}")
        if new is None:
            stop = f"solver failure at RPT EFC {target:.0f}"
            break
        sol = new
    log(f"{cond.label}: stop reason: {stop}")
    return sol


def _downsample_index(x, n):
    grid = np.linspace(x[0], x[-1], n)
    return np.unique(np.clip(np.searchsorted(x, grid), 0, x.size - 1))


def reduce(sol, cond, efc_basis_ah=EFC_BASIS_AH, n_series=20000):
    """Compact results: per-cycle and per-RPT quantities in full; the life
    time series (LAM, LLI, crack diagnostics; whole run including the final
    RPT) on n_series points evenly spaced in throughput. Everything the
    scores and figures need, nothing else."""
    age_cap, age_thr, age_exp, rpt_cap, rpt_thr, rpt_k, curves = [], [], [], [], [], [], []
    cyc_thr, cyc_sig = [], []
    for cyc in sol.cycles:
        thr_end = float(cyc["Throughput capacity [A.h]"].entries[-1])
        try:
            cyc_sig.append(float(np.max(np.abs(cyc[CRACK_VARS[2]].entries))) / 1e6)
            cyc_thr.append(thr_end)
        except Exception:  # noqa: BLE001
            pass
        cap, end_thr = ageing_discharge(cyc)
        if cap > 0 and end_thr is not None:
            age_cap.append(cap)
            age_thr.append(thr_end)
            d = cyc["Cell thickness change [m]"].entries
            age_exp.append(float(d.max() - d.min()) * 1e6 if d.size >= 2 else np.nan)
            continue
        rcap, rend = rpt_discharge(cyc)
        if rcap > 0 and rend is not None:
            rpt_cap.append(rcap)
            rpt_thr.append(thr_end)
            k = cyc["Negative electrode transfer ratio k"].entries
            rpt_k.append(float(np.max(np.abs(k))) if k.size else np.nan)
            q, v = rpt_discharge_curve(cyc)
            if q is not None:
                curves.append((thr_end, q, v))

    qt_full = sol["Throughput capacity [A.h]"].entries
    t = sol["Time [s]"].entries
    idx = _downsample_index(qt_full, n_series)

    def series(name, scale=1.0):
        return np.asarray(sol[name].entries, dtype=float)[idx] * scale

    def cumulative(name):  # time integral, at full resolution then sampled
        y = np.asarray(sol[name].entries, dtype=float)
        c = np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(t))])
        return c[idx]

    to_ah = F / 3600.0
    res = dict(
        age_cap=np.array(age_cap), age_thr=np.array(age_thr), age_expansion_um=np.array(age_exp),
        rpt_cap=np.array(rpt_cap), rpt_thr=np.array(rpt_thr), rpt_k=np.array(rpt_k),
        cycle_thr=np.array(cyc_thr), cycle_peak_stress_MPa=np.array(cyc_sig),
        thr=qt_full[idx],
        LAM_gr=series("Loss of active material in primary phase in negative electrode [%]"),
        LAM_si=series("Loss of active material in secondary phase in negative electrode [%]"),
        LAM_pos=series("Loss of active material in positive electrode [%]"),
        Q_side=series("Total capacity lost to side reactions [A.h]"),
        Q_lli_lam_trap=sum(series(f"Loss of lithium due to loss of {p}active material in {d} electrode [mol]", to_ah)
                           for p, d in (("primary ", "negative"), ("secondary ", "negative"), ("", "positive"))),
        Q_lli_particles=series("Total lithium lost from particles [mol]", to_ah),
        crack_length=series(CRACK_VARS[0]), roughness=series(CRACK_VARS[1]), porosity=series(CRACK_VARS[3]),
        Q_si_sei=series(CRACK_VARS[4], to_ah), Q_si_sei_cracks=series(CRACK_VARS[5], to_ah),
        strain_E=cumulative(STRAIN_VARS[0]), dl_strain=cumulative(STRAIN_VARS[1]),
        dl_total=cumulative(STRAIN_VARS[2]),
        soh_ref_index=np.array(1 if cond.partial else -1),
        efc_basis_ah=np.array(efc_basis_ah),
    )
    res["dl_paris"] = res["dl_total"] - res["dl_strain"]
    rpt_curves = [dict(thr=th, q=q, v=v) for th, q, v in curves]
    return res, rpt_curves
