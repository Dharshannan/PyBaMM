"""Run one condition: formation, then C/3 ageing in chunks with a C/20 RPT
placed at each real RPT's EFC (aligned RPTs), collecting compact results
chunk by chunk.

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
from .parameters import MODEL_OPTIONS, V6, VAR_PTS, parameter_values

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
    "Loss of lithium to negative secondary SEI [mol]",
    "Loss of lithium to negative secondary SEI on cracks [mol]",
]
# Negative electrode porosity (x-averaged). It oscillates within each cycle
# (Si swelling buffered into the pores); its per-cycle minimum (charged end)
# decreases smoothly from the initial porosity to the porosity floor.
POROSITY = "X-averaged negative electrode porosity"
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
    POROSITY,
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


class _Collector:
    """Builds the compact results chunk by chunk, so only the current chunk's
    solution is ever held in memory. Per-cycle and per-RPT quantities are kept
    in full; the life time series (LAM, LLI, crack diagnostics) are sampled
    every dthr A.h of throughput; time integrals are accumulated at full
    resolution across chunks."""

    SERIES = ("thr", "LAM_gr", "LAM_si", "LAM_pos", "Q_side", "Q_lli_lam_trap", "Q_lli_particles",
              "crack_length", "roughness", "Q_si_sei", "Q_si_sei_cracks", "strain_E", "dl_strain", "dl_total")

    def __init__(self, cond, efc_basis_ah, dthr=0.1):
        self.cond, self.efc_basis_ah, self.dthr = cond, efc_basis_ah, dthr
        self.age_cap, self.age_thr, self.age_porosity, self.age_cycle_thr = [], [], [], []
        self.rpt_cap, self.rpt_thr, self.curves, self.cyc_thr, self.cyc_sig = [], [], [], [], []
        self.series = {k: [] for k in self.SERIES}
        self.integral = {name: 0.0 for name in STRAIN_VARS}
        self.next_thr = None
        self.initial_porosity = None

    def add(self, sol):
        for cyc in sol.cycles:
            if cyc["Time [s]"].entries.size < 2:  # the one-point cycle a chunk starts from
                continue
            thr = cyc["Throughput capacity [A.h]"].entries
            thr_end = float(thr[-1])
            self.cyc_sig.append(float(np.max(np.abs(cyc[CRACK_VARS[2]].entries))) / 1e6)
            self.cyc_thr.append(thr_end)
            cap, end_thr = ageing_discharge(cyc)
            if cap > 0 and end_thr is not None:
                self.age_cap.append(cap)
                self.age_thr.append(thr_end)
                self.age_porosity.append(float(np.min(cyc[POROSITY].entries)))
                self.age_cycle_thr.append(float(thr[-1] - thr[0]))
                continue
            rcap, rend = rpt_discharge(cyc)
            if rcap > 0 and rend is not None:
                self.rpt_cap.append(rcap)
                self.rpt_thr.append(thr_end)
                q, v = rpt_discharge_curve(cyc)
                if q is not None:
                    self.curves.append(dict(thr=thr_end, q=q, v=v))

        thr = np.asarray(sol["Throughput capacity [A.h]"].entries, dtype=float)
        t = np.asarray(sol["Time [s]"].entries, dtype=float)
        if self.initial_porosity is None:
            self.initial_porosity = float(sol[POROSITY].entries[0])
        cum = {}
        for name in STRAIN_VARS:  # running time integrals, at full resolution
            y = np.asarray(sol[name].entries, dtype=float)
            cum[name] = self.integral[name] + np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(t))])
            self.integral[name] = float(cum[name][-1])
        if self.next_thr is None:
            self.next_thr = thr[0]
        grid = np.arange(self.next_thr, thr[-1] + 1e-12, self.dthr)
        if grid.size == 0:
            return
        idx = np.unique(np.clip(np.searchsorted(thr, grid), 0, thr.size - 1))
        self.next_thr = grid[-1] + self.dthr
        to_ah = F / 3600.0

        def col(name, scale=1.0):
            return np.asarray(sol[name].entries, dtype=float)[idx] * scale

        new = dict(
            thr=thr[idx],
            LAM_gr=col("Loss of active material in primary phase in negative electrode [%]"),
            LAM_si=col("Loss of active material in secondary phase in negative electrode [%]"),
            LAM_pos=col("Loss of active material in positive electrode [%]"),
            Q_side=col("Total capacity lost to side reactions [A.h]"),
            Q_lli_lam_trap=sum(col(f"Loss of lithium due to loss of {p}active material in {d} electrode [mol]", to_ah)
                               for p, d in (("primary ", "negative"), ("secondary ", "negative"), ("", "positive"))),
            Q_lli_particles=col("Total lithium lost from particles [mol]", to_ah),
            crack_length=col(CRACK_VARS[0]), roughness=col(CRACK_VARS[1]),
            Q_si_sei=col(CRACK_VARS[3], to_ah), Q_si_sei_cracks=col(CRACK_VARS[4], to_ah),
            strain_E=cum[STRAIN_VARS[0]][idx], dl_strain=cum[STRAIN_VARS[1]][idx], dl_total=cum[STRAIN_VARS[2]][idx],
        )
        for k, v in new.items():
            self.series[k].append(v)

    def efc_per_cycle(self):
        """EFC of the last ageing cycle (None before the first one)."""
        return self.age_cycle_thr[-1] / (2 * self.efc_basis_ah) if self.age_cycle_thr else None

    def finish(self):
        res = {k: np.concatenate(v) for k, v in self.series.items()}
        res.update(
            age_cap=np.array(self.age_cap), age_thr=np.array(self.age_thr),
            age_min_porosity=np.array(self.age_porosity), initial_porosity=np.array(self.initial_porosity),
            rpt_cap=np.array(self.rpt_cap), rpt_thr=np.array(self.rpt_thr),
            cycle_thr=np.array(self.cyc_thr), cycle_peak_stress_MPa=np.array(self.cyc_sig),
            soh_ref_index=np.array(1 if self.cond.partial else -1), efc_basis_ah=np.array(self.efc_basis_ah),
        )
        res["dl_paris"] = res["dl_total"] - res["dl_strain"]
        return res, self.curves


def run(cond, data_dir, rate_basis_ah=RATE_BASIS_AH, efc_basis_ah=EFC_BASIS_AH, log=print, options=None, params=None):
    """Simulate one condition; returns (results, rpt_curves) (see _Collector).
    Each chunk starts from the previous chunk's last state only (which is all
    PyBaMM uses to continue a run), so memory use stays flat over the run.
    options / params: overrides of MODEL_OPTIONS / V6 (e.g. for ablations)."""
    model = pybamm.lithium_ion.DFN({**MODEL_OPTIONS, **(options or {})})
    param = parameter_values(cond, data_dir, rate_basis_ah, p={**V6, **(params or {})})
    solver = pybamm.IDAKLUSolver(root_tol=SOLVER_TOL, atol=SOLVER_TOL, rtol=SOLVER_TOL,
                                 output_variables=OUTPUT_VARIABLES)

    def make_sim(experiment, this_solver):
        return pybamm.Simulation(model, parameter_values=param, experiment=experiment,
                                 solver=this_solver, var_pts=VAR_PTS)

    def solve_with_retry(experiment, state, label):
        last_exc = None
        for attempt, mult in enumerate([1.0] + RETRY_TOL_MULTIPLIERS):
            tol = SOLVER_TOL * mult
            this_solver = solver if mult == 1.0 else pybamm.IDAKLUSolver(
                root_tol=tol, atol=tol, rtol=tol, output_variables=OUTPUT_VARIABLES)
            try:
                sol = make_sim(experiment, this_solver).solve(starting_solution=state)
                if attempt:
                    log(f"[{label}] recovered on retry {attempt} (tol={tol:.0e})")
                return sol
            except Exception as exc:  # noqa: BLE001
                last_exc = exc
        log(f"[{label}] solver failure, retries exhausted: {type(last_exc).__name__}: {str(last_exc)[:200]}")
        return None

    def efc_now(state):
        return float(state["Throughput capacity [A.h]"].entries[-1]) / (2 * efc_basis_ah)

    age = ageing_cycle(cond, rate_basis_ah)
    rpt = (f"Discharge at {RPT_RATE} until {LOWER_CUTOFF_V} V", f"Charge at C/3 until {UPPER_CUTOFF_V} V",
           f"Hold at {UPPER_CUTOFF_V} V until C/50")
    if cond.partial:
        top = float(age[2].split("Hold at ")[1].split(" V")[0])
        cc = () if UPPER_CUTOFF_V - top < 0.1 else (f"Charge at C/3 until {UPPER_CUTOFF_V} V",)
        rpt = ("Rest for 10 minutes",) + cc + (f"Hold at {UPPER_CUTOFF_V} V until C/50",) + rpt
    log(f"{cond.label}: ageing cycle: " + " | ".join(age))
    log(f"{cond.label}: RPT cycle: " + " | ".join(rpt))

    out = _Collector(cond, efc_basis_ah)
    # Three one-step cycles, so the formation RPT is placed at the end of its discharge.
    formation = pybamm.Experiment([f"Discharge at {RPT_RATE} until {LOWER_CUTOFF_V} V",
                                   f"Charge at C/3 until {UPPER_CUTOFF_V} V",
                                   f"Hold at {UPPER_CUTOFF_V} V until C/50"])
    sol = make_sim(formation, solver).solve(initial_soc=1.0)
    q = sol.cycles[0]["Discharge capacity [A.h]"].entries
    epc_formation = float(q.max() - q.min()) / MODEL_CAPACITY_AH  # EFC per cycle before any ageing
    out.add(sol)
    state = sol.last_state
    del sol
    if cond.partial:  # conditioned reference RPT (the SoH reference)
        new = solve_with_retry(pybamm.Experiment([rpt[-3:]]), state, "RPT@ref")
        if new is not None:
            out.add(new)
            state = new.last_state

    n_age, stop = 0, "all real RPT EFCs reached"
    for target in data.rpt_targets(cond, data_dir):
        while True:
            e, epc = efc_now(state), out.efc_per_cycle() or epc_formation
            if target - e <= 0.5 * epc:
                break
            n = int(min(BATCH_SIZE, max(1, np.floor((target - e) / epc + 0.5))))
            if n_age + n > MAX_TOTAL_CYCLES:
                stop = "reached MAX_TOTAL_CYCLES"
                break
            new = solve_with_retry(pybamm.Experiment([age] * n), state, f"AGE x{n}")
            if new is None:
                stop = f"solver failure before RPT at EFC {target:.0f}"
                break
            out.add(new)
            state, n_age = new.last_state, n_age + n
            del new
            soh = 100 * out.age_cap[-1] / out.age_cap[0]
            log(f"[AGE] +{n} cycles ({n_age} total), EFC={efc_now(state):.1f} -> {target:.1f}, SoH={soh:.2f}%")
            if soh <= SOH_FLOOR_PERCENT:
                stop = f"reached {SOH_FLOOR_PERCENT:.0f}% SoH floor"
                break
        if stop != "all real RPT EFCs reached":
            break
        new = solve_with_retry(pybamm.Experiment([rpt]), state, f"RPT@{target:.1f}")
        if new is None:
            stop = f"solver failure at RPT EFC {target:.0f}"
            break
        out.add(new)
        state = new.last_state
        del new
        log(f"[RPT] at EFC {efc_now(state):.1f} (real RPT at {target:.1f})")
    log(f"{cond.label}: stop reason: {stop}")
    return out.finish()
