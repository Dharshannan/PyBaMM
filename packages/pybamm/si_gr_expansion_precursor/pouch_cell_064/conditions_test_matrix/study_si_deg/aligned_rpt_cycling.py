# -*- coding: utf-8 -*-
"""aligned_rpt_cycling: a version of the fork's cycling driver
(high_temp_45C/cell064_degradation_fit_45C.py run_degradation) that puts the
model's C/20 RPTs at the REAL cell's RPT EFCs instead of every 50 cycles.

Why: the fork runs an RPT at the end of every 50-cycle batch, so after the
knee (where EFC per cycle falls with capacity) the model's RPTs land 10-25
EFC away from the real ones. The 25 degC late-RPT voltage gap changes by
~50 mV per model RPT, so nearest-EFC matching made the voltage comparison
depend on where the model's RPTs happened to fall.

How: ageing C/3 cycles (the same 3 steps as the fork's batches) run in
chunks. Before each chunk, the EFC per cycle is estimated from the last
ageing cycle's throughput, and the chunk is sized to land on the next real
RPT EFC. When the next target is within half a cycle, one RPT cycle runs
(C/20 discharge, C/3 charge, CV hold, the same as the fork's RPT cycle). The
RPT discharge therefore starts within +-0.5 ageing cycles of the real RPT
EFC (about +-0.5 EFC before the knee, +-0.3 after). Chunks are capped at
BATCH_SIZE cycles and use the fork's tolerance-retry logic. The run stops
after the last real RPT, or at the fork's SoH floor / cycle budget.

The formation discharge (C/20) is the RPT at EFC 0, as in the fork.

Use: import after the fork and call run_degradation_aligned(m, targets).
rpt_soc_plots.py --aligned does this with the per-temperature real RPT EFCs.
"""
import numpy as np
import pybamm


def real_rpt_efcs(m, t_arg):
    """Real RPT EFCs in the model's frame (EFC since RPT1), excluding 0,
    which is the formation discharge. C64_REAL_CELL_DIR (fork) takes
    precedence over the temperature default (CELL064 / CELL017)."""
    if m.REAL_CELL_DIR_NAME:
        cap = m.load_real_cell("capacity_fade")
    else:
        cap = m.load_experimental_capacity_fade() if t_arg == "25" else m.load_cell017("capacity_fade")
    targets = sorted(float(e) for e in cap["efc"] if e > 0.5)
    # Optional extension past the last real RPT (e.g. to see a delayed knee):
    # RPT_EXTEND_TO_EFC=800 adds model-only RPTs every RPT_EXTEND_STEP (default
    # 50) EFC after the last real one. The run still stops at the SoH floor.
    import os
    ext = float(os.environ.get("RPT_EXTEND_TO_EFC", 0) or 0)
    step = float(os.environ.get("RPT_EXTEND_STEP", 50))
    e = targets[-1] if targets else 0.0
    while ext and e + step <= ext:
        e += step
        targets.append(e)
    return targets


def run_degradation_aligned(m, targets, ageing_cycle=None):
    """ageing_cycle: optional tuple of pybamm experiment step strings for one
    ageing cycle (e.g. a partial-SoC window); default = the fork's full
    C/3 cycle (discharge to the lower cut-off, CC-CV charge to the upper)."""
    model = pybamm.lithium_ion.DFN(dict(m.MODEL_OPTIONS_BASE))
    param = m.build_parameter_values()

    def make_sim(experiment, this_solver):
        return pybamm.Simulation(model, parameter_values=param, experiment=experiment,
                                 solver=this_solver, var_pts=m.VAR_PTS)

    partial = ageing_cycle is not None
    ageing_cycle = ageing_cycle or (
        f"Discharge at C/3 until {m.LOWER_CUTOFF_V} V",
        f"Charge at C/3 until {m.UPPER_CUTOFF_V} V",
        f"Hold at {m.UPPER_CUTOFF_V} V until C/50",
    )
    print("Ageing cycle:", " | ".join(ageing_cycle), flush=True)
    rpt_cycle = (
        f"Discharge at {m.RPT_RATE} until {m.LOWER_CUTOFF_V} V",
        f"Charge at C/3 until {m.UPPER_CUTOFF_V} V",
        f"Hold at {m.UPPER_CUTOFF_V} V until C/50",
    )
    if partial:
        # A partial-window ageing cycle ends at the window top, not at 100%
        # SoC, so the RPT starts with a full charge (as the real test does
        # before every C/20 RPT): a 10 min rest, then CC-CV to the upper
        # cut-off. If the window top is within 0.1 V of it (CELL009: 4.15 V),
        # the CC step is dropped. There, even C/5 charging overpotential
        # exceeds the cut-off on the first sample, and that zero-length step
        # hits an off-by-one in pybamm Solution.__add__ with restricted
        # output variables (the solution's t is one shorter than its
        # variables). The full-window default already ends fully charged, so
        # it is left exactly as before.
        tops = [float(st.split("Hold at ")[1].split(" V")[0]) for st in ageing_cycle if st.startswith("Hold at ")]
        cc = () if tops and m.UPPER_CUTOFF_V - tops[0] < 0.1 else (f"Charge at C/3 until {m.UPPER_CUTOFF_V} V",)
        rpt_cycle = ("Rest for 10 minutes",) + cc + (f"Hold at {m.UPPER_CUTOFF_V} V until C/50",) + rpt_cycle
    print("RPT cycle:", " | ".join(rpt_cycle), flush=True)

    def solve_with_retry(experiment, last_sol, label):
        last_exc = None
        for attempt, tol_mult in enumerate([1.0] + m.RETRY_TOL_MULTIPLIERS):
            tol = m.SOLVER_TOL * tol_mult
            this_solver = (m.solver if tol_mult == 1.0 else pybamm.IDAKLUSolver(
                root_tol=tol, atol=tol, rtol=tol,
                output_variables=m._CORE_OUTPUT_VARIABLES if m._RESTRICT_OUTPUT_VARIABLES else None))
            try:
                sol = make_sim(experiment, this_solver).solve(starting_solution=last_sol)
                if attempt > 0:
                    print(f"[{label}] recovered on retry {attempt} (tol={tol:.0e})", flush=True)
                return sol
            except Exception as exc:  # noqa: BLE001
                last_exc = exc
        print(f"[{label}] solver failure, all retries exhausted: {type(last_exc).__name__}: "
              f"{str(last_exc)[:200]}", flush=True)
        return None

    def efc_now(sol):
        return float(m.efc_from_throughput(sol["Throughput capacity [A.h]"].entries[-1]))

    def efc_per_cycle(sol):
        """EFC of the last ageing cycle; before the first one, the formation
        discharge capacity stands in (EFC per C/3 cycle ~ capacity / Q_nom)."""
        for cyc in reversed(sol.cycles):
            cap, rate, _, _ = m.cycle_ageing_leg(cyc)
            if cap > 0 and rate > 0.5:
                thr = cyc["Throughput capacity [A.h]"].entries
                return float(m.efc_from_throughput(thr[-1] - thr[0]))
        q = sol.cycles[0]["Discharge capacity [A.h]"].entries
        return float(q.max() - q.min()) / m.NOMINAL_CAP_AH

    last_sol = make_sim(m.formation_exp, m.solver).solve(initial_soc=1.0)
    print("Formation cycle solved (RPT at EFC 0).", flush=True)
    if partial:
        # Conditioned reference RPT: formation ends fully CC-CV charged, so a
        # plain C/20 discharge + recharge here is a C/20 capacity measured the
        # same way as every later RPT (and as the real RPT1). It becomes the
        # SoH reference (m.RPT_SOH_REF_INDEX = 1); the formation discharge
        # from initial_soc=1.0 sits ~1.7% higher.
        ref_cycle = (f"Discharge at {m.RPT_RATE} until {m.LOWER_CUTOFF_V} V",
                     f"Charge at C/3 until {m.UPPER_CUTOFF_V} V",
                     f"Hold at {m.UPPER_CUTOFF_V} V until C/50")
        new = solve_with_retry(pybamm.Experiment([ref_cycle]), last_sol, "RPT@ref")
        if new is not None:
            last_sol = new
            m.RPT_SOH_REF_INDEX = 1
            print(f"[RPT] conditioned reference RPT at EFC {efc_now(last_sol):.1f}", flush=True)

    n_age = 0
    stop_reason = "all real RPT EFCs reached"
    first_cap = None
    for target in targets:
        while True:
            e, epc = efc_now(last_sol), efc_per_cycle(last_sol)
            remaining = target - e
            if remaining <= 0.5 * epc:
                break
            n = int(min(m.BATCH_SIZE, max(1, np.floor(remaining / epc + 0.5))))
            if n_age + n > m.MAX_TOTAL_CYCLES:
                stop_reason = "reached MAX_TOTAL_CYCLES"
                break
            new = solve_with_retry(pybamm.Experiment([ageing_cycle] * n), last_sol, f"AGE x{n}")
            if new is None:
                stop_reason = f"solver_failure_before_rpt_efc_{target:.0f}"
                break
            last_sol, n_age = new, n_age + n
            caps = [c for c, r, _, _ in (m.cycle_ageing_leg(cyc) for cyc in last_sol.cycles) if c > 0 and r > 0.5]
            first_cap = first_cap or caps[0]
            soh = 100 * caps[-1] / first_cap
            print(f"[AGE] +{n} cycles ({n_age} total), EFC={efc_now(last_sol):.1f} -> target {target:.1f}, "
                  f"SoH={soh:.2f}%", flush=True)
            if soh <= m.SOH_TERMINATION_PERCENT:
                stop_reason = f"reached_{m.SOH_TERMINATION_PERCENT:.0f}pct_floor"
                break
        if stop_reason != "all real RPT EFCs reached":
            break
        start = efc_now(last_sol)
        new = solve_with_retry(pybamm.Experiment([rpt_cycle]), last_sol, f"RPT@{target:.1f}")
        if new is None:
            stop_reason = f"solver_failure_at_rpt_efc_{target:.0f}"
            break
        last_sol = new
        print(f"[RPT] real EFC {target:.1f}: model RPT discharge starts at EFC {start:.1f} "
              f"(offset {start - target:+.2f})", flush=True)

    print(f"Stop reason: {stop_reason}", flush=True)
    return last_sol
