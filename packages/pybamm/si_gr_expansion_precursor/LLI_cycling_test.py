import matplotlib.pyplot as plt
import numpy as np
import pybamm

"""
NOTICE: This script is to test LLI modes to simulate knee due to porosity failure.
Parameter set: Mayur2024 (Newly added).
"""

pybamm.set_logging_level('NOTICE')

# ============================================================================
# CONFIG
# ============================================================================
N = 50  # ageing cycles per block (RPT every N cycles)
M_MAX = 40  # max number of ageing+RPT blocks (hard cap)
SOH_STOP = 0.70  # stop when discharge-capacity SoH falls to this fraction
VMIN = 2.5  # discharge cutoff [V]
VMAX = 4.2  # charge cutoff [V]

# ============================================================================
# MODEL / PARAMETERS  (same composite set-up as the base script)
# ============================================================================
model = pybamm.lithium_ion.DFN(
    {
        "particle phases": ("2", "1"),
        "open-circuit potential": (("single", "current sigmoid"), "single"),
        "SEI": "solvent-diffusion limited",
        "SEI porosity change": "true",
        "lithium plating": "partially reversible",
        "lithium plating porosity change": "true",
        "particle mechanics": ("swelling only", "swelling only"),
        # "SEI on cracks": "false",
        # "loss of active material": "stress-driven",
    }
)

param = pybamm.ParameterValues("Mayur2024")
# Update for more severe degradation
param.update(
    {
        # --- SEI growth (gradual LLI) ---
        "Primary: SEI solvent diffusivity [m2.s-1]": 2.5e-22,    # was 2.5e-22
        "Secondary: SEI solvent diffusivity [m2.s-1]": 2.5e-22,  # was 2.5e-22
        "Primary: SEI kinetic rate constant [m.s-1]": 1e-12,     # was 1e-12
        "Secondary: SEI kinetic rate constant [m.s-1]": 1e-12,   # was 1e-12
        # --- lithium plating (accelerating LLI -> knee lever) ---
        "Lithium plating kinetic rate constant [m.s-1]": 1e-08 * 0.5,  # was 1e-09
        "Dead lithium decay constant [s-1]": 1e-05 * 0.5,              # was 1e-06
        # --- silicon cracking ---
        # "Secondary: Negative electrode cracking rate": 3.9e-20 / 100,
    }
)

var_pts = {
    "x_n": 5, "x_s": 5, "x_p": 5,
    "r_n_prim": 15, "r_n_sec": 15, "r_n": 15, "r_p": 15,
}
solver = pybamm.IDAKLUSolver(root_tol=1e-04, rtol=1e-04, atol=1e-04)


def make_sim(experiment):
    return pybamm.Simulation(
        model, parameter_values=param, experiment=experiment,
        solver=solver, var_pts=var_pts,
    )


# ============================================================================
# EXPERIMENTS
# ============================================================================
# Ageing block: standard 1C cycle (ends discharged + rested).
ageing_experiment = pybamm.Experiment(
    [
        (
            f"Charge at 2C until {VMAX}V",
            f"Hold at {VMAX}V until C/100",
            f"Discharge at 1C until {VMIN}V",
            "Rest for 1 hour",
        )
    ]
    * N
)
# Bring the cell to a well-defined full state before each RPT.
charge_experiment = pybamm.Experiment(
    [(f"Charge at 0.3C until {VMAX}V", f"Hold at {VMAX}V until C/100")]
)
# RPT: slow C/3 discharge to measure capacity (single cycle, single step).
rpt_experiment = pybamm.Experiment([(f"Discharge at C/3 until {VMIN}V",)])


# ============================================================================
# HELPERS
# ============================================================================
def rpt_discharge_step(rpt_sol):
    """The C/3 discharge step of the most recent RPT."""
    return rpt_sol.cycles[-1].steps[0]


def step_discharge_capacity(step):
    q = step["Discharge capacity [A.h]"].entries
    return float(q[-1] - q[0])


# ============================================================================
# RUN: initial RPT, then ageing + RPT blocks
# ============================================================================
rpt_records = []  # list of dicts: cycle, throughput_Ah, cap_Ah, soh, Q(Ah), V

print("=" * 60)
print("Initial capacity check (RPT at BoL)")
print("=" * 60)
sim = make_sim(rpt_experiment)
rpt_sol = sim.solve(initial_soc=1.0)  # discharge from full -> C0
step0 = rpt_discharge_step(rpt_sol)
C0 = step_discharge_capacity(step0)
Qd0 = step0["Discharge capacity [A.h]"].entries
rpt_records.append({
    "cycle": 0,
    "throughput_Ah": float(rpt_sol["Throughput capacity [A.h]"].entries[-1]),
    "cap_Ah": C0, "soh": 1.0,
    "Q": Qd0 - Qd0[0], "V": step0["Voltage [V]"].entries,
})
print(f"  Initial capacity C0 = {C0:.4f} A.h")

last_sol = rpt_sol
stop_reason = "reached M_MAX"

for i in range(M_MAX):
    n_cycles = (i + 1) * N
    print("=" * 60)
    print(f"Block {i + 1}/{M_MAX}:  ageing cycles {i * N + 1}–{n_cycles}, then RPT")
    print("=" * 60)
    try:
        # --- ageing block ---
        sim = make_sim(ageing_experiment)
        ageing_sol = sim.solve(starting_solution=last_sol)
        # --- full charge before RPT ---
        sim = make_sim(charge_experiment)
        charge_sol = sim.solve(starting_solution=ageing_sol)
        # --- RPT ---
        sim = make_sim(rpt_experiment)
        rpt_sol = sim.solve(starting_solution=charge_sol)
    except Exception as exc:  # noqa: BLE001
        print(f"  [WARN] solve failed in block {i + 1}: {exc}")
        stop_reason = f"solver_failure_block_{i + 1}"
        break

    step = rpt_discharge_step(rpt_sol)
    cap = step_discharge_capacity(step)
    soh = cap / C0
    Qd = step["Discharge capacity [A.h]"].entries
    rpt_records.append({
        "cycle": n_cycles,
        "throughput_Ah": float(rpt_sol["Throughput capacity [A.h]"].entries[-1]),
        "cap_Ah": cap, "soh": soh,
        "Q": Qd - Qd[0], "V": step["Voltage [V]"].entries,
    })
    print(f"  RPT capacity = {cap:.4f} A.h   SoH = {soh * 100:.2f}%")

    last_sol = rpt_sol
    if soh <= SOH_STOP:
        stop_reason = f"SoH<={SOH_STOP * 100:.0f}% at cycle {n_cycles}"
        break

final_sol = rpt_sol
print("\nStopped:", stop_reason)
print(f"RPTs recorded: {len(rpt_records)}  "
      f"(final SoH = {rpt_records[-1]['soh'] * 100:.2f}%)")

# ============================================================================
# PLOTS — continuous degradation vs throughput (from the cumulative solution)
# ============================================================================
Qt = final_sol["Throughput capacity [A.h]"].entries

# Plot 1: capacity-loss breakdown
Q_SEI = final_sol["Loss of capacity to negative SEI [A.h]"].entries
Q_plating = final_sol["Loss of capacity to negative lithium plating [A.h]"].entries
Q_side = final_sol["Total capacity lost to side reactions [A.h]"].entries
Q_LLI = final_sol["Total lithium lost [mol]"].entries * 96485.3 / 3600
plt.figure()
plt.plot(Qt, Q_SEI, label="SEI", linestyle="dashed")
plt.plot(Qt, Q_plating, label="Li plating", linestyle="dotted")
plt.plot(Qt, Q_side, label="All side reactions", linestyle=(0, (6, 1)))
plt.plot(Qt, Q_LLI, label="All LLI")
plt.xlabel("Throughput capacity [A.h]")
plt.ylabel("Capacity loss [A.h]")
plt.title("Capacity-loss breakdown")
plt.legend()
plt.tight_layout()

# Plot 2: degradation modes
LLI = final_sol["Loss of lithium inventory [%]"].entries
LAM_neg = final_sol["Loss of active material in negative electrode [%]"].entries
LAM_pos = final_sol["Loss of active material in positive electrode [%]"].entries
plt.figure()
plt.plot(Qt, LLI, label="LLI")
plt.plot(Qt, LAM_neg, label="LAM (negative)")
plt.plot(Qt, LAM_pos, label="LAM (positive)")
plt.xlabel("Throughput capacity [A.h]")
plt.ylabel("Degradation modes [%]")
plt.title("Degradation modes")
plt.legend()
plt.tight_layout()

# Plot 2b: negative-electrode LAM by phase
LAM_gr = final_sol["Loss of active material in primary phase in negative electrode [%]"].entries
LAM_si = final_sol["Loss of active material in secondary phase in negative electrode [%]"].entries
plt.figure()
plt.plot(Qt, LAM_gr, label="LAM graphite (primary)")
plt.plot(Qt, LAM_si, label="LAM silicon (secondary)", linestyle="dashed")
plt.xlabel("Throughput capacity [A.h]")
plt.ylabel("Loss of active material [%]")
plt.title("Negative-electrode LAM by phase")
plt.legend()
plt.tight_layout()

# Plot 3: negative electrode porosity
eps_neg_avg = final_sol["X-averaged negative electrode porosity"].entries
eps_neg_sep = final_sol["Negative electrode porosity"].entries[-1, :]
eps_neg_CC = final_sol["Negative electrode porosity"].entries[0, :]
plt.figure()
plt.plot(Qt, eps_neg_avg, label="Average")
plt.plot(Qt, eps_neg_sep, label="Separator side", linestyle="dotted")
plt.plot(Qt, eps_neg_CC, label="Current-collector side", linestyle="dashed")
plt.xlabel("Throughput capacity [A.h]")
plt.ylabel("Negative electrode porosity")
plt.title("Negative electrode porosity")
plt.legend()
plt.tight_layout()

# ============================================================================
# NEW Plot 4: RPT discharge voltage curves over life
# ============================================================================
plt.figure()
cmap = plt.cm.viridis
nrpt = len(rpt_records)
for k, rec in enumerate(rpt_records):
    c = cmap(k / max(nrpt - 1, 1))
    plt.plot(rec["Q"], rec["V"], color=c, lw=1.4,
             label=f"cycle {rec['cycle']} (SoH {rec['soh'] * 100:.0f}%)")
plt.xlabel("Discharge capacity [A.h]")
plt.ylabel("Voltage [V]")
plt.title("RPT discharge voltage curves (C/3) over life")
plt.legend(fontsize=7, loc="lower left")
plt.tight_layout()

# ============================================================================
# NEW Plot 5: SoH (from RPT discharge capacity) vs cycle number  — knee check
# ============================================================================
rpt_cyc = np.array([r["cycle"] for r in rpt_records], float)
rpt_soh = np.array([r["soh"] for r in rpt_records], float) * 100.0
rpt_thr = np.array([r["throughput_Ah"] for r in rpt_records], float)

fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
ax[0].plot(rpt_cyc, rpt_soh, "o-", lw=1.8, ms=5)
ax[0].axhline(80, color="orange", ls=":", lw=1, label="80% SoH")
ax[0].axhline(SOH_STOP * 100, color="red", ls="--", lw=1, label=f"{SOH_STOP * 100:.0f}% SoH")
ax[0].set_xlabel("Cycle number")
ax[0].set_ylabel("SoH [%]  (RPT discharge capacity)")
ax[0].set_title("SoH vs cycle — knee check")
ax[0].legend(fontsize=8)
ax[0].grid(True, alpha=0.3)

ax[1].plot(rpt_thr, rpt_soh, "s-", color="teal", lw=1.8, ms=5)
ax[1].axhline(80, color="orange", ls=":", lw=1)
ax[1].axhline(SOH_STOP * 100, color="red", ls="--", lw=1)
ax[1].set_xlabel("Throughput capacity [A.h]")
ax[1].set_ylabel("SoH [%]")
ax[1].set_title("SoH vs throughput")
ax[1].grid(True, alpha=0.3)
fig.tight_layout()

plt.show()


# ============================================================================
# NEW Plot 6: EVERY cycle's discharge capacity (not just the RPTs)
# ----------------------------------------------------------------------------
# Walks the full cumulative solution and pulls the discharge capacity of every
# cycle so you can see how fast the collapse actually happens at cycle
# resolution (the RPTs only sample it every N cycles).
# ============================================================================
def _valid_steps(cyc):
    """Non-empty, subscriptable steps of a cycle (skip EmptySolution / failures)."""
    good = []
    for step in getattr(cyc, "steps", []):
        try:
            _ = step["Current [A]"].entries  # probe: raises for EmptySolution
        except (KeyError, TypeError, AttributeError):
            continue
        good.append(step)
    return good


def discharge_capacity_of_steps(steps):
    """Sum discharge capacity over the discharging step(s) [A.h]."""
    total = 0.0
    for step in steps:
        try:
            I = step["Current [A]"].entries
            q = step["Discharge capacity [A.h]"].entries
        except (KeyError, TypeError, AttributeError):
            continue
        if q.size >= 2 and np.mean(I) > 0:  # discharging (PyBaMM: I>0 = discharge)
            total += float(q[-1] - q[0])
    return total


ageing_x, ageing_cap = [], []
rpt_x2, rpt_cap2 = [], []
ageing_count = 0
for cyc in final_sol.cycles:
    steps = _valid_steps(cyc)
    nsteps = len(steps)
    if nsteps == 0:
        continue  # fully empty / failed cycle -> skip
    dcap = discharge_capacity_of_steps(steps)
    if nsteps == 1:  # RPT: single C/3 discharge
        rpt_x2.append(ageing_count)
        rpt_cap2.append(dcap)
    elif nsteps >= 3:  # ageing cycle: charge/hold/discharge/rest
        ageing_count += 1
        ageing_x.append(ageing_count)
        ageing_cap.append(dcap)
    # 2-step charge-only cycles have no discharge -> skipped

ageing_x = np.array(ageing_x, float)
ageing_cap = np.array(ageing_cap, float)

fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
ax[0].plot(ageing_x, ageing_cap, "-", color="steelblue", lw=1,
           label="Ageing cycles (1C discharge)")
if rpt_cap2:
    ax[0].scatter(rpt_x2, rpt_cap2, color="crimson", s=18, zorder=3,
                  label="RPT (C/3 discharge)")
ax[0].set_xlabel("Ageing cycle number")
ax[0].set_ylabel("Discharge capacity [A.h]")
ax[0].set_title("Discharge capacity — every cycle")
ax[0].legend(fontsize=8)
ax[0].grid(True, alpha=0.3)

if ageing_cap.size:
    soh_cycle = 100 * ageing_cap / ageing_cap[0]
    ax[1].plot(ageing_x, soh_cycle, "-", color="steelblue", lw=1)
    ax[1].axhline(80, color="orange", ls=":", lw=1, label="80%")
    ax[1].axhline(SOH_STOP * 100, color="red", ls="--", lw=1, label=f"{SOH_STOP * 100:.0f}%")
    ax[1].set_xlabel("Ageing cycle number")
    ax[1].set_ylabel("SoH [%]  (per-cycle 1C discharge cap)")
    ax[1].set_title("Per-cycle SoH — collapse speed")
    ax[1].legend(fontsize=8)
    ax[1].grid(True, alpha=0.3)
fig.tight_layout()

plt.show()

# ============================================================================
# Summary table (handy for the tuning step)
# ============================================================================
print("\n  cycle   throughput[A.h]   cap[A.h]    SoH[%]")
for r in rpt_records:
    print(f"  {r['cycle']:6d}   {r['throughput_Ah']:13.2f}   "
          f"{r['cap_Ah']:8.4f}   {r['soh'] * 100:7.2f}")
