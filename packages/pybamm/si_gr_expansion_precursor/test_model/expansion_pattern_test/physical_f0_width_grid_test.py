"""
physical_f0_width_grid_test.py -- follow-up to A.0k: combines a lower f0 (BoL plateau)
with the "sweet spot" widths found there (0.12, 0.17), still on the existing, unmodified
"physical" option -- no source changes needed.

Rationale: A.0k found width=0.12 got closest to zero (-0.00029 um/A.h) but never crossed,
using the recipe's existing f0=0.7. Lower f0 makes k_cmax=(1-f0)/f0 LARGER (f0=0.4 ->
k_cmax=1.5 vs f0=0.7 -> k_cmax=0.4286), which makes f MORE sensitive to a given change in
c_pore_normalised -- i.e. amplifies whatever gradual change in c the widened width already
produces over pre-transition life. Since width=0.12/0.17 already gives c_pore_normalised a
genuine (if modest) slope from early life (A.0k), a larger k_cmax could amplify that into
a big enough f swing to finally cross the zero-slope threshold that was so close at f0=0.7.

Also changes the *starting point*: lower f0 means more of BoL swelling is buffered away,
so the absolute buffered amplitude starts lower -- worth watching in the plots, not just
the slope number, since a very low f0 could suppress the whole curve into an uninteresting
regime even if the *sign* of the slope improves.

Grid: f0 in {0.7 (reference, already known from A.0h/A.0j/A.0k), 0.6, 0.5, 0.4} x
width in {0.12, 0.17}. Same base-case recipe (SI_MULT=18, LAM_MULT=1.0), "pore buffering
transition": "physical" throughout. Each case now runs to a 50% SoH floor (matching A.0e's
extended-run convention) rather than a fixed 120-cycle window, so the full trajectory
(including well past the knee) can be inspected -- not just the early-life slope."""
import os
import matplotlib.pyplot as plt
import numpy as np
import pybamm

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

pybamm.set_logging_level("NOTICE")

MODEL_OPTIONS = {
    "particle phases": ("2", "1"),
    "open-circuit potential": (("single", "current sigmoid"), "single"),
    "SEI": "reaction limited",
    "SEI porosity change": "true",
    "particle mechanics": ("swelling and cracking", "swelling only"),
    "SEI on cracks": "true",
    "loss of active material": "stress-driven",
    "pore buffering": "true",
    "pore buffering transition": "physical",  # <-- existing option, unmodified
}

UNBUFFERED_OPTIONS = dict(MODEL_OPTIONS, **{"pore buffering": "false"})

CAPACITY_PROBE_OPTIONS = {
    "particle phases": ("2", "1"),
    "open-circuit potential": (("single", "current sigmoid"), "single"),
}

VAR_PTS = {
    "x_n": 5, "x_s": 5, "x_p": 5,
    "r_n_prim": 20, "r_n_sec": 20, "r_n": 20, "r_p": 20,
}

SI_MAX_CONC_DEFAULT = 278000.0
NOMINAL_CAP_AH = 5.0


def scaled_si_initial_conc(si_max_conc):
    return 275220.0 * (si_max_conc / SI_MAX_CONC_DEFAULT)


def measure_capacity_constants(si_max_conc):
    param = pybamm.ParameterValues("si_gr_expansion")
    param.update(
        {
            "Secondary: Maximum concentration in negative electrode [mol.m-3]": si_max_conc,
            "Secondary: Initial concentration in negative electrode [mol.m-3]": scaled_si_initial_conc(si_max_conc),
        },
        check_already_exists=False,
    )
    EPS_PROBE = 0.1
    p = param.copy()
    p.update(
        {
            "Primary: Negative electrode active material volume fraction": EPS_PROBE,
            "Secondary: Negative electrode active material volume fraction": EPS_PROBE,
        }
    )
    model = pybamm.lithium_ion.DFN(CAPACITY_PROBE_OPTIONS)
    sim = pybamm.Simulation(
        model, parameter_values=p, var_pts=VAR_PTS,
        experiment=pybamm.Experiment(["Rest for 1 second"]),
    )
    sol = sim.solve()
    cap_gr = float(np.squeeze(sol["Negative electrode primary phase capacity [A.h]"].entries[0]))
    cap_si = float(np.squeeze(sol["Negative electrode secondary phase capacity [A.h]"].entries[0]))
    return cap_gr / EPS_PROBE, cap_si / EPS_PROBE


def solve_capacity_split_at_fixed_volume(si_vol_frac, cap_gr_frac, k_gr, target_cap):
    gr_eps = target_cap * cap_gr_frac / k_gr
    si_eps = (si_vol_frac / (1 - si_vol_frac)) * gr_eps
    cap_si_target = target_cap * (1 - cap_gr_frac)
    k_si_needed = cap_si_target / si_eps
    si_max_conc_needed = SI_MAX_CONC_DEFAULT * (k_si_needed / k_si_orig_global)
    return gr_eps, si_eps, si_max_conc_needed


print("Measuring capacity constants at DEFAULT Si max concentration (pure Si)...")
k_gr_orig, k_si_orig = measure_capacity_constants(SI_MAX_CONC_DEFAULT)
k_si_orig_global = k_si_orig

CAP_GR_FRAC = 0.55
SI_VOL_FRAC = 0.20
GR_EPS, SI_EPS, SI_MAX_CONC_NEEDED = solve_capacity_split_at_fixed_volume(
    SI_VOL_FRAC, CAP_GR_FRAC, k_gr_orig, NOMINAL_CAP_AH)
print(f"Solved composition: Gr_eps={GR_EPS:.5f}  Si_eps={SI_EPS:.5f}  "
      f"Si_max_conc={SI_MAX_CONC_NEEDED:.1f} mol/m3")

# Base-case recipe -- IDENTICAL to A.0h/A.0j/A.0k's baseline point.
SI_MULT = 18.0
LAM_MULT = 1.0
SI_CRIT_STRESS = 2.2e8
SI_LAM_PROP_BASELINE = 7.5e-7
SI_LAM_EXP = 2.5
GR_DIV = 4.0
GR_CRIT_STRESS = 3.0e7
GR_LAM_PROP_BASELINE = 4.5e-6
SI_CRACK_RATE_MULT = 0.1
GR_CRACK_RATE_MULT = 0.1
BASE_CRACK_RATE = 3.9e-20
NEG_POROSITY_FLOOR = 0.01
EXPONENT_MAX_SEI = 10.0

BATCH_SIZE = 20
MAX_TOTAL_CYCLES = 400  # safety cap only -- runs normally stop at SOH_TARGET first
SOH_TARGET = 0.5  # stop once discharge capacity falls to 50% of the first cycle's

F0_SWEEP = [0.6, 0.5, 0.4]  # 0.7 skipped -- already known from A.0h/A.0j/A.0k
WIDTH_SWEEP = [0.12, 0.17]

solver = pybamm.IDAKLUSolver(root_tol=1e-06, atol=1e-06, rtol=1e-06)

formation_exp = pybamm.Experiment(
    [
        "Discharge at 0.1C until 2.5 V",
        "Charge at C/3 until 4.2 V",
        "Hold at 4.2 V until C/100",
    ]
)
ageing_cycle = (
    "Discharge at C/3 until 2.5 V",
    "Charge at C/3 until 4.2 V",
    "Hold at 4.2 V until C/100",
)
ageing_batch_exp = pybamm.Experiment([ageing_cycle for _ in range(BATCH_SIZE)])


def build_param_updates(width, f0):
    updates = {
        "Primary: Negative electrode active material volume fraction": GR_EPS,
        "Secondary: Negative electrode active material volume fraction": SI_EPS,
        "Secondary: Maximum concentration in negative electrode [mol.m-3]": SI_MAX_CONC_NEEDED,
        "Secondary: Initial concentration in negative electrode [mol.m-3]": scaled_si_initial_conc(SI_MAX_CONC_NEEDED),
        "Nominal cell capacity [A.h]": NOMINAL_CAP_AH,
        "Secondary: Negative electrode Young's modulus [Pa]": 1.0e10,
        "Secondary: SEI reaction exchange current density [A.m-2]": 1.5e-07 * SI_MULT,
        "Secondary: Negative electrode LAM constant proportional term [s-1]": SI_LAM_PROP_BASELINE * LAM_MULT,
        "Secondary: Negative electrode LAM constant exponential term": SI_LAM_EXP,
        "Secondary: Negative electrode critical stress [Pa]": SI_CRIT_STRESS,
        "Secondary: Negative electrode cracking rate": BASE_CRACK_RATE * SI_CRACK_RATE_MULT,
        "Primary: SEI reaction exchange current density [A.m-2]": 1.5e-07 / GR_DIV,
        "Primary: Negative electrode LAM constant exponential term": 2.0,
        "Primary: Negative electrode critical stress [Pa]": GR_CRIT_STRESS,
        "Primary: Negative electrode LAM constant proportional term [s-1]": GR_LAM_PROP_BASELINE * LAM_MULT,
        "Primary: Negative electrode cracking rate": BASE_CRACK_RATE * GR_CRACK_RATE_MULT,
        "Negative electrode porosity floor": NEG_POROSITY_FLOOR,
        "Primary: SEI reaction exponent cap": EXPONENT_MAX_SEI,
        "Secondary: SEI reaction exponent cap": EXPONENT_MAX_SEI,
        "Negative electrode pore buffering transition width": width,
        "Negative electrode transmitted fraction plateau": f0,
    }
    return updates


def cycle_ageing_leg(cyc):
    cap, discharge_end_thr, charge_end_thr = 0.0, None, None
    for step in getattr(cyc, "steps", []):
        try:
            I = step["Current [A]"].entries
            q = step["Discharge capacity [A.h]"].entries
            thr_step = step["Throughput capacity [A.h]"].entries
        except (KeyError, TypeError, AttributeError):
            continue
        if I.size < 2:
            continue
        mean_I = np.mean(I)
        if mean_I > 1.0 and q.size >= 2:
            c = float(q[-1] - q[0])
            if c > cap:
                cap = c
                discharge_end_thr = float(thr_step[-1])
        elif mean_I < -1e-3:
            charge_end_thr = float(thr_step[-1])
    return cap, discharge_end_thr, charge_end_thr


def current_soh(sol):
    """Ratio of the most recent valid ageing cycle's discharge capacity to the first
    valid cycle's -- reuses the same per-cycle extraction as the post-loop analysis
    below (rather than indexing sol.cycles[0]/[-1] directly, which returned cap=0 for
    every cycle when tried -- likely stale/unprocessed step data on a Solution object
    still being extended via starting_solution). Returns None if <2 valid cycles yet."""
    caps = []
    for cyc in sol.cycles:
        cap, discharge_end_thr, charge_end_thr = cycle_ageing_leg(cyc)
        if cap > 0 and discharge_end_thr is not None and charge_end_thr is not None:
            caps.append(cap)
    if len(caps) < 2:
        return None
    return caps[-1] / caps[0]


def run_case(label, model_options, width, f0):
    param = pybamm.ParameterValues("si_gr_expansion")
    param.update(build_param_updates(width, f0), check_already_exists=False)
    model = pybamm.lithium_ion.DFN(model_options)

    def make_sim(experiment):
        return pybamm.Simulation(
            model, parameter_values=param, experiment=experiment,
            solver=solver, var_pts=VAR_PTS,
        )

    tag = f"[{label}]"
    sim = make_sim(formation_exp)
    sol = sim.solve(initial_soc=1.0)
    print(f"{tag} formation solved.")

    total_cycles = 0
    stop_reason = "reached MAX_TOTAL_CYCLES safety cap"
    while total_cycles < MAX_TOTAL_CYCLES:
        sim = make_sim(ageing_batch_exp)
        try:
            sol = sim.solve(starting_solution=sol)
        except Exception as exc:  # noqa: BLE001
            print(f"{tag} solver failure after {total_cycles} cycles: "
                  f"{type(exc).__name__}: {str(exc)[:200]}")
            stop_reason = "solver_failure"
            break
        total_cycles += BATCH_SIZE
        soh_now = current_soh(sol)
        print(f"{tag} {total_cycles} cycles solved. "
              f"SoH={'n/a' if soh_now is None else f'{soh_now:.3f}'}")
        if soh_now is not None and soh_now <= SOH_TARGET:
            stop_reason = f"reached SOH_TARGET ({SOH_TARGET:.0%})"
            break

    print(f"{tag} stop reason: {stop_reason}")

    Qt_full = sol["Throughput capacity [A.h]"].entries
    tc_cell_full = sol["Cell thickness change [m]"].entries

    thr_list, cap_list, discharge_end_list, charge_end_list = [], [], [], []
    for cyc in sol.cycles:
        cap, discharge_end_thr, charge_end_thr = cycle_ageing_leg(cyc)
        if cap > 0 and discharge_end_thr is not None and charge_end_thr is not None:
            thr_list.append(float(cyc["Throughput capacity [A.h]"].entries[-1]))
            cap_list.append(cap)
            discharge_end_list.append(discharge_end_thr)
            charge_end_list.append(charge_end_thr)

    thr_arr = np.array(thr_list)
    cap_arr = np.array(cap_list)
    soh_arr = cap_arr / cap_arr[0] if len(cap_arr) else cap_arr
    discharge_end_arr = np.array(discharge_end_list)
    charge_end_arr = np.array(charge_end_list)

    tc_charge_end = np.interp(charge_end_arr, Qt_full, tc_cell_full)
    tc_discharge_end = np.interp(discharge_end_arr, Qt_full, tc_cell_full)
    cell_amplitude = tc_charge_end - tc_discharge_end

    if model_options.get("pore buffering") == "true":
        k_full = sol["Negative electrode transfer ratio k"].entries
        k_at_charge_end = np.interp(charge_end_arr, Qt_full, k_full)
    else:
        k_at_charge_end = None

    LLI_final = float(sol["Loss of lithium inventory [%]"].entries[-1])
    LAM_neg_final = float(sol["Loss of active material in negative electrode [%]"].entries[-1])

    n_half = max(3, len(thr_arr) // 2)
    if n_half >= 2:
        slope, _ = np.polyfit(thr_arr[:n_half], cell_amplitude[:n_half] * 1e6, 1)
    else:
        slope = float("nan")

    amp_um = cell_amplitude * 1e6
    peak_idx = int(np.argmax(amp_um))
    pre_peak_min = float(np.min(amp_um[: peak_idx + 1])) if peak_idx > 0 else float(amp_um[0])
    hump_rise = float(amp_um[peak_idx] - pre_peak_min)

    print(f"{tag} final: LLI={LLI_final:.2f}%  LAM_neg={LAM_neg_final:.2f}%  "
          f"early-life amplitude slope={slope:.5f} um/A.h "
          f"({'DECLINING' if slope < 0 else 'INCREASING'})  hump_rise={hump_rise:.2f} um  "
          f"BoL amplitude={amp_um[0]:.2f} um  final SoH={soh_arr[-1]:.3f}  "
          f"cycles={total_cycles}  stop_reason={stop_reason}")

    return dict(
        label=label, width=width, f0=f0,
        thr=thr_arr, cap=cap_arr, soh=soh_arr,
        cell_amplitude=cell_amplitude, k=k_at_charge_end,
        LLI_final=LLI_final, LAM_neg_final=LAM_neg_final, early_slope=slope,
        hump_rise=hump_rise, bol_amplitude=amp_um[0],
    )


print("\n=== Unbuffered reference (same recipe, pore buffering off) ===")
unbuffered_ref = run_case("Unbuffered reference", UNBUFFERED_OPTIONS, WIDTH_SWEEP[0], F0_SWEEP[0])

results = []
for width in WIDTH_SWEEP:
    for f0 in F0_SWEEP:
        print(f"\n=== width={width:g}, f0={f0:g} ===")
        results.append(run_case(f"w={width:g},f0={f0:g}", MODEL_OPTIONS, width, f0))

# Measured (observable) transfer ratio: buffered amplitude / same-recipe unbuffered
# amplitude at matching throughput -- distinct from the "internal" k variable (the raw
# model f(eps_struct) at each instant). Buffering doesn't alter degradation kinetics
# (LLI/LAM essentially unchanged across the whole grid, confirmed in A.0l), so the
# unbuffered reference's own throughput range safely covers every buffered case's range.
for r in results:
    unbuf_amp_interp = np.interp(r["thr"], unbuffered_ref["thr"], unbuffered_ref["cell_amplitude"])
    r["measured_k"] = r["cell_amplitude"] / (unbuf_amp_interp + 1e-30)

# ---------------------------------------------------------------------------
# plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(2, 2, figsize=(14, 10))

width_styles = {WIDTH_SWEEP[0]: "-", WIDTH_SWEEP[1]: "--"}
cmap_f0 = {f0: c for f0, c in zip(F0_SWEEP, plt.cm.viridis(np.linspace(0.1, 0.9, len(F0_SWEEP))))}

for r in results:
    ax[0, 0].plot(r["thr"], r["cell_amplitude"] * 1e6, width_styles[r["width"]], lw=2,
                  color=cmap_f0[r["f0"]],
                  label=f"w={r['width']:g},f0={r['f0']:g} (slope={r['early_slope']:.4f})")
ax[0, 0].plot(unbuffered_ref["thr"], unbuffered_ref["cell_amplitude"] * 1e6, "k:",
              lw=1.5, label=f"Unbuffered ref (slope={unbuffered_ref['early_slope']:.4f})")
ax[0, 0].set_title("Cell-level amplitude (solid=w0.12, dashed=w0.17)")
ax[0, 0].set_xlabel("Throughput capacity [A.h]")
ax[0, 0].set_ylabel("Cell-level within-cycle amplitude [um]")
ax[0, 0].legend(fontsize=6, ncol=2)
ax[0, 0].grid(alpha=0.3)

for r in results:
    ax[0, 1].plot(r["thr"], r["k"], width_styles[r["width"]], lw=2, color=cmap_f0[r["f0"]],
                  label=f"w={r['width']:g},f0={r['f0']:g}")
ax[0, 1].set_title("Internal transfer ratio k (= f at charge-end)")
ax[0, 1].set_xlabel("Throughput capacity [A.h]")
ax[0, 1].set_ylabel("k [-]")
ax[0, 1].legend(fontsize=6, ncol=2)
ax[0, 1].grid(alpha=0.3)

for r in results:
    ax[1, 0].plot(r["thr"], r["soh"] * 100, width_styles[r["width"]], lw=2,
                  color=cmap_f0[r["f0"]], label=f"w={r['width']:g},f0={r['f0']:g}")
ax[1, 0].plot(unbuffered_ref["thr"], unbuffered_ref["soh"] * 100, "k:", lw=1.5,
              label="Unbuffered ref (same degradation kinetics)")
ax[1, 0].axhline(50, color="gray", ls="--", lw=1)
ax[1, 0].set_xlabel("Throughput capacity [A.h]")
ax[1, 0].set_ylabel("SoH (discharge capacity) [%]")
ax[1, 0].set_title("State of Health vs. throughput")
ax[1, 0].legend(fontsize=6, ncol=2)
ax[1, 0].grid(alpha=0.3)

for r in results:
    ax[1, 1].plot(r["thr"], r["measured_k"], width_styles[r["width"]], lw=2,
                  color=cmap_f0[r["f0"]], label=f"w={r['width']:g},f0={r['f0']:g}")
ax[1, 1].set_xlabel("Throughput capacity [A.h]")
ax[1, 1].set_ylabel("Measured k = buffered amplitude / unbuffered amplitude [-]")
ax[1, 1].set_title("Measured (observable) transfer ratio -- not internal f(eps_struct)")
ax[1, 1].legend(fontsize=6, ncol=2)
ax[1, 1].grid(alpha=0.3)

fig.suptitle(
    f"'physical' option: f0 x eps_transfer_width grid (base-case recipe SI_MULT={SI_MULT:g}, "
    f"run to {SOH_TARGET:.0%} SoH)"
)
plt.tight_layout()
outpath = os.path.join(SCRIPT_DIR, "physical_f0_width_grid_test_result.png")
plt.savefig(outpath, dpi=150)
print(f"\nSaved: {outpath}")

print("\n--- summary ---")
for r in results:
    print(f"{r['label']}: early-life slope={r['early_slope']:.5f} um/A.h  "
          f"({'DECLINING' if r['early_slope'] < 0 else 'INCREASING'})  "
          f"hump_rise={r['hump_rise']:.2f} um  BoL amplitude={r['bol_amplitude']:.2f} um")
print(f"Unbuffered reference: slope={unbuffered_ref['early_slope']:.5f} um/A.h "
      f"({'DECLINING' if unbuffered_ref['early_slope'] < 0 else 'INCREASING'})")
any_positive = any(r["early_slope"] > 0 for r in results)
print(f"\nAny (width, f0) combination produced a POSITIVE early-life slope? {any_positive}")
