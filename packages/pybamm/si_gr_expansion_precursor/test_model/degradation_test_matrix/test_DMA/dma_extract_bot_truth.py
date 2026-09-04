"""
dma_extract_bot_truth.py -- one-off helper: runs ONLY formation + the first ageing batch
(reaching RPT0, identical recipe/protocol to dma_baseline_run.py, but far cheaper than the
full ~15-RPT run) and reads off the TRUE electrode stoichiometry at the start (fully
charged, full-cell SOC=1) and end (fully discharged, SOC=0) of RPT0's own C/20 discharge
leg, directly from PyBaMM's own particle-level state.

Used to anchor dma_method_b_plot.py's RPT0 fit to ground truth (see that script's
docstring for why): RPT0 has essentially zero real degradation, so its true electrode
windows are not actually ambiguous, but the DMA fit's own (x_ne_lo, x_ne_hi, x_pe_lo,
x_pe_hi, xi) had been landing in a self-consistent-but-wrong basin regardless (LAM_pos
absorbing degradation that is really the anode's). This script produces the numbers to
fix that at the one point (RPT0) where we can directly check them against simulation
ground truth.

x_pe_lo/x_pe_hi map directly onto "X-averaged positive particle surface stoichiometry"
at RPT0's discharge end/start respectively. x_ne_lo/x_ne_hi need one extra step: PyBaMM
tracks graphite and silicon's OWN stoichiometries separately (there is no single "blend
SOC" state), so the equivalent blend-model quantity is computed the same way
dma_ocp_fit.py's reconstruct_pocv would define it: z_ne_blend = xi*z_Si + (1-xi)*z_Gr,
evaluated with this recipe's OWN known xi (~0.45, from CAP_GR_FRAC=0.55 in
dma_baseline_run.py).
"""
import numpy as np
import pybamm

pybamm.set_logging_level("NOTICE")

MODEL_OPTIONS_BASE = {
    "particle phases": ("2", "1"),
    "open-circuit potential": (("single", "current sigmoid"), "single"),
    "SEI": "ec reaction limited",
    "SEI porosity change": "true",
    "particle mechanics": ("swelling and cracking", "swelling only"),
    "SEI on cracks": "true",
    "loss of active material": "stress-driven",
    "pore buffering": "true",
    "pore buffering transition": "physical",
}
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


print("Measuring capacity constants...")
k_gr_orig, k_si_orig = measure_capacity_constants(SI_MAX_CONC_DEFAULT)
k_si_orig_global = k_si_orig
CAP_GR_FRAC = 0.55
SI_VOL_FRAC = 0.20
GR_EPS, SI_EPS, SI_MAX_CONC_NEEDED = solve_capacity_split_at_fixed_volume(
    SI_VOL_FRAC, CAP_GR_FRAC, k_gr_orig, NOMINAL_CAP_AH)
XI_TRUE = 1 - CAP_GR_FRAC  # ~0.45 -- this recipe's known Si capacity share

SI_CRIT_STRESS = 2.2e8
SI_LAM_PROP_BASELINE = 7.5e-7
SI_LAM_EXP = 2.5
GR_CRIT_STRESS = 3.0e7
GR_LAM_PROP_BASELINE = 4.5e-6
SI_CRACK_RATE_MULT = 0.1
GR_CRACK_RATE_MULT = 0.1
BASE_CRACK_RATE = 3.9e-20
NEG_POROSITY_FLOOR = 0.01
EXPONENT_MAX_SEI = 10.0
F0_BASELINE = 0.7
WIDTH_BASELINE = 0.01
K_SEI_DEFAULT = 1e-12
K_SEI_ANCHOR_MULT = 0.0017

param = pybamm.ParameterValues("si_gr_expansion")
param.update(
    {
        "Primary: Negative electrode active material volume fraction": GR_EPS,
        "Secondary: Negative electrode active material volume fraction": SI_EPS,
        "Secondary: Maximum concentration in negative electrode [mol.m-3]": SI_MAX_CONC_NEEDED,
        "Secondary: Initial concentration in negative electrode [mol.m-3]": scaled_si_initial_conc(SI_MAX_CONC_NEEDED),
        "Nominal cell capacity [A.h]": NOMINAL_CAP_AH,
        "Secondary: Negative electrode Young's modulus [Pa]": 1.0e10,
        "Secondary: Negative electrode LAM constant proportional term [s-1]": SI_LAM_PROP_BASELINE,
        "Secondary: Negative electrode LAM constant exponential term": SI_LAM_EXP,
        "Secondary: Negative electrode critical stress [Pa]": SI_CRIT_STRESS,
        "Secondary: Negative electrode cracking rate": BASE_CRACK_RATE * SI_CRACK_RATE_MULT,
        "Primary: Negative electrode LAM constant exponential term": 2.0,
        "Primary: Negative electrode critical stress [Pa]": GR_CRIT_STRESS,
        "Primary: Negative electrode LAM constant proportional term [s-1]": GR_LAM_PROP_BASELINE,
        "Primary: Negative electrode cracking rate": BASE_CRACK_RATE * GR_CRACK_RATE_MULT,
        "Negative electrode porosity floor": NEG_POROSITY_FLOOR,
        "Primary: SEI reaction exponent cap": EXPONENT_MAX_SEI,
        "Secondary: SEI reaction exponent cap": EXPONENT_MAX_SEI,
        "Negative electrode pore buffering transition width": WIDTH_BASELINE,
        "Negative electrode transmitted fraction plateau": F0_BASELINE,
        "Primary: SEI kinetic rate constant [m.s-1]": K_SEI_DEFAULT * K_SEI_ANCHOR_MULT,
        "Secondary: SEI kinetic rate constant [m.s-1]": K_SEI_DEFAULT * K_SEI_ANCHOR_MULT,
    },
    check_already_exists=False,
)

solver = pybamm.IDAKLUSolver(root_tol=1e-06, atol=1e-06, rtol=1e-06)
model = pybamm.lithium_ion.DFN(MODEL_OPTIONS_BASE)

formation_exp = pybamm.Experiment(
    [
        "Discharge at 0.1C until 2.5 V",
        "Charge at C/3 until 4.2 V",
        "Hold at 4.2 V until C/100",
    ]
)
# Matches dma_baseline_run.py's batch 1 (RPT_INTERVAL=10, position 10 = RPT leg). Each
# element must be a TUPLE (one 3-step cycle), not flattened -- a flat list of strings
# makes pybamm.Experiment treat every individual STEP as its own 1-step "cycle" instead.
_ageing_cycle = ("Discharge at C/3 until 2.5 V", "Charge at C/3 until 4.2 V", "Hold at 4.2 V until C/100")
_rpt_cycle = ("Discharge at C/20 until 2.5 V", "Charge at C/3 until 4.2 V", "Hold at 4.2 V until C/100")
rpt0_batch_exp = pybamm.Experiment([_ageing_cycle] * 9 + [_rpt_cycle])

sim0 = pybamm.Simulation(model, parameter_values=param, experiment=formation_exp,
                          solver=solver, var_pts=VAR_PTS)
sol0 = sim0.solve(initial_soc=1.0)
print("formation solved.")

sim1 = pybamm.Simulation(model, parameter_values=param, experiment=rpt0_batch_exp,
                          solver=solver, var_pts=VAR_PTS)
sol1 = sim1.solve(starting_solution=sol0)
print("batch 1 (containing RPT0) solved.")

# Find the RPT0 cycle (last cycle in this batch) and its charge-end (SOC=1, RPT start)
# and discharge-end (SOC=0, RPT end) timepoints.
rpt_cycle = sol1.cycles[-1]
gr_stoich_full = sol1["X-averaged negative primary particle surface stoichiometry"].entries
si_stoich_full = sol1["X-averaged negative secondary particle surface stoichiometry"].entries
pe_stoich_full = sol1["X-averaged positive particle surface stoichiometry"].entries
t_full = sol1["Time [s]"].entries

for step in rpt_cycle.steps:
    I = step["Current [A]"].entries
    t_step = step["Time [s]"].entries
    if I.size < 2:
        continue
    mean_I = np.mean(I)
    if 0.01 < mean_I <= 1.0:  # the C/20 RPT discharge leg
        t_start, t_end = t_step[0], t_step[-1]
        gr_start = np.interp(t_start, t_full, gr_stoich_full)
        gr_end = np.interp(t_end, t_full, gr_stoich_full)
        si_start = np.interp(t_start, t_full, si_stoich_full)
        si_end = np.interp(t_end, t_full, si_stoich_full)
        pe_start = np.interp(t_start, t_full, pe_stoich_full)
        pe_end = np.interp(t_end, t_full, pe_stoich_full)
        break
else:
    raise RuntimeError("RPT0 discharge leg not found in this batch")

x_ne_hi_true = XI_TRUE * si_start + (1 - XI_TRUE) * gr_start  # SOC=1 (charged)
x_ne_lo_true = XI_TRUE * si_end + (1 - XI_TRUE) * gr_end        # SOC=0 (discharged)
x_pe_hi_true = float(pe_start)
x_pe_lo_true = float(pe_end)

print(f"\nTrue Si capacity share xi = {XI_TRUE:.4f}")
print(f"Gr stoichiometry: start(SOC=1)={gr_start:.4f}  end(SOC=0)={gr_end:.4f}")
print(f"Si stoichiometry: start(SOC=1)={si_start:.4f}  end(SOC=0)={si_end:.4f}")
print(f"PE stoichiometry: start(SOC=1)={pe_start:.4f}  end(SOC=0)={pe_end:.4f}")
print(f"\nImplied blend-model ground truth for RPT0:")
print(f"  x_ne_lo (blend, SOC=0) = {x_ne_lo_true:.4f}")
print(f"  x_ne_hi (blend, SOC=1) = {x_ne_hi_true:.4f}")
print(f"  x_pe_lo (SOC=0)        = {x_pe_lo_true:.4f}")
print(f"  x_pe_hi (SOC=1)        = {x_pe_hi_true:.4f}")
print(f"  xi                     = {XI_TRUE:.4f}")
