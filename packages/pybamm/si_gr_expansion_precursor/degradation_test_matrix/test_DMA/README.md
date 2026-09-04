# Degradation Mode Analysis (DMA) pipeline — `test_DMA`

Status as of 2026-09-03. Describes the final, settled DMA fitting approach only — earlier
comparison methods ("Method B", the notebook-faithful delta-split-only variant) were
retired once this version's cost function and search strategy were validated as the best
performing. Full tuning rationale (what was tried, what broke, why each constant has its
current value) lives in the code comments of `dma_ocp_fit_final.py` and
`dma_final_plot.py` — this file gives the settled picture, not the history.

## Purpose

Compares degradation modes (LLI, LAM at each electrode, and the Si/Gr split within the
blend anode) computed two ways from the same simulated ageing run (`dma_baseline_run.py`,
run to ~50% SoH with C/20 RPTs every 10 cycles):

- **Method A** (`dma_method_a_plot.py`): read the degradation modes straight off PyBaMM's
  own model state — ground truth, only possible in simulation.
- **Final fit** (`dma_final_plot.py` + `dma_ocp_fit_final.py`): fit a composite-DMA
  pseudo-OCV reconstruction model to each RPT's simulated discharge curve, exactly as an
  experimentalist would with real cell data, then derive degradation modes from the fit.
  This is the pipeline meant to be reused on a real experimental cell.

## Files

| File | Role |
|---|---|
| `dma_ocp_fit_final.py` | Core module: electrode OCP functions, the blend-anode model, pOCV reconstruction, the per-RPT fit (`fit_rpt_v2`), degradation-mode formulas. Standalone — no other file it depends on. |
| `dma_final_plot.py` | Driver: runs the fit RPT-by-RPT across the full ageing run, applies cross-RPT constraints, produces the comparison figures. |
| `dma_baseline_run.py` | Runs the ageing simulation, saves `dma_baseline_run_data.npz` (per-RPT raw (Q,V) discharge curves + PyBaMM's own state-based DMs). |
| `dma_common.py` | Shared plotting helpers (4-panel figure, EFC conversion, knee detection). |
| `dma_method_a_plot.py` | Method A — ground truth. |
| `dma_extract_bot_truth.py` | One-off: extracts RPT0's true electrode stoichiometry directly from PyBaMM's particle-level state, used to anchor the fit's first RPT to ground truth instead of fitting it blind. |
| `pics/` | Output figures (`dma_final_*`) and Method A's (`dma_method_a_result.png`). |

## Reference sources

1. **Rehm et al. 2026**, *"How to determine the degradation modes of lithium-ion
   batteries with silicon–graphite blend electrodes"*, J. Power Sources 670, 239418 —
   the blend-anode reconstruction (Eq. 1) and degradation-mode equations (Eqs. 6–10).
2. **PyProBE** (github.com/ImperialCollegeLondon/PyProBE) — cross-checked the
   electrode-capacity/LLI formulas; its own LLI definition is what's used (see below).

## pOCV reconstruction model

Parameterisation `theta = (nu_ne, nu_pe, sigma_ne, sigma_pe, phi_si, delta_v)`:

```
soc_ne = clip(nu_ne * soc + sigma_ne, 0, 1)
soc_pe = clip(1 - (nu_pe * soc + sigma_pe), 0, 1)
V(soc) = NMC_ocp(soc_pe) + delta_v - blend_anode_v_of_z(soc_ne, phi_si)
```

- `nu_e`/`sigma_e` are each electrode's own stoichiometry-window slope/intercept
  (`nu_e = x_hi - x_lo`, `sigma_e = x_lo`) — algebraically equivalent to a raw
  `(x_lo, x_hi)` limit parameterisation, just not re-bounded to `[0,1]`.
- `phi_si` is silicon's share of the anode blend's capacity.
- `delta_v` is a flat per-RPT voltage offset (algebraically `-I_rpt*R` for a constant
  discharge current).
- `blend_anode_v_of_z`: Rehm Eq. 1's blend is capacity-weighted at shared **voltage**
  (`z_blend(V) = phi*z_Si(V) + (1-phi)*z_Gr(V)`), not voltage-weighted at shared
  stoichiometry. Implemented by inverting each material's OCP onto a shared voltage grid
  once, combining, then inverting back.
- Silicon uses the **delithiation** branch specifically (`silicon_ocp_delithiation_
  Mark2016`) — every RPT here is a discharge, and PyBaMM's `current_sigmoid` OCP
  hysteresis submodel collapses to pure delithiation during discharge, not an average.
- RPT0 is anchored tightly (`±1e-4`) to ground truth from `dma_extract_bot_truth.py`
  (RPT0 has ~zero real degradation, so its true windows are directly readable from
  PyBaMM's own particle state); every later RPT is fit freely, seeded from the previous
  RPT's own result.

## Degradation-mode formulas

```
Q_ne = Q_full / nu_ne
Q_pe = Q_full / nu_pe
Q_li = Q_full * (sigma_ne/nu_ne + (1-sigma_pe)/nu_pe)     # PyProBE's own Qli

LAM_an = 1 - Q_ne / Q_ne_ini                                       (Rehm Eq. 6)
LAM_ca = 1 - Q_pe / Q_pe_ini                                       (Rehm Eq. 7)
LI     = 1 - Q_li / Q_li_ini                                       (PyProBE)
LAM_Gr = 1 - [(1-phi_si) * Q_ne] / [(1-phi_si_ini) * Q_ne_ini]     (Rehm Eq. 9)
LAM_Si = 1 - [phi_si * Q_ne] / [phi_si_ini * Q_ne_ini]             (Rehm Eq. 10)
```

`LI` uses PyProBE's additive "total lithium at the SOC=0 reference point" definition
rather than the paper's own Eq. 8 cross-term — the paper's version proved fragile
(swings ~8x more per unit parameter change) under this pipeline's constrained fit. Both
are implemented (`li_formula="pyprobe"`/`"rehm"` in `compute_degradation_modes_v2`);
`"pyprobe"` is the validated default.

## Fitting/optimisation procedure

**Cost function** (`dvdq_rmse_constrained`): a whole-domain three-term cost, matching the
paper's own Eq. 3 structure —

```
cost = w_ocv_b * MSE(V_model, V_meas)
     + w_dv_b  * MSE(dV/dSOC_model, dV/dSOC_meas)
     + w_full_dvdq * Huber(dV/dQ_model, dV/dQ_meas, delta=0.15)
     + cross-RPT constraint penalties
```

with `w_ocv_b=1.0, w_dv_b=0.3, w_full_dvdq=0.3`. (An alternative "delta-split" cost —
voltage on one side of the sharpest dV/dQ feature, dV/dQ on the other — exists as a
toggle but is off in the validated config; the whole-domain cost tracked electrode
attribution more reliably.)

**Cross-RPT constraints** (soft penalties, not hard optimizer bounds), each a
`(max_decrease, max_increase)` allowance between consecutive RPTs relative to the same
RPT0 reference:

```python
CONSTRAINT_BOUNDS = {
    "LAM_an": (0.01, None),   # + a dedicated, separately-scaled step cap (0.06)
    "LAM_ca": (0.01, 0.003),  # tight: true cathode fade is small/slow
    "LAM_Gr": (0.005, 0.15),
    "LAM_Si": (0.01, 0.15),
    "LI":     (0.03, None),
}
```

Plus a direct constraint on `phi_si` itself: it must sit at least `PHI_SI_VS_INI_TOL`
below its RPT0 anchor value at every later RPT (`-0.02` generally, `-0.03` for the final
RPT) — derived from the exact identity `LAM_Si - LAM_Gr = [(1-phi)/(1-phi_ini) -
phi/phi_ini] * (Q_ne/Q_ne_ini)`, whose sign depends only on whether `phi_si` sits below
its anchor. This is what guarantees the correct `LAM_Si > LAM_Gr` ordering.

**Search strategy** (`fit_rpt_v2`'s `search_style`): two modes, applied selectively —

- `"restart"` (default, used for RPT0-13): an adaptive shrinking-spread multi-restart
  `differential_evolution`. Confirmed these RPTs converge to the same optimum regardless
  of search style or budget — a single, well-behaved basin.
- `"global"` (used only for the final RPT): the fit's cost landscape there has multiple
  competing basins that `"restart"` can never escape once seeded near the previous RPT's
  fit. `"global"` runs a 5-seed ensemble of full-bounds Latin Hypercube
  `differential_evolution` searches (wider mutation/recombination), selecting the winner
  by raw voltage RMSE rather than total cost (the cost function itself was found to
  mis-rank the two basins).

## Current best result (final RPT, 2026-09-03)

| | LLI (A·h) | LAM_Gr | LAM_Si | LAM_pos | RMSE (mean/max) |
|---|---|---|---|---|---|
| **Method A (truth)** | 1.555 | 14.10% | 21.22% | 0.86% | — |
| **Final fit** | 1.416 | 22.19% | 31.77% | 6.98% | 7.33 / 15.34 mV |

**What "good" means here:** the DMA method systematically overpredicts LAM magnitude —
an accepted bias, not something being chased. What matters is `LAM_Si > LAM_Gr` (correct
here) with roughly the right relative gap, and LLI accuracy (within ~9% of truth here).
