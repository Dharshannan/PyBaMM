"""
dma_method_a_plot.py -- "Method A": degradation modes read DIRECTLY off PyBaMM's own
model-state variables (Loss of lithium inventory [%], Loss of active material in
{primary,secondary} phase in negative electrode [%], Loss of active material in positive
electrode [%]) -- the "ground truth" a simulation gives for free, with no curve-fitting
needed. Produces a CELL064-style 4-panel figure (dma_common.make_four_panel_figure) from
dma_baseline_run.py's saved trajectory, for qualitative comparison against
dma_method_b_plot.py's composite-DMA-fit version of the SAME run.
"""
import os

import numpy as np

from dma_common import efc_from_throughput, make_four_panel_figure

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

data = np.load(os.path.join(SCRIPT_DIR, "dma_baseline_run_data.npz"), allow_pickle=True)
nominal_cap = float(data["nominal_cap"])

age_efc = efc_from_throughput(data["age_thr"], nominal_cap)
rpt_efc = efc_from_throughput(data["rpt_thr"], nominal_cap)

outpath = os.path.join(SCRIPT_DIR, "pics", "dma_method_a_result.png")
efc_knee = make_four_panel_figure(
    outpath,
    "Method A: degradation modes from direct PyBaMM model-state variables",
    age_efc, data["age_amplitude"], data["age_k"], data["age_cap"],
    rpt_efc, data["rpt_cap"],
    lam_si=data["rpt_LAM_si"], lam_gr=data["rpt_LAM_gr"], lam_pos=data["rpt_LAM_pos"],
    lli=data["rpt_LLI"], lli_is_percent=True, nominal_cap=nominal_cap,
)
print(f"Saved: {outpath}  (knee at {efc_knee:.1f} EFC)")

print("\n--- final degradation modes (Method A, direct model state) ---")
print(f"LLI={data['rpt_LLI'][-1]:.2f}%  LAM_Gr={data['rpt_LAM_gr'][-1]:.2f}%  "
      f"LAM_Si={data['rpt_LAM_si'][-1]:.2f}%  LAM_pos={data['rpt_LAM_pos'][-1]:.2f}%")
