# Single-law SEI variants (rejected, Sec. 3.1 and Supplementary Figs. S5-S11)

The paper (Sec. 3.1) uses one mixed kinetic-diffusion SEI law (Safari2009, in the form of Yang2017):
the kinetic term is kept on graphite and silicon uses its diffusion limit (Bonkile2024). These runs are the
evidence for rejecting the two single-law alternatives, each swept over its one free constant on all four cells
with every other parameter of the shared set unchanged:

| variant | SEI option (graphite, silicon) | swept | values |
|---|---|---|---|
| `diffusion_both` | solvent-diffusion limited, solvent-diffusion limited | `gr_sei_solvent_diffusivity` (D_sol,Gr, m2/s) | 1.25e-23 ... 2e-21 |
| `kinetic_both` | ec reaction limited, ec reaction limited | `si_sei_k` (k_SEI,Si, m/s) | 3e-17 ... 3e-15 (8 values) |

- `sweep.py <variant> [values] [cells]` runs a sweep in parallel; runs go to `runs/<variant>/<value>/`
  (git-ignored, ~200 MB).
- `figures.py` writes, with our pyprobe-dma DMA:
  - `figures/si_<variant>.png|pdf`: Supplementary Figs. S5 (diffusion) and S6 (kinetic), one row per cell
    (C/20 SoH, graphite SEI, Si SEI, LLI; SEI = particle surface + crack faces);
  - `figures/<variant>/compare_<cell>.png|pdf` and `figures/<variant>/<value>/<cell>_overall_summary.png|pdf`;
  - `figures/low_d/`: the negligible-graphite-SEI run (D_sol,Gr = 1.25e-23), Supplementary Figs. S7-S11;
  - `scores.csv`: mean |SoH gap| at the C/20 RPTs reached, mean C/20 voltage RMSE, SEI and SoH at end of test.

The Supplementary PDFs are copied to `manuscript_overleaf/figures/paper_pics/si/`.

Result (2026-10-28; kinetic sweep extended to 8 values 2026-10-29), mean |SoH gap| over the four cells: 3.4 pp for the paper runs; 3.8-16 pp for
`diffusion_both`; 6.8-68 pp for `kinetic_both`.

- **Kinetic on both:** the silicon SEI rises abruptly once the pores close, so the knees are near-vertical. No k_SEI,Si
  gives both the 25 degC and the 45 degC knee: at 1e-16 the 45 degC cell is far too slow (86 vs 56 % at the last RPT),
  at 2e-16 its knee is ~150 EFC early. Mean gaps 16.6, 10.5, 6.8, 14.5, 17.5, 28.6, 39.6, 67.6 pp (3e-17 ... 3e-15).
- **Diffusion on both:** a D_sol,Gr large enough to supply the late graphite SEI grows it from the first cycle, which
  fills the pores early and brings the knees of the long-life cells forward. The smallest D leaves almost no graphite
  SEI: the fits are close, but the late fade is under-predicted (25 degC full window: 53.7 vs 50 % at the last RPT).

Rerun 2026-11-05 with the adopted (design-geometry) set; runs in `runs/`, the previous ones in
`runs_paper_geometry/`; `figures.py refit_local/runs/ADOPT2` compares against the adopted runs. D_sol,Gr values
are the old ones x SEI_VOLUME_FACTOR (3.25e-23 ... 5.2e-21). Mean |SoH gap|: 2.8 pp for the adopted set,
3.0-15.8 pp for `diffusion_both` (low D: 2.95 pp, but V RMSE 68 vs 76 mV, late fade under-predicted: CELL064
53.2 vs 50 % at the last RPT), 7.4-66 pp for `kinetic_both` (best 1e-16; at 1e-16 CELL017 ends at 92 vs 56 %,
at 2e-16 it is at 35 % by 334 EFC).
