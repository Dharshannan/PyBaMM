"""pyprobe-dma CELL017 RPT7 with the Si OCP fixed, refitted with the graphite fraction of the anode bounded to
0.70-0.80 (a direct bound; the PyProBE 2.6.1 fits of that RPT have 0.75-0.76). With the cathode box instead
(refit_rpt7_cpbox_v3.py) the fit shrinks the anode to the cycled window and gives 42 % graphite LAM. The anode
capacity keeps the package default bound and ends on it (1.91 Ah): the uncycled top of the anode is not
identified by the voltage. Patches
new_nodef.csv and its curves.

    python refit_rpt7_frac_v3.py
"""
import os
import sys

import pyprobe_dma as dma

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
import composite_dma as C  # noqa: E402
import pyprobe_dma_v3 as B  # noqa: E402
from refit_rpt7_cpbox_v3 import DATA, patch  # noqa: E402

GR_FRACTION = (0.70, 0.80)

B.constraints = lambda: [dma.Constraint("ne_comp1_fraction", lo=GR_FRACTION[0], hi=GR_FRACTION[1]),
                         dma.Constraint("r_params", lo=[0.0], hi=[0.3])]
hc = B.HalfCells(DATA)
vc, ic = B.curves(C.rpt_files(DATA, "CELL017")[7])
res = B.fit_best(hc, vc, ic)
row = B.summarise(res, "CELL017", 7, 1.0, 0.0)
print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()})
patch("new_nodef", row, B.fit_curve(res))
print("patched")
