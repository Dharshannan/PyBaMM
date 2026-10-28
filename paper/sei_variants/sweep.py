"""The two single-law SEI variants rejected in Sec. 3.1, each swept over its free constant on all four cells:

    diffusion_both  diffusion limit on graphite too; sweep D_sol,Gr [m2 s-1]
    kinetic_both    kinetic (EC) law on silicon too; sweep k_SEI,Si [m s-1]

Runs go to sei_variants/runs/<variant>/<value>/ (git-ignored).

    python paper/sei_variants/sweep.py diffusion_both
    python paper/sei_variants/sweep.py kinetic_both 1e-16,3e-16 25deg_full_soc
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
sys.path.insert(0, PAPER)
import sigr  # noqa: E402

KEYS = ["25deg_full_soc", "45deg_full_soc", "25deg_15-95_soc", "25deg_20-80_soc"]
VARIANTS = {
    "diffusion_both": dict(sei=("solvent-diffusion limited", "solvent-diffusion limited"),
                           param="gr_sei_solvent_diffusivity", prefix="D",
                           values=[1.25e-23, 2.5e-23, 5e-23, 1e-22, 2.5e-22, 5e-22, 1e-21, 2e-21]),
    "kinetic_both": dict(sei=("ec reaction limited", "ec reaction limited"),
                         param="si_sei_k", prefix="k",
                         values=[1e-16, 3e-16, 1e-15, 3e-15]),
}


def run_dir(variant, value):
    return os.path.join(HERE, "runs", variant, f"{VARIANTS[variant]['prefix']}{value:.1e}")


if __name__ == "__main__":
    variant = sys.argv[1]
    v = VARIANTS[variant]
    values = [float(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else v["values"]
    keys = sys.argv[3].split(",") if len(sys.argv) > 3 else KEYS
    options = {"SEI": (v["sei"], "none")}
    jobs = [dict(key=key, data_dir=os.path.join(PAPER, "data"), out_dir=run_dir(variant, x),
                 rate_basis_ah=sigr.RATE_BASIS_AH, efc_basis_ah=sigr.EFC_BASIS_AH,
                 options=options, params={v["param"]: x}, label=f"{variant} {key} {v['prefix']}={x:.1e}")
            for x in values for key in keys]
    sigr._run_jobs(jobs, report_every_s=600)
