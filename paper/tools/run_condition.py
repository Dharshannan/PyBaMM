"""Command-line run of one condition (same as the notebook's run cell).

    python paper/tools/run_condition.py 25deg_full_soc --data paper/data --out paper/outputs
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import sigr  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("key", choices=list(sigr.CONDITIONS))
    ap.add_argument("--data", default=os.path.join(os.path.dirname(__file__), "..", "data"))
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "outputs"))
    ap.add_argument("--rate-basis-ah", type=float, default=sigr.RATE_BASIS_AH)
    ap.add_argument("--efc-basis-ah", type=float, default=sigr.EFC_BASIS_AH)
    a = ap.parse_args()
    cond = sigr.CONDITIONS[a.key]
    res, curves = sigr.run_condition(a.key, a.data, a.out, a.rate_basis_ah, a.efc_basis_ah,
                                     log=lambda *x: print(*x, flush=True))
    print(sigr.analysis.summary(res, curves, cond, a.data), flush=True)
    for path in sigr.plots.save_all(res, curves, cond, a.data, a.out).values():
        print("Saved:", path, flush=True)


if __name__ == "__main__":
    main()
