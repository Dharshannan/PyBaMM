"""Write every paper figure, at its printed JES size, to paper/paper_pics.

    python paper/tools/make_paper_figures.py [--results paper/results] [--ablations paper/outputs_ablations]
"""
import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")
HERE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, HERE)
import sigr.figures  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(HERE, "data"))
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    ap.add_argument("--ablations", default=os.path.join(HERE, "outputs_ablations"))
    ap.add_argument("--out", default=os.path.join(HERE, "paper_pics"))
    a = ap.parse_args()
    sigr.figures.make_all(a.data, a.results, a.ablations, a.out, log=lambda *x: print(*x, flush=True))


if __name__ == "__main__":
    main()
