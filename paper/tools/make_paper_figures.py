"""Write every paper figure, at its printed JES size, to paper/paper_pics.

    python paper/tools/make_paper_figures.py [--results paper/results] [--ablations paper/outputs_ablations]

Figures with DMA data (overall summaries, data overview, ablations) are written once per DMA, each from a
copy of the data directory with that DMA's LAM/LLI tables (dma/tables.py), to paper_pics/<folder>:
dma_ours (our pyprobe-dma fits; the planner and the paper), dma_collaborators (their expansion-free
discharge fits) and dma_collaborators_expansion (their fits with the expansion; local only, git-ignored).
The data directory itself is given our tables.
"""
import argparse
import os
import shutil
import sys

import matplotlib

matplotlib.use("Agg")
HERE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "dma"))
import sigr.figures  # noqa: E402
import tables  # noqa: E402

DMA_FOLDERS = {"dma_ours": "ours", "dma_collaborators": "collaborators",
               "dma_collaborators_expansion": "collaborators-expansion"}


def dma_data(data_dir, work_dir):
    """{folder: data directory with that DMA's LAM/LLI tables} (copies of data_dir under work_dir)."""
    out = {}
    for folder, source in DMA_FOLDERS.items():
        d = os.path.join(work_dir, folder)
        shutil.rmtree(d, ignore_errors=True)
        shutil.copytree(data_dir, d)
        tables.write_tables(tables.load_fits(source, d), d, source)
        out[folder] = d
    tables.write_tables(tables.load_fits("ours", data_dir), data_dir, "ours")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(HERE, "data"))
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    ap.add_argument("--ablations", default=os.path.join(HERE, "outputs_ablations"))
    ap.add_argument("--out", default=os.path.join(HERE, "paper_pics"))
    ap.add_argument("--work", default=os.path.join(HERE, "outputs", "dma_data"))
    a = ap.parse_args()
    sigr.figures.make_all(a.data, a.results, a.ablations, a.out, log=lambda *x: print(*x, flush=True),
                          dma=dma_data(a.data, a.work))


if __name__ == "__main__":
    main()
