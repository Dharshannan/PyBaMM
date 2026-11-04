"""Regenerate the model-dependent figures of main.tex and supplementary.tex and copy them into the Overleaf clone.

    python tools/refresh_paper_figures.py --results DIR [--ablations DIR] [--sei-variants]
                                          [--candidate refit.json] [--install] [--list]      (from paper/)

Every \\includegraphics of main.tex and supplementary.tex is resolved against its \\graphicspath and given a
source:
    model         sigr.figures.make_all from --results (overall summaries, schematics, BOL fit, RPT selection,
                  crack comparison, partial-window Discussion figure, C/20 RPT V and dV/dSoC, crack diagnostics)
    ablations     the ablation figure, only from --ablations runs made with the same parameter set
    sei_variants  the single-law SEI sweeps (sei_variants/figures), only with --sei-variants after a rerun
    data          data and DMA only (data overview, Supplementary DMA fits): unchanged by a refit
The figures go to a staging folder (--staging); --install copies each regenerated PDF over its Overleaf copy
(manuscript_overleaf/figures/...), and --paper-pics also over paper_pics/ (the repository's copy).

--candidate: provisional mode for a refit candidate not yet adopted in sigr/ (refit_local/*.json). It applies
the pending adoption changes at run time: SoH normalised C/20 RPT to C/20 RPT (as the refit scores), and the
candidate's initial and floor porosity in the schematic and crack figures. Once the set is adopted in
sigr/parameters.py and sigr/analysis.py, run without it.

--run-ablations (with --candidate): run the ablations of sigr.ablations for the candidate through
refit_local/refit_run.py (results in refit_local/runs/<tag>_abl/<ablation>/), then use them.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
HERE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "dma"))
OVERLEAF = os.path.join(HERE, "manuscript_overleaf")
TEX = ("main.tex", "supplementary.tex")
DATA_ONLY = re.compile(r"^(data_overview|si_dma_CELL\d+)$")
SEI_VARIANTS = re.compile(r"^si_(diffusion_both|kinetic_both|low_d_.*)$")


def inventory():
    """[(tex, line, name, path relative to manuscript_overleaf/figures/paper_pics, source)]."""
    out = []
    for tex in TEX:
        text = open(os.path.join(OVERLEAF, tex), encoding="utf-8").read()
        m = re.search(r"\\graphicspath\{((?:\{[^}]*\})+)\}", text)
        paths = re.findall(r"\{([^}]*)\}", m.group(1)) if m else [""]
        for i, line in enumerate(text.splitlines(), 1):
            if line.lstrip().startswith("%"):
                continue
            for name in re.findall(r"\\includegraphics(?:\[[^]]*\])?\{([^}]*)\}", line):
                found = next((os.path.join(p, name + ".pdf") for p in paths
                              if os.path.exists(os.path.join(OVERLEAF, p, name + ".pdf"))), None)
                rel = os.path.relpath(found, "figures/paper_pics").replace("\\", "/") if found else None
                src = ("data" if DATA_ONLY.match(name) else "ablations" if name == "ablations"
                       else "sei_variants" if SEI_VARIANTS.match(name) else "model")
                out.append((tex, i, name, rel, src))
    return out


def sei_variant_file(name):
    """Staged path of a Supplementary SEI-variant figure in sei_variants/figures."""
    base = os.path.join(HERE, "sei_variants", "figures")
    if name.startswith("si_low_d_"):
        return os.path.join(base, "low_d", name[len("si_low_d_"):] + ".pdf")
    return os.path.join(base, name + ".pdf")


def provisional(cand, results_dir, data_dir):
    """Apply the not-yet-adopted changes of a refit candidate to the plotting code, and simulate its
    beginning-of-life curves (bol_fit; the refit runs do not cache them) with its geometry and parameters."""
    from sigr import analysis as an, figures, CONDITIONS
    p = cand["params"]
    sys.path.insert(0, os.path.join(HERE, "refit_local"))
    import refit_run
    sim = refit_run._setup(cand.get("geometry", "paper"), cand.get("raw"))
    sim.V7 = {**sim.V7, **p}  # bol_curves reads the module's V7
    import dataclasses
    key = "25deg_full_soc"
    cond = dataclasses.replace(CONDITIONS[key], **cand.get("cond", {}).get(key, {}))
    c = sim.bol_curves(cond, data_dir, rate_basis_ah=cand.get("rate_basis_ah", sim.RATE_BASIS_AH))
    import numpy as np
    np.savez(os.path.join(results_dir, f"{key}_bol_curves.npz"),
             **{f"{k}_{x}": v for k, (q, vv) in c.items() for x, v in (("q", q), ("v", vv))})
    figures.EPS0, figures.EPS_MIN = p.get("neg_porosity", figures.EPS0), p.get("neg_porosity_floor", figures.EPS_MIN)

    def soh_reference(res):  # C/20 RPT to C/20 RPT: the first C/20 RPT (partial windows: the conditioned one)
        i = int(res["soh_ref_index"])
        return res["rpt_cap"][i if i >= 0 else 0]

    an.soh_reference = soh_reference


def run_ablations(cand_path):
    """Ablations of the candidate through refit_run.py; returns the folder holding <ablation>/<key>_results.npz."""
    import sigr.ablations as ab
    cand = json.load(open(cand_path, encoding="utf-8"))
    root = os.path.join(HERE, "refit_local", "runs", cand["tag"] + "_abl")
    procs = []
    for name, a in ab.ABLATIONS.items():
        if all(os.path.exists(os.path.join(root, name, f"{k}_results.npz")) for k in a["keys"]):
            continue
        c = dict(cand, tag=f"{cand['tag']}_abl/{name}", params={**cand["params"], **a.get("params", {})})
        p = os.path.join(HERE, "refit_local", "cmaes", f"{cand['tag']}_abl_{name}.json")
        os.makedirs(os.path.dirname(p), exist_ok=True)
        json.dump(c, open(p, "w"), indent=2)
        procs.append(subprocess.Popen([sys.executable, "-W", "ignore", os.path.join(HERE, "refit_local", "refit_run.py"),
                                       p, *a["keys"]], cwd=HERE))
    for pr in procs:
        pr.wait()
    return root


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    ap.add_argument("--ablations", default=None, help="ablation runs made with the same parameter set")
    ap.add_argument("--sei-variants", action="store_true", help="take the (rerun) sei_variants/figures")
    ap.add_argument("--candidate", default=None, help="provisional: refit candidate json not yet adopted")
    ap.add_argument("--run-ablations", action="store_true")
    ap.add_argument("--data", default=os.path.join(HERE, "data"))
    ap.add_argument("--staging", default=os.path.join(HERE, "outputs", "paper_pics_refresh"))
    ap.add_argument("--install", action="store_true", help="copy the regenerated PDFs into the Overleaf clone")
    ap.add_argument("--paper-pics", action="store_true", help="also copy PNG + PDF into paper_pics/")
    ap.add_argument("--list", action="store_true", help="only list the figures and their sources")
    a = ap.parse_args()
    os.chdir(HERE)

    inv = inventory()
    if not a.list:
        import sigr.figures
        import tables
        if a.candidate:
            cand = json.load(open(a.candidate, encoding="utf-8"))
            provisional(cand, a.results, a.data)
            if a.run_ablations:
                a.ablations = run_ablations(a.candidate)
        tables.write_tables(tables.load_fits("ours", a.data), a.data, "ours")
        sigr.figures.make_all(a.data, a.results, a.ablations, a.staging, log=lambda *x: print(*x, flush=True))

    status = []
    for tex, line, name, rel, src in inv:
        if src == "model" or (src == "ablations" and a.ablations):
            new = os.path.join(a.staging, rel) if rel else None
        elif src == "sei_variants" and a.sei_variants:
            new = sei_variant_file(name)
        else:
            new = None
        if a.list:
            state = src
        elif src == "data":
            state = "unchanged (data only)"
        elif new is None:
            state = f"STALE: needs {src} rerun with the new set"
        elif not os.path.exists(new):
            state = "MISSING from staging"
        else:
            state = "regenerated"
            if a.install:
                shutil.copy2(new, os.path.join(OVERLEAF, "figures", "paper_pics", rel))
                state = "replaced in Overleaf"
            if a.paper_pics and src != "sei_variants":
                for ext in (".pdf", ".png"):
                    s = os.path.splitext(new)[0] + ext
                    if os.path.exists(s):
                        shutil.copy2(s, os.path.join(HERE, "paper_pics", os.path.splitext(rel)[0] + ext))
        status.append((tex, line, name, rel, state))
    w = max(len(s[2]) for s in status)
    for tex, line, name, rel, state in status:
        print(f"{tex}:{line:<5d} {name:<{w}s} {state}")


if __name__ == "__main__":
    main()
