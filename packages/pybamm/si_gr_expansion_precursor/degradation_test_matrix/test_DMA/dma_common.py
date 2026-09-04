"""
dma_common.py -- shared knee-detection and 4-panel plotting helpers for
dma_method_a_plot.py / dma_method_b_plot.py, so the two "different ways to compute LLI
and LAM" plots are laid out identically (mirroring the CELL064-style experimental figure
the user shared) and only differ in where panels (c)/(d)'s LAM/LLI numbers come from.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter


def efc_from_throughput(thr, nominal_cap):
    """Equivalent full cycles = throughput / (2 * nominal capacity), this project's
    standing convention (conditions_test_matrix_plan.md)."""
    return thr / (2.0 * nominal_cap)


def find_knee_bisector(thr, soh, window=11, polyorder=3):
    """Knee estimate via maximum curvature: smooth the (thr, soh) trajectory with a
    Savitzky-Golay filter, then locate the point of steepest downward bend (most negative
    second derivative d^2(soh)/d(thr)^2). Returns (thr_knee, soh_knee) or None if too few
    points.

    REVISION HISTORY (two fixes, in order -- both found by comparing the detected knee
    against the actual dma_baseline_run data/figures, not just trusting the number):
      v1 used a FIXED (pre_frac=0.3, post_frac=0.15) two-line bisector -- always fit a
         "pre" line to the first 30% of points and a "post" line to the last 15%,
         regardless of where the real bend is. Wrong here: this recipe's capacity stays
         close to linear until ~72 EFC (of a 0-96.5 EFC run) before a sharp fast-fade
         onset, i.e. ~75% through, not 30%. The fixed windows' extrapolated lines crossed
         at ~42 EFC, visibly nowhere near the actual kink in the capacity/expansion plots.
      v2 replaced the fixed split with a best-split segmented regression (search every
         candidate breakpoint for the two-line fit with lowest combined SSE) -- an
         improvement (62.2 EFC) but still found to sit ~10 EFC too early on visual
         inspection. Root cause: the transition here isn't a sharp instantaneous corner,
         it's a smoothly ACCELERATING bend spanning ~EFC 55-73 (the fade rate ramps up
         gradually before the sharpest drop) -- a global two-line SSE fit averages over
         that whole ramp and is pulled earlier by it, structurally unable to land at the
         point of actual peak curvature.
      v3 (this version) finds that point directly via the curve's own second derivative
         instead of fitting two lines at all -- confirmed to land at ~71.5 EFC, matching
         both the single sharpest RPT-to-RPT capacity drop in the data (EFC 72.2->72.9)
         and the point where the internal transfer ratio k saturates to ~1.0, and visibly
         aligned with the bend in the four-panel plots' capacity/expansion panels."""
    n = len(soh)
    if n < 8:
        return None
    win = min(window, n - 1 if (n - 1) % 2 == 1 else n - 2)
    if win % 2 == 0:
        win -= 1
    if win < polyorder + 2:
        return None
    soh_s = savgol_filter(soh, win, polyorder)
    d1 = np.gradient(soh_s, thr)
    d2 = np.gradient(d1, thr)
    i = int(np.argmin(d2))
    return float(thr[i]), float(soh_s[i])


def make_four_panel_figure(
    outpath, suptitle,
    age_efc, age_amplitude, age_k, age_cap,
    rpt_efc, rpt_cap,
    lam_si, lam_gr, lam_pos, lli,
    lli_is_percent, nominal_cap,
):
    """Reproduces the CELL064-style 4-panel layout:
      (1) reversible expansion (dense, C/3) + k ratio (twin axis)
      (2) discharge capacity: dense C/3 "reconstruction" line + C/20 RPT markers
      (3) LAM_Si / LAM_Gr / LAM_pos vs EFC-knee
      (4) LLI vs EFC-knee (as A.h, converting from % via nominal_cap if lli_is_percent)

    All x-axes are EFC - EFC_at_knee (knee found from the dense age_cap/age_cap[0] curve).
    """
    cap0 = age_cap[0]
    knee = find_knee_bisector(age_efc, age_cap / cap0)
    efc_knee = knee[0] if knee is not None else 0.0

    lli_ah = (lli / 100.0 * nominal_cap) if lli_is_percent else lli

    fig, ax = plt.subplots(4, 1, figsize=(9, 15), sharex=True)

    ax0b = ax[0].twinx()
    ax[0].plot(age_efc - efc_knee, age_amplitude * 1e6, "o", ms=3, color="tab:green",
               label="reversible expansion")
    ax0b.plot(age_efc - efc_knee, age_k, "o-", ms=3, color="tab:orange", label="k")
    ax[0].axvline(0, color="gray", ls="--", lw=1)
    ax[0].set_ylabel("Reversible\nexpansion [um]")
    ax0b.set_ylabel("k = transfer ratio [-]", color="tab:orange")
    ax[0].set_title("(a) Reversible expansion + internal transfer ratio k")
    lines0, labels0 = ax[0].get_legend_handles_labels()
    lines0b, labels0b = ax0b.get_legend_handles_labels()
    ax[0].legend(lines0 + lines0b, labels0 + labels0b, fontsize=8, loc="best")
    ax[0].grid(alpha=0.3)

    ax[1].plot(age_efc - efc_knee, age_cap, "-", lw=1.5, color="gray",
               label="knee reconstruction (dense C/3)")
    ax[1].plot(rpt_efc - efc_knee, rpt_cap, "o", ms=6, color="tab:blue", label=f"C/20 RPT")
    ax[1].axvline(0, color="gray", ls="--", lw=1)
    ax[1].set_ylabel("Discharge\ncapacity [A.h]")
    ax[1].set_title("(b) Capacity fade")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)

    ax[2].plot(rpt_efc - efc_knee, lam_si, "o-", ms=5, color="tab:purple", label="LAM NE Si")
    ax[2].plot(rpt_efc - efc_knee, lam_gr, "o-", ms=5, color="tab:blue", label="LAM NE Gr")
    ax[2].plot(rpt_efc - efc_knee, lam_pos, "o-", ms=5, color="tab:brown", label="LAM PE")
    ax[2].axvline(0, color="gray", ls="--", lw=1)
    ax[2].set_ylabel("LAM [%]")
    ax[2].set_title("(c) Loss of active material")
    ax[2].legend(fontsize=8)
    ax[2].grid(alpha=0.3)

    ax[3].plot(rpt_efc - efc_knee, lli_ah, "o-", ms=5, color="tab:red", label="LLI loss")
    ax[3].axvline(0, color="gray", ls="--", lw=1)
    ax[3].set_xlabel("EFC - capacity knee [n]")
    ax[3].set_ylabel("LLI loss [A.h]")
    ax[3].set_title("(d) Loss of lithium inventory")
    ax[3].legend(fontsize=8)
    ax[3].grid(alpha=0.3)

    fig.suptitle(f"{suptitle}\n(knee = {efc_knee:.1f} EFC)")
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    return efc_knee
