# Expansion Precursor: Hypotheses & Test Plan

Source: slide deck **"Expansion as an Early Precursor to Capacity Knees — Demonstration
on Si/Gr anode batteries"** (Zhiwen Wan, UM Battery Control Lab).
(`C:\Users\ds3420\Downloads\Expansion_as_Precursor_Modeling_Discussion.pdf`)

Status: **planning only**, per request — no implementation in this pass.

Three separate threads, matching the request:
- **Part A** — hypotheses for *why* the pre-hump expansion trend differs across the
  operating conditions on page 2, and whether the current model could capture each.
  **Part A.0 is a correction worth reading first, and Part A.0g (below it) is a second,
  more important correction that supersedes A.0's own "majority pattern" framing** — an
  updated, more complete experimental figure (5 conditions) shows **decline-then-hump is
  actually the majority pattern**, matching the model's existing default behaviour for
  baseline, high C-rate, high temperature, and a narrow-SoC-window condition; **only low
  pressure is a genuine outlier**, showing a continuously increasing pre-hump trend. Read
  A.0g before treating A.0–A.0f's "fix the majority pattern via LLI-vs-LAM rebalancing" as
  the live priority — that work remains a valid, informative investigation, but the premise
  motivating it (that most conditions need an increasing trend the model doesn't produce)
  has been revised. **A.0h–A.0l are a further pivot worth reading too**: the natural
  pressure-coupling extension (a `K_stack` rescale) was implemented, tested, and found to
  make things worse, not better (A.0h); that motivated a genuinely-unmodeled-mechanism
  hypothesis (particle/stack rearrangement, A.0i) — but before committing to that bigger
  lift, a cheaper shape-only fix (an asymmetric, gradual-to-sharp transition instead of the
  existing tanh/physical forms) was tried and flips the sign (A.0j), though without an
  obviously proper hump; A.0k then showed the *existing*, unmodified "physical" form can't
  fix that alone by widening its transition either; **A.0l found the resolution**: combining
  a lower `f0` with a wider `eps_transfer_width`, both within the existing "physical" form,
  produces a genuinely increasing pre-hump trend *and* a proper hump, no new functional
  form or new mechanism needed, confirmed holding to 50% SoH — currently the best-supported,
  cheapest lead for the low-pressure outlier. **A.0m checks this directly against the real
  numbers** (an updated experimental figure, `fig2_population_lead_time_v2.png`, embedded
  in A.0g, which also adds a fifth condition, cross-condition overlay/summary panels, and a
  new unexplained combined condition) and finds a genuinely good qualitative match — the
  most positive result the whole investigation has produced.
- **Part B** — a test-matrix plan reproducing page 6's 5-row mechanism matrix, run at a
  **single** operating condition.
- **Part C** — separately, a staged plan for what it would take to actually capture the
  page-2 operating-condition-*dependent* pattern differences (pressure/rate/temperature),
  which Part B deliberately does not touch. **Per A.0g, pressure-coupling (Part C.4) is
  the top-priority item** — it directly targets the one confirmed outlier — ahead of the
  rate/temperature extensions and ahead of further LLI-vs-LAM rebalancing work, which was
  motivated by a "majority pattern" premise A.0g revises. **Per A.0h/A.0i, the concrete
  mechanism has since moved on from C.4's original `K_stack(P_stack)` proposal**: a direct
  test showed a constant-rescale extension doesn't just fail to reproduce the low-pressure
  shape, it moves the wrong way — the leading candidate is now a genuinely new,
  pressure-gated particle/stack-rearrangement mechanism (A.0i), not a new parameter on the
  existing pore-buffering partition.

---

## Part A — Why does the pre-hump trend differ across conditions? (page 2)

### A.0 Cross-check against our own model's simulated behaviour — important correction

Before treating "baseline = flat, model already reproduces it" (an earlier draft of this
plan said exactly that) as settled, it's worth checking what our own model *actually*
simulates. It doesn't hold up: looking at the already-generated
`test_pore_buffering/pore_buffering_degradation_test_1000cyc_result_physical.png` (bottom
panel, cell-level and negative-electrode-level within-cycle expansion amplitude vs.
throughput), the simulated trend **declines** from BoL (~12 µm) down to a minimum around
4,500 A.h throughput (~9.8 µm), *then* rises to the hump peak (~13.8 µm) around 6,600 A.h,
before falling post-knee. That is a **decline-then-hump** shape.

Cross-referencing against Part A.1's four experimental conditions: decline-then-hump is
the shape of exactly **one** of them — **High Temperature**, the outlier/minority
condition — not baseline, low pressure, or high C-rate, all three of which show
*increasing* trends before the hump. In other words: **our model's current default
calibration reproduces the minority experimental pattern, not the majority one.** This is
a more important finding than the original per-condition gap table below suggested, and
changes the priority ordering of what to fix first (see the revised A.3 and Part C.1).

**Why, mechanistically** (directly explains the observed shape, not just describes it):
in the pre-transition regime, `f(ε_struct)` sits essentially flat at `f0` — visible
directly in the same figure's middle panel, where both the internal and observable `k`
traces are flat from throughput 0 all the way to ~3,500 A.h, only then rising sigmoidally
toward saturation by ~6,500 A.h. Since `dv_thickness ≈ f · dv_solid` (§ pore-buffering
math, `reaction_driven_porosity.py`), the measured/buffered amplitude during this whole
flat-`f` stretch is just `f0 × (unbuffered particle amplitude)(t)` — a **constant
multiple** of the underlying unbuffered swelling amplitude. And the unbuffered amplitude
itself declines over life in this recipe (already established via the sensitivity work
behind `CHANGES.md`'s "peak-height" finding — driven by LAM reducing active material, and
by the fact these are CC-CV cycles bounded by fixed voltage cutoffs, so the *total* charge
moved each cycle shrinks as capacity fades, which in this recipe evidently outweighs the
opposing effect of the shrinking active-material pool concentrating more current onto the
particles that remain). A roughly-constant multiplier can't change that declining
*shape* — it can only rescale it. So for **most** of pre-hump life, the buffered signal is
mechanically forced to just be a scaled-down copy of the declining unbuffered trend, and
only once `ε_struct` finally nears the transition band does `f`'s rise become large
enough, and fast enough, to overcome the ongoing decline and turn the curve upward into
the hump.

**Update — recalibration was tried and does NOT work; ruled out, not just risky.** A
direct 2×2 sweep (band `WIDTH` ∈ {0.005, 0.010, 0.020, 0.040} at fixed `eps_min=0.08`, and
`f0` ∈ {0.5, 0.6, 0.7, 0.85} at a fixed band) shows every single curve — all four widths,
all four `f0` values, including the widest/most permissive settings tested — still
**declines from BoL before rising**. Widening the band or raising `f0` only changes *where*
the eventual transition sits and how sharp/tall it is; it never flips the sign of the
pre-hump trend. This is not a coincidence of the specific values swept — it's structurally
guaranteed by the mechanism itself: `buffered_amplitude(t) = f(ε_struct(t)) ×
unbuffered_amplitude(t)`, and **before** `ε_struct(t)` actually reaches the transition
band, `f` is sitting at (very close to) the constant `f0`, *regardless of how wide the
band is or what `f0`'s value is* — widening the band or changing `f0` only changes what
happens once the transition starts, and can't touch the interval before it, because by
construction `f` isn't moving yet during that interval. So the pre-band trend is entirely
inherited from `unbuffered_amplitude(t)`'s own sign, which is declining, for *any* choice
of these two parameters. The only lever that could in principle change this is
repositioning `eps_min_transfer` itself close enough to the BoL `ε_struct` value that the
cell starts *inside* the transition band from cycle 1 — but that's not "recalibration"
within the mechanism's intended meaning any more, it amounts to asserting the electrode is
already near percolation closure at BoL, contradicting what `eps_min_transfer` is supposed
to represent (a genuine degradation-endpoint threshold, not a BoL state), and is exactly
the kind of "implausible parameter stretch" flagged below as a warning sign.

**Conclusion: within the single-mechanism pore-buffering model as currently structured,
there is no physically reasonable parameter choice that produces an increasing pre-hump
trend.** This promotes what was originally framed as a parallel, "worth checking" caveat
into the **leading hypothesis**: the majority (increasing) pattern most likely needs a
**different degradation mechanism to become the dominant one**, not (only) retuned
buffering-transition parameters.

**Leading candidate, and the cheapest of all to test — LLI-vs-LAM dominant-pathway
balance, using mechanisms already fully implemented.** This doesn't require any new
submodel at all, just rebalancing rate constants already present in the existing recipe.
It's the standard eSOH distinction between the two ways lithium/capacity is lost, and the
two have *opposite* effects on the remaining active material's per-cycle stoichiometry
swing (`Δx`), which is what actually sets particle-level swelling amplitude:

- **LLI-dominant** (lithium consumed irreversibly by SEI growth): the negative and
  positive electrodes' stoichiometry windows increasingly misalign relative to the fixed
  voltage cutoffs as cyclable lithium shrinks — the *utilized* `Δx` per cycle for the
  remaining active material **shrinks** over life. Declining particle-level swelling
  amplitude follows directly.
- **LAM-dominant** (active material lost to cracking, with lithium inventory relatively
  preserved): the *same* cyclable lithium now has to be accommodated by a *smaller*
  remaining active-material pool each cycle — the utilized `Δx` for the surviving
  particles **grows** over life, a well-established eSOH result (LAM concentrates the same
  lithium throughput onto fewer particles). Increasing particle-level swelling amplitude
  follows directly, right up until the pore-buffering transition takes over near the knee.

This maps cleanly onto the current recipe and gives a natural explanation for *why* our
simulated baseline currently sits in the declining regime: `si_gr_expansion`'s tuned
recipe uses `SI_MULT = 18.0` — a large, deliberate boost specifically to silicon's SEI
exchange-current density (tuned earlier this session for knee-timing/LLI-matching
purposes, not for pre-hump shape) — a strongly LLI-weighted choice. It's very plausible
this pushes the recipe's LLI/LAM balance far enough toward LLI-dominance that the
window-shrinking effect currently dominates the window-growing effect from Si/Gr cracking
-driven LAM, net declining. **Concrete, cheap test**: sweep the existing recipe's relative
weighting — dial `SI_MULT`/`GR_DIV` down (less LLI-weighted SEI growth) and/or the LAM
proportional-rate constants (`SI_LAM_PROP`/`GR_LAM_PROP`)/crack-rate multipliers up (more
LAM-weighted), and check whether the *unbuffered* particle amplitude's own trend (already
computed as `dv_solid`'s within-cycle swing in the existing scripts, independent of
buffering) flips sign from declining to increasing as LAM's contribution starts to
dominate LLI's. If it does, that's strong, direct, mechanistic support — no new model code
required, only a different point in the existing recipe's own parameter space — and it
also gives a coherent explanation for *why different operating conditions* land in
different regimes: temperature accelerates SEI reaction kinetics strongly (Arrhenius),
plausibly pushing high-T toward LLI-dominance (matching its observed declining trend and
the SI_MULT-heavy recipe's default behaviour); lower external pressure and higher C-rate
both increase mechanical stress on the particles, plausibly pushing those conditions
toward LAM/cracking-dominance instead (matching their observed increasing trends) — turning
Part C's rate/temperature/pressure work into a direct test of this same hypothesis, not a
separate track.

**If this alone doesn't fully explain the pattern** (worth checking honestly, since a
recipe rebalanced for pre-hump shape also needs to still match the recipe's other
established validation targets — knee timing, LLI/LAM magnitudes at 50% SoH, etc., not
just this one curve), the following genuinely-new-mechanism candidates remain the next
things to try, in roughly increasing order of implementation cost:

- **Lithium (soft) plating.** Not enabled anywhere in the current `si_gr_expansion`
  recipe. Plated lithium has its own volume and its own (partially reversible, stripping-
  on-discharge) swelling contribution, distinct from particle-level Si/Gr lithiation
  breathing — if plating onset grows gradually over life (very plausible in a high-silicon
  anode as porosity/kinetics shift), it would add a genuinely separate, *growing* reversible
  expansion term on top of whatever the buffering mechanism is doing, with no need for the
  buffering transition band to be doing all the work itself.
- **Gas generation.** The SEI reaction this whole session's dry-out work is built on (Li
  et al. 2022, Eq. 1: 2Li⁺ + 2EC + 2e⁻ → SEI + **C₂H₄ gas**) explicitly produces gas as a
  side product — and that paper's own derivation explicitly *ignores* the gas phase
  ("if we ignore the gas phase and the solute-volume effects...", ec_dryout's source
  material) to keep the electrolyte-volume bookkeeping tractable. In a real pouch cell
  without venting, accumulating gas is a directly measurable, physically real, and
  currently **completely unmodeled** contributor to cell thickness that has nothing to do
  with particle swelling, pore buffering, or electrolyte volume at all.
- **Evolving mechanical stiffness of the SEI/pore network itself** (not just its
  *porosity fraction*). The deferred pressure-coupling notes (`pore_buffering_
  implementation_plan.md`) reference a third quantity alongside `K_stack`/`C_pore(ε)`:
  `E_s(φ)`, the solid/SEI network's own compliance as a function of degradation state.
  Fresh SEI is soft/gel-like; aged SEI is typically stiffer/more mineralized (a
  well-documented real phenomenon) — if the pore network's *stiffness*, not just its
  *volume*, changes over life, that's a distinct physical driver of a rising transmitted
  fraction that isn't really "the same mechanism, different parameters" so much as a
  related-but-separate piece of physics the current single-porosity-driven `f(ε_struct)`
  doesn't represent at all.
- **Growing cathode-side contribution.** The composite model's cathode volume-change
  parameter was deliberately tuned down to a small, roughly fixed ~10% max contribution
  (`CHANGES.md`) based on reference data. If in reality the cathode's *own* separate
  degradation (particle cracking, surface reconstruction) causes its contribution to grow
  over life, that's a source of increasing measured cell-level expansion with nothing to
  do with the anode-side pore-buffering story this whole plan has otherwise focused on.
- **Stack-level mechanical relaxation over cycling** — distinct from the *externally
  applied* pressure setpoint (Part C.4's `K_stack(P_stack)`): the stack's own effective
  stiffness could soften over life independent of the fixed external pressure, e.g. via
  binder or current-collector degradation, separator creep, etc. — a genuinely different,
  currently entirely unmodeled degradation pathway.

None of these can be ruled in or out from the existing simulation output alone — they're
listed here as fallback candidates if the LLI-vs-LAM balance test above doesn't fully
explain the pattern, in roughly increasing order of implementation cost. **Lithium
plating** is the cheapest of the four: PyBaMM already has a working plating submodel, so
this needs a model-option change and a parameter/kinetics estimate, not new code — a much
smaller lift than gas generation, `E_s(φ)` stiffness evolution, cathode-side, or
stack-relaxation, all of which would need genuinely new model structure. Gas generation is
the next most natural candidate given this session's dry-out work already derived the
reaction stoichiometry that produces it (Eq. 1) — "just" not its volume/venting
consequences.

**Priority order overall**: (1) LLI-vs-LAM balance sweep on the existing recipe — cheapest,
uses only already-implemented mechanisms; (2) lithium plating — cheap, existing submodel,
new physics; (3) gas generation — moderate, stoichiometry already derived this session; (4)
`E_s(φ)`/cathode-side/stack-relaxation — genuinely new model structure, try only if 1–3
don't close the gap.

Pressure-coupling (Part C.4) remains valuable regardless of which alternative mechanism
turns out to be right — it's still the best candidate for *why low pressure specifically*
shows a steeper/earlier-starting increase than baseline or high-rate, whether that acts
through `f(ε_struct)` directly or through one of the alternative mechanisms above (e.g.
lower external pressure plausibly also lowers `K_stack` in a genuinely mechanical sense, or
changes SEI/plating morphology and hence gas retention). Part C.1 below reflects the
ruled-out status of recalibration and prioritises the LLI-vs-LAM balance test (A.0a below)
ahead of adding any genuinely new physics.

### A.0a Test result: LLI-vs-LAM balance hypothesis — partially confirmed, one lever cleanly supports it, the other is confounded

Tested directly (`expansion_pattern_test/lli_vs_lam_balance_test.py`): the existing
`si_gr_expansion` recipe, unbuffered (`"pore buffering": "false"`, so `"Cell thickness
change [m]"` is the raw signal), 120 cycles, two independent 1-parameter sweeps off the
tuned baseline (`SI_MULT=18.0`, `LAM_MULT=1.0`).

**Sweep 1 — `SI_MULT` down (less LLI-weighted), `LAM_MULT` fixed at 1.0 — clean, monotonic
support.** `SI_MULT` ∈ {18.0, 9.0, 4.5, 2.25} gives early-life amplitude slopes of
**−0.0053, −0.0037, −0.0025, −0.0017 µm/A·h** — monotonically approaching zero as `SI_MULT`
drops, exactly the direction the hypothesis predicts, and cleanly isolated: final LAM
stayed essentially flat (~17.0–17.0%) across all four points while final LLI dropped
20.68% → 10.18% → 7.82% → 6.47%, confirming only the intended lever moved. **The slope
never actually crossed zero within the tested range** — even at an 8× reduction it's still
declining, just three times less steeply. This is real, unambiguous evidence that reducing
LLI-weighting moves the trend the right way, but on its own (without also raising LAM) it
doesn't fully flip the sign at these `SI_MULT` values.

**Sweep 2 — `LAM_MULT` up (more LAM-weighted), `SI_MULT` fixed at baseline — moved the
WRONG way, and turned out to be confounded, not a clean refutation.** `LAM_MULT` ∈ {1.0,
2.0, 4.0} gives slopes of **−0.0053, −0.0066, −0.0111 µm/A·h** — *more* declining, opposite
to the predicted direction. But this sweep did not actually isolate "more LAM" the way it
was intended to: final LLI *also* rose substantially alongside LAM (20.68% → 23.73% →
43.72%, tracking LAM's 17.03% → 38.43% → 80.04%), even though the crack-rate parameters
themselves were held fixed. This happens because LAM and cracking-driven SEI-on-cracks LLI
share the same underlying stress-driven mechanics in this submodel — scaling up the LAM
proportional-rate constant doesn't leave crack-driven LLI untouched, it changes the whole
coupled stress trajectory. At `LAM_MULT=4.0` the cell effectively fails within the run (see
the plot: amplitude collapses to ~0 around 550–600 A·h, `LAM_neg=80%` — the electrode has
lost the large majority of its active material, well past any regime "moderate LAM
increase" was meant to represent). So sweep 2's declining result is consistent with "more
LLI still dominates when LLI rises alongside LAM," not evidence against the hypothesis
itself — it's evidence that this particular lever (scaling `SI_LAM_PROP`/`GR_LAM_PROP`
while holding nominal crack rate fixed) doesn't cleanly separate the two pathways in this
model.

**Net assessment**: genuinely supportive, not yet conclusive. The clean half of the
experiment (sweep 1) behaves exactly as predicted and rules out "the hypothesis is simply
wrong." The confounded half (sweep 2) doesn't refute it, but doesn't confirm the LAM side
either — it shows this specific implementation lever can't isolate LAM from LLI at high
`LAM_MULT`. See `expansion_pattern_test/lli_vs_lam_balance_test_result.png` for the plot.

### A.0b Test result: combined sweep — a U-shape, not monotonic progress; the LAM side is weaker than expected, possibly counterproductive

Tested the natural follow-up (`expansion_pattern_test/combined_sweep_test.py`): reduce
`SI_MULT` *and* increase `LAM_MULT` together, diagonally through five paired points from
the tuned baseline, at more moderate values than sweep 2's extremes.

| (SI_MULT, LAM_MULT) | Early-life slope [µm/A·h] | Final LLI | Final LAM (neg) |
|---|---|---|---|
| (18, 1.0) — baseline | −0.0053 | 20.7% | 17.0% |
| (9, 1.5) | −0.0041 | 12.7% | 27.0% |
| (4.5, 2.0) | **−0.0027** | 13.2% | 37.9% |
| (2.25, 2.5) | **−0.0026** | 14.7% | 52.2% |
| (1.125, 3.0) | −0.0051 | 18.7% | 69.4% |

**Result: a clear U-shape, not monotonic progress toward zero.** The slope improves
(halves, roughly) from the baseline through the middle two points, then **reverses** and
lands almost back at the baseline value by the most aggressive point — visible directly in
both the amplitude and SoH panels of the plot (the last point's SoH collapses to ~32% by
throughput 800 A·h, far faster than any other point; `LAM_Si=91.6%` there — essentially
complete silicon depletion, a cell-failure regime, not "aggressive but reasonable"
degradation). The co-reduction of `SI_MULT` DID successfully damp sweep 2's LLI-runaway
confound through the middle of the sweep (LLI stays ~13–15%, not climbing to 40%+ the way
unmitigated `LAM_MULT` alone did) — but at the extreme point, LAM itself has become so
severe that LLI creeps back up anyway (18.7%, nearly baseline), suggesting there's a
second, independent route from extreme cracking/LAM back to elevated LLI even with the
base SEI rate turned down.

**More importantly — and worth stating plainly rather than glossing over**: the best
result in this entire investigation (all three test scripts combined) is still **sweep
1's pure `SI_MULT=2.25` point, slope −0.0017** (A.0a) — *better* (closer to zero) than
this combined sweep's best point (−0.0026 at `SI_MULT=2.25, LAM_MULT=2.5`). In other
words, holding `SI_MULT` fixed at 2.25 and *adding* `LAM_MULT=2.5` on top made the slope
*worse*, not better. This is a real, unforced result, not an artefact of the LLI confound
(LLI was successfully controlled in this range) — it suggests the simple eSOH argument
("LAM concentrates the same lithium onto fewer particles, growing their swing") is
missing something. The most likely gap: that argument implicitly assumes the *active*
material fraction is what matters and lithium inventory is what's conserved, but it
ignores that **material lost to LAM stops swelling entirely** (it's no longer
participating in lithiation/delithiation at all) — so the aggregate cell-level amplitude
depends on two competing LAM effects, not one: (a) surviving particles' own stoichiometry
swing grows (the effect the hypothesis emphasised), but (b) the active *volume fraction*
still contributing to the swelling signal at all shrinks. Empirically, at the LAM
magnitudes tested here, effect (b) appears to dominate effect (a).

**Revised conclusion**: LLI reduction (A.0a's sweep 1) remains the most reliable, cleanly-
demonstrated lever toward the majority (increasing) pattern found so far — genuinely
supportive of that half of the hypothesis. The LAM side is now on weaker footing than
originally proposed: it doesn't reliably help even when its LLI-confound is controlled,
and may need to be reasoned about differently (e.g. a pure, uncracked/undamaged LAM
pathway isolated from cracking-driven surface generation entirely — not available as a
simple parameter lever in the current recipe, since LAM here is explicitly stress/crack-
driven by construction) rather than assumed. Given this, the fallback mechanisms in A.0
(lithium plating first, per the original priority order) are now relatively more
attractive next steps than pushing the LAM lever further.

### A.0c Test result: amplitude measurement convention — makes no difference, ruled out as an explanation

A separate methodological check (`expansion_pattern_test/amplitude_convention_test.py`):
everything so far measured within-cycle amplitude as "this cycle's charge-end minus this
cycle's discharge-end" (charge-referenced). Tested an alternative convention —
"this cycle's discharge-end minus the *previous* cycle's charge-end" (discharge-
referenced, using where charging initially left the cell as the baseline) — on the
tuned baseline recipe, to check whether the declining pre-hump trend found throughout A.0
is an artefact of which reference points were chosen, rather than a real property of the
underlying signal.

**Result: no meaningful difference.** Early-life slopes: −0.00530 µm/A·h (charge-
referenced) vs. −0.00532 µm/A·h (discharge-referenced) — the two trajectories are visually
indistinguishable throughout the full run (see the plot). This makes sense in hindsight
given the cycling protocol has no rest step between legs, so the cell's thickness
trajectory is essentially continuous across cycle boundaries and both conventions sample
nearly the same underlying smooth signal, just offset by one cycle. **Conclusion: the
declining pre-hump trend is a robust property of the underlying degradation physics, not
a measurement-convention artefact** — this rules out one more alternative explanation and
leaves the LLI-vs-LAM balance question (A.0a/A.0b) as the genuine open item.

### A.0d Test result: buffered confirmation — buffering rescales but doesn't flip the sign (within this window); one promising sign for longer runs

Re-ran A.0a's cleanest sweep (`SI_MULT` ∈ {18.0, 9.0, 4.5, 2.25}, `LAM_MULT=1.0`) with
`"pore buffering": "true"` enabled this time (`expansion_pattern_test/
buffered_confirmation_test.py`), to check whether the actual observable/buffered signal —
not the idealised unbuffered one used throughout A.0a–c — tells a different story.

| SI_MULT | Unbuffered slope | Buffered slope | Ratio |
|---|---|---|---|
| 18.0 | −0.00530 | −0.00364 | 0.687 |
| 9.0 | −0.00369 | −0.00249 | 0.675 |
| 4.5 | −0.00247 | −0.00164 | 0.664 |
| 2.25 | −0.00174 | −0.00114 | 0.655 |

**Result: buffering consistently helps (makes the slope less negative) but doesn't flip
the sign for any of the four within this 120-cycle window.** The ratio is close to
constant (~0.66–0.69) across all four points — this is not a coincidence: the internal
transfer ratio `k` (right panel of the plot) stays pinned flat at `f0=0.70` for the entire
120 cycles at `SI_MULT` ∈ {9, 4.5, 2.25} (their own buffering transition hasn't been
reached yet at this throughput — consistent with lower `SI_MULT` meaning slower SEI growth
meaning slower `ε_struct` decline meaning a later transition), so for those three,
`buffered_amplitude ≈ f0 × unbuffered_amplitude` is a near-exact constant rescale, and the
early-life slope scales by essentially the same `f0≈0.7` factor as a direct consequence —
not a new, independent finding beyond what the pore-buffering math (A.0) already predicts.

**`SI_MULT=18` (the tuned baseline) is the one case that actually reaches its own
transition within the window**: `k` rises from 0.70 to 0.9995 by throughput ~900 A·h, and
the buffered amplitude trace shows the full decline-then-hump shape already familiar from
`pore_buffering_degradation_test.py`'s validated results (minimum ~9.8 µm around 600–650
A·h, peak ~14.0 µm around 900 A·h). The "early-life" slope reported here is fit over only
the first half of the run, which is still mostly pre-transition for this case too, so it
still shows the same ~0.69 scaling rather than reflecting the later upturn — a reminder
that "early-life slope" is a deliberately narrow window and doesn't capture the eventual
hump for any of these points.

**One genuinely promising observation worth flagging for a longer run**: visually, the
`SI_MULT=2.25` buffered curve (lightest green, top-left panel) is not just linearly
declining — it visibly *flattens* over the back third of the 120-cycle window, more than
a simple constant-rescale of the unbuffered trend would produce on its own. That's
consistent with `ε_struct` for this much-slower-degrading recipe just beginning to
approach its own transition band by the end of the window. **Not yet confirmed**, since
120 cycles isn't enough to see whether it continues flattening into a genuine upturn —
that would need one of the reduced-`SI_MULT` recipes run for substantially longer (the
1000-cycle-scale protocol, or at least several hundred cycles) to actually observe.

**Revised bottom line**: buffering, on its own, does not rescue any of these into the
majority (increasing) pattern within a short window — it's mathematically just a rescale
wherever the transition hasn't engaged yet, exactly as A.0's own derivation predicted, so
this test is best read as a *confirmation* of that math rather than a new mechanism. The
one substantive new lead is the `SI_MULT=2.25` flattening — worth a longer follow-up run
if pursuing this further, but not yet evidence of a sign flip.

### A.0e Test result: extended SI_MULT=2.25 to 500 cycles — the flattening lead did NOT pan out; a bigger, unintended problem instead

Followed up on A.0d's lead directly (`expansion_pattern_test/si_mult_2p25_extended_test.py`):
ran `SI_MULT=2.25`, `LAM_MULT=1.0`, buffered, extended to 500 cycles (up from 120).

**The run terminated early, at 325 cycles / 2238 A·h throughput, hitting the 50% SoH
floor** (`LLI=20.2%`, `LAM_neg=57.7%` at that point). Over this full run to end-of-life:

- **`k` never meaningfully leaves the `f0=0.70` plateau**: it stays flat to 4 decimal
  places until throughput ~1700 A·h, then rises only to **0.7005** by the 50%-SoH cutoff —
  a 0.07% relative change, compared to the tuned baseline's `k` rising the *entire* way
  from 0.70 to 1.00 within its own (much shorter) life. The buffering transition
  essentially never engages for this recipe within its useful life.
- **The amplitude does NOT turn around**: it declines mildly through the region that
  looked like flattening in the 120-cycle window (roughly 1000–1300 A·h, visible as a
  slight inflection in the plot), but then **resumes declining and actually accelerates**
  as the cell approaches end-of-life, dropping from ~12 µm at BoL to ~3.7 µm by the 50%
  SoH cutoff — a much steeper decline in the back half than the front half, the opposite
  of what the earlier "flattening" observation suggested might be coming.

**What this actually shows**: the local flattening seen in the 120-cycle diagnostic was a
transient feature, not the start of a genuine turnaround — over the full lifetime it's
swamped by accelerating LLI/LAM-driven decline. More importantly, this surfaces a
previously-unrecognised **structural tension** in the "reduce `SI_MULT`" fix itself:
lowering `SI_MULT` doesn't just slow the LLI pathway (the intended effect, A.0a) — it also
slows `ε_struct`'s own decline enough that **the buffering transition (and therefore any
hump at all) may never be reached before the OTHER degradation modes (LAM, LLI-driven
capacity fade, independent of buffering) end the cell's useful life first.** A fix that
makes the pre-hump slope less negative is not useful if it also makes the hump itself
unreachable within a realistic lifetime — worth keeping in mind for any further `SI_MULT`-
based tuning: the pre-hump shape and the transition's reachability are coupled by the same
parameter, not independently adjustable.

### A.0f Test result: SEI film thickness — confirmed missing from cell thickness, but too small/too slow to explain the amplitude trend; real gap for a different metric

Checked directly (`expansion_pattern_test/sei_thickness_diagnostic.py`), prompted by the
question of whether the SEI layer's own physical thickness could be a place the model
differs from reality. First, read the source: **confirmed SEI thickness is currently NOT
included in `"Cell thickness change [m]"` anywhere** — `base_mechanics.py`'s
`_extract_mechanical_variables` builds `electrode_thickness_change` purely from particle
volume change (`v_change`, driven by the lithiation stoichiometry `t_change` function);
SEI thickness (`L_sei_k`) only ever feeds into the porosity decrease
(`delta_eps_k = -a_k * L_tot` in `reaction_driven_porosity.py`), never into the mechanical/
expansion output. This surfaced a genuine structural asymmetry: particle swelling gets a
buffered-vs-transmitted *partition* as pores close (the pore-buffering submodel), but
SEI's own volume growth does not — it's unconditionally treated as 100% pore-consuming
forever, even once pores are already closed, at which point the model just relies on the
numerical porosity floor rather than asking where that continued volume goes.

Two distinct hypotheses tested for magnitude, on the base case (SI_MULT=18, buffered):

- **Hypothesis A (naive additive — add raw SEI thickness to cell thickness)**: raw SEI
  film thickness is **5–80 nm** over the tested window — three to four orders of magnitude
  smaller than the ~12 µm amplitude scale. Even the *entire* raw SEI thickness, added
  wholesale, would barely register. **Ruled out** — SEI film thickness itself is simply
  too thin to matter at this scale.
- **Hypothesis B (structural — what if SEI volume were *fully* redirected to thickness,
  the same way saturated particle swelling is)**: converted the SEI-driven porosity term
  to an equivalent cell-level length using the *same* specific-surface-area scaling
  (`a_k`) the porosity submodel already computes internally, times the same
  `n_electrodes_parallel * L_electrode` conversion particle volume change uses. This
  upper-bound equivalent thickness is **genuinely substantial and growing**: 0.35 → 1.45
  µm over just the first 180 A·h of throughput (~20 cycles) — order-of-magnitude
  comparable to the reversible amplitude itself, and still climbing. **This is a real,
  currently-missing physical contribution, not a negligible one.**

**But — critically — its *within-cycle* oscillation is tiny**: even this fully-redirected
upper bound only swings ~50–60 nm between a cycle's charge-end and discharge-end (<0.5% of
the ~12 µm amplitude scale), because SEI growth is continuous and slow relative to a single
cycle's timescale, not something that visibly grows-then-shrinks within one charge/
discharge. So Hypothesis B's large, real magnitude shows up almost entirely as a
**growing irreversible/baseline thickness contribution**, not as a change to the
**reversible amplitude** (charge-end minus discharge-end within the same cycle) that this
whole A.0 investigation has been about.

**Conclusion**: SEI film thickness is *not* a plausible explanation for the declining-vs-
increasing pre-hump *amplitude* pattern — ruled out for that specific question. But it
*is* a genuine, quantitatively significant, currently-unmodelled contributor to the cell's
**absolute/irreversible expansion over life** (a different, real metric this project
hasn't been tracking at all) — worth fixing on its own merits as a separate follow-up
(implementing the same buffered-vs-transmitted partition for SEI volume that particle
swelling already has, redirecting it to thickness once pores saturate), just not as an
explanation for this specific pre-hump puzzle.

### A.0g PREMISE CORRECTION: updated experimental data flips the majority/minority framing — low pressure is the outlier, not the norm

An updated, more complete experimental figure supersedes the 4-condition reading in A.1
below — source file: `fig2_population_lead_time_v2.png` (same folder as this plan).

![Population lead-time figure: panel a shows discharge capacity and reversible expansion vs. EFC for 5 conditions (baseline, low pressure, high C-rate, high temperature, narrow SoC window) with lead/knee EFC annotated; panel b overlays all conditions aligned to EFC-EFC_knee; panel c summarises lead EFC by category; panel d shows a combined high-rate + narrow-SoC-window condition.](fig2_population_lead_time_v2.png)

**Panel a** adds a **fifth condition** (a narrower voltage/SoC window, 3.15–4.12 V vs. the
2.6–4.2 V baseline) to the original four, with precise lead/knee EFC now readable directly
off the figure:

| Condition | Lead EFC | Knee EFC | Pre-hump shape |
|---|---|---|---|
| Baseline (25 °C, 0.33C, 103 kPa, 2.6–4.2 V) | 121 | 145 | flat-to-mildly-declining (~52 µm), then a hump to ~64 µm right at the knee |
| **Low pressure (25 °C, 0.33C, 34 kPa, 2.6–4.2 V)** | 119 | 158 | **rises continuously from ~53 µm to a ~78 µm peak at the knee, no flat plateau** |
| High C-rate (25 °C, 2C, 103 kPa, 2.6–4.2 V) | 95 | 143 | flat-to-mildly-declining (~44 µm), then a hump to ~56 µm — same shape as baseline, compressed into the shortest window of the five |
| High temperature (45 °C, 0.33C, 103 kPa, 2.6–4.2 V) | 296 | 319 | flat/mildly declining (~50→47 µm), muted hump to only ~49 µm — same shape as baseline, stretched over by far the longest window |
| Narrow SoC window (25 °C, 0.33C, 103 kPa, 3.15–4.12 V) | 255 | 298 | flat (~37 µm), then a sharp hump to ~52 µm right at the knee — same shape as baseline |

**Panels b and c** add cross-condition structure the earlier 4-condition reading didn't
have: panel b overlays all conditions aligned to `EFC − EFC_knee`, and directly visualises
what the table says — the low-pressure squares are the clear outlier, rising earliest and
highest of any condition (peak ~88 µm, the tallest hump in the whole dataset) with no flat
lead-in, while every other condition (including the 172 kPa *high*-pressure stars, both
temperatures, both C-rates, and the SoC-window variants) stays comparatively flat through
most of the lead-EFC window before a hump concentrated right at the knee. Panel c's lead-EFC
breakdown by category (overall median 36 EFC) shows temperature has the *shortest* lead
(~20 EFC, both 25 °C and 45 °C), C-rate the *longest* (~48–54 EFC), with pressure and SoC
window in between (~24–44 EFC) — i.e. even the *timing* of when the precursor signal first
appears varies by category, on top of the *shape* difference this section focuses on.

**Panel d is a new, so-far-unexplained data point worth flagging on its own**: a
*combined* high-rate + narrow-SoC-window condition (25 °C, 2C, 103 kPa, 50–100 % SoC)
shows discharge capacity declining smoothly with **no visible knee** over the full ~800 EFC
tested, and reversible expansion declining **monotonically** from ~58 µm to ~28 µm with
**no hump at all** in that range — unlike either single-variable condition (high C-rate
alone, or any of the SoC-window variants in panel b) which each showed a hump within their
own tested windows. Not yet investigated why combining the two variables suppresses the
hump entirely (or pushes it beyond 800 EFC) — flagged here as an open item, not folded into
the reasoning below, since it doesn't obviously fit the "pressure is the only qualitative
outlier" framing that panels a–c support.

**Four of five panel-a conditions share the same qualitative shape as baseline**
(flat/declining pre-hump, hump at the knee) — only their *timing* differs (knee EFC ranges
143→319 across them). **Low pressure alone is qualitatively different**: no flat region at
all, rising continuously from cycle 1. This reverses A.0's original characterisation (which
read baseline/low-pressure/high-C-rate as the increasing majority and high-temperature as
the declining minority, based on an earlier, less complete version of the figure) — with
the fuller picture, **declining-then-hump is the majority pattern, and it's the pattern the
model's existing default calibration already reproduces** (A.0's own decline-then-hump
finding for the simulated baseline). The thing that actually needs explaining is narrower
and more specific than A.0–A.0f assumed: not "why does the model decline when most
conditions increase," but **"why does low pressure specifically increase when everything
else — including large swings in temperature, rate, and SoC window — declines like
baseline."**

**Theory: pressure is a mechanical/structural parameter; temperature, rate, and SoC window
are kinetic/rate parameters — and only a structural change can alter the pre-hump *shape*,
not just its *timing*.**

Temperature, C-rate, and SoC window all act on the same thing: how fast the underlying
chemistry consumes pore volume — i.e. how fast `ε_struct` declines toward the buffering
transition band (`reaction_driven_porosity.py`). None of them touch the *relationship*
between `ε_struct` and the transmitted fraction `f(ε_struct)` itself. Per A.0's own
derivation, the pre-transition regime is `f ≈ f0` (a plateau, essentially constant)
regardless of how quickly the cell approaches the band — so any condition that only
changes the *rate* of approach should preserve the same flat/declining-then-hump shape,
just shifted earlier or later in EFC. That is exactly what baseline, high-C-rate,
high-temperature, and the SoC-window condition show: the same shape at different knee
timings (143, 133/short-window, 319, and 298 EFC respectively vs. baseline's 145).

Pressure is categorically different: it's not a rate constant on any reaction or transport
process, it's a **mechanical boundary condition** — how much the electrode/separator stack
is externally confined. It plausibly changes `f`'s functional relationship to `ε_struct`
itself (i.e. `K_stack` in Part C.4's notation), not just how fast `ε_struct` gets there,
because it acts on the model's mechanics at *every* porosity state, including from BoL —
not on a reaction rate that only matters while it's running. Two concrete, physically
grounded mechanisms for why lower confinement specifically produces a *continuous* rise
rather than a *late, sharp* one:

- **Particle rearrangement is unconstrained at low pressure.** Under normal/higher stack
  pressure, particles are already pressed together at BoL; as SEI growth consumes pore
  volume, particles can't easily repack into the shrinking pore space, so nothing
  macroscopically moves until the pore network is nearly closed — a late, sharp `f`
  transition, exactly the shape the "physical" transmitted-fraction form (and the model's
  current default calibration) already produces. Under low pressure, particles have room
  to continuously resettle as pores shrink, so the "give" is distributed across life
  instead of reserved for one late event near the knee.
- **Less opposing force from the stack/fixture at every porosity, not just near
  saturation.** A softer, less-compressed stack transmits a larger fraction of any given
  internal volume change to the measurable cell boundary throughout life, not only once
  pores are nearly exhausted — i.e. a smaller effective `K_stack` biases `f` upward across
  the *whole* trajectory, rather than only flipping it up sharply once `ε_struct` nears
  `eps_min_transfer`.

**This directly revives, and correctly re-scopes, A.0's original "recalibration" idea.**
A.0's width-×-`f0` sweep concluded recalibrating the transition band cannot produce an
increasing pre-hump trend for the *baseline* recipe — true, and still true, since baseline
doesn't need one. But that sweep only tested a single, pressure-independent band. The
band-widening/early-starting direction that sweep explored is exactly the right shape of
fix for the *pressure-varied* condition specifically — it just needs to be applied as a
`K_stack(P_stack)`-style, pressure-dependent modulation of the transition band/compliance
(Part C.4), not a universal recalibration of the one fixed baseline recipe.

**What this means for A.0–A.0f's LLI-vs-LAM investigation.** That work is real and stays
valid as a record (A.0a's clean `SI_MULT`-down result in particular is a genuine, useful
finding about this recipe's degradation-pathway balance), but its motivating premise — that
the model's declining baseline needs fixing because most experimental conditions increase
— no longer holds. It should be read as a parallel, lower-priority investigation now, not
the leading candidate for explaining the pre-hump pattern. The **leading, best-supported
candidate going forward is Part C.4's pressure-coupling extension**, since it's the one
mechanism that (a) is structurally distinct in exactly the right way (a mechanical boundary
condition vs. a reaction/transport rate), (b) already has a natural implementation path
(`K_stack(P_stack)` replacing the existing fixed `k_cmax` in
`_transmitted_fraction_physical` — no new submodel, see C.4), and (c) makes a falsifiable,
narrow prediction: implementing it and running baseline vs. 34 kPa with everything else
held fixed should reproduce the continuously-increasing low-pressure shape *without*
needing to touch the recipe's LLI/LAM balance at all, since baseline's own shape is already
correct as-is.

**Update — prediction (c) was tested directly and falsified; see A.0h immediately below.**
A pure `K_stack` rescale (holding `eps_transfer_width`/`eps_min_transfer`/the underlying
degradation kinetics fixed) does not produce a continuously-increasing trend, and moves in
the *wrong* direction from baseline (lower `K_stack` makes the slope *more* declining, not
less). The theory in this section — pressure as a mechanical, not kinetic, parameter — is
still the best explanation for *why* low pressure is qualitatively different from
temperature/rate/SoC-window. What's wrong is the specific mechanism proposed: a constant
partition-height rescale can't do it. A.0h works out why, and what a viable pressure-
coupling mechanism needs instead.

### A.0h Test result: `K_stack` sweep — falsified, and in the wrong direction; a constant partition rescale cannot produce a rising trend

Tested directly (`expansion_pattern_test/k_stack_pressure_sweep_test.py`), enabled by the
new `"pore buffering stack compliance": "independent"` model option (decouples `K_stack`
from `f0` — see `reaction_driven_porosity.py`/`lithium_ion_parameters.py`, and the "easily
switch back" note at the end of this section). Base-case recipe (`SI_MULT=18`,
`LAM_MULT=1.0`, identical to `test_pore_buffering`'s validated baseline), `eps_min_transfer`/
`eps_transfer_width`/`f0` all held at their existing values, `K_stack` swept from the
baseline value (`(1-0.7)/0.7 = 3/7 ≈ 0.4286`, representing 103 kPa) down to a near-limiting
`1e-4`, 120-cycle diagnostic.

**Result: every value declines, none cross zero, and the direction is backwards.**

| `K_stack` | Implied BoL plateau `f0_eff=1/(1+K)` | Early-life slope [µm/A·h] |
|---|---|---|
| 0.4286 (baseline, 103 kPa) | 0.700 | **−0.00364** |
| 0.2 | 0.833 | −0.00441 |
| 0.1 | 0.909 | −0.00483 |
| 0.05 | 0.952 | −0.00506 |
| 0.01 | 0.990 | −0.00525 |
| 0.0001 (near-limit) | 0.9999 | −0.00530 |
| Unbuffered reference (no partition at all) | — | −0.00530 |

As `K_stack → 0`, the buffered trajectory converges exactly onto the unbuffered reference
(as it must — `f0_eff → 1` means essentially everything transmits, so buffering becomes a
no-op) — and the unbuffered reference is *more* declining than baseline, not less. So
**lowering `K_stack` moves the slope away from zero, not toward it**: baseline's own
`K_stack=0.4286` is already the *best* (least-declining) point in the whole sweep, and every
smaller value (i.e., every "softer stack" in the naive reading of the mechanism) makes the
trend worse.

**Why, mechanistically (worth understanding, not just observing)**: during the pre-
transition regime, `c_pore_normalised(ε_struct(t))` sits essentially flat near 1 for
*any* `K_stack` value — `K_stack` only sets the height of the resulting plateau,
`f0_eff = 1/(1+K_stack·1) = 1/(1+K_stack)`, not whether that plateau is flat. So
`buffered_amplitude(t) ≈ f0_eff × unbuffered_amplitude(t)` throughout this whole regime,
for every `K_stack` — a **constant multiplier** applied to the *already-declining*
unbuffered signal, exactly the same structural argument A.0 originally used to rule out
`f0`/width recalibration for the baseline recipe. A smaller `K_stack` just means a *larger*
`f0_eff`, i.e. a multiplier closer to 1 — which necessarily pulls the buffered curve closer
to the (more steeply declining) unbuffered curve, not further from it. There is no value of
`K_stack` alone, however extreme, that can turn a multiplicative rescale of a declining
signal into an increasing one — this is now confirmed for *both* free constants in the
"physical" transmitted-fraction form (`f0`/width via A.0's original sweep, and now
`K_stack` independently via this one), which makes it a general property of the
single-mechanism pore-buffering partition, not an artefact of which constant was swept.

**What this means for Part C.4 and the pressure-coupling theory.** The *theory* in A.0g
(pressure is a mechanical, not kinetic, parameter — the reason low pressure alone shows a
qualitatively different shape) is not undermined by this result. What's falsified is the
*specific mechanism*: pressure acting *only* through a constant rescale of the transmitted-
fraction partition cannot work, for the same structural reason recalibrating `f0`/width
couldn't. A viable pressure-coupling mechanism needs to change something that is not flat
during the pre-transition regime — two candidates, both already flagged as fallback options
elsewhere in this plan and now more clearly the *required* next step rather than an
optional Phase 2:
- **Feed `Stack pressure [Pa]` into the particle mechanical-stress boundary condition**
  (`σr(R)` in `base_mechanics.py`), not just the pore-buffering partition — this is exactly
  the "explicitly NOT included" Phase-2 extension C.4 originally deferred. Lower external
  pressure changing the particles' own stress state would change *cracking rate* itself
  (LAM), which changes the *underlying unbuffered amplitude's own trend* — not a multiplier
  on top of it — which is the only kind of change this result shows is actually capable of
  altering the pre-hump *shape*, not just its height.
- **Pressure-dependent `eps_transfer_width`/`eps_min_transfer`** (not `K_stack`): if lower
  pressure widens the transition band enough that `c_pore_normalised` starts measurably
  decreasing from much earlier in life (rather than staying pinned near 1 until close to
  `eps_min_transfer`), `f(t)` itself would no longer be flat pre-transition, which could in
  principle produce genuine curvature rather than a fixed rescale. Not yet tested — A.0's
  original width sweep varied width alongside `f0`-derived `K_stack`, not width alone at a
  fixed, independent `K_stack`, so this specific combination remains an open, cheap next
  test (reuse this script's plumbing, sweep `eps_transfer_width` instead of `K_stack`).

Either way, this re-links the pressure question back to the mechanical/LAM side of the
model rather than the pore-buffering partition in isolation — plausibly converging with the
lower-priority LLI-vs-LAM track after all (A.0's own text already speculated lower pressure
plausibly pushes toward LAM/cracking-dominance via mechanical stress), just entered through
pressure's effect on stress rather than a generic "majority pattern" rebalancing.

**Update — the `K_stack`-decoupling source change has been reverted.** The `"pore
buffering stack compliance"` model option, the `f_transmit_kmax` domain parameter, and
`k_stack_pressure_sweep_test.py` were a diagnostic patch to test the hypothesis above; now
that the test is done and answered (negatively), they've been removed from
`reaction_driven_porosity.py`/`lithium_ion_parameters.py`/`base_battery_model.py`/
`si_gr_expansion.py`, which are back to exactly their pre-test state (`k_cmax = (1-f0)/f0`,
no independent lever). `expansion_pattern_test/` has been cleaned up to just its result
images, moved to `expansion_pattern_test/pics/` — see A.0i below for what's next instead.

### A.0i Pivot: particle/stack rearrangement — a genuinely new, currently-unmodeled mechanism, not a recalibration

A.0h's result (no constant-rescale lever, on either free parameter of the existing
"physical" transmitted-fraction form, can flip the pre-hump sign) points toward a
conclusion worth stating plainly: **whatever is different about low pressure is very
likely something the model doesn't represent at all**, not a wrong calibration of
something it does. The two structural candidates A.0h ended on (stress-coupled cracking;
a wider/earlier pressure-dependent transition band) are still *extensions of the existing
buffering/cracking machinery* — both assume the right physics is present and just needs a
pressure-dependent input. A simpler, more direct explanation for why *only* low pressure
breaks the pattern: **loss of mechanical confinement lets particles/fragments physically
rearrange within the electrode as degradation proceeds** — a distinct, real, currently
entirely unrepresented volumetric process, not a modulation of anything already coded.

**Why this fits better than the constant-rescale extension did.** At normal/high external
pressure (baseline, and — critically — high-rate, high-temperature, and the narrow-SoC-
window condition, none of which vary pressure at all), particles are pressed together
tightly enough from BoL that there's no room for them to move: as SEI consumes pore
volume, the particles stay pinned in their original packing, so nothing macroscopic
happens until the pore network itself is essentially closed (matching the flat-then-sharp-
hump shape common to all four of those conditions, and matching the existing "physical"
transmitted-fraction form's own qualitative shape, however its constants are tuned — which
is exactly why A.0/A.0h could never fix this by tuning constants: the *model's structure*,
not just its calibration, assumes locked particle positions). Under insufficient external
confinement (34 kPa), particles have room to shift, rotate, or resettle as the surrounding
structure degrades — a real, physically distinct process (sometimes discussed in the
literature as stack decompaction or cyclic "ratcheting" under low preload) that could add
its own *growing*, largely irreversible contribution to measured thickness, superimposed
on top of (not instead of) the particle-swelling/pore-buffering signal the model already
computes. Unlike a rescaled `f(ε_struct)`, a genuine rearrangement/decompaction term would
naturally be small-to-negligible at BoL (structure still close to its as-built packing) and
grow progressively across life as cumulative cycling loosens the structure further — i.e.
it produces real curvature (a continuously *increasing* contribution), not a flat
multiplier on an already-declining signal, which is precisely the property A.0h showed the
existing mechanism structurally cannot provide.

**Scoping, at a "plan," not implementation, level** (matching this document's own
convention): this would need a genuinely new state, not a new parameter on an existing one
— something like an irreversible "structural looseness"/contact-loss variable that (a)
only accumulates meaningfully when local confinement is below some threshold (so it's
naturally near-zero for every fixed-pressure condition, `mult`- pressure-gated rather than
always-on), (b) grows monotonically with cumulative mechanical cycling (throughput, or
cumulative stress-reversal count, rather than being a pure function of `ε_struct` the way
the existing buffering partition is), and (c) contributes directly to cell thickness
alongside (additively to, not multiplicatively rescaling) the existing particle-swelling
and pore-buffering terms. This is a different kind of extension than C.4 originally
scoped — C.4 assumed the existing `f(ε_struct)` machinery just needed a pressure-dependent
constant; this instead needs a new, pressure-*gated* irreversible process alongside it.
Worth flagging explicitly: this is a bigger lift than C.4's original "smaller than it
looks" framing suggested, precisely because A.0h has now ruled out the cheap version.

**Relationship to the existing fallback list.** A.0's original fallback list already
included "stack-level mechanical relaxation over cycling" as a candidate, but framed it as
*independent of the fixed external pressure* (a general aging effect, not
pressure-specific). This pivot sharpens that into a more specific, better-motivated,
pressure-*gated* version of the same idea — worth reading as a refinement of that existing
fallback rather than a sixth new item.

**Update — A.0j below tested a cheaper, shape-only alternative before committing to this
bigger lift, and it works.** The particle-rearrangement mechanism above is still a good
*physical* story, but it turns out the existing single-mechanism pore-buffering framework
*can* reproduce the target shape after all — it just needed a different transition
functional form, not a new state variable. Read A.0j before treating "needs genuinely new
model structure" as settled.

### A.0j Test result: asymmetric (gradual-to-sharp) transition — the first lever in this whole investigation that actually flips the sign

A.0h ruled out constant-*height* rescales (`K_stack`, and earlier `f0`/width) because the
existing tanh/physical transmitted-fraction forms stay essentially flat throughout
pre-transition life — `f` only starts moving once `ε_struct` is already close to closure,
so `buffered_amplitude(t) ≈ (constant) × unbuffered_amplitude(t)` for nearly the
whole window, inheriting the declining unbuffered signal's sign regardless of the constant.
Tested directly whether a transition that's not flat pre-transition — one that rises
*continuously* from early life, only sharpening near closure — could do better, since
`buffered_amplitude(t) = f(t) × unbuffered_amplitude(t)` is then a genuine product of a
*rising* and a *declining* factor, which by the product rule can increase overall if `f`'s
relative growth rate outpaces the unbuffered decline rate — a possibility the flat forms
structurally excluded, not something A.0/A.0h actually tested.

Implemented a third `"pore buffering transition"` option, `"asymmetric"`
(`reaction_driven_porosity.py`'s `_transmitted_fraction_asymmetric`): progress is
normalised against the *full* BoL-to-closure porosity span (`ε_init − eps_min_transfer`),
not `eps_transfer_width`'s much narrower decay length, then raised to a shape exponent `n`
(`"{Domain} electrode pore buffering asymmetric transition exponent"`, a new parameter):
`n=1` is a pure linear ramp from `f0` at BoL, no "shoot up" at all; larger `n` concentrates
progressively more of the rise near closure. Default option stays `"tanh"`, fully additive,
no effect on any existing script (`asymmetric_transition_sweep_test.py`,
`expansion_pattern_test/`).

Tested on the same base-case recipe as A.0h (`SI_MULT=18`, unchanged degradation kinetics),
sweeping `n ∈ {1, 2, 3, 5, 10}`, 120-cycle diagnostic.

| `n` | Early-life slope [µm/A·h] |
|---|---|
| 1 (linear ramp) | **+0.00332** |
| 2 | **+0.00121** |
| 3 | −0.00071 |
| 5 | −0.00263 |
| 10 | −0.00357 |
| Unbuffered reference | −0.00530 |

**Result: `n=1` and `n=2` both cross into genuinely increasing territory** — the sign flips
between `n=2` and `n=3`. Visually
(`expansion_pattern_test/pics/asymmetric_transition_sweep_test_result.png`), `n=1`'s
amplitude rises essentially monotonically from ~12 µm at BoL to its ~14.5 µm peak with no
declining region at all, a good qualitative match to the low-pressure condition's
continuously-increasing, no-flat-plateau shape from A.0g — a genuinely different outcome
from every other lever tried in this investigation (`f0`/width in A.0, `SI_MULT`/`LAM_MULT`
in A.0a/b, `K_stack` in A.0h), none of which ever got past a less-steep decline. `n=3`
onward reverts to a U-shape (decline then hump), converging toward the existing forms'
behaviour as `n` grows, as expected. LLI/LAM at the final cycle are essentially unchanged
across the sweep (17.0–17.1% LAM, 18–20% LLI) — confirms the shape change isn't
accidentally altering the underlying degradation kinetics, only the partition.

**What this means.** The single-mechanism pore-buffering framework *can* produce an
increasing pre-hump trend — A.0's original conclusion ("no physically reasonable parameter
choice... produces an increasing pre-hump trend") was correct for the two forms tested
there (tanh, physical) but not a property of the framework in general; a third functional
form closes the gap. This is a materially cheaper fix than A.0i's new-state-variable
proposal, *if* the exponent `n` can be given a physically coherent story — and it maps
onto A.0g/A.0i's own confinement reasoning unusually well: `n` small (gradual,
early-starting transition) reads naturally as "unconfined, particles free to resettle
throughout life" (low pressure), while `n` large (flat-then-late-sharp, i.e. converging
toward the existing tanh/physical forms) reads as "confined, particles locked until pores
are nearly closed" (normal/high pressure, matching why baseline/high-rate/high-T/SoC-window
all already show the model's default shape). That reframes `n` itself as a plausible,
much cheaper `n(P_stack)` pressure-coupling lever — replacing C.4's original
`K_stack(P_stack)` idea (falsified by A.0h) with an exponent-based one instead.

**Caveats, stated plainly rather than declaring victory**:
- This is a demonstration that the *sign* can flip with the right shape, on the *base-case*
  recipe only — it is not yet a calibrated low-pressure model. No `n(P_stack)` mapping
  exists yet; `n` was swept as a free constant, not tied to 34 kPa vs. 103 kPa by anything
  but analogy.
- Only tested at 120 cycles (matching this investigation's short-diagnostic convention) —
  A.0e's extended-run result (SI_MULT reduction's "flattening" lead not surviving to 500
  cycles) is a direct reminder that short-window trends don't always hold up; not yet
  checked whether `n=1`/`n=2`'s rise continues cleanly into a proper hump at longer runs,
  or whether it saturates/reverses awkwardly given `f` still hard-caps at 1.
- Doesn't resolve which of this (existing-framework, new functional form) or A.0i
  (genuinely new state variable) is the *physically correct* explanation — both remain live
  until backed by more than the pre-hump shape alone (e.g. does either one also get the
  post-knee behaviour, or the muted-hump high-temperature condition, right without further
  tuning). This result changes the *cost* comparison (this is now clearly cheaper), not
  necessarily which one is *true*.
- **Visually, does it produce a proper hump?** Not obviously — worth checking directly,
  which is what A.0k below does, by asking a related, cheaper question first: can the
  *existing*, already-validated "physical" form (whose hump shape is well understood) be
  widened to do the same job, without introducing a new functional form at all?

### A.0k Test result: widening the existing "physical" form's `eps_transfer_width` — gets closer to zero but never crosses, and the hump vanishes before it would

A.0j's power-law `"asymmetric"` option flips the sign, but raised a fair question: does it
still produce a *proper* hump, or just a broad, gradually-rounded peak? Rather than
scrutinise the power-law's shape further, tested a related, cheaper idea first: the
*existing* `"physical"` option's own `eps_transfer_width` is already a free parameter — A.0's
original sweep only tried it up to `0.04`, narrow compared to the ~0.17 structural-porosity
range (`eps_init=0.25` to `eps_min_transfer=0.08`) actually traversed over life. No source
changes needed — `physical_wide_width_sweep_test.py` just reruns the unmodified `"physical"`
option at much larger widths.

| `eps_transfer_width` | Early-life slope [µm/A·h] | Hump prominence (peak − pre-peak min) |
|---|---|---|
| 0.04 (A.0's original max) | −0.00204 | 3.2 µm |
| 0.08 | −0.00056 | 2.3 µm |
| **0.12 (best/closest to zero)** | **−0.00029** | 1.6 µm |
| 0.17 (≈ full structural range) | −0.00048 | 1.0 µm |
| 0.25 | −0.00104 | 0.0 µm (hump gone) |
| 0.40 | −0.00194 | 0.0 µm (hump gone) |

**Result: never crosses zero, and there's a genuine trade-off, not just a tuning miss.**
The slope improves from `width=0.04` to a best point around `width=0.12`, then gets *worse*
again as width increases further — non-monotonic, with the best point still negative. Hump
prominence, meanwhile, declines *monotonically* with width and vanishes entirely by
`width=0.25`, well before the slope could ever turn positive even if the trend continued
(`physical_wide_width_sweep_test_result.png` shows this directly: the wide-width curves
become smooth, featureless declining-to-flat shoulders with no real peak left at all).

**Why**: `_transmitted_fraction_physical`'s `c_pore_normalised = 1 − exp(−headroom/width)`
has only *one* shape parameter doing two jobs at once — how gradual the early rise is, and
how sharp the late closure is — and turning that single knob trades one against the other.
Widening it enough to soften the early plateau necessarily also softens the late transition,
because both are governed by the same exponential decay length. A.0j's `"asymmetric"`
option avoids this because it uses *two* decoupled parameters instead of one: a wide
BoL-to-closure span sets the gradual part, and a separate exponent `n` independently sets
how sharp the finish is. Concretely, A.0j's `n=1` case beats every point in this sweep on
*both* axes simultaneously — slope +0.00332 (vs. this sweep's best of −0.00029, still
negative) and hump prominence ~2.5 µm (comparable to this sweep's `width=0.08` point, which
had no early rise at all).

**Conclusion (width alone)**: widening `eps_transfer_width` on its own doesn't work — not a
failed attempt at the same idea A.0j already solved, but a genuine confirmation that
varying only one shape parameter forces the early-gradual/late-sharp trade-off. **Update —
A.0l below found that adding `f0` as a second free parameter, still within this exact same
"physical" form, resolves the trade-off after all.** Read A.0l before treating "the
physical form can't do this" as the final word — it turns out it can, once two parameters
move together rather than one.

### A.0l Test result: combining lower `f0` with the A.0k widths — every combination crosses the sign, no source changes at all

A.0k varied `eps_transfer_width` alone (at the recipe's existing `f0=0.7`) and never
crossed zero. Natural next question: does adding `f0` as a second free dimension, still
within the *exact same*, unmodified `"physical"` form, escape the single-parameter
trade-off A.0k found? Lower `f0` makes `k_cmax=(1-f0)/f0` larger — e.g. `f0=0.4` gives
`k_cmax=1.5` vs. `f0=0.7`'s `0.4286` — which makes `f` more sensitive to whatever gradual
change in `c_pore_normalised` a widened `width` already produces, potentially amplifying a
near-miss into a genuine sign flip. Tested directly
(`physical_f0_width_grid_test.py`, `f0=0.7` skipped — already known from A.0h/A.0j/A.0k):
`f0 ∈ {0.6, 0.5, 0.4} × eps_transfer_width ∈ {0.12, 0.17}` (A.0k's two best widths), same
base-case recipe, 120-cycle diagnostic.

| width | f0 | Early-life slope [µm/A·h] | Hump rise [µm] | BoL amplitude [µm] |
|---|---|---|---|---|
| 0.12 | 0.6 | +0.00119 | 3.14 | 11.21 |
| 0.12 | 0.5 | +0.00247 | 5.00 | 9.25 |
| 0.12 | 0.4 | +0.00345 | 7.02 | 7.17 |
| 0.17 | 0.6 | +0.00110 | 2.35 | 12.02 |
| 0.17 | 0.5 | +0.00258 | 4.14 | 10.15 |
| 0.17 | 0.4 | **+0.00385** | 6.12 | 8.11 |

**Result: every single combination crosses into positive territory**, with both slope and
hump prominence improving *monotonically* as `f0` drops — no trace of A.0k's trade-off.
The best point (`width=0.17, f0=0.4`) beats even A.0j's best `"asymmetric"` case on slope
(+0.00385 vs. +0.00332) while keeping a substantial, well-shaped hump (6.1 µm rise, using
the physical form's genuine exponential closure, not a bare power law) —
`physical_f0_width_grid_test_result.png` shows all six curves rising smoothly from BoL into
a clear, shared peak around 850–900 A·h, with no sign of the featureless flattening A.0k's
wide-width-alone cases showed.

**A real tension worth stating plainly, not glossing over**: this *inverts* the direction
A.0g originally predicted. A.0g's argument was low pressure → less mechanical confinement →
*higher* effective `f0` (less BoL buffering headroom, since unconfined particles have less
room to be pre-compressed into place before any swelling shows up). What actually works
here is the opposite: *lower* `f0`, i.e. *more* buffering at BoL, not less — visible
directly in the BoL amplitude column, which nearly halves (11.2 → 7.2 µm) as `f0` drops
from 0.6 to 0.4. A more consistent story for what's actually happening, inferred from what
worked rather than derived first: a softer, more compliant pore network absorbs *more*
swelling per unit degradation early in life (hence low `f0`), but spends that absorption
capacity continuously across a wide band (wide `width`) rather than holding it in reserve
for a late, sharp event — so its finite capacity to keep absorbing gets used up steadily,
showing up as a smoothly *rising* transmitted fraction well before near-total pore closure,
rather than A.0g's "stays locked, then flips late" framing. Both `f0` and `width` moving
together is what makes this work — A.0k already showed `width` alone can't do it, and by
the same "single lever, one trade-off" logic demonstrated there, `f0` alone (at the
recipe's default `width=0.01`) would very likely fail too (not separately tested — `f0`
was always swept together with a widened `width` here).

**Conclusion**: the *existing*, already-validated "physical" transmitted-fraction form
*can* produce a genuinely increasing pre-hump trend with a proper hump — no new functional
form (A.0j's `"asymmetric"` option) and no new mechanism (A.0i's rearrangement idea) are
strictly necessary to get the *shape* right. What's now needed to turn this into an actual
low-pressure model, same open items as A.0j left: a `(f0, width)(P_stack)` mapping (only
two data points calibrated so far — this grid, not yet tied to 34 kPa vs. 103 kPa by
anything but the direction that works), and a check that this holds up beyond the
120-cycle/base-case-recipe diagnostic used throughout A.0h–A.0l.

**Update — re-run to full life (50% SoH floor, not just 120 cycles); result holds.**
Extended `physical_f0_width_grid_test.py` to run each of the 6 combinations (plus the
unbuffered reference) to a 50% SoH floor rather than a fixed 120-cycle window, tracking
SoH and the *measured* (observable) transfer ratio — `buffered amplitude / unbuffered
amplitude at matching throughput`, distinct from the internal `f(ε_struct)` variable —
instead of the earlier slope/hump-prominence summary panels. (Caught and fixed a real bug
first: the initial SoH-termination check indexed `sol.cycles[0]`/`[-1]` directly and always
returned `cap=0`, so every run silently burned through a 400-cycle safety cap instead of
stopping at 50% — fixed by reusing the same per-cycle extraction loop already proven
correct elsewhere in the script, rather than a new ad-hoc access pattern.)

All 6 combinations reach 50% SoH at 140–160 cycles (~1,050–1,100 A·h throughput) — a
consistent, clean rise from BoL through a shared peak around 850–900 A·h (matching the
knee), then a shared post-knee decline, for every `(f0, width)` point, confirming A.0l's
result isn't a 120-cycle artefact (recall A.0e's warning: the earlier `SI_MULT` reduction's
apparent improvement did *not* survive to a longer run — this one does).
**SoH trajectories are visually indistinguishable across all 6 buffered cases and the
unbuffered reference** — buffering redistributes swelling into thickness vs. pores, but
doesn't materially alter the degradation-kinetics clock, as expected.

**One genuine, honest nuance from the measured-transfer-ratio panel**: it doesn't just
rise monotonically to 1 — it overshoots to ~1.05–1.15 right around the knee before falling
back below 1 post-knee. That means the buffered case's hump is briefly *taller* than the
same-throughput point on the unbuffered reference curve, which looks odd for a mechanism
that's supposed to only ever route swelling *away* from measured thickness. It's a real,
self-consistent model feature, not a bug: buffering changes reported porosity
(`eps_k = eps_struct − dv_buffered`), which feeds back into the transport/resistance
equations, so the buffered and unbuffered runs aren't running identical electrochemistry
underneath — they can diverge slightly in exactly *when* the knee's stress/current
concentration peaks, not just in how much of a given swelling event reaches the surface.
Worth keeping in mind for future work on this (e.g. if a precise measured-`k`-vs-1
crossover point ever becomes load-bearing for a claim), but doesn't affect anything
concluded so far.

### A.0m Direct comparison: does A.0l's model result actually look like the real low-pressure curve?

Worth checking directly against `fig2_population_lead_time_v2.png`'s numbers (A.0g) rather
than leaving A.0l's positive result as an abstract "the sign flips" finding.

| | BoL amplitude | Peak amplitude | Peak/BoL ratio | Shape |
|---|---|---|---|---|
| **Real low-pressure data** (34 kPa) | ~53 µm | ~78 µm (at knee, EFC 158) | ~1.47× | continuous rise from cycle 1, no flat plateau |
| **Model, best A.0l point** (`width=0.17, f0=0.4`) | 8.11 µm | ~14.2 µm (at the shared knee/hump) | ~1.75× | continuous rise from BoL, no flat plateau |
| Real baseline (103 kPa), for contrast | ~52 µm | ~64 µm (at knee, EFC 145) | ~1.23× | flat, then a hump concentrated at the knee |

**The qualitative shape match is good**: both the real low-pressure curve and the model's
best A.0l point rise continuously from BoL with no flat lead-in, unlike baseline's
flat-then-hump — exactly the distinguishing feature A.0g identified as the one genuine
outlier pattern. The proportional rise (peak/BoL) is in the same ballpark too (1.47× real
vs. 1.75× model) — the model is somewhat *more* aggressive than the real data, which is a
reasonable place to land for a first parameter choice pulled from a coarse 3×2 grid (A.0l),
not a proper fit.

**What this comparison is *not***: a validated, calibrated pressure model. The absolute
magnitude scale differs by roughly 6× (tens of µm for the real cell vs. single-digit-to-
low-teens µm for this diagnostic recipe, which was never calibrated to this specific paper's
cell geometry/capacity) — so this is a shape-only comparison, matching this whole
investigation's convention of using un-rescaled diagnostic runs. And as A.0l already noted,
`(f0=0.4, width=0.17)` isn't tied to "34 kPa" by anything but the *direction* that works —
there is still no `(f0, width)(P_stack)` mapping, only a demonstration that a point in that
parameter space reproduces the right qualitative signature. Still, this is the first time
in the whole investigation (A.0a–A.0m) that a model result has been checked directly
against real numbers from the source figure, rather than just its own internal sign/shape
— worth treating as a meaningfully positive, if still preliminary, result.

### A.1 What the figure shows

**Superseded by A.0g above for the majority/minority framing — kept here as the original,
per-condition reading (still useful for the individual condition-level hypotheses in A.2),
but "increasing is the majority pattern" should be read as corrected by A.0g.**

Four conditions, each varying one parameter from baseline, all four eventually showing
the same qualitative **hump**: reversible expansion rises to a peak right around the
capacity knee, then declines. What differs is the shape **before** the hump:

| Condition | Params | Lead EFC → knee EFC | Pre-hump trend (my reading) |
|---|---|---|---|
| Baseline | 25 °C, 0.33C, 103 kPa | 127 → 145 | roughly **flat** (~53–55 µm), then a fairly sharp rise into the hump |
| Low pressure | 25 °C, 0.33C, 34 kPa | 129 → 155 | **increasing** almost continuously from cycle 1 (~48 → ~80 µm) — no distinct flat region; the "pre-hump" and "hump rise" essentially merge |
| High C-rate | 25 °C, 2C, 103 kPa | 88 → 133 | much shorter window (fewest EFC to knee of all four); mild rise, hard to distinguish from "flat" given how little pre-hump data exists |
| High temperature | 45 °C, 0.33C, 103 kPa | 292 → 319 | by far the longest pre-hump window; trend looks **decreasing** (~57 → ~47 µm), and the eventual hump is comparatively muted |

(This is my reading of the plotted curves — worth a quick visual double-check against the
source slide before treating the flat/increasing/decreasing assignment as settled, but
the qualitative claim — three distinct pre-hump shapes across the four conditions — is
robust either way and is the important part.)

### A.2 Hypotheses, grounded in the pore-buffering mechanism already implemented

Recall the mechanism (`reaction_driven_porosity.py`): particle-level swelling is
partitioned between pore-volume buffering and measurable cell/electrode thickness change
via a transmitted fraction `f(ε_struct)`, bounded between `f0` (BoL plateau, most
swelling absorbed by pores) and 1 (pores closed, all swelling transmitted). `ε_struct`
declines monotonically over life as SEI/plating/cracks consume pore volume. The **hump
itself** is well understood from earlier sensitivity work this session (`CHANGES.md`):
once `f → 1` (saturated), amplitude = the underlying *unbuffered* particle-swelling
amplitude exactly — a quantity that's typically itself declining over life (as LAM
reduces active material and hence per-particle current density/swing) — so a transition
from "flat/buffered" to "tracking a declining unbuffered signal" naturally produces a
peak right at the point buffering saturates. That part of the mechanism is condition-
agnostic and already explains the *hump* uniformly across all four conditions. What
differs here is the **pre-hump regime** — i.e. what `f(ε_struct(t))` looks like *before*
saturation, which depends on where the cell sits relative to the transition band from
the start, and on the FIRST derivative of the underlying swelling trend before buffering
saturates.

**Baseline — flat**: `ε_struct` starts comfortably above `eps_min_transfer +
eps_transfer_width`, so `f` sits at (or very near) the `f0` plateau for most of life — in
principle the partition does near-maximal, near-constant buffering, which *would* produce
a roughly flat measured trend if the underlying unbuffered amplitude were itself flat.
**Caveat, per the cross-check in A.0**: at the model's *current* calibration this doesn't
actually hold — the underlying unbuffered amplitude is declining strongly enough over this
regime that our simulated baseline is decline-then-hump, not flat-then-hump. Worse, per
A.0's follow-up sweep, this isn't fixable by retuning the transition band either — the
mechanism structurally can't flatten it out for any physically reasonable band choice, so
"flat" (or the more common "increasing") needs one of A.0's alternative degradation
mechanisms added, not a recalibration of this one.

**Low pressure — increasing continuously, no flat region**: hypothesis: lower external
stack pressure means less mechanical preload compressing the electrode/separator stack.
Physically this plausibly does one or both of: (a) leaves the cell starting closer to
the `f0` plateau's *edge* rather than deep in it — i.e. less "buffering headroom" from
BoL, so smaller changes in `ε_struct` produce visible changes in `f` immediately, instead
of `f` staying pinned at `f0` for a long flat stretch first; (b) makes the underlying pore
network itself more compliant/less rigid at a given `ε_struct` (a real stack-mechanics
effect — see Part C.4), so the transmitted fraction rises smoothly and monotonically
with declining `ε_struct` from very early in life, rather than showing a sharp late
transition. Either mechanism predicts exactly what's observed: a continuously rising
trend with no distinct flat pre-hump plateau. **This is currently NOT modelled** — see
A.3.

**High C-rate — short window, mild trend**: two effects plausibly compound. First,
timing: higher current density accelerates SEI growth and (via higher diffusion-induced
particle stress) cracking, so `ε_struct` declines faster per equivalent full cycle,
reaching the transition band much sooner (matching the much lower lead EFC, 88, the
earliest of the four). This part — earlier *timing* — is plausibly within reach of the
existing rate-dependent kinetics (Part C.2). Second, shape: because the pre-hump window
is so short in EFC terms, there's little time for a genuinely flat plateau to establish
before the transition band is reached, so almost any nonzero rate of `ε_struct` decline
would show up as a "mild rise" over such a short window even under the SAME buffering
physics as baseline — i.e. part of this condition's apparent shape difference could be a
*windowing* artefact of comparing a short pre-hump interval to a long one, not evidence
of a qualitatively different mechanism. Worth testing directly (Part C.2).

**High temperature — decreasing, muted hump**: the most striking and least obviously
explained pattern, since it's also the most *delayed* knee (lead EFC 292, by far the
latest) despite temperature normally accelerating SEI kinetics via Arrhenius terms —
naively predicting an *earlier* knee, the opposite of what's observed. Hypothesis (more
speculative, flagged explicitly as needing testing, not assumed): higher temperature
increases Li-ion and electrolyte diffusivity, which lowers concentration-polarisation-
driven mechanical stress in the particles for a given current — i.e. even though reaction
kinetics (SEI growth) speed up with temperature, the *mechanical driving force* for
cracking (which depends on concentration gradients, not just reaction rate) could
simultaneously go down, and if the mechanical/cracking pathway dominates `ε_struct`'s
practical decline rate more than pure SEI kinetics do, the net effect could be a slower
overall approach to the transition band despite faster underlying chemistry — consistent
with the much later knee. Separately, the *decreasing* early trend and *muted* hump could
reflect a genuinely different, non-buffering mechanism entirely: a "settling-in" effect
where a fresher, rougher electrode/SEI network at BoL is more compliant and gradually
consolidates over the first many cycles (an effect that would show up as a
temperature-accelerated *early* relaxation, superimposed on top of, and initially
dominating, the buffering mechanism). This second explanation is not currently represented
anywhere in the model (buffering only ever describes ε_struct *decline*, never an initial
compliance-relaxation transient) and would need new physics if it turns out to be the
right explanation, not just new temperature dependence on existing physics.

### A.3 Can the current model capture each pattern? Gap analysis

**Superseded by A.0g–A.0m above — kept here as the original per-condition record, but read
the update table immediately below first.** This table predates the premise correction
(A.0g: declining-then-hump is the *majority* pattern, matching the model's default, not a
gap needing a fix) and predates A.0h–A.0m's finding that the "Low pressure" gap is now
partially closed (shape reproducible with the *existing* model, via a specific `(f0,
width)` combination — A.0l/A.0m), not the structurally-unfixable gap this table describes.

| Condition | Capturable today? | Why / what's missing |
|---|---|---|
| Baseline (flat, and — per A.0 — really "increasing" is the more common target shape across conditions) | **No — recalibration tested directly and ruled out; needs the dominant degradation pathway to shift** | A direct width/`f0` sweep (A.0 update) confirms the pore-buffering mechanism cannot produce an increasing (or even genuinely flat) pre-hump trend for any physically reasonable parameter choice — the pre-band regime is structurally governed by the declining unbuffered amplitude regardless of transition-band tuning. Leading hypothesis (A.0): the recipe's `SI_MULT=18` makes it strongly LLI-weighted, and LLI shrinks the usable stoichiometry window (declining amplitude) while LAM/cracking grows it for surviving particles (increasing amplitude) — cheapest test is rebalancing existing `SI_MULT`/`GR_DIV`/LAM-rate constants, zero new model code. Lithium plating is the next-cheapest fallback if that doesn't close the gap. |
| Low pressure (increasing) | **No — confirmed gap on two fronts now** | (1) No external/applied stack pressure is modelled as an input anywhere in this codebase (verified by grep across `src/pybamm` and the input parameter sets — zero matches for any `"Stack pressure [Pa]"`-style parameter; the only "pressure" in the whole tree is *internal* electrolyte fluid pressure in unrelated lead-acid convection submodels inherited from upstream PyBaMM). This is not news to the project — `pore_buffering_implementation_plan.md` (§"What NOT to touch", §2.6) already explicitly deferred exactly this. (2) Per the recalibration result above, even a pressure-coupled `K_stack(P_stack)` would need to act through, or alongside, one of the alternative mechanisms — pressure-coupling the existing single-mechanism `f(ε_struct)` alone inherits the same structural limitation. Part C.4 below picks this back up concretely. |
| High C-rate (short window, mild) | **Partially, and also blocked on the same limitation** | Timing (earlier knee at higher rate) is plausible via already-rate-dependent SEI/cracking kinetics — worth a direct check (Part C.2), not yet run. But per the recalibration result, whatever pre-hump *shape* baseline ends up needing an alternative mechanism to produce, high C-rate will need the same mechanism, just compressed into a shorter window. |
| High temperature (decreasing, muted) | **Already qualitatively matched — but for reasons that need auditing** | Per A.0, this is the ONE condition the model's current, *unmodified* calibration already reproduces the right qualitative shape for (decline-then-hump) — worth confirming this isn't a coincidence (i.e. check it's for the *right* physical reasons, not just because both happen to decline). Standard Arrhenius activation-energy terms exist on SEI/diffusivity parameters already, but naively predict the *opposite* timing of what's observed (faster, not slower, degradation at higher T) — reproducing the observed *timing* reversal needs either (a) confirming the diffusion-driven-stress-reduction hypothesis actually dominates in this model's own equations at typical parameter values, or (b) new physics (a settling-in/compliance-relaxation term) not currently represented at all. Flagged as a genuine open research question — see Part C.3. |

**A.3 update, reflecting A.0g–A.0m:**

| Condition | Capturable today? | Status |
|---|---|---|
| Baseline, high C-rate, high temperature, narrow SoC window | **Yes, qualitatively, already** | A.0g: all four share the model's existing decline-then-hump default shape, just at different knee timings (already rate/kinetics-sensitive). High-temperature's *timing reversal* specifically remains a genuine open question (Part C.3) — shape yes, timing mechanism no. |
| **Low pressure** | **Shape: yes, with a specific parameter combination. Calibration: no, not yet.** | A.0l found `(f0, width)` pairs on the *existing*, unmodified `"physical"` transition (e.g. `f0=0.4, width=0.17`) reproduce a genuinely increasing pre-hump trend with a proper hump, holding to 50% SoH. A.0m checked this directly against the real low-pressure numbers (`fig2_population_lead_time_v2.png`): same qualitative shape (continuous rise, no flat plateau) and a comparable proportional rise (model 1.75× vs. real 1.47×). Still missing: any `(f0, width)(P_stack)` mapping tying this to an actual pressure value, and confirmation this is the physically *correct* mechanism rather than a shape-matching coincidence (A.0i's rearrangement idea remains an live, untested alternative explanation for the same shape). |
| Panel d (combined high-rate + narrow-SoC-window, new in the updated figure) | **Not yet assessed** | Shows no hump at all within ~800 EFC — not yet compared against any model run; flagged in A.0g as an open item requiring its own investigation before folding into this table. |

---

## Part B — Test matrix (page 6), single operating condition

**A detailed, elaborated plan for this now lives in
`degradation_test_matrix/degradation_test_matrix_plan.md`** — concrete parameter values per
row, the dry-out/composite coupling groundwork needed for rows 4–5 (building on B.3 below),
a buffered-vs-unbuffered pairing per row, and the 50%-SoH run convention. B.1–B.5 below
remain the original, higher-level record; read the dedicated plan for anything
implementation-facing.

### B.1 The slide's own "Missing Model" is the pore-buffering submodel we already built

Page 6's flow (Si fatigue cracking → excessive SEI generation → pore clogging) maps onto
existing PyBaMM submodels for the first two steps (silicon fatigue/cracking model, SEI
generation model), but explicitly calls out a **"Missing Model": "Particle to cell
expansion model considering porosity reduction — essential to create expansion hump"**
for the third. That's exactly the pore-buffering submodel (`reaction_driven_porosity.py`,
`"pore buffering"` option) already implemented and validated in this project. This means:
the test matrix's real purpose, in this codebase's terms, is to show pore-buffering's
effect *conditioned on* which of the three named degradation mechanisms are active — pore
buffering is the constant "lens" the matrix is viewed through, not a 4th independent
toggle.

### B.2 Row-by-row mapping

`"pore buffering"` requires `"SEI porosity change"` (or `"lithium plating porosity
change"`, unused here) to be `"true"` — enforced by an existing `OptionError` guard
(`base_battery_model.py:766-775`). So pore buffering is **only meaningfully engageable**
in rows where "SEI porosity loss" is On; in the rows where it's Off, `ε_struct` never
declines at all, so `f` stays pinned at `f0` and pore buffering — on or off — is a
no-op by construction. That's not a gap to work around, it's the *expected*, physically
correct outcome for those two rows (no porosity change → no buffering transition → no
hump possible, trivially) and should be reported as such, not treated as a limitation.

| Row | Si cracking | SEI porosity loss | Dry-out | Pore buffering | Concrete PyBaMM mapping |
|---|---|---|---|---|---|
| Healthy baseline | Off | Off | Off | n/a (moot) | `"SEI porosity change": "false"`; `"Primary/Secondary: Negative electrode cracking rate"` → 0 (or leave `"particle mechanics"` at `"none"`/`"swelling only"` for the negative electrode) |
| Silicon fatigue only | **On** | Off | Off | n/a (moot) | Same as above but Si's (`Secondary:`) cracking rate at its normal value; Gr's (`Primary:`) cracking rate → 0 to cleanly isolate "silicon fatigue" per the slide's own framing, rather than conflating it with graphite's separate (already-suppressed) cracking |
| Pore clogging, Si intact | Off | **On** | Off | **On** | `"SEI porosity change": "true"`, `"pore buffering": "true"`; both phases' cracking rates → 0 |
| Dry-out, Si intact | Off | **On** | **On** | **On** | as above, plus the `ec_dryout_wrapper.py` batch-restart loop engaged — see B.3, open question |
| Fully coupled | **On** | **On** | **On** | **On** | all of the above combined, cracking rates at their normal (tuned) values |

**Important nuance on "Si cracking On/Off"**: `"particle mechanics"` is a *domain*-level
option (a 2-tuple over negative/positive electrode, confirmed at
`base_battery_model.py:215-220` — "A 2-tuple can be provided for different behaviour in
negative and positive electrodes" — **not** phase-resolved). It cannot independently turn
cracking on for Si but off for Gr within the same composite negative electrode. The
correct lever for phase-resolved on/off is the phase-prefixed cracking-rate parameter
(`"Primary: Negative electrode cracking rate"` / `"Secondary: Negative electrode cracking
rate"`), set to (near) zero for "off" while leaving `"particle mechanics": "swelling and
cracking"` enabled at the domain level throughout, so the same model-option skeleton
applies to every row and only rate constants change — this also matches the established
convention in this project's own scripts (`SI_CRACK_RATE_MULT`/`GR_CRACK_RATE_MULT`).

### B.3 Open design question: dry-out needs a different SEI submodel than the tuned recipe uses

`ec_dryout_wrapper.py` requires `"SEI": "solvent-diffusion limited"` (the only PyBaMM SEI
submodel whose current-density expression depends on a bulk solvent-concentration
parameter the wrapper can meaningfully update — see `ec_dryout/implementation_plan.md`
§3.1). The existing, tuned `si_gr_expansion` degradation recipe (crit_stress, LAM rates,
crack-rate multipliers, all calibrated together this session) uses `"SEI": "reaction
limited"` instead, which has no solvent-concentration term in vanilla PyBaMM at all.
Combining pore buffering + dry-out in the same simulation (rows 4–5) therefore needs one
of:
- **(a)** Recalibrate the whole degradation recipe under `"solvent-diffusion limited"`
  SEI kinetics instead of `"reaction limited"` — significant recalibration work (the
  entire crit_stress/LAM-rate tuning this project has done was against reaction-limited
  kinetics specifically), or
- **(b)** Investigate whether `"reaction limited"`'s own rate expression has *any* hook a
  solvent-depletion effect could attach to (unlikely as-is, since it's not
  solvent-transport-limited by construction — would need its own small design pass, not
  just a parameter change).

This needs a decision before rows 4–5 can actually be run; not resolved here per "no
implementation yet."

### B.4 Proposed script/output plan

Mirror the established convention (`test_pore_buffering/pore_buffering_degradation_test.py`
/ `pore_buffering_degradation_test_1000cyc.py`): one script per matrix row (or one
parameterised script looping over the 5 configurations, given they share almost all
scaffolding), each producing:
- A short run (~130–150 cycles, matching the existing short-cycle convention) and a
  1000-cycle timescale-stretched run (matching the existing long-cycle convention).
- The same core outputs already established: SoH/knee, cell-level within-cycle expansion
  amplitude, internal transfer ratio `k`, and (rows with pore buffering on) the
  observable `k = δ_cell/δ_particle`.
- A single 5-row (or fewer, since 1–2 are moot for buffering) comparison figure — same
  spirit as the existing unbuffered-vs-buffered `_comparison_*.png` plots — showing
  capacity and expansion-amplitude trajectories side by side across all runnable rows, at
  the SAME single operating condition throughout (no pressure/rate/temperature variation
  here — that's Part C).

### B.5 Success criteria per row (matching the slide's own framing)

Page 6 states the expectation plainly: **"Expansion (reversible and irreversible) has
different behavior"** across the matrix. Concretely, per row:
- *Healthy baseline*: flat capacity, flat/no-hump expansion (nothing to consume pore
  volume).
- *Silicon fatigue only*: capacity fade from LAM (cracking-driven active-material loss),
  but expansion should stay relatively flat/monotonic — **no hump**, since without SEI
  porosity loss there's no `ε_struct` decline to trigger the buffering transition. This is
  a directly falsifiable, specific prediction from the existing pore-buffering math and
  worth explicitly checking (would be a good sanity check that the mechanism is wired
  correctly, complementary to the sensitivity-sweep validation already done this session).
- *Pore clogging, Si intact*: capacity fade from SEI/LLI, and — for the first time in the
  matrix — an expansion hump should appear, driven purely by porosity-decline-triggered
  buffering saturation, with no cracking contribution.
- *Dry-out, Si intact*: as above, plus a further capacity-fade/knee-timing shift from
  electrolyte dry-out (per `ec_dryout/implementation_plan.md` §6's validated result:
  0%-reservoir dry-out measurably worsens capacity retention vs. a protected/high-reservoir
  case) — pending B.3's resolution.
- *Fully coupled*: the full picture, all mechanisms compounding — the closest analogue to
  the real experimental Si/Gr cell behaviour and the main point of comparison against
  page 2's actual data.

---

## Part C — Plan for capturing the operating-condition-dependent pattern differences

Separate from Part B: Part B holds operating condition fixed and varies mechanisms; this
part holds mechanisms fixed (fully-coupled, or whichever combination Part B settles on as
the best match) and asks what's needed to reproduce page 2's condition-*dependent* shape
differences. Staged by how much groundwork already exists.

**A detailed, elaborated plan for this now lives in
`conditions_test_matrix/conditions_test_matrix_plan.md`** — the concrete 5-condition
matrix (mirroring `fig2_population_lead_time_v2.png`'s panel a), a temperature-dependence
audit with specific findings (cracking-rate Arrhenius activation energy is hard-set to
zero in the upstream parameter set; SEI reaction kinetics have only indirect T-sensitivity
— both relevant to C.3 below), a knee/Lead-EFC detection methodology, and the EFC-conversion
convention needed to compare model output against the real figure's x-axis. Read the
dedicated plan for anything implementation-facing; C.1–C.5 below remain the original,
higher-level staging record.

### C.1 Staging rationale

**Reordered per A.0g.** The updated experimental picture shows baseline's declining
pre-hump shape is *not* an error — it matches the majority of conditions (baseline,
high-rate, high-T, narrow-SoC-window all share the same shape, just at different knee
timings). Only low pressure is a genuine outlier. That promotes **pressure (item 3 below,
detailed in C.4)** to the top of this list, ahead of rate/temperature and ahead of further
LLI-vs-LAM rebalancing work (item 0 below), which was motivated by the now-revised
"majority needs fixing" premise. Item 0 is kept, in its original form, as a lower-priority
parallel track — it produced some genuinely useful findings (A.0a) independent of which
premise turns out to matter more — but should not be treated as the leading candidate.

0. **(Lower priority than previously stated — see above.) Baseline pre-hump shape — shift
   the dominant degradation pathway, not a recalibration.** Originally staged first under
   the belief that the model's declining baseline needed fixing because most conditions
   increase; A.0g shows baseline's shape is likely already correct, so this is now a
   secondary investigation, not a blocker for the other items. **Recalibration was tried
   directly (a width × `f0` sweep) and ruled out** for producing an increasing trend at a
   *fixed, pressure-independent* band: every combination tested still declines from BoL
   before rising — the pre-band regime is structurally governed by the underlying declining
   unbuffered amplitude, independent of transition-band width or `f0`, so no amount of
   retuning `eps_transfer_width`/`eps_min_transfer`/`f0` within physically reasonable bounds
   gets there for the baseline recipe. (A.0g reframes this same "recalibrate the band" idea
   as the right *kind* of fix, just needing to be pressure-*dependent* rather than
   universal — see C.4.) Priority order for what to try instead, if pursued (A.0 has the
   full reasoning):
   0. **LLI-vs-LAM balance sweep** — cheapest by far, zero new model code: dial the
      existing recipe's `SI_MULT`/`GR_DIV` (SEI/LLI strength) down and/or its LAM
      proportional-rate constants and crack-rate multipliers up, and check whether the
      unbuffered particle amplitude's own trend flips from declining to increasing as
      LAM-driven concentration of lithium onto fewer surviving particles starts to
      outweigh LLI-driven shrinkage of the usable stoichiometry window.
   1. **Lithium plating** — cheap, existing PyBaMM submodel, needs a model-option change
      and a kinetics estimate, not new code.
   2. **Gas generation** — moderate; this session's dry-out work already derived the
      reaction stoichiometry that produces it (Eq. 1), just not its volume/venting
      consequences.
   3. **`E_s(φ)` stiffness evolution / growing cathode contribution / stack relaxation** —
      genuinely new model structure, try only if 0–2 don't close the gap.
   Whichever explains baseline correctly changes how conditions 1–3 below ought to be
   modelled too — notably, if the LLI-vs-LAM balance hypothesis is right, temperature and
   pressure/rate (items 1 and 3 below) are themselves plausible levers on that SAME
   balance (Arrhenius-accelerated SEI kinetics pushing high-T toward LLI-dominance;
   mechanical-stress-driven cracking pushing low-pressure/high-rate toward LAM-dominance),
   so this item and the condition-specific extensions below may turn out to be the same
   piece of work, not two separate ones — worth confirming with item 0's sweep before
   committing time to the rate/temperature runs below as if they were independent.
1. **Rate** — cheapest of the condition-specific extensions to test: current-density
   dependence is already baked into every relevant submodel (SEI growth, cracking). No new
   parameters needed, just running the model — ideally already extended per item 0 above —
   at a higher C-rate and checking the *shape*, not just the timing.
2. **Temperature** — moderate: Arrhenius terms already exist on the relevant rate
   parameters, but (per A.3) the observed trend is counter-intuitive and may need a
   genuinely new mechanism, not just temperature-scaling the existing one. Needs
   hypothesis-testing before committing to an implementation approach.
3. **Pressure — now the top priority overall, per A.0g.** It's the only condition
   confirmed to need genuinely new physics (vs. rate/temperature, which may be
   timing-only), but the *design* already exists on paper (per
   `pore_buffering_implementation_plan.md`'s deferred phase) and, as it turns out, is a
   surprisingly small structural extension of what's already implemented — see C.4.

### C.2 Rate: test first, minimal new work expected

Plan: run the existing `si_gr_expansion` + pore-buffering model at the baseline recipe
but with the ageing C-rate raised (e.g. to 2C matching the slide, discharge/charge legs
adjusted accordingly), same `eps_min_transfer`/`eps_transfer_width`/`f0` as baseline (no
retuning), and directly compare the resulting pre-hump `k(EFC)`/expansion-amplitude shape
against the flat baseline trajectory. Two clean, falsifiable outcomes:
- If the higher-rate run just produces an **earlier, otherwise identically-shaped**
  flat-then-hump curve → confirms the "short-window artefact" hypothesis from A.2, no new
  physics needed, timing alone (already rate-sensitive) explains the observation.
- If it instead needs a change to `f`'s shape itself to match — e.g. the transition band
  needs to be rate-dependent, not just reached sooner — that's a genuine new finding
  worth its own follow-up, not assumed here.

### C.3 Temperature: hypothesis-testing before implementation

Plan: before touching any code, use the *existing* Arrhenius-parameterised submodels to
directly check A.2's mechanical-stress hypothesis quantitatively — i.e. compute (from
existing, already-available diagnostic variables) how much particle-level concentration-
gradient-driven stress actually changes between 25 °C and 45 °C under this recipe's own
diffusivity/activation-energy parameters, and whether that swing is large enough to
plausibly offset the accelerated SEI kinetics at 45 °C. This is a **diagnostic run**, not
new model code — reuses existing stress/crack diagnostic variables already exposed by
`crack_propagation.py`. Two outcomes:
- If the existing stress-vs-temperature relationship in the model's own equations
  already points the right direction and is large enough → the reversed trend may already
  be *latent* in the model, just never checked at this specific temperature comparison;
  worth an actual full 45 °C run to confirm.
- If not → the "settling-in"/compliance-relaxation hypothesis (A.2's second explanation)
  becomes the more likely candidate, which would need genuinely new model structure (an
  early-life relaxation term, not currently represented anywhere) — a bigger, separate
  design task, not scoped further here.

### C.4 Pressure: extend the already-implemented "physical" transition, don't start from scratch

**Status update (A.0h–A.0m): the shape-reproduction half of this section is now
substantially done, ahead of the plan below** — worth reading before the original proposal
that follows. A direct `K_stack(P_stack)` extension (a single independent constant,
A.0h) was implemented, tested, and *failed* (made the trend worse, not better — wrong
direction entirely). What worked instead (A.0l/A.0m) was jointly retuning `f0` and
`eps_transfer_width`, two parameters the `"physical"` option *already exposes*, no source
change required — `(f0=0.4, width=0.17)` reproduces a genuinely increasing pre-hump trend
with a proper hump, and checks out qualitatively against the real low-pressure numbers
(A.0m). The source-code experiments along the way (`"pore buffering stack compliance"`,
the `"asymmetric"` transition option) were deliberately reverted once superseded — kept
here as a paper trail in A.0h–A.0k, not as live code. What the original proposal below
still gets right, and remains the real open work: turning "this `(f0, width)` point
reproduces the shape" into an actual `(f0, width)(P_stack)` function calibrated against
real pressure values, which is what steps 1–3 below still describe (now reframed around
`(f0, width)` rather than a single new `K_stack` parameter).

This is the best-understood of the three gaps, because the deferred design already
exists. `pore_buffering_implementation_plan.md` (§2.6, §3.4) references the *original*
source notes' mechanistic form for `f`:

```
f = 1 / (1/K_stack) / (1/K_stack + C_pore(ε))     [PDF Eq. 24, as documented]
  = 1 / (1 + K_stack * C_pore(ε))
```

Compare this against what's **already implemented** as `_transmitted_fraction_physical`
in `reaction_driven_porosity.py`:

```python
k_cmax = (1 - f0) / f0                      # a FIXED constant, calibrated from f0 only
c_pore_normalised = 1 - exp(-headroom/width)
f = 1 / (1 + k_cmax * c_pore_normalised)
```

These are **structurally the same functional form** — `k_cmax` is already playing exactly
the role the notes call `K_stack`. The only difference is that the current implementation
treats it as a fixed calibration constant (back-solved from the BoL transfer ratio `f0`
alone) rather than a genuine function of applied stack stiffness/pressure. This means the
pressure-coupling extension is **smaller than it looks** — not a new submodel from
scratch, but replacing one already-isolated constant with a function of a new input:

1. Add a new parameter, `"Stack pressure [Pa]"` (or per-cell, matching the notes'
   original naming) — a genuinely new BoL/operating-condition input, not derived from
   anything currently in the model.
2. Define `K_stack(P_stack)` — a new, small sub-relation (start simple: linear or
   power-law in `P_stack`, in the spirit of the deferred notes' `K_stack`/`E_s(φ)`
   framing) replacing the current fixed `k_cmax`. Direction of the relation should follow
   from stack mechanics: lower external pressure → softer/more compliant stack → *more*
   of a given `ε_struct` change transmits to measurable expansion at the same porosity,
   i.e. `K_stack` should move in whatever direction makes `f` rise faster (less
   BoL-buffering headroom) as `P_stack` decreases from the baseline 103 kPa toward the
   low-pressure condition's 34 kPa.
3. Calibrate the new relation's one or two free constants against exactly the two
   available data points (baseline vs. low-pressure pre-hump shape) — the same
   BoL-constants-only philosophy already used for `eps_min_transfer`/`eps_transfer_width`/
   `f0` (§3 of `pore_buffering_implementation_plan.md`): fit once from early-life behaviour,
   never re-fit against the full trajectory.
4. Leave `_transmitted_fraction_tanh` and the negative electrode's `eps_min_transfer`/
   `eps_transfer_width` untouched — this only touches the physical/compliance-ratio branch,
   consistent with `"pore buffering transition": "physical"` already being the actively-used
   setting in the current scripts.

**What this explicitly does NOT include** (per the original notes' own scoping, §2.6):
feeding `Stack pressure [Pa]` into the particle-stress boundary condition (`σr(R)` in
`base_mechanics.py`'s `_compute_stress_and_displacement`) — that's a mechanically deeper
coupling (external pressure directly altering internal particle stress state, hence
cracking rates too) than what's needed just to reproduce the *pre-hump expansion shape*
difference. Worth flagging as a natural Phase 2 if the simpler `K_stack(P_stack)` extension
alone doesn't fully reproduce the low-pressure condition's continuously-increasing trend.

**Validation plan**: once implemented, run baseline vs. `P_stack` = 34 kPa with the SAME
`eps_min_transfer`/`eps_transfer_width`/`f0` (only `K_stack(P_stack)` differs) and check
for the qualitative signature from A.1/A.2 — continuously increasing `k`/expansion from
early life, no flat pre-hump plateau — not just a shifted knee timing (which the existing
rate-dependent kinetics could produce for the wrong reason, i.e. a false positive worth
explicitly guarding against by checking the *shape*, not just the timing, matches).

### C.5 Relationship to Part B

Parts B and C are orthogonal and composable: Part B's 5-row mechanism matrix and Part C's
rate/temperature/pressure extensions can eventually be crossed (a small mechanism × 3
extra operating-condition axis) once both are independently validated, but that full
cross-product is out of scope for this planning pass — get each axis right in isolation
first.
