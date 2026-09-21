"""CELL064 degradation-rate fit: match the composite Si/graphite model's
capacity fade to CELL064's own real aging data.

Reuses the degradation coupling from the recipe that ACTUALLY produced the
knee/LLI/LAM pattern this project validated via DMA
(../../test_model/degradation_test_matrix/test_DMA/dma_baseline_run.py) --
not the plain "reaction limited" exchange-current-density recipe from
../../test_model/test_pore_buffering/pore_buffering_degradation_test_1000cyc.py
this script started from (an earlier iteration here tried that one first;
its rate constants transfer poorly to CELL064's real composition and its
"reaction limited"/"stress-driven LAM" combination hits a sharp, numerically
fragile porosity-floor cliff well before a real knee-like shape emerges --
see git history / calibration log for that dead end). The DMA-validated
recipe instead uses:

  - `"SEI": "ec reaction limited"` (Yang et al. 2017 reaction-kinetic law,
    tunable via "Primary:"/"Secondary: SEI kinetic rate constant [m.s-1]",
    NOT "SEI reaction exchange current density") -- genuinely accelerates
    into a knee (unlike "solvent-diffusion limited") while still coupling to
    an EC-concentration parameter, per
    ../../test_model/degradation_test_matrix/degradation_test_matrix_plan.md
    section 1a's joint calibration.
  - Pore buffering ON (`"physical"` transition, f0=0.7, width=0.01) --
    included for parity with the reference recipe, though (per
    ../../submodel_stability_audit.md and the submodel's own
    get_coupled_variables) it has no feedback into the electrochemistry, so
    it cannot itself affect capacity/LLI/LAM.
  - Critical stresses, LAM proportional terms, crack-rate multipliers and
    the negative-electrode porosity floor (0.01) kept at the reference
    recipe's own validated values, unchanged -- these produced a physically
    sane, numerically robust knee for si_gr_expansion's own composition, so
    they are the right starting point here too, not something to re-derive
    from scratch.

What differs from that reference recipe: base parameters are CELL064's own
real BoL fit (cell064_parameters.py), not si_gr_expansion's illustrative,
solved composition -- Si/graphite volume fractions, capacity, geometry and
OCPs are already fixed by the real cell. The SEI kinetic rate constant
multiplier (`K_SEI_MULT`, analogous to the reference recipe's
`K_SEI_ANCHOR_MULT=0.0017`) is therefore the main knob still being tuned for
CELL064 -- see "TUNING STATUS" below.

Cycling protocol matches CELL064's own real aging schedule as closely as a
single fixed-C-rate DFN experiment can: continuous ~C/3 charge/discharge
aging ("c3d3" in the source data), with a full C/20 discharge RPT every
RPT_INTERVAL cycles (real RPTs land roughly every ~50 EFC -- see
../experimental_data/CELL064_capacity_fade.csv).

TUNING STATUS (update this section as the fit progresses):
    Current defaults (K_SEI_MULT=1.78e-3, SI_LAM_PROP=4.425e-8,
    GR_LAM_PROP=6.66e-7, NEG_POROSITY_FLOOR=0.03) give the best fit so far --
    see cell064_degradation_fit_result.png (and
    cell064_degradation_fit_result_v1_winner_backup.png for the previous,
    LAM-only-fit iteration this replaces). History, briefly:

    1. First fit (RMSE~6.4%SoH) used K_SEI_MULT=2e-5 -- deep enough in the
       "ec reaction limited" flux's reaction-limited regime
       (src/pybamm/models/submodels/interface/sei/sei_growth.py:
       `j_sei=-F*c_0*k_exp/(1+L_over_D*k_exp)`) that essentially all fade
       came from silicon LAM, with LLI ~0. A 25-point sweep confirmed
       D_ec_mult has no leverage in this regime (its effect only shows up
       once L_over_D*k_exp isn't <<1) and found SI_LAM_PROP_mult=0.10 the
       best single-mechanism fit.
    2. Fixing LLI meant raising K_SEI_MULT (into the regime where it isn't
       saturated) while suppressing BOTH LAM proportional terms (not just
       silicon's -- graphite's was still at its full, un-suppressed baseline
       the whole time) to stop pre-knee fade from being too fast. A second
       27-point sweep (K_SEI_MULT x SI_LAM_PROP_MULT x GR_LAM_PROP_MULT,
       D_ec fixed at its default -- confirmed never worse for LLI) found
       excellent pre-knee RMSE (~2.3%SoH) at K_SEI_MULT=6e-4, but NO knee
       formed within 300 cycles at any tested K_SEI_MULT -- LAM suppression
       had also removed one of the two ingredients (LAM's own contribution
       right at the knee) that produced the original run's sharp collapse.
    3. Rather than re-sweep at full (expensive) cost, decomposed the problem
       per instruction: find the RATE-INDEPENDENT ratio of (k_sei, SI_LAM,
       GR_LAM) that makes a knee trigger at the right SoH LEVEL (~90%,
       matching the real RPT148=91.4% just before the knee) using a cheap,
       fast probe (all three rates scaled up together so the knee appears
       within ~50 cycles instead of 300), then apply a uniform TIME-STRETCH
       (divide all three back down together) to slide the knee out to the
       real EFC~152 -- the same principle as
       ../../test_model/test_pore_buffering/pore_buffering_degradation_test_1000cyc.py's
       own TIMESCALE_STRETCH. Worked: 10x-scaled (k=6e-3, si=0.2, gr=0.5)
       gave a knee at 91%->43% between cycle 40-50; scaling back down by
       ~2.96x (to k=1.78e-3, si=0.059, gr=0.148) landed the collapse around
       cycle ~120-150, close to the target.
    4. That collapse was still far too violent (91%->43% within 10 cycles,
       vs. the real ~70-100 EFC-wide post-knee decline) -- diagnosed (by
       user instinct, confirmed empirically) as NEG_POROSITY_FLOOR=0.01
       being inherited unchanged from the reference recipe without ever
       being re-checked for CELL064's own geometry. A fast-probe comparison
       at floor=0.01/0.02/0.03 (same rates otherwise) showed the collapse
       violence is extremely sensitive near the low end -- 0.01 collapses
       to plateau within one batch, 0.02 is steep enough to also trigger a
       genuine SolverError (IDA_BAD_K) at the transition, 0.03 gives a
       clean, gradually-decelerating decline. 0.03 was adopted.
    5. Along the way, found the CV hold step ("Hold at 4.195V until C/100",
       inherited unchanged from the same reference recipe -- both
       pore_buffering_degradation_test_1000cyc.py and dma_baseline_run.py
       use it) was truncating batches near the knee ("final time"
       termination -- the 24h default step duration wasn't enough for the
       current to decay to C/100 once porosity crashed). Loosened to C/100
       -> C/50 throughout this script; confirmed via direct comparison this
       alone let a previously-truncating run (batch 1 at 8/25 cycles) reach
       200+ cycles instead.

    Current result (300 cycles, all fixes above): capacity fade tracks the
    real RPT points closely EFC 0-120, a knee around EFC~120 (real ~152,
    still a bit early), then a genuinely gradual, decelerating post-knee
    decline (unlike every earlier attempt's near-vertical collapse) reaching
    48.6% by ~EFC 207 -- extrapolating that trend lands very close to the
    real 50.0% at EFC 248. LLI now reaches ~60% by EFC 207 (real ~49% at
    EFC 248) -- a real, substantial contribution, up from ~2% in fit #1.
    LAM_Si reaches ~19% by EFC 200 (real ~58-80%, crude non-DMA estimate) --
    still undershooting in magnitude, though the qualitative split
    (Si >> Gr) is right. Not yet re-swept jointly with these new floor/CV
    settings -- there is likely still room to nudge the knee later (toward
    EFC~152) and boost LAM_Si's magnitude without re-breaking the now-good
    LLI/pre-knee behaviour.

    6. EFC re-anchoring (EFC_START_OFFSET, see top of file): RPT1 (EFC50.6)
       reads identical to BoL (100% SoH) -- the real record shows no
       resolvable aging before ~EFC50, so all the "knee timing"/EFC
       comparisons above, which used raw EFC, were comparing the model's
       own t=0 against a stretch of real data carrying no information.
       Re-anchoring the experimental axis to RPT1 shifts the knee target
       from raw-152 to ~101.4 in model-EFC terms -- against THAT target,
       the manually-tuned point's knee at model-EFC~120 is not early, it is
       if anything already slightly LATE. Next tuning step should account
       for this (a small reduction in the time-stretch, not an increase) --
       not yet re-tested.

    7. Finer porosity-floor bracket (0.025, 0.0275) between the working 0.03
       and the crash-prone 0.02, same rates/rest of recipe otherwise, 200
       cycles each: 0.025 hit TWO solver errors (IDA_ERR_FAIL then
       IDA_CONV_FAIL) during its collapse and never completed a full
       requested batch again afterward (91.58%->89.09%->86.34%->76.75%->
       67.16%->65.22% by cycle ~162 of 200 requested) -- i.e. 0.025 already
       reintroduces 0.02-like fragility, it is not a safe middle ground.
       0.0275 completed cleanly with NO solver errors at all (192/200
       cycles, 91.58%->89.09%->86.34%->77.98%->64.67%->58.19%) while giving
       a visibly SHARPER knee than 0.03's gradual decline. 0.0275 is the
       best of both: sharper/more knee-like than 0.03, without 0.025/0.02's
       crash risk -- worth adopting in place of 0.03 for the next iteration
       (see cell064_degradation_fit_result_floor0p025.png and
       _floor0p0275.png, run via C64_NEG_POROSITY_FLOOR/C64_OUT_TAG env
       vars, alongside the current cell064_degradation_fit_result.png at
       floor=0.03 for comparison).

    8. Attempt to push the knee-onset SoH up to ~90-91% (matching real
       RPT3=91.4%, the point just before the real knee) BEFORE doing the
       EFC-timing fix from item 6 -- per instruction, so the "shape" is
       right before re-doing the "timescale" step. Used fine-grained batch
       tracking (C64_BATCH_SIZE=C64_RPT_INTERVAL=15, 150-cycle budget) to
       locate the knee onset precisely for each candidate. Result: **every
       knob tried left the knee onset clustered at SoH~83-86%**, not
       90-91%:
         - NEG_POROSITY_FLOOR 0.03/0.032/0.035 (all with original rates):
           identical pre-knee trajectory through cycle 126 (SoH 84.92% for
           all three), knee onset unchanged: only the POST-onset slope
           differed slightly (steeper for lower floor), confirming floor
           value (at least in 0.03-0.035) sets knee SHAPE, not its SoH
           level.
         - SI_LAM_PROP/GR_LAM_PROP at 0.7x and 0.5x of the reference
           (floor=0.0275): nearly identical trajectories to EACH OTHER
           despite a 30% relative difference (e.g. both ~85.6-85.9% at
           cycle 126) -- LAM's own proportional rate contributes little to
           pre-knee fade at these settings; cutting it only delayed the
           whole curve (and the knee) by ~15-20 cycles without changing the
           SoH level at which the knee eventually triggered.
         - SI_CRIT_STRESS at 0.85x and 0.7x of the reference (floor=0.0275):
           this DID pull the whole curve down faster pre-knee (more visible
           divergence than the LAM cuts), and did make the eventual
           collapse steeper/more violent, but the actual onset SoH still
           landed at ~83.7% (0.85x) and ~83.8% (0.7x) -- same band as
           everything else, just reached a little sooner in cycle count.
         - K_SEI_MULT at 0.5x and 1.5x of the reference (floor=0.0275):
           strong, immediate pre-knee divergence (unlike LAM/crit-stress --
           k_sei clearly has real leverage over the pre-knee SoH-vs-cycle
           curve). But the eventual transition still landed close to the
           same band either way: 1.5x hit TWO solver errors (IDA_CONV_FAIL)
           and knee'd at ~86% before breaking down entirely (final readings
           post-crash are numerical garbage, e.g. SoH=2.02%, tripping the
           45%-floor stop on invalid data); 0.5x ran cleanly for far longer
           (250 cycles) with a near-PERFECTLY LINEAR decline the whole way
           (rate 0.98-1.5%/14cyc, barely accelerating) before a real knee
           finally hit at cycle ~224-238, at ~82% -- i.e. LOWER than the
           ~85-86% band, not higher: slower k_sei delays reaching the
           floor, so MORE background fade (mostly LAM) has accumulated by
           the time it's reached. So k_sei moves the onset SoH, but only
           DOWNWARD as it's reduced; raising it hits a solver-crash ceiling
           around the same ~86% the other knobs converge on.
         - SI_CRACK_RATE_MULT/GR_CRACK_RATE_MULT (the actual lever for how
           much fresh cracked area is exposed to grow SEI-on-cracks, as
           distinct from k_sei's per-area rate or LAM's active-material-loss
           rate) at 1.5x and 2x of the reference (floor=0.0275, both crit-
           stresses/LAM/k_sei at baseline): pre-knee curves nearly identical
           to baseline through cycle ~100 (tiny, slowly-growing gap), then a
           MUCH MORE VIOLENT collapse once triggered -- 1.5x dropped
           86.39%->80.95% in one batch then solver-crashed; 2x dropped
           86.06%->73.83% in one batch (worse) then also crashed. Onset SoH
           unchanged (~86% again) in both cases -- crack rate governs
           collapse VIOLENCE, same role as lowering the porosity floor, not
           knee-SoH level.

       Also tested, per instruction, whether JOINT tuning (rather than one-
       knob-at-a-time) could do better -- the mechanistic reasoning being
       that SEI-on-cracks growth should depend on CRACK RATE (how much
       fresh negative-electrode area gets exposed), not on the LAM
       proportional rate (governs active-material LOSS, a different
       quantity) or on k_sei alone (sets the reaction rate per unit area,
       not how much area exists). Working theory for why every single knob
       above converged on ~85%: the model's SEI growth has a "base" term
       (present everywhere, using k_sei, responsible for the near-perfectly
       LINEAR pre-knee decline seen in literally every probe above) and an
       "on-crack" term (only on freshly-exposed cracked area). If the BASE
       term dominates how fast the floor is approached, boosting crack rate
       alone just adds a violent kick once the base process is already
       nearly there (matches: crack rate changed violence, not onset); and
       scaling k_sei scales both the base fade rate and base porosity-
       consumption rate together, preserving the same SoH-at-floor ratio,
       just compressed/stretched in time (matches: crash ceiling and
       downward-only shift). Prediction: LOWERING k_sei (to shrink the
       base term's share) while RAISING crack rate (to make the floor-
       trigger more crack-dominated) should be able to push the onset UP.
       Tested (floor=0.0275, fine-grained tracking): k_sei=0.85x+crack=1.5x,
       and k_sei=0.7x+crack=2x. Result: **no improvement** -- both still
       knee'd at ~86% (85-90% pre-knee zone drifted through without
       accelerating there, same as every single-knob attempt), then
       solver-crashed (IDA_CONV_FAIL) shortly after, same failure mode as
       the individual crack-rate/high-k_sei probes. Also tested
       GR_CRIT_STRESS x0.7 alone (the one single knob not yet tried,
       graphite's own crack-onset threshold): reproduced baseline almost
       exactly at every batch (e.g. 84.59% vs baseline's 84.92% at cycle
       126) -- negligible individual leverage, as expected since silicon
       dominates the swelling/cracking budget in this composition.

       Conclusion after SIX single knobs (floor 0.025-0.035, LAM
       proportional rate 0.5x-1x, SI_CRIT_STRESS 0.7x-0.85x, K_SEI_MULT
       0.5x-1.5x, crack rate 1x-2x, GR_CRIT_STRESS 0.7x) AND two joint
       (k_sei+crack-rate) combinations: **the ~83-86% knee-onset SoH is
       extremely robust** -- nothing tried moves it up; several things
       (K_SEI_MULT down, and by inference the joint combos) can only push
       it DOWN. This is strong evidence the ~85% band reflects something
       closer to structural in this recipe/geometry (most likely CELL064's
       own real BoL negative-electrode porosity budget in
       cell064_parameters.py, or an emergent property of the "ec reaction
       limited"+"SEI porosity change"+stress-driven-LAM combination itself)
       rather than a rate constant waiting to be found. Untried at this
       point: a genuine multi-dimensional grid sweep over (k_sei, crack
       rate, floor) jointly; or a structural change to the BoL negative
       porosity itself -- turned out to be exactly the missing piece, see
       item 9 below (the "accept ~85%" fallback recommended at the end of
       item 8 was superseded before being acted on).

    9. BoL negative electrode porosity (NEG_POROSITY_BOL, new env-var
       override, default None = cell064_parameters.py's own 0.35
       unchanged): the breakthrough. Realised that every floor value tested
       in item 7/8 (0.025-0.035) is a tiny perturbation against a much
       larger quantity -- cell064_parameters.py's BoL "Negative electrode
       porosity" is 0.35 [SOLVED GEOMETRY], so the "headroom" porosity must
       fall through before hitting even the highest floor tested (0.035) is
       still ~0.315 -- a 0.01 change in floor is a <3% change in that
       headroom, which is almost certainly why floor showed zero effect on
       onset SoH. The real lever is the headroom itself, i.e. the BoL value,
       not the floor. Confirmed against the source project itself
       (../../../../Si_Gr_Expansion_Precursor/Pouch_Data/model_cell064/
       pybamm_parameters/CELL064_BoL_parameter_set.py): its own comment
       reads `"Negative electrode porosity": 0.350000, # [SOLVED GEOMETRY]
       = 1 - (Gr + Si active fraction)` -- i.e. porosity was DEFINED as the
       complement of the active-material fractions, with NO binder/
       conductive-additive volume subtracted out at all. Real composite
       electrodes do have a binder+carbon fraction (typically several %),
       so 0.35 is arguably an overestimate of the true void fraction, not
       just a free tuning knob -- reducing it is a legitimate correction,
       not a fudge (also, per the same source file, its own reference
       recipe already specifies "Negative electrode porosity floor": 0.08,
       "[LITERATURE - TO FIT]" -- never carried into cell064_parameters.py
       since it's a degradation-recipe knob, not a BoL one; noted here for
       completeness but NOT used in the result below, which kept floor at
       the already-validated 0.0275).

       Verified this doesn't corrupt the BoL fit before trusting it (see
       ../_porosity_bol_check.py): a standalone C/20 charge/discharge check
       comparing porosity=0.35 (baseline) against 0.25 and 0.15, with BOTH
       active-material volume fractions left untouched (capacity/OCV/dV-dQ
       depend on those, not on porosity -- porosity only enters the
       electrolyte-transport equations, negligible at C/20's very low rate).
       Result: capacity identical to 4 decimal places (2.4336 Ah) at every
       porosity tested; max voltage deviation over the full discharge curve
       was 1.02 mV at porosity=0.25 and 4.46 mV at porosity=0.15 -- both far
       below the original ~0.7%-capacity-error BoL fit's own tolerance.
       Porosity is safe to move independently.

       Fine-grained probes (floor fixed at 0.0275, all other rates at
       reference) at BoL porosity 0.30, 0.25, 0.20 (headroom to floor:
       0.2725, 0.2225, 0.1725 respectively, vs. baseline's 0.3225) all
       PUSHED THE KNEE-ONSET SOH UP, cleanly ordered by porosity cut size --
       the missing dose-response nothing else produced:
           porosity 0.35 (baseline): onset ~85-86%
           porosity 0.30:            onset ~87%
           porosity 0.25:            onset ~91-92%  <- squarely in the real
                                      target band (RPT3=91.4%)
           porosity 0.20:            onset ~91%      (same band, more violent)
       All four still hit the same post-onset IDA_CONV_FAIL solver crash
       partway down the collapse (a separate, already-understood numerical
       issue -- see item 4/7's porosity-floor discussion), but unlike every
       single-knob/joint attempt in items 7-8, porosity 0.25 and 0.30 both
       recovered and ran to completion (reached MAX_TOTAL_CYCLES cleanly),
       landing at a physically sensible final SoH (53.44% and 54.50%
       respectively) close to the real ~50.0% target at EFC248 -- not
       post-crash numerical garbage like every other "successful" knee in
       items 7-8.

       Recommended next step: pick a BoL porosity in the 0.22-0.25 range
       (0.25 is the best single point tried) as the new default, then redo
       the shape/timescale-decomposition + EFC-re-anchoring steps (items 3
       and 6) around it -- the knee-onset SoH problem this whole TUNING
       STATUS section has been chasing since item 6 is essentially solved;
       what remains is re-fitting the RATE constants (K_SEI_MULT/LAM props)
       to this new, correctly-shaped baseline so the absolute EFC timing and
       LAM/LLI magnitudes line up too.

    10. v2 confirmation (this script's current defaults): adopted
        NEG_POROSITY_BOL=0.25, NEG_POROSITY_FLOOR=0.0275 (item 7's winner,
        not item 9's 0.03 leftover default -- fixed), and re-derived the
        time-stretch for K_SEI_MULT/SI_LAM_PROP/GR_LAM_PROP: the item-3
        values, un-stretched, knee'd around model-EFC~60 at porosity=0.25
        (see item 9's probe); target (EFC_START_OFFSET-adjusted) is ~101, so
        divided all three by 101/60=1.68 (K_SEI_MULT 1.78e-3->1.06e-3,
        SI_LAM_PROP 4.425e-8->2.634e-8, GR_LAM_PROP 6.66e-7->3.964e-7) --
        the same uniform-scaling principle as item 3's own decomposition.
        Result (260-cycle budget, fine-grained batch tracking): knee onset
        at ~91% SoH (cycle ~95->114, i.e. model-EFC ~95-108) -- matches the
        real RPT3=91.4% AND lands almost exactly on the ~101 EFC target in
        the SAME run, on the first attempt -- both the "shape" (item 9) and
        "timescale" (this item) fixes composed cleanly, as the decomposition
        strategy's whole premise predicted they should. Post-knee: hit the
        same IDA_ERR_FAIL/IDA_CONV_FAIL pattern partway down (cycle ~133),
        but recovered and ran to completion at 219/260 cycles, final
        SoH=54.75% -- close to the real 50.0%-at-EFC248 target. See
        cell064_degradation_fit_result.png for the current default result.

        Also tried, per instruction: tightening the solver tolerance
        (root_tol=atol=rtol, via the new C64_SOLVER_TOL env var) from 1e-06
        to 1e-07 and 1e-08, to check whether it explains some visibly-off
        capacity points near the knee and/or lets the solver get further
        through the post-knee crash. Result: both were ~4-5x slower per
        ageing cycle than 1e-06 (~6.5s vs ~1.2-1.7s), with NO further
        slowdown between 1e-07 and 1e-08 -- the cost jumps as soon as you
        leave 1e-06 at all, it does not scale smoothly with tolerance.
        Reverted to 1e-06 per instruction (the extra cost wasn't justified
        without confirmed benefit); the "weird points" issue remains open
        and is NOT yet attributed to solver tolerance -- worth a closer,
        more targeted look (e.g. diffing the raw discharge-capacity array
        near the knee, rather than a blanket tolerance change) before
        trying this again.

        Not yet done: LAM_Si/LAM_Gr split and LLI magnitude have not been
        re-checked against experimental data at this new operating point
        (item 5's fit was against the pre-porosity-fix recipe); DMA
        cross-check still not wired (see below); the post-knee solver
        crash itself is still unresolved (a recoverable one, but a genuine
        crash-free run would be preferable for a final result).

    11. Pushed the knee later (stretch 1.68x -> 2.3x) plus a small floor
        bump (0.0275 -> 0.028), per instruction, because the LLI curve (item
        10) was overshooting real CELL064_LLI.csv at the point crashes cut
        the data off -- less elapsed time between knee and the fixed real
        checkpoints (EFC~171, ~198) means less accumulated LLI by those
        points, without touching k_sei's magnitude. Confirmed: at
        stretch=2.0x, model LLI reached ~42.7% by EFC~152 (vs. real's 41.7%
        at EFC171 -- a near-exact match despite being earlier), a big
        improvement over 1.68x's ~53% overshoot at EFC~170. At 2.3x+0.028,
        model LLI hit ~38-40% at EFC~160-163 -- matching an interpolated
        real curve at that same EFC almost exactly.

        Root-caused TWO real bugs uncovered while chasing this:
        (a) A corrupted-data bug, not a "hard" crash: after
        on_experiment_error logs a solver failure, the returned partial
        solution sometimes contains one near-instant, near-zero-capacity
        cycle (visible as a vertical dip to ~0 A.h in every earlier crash-
        affected plot). This was feeding straight into the SoH-vs-
        SOH_TERMINATION_PERCENT check, causing runs to falsely terminate
        early on garbage (e.g. a real 68%-SoH run reporting SoH=4.02% and
        stopping). Fixed with MIN_VALID_CAP_AH (10% of nominal) in
        cycle_ageing_leg -- cycles below it are now treated as missing data,
        not real data, in both the termination check and the plots (the
        vertical-dip artifact is also gone from every plot since).
        (b) A genuinely separate, unrecoverable exception-raising failure
        (distinct from the graceful, on_experiment_error-logged kind) was
        the TRUE remaining blocker after fixing (a) -- e.g. run halted
        outright with "solver_failure_batch_13" even though the physics
        were still evolving normally. Added a tolerance-retry: on this
        exception, the SAME batch is retried with the solver tolerance
        scaled by [10x, 0.1x, 100x, 0.01x] (a fresh IDAKLUSolver each time)
        before giving up, per instruction that either direction might help
        an otherwise-unpredictable stuck Newton step -- confirmed working
        via repeated "recovered on retry 1 (tol=1e-05)" log lines, i.e. the
        10x-looser retry succeeded every time it was needed in this run.

        With both fixes, three 450-cycle runs (main at the item-10 rates,
        and two with SI_LAM_PROP/GR_LAM_PROP scaled +30%/+60% on top, per
        instruction to test whether boosting LAM -- expected to hurt the
        knee fit, per instinct -- might help the LLI/LAM balance) all
        pushed FAR past every previous stopping point (172-180 cycles) to
        215-226 real cycles (EFC~184-185), reaching MUCH closer to the real
        endpoint (EFC~198) than any prior attempt:
            main (no LAM change):  EFC~185, LLI~48%, SoH~60%
            LAM proportional +30%: EFC~183, SoH~59%  (solver_failure_batch_16,
                                    genuinely exhausted all 4 retries this
                                    time -- a true stop, not a bug)
            LAM proportional +60%: EFC~185, LLI~49% (real EFC198=49.0% --
                                    an almost exact match), SoH~58% (real
                                    EFC171=57.1% -- also very close),
                                    reached MAX_TOTAL_CYCLES cleanly (no
                                    crash at all this time)
        The knee-onset SoH (~90%, from item 9) was NOT degraded by the LAM
        increase -- all three still transition at the same ~90% point per
        batch 7 of each run; LAM's effect is confined to how much EXTRA
        fade the transition batch itself carries, not the trigger SoH. So
        the instinct that "increasing LAM would hurt the knee fit" turned
        out not to bite here -- the +60% LAM variant is adopted as the new
        default (SI_LAM_PROP=3.078e-8, GR_LAM_PROP=4.634e-7, K_SEI_MULT=
        7.74e-4, NEG_POROSITY_FLOOR=0.028) since it gives the best LLI/SoH
        match at the furthest reach of any run so far. MAX_TOTAL_CYCLES
        raised to 450 as the new default (retried batches often make
        partial/zero progress, so more nominal budget is needed to reach
        the same real EFC than before the retry mechanism existed).
        LAM_Si is still dramatically undershooting real's crude estimate
        (~12% model vs. ~58-80% real at EFC170-198) -- unaddressed, and
        likely the next real target now that LLI and the knee are both in
        good shape.

    12. Noticed (per instruction) that even with item 11's fix, the model's
        C/20 RPT capacity post-knee sat noticeably ABOVE the real C/20 RPT
        points at matching EFC -- e.g. the model's own simulated RPT
        capacity was ~1.79 Ah at EFC~154, well above real RPT4's 1.40 Ah at
        EFC~171, and the model's continuous aging curve didn't reach that
        1.40 Ah value until EFC~185 -- i.e. not really a magnitude problem,
        a ~15-20 EFC LAG in how fast the post-knee collapse deepens relative
        to real, even though the knee onset and eventual final SoH both
        already matched well. Hypothesis: NEG_POROSITY_FLOOR sets the
        asymptotic severity of the collapse (via how much further porosity
        can fall before the Bruggeman-relation resistance increase
        saturates) -- a lower floor should deepen the SAME collapse window
        without needing to touch onset timing or rate constants. Retested
        floor=0.025 (previously flagged as crash-prone, but that was at the
        OLD un-stretched, non-LAM-boosted rates -- worth a fresh look now
        that the tolerance-retry mechanism (item 11) exists to absorb any
        extra instability) on top of item 11's winning config (K_SEI_MULT=
        7.74e-4, SI/GR_LAM_PROP +60%, MAX_TOTAL_CYCLES=450). Confirmed the
        hypothesis directly: the SAME transition batch that dropped
        90.40%->82.31% at floor=0.028 dropped 89.85%->77.53% at floor=0.025
        -- a visibly steeper collapse in the identical window. Final result:
        the model's SoH curve now sits almost exactly ON the real EFC171
        point (57.06%) instead of passing well above it, and the final
        reachable point (EFC~183, SoH~55.6%) is closing in on the real
        EFC198 target (49.99%) more directly. Did stall for good this time
        at 226/460 cycles (SoH=55.59%, unchanged for the last ~9 batches --
        the tolerance retries stopped helping, unlike item 11's cases,
        though MAX_TOTAL_CYCLES was reached before this became a real
        problem) -- a slightly more fragile floor than 0.028, as expected,
        but the fit quality gain was worth it. Superseded by item 13 below
        (NEG_POROSITY_FLOOR=0.022 is now the default).

    13. Root-caused (per instruction) and largely fixed the post-knee solver
        crash itself, rather than continuing to route around it. Two
        candidate mechanisms were checked:

        (a) SEI-exponent overflow in the "ec reaction limited" flux. Found a
        REAL bug: sei_growth.py already has a smooth (tanh) cap on the SEI
        reaction exponent, added specifically to fix these exact
        IDA_ERR_FAIL/IDA_CONV_FAIL failures (see its own comment) -- but
        that cap is only wired into the "reaction limited" branch, NOT the
        "ec reaction limited" branch this recipe actually uses. This
        script's own "SEI reaction exponent cap" parameter (EXPONENT_MAX_SEI)
        has therefore been a complete NO-OP the whole time. Root cause:
        `k_exp = k_sei * exp(-alpha_SEI*F_RT*eta_SEI)` is uncapped -- j_sei
        itself is mathematically bounded even as k_exp->infinity (it
        saturates at the diffusion-limited value), but k_exp can numerically
        overflow if eta_SEI swings very negative (plausible right at the
        porosity-floor singularity), turning a finite j_sei into an inf/inf
        NaN in floating point. Fixed by applying the same tanh cap to this
        branch too (src/pybamm/models/submodels/interface/sei/sei_growth.py
        -- core library code, affects any recipe using "ec reaction
        limited", not just this one). Chose EXPONENT_MAX_SEI=50 (not the
        other branch's calibrated 10) after checking the actual regime here:
        exponent_max=10 caps k_exp at only ~4% of the diffusion-limited
        threshold (D_ec/L_sei_0) for this recipe's k_sei/D_ec -- i.e. it
        would have suppressed the whole porosity-floor knee mechanism.
        exponent_max=50 leaves ~1e16x headroom above that threshold (never
        interferes with the physics) while still capping the exponent
        argument to exp() far below float64's overflow point (~709).
        Result: genuine but partial improvement -- e.g. at floor=0.022,
        final reachable point moved from 168 to 174 cycles, and the final
        error type changed (IDA_BAD_K -> "infeasible bounds at initial
        conditions"), confirming this was a real contributor but not the
        sole cause.

        (b) eps_solid (active-material volume fraction) integrating through
        zero unbounded -- ../../submodel_stability_audit.md's own top-
        flagged candidate (`eps_solid = pybamm.Variable(...)` in
        loss_active_material.py has no bounds=, unlike e.g. c_s; if it went
        through zero, a = 3*eps_solid/R -> 0 forces j = i/a -> infinity).
        Checked directly (per instruction) rather than assumed: added a
        diagnostic (C64_DIAG_EPS=1) printing min(eps_solid) for both
        negative-electrode phases every batch. Result: RULED OUT --
        eps_solid Si only dropped from BoL 0.058 to ~0.039 (~33% relative)
        over the ENTIRE run, including well past the knee and every crash
        point; Gr barely moved at all (0.592->0.570). Nowhere close to
        zero at any point. This matches item 8's much earlier finding that
        the knee is LLI-dominated, not LAM-dominated, at this operating
        point -- eps_solid was never a real risk here.

        With fix (a) applied, floor=0.022 -- the value most representative
        of the real dose-response we'd mapped out in items 11-12, but
        previously the least stable -- was retested and gave the best
        result of this entire project: the run reached the SOH_TERMINATION_
        PERCENT=45% floor NATURALLY (336 cycles, EFC~232, SoH=44.5%) instead
        of crash-stopping, the first fully natural completion at any
        aggressive floor value. Its SoH curve passes almost exactly through
        BOTH real RPT4 (EFC171, 57.1%) and RPT5 (EFC198, 50.0%) -- the best
        match to the real capacity-fade curve achieved so far. Adopted as
        the new default (NEG_POROSITY_FLOOR=0.022, EXPONENT_MAX_SEI=50,
        all item-11 rates unchanged; C64_DIAG_EPS left off by default).

        Not yet done: LAM_Si is still dramatically undershooting real's
        crude estimate (unaffected by this item's fixes) -- still the next
        open target; DMA cross-check still not wired.

    14. Two further findings while continuing to chase the remaining post-
        knee crash and the RPT-vs-RPT capacity gap (item 12's finding that
        the model's own C/20 RPT points sit above real RPT4/RPT5 even
        accounting for the C/3-vs-C/20 rate difference -- persists at a
        similar ~15-20 point magnitude across every floor/EXPONENT_MAX_SEI
        combination tried, e.g. NEG_POROSITY_FLOOR=0.021 with
        EXPONENT_MAX_SEI in {20, 50, 75, 100} all show it, suggesting it's
        not something these two knobs can fix -- more likely inherent to
        the rate constants themselves, still unaddressed):

        (a) RUN-TO-RUN NON-DETERMINISM, confirmed directly: re-ran the exact
        same config (NEG_POROSITY_FLOOR=0.022, EXPONENT_MAX_SEI=50, no env
        vars differing, verified via `env | grep C64_`) that reached
        MAX_TOTAL_CYCLES naturally at 336 cycles in item 13 -- this time it
        hit a hard solver_failure_batch_13 at only 174 cycles. Same inputs,
        different outcome. This means the "336-cycle" result was not a
        deterministic property of that config, just one possible outcome;
        likely floating-point-order sensitivity in the sparse solver's
        internals right at the marginal, stiff porosity-floor state (see
        (b) below) tips the outcome either way. Practical implication:
        single-run cycle-count/EFC-reach comparisons between configs (as
        used throughout items 11-13) are noisy signals, not clean ones --
        a config reaching further in one run isn't reliably better, only
        probably so.

        (b) Added a second diagnostic (C64_DIAG_POROSITY=1, alongside item
        13's eps_solid one) printing min("Negative electrode porosity")
        every batch, to directly check (rather than continue assuming)
        that porosity is actually pinned at NEG_POROSITY_FLOOR when the
        remaining crashes happen. CONFIRMED directly and repeatably: at
        floor=0.021, porosity reads exactly 0.021000 (to 6 decimal places)
        for multiple consecutive batches (spanning ~35 cycles) before any
        crash appears, and the eventual IDA_CONV_FAIL/IDA_ERR_FAIL/
        IDA_BAD_K failures all occur while porosity is still pinned at
        exactly this value -- confirming the porosity-floor softplus
        (reaction_driven_porosity.py, k_eps=100 sharpness, see the earlier
        discussion this project had about whether reducing k_eps would
        help) is indeed the site of the singularity, not a red herring.
        Notably, the model runs STABLY while pinned at the floor for a
        variable number of cycles before failing -- consistent with (a)'s
        non-determinism finding: this is a marginal/stiff numerical regime
        where small perturbations eventually trigger failure, not a hard
        immediate blowup the instant porosity is reached. This reframes
        the remaining crash as a genuine, not-yet-fixed numerical
        singularity at the porosity floor itself -- item 13's SEI-exponent
        fix addressed a real, separate contributing bug, but the floor's
        own softplus clamp is a second, distinct source of stiffness that
        would need its own fix (e.g. a genuinely smoother/lower-sharpness
        transition, or reformulating the floor mechanism) to eliminate
        outright, not just work around via floor value/retry-tolerance
        tuning as done so far.

    15. Tested item 14(b)'s "smoother floor transition" hypothesis directly
        (per instruction: "try the softplus clamp to be sure but I do not
        think it may help"). Temporarily made the porosity-floor softplus's
        sharpness (`reaction_driven_porosity.py`'s hardcoded k_eps=100,
        transition width ~1/k_eps) overridable via a
        PYBAMM_POROSITY_FLOOR_KEPS env var and ran three same-config
        (NEG_POROSITY_FLOOR=0.022, EXPONENT_MAX_SEI=50) probes differing
        only in k_eps: 100 (a third baseline sample), 20 (5x wider
        transition), and 5 (20x wider). Result: the hypothesis does NOT
        hold, confirming the user's instinct --

        - k_eps=100 (fresh baseline sample): natural completion at 362
          cycles, no crash -- a THIRD distinct outcome alongside item 14's
          336 (natural) and 174 (crash), reinforcing that this config's
          run-to-run variance spans roughly 170-360+ cycles regardless of
          k_eps.
        - k_eps=20 (moderately smoother): still crashed, same failure
          signature as ever (IDA_ERR_FAIL exhausting all 4 tolerance
          retries, batch 8), merely later (350 cycles) than some k_eps=100
          samples and earlier than others. Smoothing the transition did not
          remove the singularity, only nudged where in the noisy 170-360+
          range this particular run happened to land -- consistent with it
          not being the operative variable at all.
        - k_eps=5 (heavily smoother): by 350 cycles (SoH=79.6%, ageing
          cycles nearly 2.5x the real knee's ~101-152 EFC target) porosity
          had STILL not approached the floor (min=0.103, vs. floor=0.022) --
          this level of smoothing doesn't fix the singularity, it just
          defers ever reaching the region where it occurs past the
          physically relevant range, effectively disabling the whole
          floor/knee mechanism for this fit. Not a usable fix even if it
          "worked" -- it would silently break the knee timing that items
          9-13 spent most of this campaign getting right.

        Conclusion: the porosity-floor softplus sharpness is not the lever
        that controls the crash/non-determinism; item 14's characterization
        (a marginal/stiff regime right at the pinned floor, sensitive to
        floating-point path) stands, but its proposed remedy doesn't pan
        out. Reverted the temporary k_eps override in
        `reaction_driven_porosity.py` (confirmed clean via `git diff`) since
        it added an unused knob to core PyBaMM without fixing anything.
        Given the non-determinism is now confirmed robust to this knob too,
        the practical takeaway is unchanged from item 14: treat any single
        run's cycles-reached as one noisy sample of a 170-360+ cycle range
        for this config, not a reproducible number, and stop tuning against
        it as if it were.

    16. Per instruction: is there a way to give Si more LAM specifically
        AFTER the knee, purely via parameterisation (no new code)? First
        checked whether the EXISTING "stress-driven" mechanism's driving
        signal (Si particle hydrostatic stress) already varies with
        degradation state, which would let m_LAM/stress_critical alone do
        it -- diagnostic showed it does NOT: identical stress_h magnitude
        at batch 1 (pre-knee) and batch 9 (deep post-knee), since it's set
        purely by the fixed C/3 cycling amplitude, not by porosity/
        resistance state. So no coupling is available within "stress-driven"
        alone. Confirmed "stress and reaction-driven" is a pre-existing,
        published PyBaMM option (Ai2019/Reniers2019 -- not new physics,
        already implemented, just not previously enabled here) and checked
        ITS driving signal instead: Si's SEI-on-cracks volumetric current
        density (a_j_sei) jumps ~7x right at the knee (batch3->4: -48.5 ->
        -341.4 A/m3) and stays elevated -- a genuine knee-coupled signal,
        unlike stress. Switched Si (secondary) only to "stress and
        reaction-driven" (Gr/primary untouched) via SI_LAM_OPTION, added
        the new SI_BETA_LAM_SEI parameter (beta_LAM_sei, previously unused
        anywhere in this project -- every PyBaMM built-in parameter set
        ships it at 0.0, so no calibrated value existed to start from).
        Introduced score_rpt_gap() (per instruction, to stop eyeballing
        RPT-vs-RPT match quality off the plot): interpolates the model's
        own C/20 RPT trajectory onto each real RPT's EFC and reports the
        gap in percentage points -- now the standing fit-quality metric,
        printed on every run. Calibrated beta_LAM_sei via a 7-point sweep,
        all starting from the previous (item 13) best-case parameters
        unchanged:
            beta [m3/mol]   RPT4 gap   RPT5 gap   mean|gap|   LAM_Si final
            1e-6            +12.5pp    +10.1pp    8.0pp       32.6%
            2.5e-6          +11.0pp     +9.1pp    7.2pp       33.3%
            5e-6             +6.9pp     +7.4pp    5.4pp       33.4%
            7.5e-6           +4.6pp     +6.0pp    4.2pp       38.4%
            9e-6             +4.5pp     +5.1pp    3.9pp       43.1%
            10e-6            +5.4pp     +4.6pp    4.1pp       43.9%
            15e-6            +1.2pp     +2.2pp    2.0pp       54.2%
        (gap = model - real; +ve = model retains more capacity than real at
        that RPT.) Trend is monotonic overall with mild noise around 9-10e-6
        (consistent with this recipe's known run-to-run sensitivity near the
        porosity floor, see items 14-15) -- 15e-6 is unambiguously best,
        still slightly overshooting rather than crossing to undershoot.
        ADOPTED as the new default (SI_LAM_OPTION="stress and
        reaction-driven", SI_BETA_LAM_SEI=15e-6): a purely-parametric,
        zero-new-code change (one options-string entry + one PyBaMM
        parameter, both pre-existing) that improves RPT-vs-RPT match
        (2.0pp mean gap, down from item 12's ~13-20pp) AND LAM_Si magnitude
        (54.2%, now at the edge of the crude 58-80% real estimate rather
        than the v1 ceiling of ~12-33%) simultaneously, without disturbing
        the pre-knee fit (confirmed: pre-knee SoH identical to v1 to within
        0.1-0.3pp across all beta values tested, since a_j_sei is
        negligible pre-knee by construction). Not pushed further toward an
        exact zero gap (a slightly higher beta, ~17-20e-6, might get there)
        since 15e-6 was judged good enough to move on to expansion-data
        fitting (stage 3) per instruction; revisit if a tighter RPT match
        is wanted later.

    17. Stage 3 (expansion fitting) kickoff. Three new experimental-data
        panels added (reversible expansion, expansion scale k, RPT
        discharge voltage) across three new figures (plot_degradation_
        diagnostics, plot_expansion_diagnostics, plot_overall_summary --
        replacing the old single plot_comparison), plus two new standing
        fit-quality metrics (score_expansion_shape, score_k_gap) alongside
        item 12's score_rpt_gap. First idea for scaling the model's
        single-representative-electrode-pair thickness change up to the
        real pouch cell's whole-stack dilatometry magnitude -- "Number of
        electrodes connected in parallel to make a cell" -- was REJECTED
        before ever being tested: traced the parameter first (per
        instruction) and found geometric_parameters.py's
        "A_cc = L_y*L_z*n_electrodes_parallel" makes it the current-
        collector area, driving capacity/every C-rate/every current
        density -- not a thickness-only knob; scaling it ~8x to match the
        expansion magnitude would silently inflate capacity ~8x too,
        invalidating the whole degradation fit (items 1-16). No physical
        layer count is documented anywhere in the source project's data
        either (checked). ADOPTED INSTEAD (per instruction): compare
        NORMALISED reversible expansion (each series divided by its own
        early-life median) rather than absolute magnitude -- sidesteps the
        geometry coupling with no parameter change. Ran a 3x3 grid over
        f0 (transmitted-fraction plateau) x width (pore-buffering
        transition width) -- f0 in {0.7, 0.55, 0.4}, width in
        {0.01, 0.05, 0.12} -- scored via the two new metrics:
            f0     width   shape_RMSE   k_gap
            0.7    0.01    0.151        0.330
            0.7    0.05    0.096        0.226
            0.7    0.12    0.128        0.205
            0.55   0.01    0.299        0.405
            0.55   0.05    0.218        0.367
            0.55   0.12    0.110        0.315
            0.4    0.01    0.806        0.481
            0.4    0.05    0.738        0.440
            0.4    0.12    0.470        0.332
        f0=0.7 (unchanged from the inherited si_gr_expansion default)
        dominates every width tried -- lowering f0 monotonically hurt both
        metrics here, opposite of test_model's own f0-lowering sweep
        (tuned against a different recipe/target, not transferable).
        Within f0=0.7, width=0.05 and width=0.12 are close (0.05 best on
        shape by ~25% relative margin, 0.12 best on k by ~9% relative
        margin) -- ADOPTED width=0.05 as the new default since reversible
        expansion shape was the primary named target. Confirmed via the
        resulting plot: pre-knee normalised expansion now tracks the real
        near-flat trend closely (vs. the old width=0.01 default's visible
        pre-knee dip), and the knee-peak + post-knee decline both track
        well. k now matches closely pre-/near-knee (model 0.705/0.72/0.945
        vs. real 0.705/0.756/0.949 at the first three RPTs) before hitting
        a STRUCTURAL ceiling post-knee: the "physical" pore-buffering
        formulation cannot produce k>1 (dv_thickness = dv_solid -
        dv_buffered <= dv_solid whenever dv_solid>0, so k = x_average
        (dv_thickness)/x_average(dv_solid) <= 1 by construction), while
        real k rises to ~1.5 post-knee -- no (f0, width) combination can
        reach that regime; per instruction, not chased further for now
        (would need extending the pore-buffering submodel itself, a
        genuine new-physics change, not a parameter fit). Not yet done:
        the electrode-count/absolute-magnitude scaling question remains
        open (normalised comparison sidesteps it, doesn't answer it);
        no further (f0, width) refinement beyond this 3x3 grid (e.g. an
        intermediate width between 0.05/0.12) has been tried.

    18. NEXT STEPS queued (per instruction -- the current fit is already
        judged good enough for publication; these are documented for a
        future session, NOT tried yet, no runs launched for this item):
        push the capacity-fade knee a bit EARLIER and make it a bit
        SMOOTHER (currently a fairly sharp bend, see the SoH panel), via
        some combination of:
          - the degradation-rate constants (K_SEI_MULT, SI_LAM_PROP,
            GR_LAM_PROP, SI/GR_CRIT_STRESS, SI/GR_LAM_EXP -- the item-9/11
            "timescale stretch" pattern shifts EFC timing; shape-only knobs
            like SI_LAM_EXP were never swept, see item 16's stress-diagnostic
            finding that particle stress itself doesn't vary with SoH, so
            any smoothing here would have to come from the SEI/reaction
            side, not stress);
          - NEG_POROSITY_FLOOR (items 7,9,11-13: controls how deep/sharp
            the post-knee collapse is via the Bruggeman relation -- a
            higher floor could soften the bend, at the cost of the
            post-knee fade depth items 12-13 tuned);
          - EXPONENT_MAX_SEI (item 13: caps the ec-reaction-limited SEI
            exponent -- items 14-15 found non-monotonic, noisy sensitivity
            to this value near the floor-pinning regime, so any retuning
            here should re-check items 14-15's non-determinism finding
            still holds, not just the knee shape).
        Earlier timing + a smoother bend are somewhat in tension (a sharper
        transition is partly what let the knee land at the right SoH AND
        EFC simultaneously in items 9-10) -- expect this to need joint
        tuning across at least 2 of the above, not a single-knob fix.

    19. Phased joint sweep -- sweeps/cell064_joint_sweep_phase{1,2,3}.py
        (shared machinery in sweeps/cell064_sweep_lib.py; results per phase
        in sweeps/phase{1,2,3}/; heatmap via sweeps/cell064_sweep_heatmap.py
        <phaseN>). Per instruction: one campaign over BoL porosity, porosity
        floor, EXPONENT_MAX_SEI, K_SEI_MULT, a new joint LAM-proportional-
        rate multiplier C64_LAM_PROP_MULT, and SI_BETA_LAM_SEI, to improve
        the C/20 RPT capacity-fade fit AND the reversible-expansion fit AND
        the knee position at once. C64_NO_PLOT=1 (headless, scores only),
        C64_MAX_CYCLES lowered 450->250 (nothing to fit past EFC~200; the
        long tail-runs OOM'd a 3-parallel sweep), and score_knee +
        C64_TIMESCALE_MULT added this item (see their definitions).

        PHASE 1 -- 2^(6-1) res-VI fractional factorial (32 runs + baseline),
        composite = 0.45 rpt + 0.35 exp_rmse + 0.20 k_gap. Findings:
          * K_SEI_MULT is the DOMINANT knob (main effect -0.27 on composite,
            ~4-6x every other knob) -- higher SEI rate is much better,
            mostly by crushing the k-gap.
          * BoL porosity 0.28 is TOXIC on its own (RPT gap 21-26pp in
            run17/18/21/22/25/26/29/30) -- survivable ONLY when paired with
            high K_SEI_MULT (strong k_sei x porosity interaction).
          * porosity floor low (0.018) and exp_max high (100) mildly help;
            LAM_PROP_MULT high mildly helps; SI_BETA_LAM_SEI ~negligible.
          * KNEE POSITION (not in the phase-1 composite -- caught only by
            post-hoc knee extraction): high K_SEI_MULT pulls the knee
            earlier, often TOO early (run03/04/08/11/12/15/16 knee onset
            ~EFC 51-61 vs the ~101 anchor -> wrecked exp_rmse). Low
            K_SEI_MULT -> knee too late (run06 ~144, baseline ~112).
          * NO combo dominates the baseline on all of rpt+exp+k. Top by
            composite: run06 (rpt 1.2 / exp 0.125 / k 0.093, but knee ~144
            -- good RPT from a soft/late collapse, not a right knee),
            run31 (rpt 3.4 / exp 0.069 [best in sweep] / k 0.211, knee ~98
            -- right on the anchor), run19 (rpt 2.9 / exp 0.129 / k 0.076,
            knee ~100). run31 is the best "earlier knee + expansion" corner:
            poro 0.28 / floor 0.028 / exp_max 100 / k_sei 1.24e-3 (HIGH) /
            lam 0.6 / beta 7e-6.

        PHASE 2 -- focused 3-level grid around run31: K_SEI_MULT
        {1.0e-3, 1.24e-3, 1.6e-3} x BoL porosity {0.25, 0.28, 0.30} x
        porosity floor {0.024, 0.028, 0.032}, exp_max/lam/beta fixed at
        run31's values. Composite now 0.30 rpt + 0.30 exp_rmse + 0.20 k_gap
        + 0.20 knee_efc_err (score_knee, keyed off the SoH curve's
        steepest-descent EFC vs the ~101 anchor). Findings:
          * porosity FLOOR is now the dominant knob (main effect +0.108,
            ~2x k_sei) -- LOW floor (0.024) is much better than high
            (0.032). All top runs use 0.024 (phase 2's low end).
          * WINNER run01: k_sei 1.0e-3 / BoL 0.25 / floor 0.024 --
            rpt 2.7 / exp_rmse 0.052 (best of any sweep run) / k_gap 0.045
            (5x better than the default's 0.226) / knee EFC ~108 (err 6.3
            vs the default's ~32). Only regression vs the ORIGINAL default
            is rpt 2.7 vs 1.6 -- within noise, and dwarfed by the knee /
            expansion / k gains. run01 keeps BoL porosity at 0.25 (the
            original value) -- phase 1's "need 0.28" was an artifact of its
            other knob settings.
          * run25 (k_sei 1.6e-3 / BoL 0.30 / floor 0.024) ran clean to
            completion (run01 hit a solver_failure mid-run -- lower floors
            are more crash-prone, items 14-15) with an even better knee
            (err 3.5) at rpt 4.4 / exp 0.055 / k 0.023.
          * NOTE the model's SoH at its own steepest-descent point is
            ~78-86% across the sweep, below the real ~91% anchor -- i.e.
            the model still loses a bit too much capacity BEFORE the sharp
            collapse (pre-knee fade slightly too fast). Not separately
            targeted yet.

        PHASE 3 -- fine grid in phase 2's low-floor region: K_SEI_MULT
        {0.9e-3, 1.1e-3, 1.35e-3, 1.6e-3} x BoL porosity {0.25, 0.28, 0.30}
        x porosity floor {0.020, 0.024}, + run01/run25 anchors. Findings:
        floor 0.024 is the sweet spot (going to 0.020 is worse -- main
        effect flipped sign vs phase 2, so 0.023-0.024 is bracketed);
        k_sei effect now weak; BoL 0.25 still best. BUT: on inspecting the
        run logs, EVERY high-k_sei (>=1.0e-3) top combo CRASHES at batch ~5
        (~EFC 130-200) with retries eating the cycle budget -- so their
        model RPTs only span ~EFC 50-88, and score_rpt_gap / score_k_gap
        end up scored against a SINGLE pre-knee RPT (EFC 49): near-vacuous.
        The impressive "rpt 2.7 / k 0.045" on run01/run08 do NOT test the
        post-knee fit; the SoH figure shows those combos actually collapse
        ~40 EFC too fast post-knee. So the phase 2/3 composite ranking was
        misled for those two terms.

        PHASE 4 (tie-break) -- k_sei {8e-4, 9e-4, 1.0e-3} x floor {0.023,
        0.024, 0.025}, BoL 0.25 / exp_max 100 / lam 0.6 / beta 7e-6 fixed,
        C64_MAX_CYCLES RAISED to 350 so RPT4/RPT5 land inside the model's
        RPT range and rpt/k become meaningful again. Result: k_sei sets the
        knee EFC in ~discrete steps (7.74e-4 -> ~133, 0.9e-3 -> ~118,
        1.0e-3 -> ~108; floor barely moves it). The knee-to-~101 +
        best-expansion-RMSE + numerical-stability combination is NOT
        achievable with these six knobs: floor 0.024 gives run08's much
        better expansion (RMSE 0.052) but is unstable past EFC~200; floor
        0.023 is stable but its sharper post-knee collapse pushes the
        expansion RMSE back to ~0.10.

        RECOMMENDATION / ADOPTED (item 19 final): run07 --
        C64_K_SEI_MULT=1.0e-3, C64_NEG_POROSITY_FLOOR=0.023,
        C64_EXPONENT_MAX_SEI=100, C64_LAM_PROP_MULT=0.6,
        C64_SI_BETA_LAM_SEI=7e-6, C64_NEG_POROSITY_BOL=0.25. Stable (natural
        SoH-floor completion), 2 validated in-range RPTs including RPT4
        (EFC171: model 56.3% vs real 57.1%, gap -0.8pp), RPT-vs-RPT mean
        |gap| 1.7pp (= the old default's 1.6), expansion RMSE 0.102 (old
        default 0.096 -- negligibly worse), k_gap 0.20 (old 0.226), and the
        KNEE moves from EFC ~133 to ~108 (err vs the ~101 anchor: +32 ->
        +6.3). Net: ~25 EFC earlier, better-aligned knee for no real cost.
        These five defaults are now updated in this file (see each knob's
        comment). NOT resolved: the last ~6 EFC of knee lateness and the
        expansion-RMSE-vs-stability tradeoff -- these need the deferred
        porosity-floor softplus mechanism change (item 15), not more knob
        tuning. Alternative if expansion fidelity is preferred over a
        validated end-of-life fade: floor 0.024 (run08), accepting the
        ~EFC200 solver failure. All sweep artefacts: sweeps/phase{1,2,3,4}/.
        FINAL, standing-reference writeup (parameter set, fit quality, the
        k<=1 post-knee expansion finding, model-vs-real RPT discharge
        overlay): ../CELL064_final_parameterisation.md.

Cross-check via DMA (per instruction: prefer DMA-on-model-output over the
theoretical/built-in LAM/LLI variables when judging a candidate fit,
consistent with how ../../test_model/degradation_test_matrix/test_DMA/ was
itself validated against Method A). NOT yet wired -- see
../../test_model/degradation_test_matrix/test_DMA/dma_ocp_fit_final.py's
pOCV reconstruction, which needs CELL064's own three electrode OCP functions
(cell064_parameters.graphite_ocp_cell064 /
silicon_ocp_average_cell064 / nmc622_ocp_cell064) substituted for its
current Mark2016/Enertech/Chen2020 functions before it can fit CELL064's own
simulated RPT discharge curves. Do this once a capacity-fade-only fit looks
reasonable -- fitting DMA to a not-yet-converged capacity curve wastes the
effort twice over.

    20. Porosity-gated Si-LAM ("v7 isolation"): NEW physical mechanism
        added to loss_active_material.py (not just a knob on the existing
        item-1-19 recipe) -- a second, independent reaction-driven Si-LAM
        term that stays small pre-knee and grows sharply as the negative
        electrode's porosity approaches its floor, on top of the item-16
        beta_LAM_sei term. Introduced because no combination of the
        existing LAM/LLI knobs could reproduce the real DMA split's
        qualitative shape (small pre-knee Si-LAM, sharp post-knee rise
        tied to the SAME event that produces the SoH knee) without
        conflicting parameter-tuning tradeoffs.

        v7a (algebraic gate) -- SI_BETA_LAM_ISO, SI_LAM_ISO_EXPONENT (env
        C64_SI_BETA_LAM_ISO / C64_SI_LAM_ISO_EXPONENT; new PyBaMM params
        beta_LAM_iso [m3.mol-1], eta_LAM_iso; new LAM option string
        "stress and reaction-driven and porosity isolation"). Term:
        beta_LAM_iso * a_j_sei/F * isolation_gate * remaining_frac, where
        isolation_gate = (1-headroom_frac)**eta_LAM_iso is an INSTANTANEOUS
        function of current porosity (headroom_frac: 1 at as-set-up
        porosity, 0 at NEG_POROSITY_FLOOR), remaining_frac =
        min(eps_solid/eps_solid_init, 1) self-limits the term as Si is
        consumed (guards the same current-focusing singularity,
        a=3*eps_solid/R -> 0 => j=i/a -> infinity, as an earlier, fully
        reverted isolation attempt this session hit). Driving current is
        a_j_sei (Si's SEI volumetric current density), not the main
        intercalation current -- deliberately: this model's own porosity
        submodel (reaction_driven_porosity.py) already makes porosity
        decline as a multiple of SEI/plating thickness growth, so a_j_sei
        is the mechanistic cause of the pore collapse this gate responds
        to, not an arbitrary proxy (a "current-focusing" argument for
        swapping in the main current was tried and reverted -- see below).
        Achieved, for the first time this project: small pre-knee Si-LAM,
        a genuine sharp two-regime SoH knee, and a real expansion "hump",
        all as natural consequences of tying Si-LAM to the same porosity
        variable driving the knee (best point: NEG_POROSITY_BOL=0.24,
        SI_BETA_LAM_ISO=8e-4, SI_LAM_ISO_EXPONENT=3.0, K_SEI_MULT default
        1.0e-3 -- knee_err -0.1 to -1.3, LAM_Si ~80-82% at end of life,
        rpt_gap ~2-3pp). Problem: the post-knee transition was far too
        SHARP vs. the real cell (model crashes 45%-SoH-floor within ~13-15
        EFC of the knee, vs. real RPT4/RPT5 ~26 EFC apart at much lower
        SoH), and this "spread" could not be widened by any single- or
        joint-knob retune tried: eta_LAM_iso softening (spreads the ramp
        but also leaks Si-LAM earlier pre-knee, and the exponent's own
        response is non-monotonic -- exp=3.0 beat both 2.0 and 4.0);
        EXPONENT_MAX_SEI alone (non-monotonic/noisy, 13.0->12.9->10.8->17.8
        EFC spread at 70/50/30/15); jointly raising NEG_POROSITY_BOL (more
        headroom) WITH K_SEI_MULT (to re-center the now-later knee) --
        spread stayed pinned at ~13-15 EFC regardless, because
        isolation_gate is a pure function of instantaneous porosity with
        no timescale of its own, and porosity's own collapse rate through
        the critical near-floor region turned out to be dominated (~9x) by
        GRAPHITE's SEI current transiently spiking ~2.9x its BoL value
        right at the knee (untouched by anything Si-specific); a hard rate
        CAP on the isolation term (max_LAM_iso_rate [s-1], floors the term
        before -- not after -- the same singularity guard) -- negligible
        effect while loose, broke the solver (stress_h_surf-style garbage
        values) once tight enough to matter, fully reverted.

        Direct mechanism attribution (loss_active_material.py now exposes
        each LAM branch's own rate as a named variable -- "...LAM rate
        from stress/reaction (beta_LAM_sei)/porosity isolation [s-1]" --
        rather than requiring ablation runs or hand-reconstructed
        formulas): confirmed the porosity-isolation term is the SOLE
        driver of the sharp transition (3-4 orders of magnitude above the
        other two branches at every point that matters; summing all three
        against the model's own total deps_solid_dt leaves zero residual).
        The mechanical stress-driven branch is a red herring here --
        essentially zero throughout, despite a genuine, separately-fixed
        numerical singularity in it (particle-mechanics stress fields
        overflowing to ~1e35-1e38 Pa as eps_solid->0 once Si is nearly
        consumed; fixed by damping the stress INPUT by remaining_frac
        before the **m_LAM power law, not the term's output after it, so
        the power law never amplifies an already-diverging value -- kept
        as a standing numerical-hygiene fix, but it does not affect fit
        quality since this branch was never the active driver).

        v7b (relaxation-lag gate, current version) -- the fix once v7a's
        bottleneck was correctly identified as isolation_gate having no
        timescale of its own: isolation_gate is now a genuine ODE STATE
        (get_fundamental_variables' "...porosity-isolation gate state",
        one new pybamm.Variable per negative-electrode phase using
        "porosity isolation") that RELAXES toward the same algebraic
        target (isolation_gate_target, still (1-headroom_frac)**eta_LAM_iso)
        with its own time constant tau_LAM_iso [s] (env
        C64_SI_TAU_LAM_ISO; new param "...porosity-isolation gate time
        constant [s]"; d(isolation_state)/dt = (isolation_gate_target -
        isolation_state)/tau_LAM_iso, initial condition 0, matching the
        target's own BoL value). Default tau=1.0s is effectively "no lag"
        (reproduces v7a's instantaneous-gate behaviour exactly -- verified
        bit-for-bit against the v7a reference point). This is the first
        lever this project has found that decouples transition WIDTH from
        onset TIMING: bracketing at the item-19-style -0.1-knee-err
        reference point (bol=0.24, K_SEI_MULT default, exponent=3.0,
        beta_LAM_iso=8e-4) found tau<=3e5s (~28 EFC-equivalent) has
        negligible effect (spread stays ~14-15 EFC, same ceiling as every
        v7a knob), but tau=1e6s (~93 EFC-equivalent) is the first value to
        break it: spread 14.8->26.5 EFC, knee_err still -1.3 (onset
        timing preserved, as the mechanism predicts), LAM_Si 80.1% (still
        on the ~80% target), rpt_gap improved to 2.2pp (best of the whole
        session). Pushing further is NON-monotonic, not simply "more
        lag = more spread": tau=2e6/3e6s regress (model exhausts its fixed
        250-cycle budget before the now-slower isolation process can
        catch up -- stop reason flips to "reached MAX_TOTAL_CYCLES" at a
        LOWER end-EFC with LAM_Si undershooting to 54-59%), while tau=5e6s
        recovers and reaches the FURTHEST yet (end_efc~139, closest any
        config has come to real RPT4 -- 33 EFC away vs 33-115 in every
        prior config) but still undershoots LAM_Si (66%, vs the 80% seen
        at tau=1e6s) -- i.e. tau trades REACH against MAGNITUDE in a
        discontinuous, floor-crossing-vs-cycle-cap-crossing landscape, not
        a smooth one. IN PROGRESS at time of writing: a joint (tau,
        beta_LAM_iso) search (nudging beta_LAM_iso up from 8e-4 to try to
        recover magnitude at the larger, better-reaching tau values) to
        find a single point that gets both. No adopted/recommended
        defaults yet -- tau=1e6s is the best SINGLE-metric point tried so
        far (magnitude-correct, biggest single jump in spread) but not
        necessarily the final answer once the joint search lands.

    21. Follow-up to item 20's joint search: an 108+27-combo multi-
        parameter sweep (K_SEI_MULT/BOL/beta_LAM_iso, tau held at
        promising values) found the joint search itself was unreliable --
        several "best" results turned out to compare the model's own
        dead-end EFC against real RPT5 across a 90+ EFC gap (same class of
        blind spot score_rpt_gap's in-range filter has), giving misleading
        good LAM_Si/LLI matches that weren't real. Backing out to a clean,
        single-variable tau bracket (beta_LAM_iso/BOL/K_SEI_MULT all held
        at the item-20 anchor) gave the decisive, trustworthy result: knee
        EFC stays flat (err -0.8 to -2.2) across THREE-PLUS orders of
        magnitude of tau (3e4 to 1e9) -- tau is a clean, dedicated post-
        knee lever, fully confirmed. But pure-tau reach peaks at tau=3e7
        (end_efc~156, spread~57) and is NON-monotonic beyond that (cycle-
        budget exhaustion, confirmed by C64_MAX_CYCLES=250 raised to 500:
        a diagnostic run at tau=3e7/max_cycles=500 showed Q_side (LLI,
        side-reactions only, independent of Si-LAM) ALREADY exceeds its
        real target by EFC~156 (~51% model vs ~40% real interpolated) --
        LLI, governed entirely by K_SEI_MULT and untouched by anything in
        this investigation, has become the actual bottleneck crashing SoH,
        while Si-LAM (only ~37% of its 80% target at that point) is if
        anything UNDERSHOOTING, not overshooting. A follow-up 27-combo
        sweep (K_SEI_MULT lowered to {5e-4,7e-4,9e-4} paired with
        compensating NEG_POROSITY_BOL, beta_LAM_iso raised per item 20's
        original ask) confirmed the K_SEI_MULT/BOL PAIRING correctly
        preserves knee timing (knee_err flat at the SAME value across all
        3 beta_LAM_iso values tested at each pair, confirming beta doesn't
        touch knee timing either) and meaningfully improved LLI (dev
        shrinking from +12-20pp to +9.6-15.8pp) and reach (best end_efc
        166.4, spread 72.8 -- both new session bests) -- but also revealed
        a genuine, counter-intuitive trade-off: LOWERING beta_LAM_iso
        improved reach further at fixed K/BOL (a stronger isolation term
        adds to OVERALL capacity loss too, shortening life and shrinking
        the total EFC window available for everything, including Si-LAM
        itself, to develop in) while making LAM_Si's own deviation WORSE
        (-32 to -40pp undershoot even at the best-found beta for each
        pair) -- i.e. beta_LAM_iso could not simultaneously fix reach AND
        magnitude within the v7b (a_j_sei-driven) formulation.

        v7c reformulation (per instruction): root-caused why beta_LAM_iso
        alone can't fix this -- j_iso = beta_LAM_iso * a_j_sei/F *
        isolation_state * remaining_frac is a PRODUCT of the (now
        correctly tau-paced) state and a_j_sei, which is diffusion-limited
        and crashes to ~0 shortly after the knee (electrolyte starvation,
        confirmed via direct measurement much earlier this session) --
        i.e. a race between a slowly-rising state and a rapidly-dying
        current, capping how much widening tau can EVER buy regardless of
        magnitude tuning. Fix: drive the isolation term with L_sei (Si's
        SEI film THICKNESS, a monotonically non-decreasing accumulated
        STOCK -- measured directly this session: BoL ~6.9e-9 m growing
        smoothly to ~1.3e-7 m at the knee, then ACCELERATING further to
        ~3.6e-7 m by EFC~107, never collapsing) instead of a_j_sei (a
        FLOW). New term: j_iso = beta_LAM_iso * L_sei * isolation_state *
        remaining_frac (beta_LAM_iso units changed [m3.mol-1] ->
        [s-1.m-1]; parameter STRING renamed to match in both
        lithium_ion_parameters.py and this file's PARAM_UPDATES, so old
        env-var values need rescaling -- ~1.5 [s-1.m-1] matches the old
        term's peak rate scale, per direct L_sei magnitude measurement).
        Physical argument: this stays tied to the SAME underlying SEI-
        deposition process already driving porosity collapse in this
        model (reaction_driven_porosity.py) -- not a new, disconnected
        mechanism -- just measuring the ACCUMULATED deposit rather than
        its instantaneous deposition rate, which is the physically
        appropriate quantity for "how much clogging material has built up
        to isolate active material" and does not share a_j_sei's
        diffusion-starvation collapse.

        CRITICAL BUG found and fixed during v7c's first test (comparable
        config: BOL=0.24, K_SEI_MULT=1.0e-3, tau=3e7, beta_LAM_iso=1.5
        [s-1.m-1], C64_MAX_CYCLES=500): a_j_sei (v7b's driver) is always
        negative by convention, so it supplied the "this DEPLETES active
        material" sign for free; L_sei (a thickness) is always >=0, so
        that sign has to be made EXPLICIT -- the first implementation
        omitted the leading minus sign, so the isolation term was silently
        ADDING active material back instead of removing it. Symptoms this
        produced (all directly confirmed via a targeted diagnostic
        tracking eps_solid_Si/isolation_state/L_sei/j_iso per-cycle):
        eps_solid_Si growing PAST its own BoL value (LAM_Si reported as
        negative, down to -55%), a false SoH "recovery" bump post-knee,
        knee timing shifted +27.5 EFC later (breaking the previously rock-
        solid knee/tau decoupling), and eventually a genuine solver
        breakdown (repeated IDA_CONV_FAIL, step size collapsing to
        ~1e-59) once some other volume-conservation constraint was
        violated by Si's fraction growing unbounded. A SEPARATE, also-real
        bug was found and fixed in passing while investigating this (see
        remaining_frac's own doc-comment): remaining_frac was only capped
        at 1, never floored at 0, so a solver overshoot past eps_solid=0
        could flip the WHOLE term's sign the same way -- harmless with
        a_j_sei's naturally-decaying magnitude, but a latent risk with any
        stronger/non-decaying driver; fixed for both the isolation and
        stress-driven branches.

        With the sign bug fixed: knee timing is FULLY restored (knee_err
        -2.2, identical to the same config's v7b-era result) and the SoH/
        LAM_Si trajectory is genuinely monotonic (no oscillation, no false
        recovery) -- confirms isolation_state/tau/L_sei are all working
        correctly together. However, beta_LAM_iso=1.5 (calibrated by
        matching v7b's PEAK rate magnitude) reaches only end_efc~123-124
        (LAM_Si~57%) -- SHORTER than v7b's own best pure-tau result
        (end_efc~156 at the same tau=3e7) -- because L_sei ACCELERATES
        post-knee (measured directly: ~1.3e-7 m at the knee -> ~3.6-3.8e-7
        m by EFC~110-130, still climbing) rather than decaying like
        a_j_sei did, so matching a_j_sei's PEAK was too aggressive a
        starting guess. beta_LAM_iso needs its own fresh recalibration
        (likely lower) for this formulation -- peak-magnitude-matching
        against the old driver is not the right calibration method given
        the qualitatively different (growing vs decaying) time-profile.

    23. Redirect-to-LAM (real breakthrough, after a long session of
        LAM/LLI-split whack-a-mole): every attempt to reduce post-knee LLI
        by scaling SEI-kinetics parameters down (K_SEI_MULT either phase,
        EXPONENT_MAX_SEI) and re-centering the knee via a compensating
        porosity-BOL shift reproduced the SAME post-knee LLI growth rate
        relative to knee-relative time, just delayed -- proving LLI's PACE
        was never actually slowing, only starting later. Root cause: every
        prior isolation mechanism (v7a/v7b/v7c) ADDED a second, independent
        lithium-consumption channel (isolation) on TOP of an unchanged SEI
        growth rate, rather than reallocating the existing budget.

        New mechanism: "SEI reaction redirect to LAM" (new model option,
        default "false", fully backward compatible) -- a growing fraction
        of the SAME SEI reaction current (gated by the existing porosity-
        isolation gate state, tau-lagged as before) stops forming new SEI
        film (sei_growth.py's set_rhs: dcdt_sei *= (1 - isolation_gate))
        and is redirected into active-material loss instead
        (loss_active_material.py: j_iso = redirect_LAM_yield *
        isolation_gate * a_j_sei / F * remaining_frac). Conserves total
        reaction current -- no re-tuning of K_SEI_MULT needed to avoid
        double-counting. Bonus: a_j_sei's own diffusion-limited decay is
        itself slowed by the redirect (less current still building L_sei
        means slower-growing diffusion resistance), keeping the isolation
        driving signal alive longer than in the reaction-only (v7b) case.

        Verified directly (plot_redirect_lli_check.py, redirect_final_
        check.png): pure LLI (Q_side, SEI-only) decelerates sharply at the
        knee (shrinking per-EFC increments) instead of continuing to
        accelerate, while corrected LLI (Q_side + LAM-trapped) tracks
        LAM_Si's own S-curve shape post-knee -- exactly the qualitative
        behaviour predicted.

        redirect_LAM_yield tuning found to be highly non-linear: yield=1.0
        and yield=0.15 both drove LAM_Si to ~93-97% (near-total Si
        depletion) by the time the model's own SoH hit its floor --
        confirms yield only controls SPEED toward eventual full depletion,
        not final extent, since porosity stays pinned at its floor
        indefinitely once collapsed (isolation_gate always eventually
        approaches 1 given enough elapsed time/cycles). Only much smaller
        yields (~0.01-0.015) meaningfully changed the outcome. Real data's
        own steep local slope between RPT4 (58.53%) and RPT5 (79.65%, just
        26 EFC apart) means no single yield threads both exactly --
        yield=0.01 favours RPT4 (dev -2.2pp) at RPT5's expense (dev
        -8.75pp); yield=0.015 favours RPT5 (dev -0.45pp) at RPT4's expense
        (dev +8.5pp); yield=0.0125 balances both (+3.5pp / -4.6pp).

        BEST CONFIG FOUND (all env vars, on top of defaults):
        C64_SI_BETA_LAM_ISO=1.0, C64_SI_TAU_LAM_ISO=2e8,
        C64_NEG_POROSITY_BOL=0.20, C64_SI_K_SEI_MULT=6e-4,
        C64_GR_K_SEI_MULT=6e-4, C64_SI_REDIRECT_TO_LAM=1,
        C64_SI_REDIRECT_LAM_YIELD=0.0125, C64_NEG_POROSITY_FLOOR=0.03,
        C64_EXPONENT_MAX_SEI=70, C64_MAX_CYCLES=700, C64_SOH_FLOOR=30.
        Result: knee_err=+1.8, model reaches AND exceeds both real RPT4
        (171.4) and RPT5 (197.4) for the first time all session with
        isolation active. LAM_Si dev -2.7pp/-6.2pp at RPT4/RPT5. Pure LLI
        dev only +1.5pp/+1.8pp at RPT4/RPT5 (essentially exact). Capacity
        gap improved to -7.5pp/-12.0pp (from -15 to -17pp pre-floor-fix)
        but NOT fully closed -- this residual overall-capacity-too-fast
        issue is the one metric not yet fully resolved; NEG_POROSITY_FLOOR
        helped measurably but K_SEI_MULT reduction (tried 5e-4, 5.5e-4)
        did not move it further. Worth revisiting: whether floor itself
        can go a bit higher, or whether the residual gap is now dominated
        by something else entirely (e.g. positive-electrode LAM, currently
        untouched all session).

        FURTHER refined same night: raising the floor again to 0.035 (same
        everything else) improved the capacity gap further to -3.1pp/
        -7.2pp (mean 3.9pp -- best of the whole session) but pulled LAM_Si
        down to -6.2pp/-11.5pp (yield=0.0125 was calibrated against
        floor=0.03's slightly faster pace). NOTE: the automated "knee EFC"
        detector mis-fires on this floor's smoother SoH curve (reports
        wildly wrong values like +91-110 EFC error) -- the SoH-91%-crossing
        proxy is unaffected and lands within 0.1-0.2 EFC of the anchor in
        every floor=0.035 run, confirming the ACTUAL onset timing stayed
        essentially perfect; only the detector's characteristic-shape
        assumption broke, not the underlying fit. Re-bracketing yield at
        floor=0.035 (0.0125 -> 0.014 -> 0.016) traces a genuine capacity-
        vs-LAM_Si trade-off, not a further improvement in the same
        direction: 0.0125 gives the best capacity match (mean 3.9pp) at
        the cost of LAM_Si undershooting more; 0.016 gives the tightest
        LAM_Si match of the entire session (dev only -1.3pp/-3.4pp at
        RPT4/RPT5) at the cost of capacity regressing to mean 6.3pp
        (still far better than any pre-floor-fix value, ~9-17pp). 0.014
        interpolates between them (capacity 5.0pp, LAM_Si -4.2pp/-7.0pp).
        Given LAM_Si magnitude was this investigation's primary, most-
        emphasized target all session, floor=0.035/yield=0.016 is the
        recommended final candidate: knee timing essentially perfect
        (SoH-91%-crossing EFC=100.8), LAM_Si dev -1.3pp/-3.4pp, pure LLI
        dev untested at this exact point (was +1.5/+1.8pp at the
        yield=0.0125 point, expected similar), capacity gap mean 6.3pp
        (a ~2-3x improvement over every pre-floor-fix isolation-active
        config tried this session). The capacity gap is not fully closed
        at ANY point on this trade-off curve -- likely needs either an
        even higher floor (not yet tried past 0.035) or a genuinely
        different lever (positive-electrode LAM was never touched this
        whole session and remains a candidate).

        FOLLOW-UP (next morning): asked directly why pure LLI keeps
        climbing post-knee at all despite the redirect -- two real causes,
        confirmed by tracing the code: (a) graphite's own SEI growth is
        COMPLETELY UNGATED by the redirect (the option only checks for
        "porosity" in that SPECIFIC phase's own LAM option; graphite's is
        just "stress-driven", so it has no isolation_gate to reference at
        all and keeps consuming lithium at its full, normal rate the whole
        time -- "Total capacity lost to side reactions [A.h]" (Q_side) sums
        BOTH phases), and (b) even Si's own gate never reaches exactly 1
        (softplus-smoothed porosity floor + tau-lag mean a vanishingly
        small but nonzero residual current always still feeds SEI growth).
        (a) is the dominant effect given Gr and Si shared the same
        K_SEI_MULT magnitude in every config up to this point.

        Acted on directly: lowered GR_K_SEI_MULT independently (Si's own
        SI_K_SEI_MULT, which drives L_sei/isolation, left at 6e-4
        throughout) to reduce graphite's own contribution, then bumped
        redirect_LAM_yield up to compensate and pull LAM_Si closer to the
        RPT5 target given the extra headroom from lower LLI.

        NEW BEST CONFIG (supersedes the floor=0.035/yield=0.016 point
        above): C64_SI_BETA_LAM_ISO=1.0, C64_SI_TAU_LAM_ISO=2e8,
        C64_NEG_POROSITY_BOL=0.18, C64_SI_K_SEI_MULT=6e-4,
        C64_GR_K_SEI_MULT=4e-4 (lowered from 6e-4),
        C64_SI_REDIRECT_TO_LAM=1, C64_SI_REDIRECT_LAM_YIELD=0.02 (raised
        from 0.016), C64_NEG_POROSITY_FLOOR=0.035, C64_EXPONENT_MAX_SEI=70,
        C64_MAX_CYCLES=700, C64_SOH_FLOOR=30.
        Result: SoH-91%-crossing EFC=100.4 (anchor ~101, essentially
        exact), model reaches exactly both real RPT4 (171.4) and RPT5
        (197.4). LAM_Si @ RPT5 = 79.5% (real 79.65%, dev only -0.15pp --
        the best match to the primary target the whole session). LAM_Si @
        RPT4 = 65.7% (dev +7.2pp, a bit high). Pure LLI @ RPT4/RPT5 = 35.9%
        / 46.0% (real 39.94%/47.03%, dev -4.0pp/-1.0pp -- LLI now
        genuinely LOWER than target at both points, as intended). Capacity
        gap mean 7.7pp (a bit worse than the 5.1pp seen at yield=0.016 on
        this same GR=4e-4/BOL=0.18 base -- confirms the same capacity-vs-
        LAM_Si trade-off pattern as the floor bracket, now along the yield
        axis too). Plots: cell064_*_result_redirect_bestv2.png.

        Not yet explored: whether GR_K_SEI_MULT can go lower still (further
        LLI reduction) with a correspondingly adjusted BOL/yield, or
        whether RPT4's LAM_Si overshoot (+7.2pp) can be tamed independently
        (e.g. via tau, which does not affect final magnitude much but does
        affect how the rise is paced between RPT4 and RPT5) without
        disturbing the now-excellent RPT5 match.
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import pybamm
from pybamm.input.parameters.lithium_ion.si_gr_expansion import (
    silicon_LGM50_diffusivity_Bonkile2024,
    silicon_LGM50_electrolyte_exchange_current_density_Chen2020,
)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PARAMS_DIR = os.path.join(SCRIPT_DIR, "..", "parameters")
EXPERIMENTAL_DATA_DIR = os.path.join(SCRIPT_DIR, "experimental_data")
sys.path.insert(0, PARAMS_DIR)
from cell064_parameters import get_parameter_values  # noqa: E402

# RPT1 (EFC50.6) reads IDENTICAL to RPT0/BoL (100% SoH) -- nothing in the
# real capacity-fade data distinguishes EFC 0-50 from t=0 itself, i.e. the
# experimental record only starts resolving real aging from ~EFC50 onward.
# The model, by contrast, starts aging from its own EFC=0. Rather than treat
# the real x-axis as absolute, shift it back by EFC_START_OFFSET so RPT1
# lands at EFC~0 -- matching model-EFC to "EFC since aging became visible"
# instead of "EFC since cell assembly". This also re-anchors the knee target
# used elsewhere (152 -> ~152-50=102 in model-EFC terms).
EFC_START_OFFSET = 50.59298876281645  # RPT1's exact EFC (CELL064_capacity_fade.csv)

pybamm.set_logging_level("NOTICE")

# v2 (see TUNING STATUS item 16): the BUILT-IN "stress and reaction-driven"
# LAM option (Ai2019/Reniers2019, already in PyBaMM --
# loss_active_material.py's j_stress_reaction = beta_LAM_sei*a_j_sei/F term
# -- no new code) gives Si a post-knee-coupled LAM increase where plain
# "stress-driven" cannot (confirmed via diagnostic: Si's particle stress is
# IDENTICAL pre- and post-knee, so nothing in stress-driven alone tracks the
# knee, whereas Si's SEI-on-cracks current density a_j_sei jumps ~7x right
# at the knee and stays elevated -- exactly the knee-coupled signal needed).
# NOW THE DEFAULT (v2 > v1): only Si's (secondary) option changes -- Gr
# (primary) stays on plain "stress-driven", unaffected. Set
# C64_SI_LAM_OPTION="stress-driven" to fall back to v1 behaviour.
SI_LAM_OPTION = os.environ.get("C64_SI_LAM_OPTION", "stress and reaction-driven")
# Porosity-gated Si-LAM ("test_pathways" v7): a SECOND reaction-driven Si-LAM
# term (loss_active_material.py's new "porosity isolation" nested option),
# independent of SI_BETA_LAM_SEI, that ramps from ~0 while the negative
# electrode still has porosity headroom (vs. its as-set-up value,
# NEG_POROSITY_BOL) to its full rate as porosity approaches
# NEG_POROSITY_FLOOR. Motivation: every v1-v6 test this session found a
# clean Pareto frontier between matching real knee timing and matching the
# real LAM_Si/LLI split when using purely additive, time/overpotential-only
# Si-LAM levers (TIMESCALE_MULT, SI_LAM_PROP, SI_BETA_LAM_SEI, K_SEI_MULT) --
# none could give "small Si-LAM pre-knee, then a sharp jump post-knee" tied
# to the SAME event (porosity hitting its floor) that drives the SoH knee
# itself. This ties the two together directly instead of via independent
# timing knobs. SI_BETA_LAM_ISO=0.0 (default) leaves behaviour identical to
# every prior test (the new term is then always exactly 0, and SI_LAM_OPTION
# is left unmodified below).
SI_BETA_LAM_ISO = float(os.environ.get("C64_SI_BETA_LAM_ISO", 0.0))
if SI_BETA_LAM_ISO != 0.0 and "porosity" not in SI_LAM_OPTION:
    SI_LAM_OPTION = SI_LAM_OPTION + " and porosity isolation"
# Item 23: redirect-to-LAM -- see loss_active_material.py's redirect branch
# and sei_growth.py's matching gate on dcdt_sei. When enabled, a growing
# fraction of Si's SEI reaction current stops forming new SEI film (hence
# stops becoming new LLI) and is redirected into active-material loss
# instead, conserving total reaction current rather than adding a second,
# independent consumption channel on top of unchanged SEI growth. Requires
# "porosity" already in SI_LAM_OPTION (the redirect reuses that mechanism's
# isolation-gate state).
SI_REDIRECT_TO_LAM = os.environ.get("C64_SI_REDIRECT_TO_LAM", "0") == "1"
SI_REDIRECT_LAM_YIELD = float(os.environ.get("C64_SI_REDIRECT_LAM_YIELD", 1.0))
# SEI_REDIRECT_TO_LAM_OPTION is the GLOBAL model option (gates the redirect
# mechanism in both sei_growth.py and loss_active_material.py); each
# phase's OWN "porosity" membership in its LAM option decides whether THAT
# phase's SEI current actually gets gated -- see GR_REDIRECT_TO_LAM below,
# defined together with GR_LAM_OPTION.
# Sharpness of the isolation gate's ramp near the porosity floor -- default
# 4.0 matches the exponent this mechanism launched with. User feedback
# (v7t_ksei4e-4_bol0.22): the resulting knee is not just late but too
# SHARP/steep (near-vertical) vs. the real cell's transition -- lower this
# to spread the ramp over more of life (softer knee), higher to concentrate
# it (sharper knee).
SI_LAM_ISO_EXPONENT = float(os.environ.get("C64_SI_LAM_ISO_EXPONENT", 4.0))
# Relaxation time constant [s] for the porosity-isolation gate STATE (see
# loss_active_material.py's isolation-gate-lag doc-comment) -- the gate now
# relaxes toward its algebraic (porosity-driven) target with this timescale
# instead of tracking it instantaneously, giving the post-knee transition
# WIDTH a lever independent of onset timing for the first time. Default 1.0s
# is effectively "no lag" (many orders of magnitude below any cycle
# timescale), reproducing the old instantaneous-gate behaviour unless
# deliberately raised (~1e4-1e6s, i.e. ~1-100 EFC-equivalent, is the
# physically relevant test range given the ~15 EFC natural collapse window).
SI_TAU_LAM_ISO = float(os.environ.get("C64_SI_TAU_LAM_ISO", 1.0))
# Task 2 follow-up: lets the C64_GR_R_SEI_MULT=0 test (which found no
# improvement, and slightly WORSE overall fade -- see task2f) be re-tried
# cleanly. R_sei feeds BOTH (a) base_kinetics.py's per-phase throttling of
# the MAIN intercalation reaction (only active when this option != "none"),
# and (b) sei_growth.py's own SEI-growth self-limiting eta_SEI (always
# active regardless of this option) -- zeroing the R_sei PARAMETER conflates
# both effects. Overriding this OPTION to "none" instead removes only (a),
# isolating whether the main-reaction throttling specifically is what's
# starving graphite.
SEI_FILM_RESISTANCE_OPTION = os.environ.get("C64_SEI_FILM_RESISTANCE", "distributed")

# Item 24: graphite-side redirect-gating (no new LAM_Gr growth -- yield=0).
# Pure LLI (Q_side) kept climbing post-knee even with the redirect active
# for Si, because graphite's OWN SEI reaction has no gate at all (the
# redirect check in sei_growth.py is per-phase: "porosity" has to be in
# THAT PHASE's own LAM option for its isolation_gate state to even exist).
# Giving graphite "porosity isolation" too lets its SEI current be gated by
# the SAME (shared, electrode-level) porosity-collapse signal, but with its
# OWN redirect_LAM_yield=0 so the redirected current simply vanishes
# instead of manufacturing LAM_Gr growth that was never wanted -- this is
# purely a "make graphite's SEI current plateau at the knee" lever, not a
# second active-material-loss channel.
GR_LAM_OPTION = os.environ.get("C64_GR_LAM_OPTION", "stress-driven")
GR_REDIRECT_TO_LAM = os.environ.get("C64_GR_REDIRECT_TO_LAM", "0") == "1"
if GR_REDIRECT_TO_LAM and "porosity" not in GR_LAM_OPTION:
    GR_LAM_OPTION = GR_LAM_OPTION + " and porosity isolation"
GR_LAM_ISO_EXPONENT = float(os.environ.get("C64_GR_LAM_ISO_EXPONENT", 3.0))
GR_TAU_LAM_ISO = float(os.environ.get("C64_GR_TAU_LAM_ISO", 2e8))
GR_REDIRECT_LAM_YIELD = float(os.environ.get("C64_GR_REDIRECT_LAM_YIELD", 0.0))
SEI_REDIRECT_TO_LAM_OPTION = "true" if (SI_REDIRECT_TO_LAM or GR_REDIRECT_TO_LAM) else "false"

# Item 25: Si-OCP aging deformation (U_new = scale*U_base + shift), ramped by
# the SAME isolation_gate already driving Si's LAM/redirect terms -- see
# lithium_ion_parameters.py's ocp_aging_deform_scale/shift doc-comment.
# Defaults are the values derived from CELL064's own measured RPT1->RPT5
# si_ocp_scale/si_ocp_shift_V drift (Si_Gr_Expansion_Precursor/Pouch_Data
# fit_outputs), expressed as the ADDITIONAL scale/shift on top of the
# BOL-deformed curve already baked into this project's Si OCP CSV (which
# already encodes RPT1's own deformation): k_extra = s_V,RPT5/s_V,RPT1 =
# 0.8663, shift_extra = U_off,RPT5 - k_extra*U_off,RPT1 = -0.00908 V.
SI_OCP_AGING_DEFORM = os.environ.get("C64_SI_OCP_AGING_DEFORM", "0") == "1"
SI_OCP_DEFORM_SCALE = float(os.environ.get("C64_SI_OCP_DEFORM_SCALE", 0.8663))
SI_OCP_DEFORM_SHIFT = float(os.environ.get("C64_SI_OCP_DEFORM_SHIFT", -0.00908))
OCP_AGING_DEFORM_OPTION = "true" if SI_OCP_AGING_DEFORM else "false"

# Item 29: Si volume-change aging deformation (t_change exponent evolves
# with LAM fraction) -- see lithium_ion_parameters.py's
# volume_change_deform_exponent_bol/end doc-comment. Defaults are CELL064's
# own fitted Si expansion exponent (Si_Gr_Expansion_Precursor/Pouch_Data,
# joint_estimation.py's si_thickness_ratio_powerlaw, L/L0=1+3*sto^exponent)
# at RPT1 (BOL, 1.561) and RPT5 (most-aged comparison point, 2.188) -- real
# data actually peaks at RPT4 (2.880, near the knee) then relaxes by RPT5,
# a non-monotonic trend this linear-in-LAM-fraction ramp can't fully
# capture, so RPT5 is targeted as the end state (same simplification item
# 25 made for its own end scale/shift).
SI_VOLUME_CHANGE_AGING_DEFORM = os.environ.get("C64_SI_VOLUME_CHANGE_AGING_DEFORM", "0") == "1"
SI_VOLUME_CHANGE_EXP_BOL = float(os.environ.get("C64_SI_VOLUME_CHANGE_EXP_BOL", 1.561216606871758))
SI_VOLUME_CHANGE_EXP_END = float(os.environ.get("C64_SI_VOLUME_CHANGE_EXP_END", 2.18807468386046))
VOLUME_CHANGE_AGING_DEFORM_OPTION = "true" if SI_VOLUME_CHANGE_AGING_DEFORM else "false"

# Item 30 (reinstated 2026-09-19): fixed fraction of eps_s lost to isolation
# still contributes to Si's volume-change signal -- see
# base_mechanics.py/lithium_ion_parameters.py doc comments. Being re-tested
# in combination with the Si diffusivity fix (C64_SI_DIFFUSIVITY_MULT) since
# on its own (pre-diffusivity-fix) it only addressed the secondary eps_s-
# collapse factor while the (suspected diffusion-polarisation-driven)
# sto-range collapse swamped it.
SI_LAM_EXPANSION_RESIDUAL = os.environ.get("C64_SI_LAM_EXPANSION_RESIDUAL", "0") == "1"
SI_LAM_EXPANSION_RESIDUAL_FRAC = float(os.environ.get("C64_SI_LAM_EXPANSION_RESIDUAL_FRAC", 0.9))
LAM_EXPANSION_RESIDUAL_OPTION = "true" if SI_LAM_EXPANSION_RESIDUAL else "false"

# Item 26: thermal expansion. base_mechanics.py already adds
# alpha_T_cell*(T_xav - T_ref) (Ai2019 eq 13) to "Cell thickness change [m]"
# UNCONDITIONALLY -- but under the default "isothermal" thermal option T_xav
# is pinned at T_init, so this term is always exactly zero regardless of
# alpha_T_cell's value. This is the only mechanism that can push k =
# delta_rev,cell/delta_rev,particle ABOVE 1 (pore buffering can only reshape
# the particle-derived signal, never exceed it) -- real CELL064 data shows
# k rising to ~1.5 post-knee, which no (f0, width) combination could ever
# reach (see item 17's doc comment). CAUTION (flagged by user): si_gr_
# expansion.py's diffusivity/exchange-current/cracking-rate functions (Gr,
# Si, NMC) and its electrolyte diffusivity/conductivity functions ALL have
# their own Arrhenius temperature dependence -- but as hardcoded Python
# constants inside the function bodies (D_ref/E_D_s, m_ref/E_r, Eac_cr,
# E_D_c_e, E_sigma_e), NOT pybamm.Parameter() entries, so they CANNOT be
# zeroed via a parameter override the way the two SEI growth activation
# energies (also 38000 J/mol, but genuine Parameters) can. Switching to
# "lumped" therefore unavoidably reactivates ALL of these simultaneously
# the moment T departs from 298.15 K, not just the intended mechanical
# expansion term -- this needs measuring (how much does T actually rise
# post-knee?) before being treated as a safe, isolated addition.
THERMAL_OPTION = os.environ.get("C64_THERMAL_OPTION", "isothermal")

MODEL_OPTIONS_BASE = {
    "particle phases": ("2", "1"),
    "open-circuit potential": (("single", "current sigmoid"), "single"),
    "SEI": "ec reaction limited",
    "SEI porosity change": "true",
    "particle mechanics": ("swelling and cracking", "swelling only"),
    "SEI on cracks": "true",
    "loss of active material": ((GR_LAM_OPTION, SI_LAM_OPTION), "stress-driven"),
    "pore buffering": "true",
    "pore buffering transition": "physical",
    "SEI film resistance": SEI_FILM_RESISTANCE_OPTION,
    "SEI reaction redirect to LAM": SEI_REDIRECT_TO_LAM_OPTION,
    "open-circuit potential aging deformation": OCP_AGING_DEFORM_OPTION,
    "volume change aging deformation": VOLUME_CHANGE_AGING_DEFORM_OPTION,
    "active material expansion residual": LAM_EXPANSION_RESIDUAL_OPTION,
    "thermal": THERMAL_OPTION,
}

VAR_PTS = {
    "x_n": 5, "x_s": 5, "x_p": 5,
    "r_n_prim": 20, "r_n_sec": 20, "r_n": 20, "r_p": 20,
}

NOMINAL_CAP_AH = 2.5947  # cell064_parameters.get_parameter_values()'s value
UPPER_CUTOFF_V = 4.195
LOWER_CUTOFF_V = 2.5

# ---------------------------------------------------------------------------
# Degradation recipe -- same values as
# ../../test_model/degradation_test_matrix/test_DMA/dma_baseline_run.py (the
# short, un-stretched "ec reaction limited" recipe this project's DMA
# pipeline was actually validated against), kept unchanged except for
# K_SEI_MULT -- see module docstring. Env-var overrides (C64_SI_CRIT_STRESS
# etc.) let a calibration probe move any of these without editing the file;
# defaults below are the reference recipe's own values.
# ---------------------------------------------------------------------------
SI_CRIT_STRESS = float(os.environ.get("C64_SI_CRIT_STRESS", 2.2e8))   # [Pa] threshold
GR_CRIT_STRESS = float(os.environ.get("C64_GR_CRIT_STRESS", 3.0e7))   # [Pa] threshold
SI_LAM_EXP = float(os.environ.get("C64_SI_LAM_EXP", 2.5))             # shape exponent
GR_LAM_EXP = 2.0                                                      # shape exponent
NEG_POROSITY_FLOOR = float(os.environ.get("C64_NEG_POROSITY_FLOOR", 0.023))  # shape -- items 7,11,12,13; item 19 phase 4 bracketed the joint-fit sweet spot at ~0.023-0.024 (0.024 gives a smoother post-knee + better expansion RMSE but is solver-unstable past EFC~200; 0.023 is stable and validated through RPT4)
EXPONENT_MAX_SEI = float(os.environ.get("C64_EXPONENT_MAX_SEI", 100.0))  # item 19 phase 1-4: 50->100 sharpens the knee; still ~1e15x above the diffusion-limited threshold (physics-safe per item 13) and well below float64 overflow
# Was a silent no-op for this recipe's "ec reaction limited" SEI option (the
# cap was only wired into the "reaction limited" branch of sei_growth.py
# until this project's own fix -- see TUNING STATUS item 13). 10.0 (the
# other recipe's own calibrated value, for a very different k_sei/D_ec
# regime) caps k_exp at only ~4% of the diffusion-limited threshold
# (D_ec/L_sei_0) here -- i.e. it would have suppressed the porosity-floor
# knee mechanism entirely. 50.0 leaves ~1e16x headroom above that threshold
# (so the diffusion-limited saturation this whole recipe depends on is
# never touched) while still capping the exponent argument to exp() at 50,
# nowhere near float64's overflow point (~709) -- see item 13.
# Time-stretch v2 (see TUNING STATUS item 9): at NEG_POROSITY_BOL=0.25, the
# UN-stretched reference rates below knee'd around model-EFC~60 (real target,
# EFC_START_OFFSET-adjusted, is ~101) -- these three are divided by a first-
# pass stretch factor of 1.68 (101/60) from the item-3 values (K_SEI_MULT
# 1.78e-3, SI_LAM_PROP 4.425e-8, GR_LAM_PROP 6.66e-7). Not yet re-confirmed
# where the knee actually lands at this stretch -- may need a second pass.
# C64_TIMESCALE_MULT: pure KNEE-TIMING knob (item 19 phase 2+). Scales the
# three degradation-rate constants -- K_SEI_MULT, LAM_PROP_MULT and
# SI_BETA_LAM_SEI -- by ONE common factor, per item 9's decomposition:
# uniformly scaling all the rates together shifts the knee EFC while
# preserving its shape/depth. >1 => faster fade => EARLIER knee. Phase 1
# moved these three independently and so never had a clean timing lever.
TIMESCALE_MULT = float(os.environ.get("C64_TIMESCALE_MULT", 1.0))
# C64_LAM_PROP_MULT: joint multiplier on BOTH stress-driven LAM proportional
# rates below (Si and Gr together), applied on top of whatever absolute
# values C64_SI_LAM_PROP / C64_GR_LAM_PROP resolve to. Added for the item-19
# joint sweep -- lets "LAM rate" be one sweep axis without having to move
# the two absolute knobs in lockstep.
LAM_PROP_MULT = float(os.environ.get("C64_LAM_PROP_MULT", 0.6)) * TIMESCALE_MULT   # item 19 phase 2-4: 1.0->0.6 (the joint sweep's best combos all use the lower stress-driven LAM rate; the reaction-driven term now carries more of the Si LAM)
SI_LAM_PROP = float(os.environ.get("C64_SI_LAM_PROP", 3.078e-8)) * LAM_PROP_MULT   # [s-1] Si LAM proportional term -- item 11 (+60% over item 9's stretch)
GR_LAM_PROP = float(os.environ.get("C64_GR_LAM_PROP", 4.634e-7)) * LAM_PROP_MULT   # [s-1] Gr LAM proportional term -- item 11 (+60% over item 9's stretch)
SI_CRACK_RATE_MULT = float(os.environ.get("C64_SI_CRACK_MULT", 0.1))
GR_CRACK_RATE_MULT = float(os.environ.get("C64_GR_CRACK_MULT", 0.1))
# v2 (item 16) -- unused unless SI_LAM_OPTION includes "reaction". No
# calibrated non-zero value exists anywhere in PyBaMM's own parameter sets
# to anchor against (every built-in set ships this at 0.0); calibrated here
# via a 7-point beta sweep (1e-6 to 15e-6) scored against score_rpt_gap()'s
# RPT-vs-RPT mean |gap|: 8.0/7.2/5.4/4.2/3.9/4.1/2.0pp respectively -- 15e-6
# is the clear best (RPT4=+1.2pp, RPT5=+2.2pp), still slightly overshooting
# rather than crossing to undershoot, with LAM_Si reaching 54.2% (vs. v1's
# ~12-33%) -- now at the edge of the crude real estimate (58-80%).
SI_BETA_LAM_SEI = float(os.environ.get("C64_SI_BETA_LAM_SEI", 7e-6)) * TIMESCALE_MULT   # [m3.mol-1] -- item 19 phase 2-4 lowered 15e-6->7e-6; the joint sweep's near-negligible main effect for this knob meant the lower value (less aggressive post-knee Si LAM) won on stability without hurting the fit
BASE_CRACK_RATE = 3.9e-20
# Stage 3 (expansion fitting, item 17): f0 x width 3x3 grid (f0 in
# {0.7, 0.55, 0.4} x width in {0.01, 0.05, 0.12}), scored via
# score_expansion_shape()/score_k_gap() against normalised reversible
# expansion + k. f0=0.7 (unchanged from the inherited si_gr_expansion
# default) dominated every width -- lowering f0 monotonically hurt both
# metrics here (opposite of test_model's own f0-lowering sweep, which was
# tuned against a different recipe/target). Within f0=0.7: width=0.01
# (old default) shape=0.151/k_gap=0.330; width=0.05 shape=0.096/k_gap=0.226
# (best shape, close 2nd on k); width=0.12 shape=0.128/k_gap=0.205 (best k,
# close 2nd on shape). ADOPTED width=0.05 -- shape was the primary named
# target and its margin over width=0.12 (25%) is much larger than width=
# 0.12's margin over width=0.05 on k (9%). Real k's post-knee values
# (~1.5) exceed the "physical" pore-buffering formulation's structural
# ceiling of k<=1 (dv_thickness = dv_solid - dv_buffered <= dv_solid
# whenever dv_solid>0) -- no (f0, width) combination can reach that regime;
# per instruction, not chased further for now.
F0_BASELINE = float(os.environ.get("C64_F0", 0.7))
WIDTH_BASELINE = float(os.environ.get("C64_WIDTH", 0.05))
# "Closure porosity" (eps_min_transfer) -- a THIRD, independent pore-
# buffering lever never yet swept this session: sets the porosity value at
# which pore-buffering itself closes/saturates (percolation closure),
# distinct from F0 (transmitted-fraction plateau LEVEL) and WIDTH
# (transition SHARPNESS). Literature default 0.08 sits between our BOL
# (0.130-ish) and the isolation-gate floor (0.035), so it closes at a
# porosity reached BEFORE the SEI-isolation knee -- a genuinely independent
# timing channel from the redirect mechanism's own knee.
CLOSURE_POROSITY_BASELINE = float(os.environ.get("C64_CLOSURE_POROSITY", 0.08))
# Quick sensitivity test: electrolyte concentration multiplier (baseline
# 1000 mol/m3). Not tied to any specific investigation item -- just a
# direct multiplier for a fast exploratory check.
ELECTROLYTE_CONC_MULT = float(os.environ.get("C64_ELECTROLYTE_CONC_MULT", 1.0))

# Stage 3 (expansion fitting) starts here. "Number of electrodes connected
# in parallel to make a cell" was the first idea tried for scaling the
# model's single-representative-electrode-pair thickness change up to the
# real pouch cell's measured (dilatometry-style, whole-stack) expansion
# magnitude -- REJECTED before ever running it, on tracing the parameter
# (per instruction, checked whether this would touch capacity/degradation
# before testing): geometric_parameters.py's
# "A_cc = L_y * L_z * n_electrodes_parallel" makes this the CURRENT-
# COLLECTOR AREA, which drives Q_init (capacity), every C-rate, and every
# current-density term throughout the electrochemistry -- not a thickness-
# only knob. Scaling it up ~8x to match the expansion magnitude would also
# inflate the model's effective capacity ~8x, silently invalidating the
# whole degradation fit (items 1-16). No physical layer count is documented
# anywhere in the source project's data either (checked --
# CELL064_DATA_REPORT.md, the BoL parameter-set scripts, and a broad grep
# across Pouch_Data/ for "layer"/"electrode pair"/"stack count" all came up
# empty -- the source project's own geometry_solver.py deliberately solves
# a single electrode-pair's area so ITS OWN capacity alone already matches
# the real cell, an abstraction that discards physical layer-count
# information entirely). ADOPTED INSTEAD (per instruction): compare
# NORMALISED reversible expansion (each series divided by its own
# reference value) rather than absolute magnitude -- sidesteps the
# geometry coupling entirely, no parameter change needed. See
# _plot_reversible_expansion's normalise= option and TUNING STATUS.

# BoL negative electrode porosity override -- default 0.25 (was 0.35
# [SOLVED GEOMETRY] in cell064_parameters.py; see TUNING STATUS item 9).
# Porosity and the Primary/Secondary active-material volume fractions do NOT
# need to sum to 1 -- the remainder is physically the binder/conductive-
# additive fraction (inert, doesn't participate in OCV/capacity), so
# reducing porosity ALONE (leaving both active-material fractions untouched)
# does not perturb the C/20 BoL capacity/voltage/dV-dQ fit (verified in
# ../_porosity_bol_check.py: capacity identical to 4dp, max |dV| 1.02 mV at
# 0.25 -- porosity only enters the electrolyte-transport equations,
# negligible at C/20's very low rate). 0.25 was the best single point tried
# in a 0.20/0.25/0.30 bracket: pushes the porosity-floor-driven knee onset
# to ~91-92% SoH (matching real RPT3=91.4%), vs. ~85-86% at the original
# 0.35 -- see item 9. Set to None (or unset C64_NEG_POROSITY_BOL) to fall
# back to cell064_parameters.py's own unmodified 0.35.
NEG_POROSITY_BOL = os.environ.get("C64_NEG_POROSITY_BOL", "0.25")
NEG_POROSITY_BOL = float(NEG_POROSITY_BOL) if NEG_POROSITY_BOL else None

# SEI kinetic rate multiplier. The reference recipe's own anchor
# (dma_baseline_run.py's K_SEI_ANCHOR_MULT=0.0017) was ~30-50x too fast for
# CELL064's real composition (manual probes: 8x and 40x cuts both collapsed
# to ~45-55% SoH within ~70-100 cycles). 2e-5 (~85x below the anchor) was
# settled on manually before the parallel sweep (see
# ../sweeps/cell064_k_sei_lam_sweep.py) confirmed it doesn't matter much
# further: sweeping D_ec (EC diffusivity) 4 orders of magnitude at this
# k_sei value changed the capacity-fade RMSE by <0.1 percentage point --
# unlike si_gr_expansion's own composition (where D_ec had real leverage per
# ../../test_model/degradation_test_matrix/degradation_test_matrix_plan.md
# section 1a), CELL064's real geometry/composition makes k_sei the only
# lever that matters on the SEI/LLI side.
BASE_K_SEI = 1e-12
K_SEI_MULT = float(os.environ.get("C64_K_SEI_MULT", 1.0e-3)) * TIMESCALE_MULT  # item 19 phase 4: 7.74e-4->1.0e-3. This is the KNEE-TIMING lever -- 7.74e-4 knee'd at EFC~133, 1.0e-3 at EFC~108 (real anchor ~101). Higher still (>=1.24e-3) lands the knee right but crashes ~EFC200; 1.0e-3 with floor 0.023 is the stablest config that reaches the ~101 anchor and still passes through real RPT4/RPT5. C64_TIMESCALE_MULT scales this + LAM + beta together.
# Task 2 follow-up (test_pathways, task2g): per-phase override of K_SEI_MULT
# -- lets Gr's own BULK SEI kinetic rate constant be tested/reduced
# independently of Si's, complementing task2d's Gr CRACK-rate probe (which
# found no effect on graphite's stoichiometry collapse). Both default to
# K_SEI_MULT's own env value (C64_K_SEI_MULT) so behaviour is unchanged
# unless one is set explicitly.
_K_SEI_MULT_BASE = float(os.environ.get("C64_K_SEI_MULT", 1.0e-3))
GR_K_SEI_MULT = float(os.environ.get("C64_GR_K_SEI_MULT", _K_SEI_MULT_BASE)) * TIMESCALE_MULT
SI_K_SEI_MULT = float(os.environ.get("C64_SI_K_SEI_MULT", _K_SEI_MULT_BASE)) * TIMESCALE_MULT

# Task 2 follow-up (test_pathways): "SEI film resistance" is active by
# default here ("distributed", auto-enabled whenever SEI != "none") and adds
# eta_sei = -j_tot * L_sei * R_sei directly to each PHASE's OWN main
# intercalation-reaction overpotential (base_kinetics.py), not just to SEI
# growth's own self-limiting kinetics. si_gr_expansion.py sets an EQUAL
# R_sei=200000 Ohm.m for both Primary(Gr)/Secondary(Si) -- but since Gr's
# own combined SEI channels (bulk+cracks) are much larger than Si's (task2c:
# ~36% vs ~19% of nominal capacity), Gr's SEI film L_sei has likely grown
# correspondingly thicker, so this PER-PHASE resistance penalty throttles
# Gr's own charge/discharge current disproportionately even with identical
# R_sei -- a candidate explanation for task2e's finding that Gr's own
# stoichiometry window collapses to its floor early in the post-knee
# discharge (not from LAM_Gr, which stays tiny, and not fixed by suppressing
# crack_Gr alone per task2d). Independent per-phase multipliers so Gr's
# contribution can be tested in isolation from Si's.
BASE_R_SEI = 200000.0
GR_R_SEI_MULT = float(os.environ.get("C64_GR_R_SEI_MULT", 1.0))
SI_R_SEI_MULT = float(os.environ.get("C64_SI_R_SEI_MULT", 1.0))

# Si solid diffusivity multiplier (theory test, session 2026-09-19): the base
# value (silicon_LGM50_diffusivity_Bonkile2024, D_ref=3.0e-16 m2.s-1) is a
# LITERATURE placeholder never fitted against CELL064 itself -- no high-C-rate
# pouch data exists to fit it against. Hypothesis: as Si LAM shrinks eps_s,
# the remaining active Si sees higher effective current density, and if the
# real Si diffusivity is higher than this placeholder, the model
# over-estimates solid-diffusion polarisation, prematurely hitting the
# voltage cutoff and artificially truncating Si's own per-cycle stoichiometry
# excursion (the ~14.3x delta_sto collapse) well before the real accessible
# capacity is used. Default 1.0 = unchanged baseline.
SI_DIFFUSIVITY_MULT = float(os.environ.get("C64_SI_DIFFUSIVITY_MULT", 1.0))


def _silicon_diffusivity_scaled(sto, T):
    return SI_DIFFUSIVITY_MULT * silicon_LGM50_diffusivity_Bonkile2024(sto, T)


# Si exchange-current-density multiplier (theory test, session 2026-09-19,
# user's own follow-up): same "never fitted, no high-C-rate data" caveat as
# the diffusivity placeholder above -- a too-low exchange-current density
# means excessive CHARGE-TRANSFER (Butler-Volmer) overpotential at a given
# current, on top of (not instead of) the solid-diffusion overpotential the
# diffusivity fix addresses. As Si LAM shrinks the remaining active area,
# both mechanisms are driven harder by the same rising local current
# density, so this could be contributing to the same premature-cutoff /
# truncated-sto-window effect independently of diffusivity. Default 1.0 =
# unchanged baseline.
SI_EXCHANGE_CURRENT_MULT = float(os.environ.get("C64_SI_EXCHANGE_CURRENT_MULT", 1.0))


def _silicon_exchange_current_scaled(c_e, c_s_surf, c_s_max, T):
    return SI_EXCHANGE_CURRENT_MULT * silicon_LGM50_electrolyte_exchange_current_density_Chen2020(
        c_e, c_s_surf, c_s_max, T
    )


# Si particle radius multiplier (theory test, session 2026-09-19, user's own
# follow-up): the base value (1.52e-6 m, si_gr_expansion.py literature
# default) carries the same "[LITERATURE - TO FIT]" tag as the diffusivity/
# exchange-current placeholders -- never fitted to CELL064. Solid-diffusion
# time constant tau ~ R^2/D, so radius enters QUADRATICALLY -- a candidate
# for a much stronger lever per-unit-change than diffusivity (linear in D)
# if the real particles are smaller than this placeholder assumes. Smaller
# radius -> shorter diffusion length -> less surface-concentration
# polarisation for the same local current density, independent of both the
# diffusivity and exchange-current fixes above. Default 1.0 = unchanged
# baseline (1.52e-6 m).
BASE_SI_PARTICLE_RADIUS = 1.52e-06
SI_PARTICLE_RADIUS_MULT = float(os.environ.get("C64_SI_PARTICLE_RADIUS_MULT", 1.0))


PARAM_UPDATES = {
    "Initial concentration in electrolyte [mol.m-3]": 1000.0 * ELECTROLYTE_CONC_MULT,
    # Needed only for THERMAL_OPTION="lumped" (composite negative electrode
    # has no single density otherwise -- lumped thermal's heat-capacity calc
    # needs the plain, non-phase-specific key). Computed as CELL064's own
    # active-volume-fraction-weighted average of the Primary (graphite,
    # 1657 kg/m3) and Secondary (silicon, 2650 kg/m3) literature densities:
    # (0.592132*1657 + 0.0578685*2650) / (0.592132+0.0578685) = 1745.4.
    # Harmless/unused when thermal="isothermal" (the default).
    "Negative electrode density [kg.m-3]": 1745.4,
    "Secondary: Negative electrode LAM constant proportional term [s-1]": SI_LAM_PROP,
    "Secondary: Negative electrode LAM constant exponential term": SI_LAM_EXP,
    "Secondary: Negative electrode critical stress [Pa]": SI_CRIT_STRESS,
    "Secondary: Negative electrode cracking rate": BASE_CRACK_RATE * SI_CRACK_RATE_MULT,
    "Secondary: Negative electrode reaction-driven LAM factor [m3.mol-1]": SI_BETA_LAM_SEI,  # v2 default, item 16 -- calibrated 15e-6
    "Secondary: Negative electrode porosity-isolation LAM factor [s-1.m-1]": SI_BETA_LAM_ISO,  # v7c -- drives SEI THICKNESS now, not a_j_sei; only active when SI_LAM_OPTION includes "porosity isolation"
    "Secondary: Negative electrode porosity-isolation LAM exponent": SI_LAM_ISO_EXPONENT,  # v7 -- knee sharpness/width, default 4.0
    "Secondary: Negative electrode porosity-isolation gate time constant [s]": SI_TAU_LAM_ISO,  # isolation-gate-lag test -- default 1.0 (no lag)
    "Secondary: Negative electrode SEI-redirect-to-LAM yield": SI_REDIRECT_LAM_YIELD,  # item 23 -- only active when C64_SI_REDIRECT_TO_LAM=1
    "Secondary: Negative electrode OCP aging-deformation end scale": SI_OCP_DEFORM_SCALE,  # item 25 -- only active when C64_SI_OCP_AGING_DEFORM=1
    "Secondary: Negative electrode OCP aging-deformation end shift [V]": SI_OCP_DEFORM_SHIFT,  # item 25
    "Secondary: Negative electrode volume change aging-deformation BOL exponent": SI_VOLUME_CHANGE_EXP_BOL,  # item 29
    "Secondary: Negative electrode volume change aging-deformation end exponent": SI_VOLUME_CHANGE_EXP_END,  # item 29
    "Secondary: Negative electrode LAM expansion residual fraction": SI_LAM_EXPANSION_RESIDUAL_FRAC,  # item 30 -- only active when C64_SI_LAM_EXPANSION_RESIDUAL=1
    "Primary: Negative electrode porosity-isolation LAM factor [s-1.m-1]": 0.0,  # unused by the redirect branch (only fetched, never applied) but still needs a value once GR_LAM_OPTION includes "porosity isolation"
    "Primary: Negative electrode porosity-isolation LAM exponent": GR_LAM_ISO_EXPONENT,  # item 24 -- graphite's own gate sharpness, shares the same electrode-level porosity signal as Si's
    "Primary: Negative electrode porosity-isolation gate time constant [s]": GR_TAU_LAM_ISO,  # item 24
    "Primary: Negative electrode SEI-redirect-to-LAM yield": GR_REDIRECT_LAM_YIELD,  # item 24 -- default 0.0: redirected current vanishes instead of growing LAM_Gr
    "Primary: Negative electrode LAM constant proportional term [s-1]": GR_LAM_PROP,
    "Primary: Negative electrode LAM constant exponential term": GR_LAM_EXP,
    "Primary: Negative electrode critical stress [Pa]": GR_CRIT_STRESS,
    "Primary: Negative electrode cracking rate": BASE_CRACK_RATE * GR_CRACK_RATE_MULT,
    "Negative electrode porosity floor": NEG_POROSITY_FLOOR,
    "Primary: SEI reaction exponent cap": EXPONENT_MAX_SEI,
    "Secondary: SEI reaction exponent cap": EXPONENT_MAX_SEI,
    "Negative electrode pore buffering transition width": WIDTH_BASELINE,
    "Negative electrode transmitted fraction plateau": F0_BASELINE,
    "Negative electrode pore buffering closure porosity": CLOSURE_POROSITY_BASELINE,
    "Primary: SEI kinetic rate constant [m.s-1]": BASE_K_SEI * GR_K_SEI_MULT,
    "Secondary: SEI kinetic rate constant [m.s-1]": BASE_K_SEI * SI_K_SEI_MULT,
    "Primary: SEI resistivity [Ohm.m]": BASE_R_SEI * GR_R_SEI_MULT,
    "Secondary: SEI resistivity [Ohm.m]": BASE_R_SEI * SI_R_SEI_MULT,
    "Secondary: Negative particle diffusivity [m2.s-1]": _silicon_diffusivity_scaled,
    "Secondary: Negative electrode exchange-current density [A.m-2]": _silicon_exchange_current_scaled,
    "Secondary: Negative particle radius [m]": BASE_SI_PARTICLE_RADIUS * SI_PARTICLE_RADIUS_MULT,
}
if NEG_POROSITY_BOL is not None:
    PARAM_UPDATES["Negative electrode porosity"] = NEG_POROSITY_BOL

# Env-var overrides for fast calibration probes without editing the file.
BATCH_SIZE = int(os.environ.get("C64_BATCH_SIZE", 50))
MAX_TOTAL_CYCLES = int(os.environ.get("C64_MAX_CYCLES", 250))  # item 19: lowered 450->250. Nothing to fit past EFC~200 (real RPT5 is EFC~197 shifted), and the long tail-runs to 45% SoH accumulate PyBaMM solution history unboundedly -> ~24 GB/run, OOM'd a 3-parallel sweep. 250 requested cycles covers RPTs at 50/100/150/200/250 and still reaches real RPT4/RPT5's EFC range even with a few retry-stalled batches (item 11's reason for 450); raise via C64_MAX_CYCLES only for a deliberate full-life run.
RPT_INTERVAL = int(os.environ.get("C64_RPT_INTERVAL", 50))
RPT_RATE = "C/20"
assert BATCH_SIZE % RPT_INTERVAL == 0, "BATCH_SIZE must be a multiple of RPT_INTERVAL"
SOH_TERMINATION_PERCENT = float(os.environ.get("C64_SOH_FLOOR", 45.0))  # a little past the real data's lowest point (50.0% @ EFC 248); override lets a one-off extended run push past it (post-knee EFC accumulates slower than cycle count, roughly proportional to SoH, so reaching the model's own last RPT past real EFC~200 needs the floor pushed down too, not just C64_MAX_CYCLES)
EFC_LIMIT = float(os.environ.get("C64_EFC_LIMIT", "inf"))  # optional third stop condition alongside MAX_TOTAL_CYCLES/SOH_TERMINATION_PERCENT: some configs (e.g. very low GR_K_SEI_MULT) fade so slowly that neither the cycle budget nor the SoH floor triggers anywhere near the real RPT4/RPT5 EFC range, wasting a long tail of cycles past the EFC we actually score against. Checked every batch via the same "Throughput capacity [A.h]" summary variable efc_from_throughput() already uses, so it's directly comparable to the real data's EFC axis.

# Tried tightening this (1e-7, 1e-8) per instruction, to check whether the
# loose tolerance is behind the occasional visibly-off capacity points near
# the knee. Result: both were ~4-5x slower per cycle than 1e-6 (~6.5s/cycle
# vs ~1.2-1.7s) with NO further slowdown between 1e-7 and 1e-8 -- the cost
# jumps as soon as you leave 1e-6 at all, it doesn't scale smoothly with
# tolerance. Reverted to the original 1e-06 per instruction; the "weird
# points" issue is still open and not yet attributed to solver tolerance.
SOLVER_TOL = float(os.environ.get("C64_SOLVER_TOL", 1e-06))
solver = pybamm.IDAKLUSolver(root_tol=SOLVER_TOL, atol=SOLVER_TOL, rtol=SOLVER_TOL)

formation_exp = pybamm.Experiment(
    [
        # Rate-matched to RPT_RATE (not an arbitrary "0.1C"): rpt_leg()/
        # rpt_discharge_curve() classify ANY 0.05-0.5 A discharge step as an
        # "RPT-type" leg (see their mean_I window below), so this formation
        # discharge IS what extract_results() treats as the model's "RPT1"
        # for voltage_shape_rmse. It was previously "0.1C" (0.2595 A) vs.
        # every later RPT's true C/20 (0.1297 A) and the real RPT1 data's
        # own C/20 measurement -- a rate mismatch, not a degradation-physics
        # gap, that alone accounted for ~56mV of RPT1's RMSE (task2m
        # follow-up investigation, sweeps/bol_voltage_fit/): fixing just the
        # rate (no material/degradation parameter touched) cuts it to
        # ~30mV, confirmed via sweeps/bol_voltage_fit/check_formation_rate_fix.py.
        f"Discharge at {RPT_RATE} until {LOWER_CUTOFF_V} V",
        f"Charge at C/3 until {UPPER_CUTOFF_V} V",
        f"Hold at {UPPER_CUTOFF_V} V until C/50",
    ]
)


def _batch_cycle_steps(local_position):
    if local_position % RPT_INTERVAL == 0:
        return (
            f"Discharge at {RPT_RATE} until {LOWER_CUTOFF_V} V",
            f"Charge at C/3 until {UPPER_CUTOFF_V} V",
            f"Hold at {UPPER_CUTOFF_V} V until C/50",
        )
    return (
        f"Discharge at C/3 until {LOWER_CUTOFF_V} V",
        f"Charge at C/3 until {UPPER_CUTOFF_V} V",
        f"Hold at {UPPER_CUTOFF_V} V until C/50",
    )


ageing_batch_exp = pybamm.Experiment(
    [_batch_cycle_steps(i) for i in range(1, BATCH_SIZE + 1)]
)


def build_parameter_values():
    param = pybamm.ParameterValues(get_parameter_values())
    param.update(PARAM_UPDATES, check_already_exists=False)
    return param


# A post-crash cycle sometimes leaves a near-instant, near-zero-capacity
# "discharge" step in the solution (visible as a vertical dip to ~0 A.h in
# every crash-affected plot). Below this floor (~10% of nominal -- well
# under the 45% SOH_TERMINATION_PERCENT floor, so it never masks a genuine
# low-SoH reading) a cycle's capacity is treated as corrupted/incomplete
# rather than real data, so it can't falsely trip the SoH-floor termination
# check or show up as a bogus point on the plots.
MIN_VALID_CAP_AH = 0.10 * NOMINAL_CAP_AH


def cycle_ageing_leg(cyc):
    cap, rate = 0.0, 0.0
    discharge_end_thr, charge_end_thr = None, None
    for step in getattr(cyc, "steps", []):
        try:
            I = step["Current [A]"].entries
            q = step["Discharge capacity [A.h]"].entries
            thr_step = step["Throughput capacity [A.h]"].entries
        except (KeyError, TypeError, AttributeError):
            continue
        if I.size < 2:
            continue
        mean_I = np.mean(I)
        if mean_I > 0.5 and q.size >= 2:
            c = float(q[-1] - q[0])
            if c > cap:
                cap, rate = c, float(mean_I)
                discharge_end_thr = float(thr_step[-1])
        elif mean_I < -1e-3:
            charge_end_thr = float(thr_step[-1])
    if cap < MIN_VALID_CAP_AH:
        return 0.0, 0.0, None, None
    return cap, rate, discharge_end_thr, charge_end_thr


def rpt_leg(cyc):
    cap, rate = 0.0, 0.0
    discharge_end_thr, charge_end_thr = None, None
    for step in getattr(cyc, "steps", []):
        try:
            I = step["Current [A]"].entries
            q = step["Discharge capacity [A.h]"].entries
            thr_step = step["Throughput capacity [A.h]"].entries
        except (KeyError, TypeError, AttributeError):
            continue
        if I.size < 2:
            continue
        mean_I = np.mean(I)
        if 0.05 < mean_I <= 0.5 and q.size >= 2:
            c = float(q[-1] - q[0])
            if c > cap:
                cap, rate = c, mean_I
                discharge_end_thr = float(thr_step[-1])
        elif mean_I < -1e-3:
            charge_end_thr = float(thr_step[-1])
    return cap, rate, discharge_end_thr, charge_end_thr


def efc_from_throughput(thr_ah):
    return np.asarray(thr_ah) / (2.0 * NOMINAL_CAP_AH)


def cycle_expansion_ptp_um(cyc):
    """Peak-to-peak 'Cell thickness change [m]' within one cycle, in um --
    the model-side analogue of the real per-cycle reversible expansion
    amplitude (CELL064_reversible_expansion.csv), which is likewise a
    per-cycle peak-to-peak breathing amplitude, not a cumulative one."""
    try:
        d = cyc["Cell thickness change [m]"].entries
    except (KeyError, TypeError, AttributeError):
        return None
    if d.size < 2:
        return None
    return float(d.max() - d.min()) * 1e6


def cycle_k_peak(cyc):
    """Peak |Negative electrode transfer ratio k| within one cycle -- the
    model-side analogue of the real k = delta_rev,cell / delta_rev,particle
    (CELL064_k_expansion_scale.csv)."""
    try:
        k = cyc["Negative electrode transfer ratio k"].entries
    except (KeyError, TypeError, AttributeError):
        return None
    if k.size < 1:
        return None
    return float(np.max(np.abs(k)))


def rpt_discharge_curve(cyc):
    """Return (q_rel [A.h], v [V]) for cyc's RPT discharge step -- same step
    -selection rule as rpt_leg (highest-capacity step with
    0.05 < mean_I <= 0.5) -- or (None, None) if cyc has no such step."""
    best_cap = 0.0
    best = (None, None)
    for step in getattr(cyc, "steps", []):
        try:
            I = step["Current [A]"].entries
            q = step["Discharge capacity [A.h]"].entries
            v = step["Terminal voltage [V]"].entries
        except (KeyError, TypeError, AttributeError):
            continue
        if I.size < 2:
            continue
        mean_I = np.mean(I)
        if 0.05 < mean_I <= 0.5 and q.size >= 2:
            c = float(q[-1] - q[0])
            if c > best_cap:
                best_cap = c
                best = (q - q[0], v)
    return best


# On a hard (exception-raising) solver failure, retry the SAME batch with a
# different tolerance before giving up -- per instruction: sometimes a
# looser tolerance lets the solver step over a locally-stiff transient it
# can't resolve at the nominal tolerance, and sometimes a tighter one avoids
# a spurious Newton step that overshoots into an infeasible region; which
# direction helps isn't predictable in advance, so both are tried. Each
# multiplier gets a fresh Simulation (solver is bound at construction).
RETRY_TOL_MULTIPLIERS = [10.0, 0.1, 100.0, 0.01]


def run_degradation():
    options = dict(MODEL_OPTIONS_BASE)
    model = pybamm.lithium_ion.DFN(options)
    param = build_parameter_values()

    def make_sim(experiment, this_solver=solver):
        return pybamm.Simulation(
            model, parameter_values=param, experiment=experiment,
            solver=this_solver, var_pts=VAR_PTS,
        )

    sim = make_sim(formation_exp)
    last_sol = sim.solve(initial_soc=1.0)
    print("Formation cycle solved.", flush=True)

    total_cycles_requested = 0
    stop_reason = "reached MAX_TOTAL_CYCLES"
    batch_num = 0
    while total_cycles_requested < MAX_TOTAL_CYCLES:
        batch_num += 1
        new_sol = None
        last_exc = None
        for attempt, tol_mult in enumerate([1.0] + RETRY_TOL_MULTIPLIERS):
            try_tol = SOLVER_TOL * tol_mult
            this_solver = (solver if tol_mult == 1.0
                           else pybamm.IDAKLUSolver(root_tol=try_tol, atol=try_tol, rtol=try_tol))
            sim = make_sim(ageing_batch_exp, this_solver)
            try:
                new_sol = sim.solve(starting_solution=last_sol)
                if attempt > 0:
                    print(f"[BATCH {batch_num}] recovered on retry {attempt} "
                          f"(tol={try_tol:.0e})", flush=True)
                break
            except Exception as exc:  # noqa: BLE001
                last_exc = exc
                continue
        if new_sol is None:
            print(f"[BATCH {batch_num}] solver failure after {total_cycles_requested} "
                  f"successfully completed ageing cycles, all {len(RETRY_TOL_MULTIPLIERS)} "
                  f"tolerance retries exhausted: {type(last_exc).__name__}: {str(last_exc)[:200]}",
                  flush=True)
            stop_reason = f"solver_failure_batch_{batch_num}"
            break

        last_sol = new_sol
        total_cycles_requested += BATCH_SIZE

        if os.environ.get("C64_DIAG_EPS") == "1":
            # Diagnostic (per instruction): is the crash actually caused by
            # eps_solid (active-material volume fraction) integrating
            # through zero unbounded -- the mechanism flagged in
            # ../../submodel_stability_audit.md's row 1 -- rather than (or
            # in addition to) the SEI-exponent overflow fixed in item 13?
            # eps_solid has no bounds= on its Variable (loss_active_material.py),
            # so nothing stops it going negative; a = 3*eps_solid/R -> 0
            # would then force j = i/a -> infinity, a current-focusing
            # singularity distinct from the SEI-exponent one.
            for label, key in [
                ("Gr (primary)", "Negative electrode primary active material volume fraction"),
                ("Si (secondary)", "Negative electrode secondary active material volume fraction"),
            ]:
                try:
                    vals = last_sol[key].entries
                    print(f"[BATCH {batch_num}] eps_solid {label}: min={vals.min():.6f}, "
                          f"final={vals[..., -1].min():.6f} (BoL~0.59/0.058)", flush=True)
                except Exception as exc:  # noqa: BLE001
                    print(f"[BATCH {batch_num}] eps_solid {label}: unavailable ({exc})", flush=True)

        if os.environ.get("C64_DIAG_POROSITY") == "1":
            # Diagnostic: directly confirm porosity is actually pinned at/near
            # NEG_POROSITY_FLOOR right when the remaining post-item-13 crashes
            # (IDA_CONV_FAIL/IDA_ERR_FAIL/IDA_BAD_K) happen -- this has been
            # ASSUMED throughout the porosity-floor tuning (items 7,9,11-13)
            # but never directly checked the way eps_solid was in item 13.
            try:
                vals = last_sol["Negative electrode porosity"].entries
                print(f"[BATCH {batch_num}] porosity Negative: min={vals.min():.6f}, "
                      f"final={vals[..., -1].min():.6f} (floor={NEG_POROSITY_FLOOR:.4f}, BoL=0.25)",
                      flush=True)
            except Exception as exc:  # noqa: BLE001
                print(f"[BATCH {batch_num}] porosity Negative: unavailable ({exc})", flush=True)

        if os.environ.get("C64_DIAG_STRESS") == "1":
            # Diagnostic: does Si's particle hydrostatic stress -- the ONLY
            # thing driving Si's LAM rate under the current "stress-driven"
            # option (no new mechanism, just the existing j_stress_LAM =
            # -beta_LAM*(stress_tensile/stress_critical)**m_LAM term in
            # loss_active_material.py) -- already rise post-knee on its own,
            # via the existing electrochemistry (higher local current
            # density/overpotential as porosity pins at the floor)? If so,
            # sharpening the EXISTING exponent (SI_LAM_EXP/m_LAM) or lowering
            # the EXISTING threshold (SI_CRIT_STRESS) is enough to get more
            # Si LAM post-knee, with no new coupling term needed.
            try:
                st = last_sol["X-averaged negative secondary particle surface tangential stress [Pa]"].entries
                sr = last_sol["X-averaged negative secondary particle surface radial stress [Pa]"].entries
                stress_h = (sr + 2 * st) / 3
                stress_h_max_tensile = max(stress_h.max(), 0.0)
                ratio = stress_h_max_tensile / SI_CRIT_STRESS
                print(f"[BATCH {batch_num}] Si stress_h max(tensile)={stress_h_max_tensile:.3e} Pa, "
                      f"/stress_critical={ratio:.4f}, ratio**m_LAM={ratio ** SI_LAM_EXP:.3e} "
                      f"(stress_critical={SI_CRIT_STRESS:.2e}, m_LAM={SI_LAM_EXP})", flush=True)
            except Exception as exc:  # noqa: BLE001
                print(f"[BATCH {batch_num}] Si stress: unavailable ({exc})", flush=True)

        if os.environ.get("C64_DIAG_A_J_SEI") == "1":
            # Diagnostic: magnitude/trend of Si's SEI-on-cracks volumetric
            # current density (a_j_sei) across life -- the signal that would
            # drive Si's LAM if "loss of active material" were switched to
            # the built-in "stress and reaction-driven" option (see README
            # item 16/v2 test). Read here regardless of which LAM option is
            # active, purely to gauge order of magnitude for
            # C64_SI_BETA_LAM_SEI before/while calibrating it.
            try:
                a_j = last_sol[
                    "X-averaged negative electrode secondary SEI volumetric "
                    "interfacial current density [A.m-3]"
                ].entries
                print(f"[BATCH {batch_num}] Si a_j_sei: min={a_j.min():.4e}, "
                      f"max={a_j.max():.4e} A.m-3", flush=True)
            except Exception as exc:  # noqa: BLE001
                print(f"[BATCH {batch_num}] Si a_j_sei: unavailable ({exc})", flush=True)

        age_cap = []
        for cyc in last_sol.cycles:
            cp, rate, _, _ = cycle_ageing_leg(cyc)
            if cp > 0 and rate > 0.5:
                age_cap.append(cp)
        if age_cap:
            soh = 100.0 * age_cap[-1] / age_cap[0]
            print(f"[BATCH {batch_num}] {total_cycles_requested} ageing cycles requested, "
                  f"{len(age_cap)} completed, SoH={soh:.2f}%", flush=True)
            if soh <= SOH_TERMINATION_PERCENT:
                stop_reason = f"reached_{SOH_TERMINATION_PERCENT:.0f}pct_floor"
                print("Reached SoH floor -- stopping.", flush=True)
                break
        else:
            print(f"[BATCH {batch_num}] {total_cycles_requested} ageing cycles requested, "
                  "no completed ageing legs found yet", flush=True)

        if np.isfinite(EFC_LIMIT):
            current_efc = float(efc_from_throughput(
                last_sol["Throughput capacity [A.h]"].entries[-1]))
            if current_efc >= EFC_LIMIT:
                stop_reason = f"reached_efc_limit_{EFC_LIMIT:.0f}"
                print(f"[BATCH {batch_num}] EFC={current_efc:.1f} >= C64_EFC_LIMIT="
                      f"{EFC_LIMIT:.1f} -- stopping.", flush=True)
                break

    print(f"Stop reason: {stop_reason}", flush=True)
    return last_sol


def extract_results(sol):
    age_cyc, age_cap, age_thr, age_expansion_um = [], [], [], []
    rpt_cap_list, rpt_thr, rpt_k_list, rpt_voltage_curves = [], [], [], []
    for i, cyc in enumerate(sol.cycles):
        cap, rate, discharge_end_thr, charge_end_thr = cycle_ageing_leg(cyc)
        if cap > 0 and discharge_end_thr is not None:
            age_cyc.append(i)
            age_cap.append(cap)
            age_thr.append(float(cyc["Throughput capacity [A.h]"].entries[-1]))
            exp_ptp = cycle_expansion_ptp_um(cyc)
            age_expansion_um.append(exp_ptp if exp_ptp is not None else np.nan)
            continue
        rcap, rrate, rdischarge_end, rcharge_end = rpt_leg(cyc)
        if rcap > 0 and rdischarge_end is not None:
            rpt_cap_list.append(rcap)
            rpt_thr.append(float(cyc["Throughput capacity [A.h]"].entries[-1]))
            k_peak = cycle_k_peak(cyc)
            rpt_k_list.append(k_peak if k_peak is not None else np.nan)
            q_rel, v = rpt_discharge_curve(cyc)
            if q_rel is not None:
                rpt_voltage_curves.append({"thr": rpt_thr[-1], "q": q_rel, "v": v})

    age_cap = np.array(age_cap)
    age_thr = np.array(age_thr)
    age_expansion_um = np.array(age_expansion_um)
    rpt_cap_arr = np.array(rpt_cap_list)
    rpt_thr = np.array(rpt_thr)
    rpt_k_arr = np.array(rpt_k_list)

    Qt_full = sol["Throughput capacity [A.h]"].entries
    cutoff_idx = np.searchsorted(Qt_full, age_thr[-1], side="right") if age_thr.size else len(Qt_full)

    def full_res_var(name):
        return sol[name].entries[:cutoff_idx]

    Qt = Qt_full[:cutoff_idx]
    LAM_neg = full_res_var("Loss of active material in negative electrode [%]")
    LAM_pos = full_res_var("Loss of active material in positive electrode [%]")
    LAM_gr = full_res_var("Loss of active material in primary phase in negative electrode [%]")
    LAM_si = full_res_var("Loss of active material in secondary phase in negative electrode [%]")
    Q_side = full_res_var("Total capacity lost to side reactions [A.h]")

    # LAM-trapped lithium (Sulzer et al. 2021 eq 37's lli_due_to_lam channel,
    # already tracked by loss_active_material.py but not previously folded
    # into the "LLI" number plotted/compared against the real DMA target --
    # see check_corrected_lli.py's finding: this under-counts the model's
    # true comparable LLI by ~20pp). F = Faraday constant [C/mol].
    F_CONST = 96485.33212
    lli_lam_neg_primary_Ah = full_res_var(
        "Loss of lithium due to loss of primary active material in negative electrode [mol]") * F_CONST / 3600.0
    lli_lam_neg_secondary_Ah = full_res_var(
        "Loss of lithium due to loss of secondary active material in negative electrode [mol]") * F_CONST / 3600.0
    lli_lam_pos_Ah = full_res_var(
        "Loss of lithium due to loss of active material in positive electrode [mol]") * F_CONST / 3600.0
    Q_lli_lam_trap = lli_lam_neg_primary_Ah + lli_lam_neg_secondary_Ah + lli_lam_pos_Ah

    return dict(
        age_cap=age_cap, age_thr=age_thr, rpt_cap_arr=rpt_cap_arr, rpt_thr=rpt_thr,
        Qt=Qt, LAM_neg=LAM_neg, LAM_pos=LAM_pos, LAM_gr=LAM_gr, LAM_si=LAM_si,
        Q_side=Q_side, Q_lli_lam_trap=Q_lli_lam_trap, age_expansion_um=age_expansion_um,
        rpt_k_arr=rpt_k_arr, rpt_voltage_curves=rpt_voltage_curves,
    )


def _shift_efc(df):
    """Re-anchor the experimental EFC axis to EFC_START_OFFSET (see comment
    at the top of this file) and drop anything at/before the new zero --
    that region (RPT0 and, for capacity, RPT1 itself) carries no information
    once it's treated as the start rather than a real aging data point."""
    df = df.copy()
    df["efc"] = df["efc"] - EFC_START_OFFSET
    return df[df["efc"] >= 0].reset_index(drop=True)


def load_experimental_capacity_fade():
    return _shift_efc(pd.read_csv(os.path.join(EXPERIMENTAL_DATA_DIR, "CELL064_capacity_fade.csv")))


def load_experimental_lam():
    return _shift_efc(pd.read_csv(os.path.join(EXPERIMENTAL_DATA_DIR, "CELL064_LAM_derived.csv")))


def load_experimental_lli():
    return _shift_efc(pd.read_csv(os.path.join(EXPERIMENTAL_DATA_DIR, "CELL064_LLI.csv")))


def load_experimental_reversible_expansion():
    return _shift_efc(pd.read_csv(os.path.join(EXPERIMENTAL_DATA_DIR, "CELL064_reversible_expansion.csv")))


def load_experimental_k_expansion():
    return _shift_efc(pd.read_csv(os.path.join(EXPERIMENTAL_DATA_DIR, "CELL064_k_expansion_scale.csv")))


def score_rpt_gap(results, exp_cap):
    """Per explicit instruction: fit quality from now on should factor in
    the model's OWN C/20 RPT points against the real C/20 RPT points --
    NOT the continuous C/3 ageing curve (comparing that to real RPT data
    conflates the fit with a C/3-vs-C/20 rate/polarisation difference, see
    TUNING STATUS item 12's correction). Interpolates the model's RPT SoH
    trajectory onto each real RPT's EFC and reports the gap in percentage
    points (model - real; positive = model retains more capacity than the
    real cell at that point in life, negative = model is now the one aged
    ahead of the real cell). Printed on every run so this number is always
    visible when judging a fit, not just eyeballed off the plot.
    """
    sim_efc_rpt = efc_from_throughput(results["rpt_thr"]) if results["rpt_thr"].size else np.array([])
    c20 = exp_cap[exp_cap["source"] == "C/20 RPT"]
    if sim_efc_rpt.size < 2 or c20.empty:
        print("RPT-vs-RPT gap: unavailable (not enough model/real RPT points)", flush=True)
        return None
    order = np.argsort(sim_efc_rpt)
    sim_efc_sorted = sim_efc_rpt[order]
    sim_soh_sorted = 100 * results["rpt_cap_arr"][order] / results["age_cap"][0]
    # Only score real RPTs the model's own RPT trajectory actually spans --
    # extrapolating past the model's last RPT would silently invent a number.
    in_range = (c20["efc"] >= sim_efc_sorted[0]) & (c20["efc"] <= sim_efc_sorted[-1])
    print("RPT-vs-RPT gap (model C/20 RPT interpolated at each real C/20 RPT's EFC, "
          "model - real, +ve = model retains more):", flush=True)
    gaps = []
    for _, row in c20[in_range].iterrows():
        model_soh = np.interp(row["efc"], sim_efc_sorted, sim_soh_sorted)
        gap = model_soh - row["capacity_retention_pct"]
        gaps.append(gap)
        print(f"  EFC={row['efc']:.1f}: real={row['capacity_retention_pct']:.1f}%, "
              f"model={model_soh:.1f}%, gap={gap:+.1f}pp", flush=True)
    if not gaps:
        print("  (no real RPT falls within the model's own RPT range)", flush=True)
        return None
    mean_abs_gap = float(np.mean(np.abs(gaps)))
    print(f"  mean |gap| = {mean_abs_gap:.1f}pp across {len(gaps)} in-range real RPT(s)", flush=True)
    return mean_abs_gap


def score_expansion_shape(results, exp_exp, sim_efc_age, n_ref=5):
    """Stage 3 fit-quality metric (same standing-metric convention as
    score_rpt_gap): NORMALISED-shape RMSE between model and real reversible
    expansion (see _plot_reversible_expansion's normalisation rationale --
    absolute magnitude isn't comparable until/unless a proper scale factor
    is found, so this scores shape/trend only). Interpolates the model's
    own normalised trajectory onto each real (normalised) point within the
    model's own EFC range and reports the RMSE; lower is better."""
    valid = ~np.isnan(results["age_expansion_um"])
    model_efc = sim_efc_age[valid]
    model_exp = results["age_expansion_um"][valid]
    if model_exp.size < n_ref or exp_exp.empty:
        print("Expansion shape RMSE: unavailable (not enough model/real points)", flush=True)
        return None
    model_ref = np.median(model_exp[:n_ref])
    exp_vals = exp_exp["reversible_expansion_um"].to_numpy()
    exp_efc = exp_exp["efc"].to_numpy()
    exp_ref = np.median(exp_vals[:n_ref])
    model_norm = model_exp / model_ref
    exp_norm = exp_vals / exp_ref

    order = np.argsort(model_efc)
    model_efc_sorted = model_efc[order]
    model_norm_sorted = model_norm[order]
    in_range = (exp_efc >= model_efc_sorted[0]) & (exp_efc <= model_efc_sorted[-1])
    if not in_range.any():
        print("Expansion shape RMSE: unavailable (no EFC overlap)", flush=True)
        return None
    model_interp = np.interp(exp_efc[in_range], model_efc_sorted, model_norm_sorted)
    rmse = float(np.sqrt(np.mean((model_interp - exp_norm[in_range]) ** 2)))
    print(f"Expansion shape RMSE (normalised): {rmse:.4f} across {int(in_range.sum())} real points", flush=True)
    return rmse


def score_k_gap(results, exp_k, sim_efc_rpt):
    """Stage 3 fit-quality metric, same convention as score_rpt_gap but for
    the expansion scale k (already dimensionless/comparable as-is, no
    normalisation needed): interpolates the model's own RPT k trajectory
    onto each real RPT's EFC and reports the gap (model - real)."""
    if sim_efc_rpt.size < 2 or exp_k.empty:
        print("k-vs-k gap: unavailable (not enough model/real RPT points)", flush=True)
        return None
    order = np.argsort(sim_efc_rpt)
    sim_efc_sorted = sim_efc_rpt[order]
    sim_k_sorted = results["rpt_k_arr"][order]
    in_range = (exp_k["efc"] >= sim_efc_sorted[0]) & (exp_k["efc"] <= sim_efc_sorted[-1])
    print("k-vs-k gap (model interpolated at each real RPT's EFC, model - real):", flush=True)
    gaps = []
    for _, row in exp_k[in_range].iterrows():
        model_k = np.interp(row["efc"], sim_efc_sorted, sim_k_sorted)
        gap = model_k - row["k"]
        gaps.append(gap)
        print(f"  EFC={row['efc']:.1f}: real={row['k']:.3f}, model={model_k:.3f}, gap={gap:+.3f}", flush=True)
    if not gaps:
        print("  (no real RPT falls within the model's own RPT range)", flush=True)
        return None
    mean_abs_gap = float(np.mean(np.abs(gaps)))
    print(f"  mean |gap| = {mean_abs_gap:.3f} across {len(gaps)} in-range real RPT(s)", flush=True)
    return mean_abs_gap


# Real knee anchors (EFC_START_OFFSET-shifted frame): the expansion peak
# lands at absolute EFC ~152 -> ~101 shifted, and real RPT3 = 91.4% SoH at
# EFC ~98 shifted. score_knee below is judged against these.
KNEE_EFC_ANCHOR = 152.0 - EFC_START_OFFSET   # ~101
KNEE_ONSET_SOH_ANCHOR = 91.0


def score_knee(results, sim_efc_age):
    """Stage 3 (item 19+) fit-quality metric: WHERE is the model's knee?
    Phase-1 sweeps showed combos that improve the RPT gap by softening the
    post-knee slope while leaving the knee itself late -- the other score_*
    metrics don't catch that.

    Definition (robust, works on truncated runs): the model's knee EFC is
    the STEEPEST-DESCENT point of its continuous C/3 SoH curve (the sharp
    collapse / bend). Judged against KNEE_EFC_ANCHOR (~101, the real
    expansion-peak EFC). Also report the SoH-91% crossing EFC and the SoH
    at the knee for context. Composite term = |knee EFC - KNEE_EFC_ANCHOR|
    (EFC units) -- lower = better. (The bare 91% crossing is NOT used as the
    metric: on a soft/gentle fade it's reached during the pre-knee slope
    and understates how late the actual collapse is -- exactly the phase-1
    failure mode this metric exists to catch.)"""
    soh = 100.0 * results["age_cap"] / results["age_cap"][0]
    efc = np.asarray(sim_efc_age, float)
    if soh.size < 5 or efc.size != soh.size:
        print("Knee: unavailable (not enough SoH points)", flush=True)
        return None
    d = np.gradient(soh, efc)
    if d.size >= 5:
        d = np.convolve(d, np.ones(3) / 3.0, mode="same")
    knee_idx = int(np.argmin(d))
    knee_efc = float(efc[knee_idx])
    soh_at_knee = float(soh[knee_idx])
    cross91 = None
    for i in range(1, soh.size):
        if soh[i - 1] >= KNEE_ONSET_SOH_ANCHOR > soh[i]:
            f = (soh[i - 1] - KNEE_ONSET_SOH_ANCHOR) / (soh[i - 1] - soh[i])
            cross91 = float(efc[i - 1] + f * (efc[i] - efc[i - 1]))
            break
    efc_err = knee_efc - KNEE_EFC_ANCHOR
    cross_str = "n/a" if cross91 is None else f"{cross91:.1f}"
    print(f"Knee: knee EFC={knee_efc:.1f} (err {efc_err:+.1f}), "
          f"SoH at knee={soh_at_knee:.1f}%, SoH-91% crossing EFC={cross_str} "
          f"(anchor EFC~{KNEE_EFC_ANCHOR:.0f})", flush=True)
    return abs(efc_err)


# Real RPT discharge curves (Q vs. V) live in the separate Pouch_Data repo,
# not ported into experimental_data/ -- see README's "Inputs still needed".
# Only RPT1/2/4/5 have a real C/20 discharge (RPT0/3's failed QC).
REAL_RPT_DISCHARGE_DIR = ("C:/Shannan_PhD_Stuff_Local/Si_Gr_Expansion_Precursor/"
                          "Pouch_Data/cell064_data/data_timeseries/low_rate_c20")
REAL_RPT_NUMS_WITH_DISCHARGE = (1, 2, 4, 5)


def load_real_rpt_discharge(rpt_num):
    path = os.path.join(REAL_RPT_DISCHARGE_DIR,
                        f"CELL064_RPT{rpt_num:03d}_lowrate_c20_discharge.csv")
    df = pd.read_csv(path, usecols=["Q_discharge_Ah", "voltage_terminal_V"])
    q = df["Q_discharge_Ah"].to_numpy()
    v = df["voltage_terminal_V"].to_numpy()
    return q - q[0], v


def score_voltage_shape(results, exp_cap):
    """Task 2 (test_pathways) fit-quality metric: how well does the MODEL's
    own C/20 RPT discharge-curve SHAPE (not just its SoH-at-RPT single
    number) track the real RPT discharge curves across degradation? Task 1
    (sweeps/test_pathways/task1_bol_pathway_fit.py) found the model's
    simulated pathway over-weights LLI vs. Si-LAM relative to what the real
    cell's own discharge curves imply -- this metric lets a sweep optimise
    directly against the voltage-curve shape rather than only SoH/expansion/
    knee. Matches each real RPT (with a real discharge curve -- RPT1/2/4/5)
    to the model's nearest-EFC own RPT curve, computes the V(Q) RMSE over
    the overlapping capacity range for each, and returns the mean across
    all matched pairs (lower = better)."""
    model_curves = results["rpt_voltage_curves"]
    if not model_curves:
        print("Voltage shape RMSE: unavailable (no model RPT curves)", flush=True)
        return None
    model_efcs = efc_from_throughput(np.array([c["thr"] for c in model_curves]))
    c20 = exp_cap[exp_cap["source"] == "C/20 RPT"]
    real_efc_by_rpt = dict(zip(c20["rpt"], c20["efc"]))
    print("Voltage shape RMSE (model's own RPT curve vs. real, matched by nearest EFC):", flush=True)
    rmses = []
    used_js = set()
    for rpt_num in REAL_RPT_NUMS_WITH_DISCHARGE:
        if rpt_num not in real_efc_by_rpt:
            continue
        real_efc = real_efc_by_rpt[rpt_num]
        order = np.argsort(np.abs(model_efcs - real_efc))
        j = next((int(k) for k in order if int(k) not in used_js), int(order[0]))
        used_js.add(j)
        model_efc = float(model_efcs[j])
        if abs(model_efc - real_efc) > 20.0:
            print(f"  RPT{rpt_num} (EFC {real_efc:.1f}): SKIPPED, nearest model RPT is "
                  f"{abs(model_efc - real_efc):.1f} EFC away", flush=True)
            continue
        q_real, v_real = load_real_rpt_discharge(rpt_num)
        q_model, v_model = model_curves[j]["q"], model_curves[j]["v"]
        in_range = q_real <= q_model.max()
        if in_range.sum() < 5:
            continue
        v_model_interp = np.interp(q_real[in_range], q_model, v_model)
        rmse = float(np.sqrt(np.mean((v_model_interp - v_real[in_range]) ** 2)))
        rmses.append(rmse)
        print(f"  RPT{rpt_num} (EFC {real_efc:.1f}, model EFC {model_efc:.1f}): RMSE={rmse:.4f} V",
              flush=True)
    if not rmses:
        print("  (no real RPT discharge curve fell within the model's own RPT range)", flush=True)
        return None
    mean_rmse = float(np.mean(rmses))
    print(f"  mean voltage-shape RMSE = {mean_rmse:.4f} V across {len(rmses)} matched RPT(s)",
          flush=True)
    return mean_rmse


def _plot_soh(ax, results, exp_cap, sim_efc_age, sim_efc_rpt):
    ax.plot(sim_efc_age, 100 * results["age_cap"] / results["age_cap"][0],
            "o-", ms=3, color="tab:blue", label="Model (C/3 ageing)")
    if sim_efc_rpt.size:
        ax.scatter(sim_efc_rpt, 100 * results["rpt_cap_arr"] / results["age_cap"][0],
                   marker="^", color="tab:orange", s=40, label=f"Model RPT ({RPT_RATE})", zorder=5)
    c20 = exp_cap[exp_cap["source"] == "C/20 RPT"]
    hrsub = exp_cap[exp_cap["source"] != "C/20 RPT"]
    ax.scatter(c20["efc"], c20["capacity_retention_pct"], s=70, color="black",
               label="Experimental SoH (C/20 RPT)", zorder=6)
    ax.scatter(hrsub["efc"], hrsub["capacity_retention_pct"], s=70, facecolors="none",
               edgecolors="black", linewidths=1.4, label="Experimental (high-rate substitute)", zorder=6)
    knee_efc_shifted = 152.0 - EFC_START_OFFSET
    ax.axvline(knee_efc_shifted, color="gray", ls=":",
               label=f"Experimental knee EFC (~{knee_efc_shifted:.0f} since RPT1)")
    ax.set_xlabel("EFC since RPT1 [-]")
    ax.set_ylabel("SoH [%]")
    ax.set_title("SoH: model vs. experimental")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)


def _plot_lam(ax, results, exp_lam, sim_efc_full, combined_with_lli=False, exp_lli=None):
    ax.plot(sim_efc_full, results["LAM_gr"], color="saddlebrown", label="Model LAM Gr (primary)")
    ax.plot(sim_efc_full, results["LAM_si"], color="darkorange", label="Model LAM Si (secondary)")
    ax.plot(sim_efc_full, results["LAM_pos"], color="seagreen", label="Model LAM positive")
    ax.scatter(exp_lam["efc"], exp_lam["LAM_NE_graphite_pct"], color="saddlebrown", marker="s",
               label="Exp. LAM_NE graphite")
    ax.scatter(exp_lam["efc"], exp_lam["LAM_NE_silicon_pct"], color="darkorange", marker="s",
               label="Exp. LAM_NE silicon")
    ax.scatter(exp_lam["efc"], exp_lam["LAM_PE_pct"], color="seagreen", marker="s", label="Exp. LAM_PE")
    if combined_with_lli:
        LLI_pct = 100 * results["Q_side"] / NOMINAL_CAP_AH
        LLI_corrected_pct = 100 * (results["Q_side"] + results["Q_lli_lam_trap"]) / NOMINAL_CAP_AH
        ax.plot(sim_efc_full, LLI_pct, color="crimson", ls="--", label="Model LLI (side reactions only)")
        ax.plot(sim_efc_full, LLI_corrected_pct, color="crimson", ls="-",
                label="Model LLI (side reactions + LAM-trapped Li)")
        exp_lli_pct = 100 * exp_lli["charge_LLI_loss"] / exp_lli["charge_LLI"].iloc[0]
        ax.scatter(exp_lli["efc"], exp_lli_pct, color="crimson", marker="D", label="Exp. LLI")
        ax.set_ylabel("LAM / LLI [%]")
        ax.set_title("LAM by phase + LLI: model vs. experimental")
    else:
        ax.set_ylabel("LAM [%]")
        ax.set_title("LAM by phase: model vs. experimental (rough, non-DMA) estimate")
    ax.set_xlabel("EFC [-]")
    ax.legend(fontsize=6)
    ax.grid(alpha=0.3)


def _plot_lli(ax, results, exp_lli, sim_efc_full):
    LLI_pct = 100 * results["Q_side"] / NOMINAL_CAP_AH
    LLI_corrected_pct = 100 * (results["Q_side"] + results["Q_lli_lam_trap"]) / NOMINAL_CAP_AH
    ax.plot(sim_efc_full, LLI_pct, color="crimson", ls="--", label="Model LLI (side reactions only)")
    ax.plot(sim_efc_full, LLI_corrected_pct, color="crimson", ls="-",
            label="Model LLI (side reactions + LAM-trapped Li)")
    exp_lli_pct = 100 * exp_lli["charge_LLI_loss"] / exp_lli["charge_LLI"].iloc[0]
    ax.scatter(exp_lli["efc"], exp_lli_pct, color="crimson", marker="s",
               label="Exp. charge_LLI_loss (source-package fit)")
    ax.set_xlabel("EFC [-]")
    ax.set_ylabel("LLI [%]")
    ax.set_title("LLI: model vs. experimental (side-rxn-only vs. +LAM-trapped)")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)


def _plot_rpt_voltage(ax, results, exp_cap=None):
    """Model's own C/20 RPT discharge curves, coloured by EFC -- plus, when
    exp_cap is given, the real experimental RPT1/2/4/5 discharge curves
    overlaid (dashed, same colour scale, matched to the nearest model RPT by
    EFC) -- folds cell064_rpt_discharge_overlay.py's comparison into the
    standard plot_all() output so it doesn't need a separate re-run to see.
    Matching/skip-if-too-far logic mirrors score_voltage_shape()'s."""
    curves = results["rpt_voltage_curves"]
    if not curves:
        ax.set_title("Discharge voltage vs. RPT (unavailable)")
        return
    efcs = efc_from_throughput(np.array([c["thr"] for c in curves]))
    lo, hi = float(efcs.min()), float(efcs.max())
    norm = plt.Normalize(lo, hi) if hi > lo else plt.Normalize(0, 1)
    cmap = plt.cm.viridis
    for c, efc in zip(curves, efcs):
        ax.plot(c["q"], c["v"], color=cmap(norm(efc)), label=f"model RPT @ EFC~{efc:.0f}")

    if exp_cap is not None:
        c20 = exp_cap[exp_cap["source"] == "C/20 RPT"]
        real_efc_by_rpt = dict(zip(c20["rpt"], c20["efc"]))
        used_js = set()
        for rpt_num in REAL_RPT_NUMS_WITH_DISCHARGE:
            if rpt_num not in real_efc_by_rpt:
                continue
            real_efc = real_efc_by_rpt[rpt_num]
            order = np.argsort(np.abs(efcs - real_efc))
            j = next((int(k) for k in order if int(k) not in used_js), int(order[0]))
            mismatch = abs(float(efcs[j]) - real_efc)
            # Looser than score_voltage_shape's 20.0 EFC cutoff (that one
            # gates a NUMERIC composite score, where a badly-mismatched RPT
            # would be misleading) -- this is a visual overlay, so show it
            # even at a modest mismatch, but label the gap so it stays
            # honest about how close the comparison really is.
            if mismatch > 60.0:
                continue
            used_js.add(j)
            try:
                q_r, v_r = load_real_rpt_discharge(rpt_num)
            except FileNotFoundError:
                continue
            gap_note = f" (Δ{mismatch:.0f} EFC)" if mismatch > 20.0 else ""
            ax.plot(q_r, v_r, color=cmap(norm(real_efc)), ls="--", lw=1.3,
                    label=f"real RPT{rpt_num} @ EFC~{real_efc:.0f}{gap_note}")

    ax.set_xlabel(f"Discharge capacity [A.h] ({RPT_RATE})")
    ax.set_ylabel("Terminal voltage [V]")
    title = "Discharge voltage vs. RPT (solid=model, dashed=real)" if exp_cap is not None \
        else "Discharge voltage vs. RPT (model only, see README)"
    ax.set_title(title)
    ax.legend(fontsize=6.5)
    ax.grid(alpha=0.3)


def _plot_reversible_expansion(ax, results, exp_exp, sim_efc_age, n_ref=5):
    """Per instruction: absolute magnitude isn't comparable yet -- the
    model's single-representative-electrode-pair thickness change and the
    real cell's whole-stack dilatometry reading differ by an uncalibrated
    scale factor ("Number of electrodes connected in parallel" was ruled
    out as that scale knob since it also scales capacity/current via A_cc
    -- see the comment at this file's stage-3 section). Normalise each
    series by its own early-life (first n_ref-point median) value instead,
    so the comparison is of SHAPE/relative trend, not absolute magnitude."""
    valid = ~np.isnan(results["age_expansion_um"])
    model_efc = sim_efc_age[valid]
    model_exp = results["age_expansion_um"][valid]
    model_ref = np.median(model_exp[:n_ref]) if model_exp.size else np.nan
    exp_vals = exp_exp["reversible_expansion_um"].to_numpy()
    exp_ref = np.median(exp_vals[:n_ref]) if exp_vals.size else np.nan

    ax.plot(model_efc, model_exp / model_ref, "-", lw=1, color="tab:purple",
            label="Model (C/3 ageing, normalised to early-life value)")
    ax.scatter(exp_exp["efc"], exp_vals / exp_ref, s=8, alpha=0.4, color="black",
               label="Experimental (normalised to early-life value)")
    ax.set_xlabel("EFC [-]")
    ax.set_ylabel("Reversible expansion, normalised [-]")
    ax.set_title("Reversible expansion (normalised): model vs. experimental")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)


def _plot_k_ratio(ax, results, exp_k, sim_efc_rpt):
    if sim_efc_rpt.size:
        # Hollow, larger-than-the-experimental-dot, and drawn ON TOP (zorder
        # above the black dots below) -- at BOL the model and real k often
        # coincide almost exactly, and a solid marker would otherwise fully
        # hide behind (or hide) the other point at that spot. The open ring
        # stays visible around the dot even at perfect overlap.
        ax.scatter(sim_efc_rpt, results["rpt_k_arr"], marker="^", s=110,
                   facecolors="none", edgecolors="tab:orange", linewidths=1.8,
                   label=f"Model RPT ({RPT_RATE})", zorder=7)
    c20k = exp_k[exp_k["source"] == "C/20 RPT"]
    hrk = exp_k[exp_k["source"] != "C/20 RPT"]
    ax.scatter(c20k["efc"], c20k["k"], s=70, color="black", label="Experimental C/20 RPT", zorder=6)
    ax.scatter(hrk["efc"], hrk["k"], s=70, facecolors="none", edgecolors="black", linewidths=1.4,
               label="Experimental (high-rate substitute)", zorder=6)
    ax.set_xlabel("EFC [-]")
    ax.set_ylabel("k = delta_rev,cell / delta_rev,particle [-]")
    ax.set_title("Expansion scale k: model vs. experimental")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)


def _savefig(fig, base_name):
    out_tag = os.environ.get("C64_OUT_TAG", "")  # set to avoid clobbering a probe's plot, e.g. floor comparisons
    out_name = f"{base_name}{('_' + out_tag) if out_tag else ''}.png"
    outpath = os.path.join(SCRIPT_DIR, out_name)
    fig.savefig(outpath, dpi=150)
    print(f"Saved: {outpath}")


def plot_degradation_diagnostics(results, exp_cap, exp_lam, exp_lli, sim_efc_age, sim_efc_rpt, sim_efc_full):
    """Figure 1: SoH, LAM, LLI, RPT discharge voltage (each its own panel)."""
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    _plot_soh(ax[0, 0], results, exp_cap, sim_efc_age, sim_efc_rpt)
    _plot_lam(ax[0, 1], results, exp_lam, sim_efc_full)
    _plot_lli(ax[1, 0], results, exp_lli, sim_efc_full)
    _plot_rpt_voltage(ax[1, 1], results, exp_cap)
    fig.suptitle("CELL064: degradation diagnostics (SoH, LAM, LLI, RPT voltage)")
    plt.tight_layout()
    _savefig(fig, "cell064_degradation_fit_result")


def plot_expansion_diagnostics(results, exp_exp, exp_k, sim_efc_age, sim_efc_rpt):
    """Figure 2: reversible expansion, expansion scale k (each its own panel)."""
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    _plot_reversible_expansion(ax[0], results, exp_exp, sim_efc_age)
    _plot_k_ratio(ax[1], results, exp_k, sim_efc_rpt)
    fig.suptitle("CELL064: expansion diagnostics (pore-buffering params not yet retuned -- stage 3 starting point)")
    plt.tight_layout()
    _savefig(fig, "cell064_expansion_fit_result")


def plot_overall_summary(results, exp_cap, exp_lam, exp_lli, exp_exp, exp_k,
                          sim_efc_age, sim_efc_rpt, sim_efc_full):
    """Figure 3: everything above in one figure, LAM+LLI combined into one panel."""
    fig, ax = plt.subplots(2, 3, figsize=(18, 9))
    _plot_soh(ax[0, 0], results, exp_cap, sim_efc_age, sim_efc_rpt)
    _plot_lam(ax[0, 1], results, exp_lam, sim_efc_full, combined_with_lli=True, exp_lli=exp_lli)
    _plot_rpt_voltage(ax[0, 2], results, exp_cap)
    _plot_reversible_expansion(ax[1, 0], results, exp_exp, sim_efc_age)
    _plot_k_ratio(ax[1, 1], results, exp_k, sim_efc_rpt)
    ax[1, 2].axis("off")
    fig.suptitle("CELL064: overall summary (degradation + expansion)")
    plt.tight_layout()
    _savefig(fig, "cell064_overall_summary_result")


def _cycle_thermal_resistance(cyc):
    """Return (efc, T_mean [K], R_mean [Ohm], I_mean [A]) for cyc's main
    ageing discharge step (mean_I > 0.5, same step-selection rule as
    cycle_ageing_leg), or None if unavailable. Under "isothermal" (default),
    T_mean stays pinned at T_init for every cycle -- a flat line is the
    expected/correct result, not a bug; plot_thermal_diagnostics overlays a
    quasi-steady ESTIMATED temperature (from I, R, and literature h/A_cool)
    alongside it for a genuinely informative comparison."""
    for step in getattr(cyc, "steps", []):
        try:
            I = step["Current [A]"].entries
            T = step["Volume-averaged cell temperature [K]"].entries
            R = step["Local ECM resistance [Ohm]"].entries
            thr_step = step["Throughput capacity [A.h]"].entries
        except (KeyError, TypeError, AttributeError):
            continue
        if I.size < 2:
            continue
        if np.mean(I) > 0.5:
            efc = float(efc_from_throughput(np.array([float(thr_step[-1])]))[0])
            return efc, float(np.mean(T)), float(np.mean(R)), float(np.mean(np.abs(I)))
    return None


def plot_thermal_diagnostics(sol):
    """Item 26: cell temperature and internal (ECM) resistance across life,
    on the main ageing discharge leg -- a SEPARATE figure from the existing
    degradation/expansion/summary/stoichiometry ones. Purpose: check whether
    "lumped" thermal (C64_THERMAL_OPTION=lumped) produces a plausible
    temperature rise post-knee (from growing internal resistance), which is
    the physical mechanism behind base_mechanics.py's Ai2019-eq-13 thermal-
    expansion term already added to "Cell thickness change [m]" -- see
    THERMAL_OPTION's doc-comment for the caveat about si_gr_expansion.py's
    OTHER (hardcoded, non-overridable) Arrhenius temperature dependence this
    would also reactivate."""
    efcs, Ts, Rs, Is = [], [], [], []
    for cyc in sol.cycles:
        out = _cycle_thermal_resistance(cyc)
        if out is None:
            continue
        efc, T, R, I = out
        efcs.append(efc)
        Ts.append(T)
        Rs.append(R)
        Is.append(I)
    if not efcs:
        print("[THERMAL] no ageing discharge legs found -- skipping thermal diagnostics plot",
              flush=True)
        return
    efcs = np.array(efcs)
    Ts = np.array(Ts)
    Rs = np.array(Rs)
    Is = np.array(Is)

    # Quasi-steady ESTIMATED temperature rise from the model's own I, R (same
    # h/A_cool as used elsewhere -- flagged non-CELL064-specific, so this is
    # illustrative of SHAPE/relative magnitude, not a validated absolute
    # value) -- only informative overlay when the true isothermal T is flat.
    H_TOTAL, A_COOL = 10.0, 0.00170857
    T_est_C = (Is ** 2 * Rs / (H_TOTAL * A_COOL))  # degC rise above ambient

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    ax[0].plot(efcs, Ts - 273.15, "o-", color="tab:red", ms=3, label="Model T (isothermal, flat by construction)")
    if THERMAL_OPTION == "isothermal":
        ax2 = ax[0].twinx()
        ax2.plot(efcs, T_est_C, "s--", color="tab:orange", ms=3,
                  label="Quasi-steady ESTIMATED dT (from I,R; h/A_cool uncertain)")
        ax2.set_ylabel("Estimated temperature RISE above ambient [K]", color="tab:orange")
        ax2.tick_params(axis="y", labelcolor="tab:orange")
        lines1, labels1 = ax[0].get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax[0].legend(lines1 + lines2, labels1 + labels2, fontsize=7)
    else:
        ax[0].legend(fontsize=8)
    ax[0].set_xlabel("EFC [-]")
    ax[0].set_ylabel("Volume-averaged cell temperature [degC]")
    ax[0].set_title("Cell temperature across life")
    ax[0].grid(alpha=0.3)

    ax[1].plot(efcs, Rs * 1000, "o-", color="tab:blue", ms=3)
    ax[1].set_xlabel("EFC [-]")
    ax[1].set_ylabel("Local ECM resistance [mOhm]")
    ax[1].set_title("Internal resistance across life")
    ax[1].grid(alpha=0.3)

    fig.suptitle(f"CELL064: thermal diagnostics (item 26, thermal option = {THERMAL_OPTION!r})")
    plt.tight_layout()
    _savefig(fig, "cell064_thermal_diagnostics_result")
    print(f"[THERMAL] T range: {Ts.min() - 273.15:.2f} to {Ts.max() - 273.15:.2f} degC "
          f"(delta={Ts.max() - Ts.min():.2f} K); R range: {Rs.min() * 1000:.2f} to "
          f"{Rs.max() * 1000:.2f} mOhm", flush=True)


def _rpt_stoich_curve(cyc):
    """Return (q_rel [A.h], v [V], x_gr, x_si) for cyc's RPT discharge step,
    or None if cyc has no such step. Same step-selection rule as
    rpt_discharge_curve(). Folded in from the former standalone
    plot_stoichiometry_evolution.py so it reuses the SAME `sol` already
    computed for the other figures instead of re-running the simulation."""
    best_cap = 0.0
    best = None
    for step in getattr(cyc, "steps", []):
        try:
            I = step["Current [A]"].entries
            q = step["Discharge capacity [A.h]"].entries
            v = step["Terminal voltage [V]"].entries
        except (KeyError, TypeError, AttributeError):
            continue
        if I.size < 2:
            continue
        mean_I = np.mean(I)
        if 0.05 < mean_I <= 0.5 and q.size >= 2:
            c = float(q[-1] - q[0])
            if c > best_cap:
                try:
                    x_gr = step["X-averaged negative primary particle surface "
                               "stoichiometry"].entries
                    x_si = step["X-averaged negative secondary particle surface "
                               "stoichiometry"].entries
                except (KeyError, TypeError, AttributeError):
                    continue
                best_cap = c
                best = (q - q[0], v, x_gr, x_si)
    return best


def plot_stoichiometry_evolution(sol):
    """Gr's and Si's own x-averaged surface stoichiometry vs. discharge
    capacity, across EVERY RPT the model reaches -- shows whether/when each
    electrode's window collapses, and by how much, as life progresses.
    Diagnostic test of whether a config retrieves "Gr window still open
    longer, Si genuinely shrinks/empties first" (real-cell behavior) vs.
    "Gr collapses abruptly, Si carries the tail" (the model's old wrong-
    mechanism behavior). Takes the already-computed `sol` from
    run_degradation() -- no separate simulation run needed."""
    curves = []
    for cyc in sol.cycles:
        out = _rpt_stoich_curve(cyc)
        if out is None:
            continue
        q, v, x_gr, x_si = out
        try:
            thr_end = float(cyc["Throughput capacity [A.h]"].entries[-1])
        except (KeyError, TypeError, AttributeError, IndexError):
            continue
        efc = efc_from_throughput(np.array([thr_end]))[0]
        curves.append(dict(efc=efc, q=q, v=v, x_gr=x_gr, x_si=x_si))

    efc_list = [round(c["efc"], 1) for c in curves]
    print(f"[STOICH] found {len(curves)} RPT-classified discharge curves at EFC: "
          f"{efc_list}", flush=True)
    if not curves:
        print("[STOICH] no RPT curves found -- nothing to plot.", flush=True)
        return

    cmap = plt.get_cmap("viridis")
    n = len(curves)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for i, c in enumerate(curves):
        color = cmap(i / max(n - 1, 1))
        label = f"EFC={c['efc']:.0f}"
        axes[0].plot(c["q"], c["x_gr"], color=color, label=label)
        axes[1].plot(c["q"], c["x_si"], color=color, label=label)
        axes[2].plot(c["q"], c["v"], color=color, label=label)
    axes[0].set_title("Graphite (primary) surface stoichiometry")
    axes[1].set_title("Silicon (secondary) surface stoichiometry")
    axes[2].set_title("Terminal voltage")
    for ax in axes:
        ax.set_xlabel("Discharge capacity Q [A.h]")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("x (stoichiometry) [-]")
    axes[2].set_ylabel("V [V]")
    axes[0].legend(fontsize=7, ncol=2)
    fig.suptitle("Gr/Si stoichiometry-window evolution across life (all RPTs)")
    plt.tight_layout()
    _savefig(fig, "stoichiometry_evolution")


def plot_all(results, sol=None):
    exp_cap = load_experimental_capacity_fade()
    exp_lam = load_experimental_lam()
    exp_lli = load_experimental_lli()
    exp_exp = load_experimental_reversible_expansion()
    exp_k = load_experimental_k_expansion()

    sim_efc_age = efc_from_throughput(results["age_thr"])
    sim_efc_rpt = efc_from_throughput(results["rpt_thr"]) if results["rpt_thr"].size else np.array([])
    sim_efc_full = efc_from_throughput(results["Qt"])

    score_rpt_gap(results, exp_cap)
    score_expansion_shape(results, exp_exp, sim_efc_age)
    score_k_gap(results, exp_k, sim_efc_rpt)
    score_knee(results, sim_efc_age)
    score_voltage_shape(results, exp_cap)

    if os.environ.get("C64_NO_PLOT") == "1":
        # Sweep mode (sweeps/cell064_joint_sweep.py): the driver only
        # needs the three score_* lines above; skip the figures so a large
        # grid doesn't leave dozens of tagged PNGs behind.
        print("C64_NO_PLOT=1 -- skipping figure generation", flush=True)
        return

    plot_degradation_diagnostics(results, exp_cap, exp_lam, exp_lli, sim_efc_age, sim_efc_rpt, sim_efc_full)
    plot_expansion_diagnostics(results, exp_exp, exp_k, sim_efc_age, sim_efc_rpt)
    plot_overall_summary(results, exp_cap, exp_lam, exp_lli, exp_exp, exp_k, sim_efc_age, sim_efc_rpt, sim_efc_full)
    if sol is not None:
        plot_stoichiometry_evolution(sol)
        plot_thermal_diagnostics(sol)


if __name__ == "__main__":
    sol = run_degradation()
    results = extract_results(sol)
    plot_all(results, sol)
