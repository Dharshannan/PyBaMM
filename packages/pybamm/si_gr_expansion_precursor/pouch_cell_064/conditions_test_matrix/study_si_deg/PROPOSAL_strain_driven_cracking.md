# Proposal: strain-driven (volume-change fatigue) Si crack growth

**Status:** proposal only. Nothing is implemented. It is intended for a
separate branch, and it is a new mechanism, so it needs sign-off before any
code is written. Drafted 2026-10-01; revised the same day with variant B2,
the knee calibration and the three-condition refit plan.

**Notation.**
- **Route A:** the existing Paris-law cracking on the diffusion-induced
  (concentration-gradient) surface stress.
- **Route B:** the proposed cracking driven by Si volumetric strain (fatigue).

## 1. Motivation

There are four cells at 103 kPa with C/3 cycling. The table gives each
cell's knee and its cumulative reversible cell-expansion swing up to that
knee. Swings are measured per ageing cycle from the lifetime expansion data
(`extract_lifetime_reversible_expansion.py`).

| Cell | Window, T | Knee (EFC since RPT1) | Cycles to knee | Pre-knee swing per cycle | Cumulative swing to knee |
|---|---|---|---|---|---|
| CELL064 | 0–100%, 25 °C | ~101 | ~104 | 54.6 µm | **5,700 µm** |
| CELL017 | 0–100%, 45 °C | ~310 | ~316 | 51.4 µm | **16,300 µm** |
| CELL009 | 15–95%, 25 °C | ~290 | ~377 | 37.1 µm | **14,000 µm** |
| CELL026 | 20–80%, 25 °C | ~293 | ~488 | 31.7 µm | **15,500 µm** |

**Three of the four cells knee at about the same cumulative swing (and EFC),**
across two temperatures and three SoC windows. Only CELL064 (full window,
25 °C) knees about 3× earlier.

The current model has route A only. To reproduce the other cells it needs
two condition-specific corrections:

- **45 °C (CELL017):** `SI_CRACK_EAC` = 90 kJ/mol, an Arrhenius slowdown of
  cracking at high T with no strong physical basis. Faster diffusion at
  45 °C already lowers the pre-knee stress from ~43 MPa to ~16–20 MPa.
- **Partial window (CELL009):** `SI_CRACK_MULT` ×11 (1.1e4 against the shared
  1e3). During the C/3 partial cycles the model's peak Si surface stress is
  ~2 MPa, against ~35 MPa during the full-depth C/20 RPTs. The model's
  partial-window cells crack almost only during their RPTs. Plain v3 does
  knee on CELL009, but at EFC ~500–550 (FINDINGS.md), through slow pore
  clogging rather than crack growth.

**Hypothesis: two independent damage routes.**

- **A. Diffusion-stress fracture (existing).** Significant only when the
  gradient stress is high: full-window cycling at 25 °C, with deep Si
  delithiation and slow diffusion. It adds on top of B for CELL064.
- **B. Volume-change fatigue (new).** Damage grows with the cumulative Si
  volume strain traversed, i.e. how much the Si breathes, with weak or no
  temperature dependence. It sets the ~300 EFC knee common to all cells.

The literature (section 9) supports the idea. Karger, O'Kane et al. (2024)
attribute their zero volume-change exponent to constant-depth data, and say
their model is valid only for full-SoC cycling. Philipp et al. (2026) find
route A alone gives negligible cracking at intermediate SoC.

## 2. Existing law (route A)

From `models/submodels/particle_mechanics/crack_propagation.py`, with the
stress-relief cap (CHANGES.md item 2):

$$
\frac{dl}{dt}\bigg|_{A} = \frac{k_{\mathrm{cr}}(T)}{3600}
\Big(\sigma_{t,\mathrm{eff}}\, b_{\mathrm{cr}} \sqrt{\pi l}\Big)^{m_{\mathrm{cr}}},
\qquad
\sigma_{t,\mathrm{eff}} = \sigma_t \max\!\left(1 - \frac{l}{R}, 0\right)\ \ (\ge 0).
\tag{1}
$$

Downstream (unchanged by this proposal):
- roughness $\rho = 1 + 2\, l\, w_{\mathrm{cr}} \rho_{\mathrm{cr}}$;
- SEI on cracks on area $a(\rho - 1)$;
- from that: crack-SEI LLI, pore closure, the porosity-isolation gate, and
  the knee.

## 3. Proposed route B

### 3.1 Si volumetric strain and strain rate

Use the Si phase's own volume law, the same `t_change` used for thickness.
With "volume change aging deformation" on, this is the fitted
$t(\bar s) = 1 + 3\,\bar s^{\,n}$ ($n = n_{\mathrm{eff}}$). It is evaluated at
the r-averaged Si stoichiometry $\bar s$.

$$
\varepsilon_V = \ln t(\bar s),
\qquad
\dot\varepsilon_V = \frac{t'(\bar s)}{t(\bar s)}\,\frac{d\bar s}{dt},
\qquad
t'(\bar s) = 3n\,\bar s^{\,n-1},
\qquad
\frac{d\bar s}{dt} = \frac{1}{c_{\max}}\frac{\partial \bar c_s}{\partial t} = -\frac{3\, j_{\mathrm{Si}}}{F\, R_{\mathrm{Si}}\, c_{\max}}.
\tag{2}
$$

- **No new state.** $j_{\mathrm{Si}}$ is the Si interfacial current density,
  already in the model.
- **Log strain.** It is bounded and treats expansion and contraction
  symmetrically under Si's ~300% swing.

### 3.2 Crack growth with both routes

$$
\frac{dl}{dt} = \underbrace{\frac{k_{\mathrm{cr}}(T)}{3600}\big(\sigma_{t,\mathrm{eff}}\, b_{\mathrm{cr}}\sqrt{\pi l}\big)^{m_{\mathrm{cr}}}}_{A\ \text{(existing)}}
\;+\;
\underbrace{k_V(T)\; D(\dot\varepsilon_V)\; l\left(1 - \frac{l}{R}\right)}_{B\ \text{(new)}}.
\tag{3}
$$

**Strain driver $D$, two variants (one option):**

| Variant | $D(\dot\varepsilon_V)$ | Physical reading |
|---|---|---|
| **B1** | $\lvert\dot\varepsilon_V\rvert$ | damage per unit strain traversed, expansion and contraction alike |
| **B2** | $\max(-\dot\varepsilon_V,\,0)$ | contraction only, i.e. delithiation, when the Si surface is in tension. Same tensile-only logic as route A's $\sigma \ge 0$ |

**B2 is the more physically defensible; B1 is the simplest.** They differ by
a factor of ~2 in $k_V$ for symmetric cycling. They differ in shape when
lithiation and delithiation are asymmetric: CV holds, partial windows that
sit high or low, rests.

- **$k_V$** (dimensionless): fractional crack growth per unit log-strain.
  Default: no temperature dependence. Optionally
  $k_V = k_{V,0}\exp[E_V/R\,(1/T - 1/T_{\mathrm{ref}})]$, to test whether any
  is needed.
- **$l(1 - l/R)$:** self-similar growth that reuses the existing cap, so
  cracks stop at the particle radius.

### 3.3 Calibration form and knee condition

Ignoring route A and the cap, route B integrates exactly:

$$
\ln\frac{l}{l_0} = k_V\, \mathcal{E}(t),
\qquad
\mathcal{E}(t) = \int_0^t D(\dot\varepsilon_V)\, dt'.
\tag{4}
$$

$\mathcal{E}$ is the cumulative Si strain, independent of C-rate. Per cycle:

$$
\Delta\mathcal{E}_{\mathrm{cycle}} = 2\,\Delta\varepsilon_{V,\mathrm{cyc}}\ \ (\mathrm{B1}),
\qquad
\Delta\mathcal{E}_{\mathrm{cycle}} = \Delta\varepsilon_{V,\mathrm{cyc}}\ \ (\mathrm{B2}),
\qquad
\Delta\varepsilon_{V,\mathrm{cyc}} = \ln\frac{t(\bar s_{\max})}{t(\bar s_{\min})}.
$$

**Knee condition.** In the CELL009 runs with crack growth (P1c,
P3a_ref_w1e-3), the knee came at a crack growth of about
$(l/l_0)_{\mathrm{knee}} \approx 2$–2.4. So:

$$
k_V \approx \frac{\ln (l/l_0)_{\mathrm{knee}}}{\mathcal{E}_{\mathrm{knee}}} \approx \frac{0.8}{N_{\mathrm{knee}}\, \Delta\mathcal{E}_{\mathrm{cycle}}},
\tag{5}
$$

where $N_{\mathrm{knee}}$ is the number of ageing cycles to the knee.

**Order of magnitude.** A full-window Si swing ($\bar s$: 0 → 1, $t$: 1 → 4)
gives $\Delta\varepsilon_{V,\mathrm{cyc}} \approx \ln 4 = 1.39$. With
$N_{\mathrm{knee}} \approx 300$–380 for a route-B-only knee, $k_V \approx 1$e-3
(B1) or 2e-3 (B2). The real $\Delta\varepsilon_{V,\mathrm{cyc}}$ per window
must be read off the model (step 0 below). If $\bar s$ stays partly lithiated
at a window's bottom, $\Delta\varepsilon_V$ shrinks.

**Prediction.** One $k_V$ should put every route-B-dominated cell's knee at
the same $\mathcal{E}_{\mathrm{knee}}$. That is the model-side counterpart of
the common ~15,000 µm cumulative swing in the data. Because $D$ uses the Si
strain, not cell throughput, a window that mostly cycles graphite
accumulates little $\mathcal{E}$. That matches the published finding that
window position matters (Jossen group 2024; Dressler & Dahn 2025).

### 3.4 Optional Coffin–Manson-type amplitude exponent (later)

$$
D = |\dot\varepsilon_V|\left(\frac{|\varepsilon_V - \bar\varepsilon_V|}{\varepsilon_{\mathrm{ref}}}\right)^{p-1},
$$

where $\bar\varepsilon_V$ is a slow running mean (one extra relaxation state).
$p = 1$ recovers B1. Use it only if the cells can't be fitted with $p = 1$.

## 4. Expected effect per condition

| Condition | Route A (Paris, diffusion stress) | Route B (Si strain fatigue) | Knee set by |
|---|---|---|---|
| CELL064, full, 25 °C | large (σ ~40 MPa, deep delithiation) | ~same per cycle as CELL017 | A + B, earlier (~101) |
| CELL017, full, 45 °C | small (σ ~16–20 MPa) | same swing per cycle as CELL064 | B (~310) |
| CELL009, 15–95%, 25 °C | ~0 during C/3 cycling (σ ~2 MPa) | smaller swing per cycle, more cycles per EFC | B (~290) |
| CELL026, 20–80%, 25 °C | ~0 | smaller swing again | B (~293) |

**Expected outcome.** One shared set could drop both corrections:
- `SI_CRACK_EAC`: route A's own stress drop at 45 °C does the job.
- The CELL009 crack multiplier ×11: replaced by route B.

**Caveat: strain vs throughput.** Swing per EFC is ~48–56 µm for all four
cells, so cumulative swing and Ah throughput can't be separated with these
cells alone. A throughput variant ($D \propto |j_{\mathrm{Si}}|$) would fit
them equally well. Two things discriminate:
- the published result that window position matters, which favours Si strain;
- the 2C cells (CELL056/074, knees ~143–148 EFC) and the pressure series.

## 5. Alternative within the same idea: breathing SEI

Si breathing cracks the SEI shell each cycle and exposes fresh surface:

$$
\frac{\partial c_{\mathrm{SEI}}}{\partial t}\bigg|_{\mathrm{breath}}
= k_{\mathrm{br}}\; a\; D(\dot\varepsilon_V)\; \frac{L_{\mathrm{SEI,0}}}{L_{\mathrm{SEI}}},
$$

with the existing LLI and pore-closure bookkeeping.

- **Same knee prediction** as route B.
- **How to tell them apart:** pre-knee LLI vs Si LAM timing. Route B grows
  crack area first, which then drives crack SEI. The breathing-SEI variant
  gives LLI growing linearly with throughput from BoL.

## 6. Implementation sketch (separate branch)

- **Option**, per phase, default unchanged: `"particle cracking growth"`,
  with values:
  - `"Paris"`: default, route A only;
  - `"Paris + strain fatigue"`: eq. (3) with B1;
  - `"Paris + strain fatigue (contraction)"`: eq. (3) with B2.
- **New parameters:**
  - `"{Phase}: {Domain} electrode strain-fatigue cracking constant k_V [-]"`;
  - optionally `"... strain-fatigue cracking activation energy [J.mol-1]"`,
    default 0.
- **Code:**
  - `crack_propagation.py`: add route B to `dl_cr`.
  - $\dot\varepsilon_V$ from eq. (2), using the phase's interfacial current
    and `t_change` with its analytic derivative (power law), or a
    `pybamm.d_dx`-style derivative for a generic `t_change`.
- **New output variables:**
  - `"... particle strain-fatigue cracking rate [m.s-1]"`, separate from the
    Paris rate;
  - `"... cumulative volumetric strain [-]"`, i.e. $\mathcal{E}$, as an extra
    ODE state for diagnostics only.
- **Fork wiring:** `C64_SI_STRAIN_CRACK_KV` (unset = off) and
  `C64_SI_STRAIN_CRACK_MODE=abs|contraction`.

## 7. Refit plan: CELL064 (25 °C), CELL017 (45 °C), CELL009 (15–95%)

**Goal.** One shared parameter set, with only T, F0 and (for now) WIDTH per
condition, that fits all three cells:
- without `SI_CRACK_EAC` or a CELL009 crack multiplier;
- with RPTs aligned in all runs (`--aligned`);
- with CELL026 (20–80%) held out as a blind test.

**Starting point:** `crack_baseline_v3.env` (θ = 1, Gr yield 0.025 /
56.5 kJ/mol, shared `SI_CRACK_MULT` 1e3).

**Step 0: housekeeping before any fitting.**
1. **One SoH reference everywhere.** Use the conditioned RPT1 (a C/20 after a
   standard CC-CV) for the full-window runs too. Today only partial windows
   use it, and the full-window runs carry a related offset (the unexplained
   ~−2.2pp at CELL064 EFC 49).
2. **Decide the EFC basis.** Data: throughput / 5.0 Ah. Model: / 5.19 Ah,
   ~3.8% apart. The aligned RPT targets should use one basis.
3. **Diagnostic runs with B off.** Run v3 on all three cells, outputting
   $\mathcal{E}$ (B1 and B2) and $\Delta\varepsilon_{V,\mathrm{cyc}}$. This
   gives the actual Si strain per cycle in each window, and the
   $\mathcal{E}$ at each real knee, for eq. (5).

**Step 1: shared cathode LAM (removes CELL009's ×0.3).** Make the cathode's
stress-driven ("linear") LAM saturate after an early loss. The preferred
route is crack-length stress relief on the positive electrode; the fallback
is the existing power-law damping. Both are in section 10.1.
- **Fit target:** the real LAM_PE across all three cells, ~4–7% early, then
  flat.
- **Check:** SoH before the knee is unchanged.

**Step 2: fit $k_V$ where route B dominates (CELL009).**
- **Setup:** route A at the shared 1e3 (negligible there), WIDTH 1e-3,
  cathode law from step 1.
- **Fit:** $k_V$ (B2, then B1) to the knee (~290) and post-knee SoH, with
  eq. (5) as the starting value.
- **Grid:** ~4 runs per variant around the eq. (5) estimate.

**Step 3: CELL017 at 45 °C, with `SI_CRACK_EAC` = 0.**
- **Prediction:** with the same $k_V$, a knee at ~310.
- **If it knees early:** route A is too strong at 45 °C without the Arrhenius
  slowdown. Lower the shared Paris $k_{\mathrm{cr}}$ (or raise $m_{\mathrm{cr}}$,
  which suppresses the lower 45 °C stress more), then go to step 4.
- **If it knees late:** add a mild positive $E_V$, as the first sign that
  route B needs temperature dependence.

**Step 4: CELL064 at 25 °C.**
- **Target:** route A + route B knee at ~101.
- **Re-tune:** only the shared route-A parameters ($k_{\mathrm{cr}}$, and
  $m_{\mathrm{cr}}$ if step 3 moved it), never $k_V$. This is the only cell
  where A matters, so it pins A.

**Step 5: joint refinement.** A small grid over the two shared knobs,
$k_V$ (×{0.7, 1, 1.4}) × $k_{\mathrm{cr}}$ (×{0.7, 1, 1.4}): 9 points × 3 cells =
27 aligned runs.
- **Objective, per cell:**
  - mean |SoH gap|;
  - knee EFC error;
  - LLI and Si/Gr LAM at the last RPT vs DMA;
  - late-RPT voltage / dV/dSOC RMSE.
- **Pick:** the point that minimises the worst cell, not the average.
- **Then re-check** the secondary shared knobs that interact with knee
  timing: Si yield and τ (the isolation gate), and Gr yield (Gr LAM).

**Step 6: blind validation and extension.**
- **CELL026 (20–80%), coulomb-counted protocol, no tuning.** Predict the
  knee (real ~293) and the LLI / LAM split. This is the decisive test for
  Si strain over throughput: its Si strain per cycle is lower again.
- **Optional:** the 2C cells (CELL056/074) and the pressure series, to separate
  strain-driven from throughput-driven damage and to bring in rate effects
  on route A.

**Success criteria.**
- **Knees:** one set reproduces CELL064, CELL017 and CELL009 within ±15 EFC,
  with no Arrhenius slowdown on cracking and no per-cell crack multiplier.
- **SoH:** gaps within ±3pp before the knee, mean |gap| ≲ 3pp per cell.
- **LAM / LLI:** last-RPT values within ~5pp of the DMA for each cell.
- **Blind test:** CELL026's knee predicted within ±25 EFC.

**What would falsify route B.**
- One $k_V$ can't serve both CELL009 and CELL017.
- CELL009 needs a $k_V$ that knees CELL064 well before EFC 101, even with
  route A reduced.
- B1 and B2 both fail on CELL026 while a throughput variant succeeds. That
  would favour throughput-driven damage, or breathing SEI (section 5).

**Remaining per-condition parameter.** WIDTH is still 5e-3 / 3e-4 / 1e-3
(25 °C full / 45 °C / 15–95%). See section 10.2.

## 8. Diagnostics

- **Already in `rpt_soc_plots.py` crack diagnostics:** crack length,
  roughness, porosity, peak Si stress per cycle, and crack vs bulk SEI LLI.
- **To add:** the route-A vs route-B split of $dl/dt$; $\mathcal{E}$ and
  $\Delta\varepsilon_{V,\mathrm{cyc}}$ per cycle; and $(l/l_0)$ at the knee.

## 9. Literature check (2026-10-01, delegated search)

Some full texts were paywalled, so a few equation details come from
abstracts.

**Verdict: the combination doesn't appear to be published.** No continuum
Si/Gr composite model was found that combines gradient-stress Paris cracking
(A) with volumetric-strain / fatigue cracking (B).

**Route A only (Si/Gr):**
- Bonkile, O'Kane, Planella, Marinescu, Offer et al. 2024, *J. Power Sources*
  606, 234256 (PyBaMM): stress cracking, LAM and SEI on cracks per phase.
  Si loss is self-limiting; higher depth of discharge (DoD) accelerates it.
- Philipp, Köbbing, Karger, Jossen, Latz, Horstmann 2026, arXiv:2604.26545
  (PyBaMM SPMe): **negligible cracking at intermediate SoC**. That is the same
  failure v3 shows on CELL009, and it contradicts the CELL009 / CELL026 data.
- Background: O'Kane et al. 2022 (*PCCP*), Deshpande et al. 2012 (*JES*),
  Purewal et al. 2014, Ai et al. 2022, Reniers et al. 2019.

**Route B, in pieces:**
- SEI break-and-repair driven by surface expansion, graphite: Laresgoiti et
  al. 2015 (*J. Power Sources* 300); Ekström & Lindbergh 2015 (*JES* 162
  A1003), with a cracked-SEI fraction set by the expansion rate.
- Single Si particle: Pinson & Bazant 2013 (surface-area change); von
  Kolzenberg, Latz, Horstmann 2022 (breathing SEI fracture).
- Empirical: Karger, O'Kane et al. 2024 (*JES* 171 090512), Eq. 22 with a
  ΔV^γ3 crack-rate term. γ3 was fitted to 0, which the authors attribute to
  their constant-depth data (see below).
- Related: Pannala et al. 2024 (*JES* 171 010532), asymmetric tensile /
  compressive LAM (already in this code base). Zhang et al. 2023 (*Int. J.
  Fatigue*), Miner's-rule fatigue in Si.

**Experimental support:**
- Kirkaldy et al. 2022 (*ACS AEM*): 0–30% SoC at 40 °C loses 80% of Si
  capacity in ~400 EFC.
- Dressler, Ingham, Dahn 2025 (*JES*): lifetime is cycle-count driven;
  graphite-only windows last long.
- De Sutter et al. 2018 (*Energies*): reduced DoD extends Si-alloy life.
- Jossen group 2024 (*JES*): at a fixed 50% DoD, window position matters.

**Reusable forms:**
- Philipp's $\partial_t l = k\,\sigma_{t,\mathrm{surf}}^m$, with the stress
  replaced by Si strain or strain rate.
- Ekström's expansion-rate cracked-SEI fraction.
- Karger's power law.
- The existing `"current-driven"` LAM option, as a throughput baseline.

**Points the proposal addresses:**
- **Philipp 2026:** route A alone gives no intermediate-SoC cracking. That is
  the motivation for B.
- **Karger 2024 (γ3 = 0): resolved by the authors themselves.** Their
  discussion attributes the zero exponent to every check-up using the same
  cycling depth, so the volume change only varied through LLI-driven loss of
  used anode capacity. They state that their cracking model lacks a
  volume-expansion dependence and is only valid for full-SoC cycling. Route
  B fills that gap, and CELL009 / CELL026 are the varying-depth data they
  lacked.
- **Jossen 2024: window position matters.** Route B is driven by **Si**
  strain (eq. 2), not total throughput. That favours eq. (3) over a
  $|j_{\mathrm{Si}}|$ throughput variant.

## 10. Caveats for calling it a single parameter set

Route B unifies the **cracking** side: one Paris law and one $k_V$ for all
conditions, with no `SI_CRACK_EAC` and no CELL009 crack multiplier. Only
T, F0 (the measured k at RPT1) and the cycling protocol are legitimately
per-condition. Two fitted items still have to be dealt with before the
whole set can be called single.

### 10.1 Cathode LAM must be unified too

**Problem.** The real cathode LAM is early and self-limiting, the same in all
windows and at both temperatures:
- CELL064: ~7% by EFC 49, then flat (~6–7%);
- CELL017: ~3–5% throughout;
- CELL009: ~4–5% after EFC 260.

The model's stress-driven ("linear") positive LAM,
$\dot\varepsilon_{s,p} = -\beta_p\,(\sigma_{h,p}^{+}/\sigma_{\mathrm{crit},p})^{m_p}$,
grows ~linearly with cycling. It lands right at CELL064's EFC 200 (6.5%),
but reaches 18% (v3, EFC 402) or 28% (EFC 752) on CELL009. CELL009's ×0.3 is
a patch for that, not a shared parameter.

**Preferred fix: crack-length stress relief on the positive electrode.**
This is physically the early NMC secondary-particle cracking that then
stabilises.

- **Current code:**
  - The positive particle mechanics in the fork are `"swelling only"`, so
    there is no positive crack length.
  - The crack stress relief in `crack_propagation.py`, $\max(1 - l/R, 0)$,
    uses the particle radius as the cap and only acts on the crack-growth
    stress.
  - Stress-driven LAM (`loss_active_material.py`) uses its own stress and
    does not see the crack length.
- **Proposed change:**
  1. Switch the positive mechanics to `"swelling and cracking"`.
  2. Add a parameter `"Positive electrode maximum crack length [m]"`,
     $l_{\max,p}$, used as the cap instead of $R$ (default = $R$, so
     behaviour is unchanged when unset).
  3. Let the stress-driven LAM see the same relief:

$$
\dot\varepsilon_{s,p} = -\beta_p\left(\frac{\sigma_{h,p}^{+}\,\phi_p(l_p)}{\sigma_{\mathrm{crit},p}}\right)^{m_p},
\qquad
\phi_p(l_p) = \max\!\left(1 - \frac{l_p}{l_{\max,p}},\, 0\right).
$$

  Cracks grow early, during the high-stress first cycles. As $l_p \to
  l_{\max,p}$ the stress is relieved and the LAM stops. That gives an early
  loss and then a plateau, set by how fast $l_p$ reaches $l_{\max,p}$
  (positive $k_{\mathrm{cr}}$, $m_{\mathrm{cr}}$, $l_{0,p}$) and by
  $\beta_p$. Opt in through a per-phase option such as
  `"stress-driven LAM crack relief": "true"`; default off.
- **Fit:** $l_{\max,p}$, and the positive crack rate, to the ~4–7% plateau
  across all three cells. It is shared by construction, since the positive
  stress at 4.15–4.195 V tops is similar in all windows.
- **Side effects to check:**
  - Positive cracking also adds positive SEI or CEI only if "SEI on cracks"
    is enabled for the positive. Keep it off.
  - Cell thickness, through positive swelling, is unchanged.

**No-code fallback.** Use the existing `"stress-driven LAM damping": "power"`
on the positive electrode, with remaining-fraction damping
$(\varepsilon_s/\varepsilon_{s,0})^{n_p}$, a large $n_p$ (~30–60) and a
higher $\beta_p$. The rate collapses after a few percent of loss, so this
also plateaus. It is less physical, since the saturation is tied to the
amount lost rather than to crack relief, but it is a parameter-only test of
whether a plateau law fits all three cells.

### 10.2 Width remains a per-condition calibration

**Problem.** The pore-buffering transition width is 5e-3 (CELL064), 3e-4
(CELL017) and 1e-3 (CELL009), an order-of-magnitude spread. It currently
acts as a timing knob for the k / expansion transition and the C/3 step,
not as a material constant. It also nudges the knee through porosity and
the isolation gate.

**Honest framing until resolved.** Report width alongside F0 as a
per-condition calibration of the expansion / pore-buffering transition.
Don't count it as part of the shared degradation set.

**Routes to unify it (separate work):**
1. **Soften the hard headroom cap**, $\min\big((1-k)\,dv_{\mathrm{solid}},\ \mathrm{headroom}\big)$,
   in `reaction_driven_porosity.py`. It is the likely source of the abrupt
   C/3 step in partial windows. If the abrupt part comes from the cap, a
   single, larger shared width might then work.
2. **Tie the transition to a physical variable.** If the transition really
   follows a temperature- or window-dependent quantity (binder compliance at
   45 °C, the Si utilisation window), express width as a function of it
   rather than per condition.
3. **Check after the route-B refit.** Route B changes when and how fast
   porosity closes in each condition, so the fitted widths may move closer
   together on their own.
