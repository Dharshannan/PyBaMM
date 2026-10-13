# Model equations: porosity-gated isolation LAM and strain-fatigue cracking

The equations below are as implemented in the fork's core PyBaMM. Parameter
values are the current shared set **v6** (`stage3_v6.env` on top of
`../crack_baseline_v3.env`).

Code locations:

| Mechanism | File | Function |
|---|---|---|
| Isolation gate and LAM | `src/pybamm/models/submodels/active_material/loss_active_material.py` | `get_coupled_variables`, `set_rhs` |
| SEI redirect | `src/pybamm/models/submodels/interface/sei/sei_growth.py` | `set_rhs` |
| Lithium trapping sink | `src/pybamm/models/submodels/particle/fickian_diffusion.py` | `_isolation_li_trap_sink` |
| Crack growth, routes A and B | `src/pybamm/models/submodels/particle_mechanics/crack_propagation.py` | `get_coupled_variables`; `_strain_fatigue_cracking_rate` (route B, called from `set_rhs`) |

Summary: both mechanisms are **semi-empirical**. Their driving variables are
the model's own physics (porosity from SEI volume and buffered swelling;
strain from Si stoichiometry). The rate laws and gates that turn those
drivers into LAM or crack growth are phenomenological closures with fitted
constants.

---

## 1. Porosity-gated, SEI-induced isolation LAM

**Physical story:**
- SEI fills the pores: both the bulk film and the film on crack surfaces.
- Si swelling buffered into the pores fills them as well.
- As the pore volume runs out, regions of the electrode lose electrolyte
  access.
- Si (and graphite) particles there become electrically or ionically
  isolated. They turn into LAM instead of continuing to grow SEI.

### 1.1 Porosity headroom

This is local in $x$. It is 1 at the as-built porosity and 0 at the floor.

$$
h = \mathrm{clip}\!\left(\frac{\varepsilon - \varepsilon_{\min}}{\varepsilon_0 - \varepsilon_{\min}},\; 0,\; 1\right)
$$

### 1.2 Isolation gate

The gate has an algebraic target and a lagged state that relaxes toward it:

$$
g^* = (1 - h)^{\eta},
\qquad
\frac{dg}{dt} = \frac{g^* - g}{\tau}
$$

### 1.3 SEI redirect (in practice: SEI-current-driven isolation LAM)

As designed, the fraction $g$ of the SEI reaction stops forming film and
becomes LAM instead, so total reaction current is conserved:

$$
\frac{dc_{\mathrm{SEI}}}{dt} = (1 - g)\,\frac{-a\,j_{\mathrm{SEI}}}{F\,z_{\mathrm{SEI}}}
$$

$$
\left.\frac{d\varepsilon_s}{dt}\right|_{\mathrm{iso}}
= Y(T)\; g\; \frac{a\,j_{\mathrm{SEI}}}{F}\;
\mathrm{clip}\!\left(\frac{\varepsilon_s}{\varepsilon_{s,0}},\, 0,\, 1\right)
\qquad (a\,j_{\mathrm{SEI}} < 0 \;\Rightarrow\; \varepsilon_s \text{ decreases})
$$

The yield has an isothermal Arrhenius factor:

$$
Y(T) = Y_0 \exp\!\left[\frac{E_Y}{R}\left(\frac{1}{T_{\mathrm{ref}}} - \frac{1}{T}\right)\right],
\qquad T_{\mathrm{ref}} = 298.15\ \mathrm{K}
$$

**The redirect does two separate things, and only one matters:**

| Role | Equation | Matters? |
|---|---|---|
| (a) **LAM driver:** isolation LAM ∝ SEI reaction current $a\,j_{\mathrm{SEI}}$ | LAM equation above | **Essential.** It produces the knee. Without the redirect option the branch uses $\beta_{\mathrm{LAM,iso}}\,L_{\mathrm{SEI}}$ (SEI thickness) instead, which is why $\beta_{\mathrm{LAM,iso}}$ is inert in v6. |
| (b) **Current conservation:** $(1-g)$ on SEI film growth | $dc_{\mathrm{SEI}}/dt$ above | **Practically none.** $g \lesssim 10^{-3}$ during the knee and $\le$ ~0.04 over the whole test (§1.5), so SEI growth changes by ≤ 0.1–4%. That is invisible in LLI and porosity. |

**Consequences:**
- **The redirect does not divert current in practice.** The model
  effectively *adds* an isolation-LAM channel on top of essentially
  unchanged SEI growth.
- **The lithium budget still closes.** Each channel is booked separately:
  SEI LLI, and LAM-trapped Li via $\bar c$ and $\theta$ (checked by
  `lithium_budget.py`). There is no double-counting; the "redirect" label
  just over-describes what happens.
- **The fitted corner is a tiny gate with a huge yield.** For (b) to
  matter, $g$ would have to be O(0.1–1). Then $Y_0$ would have to be
  ~100–1000× smaller for the same LAM, and SEI growth would visibly slow
  during the knee.
- **The $(1-g)$ factor is optional.** It can be kept for bookkeeping
  consistency or dropped (as in §1.8); results are unchanged to plotting
  precision.
- **Not yet checked numerically.** That would need a small option in
  `sei_growth.py` to switch the factor off.
- **Description:** "SEI-current-driven, porosity-gated isolation LAM".

### 1.4 Lithium trapped in isolated particles

The standard LAM bookkeeping (Sulzer 2021, eq. 37) removes lithium at the
particle-average concentration $\bar c$. The trapping option adds the extra
lithium that a particle isolated partly lithiated would carry:

$$
\frac{dn_{\mathrm{LLI,extra}}}{dt}
= V \left\langle \theta\,(c_{\max} - \bar c)\,
\max\!\left(-\left.\frac{d\varepsilon_s}{dt}\right|_{\mathrm{iso}},\, 0\right) \right\rangle_x
$$

The same amount is removed from the surviving particles as a uniform sink in
$r$. It is guarded by $c/(c + 0.01\,c_{\max})$ so it cannot drive $c$
negative.

### 1.5 The gate integrates rather than relaxes ($\tau$ ≫ test duration)

With $\tau = 3.5\times10^{8}$ s against a ~$7\times10^{6}$ s test, $g \ll g^*$
throughout, so

$$
g(t) \approx \frac{1}{\tau}\int_0^t g^*\,dt'
\qquad\Rightarrow\qquad
\left|\left.\frac{d\varepsilon_s}{dt}\right|_{\mathrm{iso}}\right|
\approx \frac{Y}{\tau}\left(\int_0^t g^*\,dt'\right)\frac{|a\,j_{\mathrm{SEI}}|}{F}
$$

To first order only $Y(T)/\tau$ matters. The gate acts as a
damage-accumulation memory, not a physical relaxation.

**The degeneracy is between $Y_0$ and $\tau$ only.** Both are
temperature-independent constants, so only their ratio $Y_0/\tau$ enters, at
any temperature. The Arrhenius factor in $Y(T)$ lies outside the
degeneracy. $E_Y$ is separately identifiable from how the ratio changes
between 25 and 45 °C.

**$E_Y$ caveat.** The isolation-LAM rate is $Y(T)\,|a\,j_{\mathrm{SEI}}(T)|$,
and the SEI current has its own activation energy ($E_{\mathrm{SEI}}$ =
35 kJ/mol). The 45 °C data therefore constrain the product's combined
temperature sensitivity. $E_Y$ and $E_{\mathrm{SEI}}$ would trade off if
isolation were their only constraint. $E_{\mathrm{SEI}}$ is pinned
separately by the pre-knee LLI, which reduces but does not remove the
trade-off.

- **Gate size:** $g \le t/\tau \lesssim 0.04$ over the whole test. During
  the knee it is only ~$10^{-3}$ (§1.9).
- **Consequence:** the SEI redirect $(1-g)$ removes a negligible share of
  the film growth. The LAM comes from the large lumped gain $Y$ acting on
  the small $g$.
- **Wall-clock time:** the gate integrates time, not cycles. Rests, CV holds
  and slow RPTs all accumulate isolation while the pores are closed.
- **Why it ended up an integrator:** the lag was introduced as a tuning knob
  to widen and delay the knee. Its default, $\tau = 1$ s, makes the gate
  instantaneous. The fit pushed $\tau$ far past the test duration, so the
  data chose cumulative memory over an instantaneous porosity threshold.
  The physical reading (progressive isolation the longer a region stays
  electrolyte-starved) is plausible but not derived. It is a lagged-state
  control, empirical in nature.

### 1.6 Parameters

| Symbol | Parameter (code name) | v6 Si / Gr | Physical meaning | Status |
|---|---|---|---|---|
| $\varepsilon_{\min}$ | `Negative electrode porosity floor` | 0.035 | porosity at which electrolyte transport is cut off | physical in concept; value set by the fit |
| $\varepsilon_0$ | initial porosity | as built | starting porosity | physical (cell design) |
| $\eta$ | `porosity-isolation LAM exponent` | 20 / 10 | how abruptly isolation sets in as headroom closes (stand-in for a percolation threshold) | empirical shape factor |
| $\tau$ | `porosity-isolation gate time constant [s]` | 3.5e8 / 2e8 s | lag between pore closure and isolation | empirical; degenerate with $Y_0$, only $Y_0/\tau$ matters (§1.5) |
| $Y_0$ | `SEI-redirect-to-LAM yield` | 0.25 / 0.025 | Si volume fraction lost per mole of redirected SEI lithium | lumped, fitted; carries m³/mol implicitly |
| $E_Y$ | fork `C64_SI_YIELD_EAC` / `C64_GR_YIELD_EAC` | 35 / 56.5 kJ/mol | temperature dependence of the yield | fitted |
| $\theta$ | `isolation lithium trapping fraction` | 1.0 (Si) | how full the isolated particles are (1 = isolated at full lithiation) | physical meaning; the fit sits at the bound |

### 1.7 Physical or semi-empirical?

**Physical:**
- The causal chain is physical: pore filling → loss of electrolyte access →
  isolation.
- The lithium budget closes: SEI LLI and LAM-trapped Li are each booked.
  The redirect's current conservation is formally there but numerically
  negligible (§1.3), so in practice isolation LAM is an additional channel
  driven by the SEI current.
- Lithium trapping in isolated particles is a documented effect.
- The porosity it reads is the model's own: SEI volume plus buffered
  swelling.

**Empirical:**
1. **The gate form.** $(1-h)^{\eta}$ is not derived from transport or
   percolation theory, and $\eta = 20$ simply means "very sharp".
2. **$\tau$.** It is much longer than the test, so the "relaxation" is in
   practice an integrator. $Y$ and $\tau$ are not separately identifiable.
3. **$Y$ is not a stoichiometric yield.** It multiplies a molar rate
   (mol m⁻³ s⁻¹) to give a volume-fraction rate, so it carries m³/mol. Its
   size means one mole of redirected lithium isolates far more than its own
   volume of Si. That only makes sense as "a pore blockage isolates a cluster
   of particles".
4. **$E_Y$.** The yield activation energies are fit knobs with no
   independent measurement.
5. **The "redirect" is cosmetic.** With $g \lesssim 10^{-3}$–$0.04$, the
   $(1-g)$ cut on SEI growth is negligible. The physically motivated "SEI
   current becomes LAM" picture is not what the fitted parameters do (§1.3).

### 1.8 Equivalent "isolation dose" form (proposed, not implemented)

Since the gate integrates (§1.5), the same behaviour can be written honestly
as a time-integrated damage law. Define an isolation "dose" $D$: the time
spent with closed pores, weighted by $g^*$, in seconds.

$$
\frac{dD}{dt} = g^*(h) = (1-h)^{\eta},
\qquad
\left.\frac{d\varepsilon_s}{dt}\right|_{\mathrm{iso}}
= k_{\mathrm{iso}}(T)\; D\; \frac{a\,j_{\mathrm{SEI}}}{F}\;
\mathrm{clip}\!\left(\frac{\varepsilon_s}{\varepsilon_{s,0}},\, 0,\, 1\right),
\qquad
k_{\mathrm{iso}}(T) = k_{\mathrm{iso},0}\,\exp\!\left[\frac{E_Y}{R}\left(\frac{1}{T_{\mathrm{ref}}} - \frac{1}{T}\right)\right],
\quad
k_{\mathrm{iso},0} \equiv \frac{Y_0}{\tau}
$$

v6: $k_{\mathrm{iso},0} = 0.25 / 3.5\times10^{8}\ \mathrm{s} \approx 7\times10^{-10}$
(m³ mol⁻¹ s⁻¹ implicitly), with $E_Y$ = 35 kJ/mol (Si) unchanged. Two
identifiable constants ($k_{\mathrm{iso},0}$, $E_Y$) replace three
($Y_0$, $\tau$, $E_Y$).

**What it changes:**
- The pretence of a physical relaxation goes.
- The $Y_0/\tau$ degeneracy goes: one constant, $k_{\mathrm{iso},0}$,
  replaces two. $E_Y$ is kept.
- The law states openly that isolation is a time-integrated damage
  accumulation.

**The SEI redirect** becomes $(1 - k_{\mathrm{iso}} D / Y)$ on the film
growth if it is kept. It can also be dropped: it removes at most ~0.3% of
the bulk SEI current during the knee (§1.5, §1.9).

**Numerical equivalence with the current form.** The two differ only by the
lag's "leak", the factor $e^{-(t-t')/\tau}$ in

$$
g(t) = \frac{1}{\tau}\int_0^t e^{-(t-t')/\tau}\,g^*(t')\,dt'
$$

- **The leak only matters while $g^* > 0$.** $g^*$ is essentially zero
  until the porosity nears its floor (§1.9), so the leak acts only over the
  time since pore closure: $\Delta t \sim 10^{6}$ s.
- **Size:** the relative difference in $g$ is
  $\Delta t / (2\tau) \approx 0.1$–$0.3\%$. That is far below the ~4% bound
  you get from the whole test duration.
- **Effect on the knee:** LAM grows as $\propto \Delta t^{2}$ after closure,
  so a 0.3% gain change shifts the knee by ~0.15% of the closure-to-knee
  time, ≪ 1 EFC.
- **Verdict:** results identical to plotting precision. Not implemented,
  because it is a reformulation, not a physics change.

### 1.9 Why a small gate still gives a fast LAM rise (order of magnitude, 20–80% final run)

Divide the LAM rate by the starting Si volume fraction to get the gain on
the gate:

$$
\frac{1}{\varepsilon_{s,0}}\left|\frac{d\varepsilon_s}{dt}\right|_{\mathrm{iso}}
= K\, g\, \frac{\varepsilon_s}{\varepsilon_{s,0}},
\qquad
K = \frac{Y\,|a\,j_{\mathrm{SEI}}|}{F\,\varepsilon_{s,0}}
\approx \frac{0.25 \times 6\times10^{-5}}{0.058} \approx 3\times10^{-4}\ \mathrm{s^{-1}}
$$

Here $|a\,j_{\mathrm{SEI}}|/F \approx 6\times10^{-5}$ mol m⁻³ s⁻¹: bulk Si
SEI grows by +0.26% of capacity between EFC 232 and 294, over an electrode
volume of $2.9\times10^{-6}$ m³.

**1. Onset.** With $\eta = 20$, $g^*$ is ~0 until the porosity is close to
its floor:

| $h$ | $g^*$ |
|---|---|
| 0.5 | $10^{-6}$ |
| 0.1 | 0.12 |
| 0.05 | 0.36 |

With $\varepsilon_0 = 0.17$ and $\varepsilon_{\min} = 0.035$, $g^*$ switches
on below $\varepsilon \approx 0.048$. The run's porosity goes 0.068 → 0.038
between EFC 232 and 294.

**2. Acceleration.** After closure, $g \approx \Delta t/\tau$, so

$$
\mathrm{LAM} \approx \frac{K\,\Delta t^{2}}{2\tau}
\approx \frac{3\times10^{-4}\,(10^{6})^{2}}{2 \times 3.5\times10^{8}} \approx 0.4
$$

That is ~40% Si LAM within ~45 EFC of closure, while $g$ itself reaches only
~$10^{-3}$. The run shows Si LAM 0.8 → 7.6 → 27 → 48% over EFC 280–340.

**3. Saturation.** The $\varepsilon_s/\varepsilon_{s,0}$ factor rounds the
curve off, giving the S-shape.

**So the knee comes from four pieces:**
- $\eta$ sets the onset (a porosity threshold);
- the integrating gate turns it into a $\Delta t^{2}$ acceleration;
- $Y_0/\tau$ (with $E_Y$ for temperature) sets the steepness;
- Si depletion rounds it off.

---

## 2. Si crack growth: route A (Paris) + route B (strain fatigue)

### 2.1 Total crack-growth rate

$$
\frac{dl}{dt}
= \underbrace{\frac{k_{cr}(T)}{3600}\Big[\sigma_t\left(1 - \tfrac{l}{R}\right) b_{cr}\sqrt{\pi l}\,\Big]^{m_{cr}}}_{\text{A: Paris law on diffusion stress } (\sigma_t \ge 0 \text{ only})}
\; + \;
\underbrace{k_V(T)\,\max\!\left(-\dot\varepsilon_V,\, 0\right)\, l\left(1 - \tfrac{l}{R}\right)}_{\text{B: strain fatigue, contraction only}}
$$

Route B is option `"particle cracking growth": "Paris + strain fatigue
(contraction)"` (CHANGES item 15). The B1 variant (`"Paris + strain fatigue"`)
uses $|\dot\varepsilon_V|$ instead. Route B is added in
`CrackPropagation.set_rhs`.

### 2.2 Volumetric log-strain from the Si volume law

$$
\varepsilon_V = \ln\frac{V}{V_0}(\bar s),
\qquad
\dot\varepsilon_V = \frac{(V/V_0)'}{V/V_0}\,\frac{d\bar s}{dt},
\qquad
\frac{d\bar s}{dt} = \frac{1}{c_{\max}}\,\big\langle \text{particle rhs} \big\rangle_r
$$

$$
\frac{V}{V_0} = 1 + 3\,\bar s^{\,n_{\mathrm{eff}}},
\qquad
n_{\mathrm{eff}} = n_{\mathrm{BOL}} + f_{\mathrm{LAM}}\,(n_{\mathrm{end}} - n_{\mathrm{BOL}})
$$

- $\bar s$ is the r-averaged Si stoichiometry, clamped to $[10^{-6}, 1]$.
- $f_{\mathrm{LAM}} = 1 - \varepsilon_s/\varepsilon_{s,0}$ is the Si LAM
  fraction.
- $(V/V_0)'$ is a central difference with $\delta = 10^{-3}$.

### 2.3 Temperature dependence

Both routes use the same sign convention (Ai2020): a positive activation
energy means **slower when hot**.

$$
k_V(T) = k_{V,0}\exp\!\left[\frac{E_V}{R}\left(\frac{1}{T} - \frac{1}{T_{\mathrm{ref}}}\right)\right],
\qquad
k_{cr}(T) = k_{cr,0}\exp\!\left[\frac{E_{cr}}{R}\left(\frac{1}{T} - \frac{1}{T_{\mathrm{ref}}}\right)\right]
$$

At 45 °C: $k_V$ ÷45 ($E_V$ = +150 kJ/mol) and $k_{cr}$ ÷8.3
($E_{cr}$ = +90 kJ/mol).

### 2.4 Feedback into SEI, porosity and isolation

$$
\rho = 1 + 2\,l\,w_{cr}\,\rho_{cr}
\;\;\longrightarrow\;\; \text{crack SEI area}
\;\;\longrightarrow\;\; \text{LLI} \;+\; \text{pore filling}
\;\;\longrightarrow\;\; \text{isolation gate (§1)}
$$

### 2.5 Parameters (route B)

| Symbol | Parameter (code / fork knob) | v6 | Physical meaning | Status |
|---|---|---|---|---|
| $k_{V,0}$ | `strain-fatigue cracking constant` / `C64_SI_STRAIN_CRACK_KV` | 2.3e-3 | fractional crack growth per unit contraction log-strain | fitted (CELL009 knee) |
| $E_V$ | `C64_SI_STRAIN_CRACK_EA` | +150 kJ/mol | reduced fatigue of warm Si (more plasticity/creep) | fitted at 45 °C; large |
| mode | `C64_SI_STRAIN_CRACK_MODE` | contraction (B2) | cracks open only while the surface is in tension (delithiation) | physically motivated, same logic as route A's $\sigma_t \ge 0$ |
| $l(1 - l/R)$ | – | – | self-similar growth ∝ $l$, capped at the particle radius | phenomenological |
| $n_{\mathrm{BOL}}, n_{\mathrm{end}}$ | `volume change aging-deformation BOL/end exponent` | 1.56 → 2.19 | Si volume vs stoichiometry, steepening with LAM | fitted to expansion data (data-anchored) |
| $R$ | Si particle radius | – | geometric cap | physical |

Route A parameters ($k_{cr,0}$ = base × 1e3, $E_{cr}$ = +90 kJ/mol, $m_{cr}$,
$b_{cr}$) are documented in `../MODEL_EQUATIONS.md` §2.

### 2.6 Physical or semi-empirical?

**Physical:**
- Si's ~300% volume swing is real and well documented as the cause of
  particle fracture.
- The strain comes from the model's own stoichiometry, with no free driver.
- Contraction-only matches where tensile surface stress actually occurs.

**Empirical:**
1. **No amplitude effect.** Damage is linear in cumulative strain: a Miner's
   rule with exponent 1, not a true Coffin–Manson law. Many small swings do
   the same damage as one large swing of the same total. This is why route B
   over-predicts CELL026 (20–80%) relative to CELL009 (15–95%) per EFC.
2. **Mean strain only.** It uses the particle-average strain and ignores the
   core–shell mismatch that actually produces the stress. Part of what $E_V$
   absorbs is probably that missing inhomogeneity.
3. **Phenomenological growth form.** $l(1 - l/R)$ gives exponential growth
   with a cap. It is not derived from fracture mechanics: there is no
   threshold and no fracture toughness.
4. **$E_V$ is an apparent activation energy.** +150 kJ/mol lumps ductility,
   homogenisation and binder/SEI compliance. It is not the barrier of a
   single process, and no measurement supports its size.

Route A is more physical in form: the Paris law is standard fracture
mechanics, applied to the model's own diffusion stress. Its $k_{cr}$ and
$E_{cr}$ are fitted, though.

---

## 3. Bottom line

| | Driving variable | Rate law / gate | Fitted constants |
|---|---|---|---|
| Isolation LAM | porosity (model physics) | $(1-h)^{\eta}$ gate, integrating lag, lumped yield | $\eta$, $Y_0/\tau$ (only the ratio is identifiable), $E_Y$, $\theta$ |
| Strain cracking | Si log-strain (model physics) | linear damage accumulation, self-similar growth | $k_V$, $E_V$ |

Both mechanisms are physically grounded closures, not first-principles laws.

**Weakest points:**
- the isolation gate's $Y_0$–$\tau$ degeneracy, and $E_Y$ trading off
  partly with $E_{\mathrm{SEI}}$;
- route B's missing amplitude and inhomogeneity effects;
- the large $E_V$.

**Candidate upgrades (proposed, not implemented):**
- a Si charge-throughput driver, or a Coffin–Manson amplitude exponent, for
  route B (PROPOSAL §3.4, §4);
- a percolation-type isolation gate.
