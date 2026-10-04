# Cracking-driven degradation model: equations and temperature dependence

Equations for the 25 °C / 45 °C model being fitted to CELL064 (25 °C, 103
kPa) and CELL017 (45 °C, 103 kPa). It runs on the high_temp_45C fork
(`../high_temp_45C/cell064_degradation_fit_45C.py`), driven by
`crack_growth_sweep.py`.

- **Core-library changes** are in `../../../CHANGES.md` (items 1–14); they're
  referenced here as [C#].
- **Script-level choices** (parameter values, Arrhenius wrappers,
  temperature switches) live in the fork and are documented here.
- **One parameter set, both temperatures:** the only thing that differs
  between the 25 °C and 45 °C runs is the operating temperature
  $T$ (`C64_T_AMBIENT_K`, isothermal) and the per-cell BoL calibration
  $F_0$ (§6).

Notation: Si = secondary phase, Gr = primary phase of the negative
electrode. Every Arrhenius factor below uses $T_{\mathrm{ref}} = 298.15$ K,
so it is **exactly 1 at 25 °C**; activation energies only act at 45 °C.

---

## 1. SEI growth

**Si, solvent-diffusion limited** (bulk and on cracks, the same law):

$$
j_{\mathrm{SEI}} = -\frac{F\, D_{\mathrm{sol}}\, c_{\mathrm{sol}}}{L_{\mathrm{SEI}}}\;
\exp\!\left[\frac{E_{\mathrm{SEI,Si}}}{R}\left(\frac{1}{T_{\mathrm{ref}}}-\frac{1}{T}\right)\right],
\qquad D_{\mathrm{sol}} = 2.5\times10^{-22}\cdot\texttt{SOLVENT\_MULT}.
$$

**Gr, EC-reaction limited**, with $E_{\mathrm{SEI,Gr}} = 38$ kJ/mol.

- **SEI on cracks** grows on the crack area $a_{\mathrm{cr}} = a\,(\rho - 1)$,
  where $\rho$ is the roughness ratio (§2), with its own thickness
  $L_{\mathrm{SEI,cr}}$.
- **Film overpotential** uses **bulk** $L_{\mathrm{SEI}}$ only
  (`"SEI film resistance": "distributed"`), so crack SEI adds LLI and pore
  closure but not film resistance.

## 2. Si crack growth (Paris law with stress relief [C2])

$$
\frac{dl_{\mathrm{cr}}}{dt} = \frac{k_{\mathrm{cr}}(T)}{3600}\,
\Big(\sigma_{t,\mathrm{eff}}\, b_{\mathrm{cr}}\sqrt{\pi l_{\mathrm{cr}}}\Big)^{m_{\mathrm{cr}}},
\qquad
\sigma_{t,\mathrm{eff}} = \sigma_t\,\max\!\left(1 - \frac{l_{\mathrm{cr}}}{R_{\mathrm{Si}}}, 0\right),
$$

$$
\rho = 1 + 2\, l_{\mathrm{cr}}\, w_{\mathrm{cr}}\, \rho_{\mathrm{cr}},
\qquad
k_{\mathrm{cr}}(T) = k_{\mathrm{cr},0}\cdot\texttt{SI\_CRACK\_MULT}\cdot
\exp\!\left[\frac{E_{\mathrm{cr}}}{R}\left(\frac{1}{T}-\frac{1}{T_{\mathrm{ref}}}\right)\right].
$$

- **Sign convention:** $(1/T - 1/T_{\mathrm{ref}})$, so a positive
  $E_{\mathrm{cr}}$ (`SI_CRACK_EAC`) *slows* cracking at high T (lithiated Si
  is more ductile). This is the same convention as
  `silicon_cracking_rate_Ai2020`.
- **Cap:** stress relief caps $l_{\mathrm{cr}} \to R_{\mathrm{Si}}$ (~300 nm),
  so crack growth saturates.
- **Initial length:** $l_{\mathrm{cr},0}$ = 2e-8 m (literature; override
  with `C64_SI_CRACK_L0`).
- **Temperature enters twice:** directly through $E_{\mathrm{cr}}$, and
  indirectly through $\sigma_t$. Faster Si diffusion at 45 °C flattens the
  concentration gradients, so pre-knee $\sigma_t$ is ~16–20 MPa against
  ~43 MPa at 25 °C. With $m_{\mathrm{cr}} = 2.2$, that alone roughly doubles
  the knee EFC.
- **Literature $k_{\mathrm{cr},0}$ gives static cracks**
  ($l/l_0 \approx 1.005$ over life). `SI_CRACK_MULT` ~1.3e3 is needed for
  crack growth to matter.

## 3. Porosity and the isolation gate

- **Porosity closure:** structural porosity falls with
  $L_{\mathrm{SEI}} + L_{\mathrm{SEI,cr}}(\rho-1)$ (plus plating), with a
  softplus floor at $\varepsilon_{\min}$ = 0.035.
- **Gate target and relaxation** (per phase):

$$
h = \mathrm{clip}\!\left(\frac{\varepsilon - \varepsilon_{\min}}{\varepsilon_0-\varepsilon_{\min}}, 0, 1\right),
\qquad
G^{\ast} = (1-h)^{\eta},
\qquad
\frac{dG}{dt} = \frac{G^{\ast} - G}{\tau}.
$$

  $\eta$ = `SI/GR_LAM_ISO_EXPONENT` sets how sharp the gate is, and $\tau$ =
  `SI/GR_TAU_LAM_ISO` sets the knee width.

**SEI redirect to LAM.** SEI film growth is gated, and the redirected
current becomes isolation LAM:

$$
\frac{dc_{\mathrm{SEI}}}{dt} \leftarrow (1-G)\,\frac{dc_{\mathrm{SEI}}}{dt},
\qquad
\dot\varepsilon_{\mathrm{iso}} = y(T)\, G\, \frac{a\, j_{\mathrm{SEI}}}{F}\,
\mathrm{clip}\!\left(\frac{\varepsilon_s}{\varepsilon_{s,0}},0,1\right) \le 0,
$$

$$
y(T) = y_{\mathrm{ref}}\exp\!\left[\frac{E_{y}}{R}\left(\frac{1}{T_{\mathrm{ref}}}-\frac{1}{T}\right)\right]
\quad(\texttt{SI/GR\_REDIRECT\_LAM\_YIELD},\ \texttt{SI/GR\_YIELD\_EAC}).
$$

- **Evaluation:** the runs are isothermal, so $y(T)$ is a constant evaluated
  at the operating T.
- **Graphite:** $E_{y,\mathrm{Gr}} \approx 43$ kJ/mol, from requiring one base
  yield to hit both measured last-RPT Gr LAM values: 9.1% at 25 °C and
  10.5% at 45 °C.

## 4. Stress-driven Si LAM [C12]

$$
\dot\varepsilon_{\mathrm{stress}} = -\beta_{\mathrm{LAM}}(T)
\left(\frac{\sigma_h^{+}\, r^{\,n}}{\sigma_{\mathrm{crit}}}\right)^{m_{\mathrm{LAM}}},
\qquad
\beta_{\mathrm{LAM}}(T) = \beta_0 \exp\!\left[\frac{E_{\mathrm{LAM}}}{R}\left(\frac{1}{T_{\mathrm{ref}}}-\frac{1}{T}\right)\right].
$$

- **Damper exponent:** $r = \varepsilon_s/\varepsilon_{s,0}$ and $n$ =
  `SI_STRESS_LAM_DAMP_EXP` (option `"power"`; $n = 1$ is the default
  `"linear"`).
- **Activation energy:** $E_{\mathrm{LAM}}$ = `SI_LAM_EAC` (fork default
  40 kJ/mol).

## 5. LLI from LAM, with isolation lithium trapping [C13]

LAM removes lithium at the local r-averaged content (Sulzer et al. 2021,
eq. 37): $\dot n_{\mathrm{Li,LAM}} = -V\langle \bar c_s\,\dot\varepsilon_s\rangle_x$.

For Si isolation LAM only, the baseline also uses "isolation lithium
trapping", with $\theta = 1$ (`C64_SI_ISO_LI_TRAP`). Isolated Si leaves
fully lithiated:

$$
c_{\mathrm{trap}} = (1-\theta)\,\bar c_s + \theta\, c_{s,\max}.
$$

The excess $\theta(c_{s,\max}-\bar c_s)\,|\dot\varepsilon_{\mathrm{iso}}|/\varepsilon_s$
is removed uniformly in r from the remaining Si particles, and added to the
LAM-trapped LLI.
- **No current, no SEI film, no pore filling,** so the knee and Si LAM
  are unchanged.
- **Stress-driven LAM keeps eq. 37.**

## 6. Pore buffering and expansion scale k

$$
k = f(\varepsilon_{\mathrm{struct}}) = \frac{1}{1 + K\,C_{\mathrm{pore}}(\varepsilon_{\mathrm{struct}})},
\qquad
C_{\mathrm{pore}} = C_{\max}\left(1 - e^{-\,\mathrm{softplus}(\varepsilon_{\mathrm{struct}}-\varepsilon_{\min,\mathrm{tr}})/w}\right),
\qquad K C_{\max} = \frac{1-F_0}{F_0}.
$$

- **$F_0$ is the BoL k plateau,** set per cell to the measured k at its
  first RPT: 0.70 for CELL064 and 0.80 for CELL017.
- **k is monotonic by construction.** $\varepsilon_{\mathrm{struct}}$ only
  decreases, so k only rises from $F_0$ towards 1. The measured CELL017 k
  dip before the knee (0.80 → 0.72, with the Si window fully used) can't be
  reproduced. Deferred (see FINDINGS.md).

## 7. Run conventions

- **EFC cut-offs** default by temperature: 203 at 25 °C and 403 at 45 °C,
  just past the last real RPTs at 197.4 and 400.0. Each 50-cycle batch ends
  with its C/20 RPT, so the model's last RPT always lies beyond the last
  real one.
- **Scoring against real data:** at 298.15 K the fork uses CELL064; at
  318.15 K it uses CELL017 (`../../degradation_test_matrix/experimental_data/CELL017_45C/`).
- **LLI** is compared on particle inventory (`"Total lithium lost from
  particles"`), since the DMA's `charge_LLI` is an inventory quantity.

## 8. Current parameter set: crack_baseline_v3 (2026-09-30, `crack_baseline_v3.env`, = R7 / R7a)

| Parameter | Value | Effect at 45 °C |
|---|---|---|
| `SOLVENT_MULT` | 5 | – |
| `SI_ESEI` | 35 kJ/mol | Si SEI ×1.88 |
| `GR_ESEI` | 38 kJ/mol | Gr SEI ×2.62 |
| `SI_CRACK_MULT` | 1.0e3 | – |
| `SI_CRACK_EAC` | 90 kJ/mol | Si cracking ÷8.3 |
| `SI_TAU_LAM_ISO` / `SI_LAM_ISO_EXPONENT` | 3.5e8 s / 20 | – |
| `SI_REDIRECT_LAM_YIELD` / `SI_YIELD_EAC` | 0.25 / 35 kJ/mol | 0.61 |
| `SI_ISO_LI_TRAP` (θ) | 1.0 | – |
| `GR_REDIRECT_TO_LAM` / `GR_LAM_ISO_EXPONENT` | 1 / 10 | – |
| `GR_REDIRECT_LAM_YIELD` / `GR_YIELD_EAC` | 0.025 / 56.5 kJ/mol | 0.105 |
| `SI_STRESS_LAM_DAMP_EXP` | 0.25 (`"power"`) | – |
| `SI_LAM_EAC` | 40 kJ/mol | ×2.8 |
| `F0` | 0.70 (CELL064) / 0.80 (CELL017) | per-cell BoL calibration |
| `WIDTH` (pore-buffering transition width) | 5e-3 (25 °C) / 3e-4 (45 °C) | per-temperature |

- **Differences between conditions:** only T, F0 and WIDTH differ. WIDTH
  is per-temperature (FINDINGS.md, R4 width test).
- **Fitting conditions:** fitted with the item-14 core fix, and with model
  RPTs aligned to the real RPT EFCs (`rpt_soc_plots.py 25|45 --aligned`).
- **Previous set:** crack_baseline_v2 (R3c) is v3 without θ, with Gr yield
  0.018 and `GR_YIELD_EAC` 43 kJ/mol.
