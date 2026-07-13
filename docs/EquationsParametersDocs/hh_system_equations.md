# Spiking HH network — governing equations

Biophysical (spiking) analog of the rate model `SRNNModelCellTypes`
(FractionalResevoir). A randomly connected network of cortical Hodgkin–Huxley
neurons of `K` cell types (default `pyr, pvalb, sst, vip`, from Campagnola 2022)
carrying spike-frequency adaptation (SFA), short-term synaptic depression (STD),
and short-term facilitation (STF). Implemented in `src/SRNNModelHH.m`.

**Units:** ms, mV, µF/cm², mS/cm², µA/cm². Campagnola time constants (s) are
converted ×1000 to ms on load.

Indices: neuron `i, j ∈ {1..N}`; cell type `θᵢ ∈ {1..K}`; SFA timescale `ℓ`.
Synaptic resource variables carry a `(pre-neuron j, post-type q)` index.
Weight orientation `W(j,i)` = |conductance| from pre `j` to post `i` (pre × post).

## Hybrid system: smooth flow + spike jumps

The system is a **hybrid** dynamical system. Between spikes the full state
evolves by the smooth ODEs below (`SRNNModelHH.dynamics_hh`). A **spike** of
neuron `j` is an upward crossing of `V_j` through `V_th` (with `V_reset`
hysteresis); at each spike, discrete jumps are applied (`integrate_hh_hybrid`).

### HH membrane (Traub–Miles kinetics)
$$C_m \frac{dV_i}{dt} = I^{ext}_i - I^{Na}_i - I^{K}_i - I^{L}_i - I^{SFA}_i - I^{syn}_i$$
$$I^{Na}_i = \bar g_{Na}\, m_i^3 h_i (V_i - E_{Na}), \quad
  I^{K}_i = \bar g_{K}\, n_i^4 (V_i - E_{K}), \quad
  I^{L}_i = g_L (V_i - E_L)$$
$$\dot m_i = \alpha_m(V_i)(1-m_i) - \beta_m(V_i)m_i \quad\text{(same for } h,n)$$
`α/β` are the Traub–Miles rates with removable-singularity patches
(`hh_gating_rates.m`), reused from `CorticalSpreadDepolarizationModel/model_int.m`.
Fixed reversals: `E_Na=+50, E_K=-100, E_L=-67` (ionic gradients constant).

### SFA — spike-triggered K-adaptation current (adapting types only)
$$I^{SFA}_i = c_i \Big(\sum_{\ell} a_{i\ell}\Big)(V_i - E_K), \qquad
  \frac{da_{i\ell}}{dt} = -\frac{a_{i\ell}}{\tau^a_\ell}
  \;\;\text{(between spikes)}, \qquad
  a_{i\ell} \leftarrow a_{i\ell} + \Delta_a \;\;\text{(own spike)}$$
Spiking analog of the rate model's linear negative feedback `ȧ = (-a + r)/τ_a`.
`c_i = c_{gain}\,\text{adapt\_index}(\theta_i)`; non-adapting types
(`adapt_index < sfa_min_index`, e.g. fast-spiking pvalb) carry no `a` state.

### STD + STF — event-based Tsodyks–Markram (per pre-neuron j, post-type q)
Availability `b_{j,q}` (rest 1) and release probability `p_{j,q}` (rest
`p^0_{θⱼ,q}`):
$$\frac{db_{j,q}}{dt} = \frac{1 - b_{j,q}}{\tau^{rec}_{θⱼ,q}},\qquad
  \frac{dp_{j,q}}{dt} = \frac{p^0_{θⱼ,q} - p_{j,q}}{\tau^{f}_{θⱼ,q}}
  \qquad\text{(between spikes)}$$
At a spike of `j`, for every post-type `q`:
$$\text{release} = p_{j,q}\, b_{j,q}, \qquad
  b_{j,q} \leftarrow b_{j,q} - \text{release}, \qquad
  p_{j,q} \leftarrow p_{j,q} + \kappa_{θⱼ,q}(1 - p_{j,q})$$
Release `= p·b` couples the two (facilitated synapses drain faster), mirroring the
rate model. STF off ⇒ `p ≡ p⁰` ⇒ pure depression; `p⁰=1` recovers plain STD.

### Synaptic conductance (per post-neuron i, pre-type P)
Sign is realized by the reversal potential (excitatory pre → `E_exc`; inhibitory
pre → `E_inh`), not a signed weight:
$$\frac{dg_{i,P}}{dt} = -\frac{g_{i,P}}{\tau^{syn}_P}
  \quad\text{(between spikes)}, \qquad
  g_{i,P} \leftarrow g_{i,P} + |W_{j,i}|\,\text{release}_{j,θᵢ}
  \;\;\text{(spike of pre } j,\ θⱼ=P)$$
$$I^{syn}_i = \sum_P g_{i,P}(V_i - E^{syn}_P)$$

## State layout
`S = [V(N); m(N); h(N); n(N); a(n_ad·n_a); b(N·K); p(N·K); g(N·K)]`, cursor-walked
with guarded blocks: `a` present iff `n_a>0` (adapting neurons only), `b` iff
`n_b>0`, `p` iff `n_u>0`; conductance `g` always present.
`N_sys = 4N + [n_a>0]·n_ad·n_a + [n_b>0]·N·K + [n_u>0]·N·K + N·K`.

## Largest Lyapunov exponent (Benettin)
Because HH neurons have no hard reset (spikes are smooth fast excursions) and the
only discontinuities are the synaptic/SFA jumps, the LLE is computed by Benettin's
**full nonlinear reshoot**: a perturbed copy is re-integrated over each
renormalisation interval `lya_dt` with the *same* event-aware integrator (so it
detects its own spike times), then the separation is renormalised to `d0`. This
needs no saltation matrices. The Jacobian/QR-spectrum path is deferred (it would
require jump matrices for the hybrid flow); `lya_method='benettin'` only.

## Parameter provenance (Campagnola 2022)
`conn_prob ← conn_prob_adj`; `|W| ← |psp_amplitude|` (normalised, scaled by
`g_syn_scale·level_of_chaos`); `τ_rec ← ml_depression_tau`;
`p⁰ ← ml_release_prob`; `τ_f ← ml_facilitation_tau`;
`κ ← ml_facilitation_amount`; `c ← c_gain·sfa_adaptation_index`;
`τ_a ← sfa_tau`. Loaded by `src/connectivity/load_campagnola_matrices.m`
(4×4 pre×post tables); NaNs → defaults, then clamped.

## Reductions
- `n_a=0` → no SFA; `n_b=0` → no depression (`b≡1`); `n_u=0` → no facilitation
  (`p≡p⁰`). Each removes its state block from `S`.
- `level_of_chaos` scales all weights — the edge-of-chaos scan knob (analog of the
  rate model's abscissa-normalised gain).
