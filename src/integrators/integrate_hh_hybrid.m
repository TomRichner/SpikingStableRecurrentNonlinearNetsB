function [t, Y] = integrate_hh_hybrid(odefun, tspan, y0, opts, jp) %#ok<INUSL>
%INTEGRATE_HH_HYBRID  Event-aware fixed-step RK4 for SRNNModelHH (Campagnola).
%
%   [t, Y] = INTEGRATE_HH_HYBRID(odefun, tspan, y0, opts, jp)
%
%   Thin wrapper that drives the generic event integrator (integrate_hh_events)
%   with the SRNNModelHH per-(pre-neuron, post-type) Tsodyks-Markram / SFA /
%   conductance jump logic (apply_hh_jumps). @ode45-compatible signature so it
%   slots into SRNNModelBase.run() and the Benettin reshoot; opts is accepted
%   for parity and IGNORED (fixed step).
%
%   jp is the jump-parameter struct built by SRNNModelHH.get_params (field
%   .jump), with fields:
%     N, K            - neuron count, cell-type count
%     n_a, n_ad       - SFA timescales, number of adapting neurons
%     ad_idx          - (n_ad x 1) neuron indices of adapting neurons
%     type_of         - (N x 1) presynaptic type label per neuron
%     Wabs            - (N x N) effective |weight| (pre x post), already scaled
%     kappa           - (K x K) facilitation increment coeff (pre x post type)
%     p0_mat          - (N x K) baseline release prob (used when STF is off)
%     a_incr          - scalar SFA increment per own spike
%     V_th, V_reset   - spike-detect threshold and re-arm (hysteresis) level (mV)
%     has_a, has_b, has_p - mechanism-present flags (g conductance always present)
%
%   See also: integrate_hh_events, apply_hh_jumps, SRNNModelHH.

    jump_fn = @(y, spiked) apply_hh_jumps(y, spiked, jp);
    [t, Y] = integrate_hh_events(odefun, tspan, y0, jp.N, jump_fn, jp.V_th, jp.V_reset);
end
