function [m_inf, h_inf, n_inf] = hh_gating_inf(v)
%HH_GATING_INF  Steady-state HH gating values at membrane potential v (mV).
%
%   [m_inf, h_inf, n_inf] = HH_GATING_INF(v) returns the steady-state
%   activation/inactivation of the m, h, n gates, x_inf = alpha/(alpha+beta),
%   using the Traub-Miles rates from HH_GATING_RATES. Used to initialise the
%   gating variables of a neuron resting at v so it starts on its own
%   nullcline rather than transiently spiking at t = 0.
%
%   See also: hh_gating_rates, SRNNModelHH.

    [am, bm, ah, bh, an, bn] = hh_gating_rates(v);
    m_inf = am ./ (am + bm);
    h_inf = ah ./ (ah + bh);
    n_inf = an ./ (an + bn);
end
