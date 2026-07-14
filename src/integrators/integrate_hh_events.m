function [t, Y] = integrate_hh_events(odefun, tspan, y0, N, jump_fn, V_th, V_reset)
%INTEGRATE_HH_EVENTS  Generic event-aware fixed-step RK4 for hybrid HH networks.
%
%   [t, Y] = INTEGRATE_HH_EVENTS(odefun, tspan, y0, N, jump_fn, V_th, V_reset)
%
%   Integrates a hybrid HH network: RK4 on the smooth between-spike flow
%   (odefun) at the native spacing of tspan, plus model-specific discrete jumps
%   applied whenever one of the N membrane potentials (the first N state
%   entries) crosses V_th upward, with V_reset hysteresis so each spike is
%   detected once. Because both the fiducial run and the Benettin reshoot call
%   this with the same jump_fn, spike-time desynchronisation contributes to the
%   divergence, giving a valid largest Lyapunov exponent for the hybrid flow.
%
%   Arguments:
%     odefun   - smooth RHS, dS/dt = odefun(t, S)
%     tspan    - full uniform time grid; solution is returned at these points
%     N        - number of neurons; the membrane potentials V are y(1:N)
%     jump_fn  - y = jump_fn(y, spiked): applies the model's discrete jumps for
%                the neurons flagged in the logical vector `spiked` (length N)
%     V_th     - spike-detection threshold (mV)
%     V_reset  - re-arm (hysteresis) level (mV)
%
%   The caller wraps its layout-specific jump logic in jump_fn (e.g.
%   apply_hh_jumps or apply_hhei_jumps). Fixed-step; no tolerances/Jacobian.
%
%   See also: apply_hh_jumps, apply_hhei_jumps, integrate_hh_hybrid.

    t  = tspan(:);
    nt = numel(t);
    y  = y0(:);
    Y  = zeros(nt, numel(y));
    Y(1, :) = y.';

    iV = 1:N;
    armed = y(iV) < V_th;

    for k = 1:nt-1
        h  = t(k+1) - t(k);
        tk = t(k);

        k1 = odefun(tk,       y);
        k2 = odefun(tk + h/2, y + (h/2) * k1);
        k3 = odefun(tk + h/2, y + (h/2) * k2);
        k4 = odefun(tk + h,   y + h * k3);

        Vprev = y(iV);
        y = y + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
        Vnew = y(iV);

        spiked = armed & (Vprev < V_th) & (Vnew >= V_th);
        armed(spiked)         = false;
        armed(Vnew < V_reset) = true;

        if any(spiked)
            y = jump_fn(y, spiked);
        end

        Y(k+1, :) = y.';
    end
end
