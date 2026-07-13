function [t, Y] = integrate_hh_hybrid(odefun, tspan, y0, opts, jp) %#ok<INUSL>
%INTEGRATE_HH_HYBRID  Event-aware fixed-step RK4 for the spiking HH network.
%
%   [t, Y] = INTEGRATE_HH_HYBRID(odefun, tspan, y0, opts, jp)
%
%   Integrates the HH network as a HYBRID system: RK4 on the smooth
%   between-spike flow (odefun), plus discrete Tsodyks-Markram / SFA / synaptic
%   jumps applied whenever a neuron's membrane potential crosses threshold
%   upward. The @ode45-compatible signature ([t,Y]=solver(rhs,tspan,y0,opts))
%   lets it slot into SRNNModelBase.run() and, crucially, into the Benettin
%   reshoot: the perturbed copy detects its OWN spike times, so spike-time
%   desynchronisation contributes to the divergence, giving a valid largest
%   Lyapunov exponent for the hybrid flow without saltation matrices.
%
%   opts (odeset options) is accepted for signature compatibility and IGNORED
%   (fixed-step; no tolerances/Jacobian). tspan must be the full uniform time
%   grid (numel >= 2); a 2-point [t0 tf] span is integrated as many steps only
%   if it has >=2 points, but callers should pass the native fs grid.
%
%   jp is the jump-parameter struct (built by SRNNModelHH.get_params, field
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
%   See also: ode_rk4 (FractionalResevoir), SRNNModelHH, SRNNModelBase.

    t  = tspan(:);
    nt = numel(t);
    y  = y0(:);
    Y  = zeros(nt, numel(y));
    Y(1, :) = y.';

    N = jp.N; K = jp.K;

    % --- state-block offsets (S = [V; m; h; n; a; b; p; g]) ----------------
    iV     = 1:N;
    base4  = 4 * N;
    len_a  = jp.has_a * jp.n_ad * jp.n_a;
    off_a  = base4;
    len_b  = jp.has_b * N * K;
    off_b  = off_a + len_a;
    len_p  = jp.has_p * N * K;
    off_p  = off_b + len_b;
    len_g  = N * K;
    off_g  = off_p + len_p;

    V_th = jp.V_th; V_reset = jp.V_reset;
    type_of = jp.type_of(:);
    Wabs = jp.Wabs;

    % A neuron is "armed" to spike once its V has fallen below the reset level.
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

        % Upward threshold crossing among armed neurons -> spikes this step.
        spiked = armed & (Vprev < V_th) & (Vnew >= V_th);
        armed(spiked)        = false;
        armed(Vnew < V_reset) = true;

        if any(spiked)
            sp_idx = find(spiked);

            % --- SFA: increment own adaptation state on own spike ----------
            if jp.has_a
                a = reshape(y(off_a + (1:len_a)), jp.n_ad, jp.n_a);
                sp_ad = spiked(jp.ad_idx);
                if any(sp_ad)
                    a(sp_ad, :) = a(sp_ad, :) + jp.a_incr;
                    y(off_a + (1:len_a)) = a(:);
                end
            end

            % --- release = p .* b (per pre-neuron j, post-type q) ----------
            if jp.has_b, b = reshape(y(off_b + (1:len_b)), N, K); else, b = ones(N, K); end
            if jp.has_p, p = reshape(y(off_p + (1:len_p)), N, K); else, p = jp.p0_mat; end
            rel = p .* b;

            % --- conductance bump at postsynaptic targets ------------------
            % g(:,P) gathers inputs from presynaptic type P; a spike of pre-
            % neuron j (type P) adds Wabs(j,i)*rel(j, type_of(i)) to target i.
            g = reshape(y(off_g + (1:len_g)), N, K);
            for P = 1:K
                pre = sp_idx(type_of(sp_idx) == P);
                if isempty(pre), continue; end
                Wp = Wabs(pre, :);                 % npre x N
                rel_exp = rel(pre, type_of);       % npre x N (target uses its post-type)
                g(:, P) = g(:, P) + sum(Wp .* rel_exp, 1).';
            end
            y(off_g + (1:len_g)) = g(:);

            % --- STD depression (deplete available resource) ---------------
            if jp.has_b
                b(spiked, :) = max(b(spiked, :) - rel(spiked, :), 0);
                y(off_b + (1:len_b)) = b(:);
            end

            % --- STF facilitation (raise release prob toward 1) ------------
            if jp.has_p
                krow = jp.kappa(type_of(spiked), :);
                p(spiked, :) = min(p(spiked, :) + krow .* (1 - p(spiked, :)), 1);
                y(off_p + (1:len_p)) = p(:);
            end
        end

        Y(k+1, :) = y.';
    end
end
