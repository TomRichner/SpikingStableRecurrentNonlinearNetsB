function y = apply_hh_jumps(y, spiked, jp)
%APPLY_HH_JUMPS  Tsodyks-Markram / SFA / conductance jumps for SRNNModelHH.
%
%   y = APPLY_HH_JUMPS(y, spiked, jp) applies, in place, the discrete jumps
%   triggered by the neurons flagged in the logical vector `spiked`, for the
%   per-(pre-neuron, post-type) STD/STF layout of SRNNModelHH:
%     State S = [V; m; h; n; a(n_ad*n_a); b(N*K); p(N*K); g(N*K)]
%     - SFA: own-spike increment a += a_incr (adapting neurons only)
%     - release = p.*b (per pre-neuron j, post-type q)
%     - conductance bump g(i,P) += Wabs(j,i)*release(j, type_of(i))
%     - STD depression b -= release ; STF facilitation p += kappa*(1-p)
%   jp is the jump-parameter struct (SRNNModelHH.get_params .jump). See
%   integrate_hh_hybrid for field docs.
%
%   Factored out of integrate_hh_hybrid so the generic event integrator
%   (integrate_hh_events) can drive either synaptic layout via a jump closure.
%
%   See also: integrate_hh_events, integrate_hh_hybrid, apply_hhei_jumps.

    N = jp.N; K = jp.K;
    base4 = 4 * N;
    len_a = jp.has_a * jp.n_ad * jp.n_a; off_a = base4;
    len_b = jp.has_b * N * K;            off_b = off_a + len_a;
    len_p = jp.has_p * N * K;            off_p = off_b + len_b;
    len_g = N * K;                       off_g = off_p + len_p;

    type_of = jp.type_of(:);
    Wabs = jp.Wabs;
    sp_idx = find(spiked);
    if isempty(sp_idx), return; end

    % --- SFA: increment own adaptation state on own spike ------------------
    if jp.has_a
        a = reshape(y(off_a + (1:len_a)), jp.n_ad, jp.n_a);
        sp_ad = spiked(jp.ad_idx);
        if any(sp_ad)
            a(sp_ad, :) = a(sp_ad, :) + jp.a_incr;
            y(off_a + (1:len_a)) = a(:);
        end
    end

    % --- release = p .* b (per pre-neuron j, post-type q) ------------------
    if jp.has_b, b = reshape(y(off_b + (1:len_b)), N, K); else, b = ones(N, K); end
    if jp.has_p, p = reshape(y(off_p + (1:len_p)), N, K); else, p = jp.p0_mat; end
    rel = p .* b;

    % --- conductance bump at postsynaptic targets -------------------------
    g = reshape(y(off_g + (1:len_g)), N, K);
    for P = 1:K
        pre = sp_idx(type_of(sp_idx) == P);
        if isempty(pre), continue; end
        Wp = Wabs(pre, :);                 % npre x N
        rel_exp = rel(pre, type_of);       % npre x N (target uses its post-type)
        g(:, P) = g(:, P) + sum(Wp .* rel_exp, 1).';
    end
    y(off_g + (1:len_g)) = g(:);

    % --- STD depression ---------------------------------------------------
    if jp.has_b
        b(spiked, :) = max(b(spiked, :) - rel(spiked, :), 0);
        y(off_b + (1:len_b)) = b(:);
    end

    % --- STF facilitation -------------------------------------------------
    if jp.has_p
        krow = jp.kappa(type_of(spiked), :);
        p(spiked, :) = min(p(spiked, :) + krow .* (1 - p(spiked, :)), 1);
        y(off_p + (1:len_p)) = p(:);
    end
end
