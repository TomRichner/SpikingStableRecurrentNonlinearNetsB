% test_hhei_sfa_timescales.m
% DC-gain balance of the multi-timescale SFA K-current: with increment
% Delta_l = a_incr0/tau_l, each logspaced pool's steady-state mean <a_l> is equal
% (independent of tau_l), matching the rate model's a_l* = r. Verified on a single
% regularly-firing neuron with SFA feedback disabled (c=0) so it fires steadily.
% Prints HHEI_SFA_PASS / HHEI_SFA_FAIL.
setup_paths();
ok = true;
report = @(name, pass) fprintf('  [%s] %s\n', tern(pass, 'ok', 'FAIL'), name);

m = SRNNModelHHEI('n', 1, 'type_names', {'exc'}, 'type_fractions', 1, ...
    'exc_type_names', {'exc'}, 'indegree', 1, 'lya_method', 'none', ...
    'n_a_vec', 3, 'n_b_vec', 0, 'c_type', 0, ...        % c=0: pools recorded, no feedback
    'tau_a_range', [50 500], 'store_full_state', true, ...
    'input_config', struct('bias', 12, 'drive_types', [], 'drive_window', [], 'drive_amp', 0), ...
    'T_range', [0 1500]);
m.build(); m.run();

L = m.cached_params.layout;
assert(numel(L.sfa) == 1 && L.sfa.n_a == 3, 'expected one 3-pool SFA block');
tau = L.sfa.tau;                         % [50 .. 500]
a_pools = m.S_out(:, L.sfa.off + (1:3)); % columns = pools (fast..slow)

n_sp = size(m.plot_data.spikes, 1);
fires_ok = n_sp >= 10;
ok = ok && fires_ok; report(sprintf('neuron fires steadily (%d spikes)', n_sp), fires_ok);

% Steady-state pool means over the late window (after ~3x slowest tau).
late = m.t_out >= 1000;
mean_a = mean(a_pools(late, :), 1);
balance = max(mean_a) / min(mean_a);
bal_ok = balance < 1.20;                 % pools within 20% -> DC-gain balanced
ok = ok && bal_ok;
report(sprintf('pool means balanced: [%.3f %.3f %.3f], ratio=%.2f', mean_a, balance), bal_ok);

% Timescales are actually distinct and increasing.
tau_ok = tau(1) < tau(2) && tau(2) < tau(3) && abs(tau(1)-50) < 1e-6 && abs(tau(3)-500) < 1e-6;
ok = ok && tau_ok; report(sprintf('distinct logspaced tau [%.0f %.0f %.0f] ms', tau), tau_ok);

if ok, disp('HHEI_SFA_PASS'); else, disp('HHEI_SFA_FAIL'); end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end
