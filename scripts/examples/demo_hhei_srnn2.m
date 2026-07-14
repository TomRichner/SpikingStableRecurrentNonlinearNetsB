% demo_hhei_srnn2.m
% Spiking HH E/I network configured to parallel the rate model
% SRNNModel2('n_a_E',3,'n_b_E',1): two populations (E/I) with RMT connectivity,
% 3 logspaced SFA timescales on E, single-timescale STD on E, no STF. Build, run,
% plot, and report the Benettin LLE (in 1/s for comparison with the rate model).
%
% NOTE: SFA timescales match SRNNModel2 (0.25-10 s), but the run is trimmed to
% 10 s with a 5 s Benettin transient for speed. CAVEAT: 10 s < a few x the 10 s
% SFA pool, so the slowest pool is NOT fully settled -- the LLE here is
% indicative, not a converged attractor value (lengthen T_range / lya_transient
% for a rigorous number; see CLAUDE.md on run-length cost). The operating rate /
% p0 / c / g_syn_scale are hand-tuning knobs; values here are illustrative.
setup_paths();

model = SRNNModelHHEI( ...
    'n', 300, ...
    'indegree', 100, ...              % alpha = 1/3 (as SRNNModel2)
    'type_fractions', [0.5 0.5], ...  % E/I, exc first
    'T_range', [0 3000], ...         % ms (10 s run)
    'level_of_chaos', 10.0, ...        % edge-of-chaos scan knob
    'n_a_vec', [3 0], ...             % 3 SFA timescales on E, none on I
    'n_b_vec', [1 0], ...             % single STD on E, none on I; STF absent
    'tau_a_range', [250 10000], ...   % ms = SRNNModel2's SFA range (0.25-10 s), 3 logspaced pools
    'tau_rec_default', 1000, ...      % ms (= SRNNModel2 tau_rec = 1 s)
    'c_type', [0.15 0], ...           % SFA strength (c_eff = 0.15/3 = 0.05, as SRNNModel2 c_E)
    'p0_type', [0.25 0], ...         % release/depletion ~ tau_rel=0.25 s at 5 Hz
    'target_rate', 1, ...             % Hz (sets DC-gain SFA increment)
    'lya_method', 'benettin', 'lya_dt', 20, 'lya_transient', 1500, ...
    'store_full_state', true);

model.input_config.bias = 0.35;          % tonic drive so the network is active (hand-tune)

model.build();
model.run();
model.plot();

fprintf('\nLargest Lyapunov exponent: %.4g /ms  (= %.4g /s)\n', ...
    model.lya_results.LLE, model.lya_results.LLE_per_s);
