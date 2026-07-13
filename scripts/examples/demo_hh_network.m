% demo_hh_network.m
% Build, run, and plot a small spiking HH network with SFA + STD + STF and a
% Benettin largest-Lyapunov-exponent estimate. A minimal runnable example of
% the SRNNModelHH workflow.
setup_paths();

model = SRNNModelHH( ...
    'n', 300, ...
    'T_range', [0 1000], ...       % ms
    'level_of_chaos', 1.5, ...     % edge-of-chaos scan knob on the weight matrix
    'n_a', 1, 'n_b', 1, 'n_u', 1, ...% SFA + STD + STF all on
    'lya_method', 'benettin', ...
    'lya_dt', 20, ...
    'store_full_state', true);

model.build();
model.run();
model.plot();

fprintf('\nLargest Lyapunov exponent: %.4g /ms\n', model.lya_results.LLE);
