% test_hhei_run.m
% End-to-end: a small E/I RMT HH network with 3-timescale SFA + STD builds, runs,
% and produces a finite Benettin LLE. Prints HHEI_RUN_PASS / HHEI_RUN_FAIL.
setup_paths();
ok = true;
report = @(name, pass) fprintf('  [%s] %s\n', tern(pass, 'ok', 'FAIL'), name);

m = SRNNModelHHEI('n', 40, 'indegree', 12, 'T_range', [0 600], ...
    'level_of_chaos', 1.5, 'n_a_vec', [3 0], 'n_b_vec', [1 0], ...
    'tau_a_range', [50 500], 'lya_method', 'benettin', 'lya_dt', 20, 'lya_transient', 200, ...
    'input_config', struct('bias', 8, 'drive_types', [], 'drive_window', [], 'drive_amp', 0), ...
    'store_full_state', true);
m.build(); m.run();

fin_ok = all(isfinite(m.S_out(:)));
ok = ok && fin_ok; report('trajectory finite', fin_ok);

n_sp = size(m.plot_data.spikes, 1);
spk_ok = n_sp > 0;
ok = ok && spk_ok; report(sprintf('network spikes (%d)', n_sp), spk_ok);

lle_ok = isfinite(m.lya_results.LLE) && isfinite(m.lya_results.LLE_per_s);
ok = ok && lle_ok; report(sprintf('finite LLE = %.4g /ms (= %.4g /s)', ...
    m.lya_results.LLE, m.lya_results.LLE_per_s), lle_ok);

try
    m.plot(); close(gcf); plot_ok = true;
catch ME
    plot_ok = false; fprintf('    plot error: %s\n', ME.message);
end
ok = ok && plot_ok; report('plot() runs', plot_ok);

if ok, disp('HHEI_RUN_PASS'); else, disp('HHEI_RUN_FAIL'); end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end
