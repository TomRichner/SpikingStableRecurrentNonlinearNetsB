% test_hh_run.m
% A small Campagnola-based HH network builds and runs end to end with SFA+STD+STF
% enabled: finite trajectory, some spiking, and the mechanism toggles change the
% state dimension as expected. Prints HH_RUN_PASS / HH_RUN_FAIL.
setup_paths();
ok = true;
report = @(name, pass) fprintf('  [%s] %s\n', tern(pass, 'ok', 'FAIL'), name);

m = SRNNModelHH('n', 40, 'indegree', 8, 'T_range', [0 400], ...
                'lya_method', 'none', 'store_full_state', true, ...
                'n_a', 1, 'n_b', 1, 'n_u', 1);
m.build(); m.run();

fin_ok = all(isfinite(m.S_out(:)));
ok = ok && fin_ok; report('trajectory finite', fin_ok);

n_sp = size(m.plot_data.spikes, 1);
spk_ok = n_sp > 0;
ok = ok && spk_ok; report(sprintf('network spikes (%d)', n_sp), spk_ok);

% State dimension bookkeeping: N_sys_eqs matches actual columns.
p = m.cached_params;
dim_ok = size(m.S_out, 2) == p.N_sys_eqs;
ok = ok && dim_ok; report('N_sys_eqs matches state width', dim_ok);

% Disabling mechanisms removes their state blocks.
m2 = SRNNModelHH('n', 40, 'indegree', 8, 'T_range', [0 50], 'lya_method', 'none', ...
                 'n_a', 0, 'n_b', 0, 'n_u', 0, 'store_full_state', true);
m2.build();
expected2 = 4*40 + 40*4;   % V,m,h,ng + g only
dim2_ok = m2.cached_params.N_sys_eqs == expected2;
ok = ok && dim2_ok; report('mechanism toggles drop state blocks', dim2_ok);

% b (STD availability) stays in [0,1]; p (STF) stays in [0,1].
st = SRNNModelHH.unpack_states_hh(m.S_out, m.cached_params);
bnd_ok = all(st.b(:) >= -1e-9 & st.b(:) <= 1 + 1e-9) && ...
         all(st.p(:) >= -1e-9 & st.p(:) <= 1 + 1e-9);
ok = ok && bnd_ok; report('b,p within [0,1]', bnd_ok);

try
    m.plot(); close(gcf); plot_ok = true;
catch ME
    plot_ok = false; fprintf('    plot error: %s\n', ME.message);
end
ok = ok && plot_ok; report('plot() runs', plot_ok);

if ok, disp('HH_RUN_PASS'); else, disp('HH_RUN_FAIL'); end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end
