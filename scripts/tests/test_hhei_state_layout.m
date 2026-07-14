% test_hhei_state_layout.m
% Per-population ragged state layout for SRNNModelHHEI: N_sys_eqs matches the
% actual state width across several n_a_vec/n_b_vec, and blocks unpack correctly.
% Prints HHEI_LAYOUT_PASS / HHEI_LAYOUT_FAIL.
setup_paths();
ok = true;
report = @(name, pass) fprintf('  [%s] %s\n', tern(pass, 'ok', 'FAIL'), name);

configs = {[3 0] [1 0]; [3 2] [1 1]; [0 0] [0 0]; [2 0] [0 0]};
for c = 1:size(configs, 1)
    na = configs{c, 1}; nb = configs{c, 2};
    m = SRNNModelHHEI('n', 50, 'indegree', 10, 'lya_method', 'none', ...
                      'n_a_vec', na, 'n_b_vec', nb, 'store_full_state', true, ...
                      'T_range', [0 30]);
    m.build();
    p = m.cached_params; L = p.layout;
    nE = nnz(m.exc_type(m.type_of)); nI = m.n - nE;
    exp_a = nE*na(1) + nI*na(2);
    exp_b = nE*nb(1) + nI*nb(2);
    exp_sys = 4*m.n + exp_a + exp_b + m.n*m.n_types;
    m.run();
    width_ok = size(m.S_out, 2) == p.N_sys_eqs && p.N_sys_eqs == exp_sys && ...
               L.len_a == exp_a && L.len_b == exp_b;
    ok = ok && width_ok;
    report(sprintf('na=[%s] nb=[%s]: N_sys=%d (exp %d), width match', ...
        num2str(na), num2str(nb), p.N_sys_eqs, exp_sys), width_ok);
end

% Unpack round-trip: a_sum and b_prod have the right shape and sane ranges.
m = SRNNModelHHEI('n', 40, 'indegree', 10, 'lya_method', 'none', ...
                  'n_a_vec', [3 0], 'n_b_vec', [1 0], 'store_full_state', true, ...
                  'input_config', struct('bias', 6, 'drive_types', [], 'drive_window', [], 'drive_amp', 0), ...
                  'T_range', [0 150]);
m.build(); m.run();
st = SRNNModelHHEI.unpack_states_hhei(m.S_out, m.cached_params);
unpack_ok = isequal(size(st.a_sum), [40 size(m.S_out,1)]) && ...
            all(st.b_prod(:) >= -1e-9 & st.b_prod(:) <= 1 + 1e-9) && ...
            all(st.a_sum(~m.exc_type(m.type_of), :) == 0, 'all');   % inh carries no SFA
ok = ok && unpack_ok; report('unpack shapes + inh has no SFA + b_prod in [0,1]', unpack_ok);

if ok, disp('HHEI_LAYOUT_PASS'); else, disp('HHEI_LAYOUT_FAIL'); end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end
