% test_hhei_connectivity.m
% RMT E/I connectivity for SRNNModelHHEI: magnitude weights + Dale-by-reversal,
% wrong-sign entries zeroed, indegree ~ alpha*n, deterministic, and the build
% transform reproduces exactly. Prints HHEI_CONN_PASS / HHEI_CONN_FAIL.
setup_paths();
ok = true;
report = @(name, pass) fprintf('  [%s] %s\n', tern(pass, 'ok', 'FAIL'), name);

n = 60; indeg = 20;
m = SRNNModelHHEI('n', n, 'indegree', indeg, 'lya_method', 'none');
m.build();
p = m.cached_params; W = m.W;

inv_ok = all(W(:) >= 0) && all(diag(W) == 0);
ok = ok && inv_ok; report('W >= 0 (magnitude), zero diagonal', inv_ok);

alpha = indeg / n; expd = alpha * n;
meanindeg = mean(sum(W > 0, 1));
indeg_ok = abs(meanindeg - expd) < 0.15 * expd;
ok = ok && indeg_ok; report(sprintf('mean in-degree %.1f ~ alpha*n=%.1f', meanindeg, expd), indeg_ok);

excT = find(m.exc_type); inhT = find(~m.exc_type);
dale_ok = all(p.E_syn_vec(excT) == m.E_exc) && all(p.E_syn_vec(inhT) == m.E_inh);
ok = ok && dale_ok; report('reversal potentials Dale-mapped (exc->E_exc, inh->E_inh)', dale_ok);

m2 = SRNNModelHHEI('n', n, 'indegree', indeg, 'lya_method', 'none'); m2.build();
det_ok = isequal(m.W, m2.W);
ok = ok && det_ok; report('deterministic W across builds (fixed seed)', det_ok);

% Reproduce the build transform (validates wrong-sign zeroing).
rng(m.rng_seeds(1));
F = 1 / sqrt(n * alpha * (2 - alpha));
rmt = RMTMatrix(n);
rmt.f = nnz(m.is_exc) / n;
rmt.mu_tilde_e = 3*F; rmt.mu_tilde_i = -4*F;
rmt.sigma_tilde_e = F; rmt.sigma_tilde_i = F; rmt.zrs_mode = 'none';
rmt.alpha = alpha;
Wpre = rmt.W.';
exc = m.is_exc;
Wpre(exc, :)  = max(Wpre(exc, :), 0);
Wpre(~exc, :) = min(Wpre(~exc, :), 0);
Wexp = abs(Wpre); Wexp(1:n+1:end) = 0;
Wexp = m.level_of_chaos * m.g_syn_scale * Wexp;
repro_ok = max(abs(Wexp(:) - m.W(:))) < 1e-12;
ok = ok && repro_ok; report('build transform reproduces (RMT + zero wrong-sign)', repro_ok);

% Tilde defaults materialized as 3F/-4F/F.
tilde_ok = abs(m.mu_E_tilde - 3*F) < 1e-12 && abs(m.mu_I_tilde + 4*F) < 1e-12 && ...
           abs(m.sigma_E_tilde - F) < 1e-12;
ok = ok && tilde_ok; report('RMT tilde defaults 3F/-4F/F', tilde_ok);

if ok, disp('HHEI_CONN_PASS'); else, disp('HHEI_CONN_FAIL'); end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end
