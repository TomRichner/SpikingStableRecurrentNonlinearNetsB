% test_benettin.m
% Benettin largest-Lyapunov-exponent sanity checks for the spiking network:
%   (a) a finite LLE is produced and the pipeline is deterministic across
%       identical constructions (fixed seeds);
%   (b) increasing recurrent coupling raises the LLE (weak/uncoupled network is
%       not more chaotic than a strongly recurrent one).
% Prints HH_BENETTIN_PASS / HH_BENETTIN_FAIL.
setup_paths();
ok = true;
report = @(name, pass) fprintf('  [%s] %s\n', tern(pass, 'ok', 'FAIL'), name);

common = {'n', 40, 'indegree', 8, 'T_range', [0 600], 'lya_method', 'benettin', ...
          'lya_dt', 20, 'n_a', 1, 'n_b', 1, 'n_u', 1};

% (a) determinism: same params/seeds -> identical LLE.
mA = SRNNModelHH(common{:}, 'level_of_chaos', 1.0); mA.build(); mA.run();
mB = SRNNModelHH(common{:}, 'level_of_chaos', 1.0); mB.build(); mB.run();
det_ok = isfinite(mA.lya_results.LLE) && ...
         abs(mA.lya_results.LLE - mB.lya_results.LLE) < 1e-9;
ok = ok && det_ok;
report(sprintf('deterministic finite LLE (%.4g /ms)', mA.lya_results.LLE), det_ok);

% (b) monotone trend: weak coupling LLE <= strong coupling LLE.
mWeak   = SRNNModelHH(common{:}, 'level_of_chaos', 0.0); mWeak.build(); mWeak.run();
mStrong = SRNNModelHH(common{:}, 'level_of_chaos', 3.0); mStrong.build(); mStrong.run();
lle_weak = mWeak.lya_results.LLE;
lle_strong = mStrong.lya_results.LLE;
trend_ok = isfinite(lle_weak) && isfinite(lle_strong) && (lle_strong >= lle_weak - 1e-6);
ok = ok && trend_ok;
report(sprintf('LLE(weak)=%.4g <= LLE(strong)=%.4g', lle_weak, lle_strong), trend_ok);

if ok, disp('HH_BENETTIN_PASS'); else, disp('HH_BENETTIN_FAIL'); end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end
