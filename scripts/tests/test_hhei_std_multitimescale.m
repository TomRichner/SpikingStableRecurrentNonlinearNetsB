% test_hhei_std_multitimescale.m
% Multi-timescale STD jump arithmetic for SRNNModelHHEI: efficacy = product over
% a presynaptic neuron's pools, conductance bump = Wabs*efficacy, and each pool
% depletes b_m -= p0*b_m with a distinct per-pool recovery tau. Deterministic
% direct test of apply_hhei_jumps. Prints HHEI_STD_PASS / HHEI_STD_FAIL.
setup_paths();
ok = true;
report = @(name, pass) fprintf('  [%s] %s\n', tern(pass, 'ok', 'FAIL'), name);

m = SRNNModelHHEI('n', 2, 'type_names', {'exc'}, 'type_fractions', 1, ...
    'exc_type_names', {'exc'}, 'indegree', 2, 'lya_method', 'none', ...
    'n_a_vec', 0, 'n_b_vec', 2, 'p0_type', 0.2, 'tau_rec_cell', {[100 1000]});
m.build();
p = m.cached_params; L = p.layout; Wabs = p.Wabs;

% State layout: [V(2);m(2);h(2);ng(2); b(2x2); g(2x1)]  -> off_b=8, off_g=12
y = zeros(p.N_sys_eqs, 1);
% b: neuron 1 = [0.8 0.6], neuron 2 = [1 1]; column-major reshape(2,2)
y(L.off_b + 1) = 0.8; y(L.off_b + 2) = 1.0;   % pool 1: b(1,1)=0.8, b(2,1)=1
y(L.off_b + 3) = 0.6; y(L.off_b + 4) = 1.0;   % pool 2: b(1,2)=0.6, b(2,2)=1

spiked = [true; false];
y2 = SRNNModelHHEI.apply_hhei_jumps(y, spiked, p);

p0 = 0.2; eff1 = 0.8 * 0.6;                    % efficacy = product (pre-depletion)
% depletion: b(1,:) *= (1-p0)
dep_ok = abs(y2(L.off_b + 1) - 0.8*(1-p0)) < 1e-12 && ...
         abs(y2(L.off_b + 3) - 0.6*(1-p0)) < 1e-12 && ...
         y2(L.off_b + 2) == 1 && y2(L.off_b + 4) == 1;
ok = ok && dep_ok; report('per-pool depletion b_m -= p0*b_m (neuron 2 untouched)', dep_ok);

% conductance bump at neuron 2 from neuron 1 = Wabs(1,2)*efficacy; neuron 1 gets none
g_ok = abs(y2(L.off_g + 2) - Wabs(1,2)*eff1) < 1e-12 && abs(y2(L.off_g + 1) - 0) < 1e-12;
ok = ok && g_ok; report(sprintf('conductance bump g(2)=Wabs(1,2)*Pi b = %.4f', Wabs(1,2)*eff1), g_ok);

% distinct per-pool recovery timescales wired through
tau_ok = isequal(L.std.tau_rec, [100 1000]);
ok = ok && tau_ok; report('distinct per-pool tau_rec [100 1000] ms', tau_ok);

if ok, disp('HHEI_STD_PASS'); else, disp('HHEI_STD_FAIL'); end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end
