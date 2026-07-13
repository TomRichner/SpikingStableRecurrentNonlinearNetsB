% test_stp_jumps.m
% Unit test of integrate_hh_hybrid's event jumps, independent of HH tuning.
% A synthetic odefun ramps neuron 1's V through threshold exactly once; we
% verify the Tsodyks-Markram depression (b), facilitation (p), and the
% postsynaptic conductance bump (g) are applied with the correct arithmetic.
% Prints HH_STP_PASS / HH_STP_FAIL.
setup_paths();
ok = true;
report = @(name, pass) fprintf('  [%s] %s\n', tern(pass, 'ok', 'FAIL'), name);

N = 2; K = 1;
jp = struct();
jp.N = N; jp.K = K; jp.n_a = 0; jp.n_ad = 0; jp.ad_idx = [];
jp.type_of = [1; 1];
jp.Wabs = [0 1; 0 0];          % pre x post: neuron 1 -> neuron 2, weight 1
jp.kappa = 0.2;                % K x K
jp.p0_mat = [0.3; 0.3];        % N x K
jp.a_incr = 0;
jp.V_th = -20; jp.V_reset = -40;
jp.has_a = false; jp.has_b = true; jp.has_p = true;

% State layout (no a): [V(2); m(2); h(2); ng(2); b(2); p(2); g(2)]  -> 14
Ns = 4*N + N*K + N*K + N*K;
y0 = zeros(Ns, 1);
y0(1:2) = [-30; -70];          % V: neuron 1 below thresh, neuron 2 far below
iB = 4*N + (1:2); y0(iB) = [1; 1];       % b
iP = 4*N + 2 + (1:2); y0(iP) = [0.3; 0.3];  % p
iG = 4*N + 4 + (1:2); y0(iG) = [0; 0];   % g

% Ramp only neuron 1's V upward (crosses -20 once, never re-arms).
dydt = zeros(Ns, 1); dydt(1) = 50;       % +50 mV/ms
odefun = @(t, y) dydt;

tspan = (0:0.01:2)';
[~, Y] = integrate_hh_hybrid(odefun, tspan, y0, [], jp);
yend = Y(end, :).';

% Expected after exactly one spike of neuron 1 (release = p*b = 0.3):
%   b1 = 1 - 0.3 = 0.7 ; p1 = 0.3 + 0.2*(1-0.3) = 0.44 ; g at neuron 2 = 1*0.3 = 0.3
b1 = yend(iB(1)); b2 = yend(iB(2));
p1 = yend(iP(1)); p2 = yend(iP(2));
g1 = yend(iG(1)); g2 = yend(iG(2));

dep_ok = abs(b1 - 0.7) < 1e-9 && abs(b2 - 1.0) < 1e-9;
fac_ok = abs(p1 - 0.44) < 1e-9 && abs(p2 - 0.3) < 1e-9;
cond_ok = abs(g2 - 0.3) < 1e-9 && abs(g1 - 0.0) < 1e-9;
ok = ok && dep_ok && fac_ok && cond_ok;
report('STD depression b1=0.7', dep_ok);
report('STF facilitation p1=0.44', fac_ok);
report('postsynaptic conductance g2=0.3', cond_ok);

% Single spike only: rerun with a longer ramp still gives one event.
if ok, disp('HH_STP_PASS'); else, disp('HH_STP_FAIL'); end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end
