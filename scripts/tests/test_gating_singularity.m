% test_gating_singularity.m
% Exercises hh_gating_rates: finiteness across a full voltage sweep, and that
% the three singularity-patched rates equal their L'Hopital limits at the
% singular voltages. Prints HH_GATING_PASS / HH_GATING_FAIL.
setup_paths();
ok = true;
report = @(name, pass) fprintf('  [%s] %s\n', tern(pass, 'ok', 'FAIL'), name);

% Full sweep: no NaN/Inf anywhere (the whole point of the patches).
Vsweep = (-100:0.001:60)';
[am, bm, ah, bh, an, bn] = hh_gating_rates(Vsweep);
finite_ok = all(isfinite([am; bm; ah; bh; an; bn]));
ok = ok && finite_ok; report('finite over full sweep', finite_ok);

% Exact limit values at the singular voltages.
[am0, ~, ~, ~, ~, ~] = hh_gating_rates(-54);
[~, bm0, ~, ~, ~, ~] = hh_gating_rates(-27);
[~, ~, ~, ~, an0, ~] = hh_gating_rates(-52);
lim_ok = abs(am0 - 1.28) < 1e-9 && abs(bm0 - 1.40) < 1e-9 && abs(an0 - 0.16) < 1e-9;
ok = ok && lim_ok; report('L''Hopital limits at singularities', lim_ok);

% Continuity: just off the singularity should be close to the limit.
[amn, ~, ~, ~, ~, ~] = hh_gating_rates(-54 + 1e-6);
[~, ~, ~, ~, ann, ~] = hh_gating_rates(-52 - 1e-6);
cont_ok = abs(amn - 1.28) < 1e-3 && abs(ann - 0.16) < 1e-3;
ok = ok && cont_ok; report('continuity near singularities', cont_ok);

% Steady-state gating in [0,1].
[minf, hinf, ninf] = hh_gating_inf(Vsweep);
ss_ok = all(minf >= 0 & minf <= 1 & hinf >= 0 & hinf <= 1 & ninf >= 0 & ninf <= 1);
ok = ok && ss_ok; report('steady-state gating in [0,1]', ss_ok);

if ok, disp('HH_GATING_PASS'); else, disp('HH_GATING_FAIL'); end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end
