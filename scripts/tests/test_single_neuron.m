% test_single_neuron.m
% A single isolated HH neuron: rests quietly with no drive, and fires a spike
% train under supra-threshold constant current. Prints HH_SINGLE_PASS/FAIL.
setup_paths();
ok = true;
report = @(name, pass) fprintf('  [%s] %s\n', tern(pass, 'ok', 'FAIL'), name);

base = {'type_names', {'pyr'}, 'type_fractions', 1, 'exc_type_names', {'pyr'}, ...
        'use_campagnola_data', false, 'n', 1, 'indegree', 0, ...
        'n_a', 0, 'n_b', 0, 'n_u', 0, 'lya_method', 'none', ...
        'T_range', [0 200], 'store_full_state', true};

% (1) Resting: no external current -> no spikes, V stays bounded near rest.
m0 = SRNNModelHH(base{:});
m0.input_config.amp = 0; m0.input_config.bias = 0;
m0.build(); m0.run();
n_rest = size(m0.plot_data.spikes, 1);
rest_ok = (n_rest == 0) && all(isfinite(m0.plot_data.V(:))) && all(m0.plot_data.V(:) < 0);
ok = ok && rest_ok; report(sprintf('rests quietly (%d spikes)', n_rest), rest_ok);

% (2) Driven: supra-threshold constant current -> repetitive spiking.
m1 = SRNNModelHH(base{:});
m1.input_config.amp = 0; m1.input_config.bias = 12;   % uA/cm^2
m1.build(); m1.run();
n_sp = size(m1.plot_data.spikes, 1);
drive_ok = (n_sp >= 3) && all(isfinite(m1.plot_data.V(:)));
ok = ok && drive_ok; report(sprintf('fires under drive (%d spikes)', n_sp), drive_ok);

% (3) Peak amplitude is spike-like (overshoots 0 mV).
peak_ok = max(m1.plot_data.V(:)) > 0;
ok = ok && peak_ok; report('action potential overshoots 0 mV', peak_ok);

if ok, disp('HH_SINGLE_PASS'); else, disp('HH_SINGLE_FAIL'); end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end
