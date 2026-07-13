% demo_canonical_4neuron.m
% Minimal, fully transparent circuit for visually inspecting SFA, STD, and STF:
% one neuron of each cell type (pyr, pvalb, sst, vip) wired with the CANONICAL
% Campagnola connectivity (every ordered inter-type pair present once, at its
% mean |psp| weight). Protocol: settle 250 ms -> drive the pyr neuron for 250 ms
% -> settle 250 ms. Then plot, per neuron and per connection, the mechanism state
% variables with presynaptic spike times marked so each jump is visible.
setup_paths();

% --- build the canonical 4-neuron network ------------------------------------
model = SRNNModelHH( ...
    'n', 4, ...
    'type_fractions', [1 1 1 1], ...       % exactly one neuron per type
    'use_campagnola_data', true, ...       % canonical mean weights from Campagnola
    'conn_prob', ones(4), ...              % deterministic full inter-type connectivity
    'w_cv', 0, ...                         % no weight heterogeneity
    'n_a', 1, 'n_b', 1, 'n_u', 1, ...      % SFA + STD + STF all on
    'T_range', [0 750], ...                % ms
    'lya_method', 'none', ...              % inspection demo -- skip Lyapunov
    'g_syn_scale', 0.30, ...               % a bit stronger so pyr can drive its targets
    'store_full_state', true);

pyr = find(strcmpi(model.type_names, 'pyr'));   % cell-type index to drive
model.input_config.drive_types  = pyr;
model.input_config.drive_window = [250 500];    % ms
model.input_config.drive_amp    = 12;           % uA/cm^2 (supra-threshold for pyr)
model.input_config.bias         = 0;

model.build();
model.run();

% --- full-resolution state + spike raster ------------------------------------
st = SRNNModelHH.unpack_states_hh(model.S_out, model.cached_params);   % V, a, b, p, g
spikes = SRNNModelHH.detect_spikes(model.t_out, model.S_out(:, 1:model.n), ...
                                   model.V_th, model.V_reset);
t = model.t_out;
type_of = model.type_of; names = model.type_names; K = model.n_types;
cmap = lines(K);
win = model.input_config.drive_window;
spk_times = unique(spikes(:, 1));                 % all spike times (few neurons fire)

% --- figure ------------------------------------------------------------------
fig = figure('Color', 'w', 'Name', 'Canonical 4-neuron: SFA / STD / STF', ...
             'Position', [80 80 900 950]);
tl = tiledlayout(fig, 5, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = gobjects(5, 1);

% (1) membrane potential ------------------------------------------------------
ax(1) = nexttile(tl); hold(ax(1), 'on');
for i = 1:model.n
    plot(ax(1), t, st.V(i, :), 'Color', cmap(type_of(i), :), 'LineWidth', 1);
end
shade_window(ax(1), win); ylabel(ax(1), 'V (mV)');
legend(ax(1), names, 'Location', 'eastoutside'); title(ax(1), 'membrane potential');

% (2) SFA a (per neuron; sum over timescales) ---------------------------------
ax(2) = nexttile(tl); hold(ax(2), 'on');
aa = reshape(sum(st.a, 2), model.n, []);          % n x nt
h2 = gobjects(model.n, 1); lab2 = {};
for i = 1:model.n
    h2(i) = plot(ax(2), t, aa(i, :), 'Color', cmap(type_of(i), :), 'LineWidth', 1.25);
    lab2{end+1} = names{type_of(i)}; %#ok<AGROW>
end
mark_spikes(ax(2), spk_times); shade_window(ax(2), win);
ylabel(ax(2), 'SFA  a'); legend(ax(2), h2, lab2, 'Location', 'eastoutside');
title(ax(2), 'spike-frequency adaptation (own-spike driven)');

% Which presynaptic neurons actually fired -> only plot their (active) synapses.
active_pre = unique(spikes(:, 2))';

% (3) STD b_{j,q}  and (4) STF p_{j,q}  for active pre neurons ----------------
ax(3) = nexttile(tl); hold(ax(3), 'on');
plot_stp_traces(ax(3), t, st.b, active_pre, type_of, names, cmap);
mark_spikes(ax(3), spk_times); shade_window(ax(3), win);
ylabel(ax(3), 'STD  b'); ylim(ax(3), [0 1.05]);
title(ax(3), 'short-term depression (available resource; drops at each pre spike)');

ax(4) = nexttile(tl); hold(ax(4), 'on');
plot_stp_traces(ax(4), t, st.p, active_pre, type_of, names, cmap);
mark_spikes(ax(4), spk_times); shade_window(ax(4), win);
ylabel(ax(4), 'STF  p'); ylim(ax(4), [0 1.05]);
title(ax(4), 'short-term facilitation (release prob; jumps up at each pre spike)');

% (5) synaptic conductance g_{i,P} (post-neuron i from pre-type P) -------------
ax(5) = nexttile(tl); hold(ax(5), 'on');
h5 = gobjects(0); lab5 = {};
for i = 1:model.n
    for P = 1:K
        gi = squeeze(st.g(i, P, :))';
        if max(gi) < 1e-6, continue; end          % skip inactive synapse channels
        h5(end+1) = plot(ax(5), t, gi, 'LineWidth', 1); %#ok<AGROW>
        lab5{end+1} = sprintf('%s <- %s', names{type_of(i)}, names{P}); %#ok<AGROW>
    end
end
mark_spikes(ax(5), spk_times); shade_window(ax(5), win);
ylabel(ax(5), 'g_{syn}'); xlabel(ax(5), 'time (ms)');
if ~isempty(h5), legend(ax(5), h5, lab5, 'Location', 'eastoutside'); end
title(ax(5), 'synaptic conductance (STP-shaped EPSGs onto pyr targets)');

linkaxes(ax, 'x'); xlim(ax(1), [t(1) t(end)]);

fprintf('\nDriven pyr neuron %d; %d total spikes across the network.\n', pyr, size(spikes, 1));

% ============================================================================
% local helpers
% ============================================================================
function plot_stp_traces(ax, t, X, active_pre, type_of, names, cmap)
    % Plot X(j,q,:) for each active presynaptic neuron j and post-type q ~= type(j)
    % (the connections that actually exist), colored by presynaptic type.
    K = size(X, 2);
    h = gobjects(0); lab = {};
    ls = {'-', '--', ':', '-.'};
    for j = active_pre
        for q = 1:K
            if q == type_of(j), continue; end          % no same-type synapse here
            tr = squeeze(X(j, q, :))';
            h(end+1) = plot(ax, t, tr, 'Color', cmap(type_of(j), :), ...
                'LineStyle', ls{mod(q-1, numel(ls)) + 1}, 'LineWidth', 1.25); %#ok<AGROW>
            lab{end+1} = sprintf('%s\\rightarrow%s', names{type_of(j)}, names{q}); %#ok<AGROW>
        end
    end
    if ~isempty(h), legend(ax, h, lab, 'Location', 'eastoutside'); end
end

function mark_spikes(ax, spk_times)
    % Faint vertical lines at spike times (drawn behind the data).
    for s = spk_times(:)'
        xline(ax, s, 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5, ...
            'HandleVisibility', 'off');
    end
end

function shade_window(ax, win)
    if isempty(win), return; end
    xregion(ax, win(1), win(2), 'FaceColor', [0.30 0.55 0.90], ...
        'FaceAlpha', 0.08, 'HandleVisibility', 'off');
end
