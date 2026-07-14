classdef SRNNModelHHEI < SRNNModelBase
    % SRNNMODELHHEI  Spiking HH E/I network comparable to the rate model SRNNModel2.
    %
    % Biophysical (spiking) analog of the hand-tuned E/I RMT rate model
    % SRNNModel2 (FractionalResevoir): two populations (excitatory / inhibitory)
    % wired with Harris-2023 Random-Matrix-Theory connectivity, carrying
    %   - multi-timescale spike-frequency adaptation (SFA) as a K-current, with
    %     PER-POPULATION timescale counts, and
    %   - multi-timescale short-term depression (STD), per presynaptic neuron,
    %     with per-population timescale counts. NO facilitation (STF).
    % so edge-of-chaos (Benettin LLE) can be compared between rate and spiking
    % descriptions.
    %
    % Differences from SRNNModelHH (the Campagnola cell-type model):
    %   - Connectivity is RMT (RMTMatrix), not Bernoulli+Campagnola. Wrong-sign
    %     entries are ZEROED; sign is realized by synaptic reversal potential
    %     (hard Dale); indegree/alpha is the sparsity knob.
    %   - SFA/STD timescale COUNTS are per population (n_a_vec, n_b_vec), giving a
    %     ragged per-type state layout that mirrors SRNNModel2's [a_E;a_I;b_E;b_I].
    %   - SFA pools are DC-gain balanced (increment Delta_l = a_incr0/tau_l) so each
    %     logspaced timescale contributes equally at steady state (matching the rate
    %     model's a_l* = r), rather than the slow pool dominating.
    %   - STD efficacy is the PRODUCT over a presynaptic neuron's pools (Prod_m b_m),
    %     multiplying the conductance bump; each pool depletes b_m -= p0*b_m per spike
    %     and recovers at its own tau_rec,m. STF is absent.
    %
    % UNITS: ms, mV, uF/cm^2, mS/cm^2, uA/cm^2. All time constants are set by hand
    % in ms (= the rate model's seconds x 1000).
    %
    % STATE LAYOUT (per-type ragged blocks; types contiguous, exc first):
    %   S = [V(N); m(N); h(N); ng(N);
    %        a-blocks (per type t: count_t x n_a(t));
    %        b-blocks (per type t: count_t x n_b(t));
    %        g(N x K)]
    %
    % WEIGHT ORIENTATION: obj.W is (pre x post) magnitude; sign via reversal.
    %
    % Lifecycle: model = SRNNModelHHEI(...); model.build(); model.run(); model.plot();
    %
    % See also: SRNNModelBase, SRNNModelHH, RMTMatrix, integrate_hh_events.

    %% Populations
    properties
        type_names = {'exc', 'inh'}      % type 1 = excitatory (first contiguous block)
        type_fractions = [0.5, 0.5]      % E/I fractions (f = type_fractions(1))
        exc_type_names = {'exc'}
    end

    %% RMT connectivity (tilde notation, Harris 2023); [] -> filled at build
    properties
        mu_E_tilde       % default  3F
        mu_I_tilde       % default -4F
        sigma_E_tilde    % default  F
        sigma_I_tilde    % default  F
        zrs_mode = 'none'
    end

    %% SFA (K-current, per-population multi-timescale)
    properties
        n_a_vec = [3, 0]                 % SFA timescale count per type
        tau_a_range = [250, 10000]       % logspace range (ms) when tau_a_cell empty
        tau_a_cell = {}                  % optional explicit per-type tau vectors (ms)
        c_type = [0.15, 0]               % per-type SFA total strength (c_eff = c_type/n_a)
    end

    %% STD (per presynaptic neuron, per-population multi-timescale). STF OFF.
    properties
        n_b_vec = [1, 0]                 % STD timescale count per type
        tau_rec_range = [200, 2000]      % logspace range (ms) for multi-timescale STD
        tau_rec_default = 1000           % single-timescale recovery (ms)
        tau_rec_cell = {}                % optional explicit per-type tau_rec vectors (ms)
        p0_type = [0.2, 0.2]             % per-type release/depletion fraction per spike
    end

    %% HH biophysics
    properties
        C_m  = 1.0
        gNa  = 30
        gK   = 25
        gL   = 0.1
        E_Na = 50
        E_K  = -100
        E_L  = -67
    end

    %% Synapse / spike detection
    properties
        g_syn_scale  = 0.15
        E_exc = 0
        E_inh = -75
        tau_syn_exc = 2
        tau_syn_inh = 6
        V_th = -20
        V_reset = -40
    end

    %% Calibration
    properties
        target_rate = 5                  % Hz; sets the DC-gain SFA increment a_incr0 = 1000/target_rate
    end

    %% Computed (SetAccess = protected)
    properties (SetAccess = protected)
        type_of
        is_exc
        exc_type
        n_types
    end

    %% Constructor
    methods
        function obj = SRNNModelHHEI(varargin)
            obj@SRNNModelBase(varargin{:});
        end

        function [t, Y] = run_integrator(obj, odefun, tspan, y0, ~)
            % RUN_INTEGRATOR Event-aware integration seam (set as obj.ode_solver).
            p = obj.cached_params;
            jump_fn = @(y, spiked) SRNNModelHHEI.apply_hhei_jumps(y, spiked, p);
            [t, Y] = integrate_hh_events(odefun, tspan, y0, p.N, jump_fn, obj.V_th, obj.V_reset);
        end
    end

    %% Subclass hooks
    methods (Access = protected)
        function set_defaults(obj)
            obj.fs        = 40;            % dt = 0.025 ms
            obj.T_range   = [0, 2000];     % ms
            obj.n         = 300;
            obj.indegree  = 100;           % -> alpha = 1/3 (RMT sparsity)
            obj.x0_std    = 1.0;
            obj.ode_solver = @obj.run_integrator;

            obj.lya_method    = 'benettin';
            obj.lya_dt        = 20;        % ms
            obj.lya_transient = 500;       % ms (raise for slow SFA pools)
            obj.lya_d0        = 1e-3;
            obj.time_units_per_second = 1000;

            obj.input_config = struct();
            obj.input_config.bias           = 0.0;   % tonic background current (uA/cm^2)
            obj.input_config.drive_types    = [];
            obj.input_config.drive_window   = [];
            obj.input_config.drive_amp      = 0.0;
        end

        function build_network(obj)
            % BUILD_NETWORK RMT connectivity, wrong-sign zeroed, magnitude + Dale sign.
            rng(obj.rng_seeds(1));
            obj.n_types = numel(obj.type_names);
            obj.assign_types();

            n = obj.n;
            alpha = obj.indegree / n;
            F = 1 / sqrt(n * alpha * (2 - alpha));
            if isempty(obj.mu_E_tilde),    obj.mu_E_tilde    =  3 * F; end
            if isempty(obj.mu_I_tilde),    obj.mu_I_tilde    = -4 * F; end
            if isempty(obj.sigma_E_tilde), obj.sigma_E_tilde =  1 * F; end
            if isempty(obj.sigma_I_tilde), obj.sigma_I_tilde =  1 * F; end

            rmt = RMTMatrix(n);
            rmt.f = nnz(obj.is_exc) / n;                 % exact: E = first block
            rmt.mu_tilde_e = obj.mu_E_tilde;
            rmt.mu_tilde_i = obj.mu_I_tilde;
            rmt.sigma_tilde_e = obj.sigma_E_tilde;
            rmt.sigma_tilde_i = obj.sigma_I_tilde;
            rmt.zrs_mode = obj.zrs_mode;
            rmt.alpha = alpha;                           % setter re-draws sparsity mask

            W_rmt = rmt.W;                               % (post x pre)
            Wpre  = W_rmt.';                             % (pre  x post)

            % Zero wrong-sign entries by presynaptic row (hard Dale).
            exc_rows = obj.is_exc;
            Wpre(exc_rows, :)  = max(Wpre(exc_rows, :), 0);    % exc pre: keep positive
            Wpre(~exc_rows, :) = min(Wpre(~exc_rows, :), 0);   % inh pre: keep negative

            Wabs = abs(Wpre);
            Wabs(1:n+1:end) = 0;                         % no autapses
            obj.W = obj.level_of_chaos * obj.g_syn_scale * Wabs;

            fprintf(['HHEI RMT network: n=%d, alpha=%.3f, F=%.4f, mu_E=%.3f mu_I=%.3f, ', ...
                'mean in-degree=%.1f\n'], n, alpha, F, obj.mu_E_tilde, obj.mu_I_tilde, ...
                mean(sum(obj.W > 0, 1)));
        end

        function build_stimulus(obj)
            dt = 1 / obj.fs;
            T  = obj.T_range(2);
            t_stim = (obj.T_range(1):dt:T)';
            nt = numel(t_stim);
            n = obj.n;

            rng(obj.rng_seeds(2));
            ic = obj.input_config;
            u_stim = ic.bias * ones(n, nt);
            if isfield(ic, 'drive_types') && ~isempty(ic.drive_types)
                driven = ismember(obj.type_of, ic.drive_types);
                if isempty(ic.drive_window)
                    win = true(1, nt);
                else
                    win = (t_stim' >= ic.drive_window(1)) & (t_stim' < ic.drive_window(2));
                end
                u_stim(driven, win) = u_stim(driven, win) + ic.drive_amp;
            end

            obj.t_ex = t_stim;
            obj.u_ex = u_stim * obj.u_ex_scale;
            obj.u_interpolant = griddedInterpolant(obj.t_ex, obj.u_ex', 'linear', 'nearest');
            obj.S0 = obj.initialize_state(obj.get_params());
            fprintf('HHEI stimulus: %d time points (%g-%g ms), %d neurons\n', nt, obj.T_range(1), T, n);
        end

        function validate(obj)
            K = numel(obj.type_names);
            if obj.n < K
                error('SRNNModelHHEI:TooFewNeurons', 'n (%d) must be >= K (%d).', obj.n, K);
            end
            if numel(obj.n_a_vec) ~= K || numel(obj.n_b_vec) ~= K
                error('SRNNModelHHEI:BadCounts', 'n_a_vec and n_b_vec must have length K=%d.', K);
            end
            if any(obj.n_a_vec < 0) || any(obj.n_b_vec < 0) || ...
               any(mod(obj.n_a_vec, 1) ~= 0) || any(mod(obj.n_b_vec, 1) ~= 0)
                error('SRNNModelHHEI:BadCounts', 'n_a_vec and n_b_vec must be non-negative integers.');
            end
            if obj.T_range(2) <= obj.T_range(1)
                error('SRNNModelHHEI:BadT', 'T_range(2) must be > T_range(1).');
            end
            if obj.V_reset >= obj.V_th
                error('SRNNModelHHEI:BadThreshold', 'V_reset (%g) must be < V_th (%g).', obj.V_reset, obj.V_th);
            end
        end

        function dS_dt = eval_dynamics(~, t, S, params)
            dS_dt = SRNNModelHHEI.dynamics_hhei(t, S, params);
        end

        function J = eval_jacobian(~, ~, ~)
            error('SRNNModelHHEI:NoJacobian', ...
                'Analytic Jacobian / QR is unsupported for the event-based HH model. Use benettin.');
        end

        function decimate_and_unpack(obj)
            params = obj.cached_params;
            deci = obj.plot_deci;
            idx = 1:deci:size(obj.S_out, 1);
            st = SRNNModelHHEI.unpack_states_hhei(obj.S_out(idx, :), params);
            pd = struct();
            pd.t = obj.t_out(idx);
            pd.V = st.V;
            pd.a_sum = st.a_sum;
            pd.b_prod = st.b_prod;
            pd.g = st.g;
            pd.spikes = SRNNModelHHEI.detect_spikes(obj.t_out, obj.S_out(:, 1:obj.n), obj.V_th, obj.V_reset);
            pd.type_of = obj.type_of;
            pd.type_names = obj.type_names;
            obj.plot_data = pd;
        end

        function S0 = initialize_state(obj, params)
            n = params.N; L = params.layout;
            V0 = obj.E_L * ones(n, 1) + obj.x0_std * randn(n, 1);
            [minf, hinf, ninf] = hh_gating_inf(V0);
            S0 = [V0; minf; hinf; ninf; zeros(L.len_a, 1); ones(L.len_b, 1); zeros(L.len_g, 1)];
        end

        function assign_types(obj)
            K = obj.n_types; n = obj.n;
            fr = obj.type_fractions(:)'; fr = fr / sum(fr);
            counts = floor(fr * n);
            counts(1) = counts(1) + (n - sum(counts));
            tof = zeros(n, 1); c = 0;
            for k = 1:K
                tof(c + (1:counts(k))) = k;
                c = c + counts(k);
            end
            obj.type_of = tof;
            obj.exc_type = ismember(lower(obj.type_names(:)), lower(obj.exc_type_names));
            obj.is_exc = obj.exc_type(tof);
        end

        % ---- per-type parameter resolvers ----
        function tau = resolve_tau_a(obj, t)
            na = obj.n_a_vec(t);
            if numel(obj.tau_a_cell) >= t && ~isempty(obj.tau_a_cell{t})
                tau = obj.tau_a_cell{t}(:)';
            elseif na == 1
                tau = sqrt(prod(obj.tau_a_range));
            else
                tau = logspace(log10(obj.tau_a_range(1)), log10(obj.tau_a_range(2)), na);
            end
        end

        function tau = resolve_tau_rec(obj, t)
            nb = obj.n_b_vec(t);
            if numel(obj.tau_rec_cell) >= t && ~isempty(obj.tau_rec_cell{t})
                tau = obj.tau_rec_cell{t}(:)';
            elseif nb == 1
                tau = obj.tau_rec_default;
            else
                tau = logspace(log10(obj.tau_rec_range(1)), log10(obj.tau_rec_range(2)), nb);
            end
        end
    end

    %% Public
    methods
        function params = get_params(obj)
            N = obj.n; K = obj.n_types;
            params = struct();
            params.N = N; params.K = K; params.n = N;
            params.type_of = obj.type_of;

            params.C_m = obj.C_m; params.gNa = obj.gNa; params.gK = obj.gK; params.gL = obj.gL;
            params.E_Na = obj.E_Na; params.E_K = obj.E_K; params.E_L = obj.E_L;

            E_syn_type = obj.E_inh * ones(K, 1);  E_syn_type(obj.exc_type)  = obj.E_exc;
            tau_syn_type = obj.tau_syn_inh * ones(K, 1); tau_syn_type(obj.exc_type) = obj.tau_syn_exc;
            params.E_syn_vec  = E_syn_type;
            params.tau_syn_row = tau_syn_type(:)';

            % DC-gain SFA increment: <a_l> = a_incr0 * nu (per-ms rate) per pool,
            % independent of tau_l when increment = a_incr0/tau_l.
            a_incr0 = 1000 / obj.target_rate;

            % ---- build the ragged per-type layout ----
            L = struct(); L.N = N; L.K = K; L.a_incr0 = a_incr0;
            off = 4 * N;
            L.off_a = off;
            sfa = struct('type', {}, 'idx', {}, 'count', {}, 'n_a', {}, 'tau', {}, 'c', {}, 'off', {}, 'len', {});
            for t = 1:K
                if obj.n_a_vec(t) > 0
                    idx = find(obj.type_of == t);
                    s = struct();
                    s.type = t; s.idx = idx; s.count = numel(idx); s.n_a = obj.n_a_vec(t);
                    s.tau = obj.resolve_tau_a(t);
                    s.c = obj.c_type(t) / max(obj.n_a_vec(t), 1);   % c_eff = c_total / n_a
                    s.off = off; s.len = s.count * s.n_a; off = off + s.len;
                    sfa(end+1) = s; %#ok<AGROW>
                end
            end
            L.len_a = off - L.off_a;
            L.off_b = off;
            std = struct('type', {}, 'idx', {}, 'count', {}, 'n_b', {}, 'tau_rec', {}, 'p0', {}, 'off', {}, 'len', {});
            for t = 1:K
                if obj.n_b_vec(t) > 0
                    idx = find(obj.type_of == t);
                    s = struct();
                    s.type = t; s.idx = idx; s.count = numel(idx); s.n_b = obj.n_b_vec(t);
                    s.tau_rec = obj.resolve_tau_rec(t);
                    s.p0 = obj.p0_type(t);
                    s.off = off; s.len = s.count * s.n_b; off = off + s.len;
                    std(end+1) = s; %#ok<AGROW>
                end
            end
            L.len_b = off - L.off_b;
            L.off_g = off; L.len_g = N * K;
            L.sfa = sfa; L.std = std;
            params.layout = L;
            params.N_sys_eqs = 4 * N + L.len_a + L.len_b + L.len_g;

            % Jump needs
            params.Wabs = obj.W;
            params.V_th = obj.V_th; params.V_reset = obj.V_reset;
        end

        function [fig, ax] = plot(obj)
            if ~obj.has_run || isempty(obj.plot_data)
                error('SRNNModelHHEI:NotRun', 'Run the model before plotting.');
            end
            pd = obj.plot_data;
            has_lya = ~strcmpi(obj.lya_method, 'none') && ~isempty(obj.lya_results);
            np = 4 + double(has_lya);
            fig = figure('Color', 'w', 'Name', 'SRNNModelHHEI'); ax = gobjects(np, 1);
            tl = tiledlayout(fig, np, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
            % Excitatory = dark red, inhibitory = blue (consistent across panels).
            c_exc = [0.55 0.08 0.08];
            c_inh = [0.15 0.40 0.75];
            cmap = repmat(c_inh, obj.n_types, 1);
            cmap(obj.exc_type, :) = repmat(c_exc, nnz(obj.exc_type), 1);
            is_e = obj.exc_type(obj.type_of);

            ax(1) = nexttile(tl); hold(ax(1), 'on');
            if ~isempty(pd.spikes)
                sp_type = pd.type_of(pd.spikes(:, 2));      % type of each spiking neuron
                scatter(ax(1), pd.spikes(:, 1), pd.spikes(:, 2), 4, cmap(sp_type, :), 'filled');
            end
            ylabel(ax(1), 'neuron'); ylim(ax(1), [0 obj.n + 1]);
            title(ax(1), 'raster (excitatory = dark red, inhibitory = blue)');

            ax(2) = nexttile(tl); hold(ax(2), 'on');
            ie = find(is_e, 1); ii = find(~is_e, 1);
            if ~isempty(ie), plot(ax(2), pd.t, pd.V(ie, :), 'Color', cmap(1, :)); end
            if ~isempty(ii), plot(ax(2), pd.t, pd.V(ii, :), 'Color', cmap(2, :)); end
            ylabel(ax(2), 'V (mV)'); legend(ax(2), obj.type_names, 'Location', 'eastoutside');

            ax(3) = nexttile(tl); hold(ax(3), 'on');
            plot(ax(3), pd.t, mean(pd.a_sum(is_e, :), 1), 'Color', cmap(1, :));
            if any(~is_e), plot(ax(3), pd.t, mean(pd.a_sum(~is_e, :), 1), 'Color', cmap(2, :)); end
            ylabel(ax(3), 'SFA \Sigma a'); title(ax(3), 'spike-frequency adaptation (pop. mean)');

            ax(4) = nexttile(tl); hold(ax(4), 'on');
            plot(ax(4), pd.t, mean(pd.b_prod(is_e, :), 1), 'Color', cmap(1, :));
            if any(~is_e), plot(ax(4), pd.t, mean(pd.b_prod(~is_e, :), 1), 'Color', cmap(2, :)); end
            ylabel(ax(4), 'STD \Pi b'); ylim(ax(4), [0 1.05]); title(ax(4), 'short-term depression (pop. mean)');

            if has_lya
                ax(5) = nexttile(tl); hold(ax(5), 'on');
                lr = obj.lya_results;
                plot(ax(5), lr.t_lya, lr.finite_lya * obj.time_units_per_second, 'LineWidth', 1.25);
                yline(ax(5), 0, 'k:'); ylabel(ax(5), '\lambda_1 (1/s)');
                title(ax(5), sprintf('finite-time LLE -> %.4g /s', lr.LLE_per_s));
            end
            xlabel(ax(end), 'time (ms)'); linkaxes(ax, 'x');
        end
    end

    %% Static numeric kernels
    methods (Static)
        function dS_dt = dynamics_hhei(t, S, params)
            L = params.layout; N = L.N; K = L.K;
            V  = S(1:N);
            m  = S(N + (1:N));
            h  = S(2*N + (1:N));
            ng = S(3*N + (1:N));

            [am, bm, ah, bh, an, bn] = hh_gating_rates(V);
            INa = params.gNa .* m.^3 .* h .* (V - params.E_Na);
            IK  = params.gK  .* ng.^4       .* (V - params.E_K);
            IL  = params.gL  .* (V - params.E_L);

            % SFA K-current + decay (per type block)
            I_SFA = zeros(N, 1);
            da_all = zeros(L.len_a, 1); pos = 0;
            for s = L.sfa
                a_block = reshape(S(s.off + (1:s.len)), s.count, s.n_a);
                I_SFA(s.idx) = I_SFA(s.idx) + s.c .* sum(a_block, 2) .* (V(s.idx) - params.E_K);
                da = -a_block ./ s.tau;               % row-broadcast per pool
                da_all(pos + (1:s.len)) = da(:); pos = pos + s.len;
            end

            % STD recovery (per type block)
            db_all = zeros(L.len_b, 1); pos = 0;
            for s = L.std
                b_block = reshape(S(s.off + (1:s.len)), s.count, s.n_b);
                db = (1 - b_block) ./ s.tau_rec;
                db_all(pos + (1:s.len)) = db(:); pos = pos + s.len;
            end

            % synaptic conductance current + decay
            g = reshape(S(L.off_g + (1:L.len_g)), N, K);
            I_syn = zeros(N, 1);
            for P = 1:K
                I_syn = I_syn + g(:, P) .* (V - params.E_syn_vec(P));
            end
            dg = -g ./ params.tau_syn_row;

            u = params.u_interpolant(t).';
            dV = (u - INa - IK - IL - I_SFA - I_syn) ./ params.C_m;
            dm = am .* (1 - m) - bm .* m;
            dh = ah .* (1 - h) - bh .* h;
            dn = an .* (1 - ng) - bn .* ng;

            dS_dt = [dV; dm; dh; dn; da_all; db_all; dg(:)];
        end

        function y = apply_hhei_jumps(y, spiked, params)
            % Discrete jumps for the E/I per-type ragged layout: DC-gain SFA
            % increment, product-efficacy conductance bump, per-pool STD depletion.
            L = params.layout; N = L.N; K = L.K;
            Wabs = params.Wabs;
            g = reshape(y(L.off_g + (1:L.len_g)), N, K);

            % SFA increments (Delta_l = a_incr0 / tau_l), per type block
            for s = L.sfa
                loc = spiked(s.idx);
                if any(loc)
                    a_block = reshape(y(s.off + (1:s.len)), s.count, s.n_a);
                    a_block(loc, :) = a_block(loc, :) + L.a_incr0 ./ s.tau;
                    y(s.off + (1:s.len)) = a_block(:);
                end
            end

            % Per-type: efficacy (product of STD pools) + conductance bump + depletion.
            % Efficacy defaults to 1 for types with no STD.
            eff = ones(N, 1);
            for s = L.std
                b_block = reshape(y(s.off + (1:s.len)), s.count, s.n_b);
                eff(s.idx) = prod(b_block, 2);
                loc = spiked(s.idx);
                if any(loc)
                    b_block(loc, :) = b_block(loc, :) .* (1 - s.p0);   % deplete
                    y(s.off + (1:s.len)) = b_block(:);
                end
            end

            for t = 1:K
                pre = find(spiked & (params.type_of == t));
                if isempty(pre), continue; end
                g(:, t) = g(:, t) + Wabs(pre, :).' * eff(pre);   % N x 1
            end
            y(L.off_g + (1:L.len_g)) = g(:);
        end

        function st = unpack_states_hhei(S, params)
            L = params.layout; N = L.N; K = L.K; nt = size(S, 1);
            st.V = S(:, 1:N).';
            st.a_sum = zeros(N, nt);
            for s = L.sfa
                blk = reshape(S(:, s.off + (1:s.len)).', s.count, s.n_a, nt);
                st.a_sum(s.idx, :) = reshape(sum(blk, 2), s.count, nt);
            end
            st.b_prod = ones(N, nt);
            for s = L.std
                blk = reshape(S(:, s.off + (1:s.len)).', s.count, s.n_b, nt);
                st.b_prod(s.idx, :) = reshape(prod(blk, 2), s.count, nt);
            end
            st.g = reshape(S(:, L.off_g + (1:L.len_g)).', N, K, nt);
        end

        function spikes = detect_spikes(t, V, V_th, V_reset)
            [nt, n] = size(V);
            armed = V(1, :) < V_th;
            rows = [];
            for k = 2:nt
                sp = armed & (V(k-1, :) < V_th) & (V(k, :) >= V_th);
                if any(sp)
                    j = find(sp);
                    rows = [rows; [repmat(t(k), numel(j), 1), j(:)]]; %#ok<AGROW>
                end
                armed(sp) = false;
                armed(V(k, :) < V_reset) = true;
            end
            spikes = rows;
        end
    end
end
