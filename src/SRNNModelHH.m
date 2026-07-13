classdef SRNNModelHH < SRNNModelBase
    % SRNNMODELHH  Spiking Hodgkin-Huxley recurrent network with SFA, STD, STF.
    %
    % Biophysical (spiking) analog of the rate model SRNNModelCellTypes: a
    % randomly connected network of cortical Traub-Miles HH neurons of K cell
    % types (default pyr / pvalb / sst / vip, from Campagnola 2022), carrying
    %   - spike-frequency adaptation (SFA): a spike-triggered K-adaptation
    %     current per adapting neuron (linear negative feedback), and
    %   - short-term synaptic depression (STD) + facilitation (STF): event-based
    %     Tsodyks-Markram release resources per (presynaptic neuron, post-type).
    %
    % The network is a HYBRID dynamical system: HH neurons + synaptic/SFA
    % variables evolve as a smooth ODE (eval_dynamics), and presynaptic spikes
    % (upward V-threshold crossings) trigger discrete jumps. Integration and the
    % Benettin reshoot both use the event-aware integrate_hh_hybrid, so the
    % largest Lyapunov exponent of the spiking network can be compared to the
    % rate model's.
    %
    % UNITS: ms, mV, uF/cm^2, mS/cm^2, uA/cm^2. Campagnola time constants (s)
    % are converted to ms on load.
    %
    % STATE LAYOUT (cursor-walked, guarded blocks):
    %   S = [V(n); m(n); h(n); ng(n); a(n_ad*n_a); b(n*K); p(n*K); g(n*K)]
    % a present iff n_a>0 (adapting neurons only); b iff n_b>0; p iff n_u>0;
    % g (synaptic conductance per post-neuron x pre-type) always present.
    %
    % WEIGHT ORIENTATION: obj.W is (pre x post): W(j,i) = |weight| from pre j to
    % post i. Sign is realized by synaptic reversal potential (exc pre -> E_exc,
    % inh pre -> E_inh), not by a signed weight.
    %
    % Lifecycle: model = SRNNModelHH(...); model.build(); model.run(); model.plot();
    %
    % See also: SRNNModelBase, integrate_hh_hybrid, hh_gating_rates,
    %           load_campagnola_matrices.

    %% Cell-type configuration
    properties
        type_names = {'pyr', 'pvalb', 'sst', 'vip'}   % K cell types; type 1 = excitatory
        type_fractions = [0.80, 0.08, 0.07, 0.05]     % renormalized at build
        exc_type_names = {'pyr'}                       % which type_names are excitatory
        use_campagnola_data = true
    end

    %% Mechanism configuration
    properties
        n_a = 1                 % # SFA timescales (0 disables SFA)
        n_b = 1                 % STD present? (0 or 1)
        n_u = 1                 % STF present? (0 or 1)
        c_gain = 0.7            % maps adapt_index -> SFA conductance strength c (mS/cm^2)
        a_incr = 1.0            % SFA state increment per own spike
        tau_a = 100             % fallback SFA time constant (ms)
        sfa_min_index = 0.01    % types below this are non-adapting (carry no SFA state)
        w_cv = 1.0              % per-edge weight heterogeneity (std/|mean|)
        tau_b_rel_ref = 500     % reference STD release tau (ms), provisional heuristic
    end

    %% HH biophysics (fixed ionic gradients)
    properties
        C_m  = 1.0              % membrane capacitance (uF/cm^2)
        gNa  = 30               % max Na conductance (mS/cm^2)
        gK   = 25               % max K conductance (mS/cm^2)
        gL   = 0.1              % leak conductance (mS/cm^2)
        E_Na = 50               % Na reversal (mV)
        E_K  = -100             % K reversal (mV)
        E_L  = -67              % leak reversal (mV)
    end

    %% Synapse / spike-detection
    properties
        g_syn_scale  = 0.15     % overall synaptic conductance scale (mS/cm^2)
        E_exc = 0               % excitatory synaptic reversal (mV)
        E_inh = -75             % inhibitory synaptic reversal (mV)
        tau_syn_exc = 2         % excitatory synaptic decay (ms)
        tau_syn_inh = 6         % inhibitory synaptic decay (ms)
        V_th = -20              % spike-detect threshold (mV)
        V_reset = -40           % re-arm (hysteresis) level (mV)
    end

    %% Per-type parameter tables (K x K unless noted; filled at build if empty)
    properties
        conn_prob      % K x K connection probability (pre -> post)
        psp_amp        % K x K signed PSP amplitude (magnitude used for |W|)
        dep_tau        % K x K STD recovery tau (ms)
        dep_amount     % K x K STD depression amount
        rel_prob       % K x K baseline release probability p0
        fac_tau        % K x K STF facilitation tau (ms)
        kappa          % K x K STF facilitation increment coefficient
        adapt_index    % K x 1 SFA adaptation index
        tau_a_type     % K x 1 per-type SFA tau (ms)
    end

    %% Computed (SetAccess = protected)
    properties (SetAccess = protected)
        type_of        % n x 1 integer type label
        is_exc         % n x 1 logical (excitatory neuron)
        exc_type       % K x 1 logical (excitatory type)
        n_types        % K
        adapting       % n x 1 logical (carries SFA state)
        ad_idx         % indices of adapting neurons
        n_ad           % count of adapting neurons
    end

    %% Constructor
    methods
        function obj = SRNNModelHH(varargin)
            obj@SRNNModelBase(varargin{:});
        end

        function [t, Y] = run_integrator(obj, odefun, tspan, y0, opts)
            % RUN_INTEGRATOR Event-aware integration seam (set as obj.ode_solver).
            % Bundles the jump-parameter struct (cached at build) into the
            % standalone integrate_hh_hybrid so both run() and the Benettin
            % reshoot apply identical spike jumps.
            [t, Y] = integrate_hh_hybrid(odefun, tspan, y0, opts, obj.cached_params.jump);
        end
    end

    %% Subclass hooks
    methods (Access = protected)
        function set_defaults(obj)
            % SET_DEFAULTS Spiking / ms-scale defaults.
            obj.fs        = 40;             % dt = 0.025 ms
            obj.T_range   = [0, 1000];      % ms
            obj.n         = 100;
            obj.indegree  = 20;
            obj.x0_std    = 1.0;            % mV jitter on initial V
            obj.ode_solver = @obj.run_integrator;

            % Lyapunov (ms units): 20 ms renormalisation interval, skip 200 ms transient.
            obj.lya_method    = 'benettin';
            obj.lya_dt        = 20;
            obj.lya_transient = 200;
            obj.lya_d0        = 1e-3;
            obj.time_units_per_second = 1000;   % model time is ms -> report LLE in 1/s

            % Stimulus: sparse current steps + optional tonic bias (uA/cm^2).
            obj.input_config = struct();
            obj.input_config.n_steps        = 3;
            obj.input_config.amp            = 8.0;    % step current (uA/cm^2)
            obj.input_config.step_density   = 0.3;    % fraction of neurons driven per step
            obj.input_config.no_stim_pattern = false(1, 3);
            obj.input_config.bias           = 0.0;    % constant background current
            obj.input_config.positive_only  = true;
            % Optional deterministic targeted drive (overrides random steps when
            % drive_types is non-empty): a rectangular current pulse of drive_amp
            % applied to all neurons of the given types within drive_window [t_on t_off] ms.
            obj.input_config.drive_types    = [];     % cell-type indices to drive ([] -> random steps)
            obj.input_config.drive_window   = [];     % [t_on t_off] ms ([] -> whole run)
            obj.input_config.drive_amp      = 0.0;    % pulse amplitude (uA/cm^2)
        end

        function build_network(obj)
            % BUILD_NETWORK Assign types, gather per-type tables, build |W| (pre x post).
            rng(obj.rng_seeds(1));
            obj.n_types = numel(obj.type_names);
            obj.assign_types();
            obj.load_parameter_tables();

            % Adapting neurons: adaptation_index >= threshold AND SFA enabled.
            type_adapts = (obj.adapt_index(:) >= obj.sfa_min_index);
            obj.adapting = (obj.n_a > 0) & type_adapts(obj.type_of);
            obj.ad_idx = find(obj.adapting);
            obj.n_ad = numel(obj.ad_idx);

            n = obj.n; K = obj.n_types; t = obj.type_of;

            % Block-structured NON-NEGATIVE conductance magnitudes, W(pre, post).
            W = zeros(n, n);
            for tp = 1:K                                  % presynaptic type -> ROWS
                pre = find(t == tp);
                for ts = 1:K                              % postsynaptic type -> COLS
                    post = find(t == ts);
                    if isempty(pre) || isempty(post), continue; end
                    p  = min(max(obj.conn_prob(tp, ts), 0), 1);
                    mu = abs(obj.psp_amp(tp, ts));
                    sig = obj.w_cv * mu;
                    npre = numel(pre); npost = numel(post);
                    blk = max(mu + sig * randn(npre, npost), 0) .* (rand(npre, npost) < p);
                    W(pre, post) = blk;
                end
            end
            W(1:n+1:end) = 0;                             % no autapses

            % Normalise typical nonzero magnitude to 1, then apply conductance scale.
            nz = W(W > 0);
            if ~isempty(nz)
                W = W / median(nz);
            end
            obj.W = obj.level_of_chaos * obj.g_syn_scale * W;

            fprintf('HH network: n=%d, K=%d, mean in-degree=%.1f, max |g_syn| bump=%.3f mS/cm^2\n', ...
                n, K, mean(sum(obj.W > 0, 1)), max(obj.W(:)));
        end

        function build_stimulus(obj)
            % BUILD_STIMULUS Sparse current-step input + tonic bias, interpolant, S0.
            dt = 1 / obj.fs;
            T  = obj.T_range(2);
            t_stim = (obj.T_range(1):dt:T)';
            nt = numel(t_stim);
            n = obj.n;

            rng(obj.rng_seeds(2));
            ic = obj.input_config;
            u_stim = ic.bias * ones(n, nt);

            if isfield(ic, 'drive_types') && ~isempty(ic.drive_types)
                % Deterministic targeted drive: a rectangular pulse of drive_amp on
                % all neurons of the requested types within drive_window.
                driven = ismember(obj.type_of, ic.drive_types);
                if isempty(ic.drive_window)
                    win = true(1, nt);
                else
                    win = (t_stim' >= ic.drive_window(1)) & (t_stim' < ic.drive_window(2));
                end
                u_stim(driven, win) = u_stim(driven, win) + ic.drive_amp;
            else
                % Default: sparse random current steps.
                step_len = max(1, round((T - obj.T_range(1)) / ic.n_steps * obj.fs));
                if ic.positive_only
                    steps = ic.amp * abs(randn(n, ic.n_steps));
                else
                    steps = ic.amp * randn(n, ic.n_steps);
                end
                steps = steps .* (rand(n, ic.n_steps) < ic.step_density);
                steps(:, ic.no_stim_pattern) = 0;
                for s = 1:ic.n_steps
                    a0 = (s - 1) * step_len + 1;
                    b0 = min(s * step_len, nt);
                    if a0 > nt, break; end
                    u_stim(:, a0:b0) = u_stim(:, a0:b0) + repmat(steps(:, s), 1, b0 - a0 + 1);
                end
            end

            obj.t_ex = t_stim;
            obj.u_ex = u_stim * obj.u_ex_scale;
            obj.u_interpolant = griddedInterpolant(obj.t_ex, obj.u_ex', 'linear', 'nearest');
            obj.S0 = obj.initialize_state(obj.get_params());
            fprintf('HH stimulus: %d time points (%g-%g ms), %d neurons\n', nt, obj.T_range(1), T, n);
        end

        function validate(obj)
            if obj.n < obj.n_types
                error('SRNNModelHH:TooFewNeurons', 'n (%d) must be >= K (%d).', obj.n, obj.n_types);
            end
            if ~ismember(obj.n_b, [0 1]) || ~ismember(obj.n_u, [0 1])
                error('SRNNModelHH:BadTimescales', 'n_b and n_u must be 0 or 1.');
            end
            if obj.n_a < 0
                error('SRNNModelHH:BadTimescales', 'n_a must be >= 0.');
            end
            if obj.T_range(2) <= obj.T_range(1)
                error('SRNNModelHH:BadT', 'T_range(2) must be > T_range(1).');
            end
            if obj.V_reset >= obj.V_th
                error('SRNNModelHH:BadThreshold', 'V_reset (%g) must be < V_th (%g).', obj.V_reset, obj.V_th);
            end
        end

        function dS_dt = eval_dynamics(~, t, S, params)
            dS_dt = SRNNModelHH.dynamics_hh(t, S, params);
        end

        function J = eval_jacobian(~, ~, ~)
            % Not supported this pass: the event-based (hybrid) system needs
            % jump/saltation matrices for a correct variational flow. Benettin
            % (full nonlinear reshoot) does not use the Jacobian.
            error('SRNNModelHH:NoJacobian', ...
                ['Analytic Jacobian / QR spectrum is not supported for the ', ...
                 'event-based HH model. Use lya_method=''benettin''.']);
        end

        function decimate_and_unpack(obj)
            params = obj.cached_params;
            deci = obj.plot_deci;
            idx = 1:deci:size(obj.S_out, 1);
            t_plot = obj.t_out(idx);
            st = SRNNModelHH.unpack_states_hh(obj.S_out(idx, :), params);

            pd = struct();
            pd.t = t_plot;
            pd.u_ext = obj.u_ex(:, idx);
            pd.V = st.V;                 % n x nt
            pd.spikes = SRNNModelHH.detect_spikes(obj.t_out, obj.S_out(:, 1:obj.n), obj.V_th, obj.V_reset);
            pd.b = st.b;                 % n x K x nt
            pd.p = st.p;                 % n x K x nt
            pd.g = st.g;                 % n x K x nt
            if ~isempty(st.a), pd.a = st.a; else, pd.a = []; end
            pd.type_of = obj.type_of;
            pd.type_names = obj.type_names;
            obj.plot_data = pd;
        end

        function S0 = initialize_state(obj, params)
            % INITIALIZE_STATE V near E_L (steady-state gating), a=0, b=1, p=p0, g=0.
            n = params.n; K = params.K;
            V0 = obj.E_L * ones(n, 1) + obj.x0_std * randn(n, 1);
            [minf, hinf, ninf] = hh_gating_inf(V0);
            a0 = zeros(params.n_ad * params.n_a, 1);
            if params.n_b > 0, b0 = ones(n * K, 1); else, b0 = []; end
            if params.n_u > 0, p0 = params.p0_mat(:); else, p0 = []; end
            g0 = zeros(n * K, 1);
            S0 = [V0; minf; hinf; ninf; a0; b0; p0; g0];
        end

        function assign_types(obj)
            % ASSIGN_TYPES Contiguous per-type blocks sized by type_fractions.
            K = obj.n_types; n = obj.n;
            fr = obj.type_fractions(:)'; fr = fr / sum(fr);
            counts = floor(fr * n);
            counts(1) = counts(1) + (n - sum(counts));    % remainder to type 1
            tof = zeros(n, 1); c = 0;
            for k = 1:K
                tof(c + (1:counts(k))) = k;
                c = c + counts(k);
            end
            obj.type_of = tof;
            obj.exc_type = ismember(lower(obj.type_names(:)), lower(obj.exc_type_names));
            obj.is_exc = obj.exc_type(tof);
        end

        function load_parameter_tables(obj)
            % LOAD_PARAMETER_TABLES Fill empty K x K / K x 1 tables from Campagnola
            % data (with s->ms conversion) or from defaults. NaNs -> defaults; clamps.
            K = obj.n_types;
            def = SRNNModelHH.default_tables(K);

            if obj.use_campagnola_data
                C = load_campagnola_matrices();
                [tf, ord] = ismember(lower(obj.type_names), lower(C.types));
                if ~all(tf)
                    error('SRNNModelHH:UnknownType', ...
                        'type_names must be a subset of Campagnola types {pyr,pvalb,sst,vip} when use_campagnola_data=true.');
                end
                data.conn_prob   = SRNNModelHH.pick(C.conn_prob_adj(ord, ord), C.conn_prob(ord, ord));
                data.psp_amp     = C.psp_amplitude(ord, ord);
                data.dep_tau     = C.ml_depression_tau(ord, ord)    * 1000;   % s -> ms
                data.dep_amount  = C.ml_depression_amount(ord, ord);
                data.rel_prob    = C.ml_release_prob(ord, ord);
                data.fac_tau     = C.ml_facilitation_tau(ord, ord)  * 1000;   % s -> ms
                data.kappa       = C.ml_facilitation_amount(ord, ord);
                data.adapt_index = C.sfa_adaptation_index(ord);
                data.tau_a_type  = C.sfa_tau(ord) * 1000;                     % s -> ms
            else
                data = def;
            end

            % Fill each table if the user has not set it, replacing NaNs with defaults.
            obj.conn_prob   = SRNNModelHH.fill(obj.conn_prob,   data.conn_prob,   def.conn_prob);
            obj.psp_amp     = SRNNModelHH.fill(obj.psp_amp,     data.psp_amp,     def.psp_amp);
            obj.dep_tau     = SRNNModelHH.fill(obj.dep_tau,     data.dep_tau,     def.dep_tau);
            obj.dep_amount  = SRNNModelHH.fill(obj.dep_amount,  data.dep_amount,  def.dep_amount);
            obj.rel_prob    = SRNNModelHH.fill(obj.rel_prob,    data.rel_prob,    def.rel_prob);
            obj.fac_tau     = SRNNModelHH.fill(obj.fac_tau,     data.fac_tau,     def.fac_tau);
            obj.kappa       = SRNNModelHH.fill(obj.kappa,       data.kappa,       def.kappa);
            obj.adapt_index = SRNNModelHH.fill(obj.adapt_index, data.adapt_index, def.adapt_index);
            obj.tau_a_type  = SRNNModelHH.fill(obj.tau_a_type,  data.tau_a_type,  def.tau_a_type);

            % Clamps to keep dynamics well-posed.
            obj.conn_prob  = min(max(obj.conn_prob, 0), 1);
            obj.dep_tau    = min(max(obj.dep_tau, 10), 5000);
            obj.rel_prob   = min(max(obj.rel_prob, 0.05), 0.95);
            obj.fac_tau    = min(max(obj.fac_tau, 10), 5000);
            obj.kappa      = max(obj.kappa, 0);
            obj.tau_a_type = min(max(obj.tau_a_type, 10), 5000);
        end
    end

    %% Public
    methods
        function params = get_params(obj)
            % GET_PARAMS Pack the smooth-RHS params and the integrator jump struct.
            n = obj.n; K = obj.n_types;
            params = struct();
            params.n = n; params.K = K;
            params.n_a = obj.n_a; params.n_b = obj.n_b; params.n_u = obj.n_u;
            params.n_ad = obj.n_ad; params.ad_idx = obj.ad_idx; params.adapting = obj.adapting;
            params.type_of = obj.type_of;
            params.N_sys_eqs = 4 * n + (obj.n_a > 0) * obj.n_ad * obj.n_a + ...
                (obj.n_b > 0) * n * K + (obj.n_u > 0) * n * K + n * K;

            % HH constants
            params.C_m = obj.C_m; params.gNa = obj.gNa; params.gK = obj.gK; params.gL = obj.gL;
            params.E_Na = obj.E_Na; params.E_K = obj.E_K; params.E_L = obj.E_L;

            % SFA per-neuron
            if obj.n_a > 0
                params.tau_a = repmat(obj.tau_a_type(obj.type_of), 1, obj.n_a);   % n x n_a
            else
                params.tau_a = repmat(obj.tau_a, n, max(obj.n_a, 1));
            end
            params.c_vec = obj.c_gain * obj.adapt_index(obj.type_of);             % n x 1

            % STD/STF per (pre-neuron j, post-type q): gather K x K -> n x K by pre-type.
            t = obj.type_of;
            params.tau_b_rec_mat = obj.dep_tau(t, :);                             % n x K
            params.p0_mat        = obj.rel_prob(t, :);                            % n x K
            params.tau_f_mat     = obj.fac_tau(t, :);                             % n x K

            % Synaptic reversal / decay per PRE type -> per column.
            E_syn_type   = obj.E_inh * ones(K, 1);  E_syn_type(obj.exc_type)  = obj.E_exc;
            tau_syn_type = obj.tau_syn_inh * ones(K, 1); tau_syn_type(obj.exc_type) = obj.tau_syn_exc;
            params.E_syn_vec  = E_syn_type;          % K x 1 (indexed by pre-type)
            params.tau_syn_row = tau_syn_type(:)';   % 1 x K

            params.W = obj.W;

            % Jump-parameter struct consumed by integrate_hh_hybrid.
            jp = struct();
            jp.N = n; jp.K = K; jp.n_a = obj.n_a; jp.n_ad = obj.n_ad;
            jp.ad_idx = obj.ad_idx; jp.type_of = obj.type_of;
            jp.Wabs = obj.W;                          % already |.| and scaled (pre x post)
            jp.kappa = obj.kappa;                     % K x K (pre x post)
            jp.p0_mat = params.p0_mat;                % n x K
            jp.a_incr = obj.a_incr;
            jp.V_th = obj.V_th; jp.V_reset = obj.V_reset;
            jp.has_a = obj.n_a > 0; jp.has_b = obj.n_b > 0; jp.has_p = obj.n_u > 0;
            params.jump = jp;
        end

        function [fig, ax] = plot(obj)
            % PLOT Spike raster + example V trace + population mechanism means.
            if ~obj.has_run || isempty(obj.plot_data)
                error('SRNNModelHH:NotRun', 'Run the model before plotting.');
            end
            pd = obj.plot_data;
            has_lya = ~strcmpi(obj.lya_method, 'none') && ~isempty(obj.lya_results);
            n_panels = 3 + double(has_lya);
            fig = figure('Color', 'w', 'Name', 'SRNNModelHH');
            tl = tiledlayout(fig, n_panels, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
            ax = gobjects(n_panels, 1);

            % 1: raster
            ax(1) = nexttile(tl); hold(ax(1), 'on');
            if ~isempty(pd.spikes)
                scatter(ax(1), pd.spikes(:, 1), pd.spikes(:, 2), 4, 'k', 'filled');
            end
            ylabel(ax(1), 'neuron'); ylim(ax(1), [0, obj.n + 1]); box(ax(1), 'off');
            title(ax(1), sprintf('HH network raster (n=%d)', obj.n));

            % 2: example membrane potentials (first neuron of each type)
            ax(2) = nexttile(tl); hold(ax(2), 'on');
            for k = 1:obj.n_types
                i = find(obj.type_of == k, 1);
                if ~isempty(i), plot(ax(2), pd.t, pd.V(i, :)); end
            end
            ylabel(ax(2), 'V (mV)'); box(ax(2), 'off');
            legend(ax(2), obj.type_names, 'Location', 'eastoutside');

            % 3: population mechanism means
            ax(3) = nexttile(tl); hold(ax(3), 'on');
            leg = {};
            if obj.n_a > 0 && ~isempty(pd.a)
                plot(ax(3), pd.t, squeeze(mean(sum(pd.a, 2), 1))); leg{end+1} = 'SFA a';
            end
            if obj.n_b > 0
                plot(ax(3), pd.t, squeeze(mean(mean(pd.b, 2), 1))); leg{end+1} = 'STD b';
            end
            if obj.n_u > 0
                plot(ax(3), pd.t, squeeze(mean(mean(pd.p, 2), 1))); leg{end+1} = 'STF p';
            end
            ylabel(ax(3), 'pop. mean'); box(ax(3), 'off');
            if ~isempty(leg), legend(ax(3), leg, 'Location', 'eastoutside'); end

            % 4: local Lyapunov exponent
            if has_lya
                ax(4) = nexttile(tl); hold(ax(4), 'on');
                lr = obj.lya_results;
                plot(ax(4), lr.t_lya, lr.finite_lya * obj.time_units_per_second, 'LineWidth', 1.25);
                yline(ax(4), 0, 'k:');
                ylabel(ax(4), '\lambda_1 (1/s)'); box(ax(4), 'off');
                title(ax(4), sprintf('finite-time LLE -> %.4g /s', lr.LLE_per_s));
            end
            xlabel(ax(end), 'time (ms)');
        end
    end

    %% Static numeric kernels / helpers
    methods (Static)
        function dS_dt = dynamics_hh(t, S, params)
            % DYNAMICS_HH Smooth (between-spike) RHS of the hybrid HH network.
            n = params.n; K = params.K; n_a = params.n_a;

            idx = 0;
            V  = S(1:n);            idx = n;
            m  = S(idx + (1:n));    idx = idx + n;
            h  = S(idx + (1:n));    idx = idx + n;
            ng = S(idx + (1:n));    idx = idx + n;

            len_a = params.n_ad * n_a * (n_a > 0);
            if len_a > 0, a_ad = reshape(S(idx + (1:len_a)), params.n_ad, n_a); else, a_ad = []; end
            idx = idx + len_a;

            len_b = n * K * (params.n_b > 0);
            if len_b > 0, b = reshape(S(idx + (1:len_b)), n, K); else, b = ones(n, K); end
            idx = idx + len_b;

            len_p = n * K * (params.n_u > 0);
            if len_p > 0, p = reshape(S(idx + (1:len_p)), n, K); else, p = params.p0_mat; end
            idx = idx + len_p;

            len_g = n * K;
            g = reshape(S(idx + (1:len_g)), n, K);

            % --- HH channel currents ---
            [am, bm, ah, bh, an, bn] = hh_gating_rates(V);
            INa = params.gNa .* m.^3 .* h .* (V - params.E_Na);
            IK  = params.gK  .* ng.^4       .* (V - params.E_K);
            IL  = params.gL  .* (V - params.E_L);

            % --- SFA adaptation current (scatter a to full n, sum timescales) ---
            if len_a > 0
                a_full = zeros(n, n_a); a_full(params.ad_idx, :) = a_ad;
                I_SFA = params.c_vec .* sum(a_full, 2) .* (V - params.E_K);
            else
                I_SFA = zeros(n, 1);
            end

            % --- synaptic current: sum over presynaptic type P ---
            I_syn = zeros(n, 1);
            for P = 1:K
                I_syn = I_syn + g(:, P) .* (V - params.E_syn_vec(P));
            end

            % --- external drive ---
            u = params.u_interpolant(t).';

            dV = (u - INa - IK - IL - I_SFA - I_syn) ./ params.C_m;
            dm = am .* (1 - m) - bm .* m;
            dh = ah .* (1 - h) - bh .* h;
            dn = an .* (1 - ng) - bn .* ng;

            if len_a > 0, da = -a_ad ./ params.tau_a(params.ad_idx, :); else, da = []; end
            if len_b > 0, db = (1 - b) ./ params.tau_b_rec_mat;         else, db = []; end
            if len_p > 0, dp = (params.p0_mat - p) ./ params.tau_f_mat; else, dp = []; end
            dg = -g ./ params.tau_syn_row;

            dS_dt = [dV; dm; dh; dn; da(:); db(:); dp(:); dg(:)];
        end

        function st = unpack_states_hh(S, params)
            % UNPACK_STATES_HH Reshape a (nt x N_sys) trajectory into named fields.
            % Returns V,m,h,ng (n x nt) and a (n x n_a x nt), b/p/g (n x K x nt).
            n = params.n; K = params.K; n_a = params.n_a;
            nt = size(S, 1);
            idx = 0;
            V  = S(:, 1:n).';           idx = n;
            m  = S(:, idx + (1:n)).';   idx = idx + n;
            h  = S(:, idx + (1:n)).';   idx = idx + n;
            ng = S(:, idx + (1:n)).';   idx = idx + n;

            len_a = params.n_ad * n_a * (n_a > 0);
            if len_a > 0
                a_ad = reshape(S(:, idx + (1:len_a)).', params.n_ad, n_a, nt);
                a = zeros(n, n_a, nt); a(params.ad_idx, :, :) = a_ad;
            else
                a = [];
            end
            idx = idx + len_a;

            len_b = n * K * (params.n_b > 0);
            if len_b > 0, b = reshape(S(:, idx + (1:len_b)).', n, K, nt); else, b = ones(n, K, nt); end
            idx = idx + len_b;

            len_p = n * K * (params.n_u > 0);
            if len_p > 0, p = reshape(S(:, idx + (1:len_p)).', n, K, nt); else, p = repmat(params.p0_mat, 1, 1, nt); end
            idx = idx + len_p;

            len_g = n * K;
            g = reshape(S(:, idx + (1:len_g)).', n, K, nt);

            st = struct('V', V, 'm', m, 'h', h, 'ng', ng, 'a', a, 'b', b, 'p', p, 'g', g);
        end

        function spikes = detect_spikes(t, V, V_th, V_reset)
            % DETECT_SPIKES Post-hoc raster from a (nt x n) V trajectory.
            % Returns [time, neuron] rows at upward V_th crossings (hysteresis).
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

        function def = default_tables(K)
            % DEFAULT_TABLES Reasonable K x K / K x 1 defaults (ms units) used when
            % Campagnola data is disabled or to fill NaNs. Excitatory = type 1.
            conn = 0.12 * ones(K);
            psp  = 2e-4 * ones(K);               % magnitude; sign irrelevant (|.| used)
            def = struct();
            def.conn_prob   = conn;
            def.psp_amp     = psp;
            def.dep_tau     = 300 * ones(K);     % ms
            def.dep_amount  = 0.3 * ones(K);
            def.rel_prob    = 0.25 * ones(K);
            def.fac_tau     = 200 * ones(K);     % ms
            def.kappa       = 0.1 * ones(K);
            ai = 0.05 * ones(K, 1); ai(1) = 0.07;
            def.adapt_index = ai;
            def.tau_a_type  = 100 * ones(K, 1);  % ms
        end

        function A = pick(A, B)
            % PICK Return A but substitute B where A is NaN (data preference chain).
            m = isnan(A); A(m) = B(m);
        end

        function out = fill(user_val, data_val, def_val)
            % FILL Choose user_val if non-empty, else data_val; replace NaNs with def_val.
            if ~isempty(user_val)
                out = user_val;
            else
                out = data_val;
            end
            m = isnan(out);
            out(m) = def_val(m);
        end
    end
end
