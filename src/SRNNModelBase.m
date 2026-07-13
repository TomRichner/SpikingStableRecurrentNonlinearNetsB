classdef SRNNModelBase < handle
    % SRNNMODELBASE Abstract base for spiking recurrent-network models.
    %
    % Holds the model-agnostic machinery shared by concrete spiking models:
    % construction / name-value parsing, the build() template, ODE integration
    % (run), Lyapunov analysis (Benettin LLE + QR spectrum), local-exponent
    % filtering, and decimation.
    %
    % Adapted from FractionalResevoir/src/SRNNModelBase.m: the lifecycle
    % (build -> run -> plot), the reflective name-value constructor, and the
    % Lyapunov static methods (benettin_algorithm_internal,
    % lyapunov_spectrum_qr_internal, compute_kaplan_yorke_dimension_internal)
    % are reused nearly verbatim so Lyapunov results are numerically comparable
    % to the rate model. RMT-connectivity dependent properties, the rate
    % activation-function library, and eigenvalue plotting were dropped.
    %
    % Model-specific behavior is supplied by concrete subclasses through the
    % abstract hooks:
    %   build_network, build_stimulus, validate, get_params,
    %   decimate_and_unpack, eval_dynamics, eval_jacobian
    %
    % eval_dynamics / eval_jacobian are the numeric "seam": run() and
    % compute_lyapunov() dispatch the RHS and Jacobian through them so each
    % subclass plugs in its own dynamics.
    %
    % IMPORTANT (time units): everything here is unit-agnostic as long as fs,
    % T_range, and lya_dt share a consistent time unit. The HH subclass works
    % in MILLISECONDS, so fs is samples/ms, T_range is in ms, and lya_dt (the
    % Benettin renormalisation interval) is in ms. Set obj.lya_dt explicitly
    % for spiking; the rate-model default (0.02) is only appropriate in seconds.
    %
    % Concrete subclass: SRNNModelHH.

    %% Network Architecture Properties
    properties
        n = 100                     % Total number of neurons
        indegree = 20               % Expected in-degree (Bernoulli target)
        level_of_chaos = 1.0        % Global scaling knob on the weight matrix W
    end

    %% Simulation Settings Properties
    properties
        fs = 40                     % Sampling frequency (samples per time unit; ms for HH)
        T_range = [0, 1000]         % Simulation interval [start, end] (ms for HH)
        ode_solver                  % Integrator handle [t,Y]=solver(rhs,tspan,y0,opts)
        ode_opts                    % Integrator options struct (ignored by fixed-step integrators)
        x0_std = 0.0                % Std of random component of the initial state
    end

    %% Input Configuration Properties
    properties
        input_config                % Struct with stimulus parameters
        u_ex_scale = 1.0            % Scaling factor for external input
        rng_seeds = [1 2]           % RNG seeds [network, stimulus]
    end

    %% Lyapunov Settings Properties
    properties
        lya_method = 'benettin'     % 'benettin', 'qr', or 'none'
        lya_dt                      % Renormalisation interval (same time unit as fs/T_range). [] -> method default
        time_units_per_second = 1   % Model time units in one second (ms model: 1000). Used to report LLE in 1/s.
        lya_T_interval              % [t_start t_end] window for accumulation ([] -> auto, skip transient)
        lya_transient = 200         % Transient (in T_range units) skipped before accumulating
        lya_d0 = 1e-3               % Benettin perturbation magnitude
        filter_local_lya = false    % Lowpass-filter the local exponent before decimation
        lya_filter_order = 2        % Butterworth order
        lya_filter_cutoff = 0.25    % Normalised cutoff (fraction of Nyquist)
    end

    %% Storage Options Properties
    properties
        store_full_state = false    % Keep the full S_out trajectory in memory
        store_decimated_state = true% Keep decimated plot_data
        plot_deci                   % Decimation factor for plotting (auto from fs/plot_freq)
        plot_freq = 1               % Target plotting frequency (in fs units)
    end

    %% Computed Properties (SetAccess = protected for subclass access)
    properties (SetAccess = protected)
        W                           % Connection weight matrix (n x n, pre x post)
        is_built = false            % Whether build() has run

        t_ex                        % Time grid for external input / integration
        u_ex                        % External input matrix (n x nt)
        u_interpolant               % griddedInterpolant for external input
        S0                          % Initial state vector
        cached_params               % Cached params struct (set by build)
    end

    %% Results Properties (conditionally stored)
    properties (SetAccess = protected)
        t_out                       % Time vector from the integrator
        S_out                       % State trajectory (nt x N_sys_eqs)
        plot_data                   % Struct with decimated data for plotting
        lya_results                 % Lyapunov analysis results struct
        has_run = false             % Whether run() has completed
    end

    %% Constructor
    methods
        function obj = SRNNModelBase(varargin)
            % SRNNMODELBASE Constructor with name-value pairs (shared by all subclasses).
            %
            % Applies defaults (set_defaults, dispatched to the most-derived
            % override), parses name-value pairs reflectively, derives plot_deci.

            obj.set_defaults();

            for i = 1:2:length(varargin)
                if isprop(obj, varargin{i})
                    obj.(varargin{i}) = varargin{i+1};
                else
                    warning('SRNNModel:UnknownProperty', 'Unknown property: %s', varargin{i});
                end
            end

            if isempty(obj.plot_deci)
                obj.plot_deci = max(1, round(obj.fs / obj.plot_freq));
            end
        end
    end

    %% Abstract hooks — implemented by concrete subclasses
    methods (Abstract, Access = protected)
        set_defaults(obj)                        % set subclass default property values
        build_network(obj)                       % construct connectivity matrix obj.W
        build_stimulus(obj)                      % generate stimulus, interpolant, and S0
        validate(obj)                            % check parameter consistency
        decimate_and_unpack(obj)                 % decimate S_out into obj.plot_data
        dS_dt = eval_dynamics(obj, t, S, params) % smooth (between-spike) RHS seam
        J = eval_jacobian(obj, S, params)        % Jacobian seam (QR only; may be unsupported)
    end

    methods (Abstract)
        params = get_params(obj)                 % pack the params struct read by the numeric kernels
    end

    %% Public lifecycle
    methods
        function build(obj)
            % BUILD Initialise the network: create W, generate stimulus, cache params.
            obj.build_network();
            obj.build_stimulus();
            obj.finalize_build();
        end

        function run(obj)
            % RUN Integrate the equations and optionally compute Lyapunov exponents.
            if ~obj.is_built
                error('SRNNModel:NotBuilt', 'Model must be built before running. Call build() first.');
            end

            params = obj.cached_params;
            dt = 1 / obj.fs;

            if isempty(obj.ode_opts)
                % Fixed-step integrators ignore these; kept for signature parity
                % and for the QR variational path (which uses ode45).
                obj.ode_opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'MaxStep', dt);
            end

            params.u_interpolant = obj.u_interpolant;
            rhs = @(t, S) obj.eval_dynamics(t, S, params);

            fprintf('Integrating equations (n=%d, %g %s, dt=%.4g)\n', ...
                obj.n, obj.T_range(2) - obj.T_range(1), 'time units', dt);
            tic
            [t_raw, S_raw] = obj.ode_solver(rhs, obj.t_ex, obj.S0, obj.ode_opts);
            fprintf('Integration complete in %.2f s.\n', toc);

            if length(t_raw) ~= length(obj.t_ex) || max(abs(t_raw(:) - obj.t_ex(:))) > 1e-6
                error('SRNNModel:TimeMismatch', ...
                    'Integrator output times do not match input grid. Max diff: %.2e', ...
                    max(abs(t_raw(:) - obj.t_ex(:))));
            end

            obj.t_out = t_raw;
            obj.S_out = S_raw;

            if ~strcmpi(obj.lya_method, 'none')
                obj.compute_lyapunov();
                if obj.filter_local_lya
                    obj.filter_lyapunov();
                end
            end

            if obj.store_decimated_state
                obj.decimate_and_unpack();
            end

            if ~obj.store_full_state
                obj.S_out = [];
            end

            obj.has_run = true;
            fprintf('Simulation complete.\n');
        end

        function compute_lyapunov(obj)
            % COMPUTE_LYAPUNOV Compute Lyapunov exponents based on lya_method.
            if isempty(obj.S_out)
                error('SRNNModel:NoStateData', ...
                    'State data not available. Set store_full_state=true or call before clearing.');
            end

            dt = 1 / obj.fs;
            params = obj.cached_params;

            if isempty(obj.lya_T_interval)
                span = obj.T_range(2) - obj.T_range(1);
                if span > obj.lya_transient
                    obj.lya_T_interval = [obj.T_range(1) + obj.lya_transient, obj.T_range(2)];
                else
                    obj.lya_T_interval = [obj.T_range(1), obj.T_range(2)];
                end
            end

            params.u_interpolant = obj.u_interpolant;
            rhs = @(t, S) obj.eval_dynamics(t, S, params);
            jac = @(tt, S, p) obj.eval_jacobian(S, p);

            fprintf('Computing Lyapunov exponents using %s method\n', obj.lya_method);
            obj.lya_results = SRNNModelBase.compute_lyapunov_exponents_internal( ...
                obj.lya_method, obj.S_out, obj.t_out, dt, obj.fs, obj.lya_T_interval, ...
                params, obj.ode_opts, obj.ode_solver, rhs, obj.t_ex, obj.u_ex, jac, ...
                obj.lya_dt, obj.lya_d0);

            % Convenience: LLE expressed in 1/s (rate-model units) regardless of
            % the model's internal time unit. For the ms-based HH model
            % time_units_per_second = 1000, so LLE_per_s = LLE * 1000.
            if isfield(obj.lya_results, 'LLE')
                obj.lya_results.LLE_per_s = obj.lya_results.LLE * obj.time_units_per_second;
                fprintf('Largest Lyapunov Exponent: %.5g /time-unit  (= %.5g /s)\n', ...
                    obj.lya_results.LLE, obj.lya_results.LLE_per_s);
            end
        end

        function filter_lyapunov(obj)
            % FILTER_LYAPUNOV Lowpass-filter the local Lyapunov exponent (before decimation).
            if isempty(obj.lya_results)
                return;
            end
            Wn = obj.lya_filter_cutoff / (obj.lya_results.lya_fs / 2);
            [b, a] = butter(obj.lya_filter_order, Wn, 'low');
            if isfield(obj.lya_results, 'local_lya') && ~isempty(obj.lya_results.local_lya)
                obj.lya_results.local_lya = filtfilt(b, a, obj.lya_results.local_lya);
            end
            if isfield(obj.lya_results, 'local_LE_spectrum_t') && ~isempty(obj.lya_results.local_LE_spectrum_t)
                for col = 1:size(obj.lya_results.local_LE_spectrum_t, 2)
                    obj.lya_results.local_LE_spectrum_t(:, col) = ...
                        filtfilt(b, a, obj.lya_results.local_LE_spectrum_t(:, col));
                end
            end
        end
    end

    %% Protected finalize
    methods (Access = protected)
        function finalize_build(obj)
            % FINALIZE_BUILD Validate parameters and cache the params struct.
            obj.validate();
            obj.cached_params = obj.get_params();
            obj.is_built = true;
        end
    end

    %% Lyapunov numerics (reused nearly verbatim from FractionalResevoir)
    methods (Static)
        function lya_results = compute_lyapunov_exponents_internal(Lya_method, S_out, t_out, dt, fs, T_interval, params, opts, ode_solver, rhs_func, t_ex, u_ex, jac_func, lya_dt_override, d0_override) %#ok<INUSD>
            % Compute Lyapunov exponents using the Benettin or QR method.
            %
            % lya_dt_override lets a caller supply the renormalisation interval in
            % the model's own time unit (required for the ms-scale HH model);
            % [] falls back to the historical per-method default.

            lya_results = struct();
            if strcmpi(Lya_method, 'none')
                return;
            end

            if nargin >= 14 && ~isempty(lya_dt_override)
                lya_dt = lya_dt_override;
            elseif strcmpi(Lya_method, 'qr')
                lya_dt = 0.1;
            elseif strcmpi(Lya_method, 'benettin')
                lya_dt = 0.02;
            else
                lya_dt = 0.1;
            end
            if nargin >= 15 && ~isempty(d0_override)
                d0 = d0_override;
            else
                d0 = 1e-3;
            end

            lya_dt_vs_dt_factor = lya_dt / dt;
            if abs(round(lya_dt_vs_dt_factor) - lya_dt_vs_dt_factor) > 1e-9
                error('lya_dt must be a multiple of dt. lya_dt/dt = %g', lya_dt_vs_dt_factor);
            end
            if lya_dt_vs_dt_factor < 3
                error('lya_dt must be >= 3*dt. lya_dt/dt = %g. Increase fs or lya_dt.', lya_dt_vs_dt_factor);
            end
            lya_fs = 1 / lya_dt;

            switch lower(Lya_method)
                case 'benettin'
                    fprintf('Computing largest Lyapunov exponent using Benettin''s algorithm...\n');
                    tic
                    [LLE, local_lya, finite_lya, t_lya] = SRNNModelBase.benettin_algorithm_internal( ...
                        S_out, t_out, dt, fs, d0, T_interval, lya_dt, params, opts, rhs_func, t_ex, u_ex, ode_solver);
                    toc
                    lya_results.LLE = LLE;
                    lya_results.local_lya = local_lya;
                    lya_results.finite_lya = finite_lya;
                    lya_results.t_lya = t_lya;
                    lya_results.lya_dt = lya_dt;
                    lya_results.lya_fs = lya_fs;

                case 'qr'
                    fprintf('Computing full Lyapunov spectrum using QR decomposition method...\n');
                    tic
                    [LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya] = ...
                        SRNNModelBase.lyapunov_spectrum_qr_internal(S_out, t_out, lya_dt, params, ...
                        ode_solver, opts, jac_func, T_interval, params.N_sys_eqs, fs);
                    toc
                    fprintf('Lyapunov Dimension: %.2f\n', ...
                        SRNNModelBase.compute_kaplan_yorke_dimension_internal(LE_spectrum));
                    [sorted_LE, sort_idx] = sort(real(LE_spectrum), 'descend');
                    lya_results.LE_spectrum = sorted_LE;
                    lya_results.local_LE_spectrum_t = local_LE_spectrum_t(:, sort_idx);
                    lya_results.finite_LE_spectrum_t = finite_LE_spectrum_t(:, sort_idx);
                    lya_results.sort_idx = sort_idx;
                    lya_results.t_lya = t_lya;
                    lya_results.lya_dt = lya_dt;
                    lya_results.lya_fs = lya_fs;
                    lya_results.LLE = lya_results.LE_spectrum(1);
                    fprintf('Largest Lyapunov Exponent (sorted): %.5g\n', lya_results.LE_spectrum(1));

                otherwise
                    error('Unknown Lyapunov method: %s', Lya_method);
            end
        end

        function [LLE, local_lya, finite_lya, t_lya] = benettin_algorithm_internal(X, t, dt, fs, d0, T, lya_dt, params, ode_options, dynamics_func, t_ex, u_ex, ode_solver) %#ok<INUSD>
            % Benettin's algorithm for the largest Lyapunov exponent.
            % Reused from FractionalResevoir/src/SRNNModelBase.m. The perturbed
            % copy is re-integrated on each interval with the SAME ode_solver
            % (for the HH model that is the event-aware hybrid integrator, so the
            % reshoot detects its own spike times), then renormalised.

            if ~isscalar(lya_dt) || ~isnumeric(lya_dt) || lya_dt <= 0
                error('lya_dt must be a positive scalar.');
            end

            deci_lya = round(lya_dt * fs);
            if deci_lya < 1
                error('lya_dt * fs must give at least 1 sample per interval.');
            end

            tau_lya = dt * deci_lya;
            t_lya = t(1:deci_lya:end);
            if t_lya(end) + tau_lya > T(2)
                t_lya(end) = [];
            end
            nt_lya = numel(t_lya);

            local_lya = zeros(nt_lya, 1);
            finite_lya = nan(nt_lya, 1);
            sum_log_stretching_factors = 0;

            n_state = size(X, 2);
            rnd_IC = randn(n_state, 1);
            pert = (rnd_IC ./ norm(rnd_IC)) .* d0;

            min_max_range = SRNNModelBase.get_minMaxRange_internal(params);
            min_bnds = min_max_range(:, 1);
            max_bnds = min_max_range(:, 2);

            % Accumulation window start (skip transient)
            t_accum_start = T(1);

            for k = 1:nt_lya
                idx_start = (k - 1) * deci_lya + 1;
                idx_end = idx_start + deci_lya;

                X_start = X(idx_start, :).';
                X_k_pert = X_start + pert;

                idx_violates_min = ~isnan(min_bnds) & (X_k_pert < min_bnds);
                X_k_pert(idx_violates_min) = min_bnds(idx_violates_min);
                idx_violates_max = ~isnan(max_bnds) & (X_k_pert > max_bnds);
                X_k_pert(idx_violates_max) = max_bnds(idx_violates_max);

                t_seg_detailed = t_lya(k) + (0:dt:tau_lya);

                [~, X_pert_output_all_steps] = ode_solver(dynamics_func, t_seg_detailed, X_k_pert, ode_options);

                X_pert_end = X_pert_output_all_steps(end, :).';
                X_end = X(idx_end, :).';

                delta = X_pert_end - X_end;
                d_k = norm(delta);
                local_lya(k) = log(d_k / d0) / tau_lya;

                if ~isfinite(local_lya(k))
                    warning('System diverged at t=%g. Truncating results.', t_lya(k));
                    last_valid = finite_lya(~isnan(finite_lya));
                    if ~isempty(last_valid), LLE = last_valid(end); else, LLE = 0; end
                    local_lya(k:end) = [];
                    finite_lya(k:end) = [];
                    t_lya(k:end) = [];
                    return;
                end

                pert = (delta ./ d_k) .* d0;

                if t_lya(k) >= t_accum_start
                    sum_log_stretching_factors = sum_log_stretching_factors + log(d_k / d0);
                    finite_lya(k, 1) = sum_log_stretching_factors / ...
                        max((t_lya(k) + tau_lya) - t_accum_start, eps);
                end
            end

            last_valid = finite_lya(~isnan(finite_lya));
            if ~isempty(last_valid), LLE = last_valid(end); else, LLE = 0; end
        end

        function [LE_spectrum, local_LE_spectrum_t, finite_LE_spectrum_t, t_lya_vec] = lyapunov_spectrum_qr_internal(X_fid_traj, t_fid_traj, lya_dt_interval, params, ode_solver, ode_options_main, jacobian_func_handle, T_full_interval, N_states_sys, fs_fid) %#ok<INUSD>
            % QR method for the full Lyapunov spectrum (variational eqs + periodic
            % reorthonormalisation). Reused from FractionalResevoir. NOTE: for the
            % event-based HH model this requires jump/saltation matrices that the
            % smooth Jacobian omits; it is retained for completeness but the HH
            % subclass supports only 'benettin' in this pass.

            fiducial_interpolants = cell(N_states_sys, 1);
            for i = 1:N_states_sys
                fiducial_interpolants{i} = griddedInterpolant(t_fid_traj, X_fid_traj(:, i), 'pchip');
            end

            dt_fid = 1 / fs_fid;
            deci_lya = round(lya_dt_interval / dt_fid);
            if deci_lya == 0
                error('lya_dt_interval is too small.');
            end
            tau_lya = dt_fid * deci_lya;

            t_lya_indices = 1:deci_lya:length(t_fid_traj);
            t_lya_vec = t_fid_traj(t_lya_indices);
            if ~isempty(t_lya_vec) && (t_lya_vec(end) + tau_lya > t_fid_traj(end) + eps(t_fid_traj(end)))
                t_lya_vec(end) = [];
                t_lya_indices(end) = [];
            end

            nt_lya = numel(t_lya_vec);
            if nt_lya == 0
                warning('No Lyapunov intervals could be formed.');
                LE_spectrum = nan(N_states_sys, 1);
                local_LE_spectrum_t = [];
                finite_LE_spectrum_t = [];
                return;
            end

            Q_current = eye(N_states_sys);
            sum_log_R_diag = zeros(N_states_sys, 1);
            local_LE_spectrum_t = zeros(nt_lya, N_states_sys);
            finite_LE_spectrum_t = nan(nt_lya, N_states_sys);
            total_positive_time_accumulated = 0;
            ode_options_var = ode_options_main;

            for k = 1:nt_lya
                t_start_segment = t_lya_vec(k);
                t_end_segment = min(t_start_segment + tau_lya, t_fid_traj(end));
                current_segment_duration = t_end_segment - t_start_segment;
                if current_segment_duration <= eps
                    if k > 1
                        local_LE_spectrum_t(k, :) = local_LE_spectrum_t(k-1, :);
                        finite_LE_spectrum_t(k, :) = finite_LE_spectrum_t(k-1, :);
                    else
                        local_LE_spectrum_t(k, :) = NaN;
                        finite_LE_spectrum_t(k, :) = NaN;
                    end
                    continue;
                end

                t_span_ode = [t_start_segment, t_end_segment];
                Psi0_vec = reshape(Q_current, [], 1);
                variational_eqs = @(tt, current_Psi_vec) SRNNModelBase.variational_eqs_ode_internal( ...
                    tt, current_Psi_vec, fiducial_interpolants, N_states_sys, jacobian_func_handle, params);
                [~, Psi_solution_vec] = ode_solver(variational_eqs, t_span_ode, Psi0_vec, ode_options_var);
                Psi_evolved_matrix = reshape(Psi_solution_vec(end, :)', [N_states_sys, N_states_sys]);

                if any(~isfinite(Psi_evolved_matrix(:)))
                    warning('System diverged at t=%g. Truncating results.', t_start_segment);
                    if total_positive_time_accumulated > eps
                        LE_spectrum = sum_log_R_diag / total_positive_time_accumulated;
                    else
                        LE_spectrum = nan(N_states_sys, 1);
                    end
                    t_lya_vec(k:end) = [];
                    local_LE_spectrum_t(k:end, :) = [];
                    finite_LE_spectrum_t(k:end, :) = [];
                    return;
                end

                [Q_new, R_segment] = qr(Psi_evolved_matrix);
                diag_R = diag(R_segment);
                log_abs_diag_R = log(abs(diag_R));
                valid_diag_R = abs(diag_R) > eps;

                current_local_LEs = zeros(N_states_sys, 1);
                current_local_LEs(valid_diag_R) = log_abs_diag_R(valid_diag_R) / current_segment_duration;
                current_local_LEs(~valid_diag_R) = -Inf;
                local_LE_spectrum_t(k, :) = current_local_LEs';

                if t_start_segment >= T_full_interval(1) - eps(0)
                    sum_log_R_diag(valid_diag_R) = sum_log_R_diag(valid_diag_R) + log_abs_diag_R(valid_diag_R);
                    total_positive_time_accumulated = total_positive_time_accumulated + current_segment_duration;
                end

                if total_positive_time_accumulated > eps
                    finite_LE_spectrum_t(k, :) = (sum_log_R_diag / total_positive_time_accumulated)';
                elseif k > 1
                    finite_LE_spectrum_t(k, :) = finite_LE_spectrum_t(k-1, :);
                else
                    finite_LE_spectrum_t(k, :) = NaN;
                end

                Q_current = Q_new;
            end

            if total_positive_time_accumulated > eps
                LE_spectrum = sum_log_R_diag / total_positive_time_accumulated;
            else
                warning('No accumulation over positive time for global LEs.');
                LE_spectrum = nan(N_states_sys, 1);
            end
        end

        function dPsi_vec_dt = variational_eqs_ode_internal(tt, current_Psi_vec, fiducial_interpolants, N_states_sys, jacobian_func_handle, params)
            % Variational ODE for the QR method.
            X_fid_at_tt = zeros(N_states_sys, 1);
            for s = 1:N_states_sys
                X_fid_at_tt(s) = fiducial_interpolants{s}(tt);
            end
            J_matrix = jacobian_func_handle(tt, X_fid_at_tt, params);
            Psi_matrix = reshape(current_Psi_vec, [N_states_sys, N_states_sys]);
            dPsi_vec_dt = reshape(J_matrix * Psi_matrix, [], 1);
        end

        function min_max_range = get_minMaxRange_internal(params)
            % Per-state bounds for Benettin perturbation clipping. Default: none
            % (NaN). Length matches the state dimension via params.N_sys_eqs.
            min_max_range = nan(params.N_sys_eqs, 2);
        end

        function D_KY = compute_kaplan_yorke_dimension_internal(lambda)
            % Kaplan-Yorke (Lyapunov) dimension from a spectrum.
            lambda = sort(lambda, 'descend');
            cumsum_lambda = cumsum(lambda);
            j = find(cumsum_lambda >= 0, 1, 'last');
            if isempty(j)
                D_KY = 0;
            elseif j == length(lambda)
                D_KY = length(lambda);
            else
                D_KY = j + cumsum_lambda(j) / abs(lambda(j + 1));
            end
        end
    end
end
