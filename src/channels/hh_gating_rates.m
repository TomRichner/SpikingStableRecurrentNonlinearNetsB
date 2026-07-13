function [alpha_m, beta_m, alpha_h, beta_h, alpha_n, beta_n] = hh_gating_rates(v)
%HH_GATING_RATES  Traub-Miles / Wei HH gating rate constants (vectorized, mV, 1/ms).
%
%   [am, bm, ah, bh, an, bn] = HH_GATING_RATES(v) returns the alpha/beta
%   opening/closing rates for the m, h and n gating variables of a cortical
%   Hodgkin-Huxley neuron, evaluated elementwise for a column (or any shape)
%   of membrane potentials v (mV). Rates are in 1/ms.
%
%   The rate expressions are the Traub-Miles cortical constants reused from
%   CorticalSpreadDepolarizationModel/model_int.m. Three of the six rates
%   (alpha_m, beta_m, alpha_n) have a removable 0/0 singularity where both the
%   numerator (v - v*) and the (1 - exp(...)) / (exp(...) - 1) denominator
%   vanish. Inside a tiny window around each singular voltage the raw value is
%   OVERWRITTEN (logical-index assignment) with the analytic L'Hopital limit:
%       alpha_m: v* = -54, limit = 0.32*4 = 1.28
%       beta_m : v* = -27, limit = 0.28*5 = 1.40
%       alpha_n: v* = -52, limit = 0.032*5 = 0.16
%   Direct assignment (not the multiply-add idiom in model_int.m) is used so
%   that a voltage landing EXACTLY on the singularity, where the raw expression
%   is NaN, is still patched correctly -- NaN.*0 would be NaN, not 0.
%   alpha_h and beta_h are smooth exponentials with no singularity.
%
%   See also: hh_gating_inf, SRNNModelHH.

    tol = 1.5e-4;   % half-width (mV) of the singularity-patch window

    % --- m gate -----------------------------------------------------------
    alpha_m = 0.32 .* (54 + v) ./ (1 - exp(-(v + 54) ./ 4));
    alpha_m(abs(-54 - v) < tol) = 1.28;

    beta_m = 0.28 .* (v + 27) ./ (exp((v + 27) ./ 5) - 1);
    beta_m(abs(-27 - v) < tol) = 1.4;

    % --- h gate (no singularity) -----------------------------------------
    alpha_h = 0.128 .* exp(-(50 + v) ./ 18);
    beta_h  = 4 ./ (1 + exp(-(v + 27) ./ 5));

    % --- n gate -----------------------------------------------------------
    alpha_n = 0.032 .* (v + 52) ./ (1 - exp(-(v + 52) ./ 5));
    alpha_n(abs(-52 - v) < tol) = 0.16;

    beta_n = 0.5 .* exp(-(v + 57) ./ 40);
end
