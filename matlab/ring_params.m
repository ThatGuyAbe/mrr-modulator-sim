function p = ring_params(P, Q_target, alpha_scale, tau_c_scale)
% normalized coefficients of one ring (time in units of 1/gamma_0).
%   the coupling time comes from the target loaded Q with the nominal loss, then the
%   loss and the coupling time are scaled by alpha_scale and tau_c_scale (1 and 1 for
%   the nominal ring).
    if nargin < 3, alpha_scale = 1; end
    if nargin < 4, tau_c_scale = 1; end
    alpha_Npm_nom = P.alpha_dBcm_nom * 100 * log(10)/10;
    gamma_lin_nom = P.c_light * alpha_Npm_nom / P.n_Si;
    gamma_0_target    = P.omega_0 / (2 * Q_target);
    gamma_coup_target = gamma_0_target - gamma_lin_nom;
    if gamma_coup_target <= 0, error('Q = %d needs a negative coupling loss', Q_target); end
    tau_c_nom = 2 / gamma_coup_target;

    alpha_dBcm = P.alpha_dBcm_nom * alpha_scale;
    tau_c      = tau_c_nom * tau_c_scale;
    alpha_Npm  = alpha_dBcm * 100 * log(10)/10;
    gamma_lin  = P.c_light * alpha_Npm / P.n_Si;
    gamma_coup = 2 / tau_c;
    gamma_0    = gamma_lin + gamma_coup;
    mu_raw     = sqrt(2 / tau_c);

    p.gamma_0     = gamma_0;
    p.mu_raw      = mu_raw;
    p.mu_n        = mu_raw / gamma_0;
    p.glin_n      = gamma_lin / gamma_0;
    p.cTPA_n      = (P.beta_TPA * P.c_light^2 / (P.n_Si^2 * P.V_TPA)) / gamma_0;
    p.cFCA_n      = (P.Gamma_FCA * P.sigma_FCA * P.c_light / (2*P.n_Si)) / gamma_0;
    p.cTh_n       = P.Gamma_th / (P.m_ring * P.c_p * gamma_0);
    p.iTFC_n      = (1/P.tau_FC) / gamma_0;
    p.iTTH_n      = (1/P.tau_th) / gamma_0;
    p.dw_NL_scale = P.omega_0 / (P.n_Si * gamma_0);
    p.dw_mod_n    = -P.dW0_dVpn / gamma_0;
    p.iTRC_n      = (1/P.tau_RC) / gamma_0;
    % for reference, not used by the solver
    p.tau_c   = tau_c;
    p.Q_ring  = P.omega_0 / (2*gamma_0);
    p.hw_GHz  = gamma_0 / (2*pi*1e9);
    p.gamma_lin = gamma_lin;
end
