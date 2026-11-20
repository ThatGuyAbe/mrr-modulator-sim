function [P, S] = mrr_params()
% model parameters (P) and run settings (S), SI units.
%
%   every script takes its values from here, so they live in one place.
%   ring mass from the geometry: 2*pi*7.5 um * 500 nm * 220 nm * 2329 kg/m^3 = 1.2e-14 kg.
%   set the environment variable MRR_SMOKE to 1 for a quick test run (coarse grid,
%   12 rings, output in out_smoke/). never use those numbers.

    % physical constants
    P.c_light = 2.99792458e8;          % m/s
    P.hbar    = 1.054571817e-34;       % J s

    % cavity and material (giron castro 2024 unless noted)
    P.m_ring     = 1.2e-14;            % kg, from the ring geometry
    P.n_Si       = 3.485;
    P.lambda_0   = 1553.49e-9;         % m
    P.c_p        = 700;                % J/(kg K)
    P.V_FCA      = 2.36e-18;           % m^3
    P.V_TPA      = 2.59e-18;           % m^3
    P.beta_TPA   = 8.4e-12;            % m/W, de cea 2019
    P.Gamma_FCA  = 0.9996;
    P.Gamma_th   = 0.9355;
    P.dn_dT      = 1.86e-4;            % 1/K
    P.dn_dN      = -1.73e-27;          % m^3
    P.sigma_FCA  = 1.0e-21;            % m^2
    P.tau_th     = 50e-9;              % s
    P.tau_FC     = 0.1e-9;             % s, de cea 2019
    P.alpha_dBcm_nom = 0.8;            % dB/cm

    % ring geometry, used for the mass
    P.ring_radius = 7.5e-6;            % m
    P.wg_width    = 500e-9;            % m
    P.wg_height   = 220e-9;            % m
    P.rho_Si      = 2329;              % kg/m^3

    % modulator (de cea 2019)
    P.dW0_dVpn = 1.57e10;              % rad/(s V)
    P.tau_RC   = 10e-12;               % s
    P.V_bias   = -2.5;                 % V
    P.V_pp     = 4.0;                  % V

    % derived
    P.omega_0 = 2*pi * P.c_light / P.lambda_0;
    P.V_high  = P.V_bias + P.V_pp/2;   % zero bit, -0.5 V
    P.V_low   = P.V_bias - P.V_pp/2;   % one bit, -4.5 V
    P.label   = 'corrected (m 1.2e-14 kg, beta 8.4e-12 m/W)';

    % run settings
    S.N_sym         = 1000;            % data symbols (prbs-10)
    S.prbs_order    = 10;
    S.dt_real       = 1e-12;           % s, rk4 step
    S.sigma_rx_frac = 0.05;            % receiver noise floor, fraction of P_in
    S.n_tau_ramp    = 5;               % warm-up ramp from the blue side, in tau_th (250 ns)
    S.n_tau_hold    = 5;               % then a hold at the target, in tau_th (250 ns)
    S.df_grid       = -80:0.25:80;     % GHz, detuning grid of every sweep
    S.Q_values      = [3000 6000 9000 12000 15000];
    S.powers        = [2 3 4];         % mW
    S.datarates     = [2.5 15];        % Gb/s
    S.P_lin_mW      = 0.01;            % low-power reference for dw_lin
    S.rate_lin      = 2.5;             % Gb/s
    S.N_rings       = 2000;
    S.n_asf         = 500;             % as-fabricated case on rings 1..n_asf only
    S.sigma_lam     = 0.65e-9;         % m, resonance offset spread (about 81 GHz)
    S.sigma_Q_rel   = 0.10;            % loss and coupling spread
    S.seed          = 42;
    S.out_root      = 'out';
    S.cw_lo         = -50;             % GHz, red end of the cw sweeps in step 4
    S.smoke         = strcmp(getenv('MRR_SMOKE'), '1');
    if S.smoke
        S.df_grid   = -20:0.5:12;
        S.N_rings   = 12;
        S.n_asf     = 4;
        S.out_root  = 'out_smoke';
    end
    S.P_low = min(S.powers);           % mW
    S.P_ref = max(S.powers);           % mW
end
