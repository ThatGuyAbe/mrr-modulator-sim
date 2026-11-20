function R = ring_population(S, P)
% the 2000-ring draw (rng(42)): resonance offset, then loss scale, then coupling-time
% scale, scales clipped at 0.5. z_res are the normal draws for the tuning error
% (step 7, seed 2026). in matlab the draw is checked against ring_population_seed42.mat
% and the script stops if they differ. octave has a different generator, so it uses
% the saved draw.
    N = S.N_rings;
    F = load(fullfile(fileparts(mfilename('fullpath')), 'ring_population_seed42.mat'));
    is_octave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if ~is_octave
        rng(S.seed);
        dlam_fab    = S.sigma_lam * randn(2000, 1);
        df_fab_Hz   = (P.c_light / P.lambda_0^2) * dlam_fab;
        df_fab_GHz  = df_fab_Hz / 1e9;
        alpha_scale = 1 + S.sigma_Q_rel * randn(2000, 1);
        alpha_scale(alpha_scale < 0.5) = 0.5;
        tau_c_scale = 1 + S.sigma_Q_rel * randn(2000, 1);
        tau_c_scale(tau_c_scale < 0.5) = 0.5;
        rs = RandStream('mt19937ar', 'Seed', 2026);
        z_res = randn(rs, 2000, 1);
        d = max([max(abs(df_fab_GHz - F.df_fab_GHz)), max(abs(alpha_scale - F.alpha_scale)), ...
                 max(abs(tau_c_scale - F.tau_c_scale)), max(abs(z_res - F.z_res))]);
        if d > 1e-12
            error('ring population differs from the saved draw (max difference %.2e)', d);
        end
        R.check = d;
    else
        df_fab_GHz = F.df_fab_GHz;  alpha_scale = F.alpha_scale;  tau_c_scale = F.tau_c_scale;
        z_res = F.z_res;  R.check = NaN;
    end
    R.df_fab_GHz  = df_fab_GHz(1:N);
    R.alpha_scale = alpha_scale(1:N);
    R.tau_c_scale = tau_c_scale(1:N);
    R.z_res       = z_res(1:N);
    R.N           = N;
end
