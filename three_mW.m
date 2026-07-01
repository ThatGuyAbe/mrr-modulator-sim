% three_mW.m
%
% Full sweep: 3 powers x 4 data rates = 12 configs.
% For each (P_in, data rate):
%   1. sweep detuning for the nominal OMA/ER/IL curves
%   2. run the 2000-ring experiment (as-fabricated vs corrected at delta_w=0)
%   3. plot OMA vs ring index and IL vs ring index
%
% Both the resonance offset and the Q factor vary ring to ring, each drawn
% from a Gaussian.
%
% Needs sim_mrr_modulator.m, mod_rk4_step.m, gen_prbs7.m

clear; clc; close all;

%% output folder for figures
if ~exist('paper', 'dir'), mkdir('paper'); end

%% Physical constants
c_light = 2.99792458e8;
hbar    = 1.054571817e-34;

%% MRR parameters (Giron Castro 2024, Table 1)
m_ring     = 1.2e-11;
tau_c_nom  = 54.7e-12;
n_Si       = 3.485;
lambda_0   = 1553.49e-9;
c_p        = 700;
V_FCA      = 2.36e-18;
V_TPA      = 2.59e-18;
beta_TPA   = 8.4e-11;
Gamma_FCA  = 0.9996;
Gamma_th   = 0.9355;
dn_dT      = 1.86e-4;
dn_dN      = -1.73e-27;
sigma_FCA  = 1.0e-21;
tau_th     = 50e-9;
tau_FC     = 0.1e-9;
alpha_dBcm_nom = 0.8;
omega_0    = 2*pi * c_light / lambda_0;

%% Modulator parameters (de Cea 2019)
dW0_dVpn  = 1.57e10;
tau_RC    = 10e-12;
V_bias    = -2.5;
V_pp      = 4.0;

%% Fabrication variation parameters
sigma_lam   = 0.65e-9;    % resonance offset std dev (Parsons 2025)
sigma_Q_rel = 0.10;       % Q variation: 10% relative std dev
N_rings     = 2000;

%% Build all 12 configurations: 3 powers x 4 data rates
powers    = [2.00, 3.00, 4.00];                    % mW
datarates = [2.5, 5.0, 10.0, 15.0];                % Gbps

configs = zeros(length(powers)*length(datarates), 2);
idx = 0;
for ip = 1:length(powers)
    for id = 1:length(datarates)
        idx = idx + 1;
        configs(idx,:) = [powers(ip), datarates(id)];
    end
end
N_configs = size(configs, 1);
fprintf('Total configurations: %d\n', N_configs);

%% Generate fabrication offsets (same for all configs)
rng(42);

% resonance offset
dlam_fab  = sigma_lam * randn(N_rings, 1);
df_fab_Hz = -(c_light / lambda_0^2) * dlam_fab;

% Q variation: scale alpha and tau_c independently
alpha_scale = 1 + sigma_Q_rel * randn(N_rings, 1);
alpha_scale(alpha_scale < 0.5) = 0.5;
tau_c_scale = 1 + sigma_Q_rel * randn(N_rings, 1);
tau_c_scale(tau_c_scale < 0.5) = 0.5;

fprintf('Generated %d rings\n', N_rings);
fprintf('Resonance offset std: %.1f GHz\n', std(df_fab_Hz)/1e9);
fprintf('Alpha scale: mean=%.2f, std=%.2f\n', mean(alpha_scale), std(alpha_scale));
fprintf('Tau_c scale: mean=%.2f, std=%.2f\n', mean(tau_c_scale), std(tau_c_scale));

%% Compute nominal normalized parameters
pn = compute_params(alpha_dBcm_nom, tau_c_nom, omega_0, c_light, n_Si, ...
    beta_TPA, V_TPA, Gamma_FCA, sigma_FCA, m_ring, c_p, ...
    Gamma_th, tau_FC, tau_th, dW0_dVpn, tau_RC);

linewidth_GHz = pn.gamma_0 / (2*pi*1e9);
Q_nom = omega_0 / (2*pn.gamma_0);
fprintf('Nominal: linewidth=%.1f GHz, Q=%.0f\n\n', linewidth_GHz, Q_nom);

%% Pre-compute per-ring parameters
ring_alpha = alpha_dBcm_nom * alpha_scale;
ring_tau_c = tau_c_nom * tau_c_scale;

ring_gamma_0     = zeros(N_rings,1);
ring_mu_raw      = zeros(N_rings,1);
ring_mu_n        = zeros(N_rings,1);
ring_glin_n      = zeros(N_rings,1);
ring_cTPA_n      = zeros(N_rings,1);
ring_cFCA_n      = zeros(N_rings,1);
ring_cTh_n       = zeros(N_rings,1);
ring_iTFC_n      = zeros(N_rings,1);
ring_iTTH_n      = zeros(N_rings,1);
ring_dw_NL_scale = zeros(N_rings,1);
ring_dw_mod_n    = zeros(N_rings,1);
ring_iTRC_n      = zeros(N_rings,1);

for r = 1:N_rings
    pr = compute_params(ring_alpha(r), ring_tau_c(r), omega_0, c_light, n_Si, ...
        beta_TPA, V_TPA, Gamma_FCA, sigma_FCA, m_ring, c_p, ...
        Gamma_th, tau_FC, tau_th, dW0_dVpn, tau_RC);
    ring_gamma_0(r)     = pr.gamma_0;
    ring_mu_raw(r)      = pr.mu_raw;
    ring_mu_n(r)        = pr.mu_n;
    ring_glin_n(r)      = pr.glin_n;
    ring_cTPA_n(r)      = pr.cTPA_n;
    ring_cFCA_n(r)      = pr.cFCA_n;
    ring_cTh_n(r)       = pr.cTh_n;
    ring_iTFC_n(r)      = pr.iTFC_n;
    ring_iTTH_n(r)      = pr.iTTH_n;
    ring_dw_NL_scale(r) = pr.dw_NL_scale;
    ring_dw_mod_n(r)    = pr.dw_mod_n;
    ring_iTRC_n(r)      = pr.iTRC_n;
end

ring_Q = omega_0 ./ (2*ring_gamma_0);
fprintf('Ring Q stats: mean=%.0f, std=%.0f (%.1f%% rel)\n', ...
    mean(ring_Q), std(ring_Q), 100*std(ring_Q)/mean(ring_Q));

%% Start parallel pool
pool = gcp('nocreate');
if isempty(pool), pool = parpool('local'); end
fprintf('Workers: %d\n\n', pool.NumWorkers);

%% Loop over configurations
delta_f_GHz = -80:1:80;
N_f = numel(delta_f_GHz);
prbs_bits = gen_prbs7();
N_bits = length(prbs_bits);
N_warmup = 200;

for ic = 1:N_configs
    P_in_mW   = configs(ic, 1);
    dr_Gbps   = configs(ic, 2);
    data_rate  = dr_Gbps * 1e9;
    P_in_W     = P_in_mW * 1e-3;
    E_in       = sqrt(P_in_W);

    fprintf('============================================================\n');
    fprintf('  Config %d/%d: P_in = %.2f mW, data rate = %.1f Gbps\n', ...
        ic, N_configs, P_in_mW, dr_Gbps);
    fprintf('============================================================\n');

    % build the voltage waveform for this data rate
    T_bit       = 1 / data_rate;
    dt_real     = 1e-12;
    M_bit       = round(T_bit / dt_real);
    dtau_nom    = pn.gamma_0 * dt_real;
    total_bits  = N_warmup + N_bits;
    total_steps = total_bits * M_bit;

    V_high = V_bias + V_pp/2;
    V_low  = V_bias - V_pp/2;

    V_wave = zeros(total_steps, 1);
    V_wave(1 : N_warmup*M_bit) = V_bias;
    for b = 1:N_bits
        if prbs_bits(b) == 1, v = V_low; else, v = V_high; end
        i1 = N_warmup*M_bit + (b-1)*M_bit + 1;
        i2 = N_warmup*M_bit + b*M_bit;
        V_wave(i1:i2) = v;
    end

    fprintf('  Bit period: %.0f ps, steps/bit: %d, total steps: %d\n', ...
        T_bit*1e12, M_bit, total_steps);

    %% Part 1: Nominal detuning sweep (no variation)
    OMA_t = zeros(N_f,1); ER_t = zeros(N_f,1); IL_t = zeros(N_f,1);

    fprintf('  Nominal sweep (%d points)...\n', N_f);
    t0 = tic;

    parfor iw = 1:N_f
        dw_ext  = 2*pi * delta_f_GHz(iw) * 1e9;
        dw_n    = dw_ext / pn.gamma_0;
        omega_p = omega_0 + dw_ext;
        cFCg_l  = (Gamma_FCA * c_light^2 * beta_TPA / ...
                   (2 * hbar * omega_p * V_FCA^2 * n_Si^2)) / pn.gamma_0;

        [~,~,~, ot,et,it] = sim_mrr_modulator(E_in, V_wave, dw_n, dtau_nom, ...
            total_steps, M_bit, N_warmup, N_bits, prbs_bits, ...
            pn.mu_n, pn.glin_n, pn.cTPA_n, pn.cFCA_n, cFCg_l, ...
            pn.cTh_n, pn.gamma_0, pn.iTFC_n, pn.iTTH_n, dn_dN, dn_dT, ...
            pn.dw_NL_scale, pn.dw_mod_n, pn.iTRC_n, P_in_W, pn.mu_raw);

        OMA_t(iw) = ot; ER_t(iw) = et; IL_t(iw) = it;
    end
    fprintf('  Done in %.1f s\n', toc(t0));

    % nominal detuning sweep plot
    fig1 = figure('Position',[40 80 900 500],'Color','w');
    yyaxis left
    plot(delta_f_GHz, ER_t, 'k-','LineWidth',1.5); hold on;
    plot(delta_f_GHz, IL_t, 'r-','LineWidth',1.5);
    ylabel('ER and IL  [dB]','FontSize',16,'FontWeight','bold');
    ax = gca; ax.YAxis(1).Color = [0.8 0 0];
    ylMax = max([max(ER_t),max(IL_t)])*1.3; if ylMax<1,ylMax=5;end
    ylim([0 ylMax]);
    yyaxis right
    plot(delta_f_GHz, OMA_t, 'b-','LineWidth',2);
    ylabel('Normalised OMA','FontSize',16,'FontWeight','bold');
    ax.YAxis(2).Color = [0 0 0.8];
    ylim([0 max(OMA_t)*1.3+1e-4]);
    xlabel('\Delta\omega / 2\pi   [GHz]','FontSize',16,'FontWeight','bold');
    leg = legend('ER','IL','Norm. OMA','Location','best');
    set(leg,'FontSize',13,'FontWeight','bold');
    set(gca,'FontSize',14,'FontWeight','bold');
    ax.XAxis.FontWeight = 'bold';
    ax.YAxis(1).FontWeight = 'bold';
    ax.YAxis(2).FontWeight = 'bold';
    grid on;
    fname1 = sprintf('paper/nominal_Pin%.0f_%.0fGbps.png', P_in_mW*1000, dr_Gbps*10);
    print(fig1, fname1, '-dpng', '-r300');
    fprintf('  Saved %s\n', fname1);
    close(fig1);

    %% Part 2: 2000-ring experiment with Q variation

    % corrected: delta_w=0, but each ring has its own Q
    OMA_corrected = zeros(N_rings, 1);
    ER_corrected  = zeros(N_rings, 1);
    IL_corrected  = zeros(N_rings, 1);

    fprintf('  Running %d corrected rings (delta_w=0, Q varied)...\n', N_rings);
    t0 = tic;

    parfor r = 1:N_rings
        dtau_r  = ring_gamma_0(r) * dt_real;
        omega_p = omega_0;
        cFCg_r  = (Gamma_FCA * c_light^2 * beta_TPA / ...
                   (2 * hbar * omega_p * V_FCA^2 * n_Si^2)) / ring_gamma_0(r);

        [~,~,~, ot,et,it] = sim_mrr_modulator(E_in, V_wave, 0, dtau_r, ...
            total_steps, M_bit, N_warmup, N_bits, prbs_bits, ...
            ring_mu_n(r), ring_glin_n(r), ring_cTPA_n(r), ring_cFCA_n(r), cFCg_r, ...
            ring_cTh_n(r), ring_gamma_0(r), ring_iTFC_n(r), ring_iTTH_n(r), dn_dN, dn_dT, ...
            ring_dw_NL_scale(r), ring_dw_mod_n(r), ring_iTRC_n(r), P_in_W, ring_mu_raw(r));

        OMA_corrected(r) = ot;
        ER_corrected(r)  = et;
        IL_corrected(r)  = it;
    end
    fprintf('  Done in %.1f s\n', toc(t0));

    % as-fabricated: each ring has its own delta_w AND its own Q
    OMA_asfab = zeros(N_rings, 1);
    ER_asfab  = zeros(N_rings, 1);
    IL_asfab  = zeros(N_rings, 1);

    fprintf('  Running %d as-fabricated rings (delta_w + Q varied)...\n', N_rings);
    t0 = tic;

    parfor r = 1:N_rings
        dtau_r  = ring_gamma_0(r) * dt_real;
        dw_ext  = 2*pi * df_fab_Hz(r);
        dw_n    = dw_ext / ring_gamma_0(r);
        omega_p = omega_0 + dw_ext;
        cFCg_r  = (Gamma_FCA * c_light^2 * beta_TPA / ...
                   (2 * hbar * omega_p * V_FCA^2 * n_Si^2)) / ring_gamma_0(r);

        [~,~,~, ot,et,it] = sim_mrr_modulator(E_in, V_wave, dw_n, dtau_r, ...
            total_steps, M_bit, N_warmup, N_bits, prbs_bits, ...
            ring_mu_n(r), ring_glin_n(r), ring_cTPA_n(r), ring_cFCA_n(r), cFCg_r, ...
            ring_cTh_n(r), ring_gamma_0(r), ring_iTFC_n(r), ring_iTTH_n(r), dn_dN, dn_dT, ...
            ring_dw_NL_scale(r), ring_dw_mod_n(r), ring_iTRC_n(r), P_in_W, ring_mu_raw(r));

        OMA_asfab(r) = ot;
        ER_asfab(r)  = et;
        IL_asfab(r)  = it;
    end
    fprintf('  Done in %.1f s\n', toc(t0));

    % stats
    n_oma_improved = sum(OMA_corrected > OMA_asfab);
    n_il_hurts     = sum(IL_corrected > IL_asfab);
    n_er_improved  = sum(ER_corrected > ER_asfab);
    fprintf('\n  Results:\n');
    fprintf('  Corrected mean OMA:    %.4f\n', mean(OMA_corrected));
    fprintf('  Corrected median OMA:  %.4f\n', median(OMA_corrected));
    fprintf('  Corrected mean ER:     %.4f dB\n', mean(ER_corrected));
    fprintf('  Corrected mean IL:     %.4f dB\n', mean(IL_corrected));
    fprintf('  As-fab mean OMA:       %.4f\n', mean(OMA_asfab));
    fprintf('  As-fab median OMA:     %.4f\n', median(OMA_asfab));
    fprintf('  As-fab mean ER:        %.4f dB\n', mean(ER_asfab));
    fprintf('  As-fab mean IL:        %.4f dB\n', mean(IL_asfab));
    fprintf('  OMA improved by corr:  %d / %d (%.1f%%)\n', ...
        n_oma_improved, N_rings, 100*n_oma_improved/N_rings);
    fprintf('  IL worse after corr:   %d / %d (%.1f%%)\n', ...
        n_il_hurts, N_rings, 100*n_il_hurts/N_rings);
    fprintf('  ER improved by corr:   %d / %d (%.1f%%)\n', ...
        n_er_improved, N_rings, 100*n_er_improved/N_rings);

    % OMA vs ring index
    fig2 = figure('Position', [80, 80, 1800, 700], 'Color', 'w');
    plot(1:N_rings, OMA_asfab, 'r-', 'LineWidth', 0.5); hold on;
    plot(1:N_rings, OMA_corrected, 'b-', 'LineWidth', 0.8);
    xlabel('Ring Index', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Normalised OMA', 'FontSize', 16, 'FontWeight', 'bold');
    leg = legend('As-fabricated (\Delta\omega + Q varied)', ...
           sprintf('Corrected (\\delta\\omega=0, Q varied) mean=%.3f', mean(OMA_corrected)), ...
           'Location', 'best', 'FontSize', 13);
    set(leg, 'FontWeight', 'bold');
    ylim([0 max([max(OMA_asfab), max(OMA_corrected)])*1.2]);
    ax2 = gca(fig2);
    ax2.FontSize = 14;
    ax2.FontWeight = 'bold';
    ax2.XAxis.FontWeight = 'bold';
    ax2.YAxis.FontWeight = 'bold';
    grid on;
    fname2 = sprintf('paper/rings_OMA_Pin%.0f_%.0fGbps.png', P_in_mW*1000, dr_Gbps*10);
    print(fig2, fname2, '-dpng', '-r300');
    fprintf('  Saved %s\n', fname2);
    close(fig2);

    % IL vs ring index
    fig3 = figure('Position', [80, 80, 1100, 500], 'Color', 'w');
    plot(1:N_rings, IL_asfab, 'r-', 'LineWidth', 0.5); hold on;
    plot(1:N_rings, IL_corrected, 'b-', 'LineWidth', 0.8);
    xlabel('Ring Index', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Insertion Loss  [dB]', 'FontSize', 16, 'FontWeight', 'bold');
    legend('As-fabricated (\Delta\omega + Q varied)', ...
           sprintf('Corrected (\\delta\\omega=0, Q varied) mean=%.2f dB', mean(IL_corrected)), ...
           'Location', 'best', 'FontSize', 11);
    ylim([0 max([max(IL_asfab), max(IL_corrected)])*1.2 + 0.1]);
    set(gca(fig3), 'FontSize', 14); grid on;
    fname3 = sprintf('paper/rings_IL_Pin%.0f_%.0fGbps.png', P_in_mW*1000, dr_Gbps*10);
    print(fig3, fname3, '-dpng', '-r300');
    fprintf('  Saved %s\n\n', fname3);
    close(fig3);

    % save this config
    results(ic).P_in_mW        = P_in_mW;
    results(ic).dr_Gbps        = dr_Gbps;
    results(ic).delta_f_GHz    = delta_f_GHz;
    results(ic).OMA_t          = OMA_t;
    results(ic).ER_t           = ER_t;
    results(ic).IL_t           = IL_t;
    results(ic).OMA_corrected  = OMA_corrected;
    results(ic).ER_corrected   = ER_corrected;
    results(ic).IL_corrected   = IL_corrected;
    results(ic).OMA_asfab      = OMA_asfab;
    results(ic).ER_asfab       = ER_asfab;
    results(ic).IL_asfab       = IL_asfab;
    results(ic).n_oma_improved = n_oma_improved;
    results(ic).n_il_hurts     = n_il_hurts;
    results(ic).n_er_improved  = n_er_improved;
end

%% Summary table
fprintf('\n\n');
sep = repmat('=', 1, 140);
div = repmat('-', 1, 140);
fprintf('%s\n', sep);
fprintf('  SUMMARY TABLE — %d powers x %d rates = %d configs\n', ...
    length(powers), length(datarates), N_configs);
fprintf('%s\n', sep);
fprintf('%-6s %-7s | %-10s %-10s %-10s | %-10s %-10s %-10s | %-18s %-18s %-18s\n', ...
    'Pin', 'Rate', 'As-fab', 'As-fab', 'As-fab', ...
    'Corr', 'Corr', 'Corr', ...
    'OMA impr', 'IL worse', 'ER impr');
fprintf('%-6s %-7s | %-10s %-10s %-10s | %-10s %-10s %-10s | %-18s %-18s %-18s\n', ...
    '(mW)', '(Gbps)', 'mean OMA', 'mean IL', 'mean ER', ...
    'mean OMA', 'mean IL', 'mean ER', ...
    'by corr', 'after corr', 'by corr');
fprintf('%s\n', div);
for ic = 1:N_configs
    fprintf('%-6.2f %-7.1f | %-10.4f %-10.2f %-10.2f | %-10.4f %-10.2f %-10.2f | %4d/%d (%5.1f%%)   %4d/%d (%5.1f%%)   %4d/%d (%5.1f%%)\n', ...
        results(ic).P_in_mW, results(ic).dr_Gbps, ...
        mean(results(ic).OMA_asfab), mean(results(ic).IL_asfab), mean(results(ic).ER_asfab), ...
        mean(results(ic).OMA_corrected), mean(results(ic).IL_corrected), mean(results(ic).ER_corrected), ...
        results(ic).n_oma_improved, N_rings, 100*results(ic).n_oma_improved/N_rings, ...
        results(ic).n_il_hurts, N_rings, 100*results(ic).n_il_hurts/N_rings, ...
        results(ic).n_er_improved, N_rings, 100*results(ic).n_er_improved/N_rings);
    if mod(ic, length(datarates)) == 0 && ic < N_configs
        fprintf('%s\n', div);
    end
end
fprintf('%s\n', sep);

%% Save everything
save('multi_config_results_full.mat', 'results', 'configs', ...
    'powers', 'datarates', ...
    'df_fab_Hz', 'dlam_fab', 'ring_Q', 'alpha_scale', 'tau_c_scale', ...
    'N_rings', 'sigma_lam', 'sigma_Q_rel');

fprintf('\nAll results saved to multi_config_results_full.mat\n');
fprintf('Done.\n');

%% local function (has to sit at the end in MATLAB)
function p = compute_params(alpha_dBcm, tau_c, omega_0, c_light, n_Si, ...
        beta_TPA, V_TPA, Gamma_FCA, sigma_FCA, m_ring, c_p, ...
        Gamma_th, tau_FC, tau_th, dW0_dVpn, tau_RC)
    alpha_Npm  = alpha_dBcm * 100 * log(10)/10;
    gamma_lin  = c_light * alpha_Npm / n_Si;
    gamma_coup = 2 / tau_c;
    gamma_0    = gamma_lin + gamma_coup;
    mu_raw     = sqrt(2 / tau_c);

    p.gamma_0     = gamma_0;
    p.mu_raw      = mu_raw;
    p.mu_n        = mu_raw / gamma_0;
    p.glin_n      = gamma_lin / gamma_0;
    p.cTPA_n      = (beta_TPA * c_light^2 / (n_Si^2 * V_TPA)) / gamma_0;
    p.cFCA_n      = (Gamma_FCA * sigma_FCA * c_light / (2*n_Si)) / gamma_0;
    p.cTh_n       = Gamma_th / (m_ring * c_p * gamma_0);
    p.iTFC_n      = (1/tau_FC) / gamma_0;
    p.iTTH_n      = (1/tau_th) / gamma_0;
    p.dw_NL_scale = omega_0 / (n_Si * gamma_0);
    p.dw_mod_n    = -dW0_dVpn / gamma_0;
    p.iTRC_n      = (1/tau_RC) / gamma_0;
end