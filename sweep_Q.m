% sweep_Q.m
%
% Outer sweep over Q = 3000:500:15000.
% Inner sweep: 3 powers x 6 data rates = 18 configs per Q.
%
% Each Q value gets its own Q_XXXX/ folder with:
%   - nominal OMA/IL/ER detuning-sweep figures
%   - as-fabricated vs corrected ring-index figures for OMA/IL/ER
%   - a summary table, printed and saved as .txt
%
% Needs sim_mrr_modulator.m, mod_rk4_step.m, gen_prbs7.m

clear; clc; close all;

%% Physical constants
c_light = 2.99792458e8;
hbar    = 1.054571817e-34;

%% MRR parameters (Giron Castro 2024, Table 1)
m_ring     = 1.2e-11;
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
sigma_lam   = 0.65e-9;
sigma_Q_rel = 0.10;
N_rings     = 2000;

%% Sweep parameters
Q_values  = 3000:500:15000;
powers    = [2.00, 3.00, 4.00];
datarates = [2.5, 5.0, 7.5, 10.0, 12.5, 15.0];

%% Build config list
configs = zeros(length(powers)*length(datarates), 2);
idx = 0;
for ip = 1:length(powers)
    for id = 1:length(datarates)
        idx = idx + 1;
        configs(idx,:) = [powers(ip), datarates(id)];
    end
end
N_configs = size(configs, 1);
fprintf('Configs per Q: %d (%d powers x %d rates)\n', N_configs, length(powers), length(datarates));
fprintf('Q values: %d to %d in steps of 500 (%d total)\n', Q_values(1), Q_values(end), length(Q_values));
fprintf('Total simulations: %d Q x %d configs = %d\n\n', length(Q_values), N_configs, length(Q_values)*N_configs);

%% Generate fabrication offsets (same across all Q values)
rng(42);
dlam_fab  = sigma_lam * randn(N_rings, 1);
df_fab_Hz = -(c_light / lambda_0^2) * dlam_fab;

alpha_scale = 1 + sigma_Q_rel * randn(N_rings, 1);
alpha_scale(alpha_scale < 0.5) = 0.5;
tau_c_scale = 1 + sigma_Q_rel * randn(N_rings, 1);
tau_c_scale(tau_c_scale < 0.5) = 0.5;

%% PRBS and detuning sweep setup
prbs_bits = gen_prbs7();
N_bits    = length(prbs_bits);
N_warmup  = 200;
delta_f_GHz = -80:1:80;
N_f = numel(delta_f_GHz);
dt_real = 1e-12;

V_high = V_bias + V_pp/2;
V_low  = V_bias - V_pp/2;

%% Start parallel pool
pool = gcp('nocreate');
if isempty(pool), pool = parpool('local'); end
fprintf('Workers: %d\n\n', pool.NumWorkers);

%% Create top-level output directory
if ~exist('Q_sweep_results', 'dir'), mkdir('Q_sweep_results'); end

%% Grand summary storage
grand_results = struct();
grand_idx = 0;

%% outer loop over Q values
for iQ = 1:length(Q_values)
    Q_target = Q_values(iQ);

    fprintf('\n################################################################\n');
    fprintf('  Q = %d  (%d / %d)\n', Q_target, iQ, length(Q_values));
    fprintf('################################################################\n');

    % back out tau_c from the target Q.
    % Q = omega_0/(2*gamma_0), gamma_0 = gamma_lin + gamma_coup, and
    % gamma_lin comes from the nominal alpha.
    alpha_Npm_nom = alpha_dBcm_nom * 100 * log(10)/10;
    gamma_lin_nom = c_light * alpha_Npm_nom / n_Si;
    % gamma_0 = omega_0/(2*Q)
    gamma_0_target = omega_0 / (2 * Q_target);
    gamma_coup_target = gamma_0_target - gamma_lin_nom;
    if gamma_coup_target <= 0
        fprintf('  WARNING: Q=%d requires negative coupling loss, skipping.\n', Q_target);
        continue;
    end
    tau_c_nom = 2 / gamma_coup_target;

    fprintf('  tau_c = %.2f ps, gamma_0 = %.3e rad/s\n', tau_c_nom*1e12, gamma_0_target);

    % make the Q folder
    Q_folder = sprintf('Q_sweep_results/Q_%05d', Q_target);
    if ~exist(Q_folder, 'dir'), mkdir(Q_folder); end

    %% Compute nominal normalized parameters for this Q
    pn = compute_params(alpha_dBcm_nom, tau_c_nom, omega_0, c_light, n_Si, ...
        beta_TPA, V_TPA, Gamma_FCA, sigma_FCA, m_ring, c_p, ...
        Gamma_th, tau_FC, tau_th, dW0_dVpn, tau_RC);

    Q_nom_check = omega_0 / (2*pn.gamma_0);
    fprintf('  Nominal Q check: %.0f\n', Q_nom_check);

    %% Pre-compute per-ring parameters (Q variation around this nominal)
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
    fprintf('  Ring Q stats: mean=%.0f, std=%.0f\n', mean(ring_Q), std(ring_Q));

    %% Local results storage for this Q
    results_Q = struct();

    %% inner loop over configs
    for ic = 1:N_configs
        P_in_mW   = configs(ic, 1);
        dr_Gbps   = configs(ic, 2);
        data_rate  = dr_Gbps * 1e9;
        P_in_W     = P_in_mW * 1e-3;
        E_in       = sqrt(P_in_W);

        fprintf('  Config %d/%d: P_in=%.2f mW, rate=%.1f Gbps\n', ic, N_configs, P_in_mW, dr_Gbps);

        % voltage waveform
        T_bit       = 1 / data_rate;
        M_bit       = round(T_bit / dt_real);
        dtau_nom    = pn.gamma_0 * dt_real;
        total_bits  = N_warmup + N_bits;
        total_steps = total_bits * M_bit;

        V_wave = zeros(total_steps, 1);
        V_wave(1 : N_warmup*M_bit) = V_bias;
        for b = 1:N_bits
            if prbs_bits(b) == 1, v = V_low; else, v = V_high; end
            i1 = N_warmup*M_bit + (b-1)*M_bit + 1;
            i2 = N_warmup*M_bit + b*M_bit;
            V_wave(i1:i2) = v;
        end

        %% Part 1: Nominal detuning sweep
        OMA_t = zeros(N_f,1); ER_t = zeros(N_f,1); IL_t = zeros(N_f,1);

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

        %% Nominal detuning sweep figure: ER/IL (left) + OMA (right)
        fig = figure('Position',[40 80 900 500],'Color','w','Visible','off');
        yyaxis left
        plot(delta_f_GHz, ER_t, 'k-', 'LineWidth', 1.5); hold on;
        plot(delta_f_GHz, IL_t, 'r-', 'LineWidth', 1.5);
        ylabel('ER and IL  [dB]', 'FontSize', 18, 'FontWeight', 'bold');
        ax = gca; ax.YAxis(1).Color = [0.8 0 0];
        ylMax = max([max(ER_t), max(IL_t)])*1.3; if ylMax < 1, ylMax = 5; end
        ylim([0 ylMax]);

        yyaxis right
        plot(delta_f_GHz, OMA_t, 'b-', 'LineWidth', 2);
        ylabel('Normalised OMA', 'FontSize', 18, 'FontWeight', 'bold');
        ax.YAxis(2).Color = [0 0 0.8];
        ylim([0 max(OMA_t)*1.3 + 1e-4]);

        xlabel('\Delta\omega / 2\pi   [GHz]', 'FontSize', 18, 'FontWeight', 'bold');
        title(sprintf('Through port — Q=%d, P_{in}=%.2fmW, %.1fGbps', Q_target, P_in_mW, dr_Gbps), ...
            'FontSize', 16, 'FontWeight', 'bold');
        leg = legend('ER', 'IL', 'Norm. OMA', 'Location', 'best');
        set(leg, 'FontSize', 14, 'FontWeight', 'bold');
        set(gca, 'FontSize', 16, 'FontWeight', 'bold');
        ax.XAxis.FontSize = 16; ax.XAxis.FontWeight = 'bold';
        ax.YAxis(1).FontSize = 16; ax.YAxis(1).FontWeight = 'bold';
        ax.YAxis(2).FontSize = 16; ax.YAxis(2).FontWeight = 'bold';
        grid on;
        fname = sprintf('%s/nominal_Pin%.0f_%.0fGbps.png', Q_folder, P_in_mW*1000, dr_Gbps*10);
        print(fig, fname, '-dpng', '-r300'); close(fig);

        %% Part 2: 2000-ring experiment

        % corrected rings: delta_w = 0, Q varied
        OMA_corrected = zeros(N_rings, 1);
        ER_corrected  = zeros(N_rings, 1);
        IL_corrected  = zeros(N_rings, 1);

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

        % as-fabricated rings: delta_w and Q both varied
        OMA_asfab = zeros(N_rings, 1);
        ER_asfab  = zeros(N_rings, 1);
        IL_asfab  = zeros(N_rings, 1);

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

        % stats
        n_oma_improved = sum(OMA_corrected > OMA_asfab);
        n_il_hurts     = sum(IL_corrected > IL_asfab);
        n_er_improved  = sum(ER_corrected > ER_asfab);

        %% Ring-index figure: OMA as-fab vs corrected
        fig = figure('Position', [80, 80, 1800, 700], 'Color', 'w', 'Visible', 'off');
        plot(1:N_rings, OMA_asfab, 'r-', 'LineWidth', 0.5); hold on;
        plot(1:N_rings, OMA_corrected, 'b-', 'LineWidth', 0.8);
        xlabel('Ring Index', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('Normalised OMA', 'FontSize', 18, 'FontWeight', 'bold');
        title(sprintf('Q=%d, P_{in}=%.1f mW, %.1f Gbps', Q_target, P_in_mW, dr_Gbps), ...
            'FontSize', 16, 'FontWeight', 'bold');
        leg = legend('As-fabricated (\Delta\omega + Q varied)', ...
               sprintf('Corrected (\\delta\\omega=0, Q varied) mean=%.3f', mean(OMA_corrected)), ...
               'Location', 'best', 'FontSize', 14);
        set(leg, 'FontWeight', 'bold');
        ylim([0 max([max(OMA_asfab), max(OMA_corrected)])*1.2]);
        set(gca, 'FontSize', 16, 'FontWeight', 'bold');
        ax2 = gca;
        ax2.XAxis.FontSize = 16; ax2.XAxis.FontWeight = 'bold';
        ax2.YAxis.FontSize = 16; ax2.YAxis.FontWeight = 'bold';
        grid on;
        fname = sprintf('%s/rings_OMA_Pin%.0f_%.0fGbps.png', Q_folder, P_in_mW*1000, dr_Gbps*10);
        print(fig, fname, '-dpng', '-r300'); close(fig);

        %% Ring-index figure: ER as-fab vs corrected
        fig = figure('Position', [80, 80, 1800, 700], 'Color', 'w', 'Visible', 'off');
        plot(1:N_rings, ER_asfab, 'r-', 'LineWidth', 0.5); hold on;
        plot(1:N_rings, ER_corrected, 'b-', 'LineWidth', 0.8);
        xlabel('Ring Index', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('Extinction Ratio [dB]', 'FontSize', 18, 'FontWeight', 'bold');
        title(sprintf('Q=%d, P_{in}=%.1f mW, %.1f Gbps', Q_target, P_in_mW, dr_Gbps), ...
            'FontSize', 16, 'FontWeight', 'bold');
        legend('As-fabricated (\Delta\omega + Q varied)', ...
               sprintf('Corrected (\\delta\\omega=0, Q varied) mean=%.2f dB', mean(ER_corrected)), ...
               'Location', 'best', 'FontSize', 14, 'FontWeight', 'bold');
        ylim([0 max([max(ER_asfab), max(ER_corrected)])*1.2 + 0.1]);
        set(gca, 'FontSize', 16, 'FontWeight', 'bold');
        ax3 = gca;
        ax3.XAxis.FontSize = 16; ax3.XAxis.FontWeight = 'bold';
        ax3.YAxis.FontSize = 16; ax3.YAxis.FontWeight = 'bold';
        grid on;
        fname = sprintf('%s/rings_ER_Pin%.0f_%.0fGbps.png', Q_folder, P_in_mW*1000, dr_Gbps*10);
        print(fig, fname, '-dpng', '-r300'); close(fig);

        %% Ring-index figure: IL as-fab vs corrected
        fig = figure('Position', [80, 80, 1800, 700], 'Color', 'w', 'Visible', 'off');
        plot(1:N_rings, IL_asfab, 'r-', 'LineWidth', 0.5); hold on;
        plot(1:N_rings, IL_corrected, 'b-', 'LineWidth', 0.8);
        xlabel('Ring Index', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('Insertion Loss [dB]', 'FontSize', 18, 'FontWeight', 'bold');
        title(sprintf('Q=%d, P_{in}=%.1f mW, %.1f Gbps', Q_target, P_in_mW, dr_Gbps), ...
            'FontSize', 16, 'FontWeight', 'bold');
        legend('As-fabricated (\Delta\omega + Q varied)', ...
               sprintf('Corrected (\\delta\\omega=0, Q varied) mean=%.2f dB', mean(IL_corrected)), ...
               'Location', 'best', 'FontSize', 14, 'FontWeight', 'bold');
        ylim([0 max([max(IL_asfab), max(IL_corrected)])*1.2 + 0.1]);
        set(gca, 'FontSize', 16, 'FontWeight', 'bold');
        ax4 = gca;
        ax4.XAxis.FontSize = 16; ax4.XAxis.FontWeight = 'bold';
        ax4.YAxis.FontSize = 16; ax4.YAxis.FontWeight = 'bold';
        grid on;
        fname = sprintf('%s/rings_IL_Pin%.0f_%.0fGbps.png', Q_folder, P_in_mW*1000, dr_Gbps*10);
        print(fig, fname, '-dpng', '-r300'); close(fig);

        %% Store results
        results_Q(ic).P_in_mW       = P_in_mW;
        results_Q(ic).dr_Gbps       = dr_Gbps;
        results_Q(ic).OMA_corrected = OMA_corrected;
        results_Q(ic).ER_corrected  = ER_corrected;
        results_Q(ic).IL_corrected  = IL_corrected;
        results_Q(ic).OMA_asfab     = OMA_asfab;
        results_Q(ic).ER_asfab      = ER_asfab;
        results_Q(ic).IL_asfab      = IL_asfab;
        results_Q(ic).n_oma_improved = n_oma_improved;
        results_Q(ic).n_il_hurts     = n_il_hurts;
        results_Q(ic).n_er_improved  = n_er_improved;

        grand_idx = grand_idx + 1;
        grand_results(grand_idx).Q_target      = Q_target;
        grand_results(grand_idx).P_in_mW       = P_in_mW;
        grand_results(grand_idx).dr_Gbps       = dr_Gbps;
        grand_results(grand_idx).mean_OMA_asfab     = mean(OMA_asfab);
        grand_results(grand_idx).mean_IL_asfab      = mean(IL_asfab);
        grand_results(grand_idx).mean_ER_asfab      = mean(ER_asfab);
        grand_results(grand_idx).mean_OMA_corrected = mean(OMA_corrected);
        grand_results(grand_idx).mean_IL_corrected  = mean(IL_corrected);
        grand_results(grand_idx).mean_ER_corrected  = mean(ER_corrected);
        grand_results(grand_idx).n_oma_improved     = n_oma_improved;
        grand_results(grand_idx).n_il_hurts         = n_il_hurts;
        grand_results(grand_idx).n_er_improved      = n_er_improved;
    end

    %% Per-Q summary table (console + file)
    table_file = sprintf('%s/summary_table_Q%05d.txt', Q_folder, Q_target);
    fid = fopen(table_file, 'w');

    hdr = sprintf('SUMMARY TABLE — Q = %d | %d powers x %d rates = %d configs', ...
        Q_target, length(powers), length(datarates), N_configs);

    sep = repmat('=', 1, 160);
    dash = repmat('-', 1, 160);

    for out = [1, fid]
        fprintf(out, '\n%s\n', sep);
        fprintf(out, '  %s\n', hdr);
        fprintf(out, '%s\n', sep);
        fprintf(out, '%-6s %-6s | %-10s %-10s %-10s | %-10s %-10s %-10s | %-18s %-18s %-18s\n', ...
            'Pin', 'Rate', 'As-fab', 'As-fab', 'As-fab', 'Corr', 'Corr', 'Corr', 'OMA improved', 'IL worse', 'ER improved');
        fprintf(out, '%-6s %-6s | %-10s %-10s %-10s | %-10s %-10s %-10s | %-18s %-18s %-18s\n', ...
            '(mW)', '(Gbps)', 'mean OMA', 'mean IL', 'mean ER', 'mean OMA', 'mean IL', 'mean ER', 'by corr', 'after corr', 'by corr');
        fprintf(out, '%s\n', dash);
        for ic = 1:N_configs
            fprintf(out, '%-6.2f %-6.1f | %-10.4f %-10.2f %-10.2f | %-10.4f %-10.2f %-10.2f | %4d/%d (%5.1f%%) %4d/%d (%5.1f%%) %4d/%d (%5.1f%%)\n', ...
                results_Q(ic).P_in_mW, results_Q(ic).dr_Gbps, ...
                mean(results_Q(ic).OMA_asfab), mean(results_Q(ic).IL_asfab), mean(results_Q(ic).ER_asfab), ...
                mean(results_Q(ic).OMA_corrected), mean(results_Q(ic).IL_corrected), mean(results_Q(ic).ER_corrected), ...
                results_Q(ic).n_oma_improved, N_rings, 100*results_Q(ic).n_oma_improved/N_rings, ...
                results_Q(ic).n_il_hurts, N_rings, 100*results_Q(ic).n_il_hurts/N_rings, ...
                results_Q(ic).n_er_improved, N_rings, 100*results_Q(ic).n_er_improved/N_rings);
            if mod(ic, length(datarates)) == 0 && ic < N_configs
                fprintf(out, '%s\n', dash);
            end
        end
        fprintf(out, '%s\n', sep);
    end
    fclose(fid);
    fprintf('  Saved %s\n', table_file);

    % save per-Q .mat
    mat_file = sprintf('%s/results_Q%05d.mat', Q_folder, Q_target);
    save(mat_file, 'results_Q', 'Q_target', 'configs', 'powers', 'datarates', ...
        'ring_Q', 'N_rings');
    fprintf('  Saved %s\n', mat_file);
end

%% Grand summary table across all Q values
fprintf('\n\n');
sep = repmat('=', 1, 180);
dash = repmat('-', 1, 180);
fprintf('%s\n', sep);
fprintf('  GRAND SUMMARY — All Q values\n');
fprintf('%s\n', sep);
fprintf('%-7s %-6s %-6s | %-10s %-10s %-10s | %-10s %-10s %-10s | %-16s %-16s %-16s\n', ...
    'Q', 'Pin', 'Rate', 'As-fab', 'As-fab', 'As-fab', 'Corr', 'Corr', 'Corr', 'OMA impr', 'IL worse', 'ER impr');
fprintf('%-7s %-6s %-6s | %-10s %-10s %-10s | %-10s %-10s %-10s | %-16s %-16s %-16s\n', ...
    '', '(mW)', '(Gbps)', 'mean OMA', 'mean IL', 'mean ER', 'mean OMA', 'mean IL', 'mean ER', 'by corr', 'after corr', 'by corr');
fprintf('%s\n', dash);

prev_Q = 0;
for gi = 1:length(grand_results)
    if grand_results(gi).Q_target ~= prev_Q && prev_Q ~= 0
        fprintf('%s\n', dash);
    end
    prev_Q = grand_results(gi).Q_target;

    fprintf('%-7d %-6.2f %-6.1f | %-10.4f %-10.2f %-10.2f | %-10.4f %-10.2f %-10.2f | %4d/%d (%5.1f%%) %4d/%d (%5.1f%%) %4d/%d (%5.1f%%)\n', ...
        grand_results(gi).Q_target, grand_results(gi).P_in_mW, grand_results(gi).dr_Gbps, ...
        grand_results(gi).mean_OMA_asfab, grand_results(gi).mean_IL_asfab, grand_results(gi).mean_ER_asfab, ...
        grand_results(gi).mean_OMA_corrected, grand_results(gi).mean_IL_corrected, grand_results(gi).mean_ER_corrected, ...
        grand_results(gi).n_oma_improved, N_rings, 100*grand_results(gi).n_oma_improved/N_rings, ...
        grand_results(gi).n_il_hurts, N_rings, 100*grand_results(gi).n_il_hurts/N_rings, ...
        grand_results(gi).n_er_improved, N_rings, 100*grand_results(gi).n_er_improved/N_rings);
end
fprintf('%s\n', sep);

%% Save grand results
save('Q_sweep_results/grand_results.mat', 'grand_results', 'Q_values', ...
    'powers', 'datarates', 'N_rings', 'sigma_lam', 'sigma_Q_rel', ...
    'df_fab_Hz', 'dlam_fab', 'alpha_scale', 'tau_c_scale');

% also save grand table to file
fid = fopen('Q_sweep_results/grand_summary_table.txt', 'w');
fprintf(fid, '%s\n', sep);
fprintf(fid, '  GRAND SUMMARY — All Q values\n');
fprintf(fid, '%s\n', sep);
fprintf(fid, '%-7s %-6s %-6s | %-10s %-10s %-10s | %-10s %-10s %-10s | %-16s %-16s %-16s\n', ...
    'Q', 'Pin', 'Rate', 'As-fab', 'As-fab', 'As-fab', 'Corr', 'Corr', 'Corr', 'OMA impr', 'IL worse', 'ER impr');
fprintf(fid, '%-7s %-6s %-6s | %-10s %-10s %-10s | %-10s %-10s %-10s | %-16s %-16s %-16s\n', ...
    '', '(mW)', '(Gbps)', 'mean OMA', 'mean IL', 'mean ER', 'mean OMA', 'mean IL', 'mean ER', 'by corr', 'after corr', 'by corr');
fprintf(fid, '%s\n', dash);
prev_Q = 0;
for gi = 1:length(grand_results)
    if grand_results(gi).Q_target ~= prev_Q && prev_Q ~= 0
        fprintf(fid, '%s\n', dash);
    end
    prev_Q = grand_results(gi).Q_target;
    fprintf(fid, '%-7d %-6.2f %-6.1f | %-10.4f %-10.2f %-10.2f | %-10.4f %-10.2f %-10.2f | %4d/%d (%5.1f%%) %4d/%d (%5.1f%%) %4d/%d (%5.1f%%)\n', ...
        grand_results(gi).Q_target, grand_results(gi).P_in_mW, grand_results(gi).dr_Gbps, ...
        grand_results(gi).mean_OMA_asfab, grand_results(gi).mean_IL_asfab, grand_results(gi).mean_ER_asfab, ...
        grand_results(gi).mean_OMA_corrected, grand_results(gi).mean_IL_corrected, grand_results(gi).mean_ER_corrected, ...
        grand_results(gi).n_oma_improved, N_rings, 100*grand_results(gi).n_oma_improved/N_rings, ...
        grand_results(gi).n_il_hurts, N_rings, 100*grand_results(gi).n_il_hurts/N_rings, ...
        grand_results(gi).n_er_improved, N_rings, 100*grand_results(gi).n_er_improved/N_rings);
end
fprintf(fid, '%s\n', sep);
fclose(fid);

fprintf('\nGrand results saved to Q_sweep_results/grand_results.mat\n');
fprintf('Grand summary table saved to Q_sweep_results/grand_summary_table.txt\n');
fprintf('Total Q folders: %d\n', length(Q_values));
fprintf('Total figures generated: %d\n', length(Q_values) * N_configs * 4);
fprintf('Done.\n');

%% Local function
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