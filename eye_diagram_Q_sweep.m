% eye_diagram_Q_sweep.m
%
% Same through-port eye diagrams as eye_diagram.m, but swept over a set of
% target Q values. Through port is what the OIP 2026 metrics use.
%
% How the eye gets built:
%   - drive with a 1000-symbol PRBS
%   - grab the through-port response over the data window
%   - chop it into 250 non-overlapping groups of 4 symbols
%   - overlay the 250 four-symbol traces
%
% Ring model (same as three_mW.m):
%   the ring has a residual Q spread (loss/coupling) you can't tune out and a
%   resonance offset you can. "corrected" = offset tuned out, residual Q left
%   in. "as-fabricated" = offset still there. Same ring either way, they only
%   differ by the offset.
%
% Per config I pick one laser detuning off the corrected ring (the OMA
% optimum), then run both rings at that detuning and at zero. Four cases per
% config: {as-fab, corrected} x {zero, optimal detuning}.
%
% Q sweep: for each target Q I back out the coupling time tau_c from the
% nominal loss (same as sweep_Q.m). Every Q gets its own subfolder under
% "Q sweep eye diagram" with the four eye grids, a table, and a .mat.
%
% Needs mod_rk4_step.m and sim_mrr_eye.m. run_eye_case is a local function below.

clear; clc; close all;
out_root = 'Q sweep eye diagram';
if ~exist(out_root, 'dir'), mkdir(out_root); end

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

%% Fabrication variation (Parsons 2025; same model as three_mW.m)
% residual Q spread you can't correct, present in both rings:
res_alpha_scale = 1.10;    % +1 sigma residual loss (sigma_Q_rel = 0.10)
res_tau_c_scale = 1.10;    % +1 sigma residual coupling time
% resonance offset, only present in the as-fab ring:
fab_offset_GHz  = 20;      % representative offset [GHz]. 1 sigma is ~81 GHz;
                           % sign sets which flank; magnitude sets severity.

%% Q sweep
Q_values = [3000 6000 9000 12000 15000];

% nominal linear loss rate, needed to back out the coupling rate from Q
alpha_Npm_nom = alpha_dBcm_nom * 100 * log(10)/10;
gamma_lin_nom = c_light * alpha_Npm_nom / n_Si;

%% Eye configuration
powers      = [2.00, 3.00, 4.00];          % mW
datarates   = [2.5, 5.0, 10.0, 15.0];      % Gbps
df_step     = 2;        % laser detuning sweep step [GHz]
df_max_GHz  = 80;       % sweep span; covers the widest linewidth (Q=3000 ~ 32 GHz)
df_grid     = 0:df_step:df_max_GHz;   % fixed detuning search, same for every Q
N_sym       = 1000;     % PRBS symbols
prbs_order  = 10;       % PRBS order (10 gives period 1023, longer than the 1000 we use)
group_size  = 4;        % sequential symbols per superimposed group (per the recipe)
n_groups    = N_sym / group_size;          % 250 non-overlapping groups
n_tau_warm  = 5;        % warmup length in multiples of tau_th, toggling data
dt_real     = 1e-12;    % fixed 1 ps time step
save_individual = false;

assert(mod(N_sym, group_size) == 0, 'N_sym must be divisible by group_size');

%% PRBS drive (Q-independent)
prbs_bits = gen_prbs_len(N_sym, prbs_order);
V_high = V_bias + V_pp/2;
V_low  = V_bias - V_pp/2;

%% constants struct for run_eye_case
consts = struct('dt_real', dt_real, 'Gamma_FCA', Gamma_FCA, 'c_light', c_light, ...
    'beta_TPA', beta_TPA, 'hbar', hbar, 'V_FCA', V_FCA, 'n_Si', n_Si, ...
    'dn_dN', dn_dN, 'dn_dT', dn_dT);

%% Config list (power outer, rate inner)
nP = numel(powers);  nR = numel(datarates);
configs = zeros(nP*nR, 2);
ic = 0;
for ip = 1:nP
    for id = 1:nR
        ic = ic + 1;
        configs(ic,:) = [powers(ip), datarates(id)];
    end
end
N_configs = size(configs, 1);

%% Parallel pool
pool = gcp('nocreate');
if isempty(pool), pool = parpool('local'); end
addAttachedFiles(pool, {'sim_mrr_eye.m', 'mod_rk4_step.m'});
fprintf('Workers: %d\n', pool.NumWorkers);

%% outer loop over target Q
for iQ = 1:numel(Q_values)
    Q_target = Q_values(iQ);

    % derive tau_c from target Q using the nominal loss (Q = omega_0/(2 gamma_0))
    gamma_0_target    = omega_0 / (2 * Q_target);
    gamma_coup_target = gamma_0_target - gamma_lin_nom;
    if gamma_coup_target <= 0
        fprintf('Q=%d requires negative coupling loss, skipping.\n', Q_target);
        continue;
    end
    tau_c_nom = 2 / gamma_coup_target;

    % per-Q output folder
    Q_folder = sprintf('%s/Q %d', out_root, Q_target);
    if ~exist(Q_folder, 'dir'), mkdir(Q_folder); end

    % normalized params: nominal reference and the residual-Q ring at this Q
    pn = compute_params(alpha_dBcm_nom, tau_c_nom, omega_0, c_light, n_Si, ...
        beta_TPA, V_TPA, Gamma_FCA, sigma_FCA, m_ring, c_p, ...
        Gamma_th, tau_FC, tau_th, dW0_dVpn, tau_RC);
    pn_r = compute_params(alpha_dBcm_nom*res_alpha_scale, tau_c_nom*res_tau_c_scale, ...
        omega_0, c_light, n_Si, beta_TPA, V_TPA, Gamma_FCA, sigma_FCA, m_ring, c_p, ...
        Gamma_th, tau_FC, tau_th, dW0_dVpn, tau_RC);

    % detuning search is the fixed df_grid set above, same window for every Q
    lw_GHz = pn.gamma_0/(2*pi*1e9);

    fprintf('\n=== Q=%d:  tau_c=%.2f ps, linewidth=%.1f GHz, ring Q=%.0f, df sweep 0..%g GHz ===\n', ...
        Q_target, tau_c_nom*1e12, lw_GHz, omega_0/(2*pn_r.gamma_0), df_max_GHz);

    %% Storage: 4 cases (all on the residual-Q ring pn_r)
    eye_a0   = cell(N_configs,1);  % as-fab, zero
    eye_c0   = cell(N_configs,1);  % corrected, zero
    eye_aopt = cell(N_configs,1);  % as-fab, optimal
    eye_copt = cell(N_configs,1);  % corrected, optimal
    Mbit_arr = zeros(N_configs,1);
    df_opt   = zeros(N_configs,1); % design-optimal laser detuning (from corrected ring)
    oma_a0 = zeros(N_configs,1); oma_c0 = zeros(N_configs,1);
    oma_ao = zeros(N_configs,1); oma_co = zeros(N_configs,1);
    eo_a0  = zeros(N_configs,1); eo_c0  = zeros(N_configs,1);
    eo_ao  = zeros(N_configs,1); eo_co  = zeros(N_configs,1);

    parfor ci = 1:N_configs
        P_in_mW = configs(ci,1);
        dr_Gbps = configs(ci,2);
        P_in_W  = P_in_mW*1e-3;
        E_in    = sqrt(P_in_W);

        T_bit   = 1/(dr_Gbps*1e9);
        M_bit   = round(T_bit/dt_real);

        % thermal-conditioned warmup + PRBS data, both toggling
        N_warmup    = ceil(n_tau_warm * tau_th / T_bit);
        warm_bits   = prbs_bits(mod((0:N_warmup-1), N_sym) + 1);
        all_bits    = [warm_bits; prbs_bits];
        total_bits  = N_warmup + N_sym;
        total_steps = total_bits * M_bit;
        V_wave = zeros(total_steps,1);
        for b = 1:total_bits
            if all_bits(b)==1, v=V_low; else, v=V_high; end
            V_wave((b-1)*M_bit+1 : b*M_bit) = v;
        end

        % corrected, zero detuning (laser phi = 0, ring on resonance: delta = 0)
        [Pc0, oc0, ec0] = run_eye_case(E_in, V_wave, 0, total_steps, ...
            M_bit, N_warmup, N_sym, prbs_bits, P_in_W, pn_r, omega_0, consts);

        % design optimum from the corrected ring (sweep laser phi; delta = phi)
        bc = -inf; Pco = Pc0; dco = 0; eco = ec0;
        for q = 1:numel(df_grid)
            [Pt, om, eo] = run_eye_case(E_in, V_wave, df_grid(q), total_steps, ...
                M_bit, N_warmup, N_sym, prbs_bits, P_in_W, pn_r, omega_0, consts);
            if om > bc, bc = om; Pco = Pt; dco = df_grid(q); eco = eo; end
        end

        % as-fabricated: same ring, uncorrected offset shifts delta = phi - offset
        [Pa0, oa0, ea0] = run_eye_case(E_in, V_wave, 0 - fab_offset_GHz, total_steps, ...
            M_bit, N_warmup, N_sym, prbs_bits, P_in_W, pn_r, omega_0, consts);
        [Pao, oao, eao] = run_eye_case(E_in, V_wave, dco - fab_offset_GHz, total_steps, ...
            M_bit, N_warmup, N_sym, prbs_bits, P_in_W, pn_r, omega_0, consts);

        eye_c0{ci}=Pc0; eye_copt{ci}=Pco; eye_a0{ci}=Pa0; eye_aopt{ci}=Pao;
        Mbit_arr(ci)=M_bit; df_opt(ci)=dco;
        oma_c0(ci)=oc0; oma_co(ci)=bc; oma_a0(ci)=oa0; oma_ao(ci)=oao;
        eo_c0(ci)=ec0;  eo_co(ci)=eco; eo_a0(ci)=ea0;  eo_ao(ci)=eao;

        fprintf('  done: %.0f mW, %.1f Gbps -> design opt df = %g GHz\n', P_in_mW, dr_Gbps, dco);
    end

    %% Shared y-limits per power across all four cases
    ymax_cfg = zeros(N_configs,1);
    for ip = 1:nP
        sel = find(configs(:,1)==powers(ip));
        ym = 0;
        for k = sel.'
            ym = max([ym, max(eye_a0{k}), max(eye_c0{k}), max(eye_aopt{k}), max(eye_copt{k})]);
        end
        ymax_cfg(sel) = 1.05 * ym * 1e3;
    end

    %% Four grid figures (panels labeled by laser detuning), saved into Q_folder
    z = zeros(N_configs,1);
    plot_eye_grid(eye_c0,   configs, Mbit_arr, z,      group_size, n_groups, ymax_cfg, nP, nR, ...
        sprintf('Corrected, \\Delta\\omega = 0 GHz   (Q=%d)', Q_target), ...
        sprintf('%s/eye_grid_corr_df0.png', Q_folder));
    plot_eye_grid(eye_a0,   configs, Mbit_arr, z,      group_size, n_groups, ymax_cfg, nP, nR, ...
        sprintf('As-fabricated, \\Delta\\omega = 0 GHz   (Q=%d, uncorrected offset %g GHz)', Q_target, fab_offset_GHz), ...
        sprintf('%s/eye_grid_asfab_df0.png', Q_folder));
    plot_eye_grid(eye_copt, configs, Mbit_arr, df_opt, group_size, n_groups, ymax_cfg, nP, nR, ...
        sprintf('Corrected, OMA-optimal \\Delta\\omega per config   (Q=%d)', Q_target), ...
        sprintf('%s/eye_grid_corr_opt.png', Q_folder));
    plot_eye_grid(eye_aopt, configs, Mbit_arr, df_opt, group_size, n_groups, ymax_cfg, nP, nR, ...
        sprintf('As-fabricated at the optimal setting   (Q=%d, offset %g GHz uncorrected)', Q_target, fab_offset_GHz), ...
        sprintf('%s/eye_grid_asfab_opt.png', Q_folder));

    if save_individual
        save_eye_pngs(eye_c0,   configs, Mbit_arr, z,      group_size, n_groups, Q_folder, 'corr_df0');
        save_eye_pngs(eye_a0,   configs, Mbit_arr, z,      group_size, n_groups, Q_folder, 'asfab_df0');
        save_eye_pngs(eye_copt, configs, Mbit_arr, df_opt, group_size, n_groups, Q_folder, 'corr_opt');
        save_eye_pngs(eye_aopt, configs, Mbit_arr, df_opt, group_size, n_groups, Q_folder, 'asfab_opt');
    end

    %% Table: signed mean OMA and vertical eye opening (console + Q_folder txt)
    table_file = sprintf('%s/eye_table_Q%d.txt', Q_folder, Q_target);
    fid = fopen(table_file, 'w');
    for outid = [1, fid]
        fprintf(outid, '\n  Q = %d   (offset %g GHz)   negative = inverted/closed eye\n', Q_target, fab_offset_GHz);
        fprintf(outid, '                       signed mean OMA                   vertical eye opening\n');
        fprintf(outid, '  Pin   Rate  | asf0   cor0   asfOpt corOpt | asf0   cor0   asfOpt corOpt\n');
        for ci = 1:N_configs
            fprintf(outid, '  %-4.0f  %-5.1f | %5.3f  %5.3f  %5.3f  %5.3f | %6.3f %6.3f %6.3f %6.3f\n', ...
                configs(ci,1), configs(ci,2), ...
                oma_a0(ci), oma_c0(ci), oma_ao(ci), oma_co(ci), ...
                eo_a0(ci),  eo_c0(ci),  eo_ao(ci),  eo_co(ci));
        end
    end
    fclose(fid);
    fprintf('Saved %s\n', table_file);

    %% Save results for this Q
    save(sprintf('%s/eye_results_Q%d.mat', Q_folder, Q_target), ...
        'Q_target', 'tau_c_nom', 'configs', 'powers', 'datarates', 'group_size', 'n_groups', ...
        'df_opt', 'oma_a0', 'oma_c0', 'oma_ao', 'oma_co', ...
        'eo_a0', 'eo_c0', 'eo_ao', 'eo_co', ...
        'fab_offset_GHz', 'res_alpha_scale', 'res_tau_c_scale', 'prbs_order');
end

fprintf('\nDone. Outputs in "%s".\n', out_root);

%% local functions

function [Pt, oma, eopen] = run_eye_case(E_in, V_wave, dw_GHz, total_steps, ...
        M_bit, N_warmup, N_sym, prbs_bits, P_in_W, params, omega_0, consts)
%RUN_EYE_CASE  one eye case for a given ring and laser detuning.
    dtau   = params.gamma_0 * consts.dt_real;
    dw_ext = 2*pi * dw_GHz * 1e9;
    dw_n   = dw_ext / params.gamma_0;
    wp     = omega_0 + dw_ext;
    cFCg   = (consts.Gamma_FCA * consts.c_light^2 * consts.beta_TPA / ...
              (2 * consts.hbar * wp * consts.V_FCA^2 * consts.n_Si^2)) / params.gamma_0;

    Pt = sim_mrr_eye(E_in, V_wave, dw_n, dtau, total_steps, ...
        M_bit, N_warmup, N_sym, ...
        params.mu_n, params.glin_n, params.cTPA_n, params.cFCA_n, cFCg, ...
        params.cTh_n, params.gamma_0, params.iTFC_n, params.iTTH_n, ...
        consts.dn_dN, consts.dn_dT, ...
        params.dw_NL_scale, params.dw_mod_n, params.iTRC_n, params.mu_raw);

    ctr = max(round(M_bit/2), 1);
    Mx  = reshape(Pt, M_bit, N_sym);
    cs  = Mx(ctr, :).';
    hi  = cs(prbs_bits == 1);
    lo  = cs(prbs_bits == 0);
    oma   = (mean(hi) - mean(lo)) / P_in_W;   % signed OMA: negative = inverted polarity
    eopen = (min(hi) - max(lo)) / P_in_W;     % signed vertical eye opening (negative = closed/inverted)
end

function [E, t] = build_eye(W, M_bit, group_size, n_groups)
    spt = group_size * M_bit;
    E   = reshape(W, spt, n_groups);
    t   = (0:spt-1)' / M_bit;
end

function plot_eye_grid(eye_cell, configs, Mbit_arr, df_lab, group_size, n_groups, ymax_cfg, nP, nR, suptxt, fname)
    fig = figure('Position', [40 40 1600 1000], 'Color', 'w');
    movegui(fig, 'center');
    tl = tiledlayout(nP, nR, 'Padding', 'compact', 'TileSpacing', 'compact');
    for ci = 1:size(configs,1)
        P = configs(ci,1);  R = configs(ci,2);  M = Mbit_arr(ci);
        [E, t] = build_eye(eye_cell{ci}, M, group_size, n_groups);
        E = E * 1e3;
        nexttile(ci);
        hl = plot(t, E, 'LineWidth', 0.4);
        set(hl, 'Color', [0.10 0.35 0.75 0.08]);
        xlim([0 group_size]); ylim([0 ymax_cfg(ci)]); grid on; box on;
        title(sprintf('%.0f mW, %.1f Gbps  (\\Delta\\omega=%g GHz)', P, R, df_lab(ci)), 'FontWeight', 'bold');
        if mod(ci-1, nR) == 0, ylabel('P_{thru}  [mW]', 'FontWeight', 'bold'); end
        if ci > (nP-1)*nR,     xlabel('Time  [UI]',     'FontWeight', 'bold'); end
        set(gca, 'FontSize', 11);
    end
    title(tl, suptxt, 'FontWeight', 'bold', 'FontSize', 14);
    print(fig, fname, '-dpng', '-r300');
    fprintf('Saved %s\n', fname);
end

function save_eye_pngs(eye_cell, configs, Mbit_arr, df_lab, group_size, n_groups, folder, suffix)
    for ci = 1:size(configs,1)
        P = configs(ci,1);  R = configs(ci,2);  M = Mbit_arr(ci);
        [E, t] = build_eye(eye_cell{ci}, M, group_size, n_groups);
        E = E * 1e3;
        f = figure('Position', [80 80 700 480], 'Color', 'w', 'Visible', 'off');
        hl = plot(t, E, 'LineWidth', 0.5);
        set(hl, 'Color', [0.10 0.35 0.75 0.08]);
        xlim([0 group_size]); grid on; box on;
        xlabel('Time  [UI]', 'FontSize', 15, 'FontWeight', 'bold');
        ylabel('Through-port power  [mW]', 'FontSize', 15, 'FontWeight', 'bold');
        title(sprintf('%.0f mW, %.1f Gbps, \\Delta\\omega = %g GHz', P, R, df_lab(ci)), ...
            'FontSize', 14, 'FontWeight', 'bold');
        set(gca, 'FontSize', 13);
        fname = sprintf('%s/eye_Pin%.0f_%.0fGbps_%s.png', folder, P*1000, R*10, suffix);
        print(f, fname, '-dpng', '-r300');
        close(f);
    end
    fprintf('Saved individual PNGs (%s)\n', suffix);
end

function bits = gen_prbs_len(N, order)
%GEN_PRBS_LEN  maximal-length PRBS of a given order, N bits (Fibonacci LFSR).
    switch order
        case 7,  taps = [7 6];
        case 9,  taps = [9 5];
        case 10, taps = [10 7];
        case 11, taps = [11 9];
        case 15, taps = [15 14];
        otherwise, error('Unsupported PRBS order %d', order);
    end
    reg  = ones(1, order);
    bits = zeros(N, 1);
    for k = 1:N
        bits(k) = reg(order);
        fb = reg(taps(1));
        for j = 2:numel(taps), fb = xor(fb, reg(taps(j))); end
        reg = [fb, reg(1:order-1)];
    end
end

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