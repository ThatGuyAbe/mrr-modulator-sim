% step 6, population. run after step 1.
% all five Q, 2, 3, 4 mW, 2.5 and 15 Gb/s (30 configurations).
%   zdc  laser at the cold resonance of each ring (dw = 0), all 2000 rings
%   pac  laser at dw_opt from step 1, all 2000 rings
%   af   as fabricated, laser at the cold design resonance, rings 1 to 500
% the mean and spread of the one and zero levels of every ring are saved (LV_*), so
% the ber can be recomputed for any receiver floor without a new run.
% a checkpoint is saved after every (Q, configuration), so a stopped run restarts
% where it left off. set Q_run below to run some Q values only.
% output: out/step6_population/pop_results_Q<value>.mat and summary.txt
% cost: 30 configurations x (2 x 2000 + 500) = 135,000 runs, about 10 h on 4 workers.

clear; clc; close all;
[P, S] = mrr_params();
ut = utils();
od = ut.out_folder(S, 'step6_population');
nW = ut.start_pool();
T1 = ut.load_targets(S);
R  = ring_population(S, P);

Q_run     = S.Q_values;            % edit to run a subset, e.g. [9000 15000]
powers    = S.powers;
datarates = S.datarates;
N_rings   = S.N_rings;
n_asf     = S.n_asf;
asf_rings = (1:n_asf).';
if ~isequal(T1.datarates(:).', datarates(:).') || ~isequal(T1.powers(:).', powers(:).')
    error('step 1 targets were computed for other powers or bit rates. rerun step1_sweeps.m.');
end

% configurations: power outer, bit rate inner
nP = numel(powers);  nR = numel(datarates);
configs = zeros(nP*nR, 2);  ic = 0;
for ip = 1:nP, for ir = 1:nR, ic = ic + 1; configs(ic,:) = [powers(ip), datarates(ir)]; end, end
nCfg = size(configs,1);
DV = cell(nR,1);
for ir = 1:nR, DV{ir} = make_drive(datarates(ir), S, P); end

% timing
p1 = ring_params(P, Q_run(1), R.alpha_scale(1), R.tau_c_scale(1));
tt = zeros(nR,1);
for ir = 1:nR, t0 = tic; run_case('eye', S.P_ref*1e-3, 0, p1, DV{ir}, P, S); tt(ir) = toc(t0); end
n_per_cfg = 2*N_rings + n_asf;
fprintf('workers %d, %d rings (AF on %d). one run: %s s. %d configs x %d Q x %d runs, about %.1f h\n', ...
    nW, N_rings, n_asf, mat2str(round(100*tt.')/100), nCfg, numel(Q_run), n_per_cfg, ...
    n_per_cfg*nP*sum(tt)*numel(Q_run)/nW/3600);

for Q_target = Q_run
    iq = find(T1.Q_values == Q_target, 1);
    pn = ring_params(P, Q_target);
    tau_c_nom = pn.tau_c;
    pr_cell = cell(N_rings,1);
    for r = 1:N_rings
        pr_cell{r} = ring_params(P, Q_target, R.alpha_scale(r), R.tau_c_scale(r));
    end
    OMA_asf = nan(N_rings, nCfg);  OMA_cor = OMA_asf;  OMA_opt = OMA_asf;
    EO_asf  = OMA_asf;  EO_cor  = OMA_asf;  EO_opt  = OMA_asf;
    BER_asf = OMA_asf;  BER_cor = OMA_asf;  BER_opt = OMA_asf;
    LV_asf  = nan(N_rings, 4, nCfg);  LV_cor = LV_asf;  LV_opt = LV_asf;
    df_opt  = zeros(nCfg,1);  t_cfg = zeros(nCfg,1);
    oma_vs_df = zeros(numel(T1.df_grid), nCfg);  veo_vs_df = oma_vs_df;
    for ci = 1:nCfg
        P_in_mW = configs(ci,1);  rate = configs(ci,2);
        ip = find(powers == P_in_mW, 1);  ir = find(abs(datarates - rate) < 1e-9, 1);
        dco = T1.dw_opt(iq, ip, ir);
        df_opt(ci) = dco;
        oma_vs_df(:,ci) = T1.oma_sw(:, iq, ip, ir);  veo_vs_df(:,ci) = T1.veo_sw(:, iq, ip, ir);
        ck = fullfile(od, sprintf('ckpt_Q%d_cfg%d.mat', Q_target, ci));
        if exist(ck, 'file') == 2
            L = load(ck);
            if L.dco == dco && isequal(L.cfg, configs(ci,:)) && L.N_rings == N_rings && L.n_asf == n_asf ...
                    && strcmp(L.label, P.label) && L.n_tau_hold == S.n_tau_hold
                OMA_asf(:,ci) = L.o(:,1);  OMA_cor(:,ci) = L.o(:,2);  OMA_opt(:,ci) = L.o(:,3);
                EO_asf(:,ci)  = L.e(:,1);  EO_cor(:,ci)  = L.e(:,2);  EO_opt(:,ci)  = L.e(:,3);
                BER_asf(:,ci) = L.b(:,1);  BER_cor(:,ci) = L.b(:,2);  BER_opt(:,ci) = L.b(:,3);
                LV_asf(:,:,ci) = L.l(:,:,1);  LV_cor(:,:,ci) = L.l(:,:,2);  LV_opt(:,:,ci) = L.l(:,:,3);
                t_cfg(ci) = L.t;
                fprintf('Q %d config %d loaded from checkpoint\n', Q_target, ci);
                continue;
            end
        end
        dv = DV{ir};  Pw = P_in_mW*1e-3;  off = R.df_fab_GHz;
        o = nan(N_rings,3);  e = o;  b = o;  l = nan(N_rings,4,3);
        tc = tic;
        parfor r = 1:N_rings
            pr = pr_cell{r};
            [~, o2, e2, b2, l2] = run_case('eye', Pw, 0,   pr, dv, P, S);     % zdc
            [~, o3, e3, b3, l3] = run_case('eye', Pw, dco, pr, dv, P, S);     % pac
            if r <= n_asf
                [~, o1, e1, b1, l1] = run_case('eye', Pw, off(r), pr, dv, P, S);   % as fabricated
            else
                o1 = NaN;  e1 = NaN;  b1 = NaN;  l1 = nan(1,4);
            end
            o(r,:) = [o1 o2 o3];  e(r,:) = [e1 e2 e3];  b(r,:) = [b1 b2 b3];
            l(r,:,:) = reshape([l1; l2; l3].', [1 4 3]);
        end
        t = toc(tc);
        cfg = configs(ci,:);  label = P.label;  n_tau_hold = S.n_tau_hold;
        save(ck, 'o', 'e', 'b', 'l', 'dco', 'cfg', 'N_rings', 'n_asf', 'label', 'n_tau_hold', 't', '-v7');
        OMA_asf(:,ci) = o(:,1);  OMA_cor(:,ci) = o(:,2);  OMA_opt(:,ci) = o(:,3);
        EO_asf(:,ci)  = e(:,1);  EO_cor(:,ci)  = e(:,2);  EO_opt(:,ci)  = e(:,3);
        BER_asf(:,ci) = b(:,1);  BER_cor(:,ci) = b(:,2);  BER_opt(:,ci) = b(:,3);
        LV_asf(:,:,ci) = l(:,:,1);  LV_cor(:,:,ci) = l(:,:,2);  LV_opt(:,:,ci) = l(:,:,3);
        t_cfg(ci) = t;
        fprintf('Q %d, %g mW, %4.1f Gb/s, dw_opt %+.2f GHz: %.1f min | invalid ZDC %d, PAC %d, AF %d of %d\n', ...
            Q_target, P_in_mW, rate, dco, t/60, sum(e(:,2) < 0), sum(e(:,3) < 0), sum(e(asf_rings,1) < 0), n_asf);
    end
    df_grid = T1.df_grid;  df_fab_GHz = R.df_fab_GHz;  alpha_scale = R.alpha_scale;  tau_c_scale = R.tau_c_scale;
    sigma_lam = S.sigma_lam;  sigma_Q_rel = S.sigma_Q_rel;  sigma_rx_frac = S.sigma_rx_frac;
    n_tau_ramp = S.n_tau_ramp;  n_tau_hold = S.n_tau_hold;  prbs_order = S.prbs_order;
    param_label = P.label;  m_ring = P.m_ring;  beta_TPA = P.beta_TPA;  smoke = S.smoke; %#ok<NASGU>
    save(fullfile(od, sprintf('pop_results_Q%d.mat', Q_target)), ...
        'Q_target', 'tau_c_nom', 'N_rings', 'n_asf', 'asf_rings', 'configs', 'powers', 'datarates', ...
        'df_opt', 'df_grid', 'oma_vs_df', 'veo_vs_df', ...
        'OMA_asf', 'OMA_cor', 'OMA_opt', 'EO_asf', 'EO_cor', 'EO_opt', 'BER_asf', 'BER_cor', 'BER_opt', ...
        'LV_asf', 'LV_cor', 'LV_opt', 'df_fab_GHz', 'alpha_scale', 'tau_c_scale', ...
        'sigma_lam', 'sigma_Q_rel', 'sigma_rx_frac', 'n_tau_ramp', 'n_tau_hold', 'prbs_order', ...
        'param_label', 'm_ring', 'beta_TPA', 'smoke', 't_cfg', '-v7');
    fprintf('Q %d saved (%.1f h)\n', Q_target, sum(t_cfg)/3600);
end

%% summary over every Q that has a results file
fid = fopen(fullfile(od, 'summary.txt'), 'w');
ut.tee(fid, 'step 6, population, %s. %d rings (ZDC, PAC), AF on rings 1-%d. invalid = VEO < 0 (an inverted ring, OMA < 0, always has VEO < 0)\n', ...
    P.label, N_rings, n_asf);
ut.tee(fid, '\n    Q   P  rate  dw_opt |   mean OMA            |   mean VEO            |  mean log10 BER       |  invalid (neg. OMA)                 | PAC better: OMA BER inv\n');
ut.tee(fid, '                        |  AF     ZDC    PAC    |  AF     ZDC    PAC    |  AF     ZDC     PAC   |  AF/%d     ZDC/%d     PAC/%d    |\n', n_asf, N_rings, N_rings);
n_cfg_all = 0;  n_win = 0;  lb = @(x) mean(log10(max(x, 1e-30)));
for Q_target = S.Q_values
    f = fullfile(od, sprintf('pop_results_Q%d.mat', Q_target));
    if exist(f, 'file') ~= 2, continue; end
    D = load(f);
    for ci = 1:size(D.configs,1)
        a = D.asf_rings;
        mo = [mean(D.OMA_asf(a,ci)), mean(D.OMA_cor(:,ci)), mean(D.OMA_opt(:,ci))];
        me = [mean(D.EO_asf(a,ci)),  mean(D.EO_cor(:,ci)),  mean(D.EO_opt(:,ci))];
        mb = [lb(D.BER_asf(a,ci)),   lb(D.BER_cor(:,ci)),   lb(D.BER_opt(:,ci))];
        ni = [sum(D.EO_asf(a,ci) < 0), sum(D.EO_cor(:,ci) < 0), sum(D.EO_opt(:,ci) < 0)];
        nn = [sum(D.OMA_asf(a,ci) < 0), sum(D.OMA_cor(:,ci) < 0), sum(D.OMA_opt(:,ci) < 0)];
        w = [mo(3) > mo(2), mb(3) < mb(2), ni(3) <= ni(2)];
        n_cfg_all = n_cfg_all + 1;  n_win = n_win + all(w);
        ut.tee(fid, '%6d  %g  %4.1f  %+6.2f | %+6.3f %+6.3f %+6.3f | %+6.3f %+6.3f %+6.3f | %6.2f %7.2f %7.2f | %3d (%3d)  %4d (%4d)  %4d (%4d) | %s %s %s\n', ...
            Q_target, D.configs(ci,1), D.configs(ci,2), D.df_opt(ci), mo, me, mb, ni(1), nn(1), ni(2), nn(2), ...
            ni(3), nn(3), ut.yesno(w(1)), ut.yesno(w(2)), ut.yesno(w(3)));
    end
end
ut.tee(fid, '\nPAC has a higher mean OMA, a lower mean BER and no more invalid rings than ZDC in %d of %d configurations\n', ...
    n_win, n_cfg_all);
fclose(fid);
fprintf('\nsaved to "%s"\n', od);
