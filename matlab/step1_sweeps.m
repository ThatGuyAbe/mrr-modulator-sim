% step 1, single-ring sweeps. nominal ring, every detuning of the -80 to +80 GHz grid
% (0.25 GHz) is its own run with its own warm-up.
%   dw_opt  detuning of largest oma at 2, 3, 4 mW for every Q at 2.5 and 15 Gb/s
%           (the power-aware detuning, pac)
%   dw_lin  the same at 0.01 mW and 2.5 Gb/s (low-power optimum)
%   the oma and veo curves of the mechanism figure are a subset of these sweeps.
% every later step reads its targets from here. a checkpoint is saved after every Q,
% so a stopped run restarts where it left off.
% output: out/step1_sweeps/results.mat, summary.txt, quicklook.png
% cost: (5 Q x 3 powers x 2 bit rates + 5 Q at 0.01 mW) x 641 points = 22,435 runs,
% about 2 h on 4 workers. the script times one run and prints an estimate.

clear; clc; close all;
[P, S] = mrr_params();
ut = utils();
od = ut.out_folder(S, 'step1_sweeps');
nW = ut.start_pool();

Q_values  = S.Q_values;   powers = S.powers;   datarates = S.datarates;
df_grid   = S.df_grid;
nQ = numel(Q_values);  nP = numel(powers);  nR = numel(datarates);  nD = numel(df_grid);
i0 = find(abs(df_grid) < 1e-9, 1);

% one drive per bit rate (the 0.01 mW reference uses the 2.5 Gb/s drive)
DV = cell(nR,1);
for r = 1:nR, DV{r} = make_drive(datarates(r), S, P); end
dv_lin = make_drive(S.rate_lin, S, P);

% nominal rings
PN = cell(nQ,1);  hw_GHz = zeros(nQ,1);  tau_c = zeros(nQ,1);
for q = 1:nQ
    PN{q} = ring_params(P, Q_values(q));
    hw_GHz(q) = PN{q}.hw_GHz;  tau_c(q) = PN{q}.tau_c;
end

% timing
t0 = tic;  run_case('eye', 2e-3, 0, PN{1}, DV{1}, P, S);  t25 = toc(t0);
t0 = tic;  run_case('eye', 2e-3, 0, PN{1}, DV{end}, P, S);  t15 = toc(t0);
n_runs = nD*nQ*(nP*nR + 1);
fprintf('workers %d. one run %.2f s (%.1f Gb/s), %.2f s (%.1f Gb/s). %d runs, about %.1f h\n', ...
    nW, t25, datarates(1), t15, datarates(end), n_runs, ...
    nD*nQ*(nP*(t25 + t15*(nR-1)) + t25)/nW/3600);

% storage
oma_sw = nan(nD, nQ, nP, nR);  veo_sw = oma_sw;  ber_sw = oma_sw;
lv_sw  = nan(nD, 4, nQ, nP, nR);
oma_lin = nan(nD, nQ);  veo_lin = oma_lin;  ber_lin = oma_lin;
t_q = zeros(nQ,1);

for q = 1:nQ
    ck = fullfile(od, sprintf('ckpt_Q%d.mat', Q_values(q)));
    if exist(ck, 'file') == 2
        L = load(ck);
        if isequal(L.df_grid, df_grid) && isequal(L.powers, powers) && isequal(L.datarates, datarates) ...
                && isequal(L.label, P.label) && L.n_tau_hold == S.n_tau_hold
            oma_sw(:,q,:,:) = L.om;  veo_sw(:,q,:,:) = L.eo;  ber_sw(:,q,:,:) = L.bb;
            lv_sw(:,:,q,:,:) = L.lv;  oma_lin(:,q) = L.oml;  veo_lin(:,q) = L.eol;  ber_lin(:,q) = L.bbl;
            t_q(q) = L.t;
            fprintf('Q %d loaded from checkpoint\n', Q_values(q));
            continue;
        end
    end
    tq = tic;
    pn = PN{q};
    % jobs: (detuning, power, rate) plus the 0.01 mW sweep as power index 0
    [jd, jp, jr] = ndgrid(1:nD, 0:nP, 1:nR);
    keep = ~(jp == 0 & jr > 1);                        % the 0.01 mW sweep only once
    jd = jd(keep);  jp = jp(keep);  jr = jr(keep);
    nJ = numel(jd);
    o = zeros(nJ,1);  e = zeros(nJ,1);  b = zeros(nJ,1);  l = zeros(nJ,4);
    parfor j = 1:nJ
        if jp(j) == 0
            Pw = S.P_lin_mW*1e-3;  dvj = dv_lin;
        else
            Pw = powers(jp(j))*1e-3;  dvj = DV{jr(j)};
        end
        [~, o1, e1, b1, l1] = run_case('eye', Pw, df_grid(jd(j)), pn, dvj, P, S);
        o(j) = o1;  e(j) = e1;  b(j) = b1;  l(j,:) = l1;
    end
    om = nan(nD, nP, nR);  eo = om;  bb = om;  lv = nan(nD, 4, nP, nR);
    for j = 1:nJ
        if jp(j) == 0
            oma_lin(jd(j), q) = o(j);  veo_lin(jd(j), q) = e(j);  ber_lin(jd(j), q) = b(j);
        else
            om(jd(j), jp(j), jr(j)) = o(j);  eo(jd(j), jp(j), jr(j)) = e(j);  bb(jd(j), jp(j), jr(j)) = b(j);
            lv(jd(j), :, jp(j), jr(j)) = l(j,:);
        end
    end
    oma_sw(:,q,:,:) = om;  veo_sw(:,q,:,:) = eo;  ber_sw(:,q,:,:) = bb;  lv_sw(:,:,q,:,:) = lv;
    oml = oma_lin(:,q);  eol = veo_lin(:,q);  bbl = ber_lin(:,q);  t = toc(tq);  t_q(q) = t;
    label = P.label;  n_tau_hold = S.n_tau_hold;
    save(ck, 'om', 'eo', 'bb', 'lv', 'oml', 'eol', 'bbl', 'df_grid', 'powers', 'datarates', 'label', ...
        'n_tau_hold', 't', '-v7');
    fprintf('Q %d done in %.1f min\n', Q_values(q), t/60);
end

%% targets: first index of the maximum
dw_opt = zeros(nQ,nP,nR);  dw_veo = dw_opt;  oma_pk = dw_opt;  veo_at_pk = dw_opt;  veo_pk = dw_opt;
oma_0 = dw_opt;  veo_0 = dw_opt;  edge = dw_opt;
for q = 1:nQ, for ip = 1:nP, for r = 1:nR
    c = oma_sw(:,q,ip,r);  v = veo_sw(:,q,ip,r);
    [oma_pk(q,ip,r), k] = max(c);   dw_opt(q,ip,r) = df_grid(k);  veo_at_pk(q,ip,r) = v(k);
    [veo_pk(q,ip,r), kv] = max(v);  dw_veo(q,ip,r) = df_grid(kv);
    oma_0(q,ip,r) = c(i0);  veo_0(q,ip,r) = v(i0);
    % hot-branch edge: nearest grid step on the red side of the peak where the oma
    % drops by more than 0.05 from one point to the next (NaN = no such step)
    dd = -diff(c(k:-1:1));                  % drop when moving one step to the red
    ke = find(dd > 0.05, 1);
    if isempty(ke), edge(q,ip,r) = NaN; else, edge(q,ip,r) = df_grid(k - ke + 1); end
end, end, end
dw_lin = zeros(nQ,1);  oma_lin_pk = zeros(nQ,1);  dw_lin_an = zeros(nQ,1);  max_dev_an = zeros(nQ,1);
for q = 1:nQ
    [oma_lin_pk(q), k] = max(oma_lin(:,q));  dw_lin(q) = df_grid(k);
    % closed-form linear limit (no heating, no carriers), constant drive levels
    pn = PN{q};  kap = pn.mu_raw^2 / pn.gamma_0;
    x1 = (2*pi*df_grid*1e9 - P.dW0_dVpn*P.V_low ) / pn.gamma_0;
    x0 = (2*pi*df_grid*1e9 - P.dW0_dVpn*P.V_high) / pn.gamma_0;
    an = (((1-kap)^2 + x1.^2) ./ (1 + x1.^2) - ((1-kap)^2 + x0.^2) ./ (1 + x0.^2)).';
    [~, ka] = max(an);  dw_lin_an(q) = df_grid(ka);
    max_dev_an(q) = max(abs(oma_lin(:,q) - an));
end

%% summary
fid = fopen(fullfile(od, 'summary.txt'), 'w');
ut.tee(fid, 'step 1, single-ring sweeps, %s. nominal ring, grid %g to %g GHz step %g, warm-up %d-bit ramp + %d-bit hold at 2.5 Gb/s\n', ...
    P.label, df_grid(1), df_grid(end), df_grid(2)-df_grid(1), DV{1}.N_ramp, DV{1}.N_hold);
ut.tee(fid, '\nPOWER-AWARE DETUNING. dw/2pi in GHz (+ = laser on the blue side of the cold resonance, - = red side)\n');
p_txt = strjoin(arrayfun(@(x) sprintf('%g', x), powers, 'UniformOutput', false), ' / ');
hdr = '    Q    HW    dw_lin |';
for r = 1:nR, hdr = [hdr, sprintf('  dw_opt %4.1f Gb/s: %s mW   |', datarates(r), p_txt)]; end %#ok<AGROW>
ut.tee(fid, '%s\n', hdr);
for q = 1:nQ
    ut.tee(fid, '%6d  %5.2f  %+6.2f |', Q_values(q), hw_GHz(q), dw_lin(q));
    for r = 1:nR, ut.tee(fid, ['      ', repmat('%+6.2f  ', 1, nP), '      |'], dw_opt(q,:,r)); end
    ut.tee(fid, '\n');
end
ut.tee(fid, '\nmove of dw_opt from %g to %g mW (GHz, and half-widths). negative = toward the red\n', powers(1), powers(end));
for q = 1:nQ
    ut.tee(fid, '%6d |', Q_values(q));
    for r = 1:nR
        d = dw_opt(q,end,r) - dw_opt(q,1,r);
        ut.tee(fid, '  %4.1f Gb/s: %+6.2f GHz (%+5.2f HW) |', datarates(r), d, d/hw_GHz(q));
    end
    ut.tee(fid, '\n');
end
ut.tee(fid, '\nnonlinear part dw_opt - dw_lin at 2.5 Gb/s (GHz, HW), %s mW\n', p_txt);
r25 = find(abs(datarates - 2.5) < 1e-9, 1);
for q = 1:nQ
    d = squeeze(dw_opt(q,:,r25)) - dw_lin(q);
    ut.tee(fid, ['%6d |', repmat(' %+6.2f (%+5.2f) ', 1, nP), '\n'], Q_values(q), [d(:).'; d(:).'/hw_GHz(q)]);
end
ut.tee(fid, '\nper curve: OMA peak, VEO at the OMA peak, VEO peak and where, OMA and VEO at dw = 0 (ZDC on the nominal ring)\n');
ut.tee(fid, '    Q   P  rate | dw_opt  OMA pk  VEO there | dw_veo  VEO pk | OMA(0)   VEO(0)  | hot-branch edge (GHz)\n');
for q = 1:nQ, for r = 1:nR, for ip = 1:nP
    ut.tee(fid, '%6d  %g  %4.1f | %+6.2f  %6.4f  %+7.4f  | %+6.2f  %+7.4f | %+7.4f  %+7.4f | %s\n', Q_values(q), powers(ip), ...
        datarates(r), dw_opt(q,ip,r), oma_pk(q,ip,r), veo_at_pk(q,ip,r), dw_veo(q,ip,r), veo_pk(q,ip,r), ...
        oma_0(q,ip,r), veo_0(q,ip,r), num2str(edge(q,ip,r), '%+.2f'));
end, end, end
ut.tee(fid, '\nlow-power reference (%.2f mW, %.1f Gb/s): dw_lin from the simulation and from the closed-form linear limit\n', ...
    S.P_lin_mW, S.rate_lin);
for q = 1:nQ
    ut.tee(fid, '%6d | dw_lin %+6.2f, closed form %+6.2f, OMA max %.4f, max |sim - closed form| %.1e\n', ...
        Q_values(q), dw_lin(q), dw_lin_an(q), oma_lin_pk(q), max_dev_an(q));
end
n_red = sum(sum(dw_opt(:,end,:) < dw_opt(:,1,:)));  n_blue = sum(sum(dw_opt(:,end,:) > dw_opt(:,1,:)));
ut.tee(fid, '\nDIRECTION: from %g to %g mW the target moves to the red in %d and to the blue in %d of %d (Q, bit rate) cases\n', ...
    powers(1), powers(end), n_red, n_blue, nQ*nR);
ut.tee(fid, 'timing: %.1f min in total on %d workers\n', sum(t_q)/60, nW);
fclose(fid);

run_info = struct('label', P.label, 'm_ring', P.m_ring, 'beta_TPA', P.beta_TPA, 'n_tau_ramp', S.n_tau_ramp, ...
    'n_tau_hold', S.n_tau_hold, 'N_sym', S.N_sym, 'P_lin_mW', S.P_lin_mW, 'smoke', S.smoke); %#ok<NASGU>
save(fullfile(od, 'results.mat'), 'Q_values', 'powers', 'datarates', 'df_grid', 'hw_GHz', 'tau_c', ...
    'oma_sw', 'veo_sw', 'ber_sw', 'lv_sw', 'oma_lin', 'veo_lin', 'ber_lin', 'dw_opt', 'dw_veo', 'dw_lin', ...
    'oma_pk', 'veo_at_pk', 'veo_pk', 'oma_0', 'veo_0', 'edge', 'oma_lin_pk', 'dw_lin_an', 'max_dev_an', ...
    't_q', 'run_info', '-v7');

%% quick look (skipped if no graphics are available)
try
    fig = figure('Position', [40 40 1500 420], 'Color', 'w', 'Visible', 'off');
    show = find(ismember(Q_values, [3000 9000 15000]));
    pc = [0.12 0.31 0.61; 0.75 0.22 0.17; 0 0 0];
    for k = 1:numel(show)
        q = show(k);
        subplot(1, numel(show), k); hold on;
        plot(df_grid, oma_lin(:,q), '-', 'Color', [0.5 0.5 0.5]);
        for ip = 1:nP
            plot(df_grid, oma_sw(:,q,ip,r25), '-', 'Color', pc(ip,:), 'LineWidth', 1.3);
            plot(df_grid, veo_sw(:,q,ip,r25), '--', 'Color', pc(ip,:));
            plot(dw_opt(q,ip,r25), oma_pk(q,ip,r25), 'o', 'Color', pc(ip,:));
        end
        plot([0 0], [-1 1], 'k-');  plot(df_grid([1 end]), [0 0], 'k-');
        xlim(dw_lin(q) + [-4 4]*hw_GHz(q));  ylim([-0.6 0.9]);  grid on;  box on;
        title(sprintf('Q = %d, 2.5 Gb/s', Q_values(q)));  xlabel('\Delta\omega/2\pi (GHz)');
        if k == 1, ylabel('OMA (solid), VEO (dashed)'); end
    end
    ut.save_png(fig, fullfile(od, 'quicklook.png'));
catch err
    fprintf('quick-look figure skipped: %s\n', err.message);
end
fprintf('\nsaved to "%s"\n', od);
