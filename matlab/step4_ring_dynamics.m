% step 4, ring dynamics: why zdc fails. run after steps 1 and 3.
% part 1  full traces (warm-up and data) with the internal states
%           case1   nominal ring, Q 15000, 2 mW, 2.5 Gb/s, zdc
%           case2   nominal ring, Q 15000, 4 mW, 2.5 Gb/s, zdc
%           case3   nominal ring, Q 15000, 4 mW, 2.5 Gb/s, pac
%           case5a-5d  the step 3 ring (Q 9000 design, loss and coupling time x 1.10)
%                   at 4 mW, zdc and pac, 2.5 and 15 Gb/s
%         each case must give the same oma and veo as step 1 (cases 1-3) or step 3.
% part 2  self-pulsation test. nominal ring, 4 mW, cw, constant voltage -0.5, -2.5,
%         -4.5 V, at dw = 0 and at dw_opt, Q 15000 and 9000. the same 250 ns ramp,
%         then 1000 ns held. swing over the last 500 ns.
% part 3  bistability. nominal ring, 4 mW, cw, constant voltage -0.5, -2.5, -4.5 V,
%         Q 15000 and 9000. the detuning steps through the grid inside one run, blue
%         to red and back, 250 ns per point.
% output: out/step4_ring_dynamics/summary.txt, traces.mat, selfpulse.mat,
% hysteresis.mat and quick-look png files. about 20 min on 4 workers.

clear; clc; close all;
[P, S] = mrr_params();
ut = utils();
od = ut.out_folder(S, 'step4_ring_dynamics');
nW = ut.start_pool();
T1 = ut.load_targets(S);
L1 = load(fullfile(S.out_root, 'step1_sweeps', 'results.mat'));
f3 = fullfile(S.out_root, 'step3_eyes', 'results.mat');
have3 = exist(f3, 'file') == 2;
if have3, L3 = load(f3); end

run_part   = [true true true];
dec_trace  = 20;                       % keep every 20th step of the whole run
zoom_bits  = 20:39;                    % data bits kept at full resolution
dw_cw      = S.cw_lo:0.25:40;          % GHz, hysteresis grid, -50 to 40
Pl = S.P_low;  Ph = S.P_ref;  Pw_ref = Ph*1e-3;   % 2 and 4 mW
P_cw_mW = Ph;                                       % saved with the cw results
hold_ns    = 250;  rec_ns = 50;  settle0_ns = 250;
bis_tol    = 0.01;                 % P_in, difference between the two sweep directions that counts as bistable
V_cw       = [P.V_high P.V_bias P.V_low];   % -0.5, -2.5, -4.5 V
sp_hold_ns = 1000;
if S.smoke, dw_cw = -30:1:20; hold_ns = 50; sp_hold_ns = 200; end

dv25 = make_drive(2.5, S, P);
dv15 = make_drive(15,  S, P);
p15  = ring_params(P, 15000);
p9   = ring_params(P, 9000);
pfig = ring_params(P, 9000, 1.10, 1.10);
d15_25 = T1.get(15000, Ph, 2.5);  d15_15 = T1.get(15000, Ph, 15);
d9_25  = T1.get(9000, Ph, 2.5);   d9_15  = T1.get(9000, Ph, 15);
summ = {};
fid = fopen(fullfile(od, 'summary.txt'), 'w');
ut.tee(fid, 'step 4, ring dynamics, %s. half-widths: Q 15000 %.3f GHz, Q 9000 %.3f GHz, step 3 ring %.3f GHz (loaded Q %.0f)\n', ...
    P.label, p15.hw_GHz, p9.hw_GHz, pfig.hw_GHz, pfig.Q_ring);
ut.tee(fid, 'dw_opt (step 1): Q 15000, %g mW: %+.2f GHz (2.5 Gb/s), %+.2f GHz (15 Gb/s). Q 9000, %g mW: %+.2f, %+.2f GHz\n', Ph, ...
    d15_25, d15_15, Ph, d9_25, d9_15);

%% part 1, traces
if run_part(1)
    i15 = find(L1.Q_values == 15000, 1);  ir25 = find(abs(L1.datarates - 2.5) < 1e-9, 1);
    k0  = find(abs(L1.df_grid) < 1e-9, 1);  kp = find(abs(L1.df_grid - d15_25) < 1e-9, 1);
    ip2 = find(L1.powers == Pl, 1);  ip4 = find(L1.powers == Ph, 1);
    C = struct('name', {}, 'p', {}, 'P_mW', {}, 'rate', {}, 'dw', {}, 'exp_oma', {}, 'exp_veo', {}, 'note', {});
    C(end+1) = struct('name','case1','p',p15,'P_mW',Pl,'rate',2.5,'dw',0, 'exp_oma', L1.oma_sw(k0,i15,ip2,ir25), ...
        'exp_veo', L1.veo_sw(k0,i15,ip2,ir25), 'note','nominal Q 15000, ZDC');
    C(end+1) = struct('name','case2','p',p15,'P_mW',Ph,'rate',2.5,'dw',0, 'exp_oma', L1.oma_sw(k0,i15,ip4,ir25), ...
        'exp_veo', L1.veo_sw(k0,i15,ip4,ir25), 'note','nominal Q 15000, ZDC');
    C(end+1) = struct('name','case3','p',p15,'P_mW',Ph,'rate',2.5,'dw',d15_25, 'exp_oma', L1.oma_sw(kp,i15,ip4,ir25), ...
        'exp_veo', L1.veo_sw(kp,i15,ip4,ir25), 'note','nominal Q 15000, PAC');
    fig_cases = {'case5a', Ph, 2.5, 0;  'case5b', Ph, 15, 0;  'case5c', Ph, 2.5, d9_25;  'case5d', Ph, 15, d9_15};
    for k = 1:4
        eo = NaN;  ev = NaN;
        if have3
            m = find(L3.case_P_mW == fig_cases{k,2} & abs(L3.case_rate - fig_cases{k,3}) < 1e-9 & ...
                     abs(L3.case_dw - fig_cases{k,4}) < 1e-9, 1);
            if ~isempty(m), eo = L3.oma(m); ev = L3.veo(m); end
        end
        if fig_cases{k,4} == 0, nm = 'step 3 ring, ZDC'; else, nm = 'step 3 ring, PAC'; end
        C(end+1) = struct('name', fig_cases{k,1}, 'p', pfig, 'P_mW', fig_cases{k,2}, 'rate', fig_cases{k,3}, ...
            'dw', fig_cases{k,4}, 'exp_oma', eo, 'exp_veo', ev, 'note', nm); %#ok<SAGROW>
    end
    nC = numel(C);
    TRB = cell(nC,1);
    t1 = tic;
    parfor c = 1:nC
        if C(c).rate == 2.5, dv = dv25; else, dv = dv15; end
        Pw = C(c).P_mW * 1e-3;
        [T, W, D] = run_case('trace', Pw, C(c).dw, C(c).p, dv, P, S, dec_trace);
        zr = (zoom_bits(1)-1)*dv.M_bit + 1 : zoom_bits(end)*dv.M_bit;
        T.W = single(W);  T.Z = single(D(zr, :));
        T.name = C(c).name;  T.note = C(c).note;  T.P_mW = C(c).P_mW;  T.rate = C(c).rate;  T.dw = C(c).dw;
        T.exp_oma = C(c).exp_oma;  T.exp_veo = C(c).exp_veo;
        T.check = max(abs([T.oma - C(c).exp_oma, T.veo - C(c).exp_veo]));
        TRB{c} = T;
    end
    t_tr = toc(t1);
    TRB = [TRB{:}];
    dt_real = S.dt_real;  N_sym = S.N_sym;  V_high = P.V_high;  V_low = P.V_low;  prbs_bits = dv25.prbs_bits;
    save(fullfile(od, 'traces.mat'), 'TRB', 'prbs_bits', 'zoom_bits', 'dec_trace', ...
        'dt_real', 'N_sym', 'V_high', 'V_low', 't_tr', '-v7');
    ut.tee(fid, '\nPART 1. traces (%.1f min). check = difference from step 1 (cases 1-3) or step 3 (case 5), must be 0 (NaN = not available)\n', t_tr/60);
    ut.tee(fid, ' case    ring Q  P(mW) rate  dw(GHz) |  OMA      VEO    check  | dwn one  dwn zero (HW) | red one  red zero | DT (K)  thermal  carrier (GHz) | U one  U zero  mean (fJ)\n');
    for c = 1:nC
        T = TRB(c);
        ut.tee(fid, ' %-6s %6.0f  %4.0f  %4.1f  %+6.2f  | %+7.4f  %+7.4f  %.0e | %+7.3f  %+7.3f      |  %5.3f   %5.3f   | %6.3f  %+7.3f  %+7.3f      | %5.1f  %5.1f  %5.1f\n', ...
            T.name, T.Q_ring, T.P_mW, T.rate, T.dw, T.oma, T.veo, T.check, T.dwn_one, T.dwn_zero, ...
            T.red_one, T.red_zero, T.DT_mean, T.thermal_GHz, T.carrier_GHz, T.U_one*1e15, T.U_zero*1e15, T.U_mean*1e15);
    end
    ut.tee(fid, ' dwn = laser position from the hot, biased resonance at the decision instants (+ = blue side).\n');
    ut.tee(fid, ' red one / red zero = share of one / zero decision instants with the laser on the red side (dwn < 0)\n');
    ut.tee(fid, '\n does the resonance sweep through the laser inside a bit, and how deep is the dip?\n');
    ut.tee(fid, ' case   | share of bits with a sign change: 0->1  1->0  none | lowest Pt in the bit / P_in: 0->1  1->0 | level zero, one (P_in) | DT drift (HW)\n');
    for c = 1:nC
        T = TRB(c);
        ut.tee(fid, ' %-6s |                                   %5.3f  %5.3f  %5.3f |                          %5.3f  %5.3f |  %5.3f  %5.3f         | %+8.4f\n', ...
            T.name, T.cross_01, T.cross_10, T.cross_non, T.pmin_01, T.pmin_10, T.lvl_zero, T.lvl_one, T.DT_drift_HW);
    end
end

%% part 2, self-pulsation
if run_part(2)
    SPQ = [15000*ones(1,6), 9000*ones(1,6)];
    SPV = repmat(V_cw, 1, 4);
    SPD = [0 0 0, d15_25*[1 1 1], 0 0 0, d9_25*[1 1 1]];
    nS = numel(SPQ);
    SP = cell(nS,1);
    Nd_sp = round(sp_hold_ns*1e-9/(dv25.M_bit*S.dt_real));
    t2 = tic;
    parfor s = 1:nS
        if SPQ(s) == 15000, pp = p15; else, pp = p9; end
        dv = struct('M_bit', dv25.M_bit, 'N_ramp', dv25.N_ramp, 'N_warmup', dv25.N_ramp, 'N_sym', Nd_sp, ...
                    'total_steps', (dv25.N_ramp + Nd_sp)*dv25.M_bit, 'prev_bit', 0, 'prbs_bits', zeros(Nd_sp,1));
        dv.V_wave = SPV(s) * ones(dv.total_steps, 1);
        cf = run_case('cfcg', P, pp, SPD(s));
        [Ws, Ds] = sim_mrr_fulltrace(sqrt(Pw_ref), dv.V_wave, 2*pi*SPD(s)*1e9/pp.gamma_0, pp.gamma_0*S.dt_real, ...
            dv.total_steps, dv.M_bit, dv.N_ramp, dv.N_warmup, dv.N_sym, pp.mu_n, pp.glin_n, pp.cTPA_n, ...
            pp.cFCA_n, cf, pp.cTh_n, pp.gamma_0, pp.iTFC_n, pp.iTTH_n, P.dn_dN, P.dn_dT, ...
            pp.dw_NL_scale, pp.dw_mod_n, pp.iTRC_n, pp.mu_raw, 1000, true);
        n_last = min(round(0.5*sp_hold_ns*1e-9/S.dt_real), size(Ds,1));
        last = Ds(end-n_last+1:end, :);
        Rs = struct();
        Rs.Q = SPQ(s);  Rs.V = SPV(s);  Rs.dw = SPD(s);
        Rs.Pt_mean = mean(last(:,1));  Rs.Pt_pp = max(last(:,1)) - min(last(:,1));
        Rs.DN_mean = mean(last(:,3));  Rs.DN_pp = max(last(:,3)) - min(last(:,3));
        Rs.DT_mean = mean(last(:,4));  Rs.DT_pp = max(last(:,4)) - min(last(:,4));
        Rs.dwn_mean = mean(last(:,6)); Rs.dwn_pp = max(last(:,6)) - min(last(:,6));
        Rs.W = single(Ws);
        SP{s} = Rs;
    end
    t_sp = toc(t2);
    SP = [SP{:}];
    save(fullfile(od, 'selfpulse.mat'), 'SP', 'sp_hold_ns', 't_sp', 'P_cw_mW', '-v7');
    ut.tee(fid, '\nPART 2. self-pulsation, cw, constant voltage, %g mW, %d ns held after the 250 ns ramp (%.1f min)\n', Ph, sp_hold_ns, t_sp/60);
    ut.tee(fid, ' over the second half of the hold. pp = max - min. a steady ring gives pp near zero.\n');
    ut.tee(fid, '    Q    V(V)  dw(GHz) | Pt mean/Pin  Pt pp/Pin | DT mean (K)  DT pp (mK) | DN pp/DN mean | dwn mean  dwn pp (HW)\n');
    for s = 1:nS
        Rs = SP(s);
        ut.tee(fid, ' %6d  %5.1f  %+6.2f | %8.4f    %8.1e  |  %7.4f     %8.1e  |   %8.1e    | %+7.3f   %8.1e\n', ...
            Rs.Q, Rs.V, Rs.dw, Rs.Pt_mean/Pw_ref, Rs.Pt_pp/Pw_ref, Rs.DT_mean, Rs.DT_pp*1e3, ...
            Rs.DN_pp/max(Rs.DN_mean, eps), Rs.dwn_mean, Rs.dwn_pp);
    end
end

%% part 3, hysteresis
if run_part(3)
    HQ = [15000 9000];  nHQ = numel(HQ);  nV = numel(V_cw);
    [jq, jv, jd] = ndgrid(1:nHQ, 1:nV, 1:2);
    jq = jq(:);  jv = jv(:);  jd = jd(:);
    nJ = numel(jq);
    RH = cell(nJ,1);
    hold_steps = round(hold_ns*1e-9/S.dt_real);  rec_steps = round(rec_ns*1e-9/S.dt_real);
    settle0_steps = round(settle0_ns*1e-9/S.dt_real);
    t3 = tic;
    parfor j = 1:nJ
        if HQ(jq(j)) == 15000, pp = p15; else, pp = p9; end
        if jd(j) == 1, dl = fliplr(dw_cw); else, dl = dw_cw; end      % 1 = blue to red, 2 = red to blue
        cf = run_case('cfcg', P, pp, dl);
        Rw = sim_mrr_cw_sweep(sqrt(Pw_ref), V_cw(jv(j)), 2*pi*dl*1e9/pp.gamma_0, cf, pp.gamma_0*S.dt_real, ...
            hold_steps, rec_steps, settle0_steps, ...
            pp.mu_n, pp.glin_n, pp.cTPA_n, pp.cFCA_n, pp.cTh_n, pp.gamma_0, pp.iTFC_n, pp.iTTH_n, ...
            P.dn_dN, P.dn_dT, pp.dw_NL_scale, pp.dw_mod_n, pp.iTRC_n, pp.mu_raw);
        if jd(j) == 1, Rw = flipud(Rw); end
        RH{j} = Rw;
    end
    t_hy = toc(t3);
    HYS = zeros(numel(dw_cw), 10, nHQ, nV, 2);
    for j = 1:nJ, HYS(:,:,jq(j),jv(j),jd(j)) = RH{j}; end
    hw_hys = [p15.hw_GHz p9.hw_GHz];
    dW0_dVpn = P.dW0_dVpn;
    save(fullfile(od, 'hysteresis.mat'), 'HYS', 'HQ', 'V_cw', 'dw_cw', 'hold_ns', 'rec_ns', 'P_cw_mW', ...
        'settle0_ns', 'hw_hys', 'dW0_dVpn', 'd15_25', 'd15_15', 'd9_25', 'd9_15', 't_hy', '-v7');
    ut.tee(fid, '\nPART 3. hysteresis, cw, %g mW, %g ns per point, grid %g to %g GHz step %g (%.1f min)\n', ...
        Ph, hold_ns, dw_cw(1), dw_cw(end), dw_cw(2)-dw_cw(1), t_hy/60);
    ut.tee(fid, ' bistable = points where the blue-to-red and red-to-blue sweeps differ by more than %g P_in\n', bis_tol);
    [~, i0] = min(abs(dw_cw));
    for iq = 1:nHQ
        if HQ(iq) == 15000, dop = [d15_25 d15_15]; else, dop = [d9_25 d9_15]; end
        for iv = 1:nV
            Pb = HYS(:,1,iq,iv,1)/Pw_ref;  Pr = HYS(:,1,iq,iv,2)/Pw_ref;
            bis = abs(Pb - Pr) > bis_tol;
            if any(bis)
                rng_txt = sprintf('%+g to %+g GHz (%d points), largest difference %.3f P_in', min(dw_cw(bis)), ...
                    max(dw_cw(bis)), sum(bis), max(abs(Pb - Pr)));
            else
                rng_txt = sprintf('none (largest difference %.1e P_in)', max(abs(Pb - Pr)));
            end
            [~, ip1] = min(abs(dw_cw - dop(1)));  [~, ip2] = min(abs(dw_cw - dop(2)));
            ut.tee(fid, ' Q %5d, V %4.1f: bistable %s. peak DT %.3f K\n', HQ(iq), V_cw(iv), rng_txt, max(HYS(:,4,iq,iv,1)));
            ut.tee(fid, '     Pt/Pin [blue-to-red, red-to-blue]: dw 0 [%.4f %.4f], dw_opt %+.2f [%.4f %.4f], dw_opt %+.2f [%.4f %.4f]\n', ...
                Pb(i0), Pr(i0), dop(1), Pb(ip1), Pr(ip1), dop(2), Pb(ip2), Pr(ip2));
        end
    end
end
fclose(fid);

%% quick-look figures (skipped if no graphics are available)
try
    if run_part(1)
        for c = 1:numel(TRB)
            T = TRB(c);  W = double(T.W);  Z = double(T.Z);
            tW = W(:,1)*S.dt_real*1e9;  tZ = (0:size(Z,1)-1).'/T.M_bit;
            lab = {'P_{thr} (mW)', 'U (fJ)', '\DeltaN (10^{23} m^{-3})', '\DeltaT (K)', 'V_{pn} (V)', '\delta\omega/\gamma_0'};
            sc  = [1e3 1e15 1e-23 1 1 1];
            fig = figure('Position', [40 40 1300 1000], 'Color', 'w', 'Visible', 'off');
            for r = 1:6
                subplot(6,2,2*r-1); plot(tW, W(:, r+1)*sc(r), 'k-'); grid on; ylabel(lab{r});
                if r == 6, xlabel('time (ns)'); end
                subplot(6,2,2*r);   plot(tZ, Z(:, r)*sc(r), 'k-'); grid on; xlim([0 numel(zoom_bits)]);
                if r == 6, xlabel(sprintf('time (UI), data bits %d to %d', zoom_bits(1), zoom_bits(end))); end
            end
            ut.save_png(fig, fullfile(od, sprintf('trace_%s.png', T.name)));
        end
    end
    if run_part(3)
        fig = figure('Position', [40 40 1500 700], 'Color', 'w', 'Visible', 'off');
        for iq = 1:2
            for iv = 1:3
                subplot(2,3,(iq-1)*3+iv); hold on;
                plot(dw_cw, HYS(:,1,iq,iv,1)/Pw_ref, 'b-', 'LineWidth', 1.4);
                plot(dw_cw, HYS(:,1,iq,iv,2)/Pw_ref, 'r--', 'LineWidth', 1.2);
                ylim([0 1.05]); grid on; box on;
                title(sprintf('Q %d, V = %.1f V (blue: blue to red, red: red to blue)', HQ(iq), V_cw(iv)));
                xlabel('\Delta\omega/2\pi (GHz)');  ylabel('P_{thr}/P_{in}');
            end
        end
        ut.save_png(fig, fullfile(od, 'hysteresis.png'));
    end
catch err
    fprintf('quick-look figures skipped: %s\n', err.message);
end
fprintf('\nsaved to "%s"\n', od);
