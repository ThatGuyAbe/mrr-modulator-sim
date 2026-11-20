% step 5, numerical checks. run after step 1.
% part a  convergence on the nominal ring
%           stiff  Q 15000, 4 mW, 15 Gb/s, zdc and pac
%           easy   Q 3000,  2 mW, 2.5 Gb/s, zdc and pac
%         variants: time step 1 (baseline), 0.5 and 0.25 ps with the bit period and
%         decision instant held at their 1 ps values, warm-up doubled (500 ns ramp +
%         500 ns hold), data window of 2000 symbols. the baseline must equal the
%         step 1 sweep exactly.
% part b  linear limit. 1 uW cw at the -2.5 V bias, every detuning of the grid its own
%         run (ramp + hold, then 50 bits averaged), all five Q, compared with the
%         closed-form lorentzian. the half-width is read off the simulated dip.
% output: out/step5_convergence/summary.txt and results.mat. about 15 min on 4 workers.

clear; clc; close all;
[P, S] = mrr_params();
ut = utils();
od = ut.out_folder(S, 'step5_convergence');
nW = ut.start_pool();
T1 = ut.load_targets(S);
L1 = load(fullfile(S.out_root, 'step1_sweeps', 'results.mat'));
fid = fopen(fullfile(od, 'summary.txt'), 'w');
ut.tee(fid, 'step 5, convergence, %s\n', P.label);

%% part a
Ph = S.P_ref;  Pl = S.P_low;                     % 4 and 2 mW
CASES = {15000, Ph, 15,  0,                       'stiff ZDC';
         15000, Ph, 15,  T1.get(15000, Ph, 15),   'stiff PAC';
          3000, Pl, 2.5, 0,                       'easy ZDC';
          3000, Pl, 2.5, T1.get(3000, Pl, 2.5),   'easy PAC'};
% variants: time step (ps), ramp (tau_th), hold (tau_th), symbols
VAR = [1    S.n_tau_ramp    S.n_tau_hold    S.N_sym;
       0.5  S.n_tau_ramp    S.n_tau_hold    S.N_sym;
       0.25 S.n_tau_ramp    S.n_tau_hold    S.N_sym;
       1    2*S.n_tau_ramp  2*S.n_tau_hold  S.N_sym;
       1    S.n_tau_ramp    S.n_tau_hold    2*S.N_sym];
var_name = {'baseline 1 ps', 'dt 0.5 ps', 'dt 0.25 ps', 'warm-up doubled', '2x symbols'};
nC = size(CASES,1);  nV = size(VAR,1);
[jc, jv] = ndgrid(1:nC, 1:nV);  jc = jc(:);  jv = jv(:);  nJ = numel(jc);
res = zeros(nJ, 3);
t0 = tic;
parfor j = 1:nJ
    q = CASES{jc(j),1};  Pm = CASES{jc(j),2};  rate = CASES{jc(j),3};  dw = CASES{jc(j),4};
    p = ring_params(P, q);
    [o1, e1, b1] = run_case('eye_dt', Pm*1e-3, dw, p, P, S, rate, VAR(jv(j),1)*1e-12, VAR(jv(j),2), ...
        VAR(jv(j),3), VAR(jv(j),4));
    res(j,:) = [o1, e1, b1];
end
t_a = toc(t0);
RA = reshape(res, nC, nV, 3);
base_chk = nan(nC,1);
for c = 1:nC
    iq = find(L1.Q_values == CASES{c,1});  ip = find(L1.powers == CASES{c,2});
    ir = find(abs(L1.datarates - CASES{c,3}) < 1e-9);  k = find(abs(L1.df_grid - CASES{c,4}) < 1e-9);
    if ~isempty(iq) && ~isempty(ip) && ~isempty(ir) && ~isempty(k)
        base_chk(c) = max(abs([RA(c,1,1) - L1.oma_sw(k,iq,ip,ir), RA(c,1,2) - L1.veo_sw(k,iq,ip,ir)]));
    end
end
ut.tee(fid, '\nPART A. convergence (%d runs, %.1f min). baseline vs the step 1 sweep: largest difference %.1e (must be 0)\n', ...
    nJ, t_a/60, max(base_chk));
worst = zeros(nV, 2);
for c = 1:nC
    ut.tee(fid, '\n %s (Q %d, %g mW, %g Gb/s, dw %+.2f GHz)\n', CASES{c,5}, CASES{c,1}, CASES{c,2}, CASES{c,3}, CASES{c,4});
    ut.tee(fid, '   variant          |  OMA       VEO      log10 BER |  dOMA (%% of |base|)   dVEO (%% of |base|)   dlogBER\n');
    for v = 1:nV
        dO = RA(c,v,1) - RA(c,1,1);  dE = RA(c,v,2) - RA(c,1,2);  dB = RA(c,v,3) - RA(c,1,3);
        ut.tee(fid, '   %-16s | %8.5f  %8.5f  %8.3f  | %+9.2e (%+6.3f)   %+9.2e (%+6.3f)   %+7.3f\n', var_name{v}, ...
            RA(c,v,1), RA(c,v,2), RA(c,v,3), dO, 100*dO/abs(RA(c,1,1)), dE, 100*dE/abs(RA(c,1,2)), dB);
        worst(v,:) = max(worst(v,:), [abs(100*dO/RA(c,1,1)), abs(100*dE/RA(c,1,2))]);
    end
end
ut.tee(fid, '\n largest relative change of OMA or VEO over the %d cases (percent):\n', nC);
for v = 2:nV
    ut.tee(fid, '   %-16s OMA %.3f, VEO %.3f\n', var_name{v}, worst(v,1), worst(v,2));
end

%% part b, linear limit at 1 uW
P_cw_W = 1e-6;  V_lin = P.V_bias;  N_win = 50;
T_bit = 400e-12;  M_bit = round(T_bit/S.dt_real);
N_ramp = ceil(S.n_tau_ramp*P.tau_th/T_bit);  N_hold = ceil(S.n_tau_hold*P.tau_th/T_bit);
dvc = struct('M_bit', M_bit, 'N_ramp', N_ramp, 'N_warmup', N_ramp + N_hold, 'N_sym', N_win, ...
    'total_steps', (N_ramp + N_hold + N_win)*M_bit);
dvc.V_wave = V_lin * ones(dvc.total_steps, 1);
Qv = S.Q_values;  nQ = numel(Qv);  dfg = S.df_grid;  nD = numel(dfg);
PNc = cell(nQ,1);
for q = 1:nQ, PNc{q} = ring_params(P, Qv(q)); end
[jd, jq] = ndgrid(1:nD, 1:nQ);  jd = jd(:);  jq = jq(:);
Tf = zeros(numel(jd),1);
t0 = tic;
parfor j = 1:numel(jd)
    p  = PNc{jq(j)};
    cf = run_case('cfcg', P, p, dfg(jd(j)));
    Pt = sim_mrr_eye(sqrt(P_cw_W), dvc.V_wave, 2*pi*dfg(jd(j))*1e9/p.gamma_0, p.gamma_0*S.dt_real, ...
        dvc.total_steps, dvc.M_bit, dvc.N_ramp, dvc.N_warmup, dvc.N_sym, p.mu_n, p.glin_n, p.cTPA_n, ...
        p.cFCA_n, cf, p.cTh_n, p.gamma_0, p.iTFC_n, p.iTTH_n, P.dn_dN, P.dn_dT, ...
        p.dw_NL_scale, p.dw_mod_n, p.iTRC_n, p.mu_raw);
    Tf(j) = mean(Pt) / P_cw_W;
end
t_b = toc(t0);
T_sim = reshape(Tf, nD, nQ);
T_an = zeros(nD, nQ);  dmax = zeros(nQ,1);  hw_sim = zeros(nQ,1);  hw_ex = zeros(nQ,1);
for q = 1:nQ
    p = PNc{q};  kap = p.mu_raw^2/p.gamma_0;
    x = (2*pi*dfg*1e9 - P.dW0_dVpn*V_lin) / p.gamma_0;
    T_an(:,q) = (((1-kap)^2 + x.^2) ./ (1 + x.^2)).';
    dmax(q) = max(abs(T_sim(:,q) - T_an(:,q)));
    [tmin, im] = min(T_sim(:,q));  lev = (1 + tmin)/2;
    xl = ut.cross_at(dfg(1:im), T_sim(1:im,q), lev);  xr = ut.cross_at(dfg(im:end), T_sim(im:end,q), lev);
    hw_sim(q) = (xr - xl)/2;  hw_ex(q) = p.hw_GHz;
end
ut.tee(fid, '\nPART B. linear limit, %.0f uW cw at %.1f V (%d runs, %.1f min)\n', P_cw_W*1e6, V_lin, numel(jd), t_b/60);
ut.tee(fid, '    Q  | max |T_sim - T_closed form| (P_in) | HW from the dip (GHz)  gamma_0/2pi (GHz)  difference (MHz)\n');
for q = 1:nQ
    ut.tee(fid, '%6d |        %8.1e                   |      %8.4f           %8.4f          %+7.2f\n', ...
        Qv(q), dmax(q), hw_sim(q), hw_ex(q), 1e3*(hw_sim(q) - hw_ex(q)));
end
ut.tee(fid, ' largest: transmission %.1e of P_in, half-width %.2f MHz\n', max(dmax), 1e3*max(abs(hw_sim - hw_ex)));
fclose(fid);
save(fullfile(od, 'results.mat'), 'CASES', 'VAR', 'var_name', 'RA', 'base_chk', 'worst', ...
    'T_sim', 'T_an', 'dmax', 'hw_sim', 'hw_ex', 'Qv', 'dfg', 't_a', 't_b', '-v7');
fprintf('\nsaved to "%s"\n', od);
