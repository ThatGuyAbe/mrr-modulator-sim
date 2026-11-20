% step 0, model checks. run first, takes a few minutes.
%   a. parameter values and derived numbers (m c_p, R_th, mass from the geometry,
%      thermo-optic shift per kelvin)
%   b. units and heat source. cw input at a constant voltage, held for 1 us, then the
%      settled DT and DN are compared with the steady state of the carrier and
%      temperature equations (time derivatives set to zero)
%   c. heating and resonance shifts of the nominal ring at Q = 9000, 12000, 15000,
%      4 mW, zero detuning, 2.5 Gb/s, and cw sweeps in both directions (bistability)
%   d. warm-up: drift of DT over the data window and the oma and veo for a 250 ns
%      ramp alone, the 250 ns ramp + 250 ns hold used everywhere, and a long
%      reference (250 ns ramp + 750 ns hold)
% output: out/step0_checks/summary.txt and results.mat

clear; clc; close all;
[P, S] = mrr_params();
ut = utils();
od = ut.out_folder(S, 'step0_checks');
fid = fopen(fullfile(od, 'summary.txt'), 'w');
nW = ut.start_pool();
t_all = tic;

%% a. parameters and derived numbers
mcp  = P.m_ring * P.c_p;                                    % J/K
Rth  = P.tau_th / mcp;                                      % K/W
V_ring = 2*pi*P.ring_radius * P.wg_width * P.wg_height;     % m^3
m_geo  = P.rho_Si * V_ring;
dfdT   = P.omega_0 * P.dn_dT / P.n_Si / (2*pi*1e9);          % GHz per K
ut.tee(fid, 'step 0, model checks. %s\n', datestr(now));
ut.tee(fid, '\nA. PARAMETERS (mrr_params.m)\n');
ut.tee(fid, '  m_ring    %.3g kg\n', P.m_ring);
ut.tee(fid, '  beta_TPA  %.3g m/W\n', P.beta_TPA);
ut.tee(fid, '  c_p       %g J/(kg K), dn/dN %.3g m^3, tau_th %g ns, sigma_FCA %.3g m^2\n', ...
    P.c_p, P.dn_dN, P.tau_th*1e9, P.sigma_FCA);
ut.tee(fid, '  ring geometry: radius %.2f um, strip %.0f nm x %.0f nm -> volume %.3f um^3, mass %.3g kg\n', ...
    P.ring_radius*1e6, P.wg_width*1e9, P.wg_height*1e9, V_ring*1e18, m_geo);
ut.tee(fid, '  m c_p = %.3g J/K\n', mcp);
ut.tee(fid, '  R_th = tau_th/(m c_p) = %.2f K/mW, with Gamma_th %.4f: %.2f K/mW\n', ...
    Rth*1e-3, P.Gamma_th, P.Gamma_th*Rth*1e-3);
ut.tee(fid, '  thermo-optic shift of the resonance: %.3f GHz per K (red)\n', dfdT);

%% b. units and heat source
ut.tee(fid, '\nB. UNITS AND HEAT SOURCE\n');
okA = (P.c_p == 700) && (P.m_ring < 1e-12);
ut.tee(fid, '  (a) c_p = %g J/(kg K) with m = %.3g kg: %s\n', P.c_p, P.m_ring, ut.yesno(okA));
ut.tee(fid, '  (b) heat source in mod_rk4_step.m: Pn = (glin_n + gTPA + gFCA)*|a|^2, dT/dtau = -iTTH_n*DT + cTh_n*g0*Pn,\n');
ut.tee(fid, '      cTh_n = Gamma_th/(m c_p gamma_0), one factor |a|^2.\n');
ut.tee(fid, '      test: cw input, constant voltage -2.5 V, approach from the blue, hold 1 us, then compare\n');
ut.tee(fid, '      the settled state with the steady state (d/dt = 0).\n');
ss_Q  = [9000 9000 15000 15000];
ss_dw = [0 -5 0 -8];                 % GHz, zero detuning and a red detuning on the hot branch
ss_P  = S.P_ref*1e-3;               % 4 mW
SS = zeros(numel(ss_Q), 8);
for k = 1:numel(ss_Q)
    p  = ring_params(P, ss_Q(k));
    hw = p.hw_GHz;
    dl = [linspace(ss_dw(k) + 15*hw, ss_dw(k), 60), ss_dw(k)*ones(1,100)];    % 60 ramp points, then 1 us held
    cf = run_case('cfcg', P, p, dl);
    Rw = sim_mrr_cw_sweep(sqrt(ss_P), P.V_bias, 2*pi*dl*1e9/p.gamma_0, cf, p.gamma_0*S.dt_real, ...
        10000, 1000, 0, p.mu_n, p.glin_n, p.cTPA_n, p.cFCA_n, p.cTh_n, p.gamma_0, ...
        p.iTFC_n, p.iTTH_n, P.dn_dN, P.dn_dT, p.dw_NL_scale, p.dw_mod_n, p.iTRC_n, p.mu_raw);
    U  = Rw(end,2);  DN = Rw(end,3);  DT = Rw(end,4);
    wp = P.omega_0 + 2*pi*ss_dw(k)*1e9;
    gFC_real  = P.Gamma_FCA * P.c_light^2 * P.beta_TPA / (2 * P.hbar * wp * P.V_FCA^2 * P.n_Si^2);
    DN_ss = P.tau_FC * gFC_real * U^2;                                   % carrier equation, d/dt = 0
    g_abs = p.gamma_lin + P.beta_TPA*P.c_light^2/(P.n_Si^2*P.V_TPA)*U + ...
            P.Gamma_FCA*P.sigma_FCA*P.c_light/(2*P.n_Si)*DN;             % c alpha/n + gTPA + gFCA
    DT_ss = P.tau_th * P.Gamma_th / (P.m_ring*P.c_p) * g_abs * U;        % temperature equation, d/dt = 0
    DT_extra = DT_ss * U;                                                % what an extra |a|^2 would give
    SS(k,:) = [ss_Q(k), ss_dw(k), U*1e15, DT, DT_ss, abs(DT - DT_ss)/DT_ss, abs(DN - DN_ss)/DN_ss, Rw(end,6)];
    ut.tee(fid, '      Q %5d, dw %+5.1f GHz: U %.2f fJ, laser at %+.3f HW | DT sim %.6f K, steady state %.6f K, rel. diff %.1e | DN rel. diff %.1e | with an extra |a|^2 DT would be %.2e K\n', ...
        ss_Q(k), ss_dw(k), U*1e15, Rw(end,6), DT, DT_ss, SS(k,6), SS(k,7), DT_extra);
end
ut.tee(fid, '  steady-state check: largest relative difference DT %.1e, DN %.1e\n', max(SS(:,6)), max(SS(:,7)));

%% c. heating and resonance shifts
ut.tee(fid, '\nC. HEATING. nominal ring, %g mW, 2.5 Gb/s, dw = 0, window means over the data window\n', S.P_ref);
ut.tee(fid, '       Q   |  DT (K)   thermal (GHz)  (HW)  | carrier (GHz)  (HW)  |  OMA     VEO\n');
dv = make_drive(2.5, S, P);
san_Q = [9000 12000 15000];
SAN = zeros(numel(san_Q), 6);
for k = 1:numel(san_Q)
    p = ring_params(P, san_Q(k));
    T = run_case('trace', S.P_ref*1e-3, 0, p, dv, P, S);
    SAN(k,:) = [T.DT_mean, T.thermal_GHz, T.carrier_GHz, T.oma, T.veo, p.hw_GHz];
    ut.tee(fid, '   %6d  | %8.4f   %+8.3f    %+6.3f  |  %+8.3f    %+6.3f  | %+6.3f  %+6.3f\n', san_Q(k), ...
        T.DT_mean, T.thermal_GHz, T.thermal_GHz/p.hw_GHz, T.carrier_GHz, T.carrier_GHz/p.hw_GHz, T.oma, T.veo);
end

% cw sweeps in both directions at -2.5 V
hys_Q = [9000 12000 15000];
dgrid = -45:0.5:30;
HY = cell(numel(hys_Q), 2);
jobs = [kron((1:numel(hys_Q)).', [1;1]), repmat([1;2], numel(hys_Q), 1)];
HYc = cell(size(jobs,1),1);
parfor j = 1:size(jobs,1)
    p = ring_params(P, hys_Q(jobs(j,1)));
    if jobs(j,2) == 1, dl = fliplr(dgrid); else, dl = dgrid; end      % 1 = blue to red
    cf = run_case('cfcg', P, p, dl);
    Rw = sim_mrr_cw_sweep(sqrt(S.P_ref*1e-3), P.V_bias, 2*pi*dl*1e9/p.gamma_0, cf, p.gamma_0*S.dt_real, ...
        100000, 10000, 250000, p.mu_n, p.glin_n, p.cTPA_n, p.cFCA_n, p.cTh_n, p.gamma_0, ...
        p.iTFC_n, p.iTTH_n, P.dn_dN, P.dn_dT, p.dw_NL_scale, p.dw_mod_n, p.iTRC_n, p.mu_raw);
    if jobs(j,2) == 1, Rw = flipud(Rw); end
    HYc{j} = Rw;
end
for j = 1:size(jobs,1), HY{jobs(j,1), jobs(j,2)} = HYc{j}; end
ut.tee(fid, '\n   cw sweeps at %g mW and -2.5 V, %g to %g GHz, 100 ns per point. bistable where the two directions differ by > 0.01 P_in\n', S.P_ref, dgrid(1), dgrid(end));
for k = 1:numel(hys_Q)
    a = HY{k,1}(:,1)/(S.P_ref*1e-3);  b = HY{k,2}(:,1)/(S.P_ref*1e-3);
    bis = abs(a - b) > 0.01;
    if any(bis)
        ut.tee(fid, '   Q %5d: bistable from %+.1f to %+.1f GHz, largest difference %.3f P_in, peak DT %.3f K\n', ...
            hys_Q(k), min(dgrid(bis)), max(dgrid(bis)), max(abs(a-b)), max(HY{k,1}(:,4)));
    else
        ut.tee(fid, '   Q %5d: no bistable region on this grid (largest difference %.1e P_in), peak DT %.3f K\n', ...
            hys_Q(k), max(abs(a-b)), max(HY{k,1}(:,4)));
    end
end

%% d. warm-up
ut.tee(fid, '\nD. WARM-UP. nominal ring, %g mW. DT drift = last quarter minus first quarter of the data window\n', S.P_ref);
wu_Q    = [9000 15000];
wu_rate = [2.5 15];
wu_dw   = [0 -5];
wu_var  = [5 0; 5 5; 5 15];          % [ramp, hold] in tau_th
wu_name = {'250 ns ramp, no hold', 'used: 250 ns ramp + 250 ns hold', 'reference: 250 ns ramp + 750 ns hold'};
[a1, a2, a3, a4] = ndgrid(1:numel(wu_Q), 1:numel(wu_rate), 1:numel(wu_dw), 1:size(wu_var,1));
WJ = [a1(:) a2(:) a3(:) a4(:)];
WU = zeros(size(WJ,1), 5);
parfor j = 1:size(WJ,1)
    q = wu_Q(WJ(j,1));  r = wu_rate(WJ(j,2));  d = wu_dw(WJ(j,3));
    dd = make_drive(r, S, P, wu_var(WJ(j,4),1), wu_var(WJ(j,4),2));
    p  = ring_params(P, q);
    T  = run_case('trace', S.P_ref*1e-3, d, p, dd, P, S);
    WU(j,:) = [T.oma, T.veo, T.DT_drift_K, T.DT_drift_HW, T.DT_mean];
end
ut.tee(fid, '    Q    rate  dw   | warm-up                                 |  OMA      VEO    | DT mean (K)  DT drift (K)  (HW)\n');
for j = 1:size(WJ,1)
    ut.tee(fid, ' %6d  %4.1f  %+3.0f  | %-39s | %+7.4f  %+7.4f |  %7.4f      %+8.4f   %+7.4f\n', ...
        wu_Q(WJ(j,1)), wu_rate(WJ(j,2)), wu_dw(WJ(j,3)), wu_name{WJ(j,4)}, WU(j,1), WU(j,2), WU(j,5), WU(j,3), WU(j,4));
end
d_old = [];  d_new = [];
for i1 = 1:numel(wu_Q), for i2 = 1:numel(wu_rate), for i3 = 1:numel(wu_dw)
    s = @(v) find(WJ(:,1)==i1 & WJ(:,2)==i2 & WJ(:,3)==i3 & WJ(:,4)==v);
    d_old(end+1) = max(abs(WU(s(1),1:2) - WU(s(3),1:2))); %#ok<SAGROW>
    d_new(end+1) = max(abs(WU(s(2),1:2) - WU(s(3),1:2))); %#ok<SAGROW>
end, end, end
ut.tee(fid, '  largest change of OMA or VEO against the long reference: ramp only %.4f, ramp + hold %.4f\n', max(d_old), max(d_new));

ut.tee(fid, '\nworkers %d, total time %.1f min\n', nW, toc(t_all)/60);
fclose(fid);
save(fullfile(od, 'results.mat'), 'SS', 'SAN', 'san_Q', 'HY', 'hys_Q', 'dgrid', ...
    'WJ', 'WU', 'wu_Q', 'wu_rate', 'wu_dw', 'wu_var', 'Rth', 'mcp', 'm_geo', 'dfdT', '-v7');
fprintf('\nsaved to "%s"\n', od);
