function varargout = run_case(kind, varargin)
% one simulated run of one ring, and its metrics.
%
%   [Pt, oma, veo, ber, lv] = run_case('eye', P_in_W, dw_GHz, p, dv, P, S)
%   [T, W, D] = run_case('trace', P_in_W, dw_GHz, p, dv, P, S, dec, snip_bits)
%   [oma, veo, lber] = run_case('eye_dt', P_in_W, dw_GHz, p, P, S, rate, dt, n_tau_ramp, n_tau_hold, N_sym)
%   cFCg = run_case('cfcg', P, p, dw_GHz)
%
% dw_GHz is the laser detuning from the cold, unbiased resonance (+ = laser on the
% blue side). metrics are taken at the center sample of every symbol:
%   oma  signed oma / P_in (negative = inverted levels)
%   veo  signed vertical eye opening / P_in (min one - max zero)
%   ber  0.5 erfc(Q_f/sqrt 2), Q_f from the level means and spreads with the
%        receiver floor sigma_rx = S.sigma_rx_frac * P_in added in quadrature
%   lv   [mean one, mean zero, std one, std zero] / P_in, so the ber can be
%        recomputed for any receiver floor
% 'trace' also records the internal states (W every dec-th step of the whole run,
% D every step of the data window) and reduces them to window and bit statistics.
% 'eye_dt' is the time-step check: the bit period and decision instant stay at
% their 1 ps values and only the integration step changes.
    switch kind
        case 'eye',    [varargout{1:nargout}] = eye_case(varargin{:});
        case 'trace',  [varargout{1:nargout}] = trace_case(varargin{:});
        case 'eye_dt', [varargout{1:nargout}] = eye_case_dt(varargin{:});
        case 'cfcg',   [varargout{1:nargout}] = cfcg_of(varargin{:});
        otherwise, error('unknown run_case kind %s', kind);
    end
end

function [Pt, oma, veo, ber, lv] = eye_case(P_in_W, dw_GHz, p, dv, P, S)
    dtau   = p.gamma_0 * S.dt_real;
    dw_ext = 2*pi * dw_GHz * 1e9;
    dw_n   = dw_ext / p.gamma_0;
    cFCg   = cfcg_of(P, p, dw_GHz);
    Pt = sim_mrr_eye(sqrt(P_in_W), dv.V_wave, dw_n, dtau, dv.total_steps, ...
        dv.M_bit, dv.N_ramp, dv.N_warmup, dv.N_sym, ...
        p.mu_n, p.glin_n, p.cTPA_n, p.cFCA_n, cFCg, ...
        p.cTh_n, p.gamma_0, p.iTFC_n, p.iTTH_n, ...
        P.dn_dN, P.dn_dT, ...
        p.dw_NL_scale, p.dw_mod_n, p.iTRC_n, p.mu_raw);
    ctr = max(round(dv.M_bit/2), 1);
    Mx  = reshape(Pt, dv.M_bit, dv.N_sym);
    cs  = Mx(ctr, :).';
    hi  = cs(dv.prbs_bits == 1);
    lo  = cs(dv.prbs_bits == 0);
    oma = (mean(hi) - mean(lo)) / P_in_W;
    veo = (min(hi) - max(lo)) / P_in_W;
    mu1 = mean(hi); mu0 = mean(lo);
    s1  = std(hi);  s0  = std(lo);
    sig = S.sigma_rx_frac * P_in_W;
    s1n = sqrt(s1^2 + sig^2);
    s0n = sqrt(s0^2 + sig^2);
    qfac = (mu1 - mu0) / (s1n + s0n);
    ber  = 0.5 * erfc(qfac / sqrt(2));
    lv   = [mu1 mu0 s1 s0] / P_in_W;
end

function [oma, veo, lber] = eye_case_dt(P_in_W, dw_GHz, p, P, S, rate_Gbps, dt, n_tau_ramp, n_tau_hold, N_sym)
    dv    = make_drive(rate_Gbps, S, P, n_tau_ramp, n_tau_hold, N_sym, dt);
    M1    = round(dv.T_bit/S.dt_real);              % samples per bit at 1 ps
    ctr1  = max(round(M1/2), 1);                    % decision sample at 1 ps
    ctr   = round(ctr1 * S.dt_real / dt);
    bits  = dv.prbs_bits;
    dw_ext = 2*pi*dw_GHz*1e9;
    cFCg   = cfcg_of(P, p, dw_GHz);
    Pt = sim_mrr_eye(sqrt(P_in_W), dv.V_wave, dw_ext/p.gamma_0, p.gamma_0*dt, dv.total_steps, ...
        dv.M_bit, dv.N_ramp, dv.N_warmup, N_sym, p.mu_n, p.glin_n, p.cTPA_n, p.cFCA_n, cFCg, ...
        p.cTh_n, p.gamma_0, p.iTFC_n, p.iTTH_n, P.dn_dN, P.dn_dT, ...
        p.dw_NL_scale, p.dw_mod_n, p.iTRC_n, p.mu_raw);
    cs = reshape(Pt, dv.M_bit, N_sym);  cs = cs(ctr, :).';
    hi = cs(bits == 1);  lo = cs(bits == 0);
    oma = (mean(hi) - mean(lo)) / P_in_W;
    veo = (min(hi) - max(lo)) / P_in_W;
    sig = S.sigma_rx_frac * P_in_W;
    qf  = (mean(hi) - mean(lo)) / (sqrt(std(hi)^2 + sig^2) + sqrt(std(lo)^2 + sig^2));
    lber = log10(max(0.5 * erfc(qf / sqrt(2)), 1e-30));
end

function [T, W, D] = trace_case(P_in_W, dw_GHz, p, dv, P, S, dec, keep_snip_bits)
% dec             keep every dec-th step of the whole run in W (Inf = none)
% keep_snip_bits  data bits kept at full resolution in T.snip, e.g. 1:20
    if nargin < 7 || isempty(dec), dec = Inf; end
    if nargin < 8, keep_snip_bits = []; end
    dtau   = p.gamma_0 * S.dt_real;
    dw_ext = 2*pi * dw_GHz * 1e9;
    dw_n   = dw_ext / p.gamma_0;
    cFCg   = cfcg_of(P, p, dw_GHz);
    if isinf(dec), dec_run = dv.total_steps + 1; else, dec_run = dec; end
    [W, D] = sim_mrr_fulltrace(sqrt(P_in_W), dv.V_wave, dw_n, dtau, dv.total_steps, ...
        dv.M_bit, dv.N_ramp, dv.N_warmup, dv.N_sym, ...
        p.mu_n, p.glin_n, p.cTPA_n, p.cFCA_n, cFCg, ...
        p.cTh_n, p.gamma_0, p.iTFC_n, p.iTTH_n, ...
        P.dn_dN, P.dn_dT, ...
        p.dw_NL_scale, p.dw_mod_n, p.iTRC_n, p.mu_raw, dec_run, true);
    T = reduce_trace(D, dv, P_in_W, P, p);
    T.dw_GHz = dw_GHz;  T.P_in_W = P_in_W;  T.gamma_0 = p.gamma_0;  T.hw_GHz = p.hw_GHz;
    T.Q_ring = p.Q_ring;  T.M_bit = dv.M_bit;  T.N_warmup = dv.N_warmup;  T.N_ramp = dv.N_ramp;
    if ~isempty(keep_snip_bits)
        zr = (keep_snip_bits(1)-1)*dv.M_bit + 1 : keep_snip_bits(end)*dv.M_bit;
        g  = 1/(2*pi*1e9);
        th = (P.omega_0/P.n_Si) * P.dn_dT * D(zr,4) * g;
        fc = (P.omega_0/P.n_Si) * P.dn_dN * D(zr,3) * g;
        % columns: Pt (W), U (J), thermal (GHz), carrier (GHz), Vpn (V), dwn (half-widths)
        T.snip = [D(zr,1), D(zr,2), th, fc, D(zr,5), D(zr,6)];
    end
end

function S = reduce_trace(D, dv, P_in_W, P, p)
% metrics and internal quantities of one data window.
% D columns: [Pt, U, DN, DT, Vpn, dwn] at every step of the data window.
% thermal and carrier values are the terms they add to the detuning, in GHz
% (+ = resonance moved red by heating, - = moved blue by carriers). dwn is in
% half-widths, + = laser blue of the hot, biased resonance.
    M_bit = dv.M_bit;  N_sym = dv.N_sym;  bits = dv.prbs_bits;
    Pt = D(:,1);  U = D(:,2);  DN = D(:,3);  DT = D(:,4);  dwn = D(:,6);
    ctr = max(round(M_bit/2), 1);
    idx = (0:N_sym-1).' * M_bit + ctr;
    cs  = Pt(idx);
    one = (bits == 1);  zer = ~one;
    S.oma = (mean(cs(one)) - mean(cs(zer))) / P_in_W;
    S.veo = (min(cs(one)) - max(cs(zer))) / P_in_W;
    S.lvl_one  = mean(cs(one)) / P_in_W;
    S.lvl_zero = mean(cs(zer)) / P_in_W;
    g   = 1/(2*pi*1e9);
    thG = (P.omega_0/P.n_Si) * P.dn_dT * DT * g;
    fcG = (P.omega_0/P.n_Si) * P.dn_dN * DN * g;
    S.DT_mean     = mean(DT);
    S.DN_mean     = mean(DN);
    S.thermal_GHz = mean(thG);
    S.carrier_GHz = mean(fcG);
    S.U_mean      = mean(U);
    S.dwn_mean    = mean(dwn);
    S.thermal_min = min(thG);  S.thermal_max = max(thG);
    S.carrier_min = min(fcG);  S.carrier_max = max(fcG);
    fs = sort(fcG);  n = numel(fs);
    S.carrier_p05 = fs(max(1, min(n, round(0.05*(n-1)) + 1)));
    S.carrier_p95 = fs(max(1, min(n, round(0.95*(n-1)) + 1)));
    fc_bit = mean(reshape(fcG, M_bit, N_sym), 1).';
    th_bit = mean(reshape(thG, M_bit, N_sym), 1).';
    DT_bit = mean(reshape(DT,  M_bit, N_sym), 1).';
    S.carrier_bitmean_min = min(fc_bit);  S.carrier_bitmean_max = max(fc_bit);
    S.thermal_bitmean_min = min(th_bit);  S.thermal_bitmean_max = max(th_bit);
    % temperature drift over the data window: last quarter minus first quarter
    q4 = max(1, floor(N_sym/4));
    S.DT_drift_K  = mean(DT_bit(end-q4+1:end)) - mean(DT_bit(1:q4));
    S.DT_drift_HW = S.DT_drift_K * (P.omega_0/P.n_Si) * P.dn_dT / p.gamma_0;
    % at the decision instants
    cd = dwn(idx);
    S.dwn_one = mean(cd(one));     S.dwn_zero = mean(cd(zer));
    S.dwn_one_sd = std(cd(one));   S.dwn_zero_sd = std(cd(zer));
    S.red_one = mean(cd(one) < 0); S.red_zero = mean(cd(zer) < 0);
    S.U_one  = mean(U(idx(one)));  S.U_zero  = mean(U(idx(zer)));
    S.fc_one = mean(fcG(idx(one)));  S.fc_zero = mean(fcG(idx(zer)));
    S.th_one = mean(thG(idx(one)));  S.th_zero = mean(thG(idx(zer)));
    % does the resonance sweep through the laser inside a bit
    Mx = reshape(dwn, M_bit, N_sym);
    dmin = min(Mx, [], 1).';  dmax = max(Mx, [], 1).';
    Px = reshape(Pt, M_bit, N_sym);
    pmin = min(Px, [], 1).';
    prev = [dv.prev_bit; bits(1:end-1)];
    t01 = (prev == 0) & one;  t10 = (prev == 1) & zer;  non = (prev == bits);
    cross = (dmin < 0) & (dmax > 0);
    S.cross_01 = mean(cross(t01));  S.cross_10 = mean(cross(t10));  S.cross_non = mean(cross(non));
    S.pmin_01 = mean(pmin(t01)) / P_in_W;  S.pmin_10 = mean(pmin(t10)) / P_in_W;
    S.c_Pt = cs;  S.c_dwn = cd;  S.c_U = U(idx);  S.c_fc = fcG(idx);  S.c_th = thG(idx);
    S.DT_bit = DT_bit;
end

function cFCg = cfcg_of(P, p, dw_GHz)
% normalized carrier generation coefficient at laser detuning dw_GHz (element-wise)
    dw_ext = 2*pi * dw_GHz * 1e9;
    wp     = P.omega_0 + dw_ext;
    cFCg   = (P.Gamma_FCA * P.c_light^2 * P.beta_TPA ./ ...
              (2 * P.hbar * wp * P.V_FCA^2 * P.n_Si^2)) / p.gamma_0;
end
