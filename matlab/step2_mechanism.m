% step 2, what moves the best detuning. run after step 1.
% nominal ring at Q = 3000, 9000, 15000 and 2.5 Gb/s, at 0.01 mW and 2, 3, 4 mW, each at
% dw = 0 (zdc) and at the power-aware detuning (dw_lin at 0.01 mW), with the internal
% states recorded. per run: window means of the thermal and carrier resonance shifts
% (GHz, + = red), their bit-to-bit range, the stored energy, where the laser sits from
% the hot resonance at the decision instants (half-widths, + = blue side), the
% temperature rise and its drift, and the first 20 data bits at full resolution.
% each run must give the same oma and veo as the step 1 sweep at that point.
% output: out/step2_mechanism/results.mat, summary.txt. 24 runs, a few minutes.

clear; clc; close all;
[P, S] = mrr_params();
ut = utils();
od = ut.out_folder(S, 'step2_mechanism');
nW = ut.start_pool();
T1 = ut.load_targets(S);
L1 = load(fullfile(S.out_root, 'step1_sweeps', 'results.mat'));

rate    = 2.5;
Q_list  = [3000 9000 15000];
P_trace = [S.P_lin_mW S.powers];
n_snip  = 20;
dv = make_drive(rate, S, P);
ir = find(abs(T1.datarates - rate) < 1e-9, 1);

% jobs: (Q, power, case) with case 1 = dw 0, case 2 = target
jobs = [];
for q = Q_list
    for pm = P_trace
        if pm < 0.1
            dwt = T1.dw_lin(T1.Q_values == q);
        else
            dwt = T1.get(q, pm, rate);
        end
        jobs = [jobs; q pm 0 1; q pm dwt 2]; %#ok<AGROW>
    end
end
nJ = size(jobs,1);
TR = cell(nJ,1);
t0 = tic;
parfor j = 1:nJ
    p = ring_params(P, jobs(j,1));
    TR{j} = run_case('trace', jobs(j,2)*1e-3, jobs(j,3), p, dv, P, S, Inf, 1:n_snip);
end
t_run = toc(t0);
TRs = [TR{:}];

% check against the step 1 sweep (same run, so the difference must be zero)
chk = nan(nJ,1);
for j = 1:nJ
    iq = find(L1.Q_values == jobs(j,1), 1);
    k  = find(abs(L1.df_grid - jobs(j,3)) < 1e-9, 1);
    if jobs(j,2) < 0.1
        if ~isempty(k), chk(j) = max(abs([TRs(j).oma - L1.oma_lin(k,iq), TRs(j).veo - L1.veo_lin(k,iq)])); end
    else
        ip = find(abs(L1.powers - jobs(j,2)) < 1e-9, 1);
        if ~isempty(k), chk(j) = max(abs([TRs(j).oma - L1.oma_sw(k,iq,ip,ir), TRs(j).veo - L1.veo_sw(k,iq,ip,ir)])); end
    end
end

fid = fopen(fullfile(od, 'summary.txt'), 'w');
ut.tee(fid, 'step 2, mechanism traces, %s. nominal ring, %.1f Gb/s. shifts are the terms in delta_omega (GHz): + = red (heating), - = blue (carriers)\n', P.label, rate);
ut.tee(fid, 'dwn = laser position from the hot, biased resonance at the decision instants, in half-widths, + = laser on the blue side\n\n');
ut.tee(fid, '    Q   P(mW)  case    dw   |  DT (K)  thermal  carrier   sum  (HW)  | carrier bit-to-bit  | U mean  U one  U zero (fJ) | dwn one  dwn zero | OMA     VEO    | DT drift (HW) | check\n');
for j = 1:nJ
    T = TRs(j);  hw = T.hw_GHz;
    if jobs(j,4) == 1, nm = 'ZDC'; elseif jobs(j,2) < 0.1, nm = 'lin'; else, nm = 'PAC'; end
    ut.tee(fid, '%6d  %5.2f  %-4s %+6.2f | %6.3f  %+7.3f  %+7.3f  %+6.3f (%+5.2f) | %+6.3f to %+6.3f  | %6.2f  %6.2f  %6.2f    | %+7.3f  %+7.3f  | %+6.4f %+6.4f | %+8.4f      | %.0e\n', ...
        jobs(j,1), jobs(j,2), nm, jobs(j,3), T.DT_mean, T.thermal_GHz, T.carrier_GHz, ...
        T.thermal_GHz + T.carrier_GHz, (T.thermal_GHz + T.carrier_GHz)/hw, ...
        T.carrier_bitmean_min, T.carrier_bitmean_max, T.U_mean*1e15, T.U_one*1e15, T.U_zero*1e15, ...
        T.dwn_one, T.dwn_zero, T.oma, T.veo, T.DT_drift_HW, chk(j));
end
ut.tee(fid, '\nlargest difference from the step 1 sweep: %.1e (must be 0)\n', max(chk));

% target move and resonance shifts per Q
ut.tee(fid, '\ntarget move and shifts:\n');
for q = Q_list
    pq = ring_params(P, q);  hw = pq.hw_GHz;
    Pl = S.P_low;  Ph = S.P_ref;                  % 2 and 4 mW
    d2 = T1.get(q, Pl, rate);  d4 = T1.get(q, Ph, rate);
    jz = find(jobs(:,1) == q & jobs(:,2) == Ph & jobs(:,4) == 1);
    jo = find(jobs(:,1) == q & jobs(:,2) == Ph & jobs(:,4) == 2);
    ut.tee(fid, ' Q %5d: dw_opt %+.2f GHz at %g mW, %+.2f GHz at %g mW, move %+.2f GHz = %+.2f HW. at %g mW: thermal %+.2f GHz (%+.2f HW), carrier %+.2f GHz (%+.2f HW) at dw_opt\n', ...
        q, d2, Pl, d4, Ph, d4 - d2, (d4 - d2)/hw, Ph, TRs(jo).thermal_GHz, TRs(jo).thermal_GHz/hw, TRs(jo).carrier_GHz, TRs(jo).carrier_GHz/hw);
    ut.tee(fid, '          stored energy at %g mW: ZDC %.2f fJ (one %.2f, zero %.2f), PAC %.2f fJ (one %.2f, zero %.2f). zero detuning sits %.2f HW from dw_opt\n', ...
        Ph, TRs(jz).U_mean*1e15, TRs(jz).U_one*1e15, TRs(jz).U_zero*1e15, TRs(jo).U_mean*1e15, TRs(jo).U_one*1e15, ...
        TRs(jo).U_zero*1e15, abs(d4)/hw);
end
ut.tee(fid, '\n%d runs in %.1f min on %d workers\n', nJ, t_run/60, nW);
fclose(fid);
save(fullfile(od, 'results.mat'), 'jobs', 'TRs', 'chk', 'rate', 'Q_list', 'P_trace', 'n_snip', ...
    'P', '-v7');
fprintf('\nsaved to "%s"\n', od);
