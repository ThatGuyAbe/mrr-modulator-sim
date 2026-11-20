% step 3, eye waveforms. run after step 1.
% ring: Q = 9000 design with loss and coupling time both 10 percent (one standard
% deviation) above nominal, a typical fabricated ring rather than the nominal ring on
% which the target is found. zdc (dw = 0) and pac (step 1 target at the same power and
% bit rate) at 2 mW / 2.5 Gb/s, 4 mW / 2.5 Gb/s and 4 mW / 15 Gb/s.
% output: out/step3_eyes/results.mat, summary.txt. six runs, about a minute.

clear; clc; close all;
[P, S] = mrr_params();
ut = utils();
od = ut.out_folder(S, 'step3_eyes');
nW = ut.start_pool();
T1 = ut.load_targets(S);

Q_target        = 9000;
res_alpha_scale = 1.10;
res_tau_c_scale = 1.10;
case_name = {'zdc', 'pac', 'zdc', 'pac', 'zdc', 'pac'};
case_P_mW = [S.P_low S.P_low S.P_ref S.P_ref S.P_ref S.P_ref];   % 2 and 4 mW
case_rate = [2.5 2.5 2.5 2.5 15 15];
case_dw   = zeros(1,6);
for c = 1:6
    if strcmp(case_name{c}, 'pac'), case_dw(c) = T1.get(Q_target, case_P_mW(c), case_rate(c)); end
end

pr = ring_params(P, Q_target, res_alpha_scale, res_tau_c_scale);
Q_ring = pr.Q_ring;
fprintf('ring: loaded Q = %.0f, half-width %.2f GHz\n', Q_ring, pr.hw_GHz);

nC = numel(case_name);
Pt_all = cell(nC,1);  Mbit_all = zeros(nC,1);
oma = zeros(nC,1);  veo = zeros(nC,1);  ber = zeros(nC,1);  lvl = zeros(nC,4);
parfor c = 1:nC
    dv = make_drive(case_rate(c), S, P);
    Pw = case_P_mW(c)*1e-3;
    [Pt, o1, e1, b1] = run_case('eye', Pw, case_dw(c), pr, dv, P, S);
    ctr = max(round(dv.M_bit/2), 1);
    cs  = Pt((0:dv.N_sym-1)*dv.M_bit + ctr) / Pw;
    hi  = cs(dv.prbs_bits == 1);  lo = cs(dv.prbs_bits == 0);
    Pt_all{c} = Pt;  Mbit_all(c) = dv.M_bit;
    oma(c) = o1;  veo(c) = e1;  ber(c) = b1;
    lvl(c,:) = [min(hi) max(hi) min(lo) max(lo)];
end
dv0 = make_drive(case_rate(1), S, P);
prbs_bits = dv0.prbs_bits;
N_sym = S.N_sym;

fid = fopen(fullfile(od, 'summary.txt'), 'w');
ut.tee(fid, 'step 3, eye runs, %s. Q %d design, loss and coupling time x %.2f, loaded Q %.0f, half-width %.2f GHz\n', ...
    P.label, Q_target, res_alpha_scale, Q_ring, pr.hw_GHz);
ut.tee(fid, ' case  P(mW)  rate  dw(GHz) |  OMA     VEO    log10 BER | ones min..max (P_in)  zeros min..max (P_in)\n');
for c = 1:nC
    ut.tee(fid, ' %-4s  %4.0f  %5.1f  %+6.2f  | %+6.3f  %+6.3f  %8.2f  |  %5.3f .. %5.3f        %5.3f .. %5.3f\n', ...
        case_name{c}, case_P_mW(c), case_rate(c), case_dw(c), oma(c), veo(c), log10(max(ber(c),1e-30)), ...
        lvl(c,1), lvl(c,2), lvl(c,3), lvl(c,4));
end
fclose(fid);
save(fullfile(od, 'results.mat'), 'case_name', 'case_P_mW', 'case_rate', 'case_dw', 'Pt_all', ...
    'Mbit_all', 'oma', 'veo', 'ber', 'lvl', 'prbs_bits', 'N_sym', 'Q_target', 'Q_ring', ...
    'res_alpha_scale', 'res_tau_c_scale', '-v7');
fprintf('\nsaved to "%s"\n', od);
