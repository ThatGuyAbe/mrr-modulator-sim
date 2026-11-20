% step 7, residual tuning error. run after step 6.
% ring i gets an extra laser detuning eps_i = sigma_res * z_i (seed 2026, saved draw),
% added to both targets:  zdc: dw = eps_i    pac: dw = dw_opt + eps_i
% eps > 0 moves the laser to the blue, eps < 0 to the red. sigma_res = 0 is the step 6
% result itself (read from pop_results, not rerun). the first n_chk rings are rerun at
% sigma_res = 0 as a check.
% configurations: 4 mW at Q 9000 and 15000, 2.5 and 15 Gb/s, sigma_res = 1 and 2 GHz.
% output: out/step7_tuning_error/summary.txt, results.mat and one checkpoint per
% configuration. cost: 4 configs x 2 sigma x 2 targets x 2000 rings = 32,000 runs,
% about 3 h on 4 workers.

clear; clc; close all;
[P, S] = mrr_params();
ut = utils();
od = ut.out_folder(S, 'step7_tuning_error');
nW = ut.start_pool();
T1 = ut.load_targets(S);
R  = ring_population(S, P);

Ph = S.P_ref;                      % 4 mW
C = [ 9000 Ph  2.5;
      9000 Ph 15;
     15000 Ph  2.5;
     15000 Ph 15];
sigma_res = [1 2];                 % GHz (sigma_res = 0 comes from pop_results)
n_chk     = 10;                    % rings rerun at sigma_res = 0 as a check
N = S.N_rings;  z = R.z_res;
nCfg = size(C,1);  nSig = numel(sigma_res);

OMA = nan(N, nSig+1, 2, nCfg);  VEO = OMA;  BER = OMA;      % (ring, sigma incl. 0, target, config)
dwo = zeros(nCfg,1);  chk0 = nan(nCfg,1);  t_cfg = zeros(nCfg,1);
for k = 1:nCfg
    q = C(k,1);  Pm = C(k,2);  rate = C(k,3);
    dwo(k) = T1.get(q, Pm, rate);
    D = load(fullfile(S.out_root, 'step6_population', sprintf('pop_results_Q%d.mat', q)));
    ci = find(D.configs(:,1) == Pm & abs(D.configs(:,2) - rate) < 1e-9, 1);
    if abs(D.df_opt(ci) - dwo(k)) > 1e-9, error('dw_opt in pop_results differs from step 1'); end
    OMA(:,1,1,k) = D.OMA_cor(:,ci);  VEO(:,1,1,k) = D.EO_cor(:,ci);  BER(:,1,1,k) = D.BER_cor(:,ci);
    OMA(:,1,2,k) = D.OMA_opt(:,ci);  VEO(:,1,2,k) = D.EO_opt(:,ci);  BER(:,1,2,k) = D.BER_opt(:,ci);
    dv = make_drive(rate, S, P);
    pc = cell(N,1);
    for r = 1:N, pc{r} = ring_params(P, q, R.alpha_scale(r), R.tau_c_scale(r)); end
    % check at sigma_res = 0
    ref0 = [D.OMA_cor(1:n_chk,ci), D.EO_cor(1:n_chk,ci), D.OMA_opt(1:n_chk,ci), D.EO_opt(1:n_chk,ci)];
    dwk  = dwo(k);
    c0 = zeros(n_chk, 1);
    parfor r = 1:n_chk
        [~, o2, e2] = run_case('eye', Pm*1e-3, 0, pc{r}, dv, P, S);
        [~, o3, e3] = run_case('eye', Pm*1e-3, dwk, pc{r}, dv, P, S);
        c0(r) = max(abs([o2 e2 o3 e3] - ref0(r,:)));
    end
    chk0(k) = max(c0);
    ck = fullfile(od, sprintf('ckpt_cfg%d.mat', k));
    if exist(ck, 'file') == 2
        L = load(ck);
        if isequal(L.cfg, C(k,:)) && isequal(L.sigma_res, sigma_res) && L.dwk == dwk && strcmp(L.label, P.label)
            OMA(:,2:end,:,k) = L.om;  VEO(:,2:end,:,k) = L.eo;  BER(:,2:end,:,k) = L.bb;  t_cfg(k) = L.t;
            fprintf('config %d loaded from checkpoint\n', k);
            continue;
        end
    end
    [jr, js, jt] = ndgrid(1:N, 1:nSig, 1:2);  jr = jr(:);  js = js(:);  jt = jt(:);
    dw_job = (jt == 2) * dwk + sigma_res(js).' .* z(jr);
    nJ = numel(jr);  om = zeros(nJ,1);  eo = om;  bb = om;
    tk = tic;
    parfor j = 1:nJ
        [~, o1, e1, b1] = run_case('eye', Pm*1e-3, dw_job(j), pc{jr(j)}, dv, P, S);
        om(j) = o1;  eo(j) = e1;  bb(j) = b1;
    end
    t = toc(tk);
    om = reshape(om, N, nSig, 2);  eo = reshape(eo, N, nSig, 2);  bb = reshape(bb, N, nSig, 2);
    cfg = C(k,:);  label = P.label;
    save(ck, 'om', 'eo', 'bb', 'cfg', 'sigma_res', 'label', 'dwk', 't', '-v7');
    OMA(:,2:end,:,k) = om;  VEO(:,2:end,:,k) = eo;  BER(:,2:end,:,k) = bb;  t_cfg(k) = t;
    fprintf('config %d of %d done in %.1f min\n', k, nCfg, t/60);
end

%% summary
sig_all = [0 sigma_res];
names = {'ZDC', 'PAC'};
fid = fopen(fullfile(od, 'summary.txt'), 'w');
ut.tee(fid, 'step 7, residual tuning error, %s. %d rings, error seed 2026 (saved draw), sigma_res = 0 (pop_results) and %s GHz\n', ...
    P.label, N, mat2str(sigma_res));
ut.tee(fid, 'check: rings 1-%d rerun at sigma_res = 0 vs pop_results, largest difference %.1e (must be 0)\n', n_chk, max(chk0));
ut.tee(fid, 'z: mean %.4f, std %.4f, %d negative (laser moved red), %d positive (blue)\n', mean(z), std(z), sum(z < 0), sum(z > 0));
for k = 1:nCfg
    pk = ring_params(P, C(k,1));
    ut.tee(fid, '\nconfig %d. Q %d, %g mW, %g Gb/s, dw_opt %+.2f GHz, half-width %.2f GHz\n', k, C(k,1), C(k,2), C(k,3), dwo(k), pk.hw_GHz);
    ut.tee(fid, '  target sigma | mean OMA  mean VEO  log10 BER | invalid  (of eps<0, of eps>0) |  VEO10  | smallest |eps| that breaks a ring valid at sigma 0: red, blue (GHz)\n');
    for t = 1:2
        for s = 1:nSig+1
            o = OMA(:,s,t,k);  e = VEO(:,s,t,k);  b = BER(:,s,t,k);
            bad = e < 0;  bad0 = VEO(:,1,t,k) < 0;  es = sort(e);
            brk = bad & ~bad0;  ep = sig_all(s) * z;
            r_red = min(abs(ep(brk & ep < 0)));  r_blu = min(abs(ep(brk & ep > 0)));
            if isempty(r_red), r_red = NaN; end
            if isempty(r_blu), r_blu = NaN; end
            ut.tee(fid, '  %-4s  %4.1f  | %+7.3f   %+7.3f   %7.2f   |  %4d    (%4d, %4d)          | %+7.3f | %6.2f  %6.2f\n', ...
                names{t}, sig_all(s), mean(o), mean(e), mean(log10(max(b, 1e-30))), sum(bad), sum(bad & z < 0), ...
                sum(bad & z > 0), es(ceil(0.1*N)), r_red, r_blu);
        end
    end
end
fclose(fid);
save(fullfile(od, 'results.mat'), 'OMA', 'VEO', 'BER', 'C', 'sigma_res', 'sig_all', 'dwo', 'z', 'chk0', 't_cfg', '-v7');
fprintf('\nsaved to "%s"\n', od);
