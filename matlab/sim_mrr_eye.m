function Pt_wave = sim_mrr_eye(E_in, V_wave, dw_n, dtau, ...
        total_steps, M_bit, N_ramp, N_warmup, N_data, ...
        mu_n, glin_n, cTPA_n, cFCA_n, cFCg, ...
        cTh_n, gamma_0, iTFC_n, iTTH_n, dn_dN, dn_dT, ...
        dw_NL_scale, dw_mod_n, iTRC_n, mu_raw) %#codegen
% through-port power of one run: warm-up, then the data window.
%   steps 1 .. N_ramp*M_bit              laser detuning ramps (quadratic) from
%                                        max(dw_n,0) + 15 half-widths on the blue side
%                                        to the target
%   up to N_warmup*M_bit                 held at the target
%   the next N_data*M_bit steps          data window, recorded
% each step is split into equal rk4 substeps when |detuning| * dtau would exceed 2.4,
% where a single 1 ps step is unstable (start of the ramp at Q = 3000, far-detuned
% as-fabricated rings). a diverged state is reset to zero.
% returns |E_in + i*mu_raw*a|^2 in W over the data window, N_data*M_bit samples.
%
% for speed in matlab this compiles with matlab coder:
%     args = repmat({0},1,24);  args{2} = coder.typeof(0,[Inf 1],[1 0]);
%     codegen sim_mrr_eye -args args
% then call sim_mrr_eye_mex in run_case.m. the results are identical.
    ar = 0;  ai = 0;  DN = 0;  DT = 0;  Vpn = V_wave(1);
    ramp_steps   = N_ramp * M_bit;
    warmup_steps = N_warmup * M_bit;
    dw_blue_far  = 15;                           % start 15 half-widths blue (normalized)
    dw_start     = max(dw_n, 0) + dw_blue_far;   % always come in from the blue side
    z_max        = 2.4;                          % largest |detuning|*step kept per RK4 step
    wmod         = abs(dw_mod_n) * max(abs(V_wave));   % largest modulation shift
    data_start   = warmup_steps + 1;
    Nd = N_data * M_bit;
    Pt_wave = zeros(Nd, 1);
    j = 0;
    for k = 1:total_steps
        Vk = V_wave(k);
        % ramp the detuning, then hold it at the target
        if k <= ramp_steps
            frac = k / ramp_steps;
            frac = frac * frac;                  % quadratic, slower at the start
            dw_k = dw_start + (dw_n - dw_start) * frac;
        else
            dw_k = dw_n;
        end
        nsub = max(1, ceil((abs(dw_k) + wmod) * dtau / z_max));
        if nsub == 1
            [ar, ai, DN, DT, Vpn] = mod_rk4_step(ar, ai, DN, DT, Vpn, ...
                E_in, Vk, dtau, ...
                dw_k, dw_NL_scale, dw_mod_n, ...
                mu_n, glin_n, cTPA_n, cFCA_n, cFCg, ...
                cTh_n, gamma_0, iTFC_n, iTTH_n, dn_dN, dn_dT, iTRC_n);
        else
            for is = 1:nsub
                [ar, ai, DN, DT, Vpn] = mod_rk4_step(ar, ai, DN, DT, Vpn, ...
                    E_in, Vk, dtau/nsub, ...
                    dw_k, dw_NL_scale, dw_mod_n, ...
                    mu_n, glin_n, cTPA_n, cFCA_n, cFCg, ...
                    cTh_n, gamma_0, iTFC_n, iTTH_n, dn_dN, dn_dT, iTRC_n);
            end
        end
        if ~isfinite(ar)
            ar = 0; ai = 0; DN = 0; DT = 0; Vpn = V_wave(k);
        end
        if k >= data_start
            j = j + 1;
            Pt_wave(j) = (E_in - mu_raw*ai)^2 + (mu_raw*ar)^2;
        end
    end
end
