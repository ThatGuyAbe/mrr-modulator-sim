function [W, D] = sim_mrr_fulltrace(E_in, V_wave, dw_n, dtau, ...
        total_steps, M_bit, N_ramp, N_warmup, N_data, ...
        mu_n, glin_n, cTPA_n, cFCA_n, cFCg, ...
        cTh_n, gamma_0, iTFC_n, iTTH_n, dn_dN, dn_dT, ...
        dw_NL_scale, dw_mod_n, iTRC_n, mu_raw, dec, keep_data) %#codegen
% the same run as sim_mrr_eye.m (same loop, substeps and reset), also recording the
% internal states.
%   W  every dec-th step of the whole run, one row per kept step:
%      [step, Pt (W), U (J), DN (m^-3), DT (K), Vpn (V), dwn]
%   D  every step of the data window if keep_data is true: [Pt, U, DN, DT, Vpn, dwn]
% dwn is the instantaneous detuning over gamma_0, in half-widths, + = laser on the
% blue side of the hot, biased resonance.
    ar = 0;  ai = 0;  DN = 0;  DT = 0;  Vpn = V_wave(1);
    ramp_steps   = N_ramp * M_bit;
    warmup_steps = N_warmup * M_bit;
    dw_blue_far  = 15;
    dw_start     = max(dw_n, 0) + dw_blue_far;
    z_max        = 2.4;                          % largest |detuning|*step kept per RK4 step
    wmod         = abs(dw_mod_n) * max(abs(V_wave));   % largest modulation shift
    data_start   = warmup_steps + 1;
    Nd = N_data * M_bit;
    nW = floor(total_steps / dec);
    W = zeros(nW, 7);
    if keep_data
        D = zeros(Nd, 6);
    else
        D = zeros(0, 6);
    end
    j = 0;  iw = 0;
    for k = 1:total_steps
        Vk = V_wave(k);
        if k <= ramp_steps
            frac = k / ramp_steps;
            frac = frac * frac;
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
        Pt  = (E_in - mu_raw*ai)^2 + (mu_raw*ar)^2;
        U   = ar^2 + ai^2;
        dwn = dw_k + dw_NL_scale * (DN*dn_dN + DT*dn_dT) + dw_mod_n * Vpn;
        if mod(k, dec) == 0
            iw = iw + 1;
            W(iw, :) = [k, Pt, U, DN, DT, Vpn, dwn];
        end
        if keep_data && k >= data_start
            j = j + 1;
            D(j, :) = [Pt, U, DN, DT, Vpn, dwn];
        end
    end
end
