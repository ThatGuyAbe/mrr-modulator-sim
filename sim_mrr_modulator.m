function [oma_drop, er_drop, il_drop, oma_thr, er_thr, il_thr] = ...
        sim_mrr_modulator(E_in, V_wave, dw_n, dtau, total_steps, ...
        M_bit, N_warmup, N_bits, prbs_bits, ...
        mu_n, glin_n, cTPA_n, cFCA_n, cFCg, ...
        cTh_n, gamma_0, iTFC_n, iTTH_n, dn_dN, dn_dT, ...
        dw_NL_scale, dw_mod_n, iTRC_n, P_in, mu_raw)
%SIM_MRR_MODULATOR  run a voltage-driven silicon MRR modulator, return the
%   OMA/ER/IL at both ports.
%
%   Adiabatic warmup: the laser detuning is ramped in from the far blue side
%   toward the target, so the cavity settles on the upper branch of the
%   thermal bistability instead of dropping to the lower one. Same trick as
%   de Cea (2019) Section 4.1.
%
%   drop port:    P_drop = (2/tau_c) * |a|^2
%   through port: P_thr  = |E_in + i*mu_raw*a|^2

    % init state
    ar = 0;  ai = 0;  DN = 0;  DT = 0;  Vpn = V_wave(1);

    % adiabatic sweep setup.
    % during warmup, ramp the detuning from far blue toward the target so the
    % cavity tracks the upper bistability branch.
    % "blue" = laser above resonance = positive dw_n.
    warmup_steps  = N_warmup * M_bit;
    dw_blue_far   = 15;         % start 15 linewidths to the blue (normalised)

    % always start on the blue side (large positive dw) and sweep toward target
    dw_start = max(dw_n, 0) + dw_blue_far;

    % average power per data bit
    Pd_bits  = zeros(N_bits, 1);
    Pt_bits  = zeros(N_bits, 1);

    data_start = warmup_steps + 1;
    acc_Pd  = 0;
    acc_Pt  = 0;
    cnt     = 0;
    cur_bit = 0;

    for k = 1:total_steps
        Vk = V_wave(k);

        % detuning: ramp during warmup, then hold fixed
        if k <= warmup_steps
            frac  = k / warmup_steps;              % 0 to 1
            frac  = frac * frac;                    % quadratic (slower start)
            dw_k  = dw_start + (dw_n - dw_start) * frac;
        else
            dw_k = dw_n;
        end

        % RK4 step
        [ar, ai, DN, DT, Vpn] = mod_rk4_step(ar, ai, DN, DT, Vpn, ...
            E_in, Vk, dtau, ...
            dw_k, dw_NL_scale, dw_mod_n, ...
            mu_n, glin_n, cTPA_n, cFCA_n, cFCg, ...
            cTh_n, gamma_0, iTFC_n, iTTH_n, dn_dN, dn_dT, iTRC_n);

        % catch divergence
        if ~isfinite(ar)
            ar = 0; ai = 0; DN = 0; DT = 0; Vpn = V_wave(k);
        end

        % accumulate power over the data window
        if k >= data_start
            % drop port: P_d = (2/tau_c)*|a|^2
            Pd = mu_raw^2 * (ar^2 + ai^2);

            % through port: P_thr = |E_in + i*mu*a|^2
            Pt = (E_in - mu_raw*ai)^2 + (mu_raw*ar)^2;

            acc_Pd = acc_Pd + Pd;
            acc_Pt = acc_Pt + Pt;
            cnt    = cnt + 1;

            if cnt == M_bit
                cur_bit = cur_bit + 1;
                Pd_bits(cur_bit) = acc_Pd / M_bit;
                Pt_bits(cur_bit) = acc_Pt / M_bit;
                acc_Pd = 0;  acc_Pt = 0;  cnt = 0;
            end
        end
    end

    % pull out OMA, ER, IL
    idx1 = (prbs_bits == 1);
    idx0 = (prbs_bits == 0);

    if sum(idx1) < 2 || sum(idx0) < 2
        oma_drop = 0; er_drop = 0; il_drop = 0;
        oma_thr  = 0; er_thr  = 0; il_thr  = 0;
        return;
    end

    % drop port
    mu1_d = mean(Pd_bits(idx1));
    mu0_d = mean(Pd_bits(idx0));
    Ph_d  = max(mu1_d, mu0_d);
    Pl_d  = min(mu1_d, mu0_d);

    oma_drop = (Ph_d - Pl_d) / P_in;
    er_drop  = safe_ratio_dB(Ph_d, Pl_d);
    il_drop  = -10*log10(max(Ph_d, 1e-30) / P_in);

    % through port
    mu1_t = mean(Pt_bits(idx1));
    mu0_t = mean(Pt_bits(idx0));
    Ph_t  = max(mu1_t, mu0_t);
    Pl_t  = min(mu1_t, mu0_t);

    oma_thr = (Ph_t - Pl_t) / P_in;
    er_thr  = safe_ratio_dB(Ph_t, Pl_t);
    il_thr  = -10*log10(max(Ph_t, 1e-30) / P_in);
end

function dB = safe_ratio_dB(Phigh, Plow)
    if Plow > 1e-30
        dB = 10*log10(Phigh / Plow);
    else
        dB = 30;
    end
    if ~isfinite(dB), dB = 0; end
end