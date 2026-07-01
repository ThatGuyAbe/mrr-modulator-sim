function [Pt_wave, Pd_wave] = sim_mrr_eye(E_in, V_wave, dw_n, dtau, ...
        total_steps, M_bit, N_warmup, N_data, ...
        mu_n, glin_n, cTPA_n, cFCA_n, cFCg, ...
        cTh_n, gamma_0, iTFC_n, iTTH_n, dn_dN, dn_dT, ...
        dw_NL_scale, dw_mod_n, iTRC_n, mu_raw)
%SIM_MRR_EYE  time-resolved through/drop port power for an eye diagram.
%
%   Same 5-state RK4 model as sim_mrr_modulator, but instead of averaging
%   power per bit it keeps the power at every step over the data window. That
%   holds onto the intra-symbol shape (edges, ringing, overshoot) that the eye
%   is built from.
%
%   Returns two column vectors of length N_data*M_bit:
%     Pt_wave  through-port power  |E_in + i*mu_raw*a|^2   [W]
%     Pd_wave  drop-port power     (2/tau_c)*|a|^2         [W]
%
%   Through port is the one the OIP 2026 metrics use.

    % init state
    ar = 0;  ai = 0;  DN = 0;  DT = 0;  Vpn = V_wave(1);

    data_start = N_warmup * M_bit + 1;
    Nd = N_data * M_bit;
    Pt_wave = zeros(Nd, 1);
    Pd_wave = zeros(Nd, 1);
    j = 0;

    for k = 1:total_steps
        Vk = V_wave(k);

        [ar, ai, DN, DT, Vpn] = mod_rk4_step(ar, ai, DN, DT, Vpn, ...
            E_in, Vk, dtau, ...
            dw_n, dw_NL_scale, dw_mod_n, ...
            mu_n, glin_n, cTPA_n, cFCA_n, cFCg, ...
            cTh_n, gamma_0, iTFC_n, iTTH_n, dn_dN, dn_dT, iTRC_n);

        % catch divergence
        if ~isfinite(ar)
            ar = 0; ai = 0; DN = 0; DT = 0; Vpn = V_wave(k);
        end

        % record only inside the data window
        if k >= data_start
            j = j + 1;
            Pt_wave(j) = (E_in - mu_raw*ai)^2 + (mu_raw*ar)^2;
            Pd_wave(j) = mu_raw^2 * (ar^2 + ai^2);
        end
    end
end
