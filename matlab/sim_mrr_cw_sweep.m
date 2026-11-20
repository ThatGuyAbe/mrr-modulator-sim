function R = sim_mrr_cw_sweep(E_in, Vc, dw_list, cFCg_list, dtau, ...
        hold_steps, rec_steps, settle0_steps, ...
        mu_n, glin_n, cTPA_n, cFCA_n, ...
        cTh_n, gamma_0, iTFC_n, iTTH_n, dn_dN, dn_dT, ...
        dw_NL_scale, dw_mod_n, iTRC_n, mu_raw) %#codegen
% quasi-static detuning sweep with a cw input and a constant junction voltage Vc.
% the laser detuning steps through dw_list (normalized, in sweep order) inside one
% run, so the ring carries its state from point to point and can stay on a branch.
% each point is held for hold_steps, the first one settle0_steps longer from a cold
% start. cFCg depends on the laser frequency, so it is passed per point.
% one row per point: [Pt, U, DN, DT, Vpn, dwn] at the end of the hold, then
% [Pt_min, Pt_max, DN_min, DN_max] over the last rec_steps of the hold.
    ar = 0;  ai = 0;  DN = 0;  DT = 0;  Vpn = Vc;
    n = numel(dw_list);
    R = zeros(n, 10);
    for i = 1:n
        dw_k = dw_list(i);
        cFCg = cFCg_list(i);
        ns = hold_steps;
        if i == 1
            ns = ns + settle0_steps;
        end
        pmin = inf;  pmax = -inf;  nmin = inf;  nmax = -inf;
        for k = 1:ns
            [ar, ai, DN, DT, Vpn] = mod_rk4_step(ar, ai, DN, DT, Vpn, ...
                E_in, Vc, dtau, ...
                dw_k, dw_NL_scale, dw_mod_n, ...
                mu_n, glin_n, cTPA_n, cFCA_n, cFCg, ...
                cTh_n, gamma_0, iTFC_n, iTTH_n, dn_dN, dn_dT, iTRC_n);
            if ~isfinite(ar)
                ar = 0; ai = 0; DN = 0; DT = 0; Vpn = Vc;
            end
            if k > ns - rec_steps
                Pt = (E_in - mu_raw*ai)^2 + (mu_raw*ar)^2;
                pmin = min(pmin, Pt);  pmax = max(pmax, Pt);
                nmin = min(nmin, DN);  nmax = max(nmax, DN);
            end
        end
        Pt  = (E_in - mu_raw*ai)^2 + (mu_raw*ar)^2;
        U   = ar^2 + ai^2;
        dwn = dw_k + dw_NL_scale * (DN*dn_dN + DT*dn_dT) + dw_mod_n * Vpn;
        R(i, :) = [Pt, U, DN, DT, Vpn, dwn, pmin, pmax, nmin, nmax];
    end
end
