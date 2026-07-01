function [ar,ai,Nn,Tn,Vn] = mod_rk4_step(ar,ai,Nv,Tv,Vv, ...
        Ein, Vext, dtau, ...
        dw_n, dw_NL_s, dw_mod_n, ...
        mu_n, glin_n, cTPA_n, cFCA_n, cFCg_n, ...
        cTh_n, g0, iTFC_n, iTTH_n, dnN, dnT, iTRC_n)
%MOD_RK4_STEP  one RK4 step of the 5-state silicon MRR modulator model.
%
%   state: [ar, ai, DN, DT, Vpn]
%     ar, ai   real/imag parts of the intracavity field a
%     DN       excess free-carrier density [m^-3]
%     DT       temperature rise above ambient [K]
%     Vpn      voltage across the PN junction [V]
%
%   Ein and Vext are held fixed over the step (sample-and-hold):
%     Ein      CW input field amplitude [sqrt(W)]
%     Vext     drive voltage [V]

    h = 0.5 * dtau;

    [r1,i1,n1,t1,v1] = mod_cavity_rhs(ar, ai, Nv, Tv, Vv, ...
        Ein, Vext, dw_n, dw_NL_s, dw_mod_n, ...
        mu_n, glin_n, cTPA_n, cFCA_n, cFCg_n, ...
        cTh_n, g0, iTFC_n, iTTH_n, dnN, dnT, iTRC_n);

    [r2,i2,n2,t2,v2] = mod_cavity_rhs( ...
        ar + h*r1, ai + h*i1, max(Nv + h*n1, 0), max(Tv + h*t1, 0), Vv + h*v1, ...
        Ein, Vext, dw_n, dw_NL_s, dw_mod_n, ...
        mu_n, glin_n, cTPA_n, cFCA_n, cFCg_n, ...
        cTh_n, g0, iTFC_n, iTTH_n, dnN, dnT, iTRC_n);

    [r3,i3,n3,t3,v3] = mod_cavity_rhs( ...
        ar + h*r2, ai + h*i2, max(Nv + h*n2, 0), max(Tv + h*t2, 0), Vv + h*v2, ...
        Ein, Vext, dw_n, dw_NL_s, dw_mod_n, ...
        mu_n, glin_n, cTPA_n, cFCA_n, cFCg_n, ...
        cTh_n, g0, iTFC_n, iTTH_n, dnN, dnT, iTRC_n);

    [r4,i4,n4,t4,v4] = mod_cavity_rhs( ...
        ar + dtau*r3, ai + dtau*i3, max(Nv + dtau*n3, 0), max(Tv + dtau*t3, 0), Vv + dtau*v3, ...
        Ein, Vext, dw_n, dw_NL_s, dw_mod_n, ...
        mu_n, glin_n, cTPA_n, cFCA_n, cFCg_n, ...
        cTh_n, g0, iTFC_n, iTTH_n, dnN, dnT, iTRC_n);

    s  = dtau / 6;
    ar = ar + s * (r1 + 2*r2 + 2*r3 + r4);
    ai = ai + s * (i1 + 2*i2 + 2*i3 + i4);
    Nn = max(Nv + s * (n1 + 2*n2 + 2*n3 + n4), 0);
    Tn = max(Tv + s * (t1 + 2*t2 + 2*t3 + t4), 0);
    Vn = Vv + s * (v1 + 2*v2 + 2*v3 + v4);
end


% the coupled ODE right-hand side, kept in its own function for speed
function [dar, dai, dN, dTv, dVpn] = mod_cavity_rhs( ...
        ar, ai, Nv, Tv, Vpn, ...
        Ein, Vext, ...
        dw_n, dw_NL_s, dw_mod_n, ...
        mu_n, glin_n, cTPA_n, cFCA_n, cFCg_n, ...
        cTh_n, g0, iTFC_n, iTTH_n, dnN, dnT, iTRC_n)
%MOD_CAVITY_RHS  the coupled ODEs in normalised time tau = gamma_0 * t.
%
%   Giron Castro (2024) Eqs. 1-6 plus de Cea (2019) Eqs. 1-3 for the
%   voltage-driven modulation.
%
%   total detuning = (omega_p - omega_0) - dW_NL - dW_mod
%     NL shift (free carriers + thermal):
%       dW_NL  = -(omega_0/n_Si)*(DN*dn_dN + DT*dn_dT)
%     modulation shift:
%       dW_mod = dW0_dVpn * Vpn
%     normalised:
%       dw_tot = dw_n + dw_NL_s*(DN*dn_dN + DT*dn_dT) + dw_mod_n*Vpn
%     with dw_mod_n = -dW0_dVpn / gamma_0.

    a2 = ar^2 + ai^2;

    % detuning
    dw_NL  = dw_NL_s * (Nv*dnN + Tv*dnT);   % free-carrier + thermal
    dw_mod = dw_mod_n * Vpn;                  % voltage modulation
    dw_tot = dw_n + dw_NL + dw_mod;

    % total cavity loss (normalised)
    gTPA = cTPA_n * a2;
    gFCA = cFCA_n * Nv;
    gtot = 1 + gTPA + gFCA;   % the 1 is (gamma_lin + gamma_coup)/gamma_0

    % cavity field da/dtau (CW input, no feedback).
    % i*mu_n*Ein drives the imag part, the real part has no direct drive.
    dar = -dw_tot * ai  -  gtot * ar;
    dai =  dw_tot * ar  -  gtot * ai  +  mu_n * Ein;

    % free-carrier density dDN/dtau
    dN = -iTFC_n * Nv  +  cFCg_n * a2^2;

    % temperature dDT/dtau
    Pn  = (glin_n + gTPA + gFCA) * a2;       % absorbed power (normalised)
    dTv = -iTTH_n * Tv  +  cTh_n * g0 * Pn;

    % junction voltage dVpn/dtau, first-order RC filter
    dVpn = iTRC_n * (Vext - Vpn);
end
