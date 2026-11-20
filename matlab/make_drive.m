function dv = make_drive(rate_Gbps, S, P, n_tau_ramp, n_tau_hold, N_sym, dt)
% nrz drive waveform for one bit rate. the prbs repeats through the warm-up (ramp,
% then hold, each rounded up to whole bits), then the N_sym-symbol data window
% follows. a one bit drives -4.5 V, a zero bit -0.5 V.
% the optional arguments override the run settings (warm-up and convergence checks).
% with dt the samples per bit are scaled from the 1 ps value, so the bit period and
% the decision instant stay where they are at 1 ps.
    if nargin < 4 || isempty(n_tau_ramp), n_tau_ramp = S.n_tau_ramp; end
    if nargin < 5 || isempty(n_tau_hold), n_tau_hold = S.n_tau_hold; end
    if nargin < 6 || isempty(N_sym),      N_sym      = S.N_sym;      end
    prbs_bits = gen_prbs(N_sym, S.prbs_order);
    T_bit = 1/(rate_Gbps*1e9);
    M_bit = round(T_bit/S.dt_real);
    if nargin >= 7 && ~isempty(dt)
        M_bit = round(M_bit * S.dt_real / dt);
    end
    N_ramp = ceil(n_tau_ramp * P.tau_th / T_bit);
    N_hold = ceil(n_tau_hold * P.tau_th / T_bit);
    if n_tau_hold == 0, N_hold = 0; end
    N_warmup    = N_ramp + N_hold;
    warm_bits   = prbs_bits(mod((0:N_warmup-1), N_sym) + 1);
    all_bits    = [warm_bits; prbs_bits];
    total_bits  = N_warmup + N_sym;
    total_steps = total_bits * M_bit;
    V_wave = zeros(total_steps,1);
    for b = 1:total_bits
        if all_bits(b)==1, v=P.V_low; else, v=P.V_high; end
        V_wave((b-1)*M_bit+1 : b*M_bit) = v;
    end
    dv = struct('V_wave', V_wave, 'M_bit', M_bit, 'N_ramp', N_ramp, 'N_hold', N_hold, ...
        'N_warmup', N_warmup, 'N_sym', N_sym, 'total_steps', total_steps, ...
        'prev_bit', warm_bits(end), 'prbs_bits', prbs_bits, 'rate', rate_Gbps, 'T_bit', T_bit);
end

function bits = gen_prbs(N, order)
% maximal-length prbs of the given order, N bits (fibonacci lfsr)
    switch order
        case 7,  taps = [7 6];
        case 9,  taps = [9 5];
        case 10, taps = [10 7];
        case 11, taps = [11 9];
        case 15, taps = [15 14];
        otherwise, error('Unsupported PRBS order %d', order);
    end
    reg  = ones(1, order);
    bits = zeros(N, 1);
    for k = 1:N
        bits(k) = reg(order);
        fb = reg(taps(1));
        for j = 2:numel(taps), fb = xor(fb, reg(taps(j))); end
        reg = [fb, reg(1:order-1)];
    end
end
