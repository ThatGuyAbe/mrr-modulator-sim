function bits = gen_prbs7()
%GEN_PRBS7  PRBS-7 sequence, 2^7 - 1 = 127 bits.
%   7-bit Fibonacci LFSR, taps at 6 and 7 (x^7 + x^6 + 1).
    reg = ones(1, 7);
    N   = 127;
    bits = zeros(N, 1);
    for k = 1:N
        bits(k) = reg(7);
        newbit  = xor(reg(6), reg(7));
        reg     = [newbit, reg(1:6)];
    end
end
