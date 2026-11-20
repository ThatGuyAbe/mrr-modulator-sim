% runs every step in order. start matlab in this folder. outputs go to out/.
%   step0_checks         model checks                                  (a few min)
%   step1_sweeps         single-ring detuning sweeps, dw_opt and dw_lin (about 2 h)
%   step2_mechanism      thermal and carrier shifts                     (a few min)
%   step3_eyes           eye waveforms                                  (1 min)
%   step4_ring_dynamics  traces, self-pulsation, bistability            (about 20 min)
%   step5_convergence    time step, warm-up, data window, linear limit  (about 15 min)
%   step6_population     2000 rings, zdc, pac and af                    (about 10 h)
%   step7_tuning_error   residual tuning error                          (about 3 h)
% times are for 4 workers. steps 1 and 6 save checkpoints, so a stopped run picks up
% where it left off.
step0_checks
step1_sweeps
step2_mechanism
step3_eyes
step4_ring_dynamics
step5_convergence
step6_population
step7_tuning_error
