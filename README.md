# mrr-modulator-sim
Nonlinear silicon microring modulator simulator in MATLAB, sweeping fabrication variation and Q to produce OMA/ER/IL metrics and eye diagrams.

# Silicon Microring Modulator Simulation

MATLAB scripts for simulating a voltage-driven silicon microring (MRR)
modulator with the full nonlinear cavity physics: two-photon absorption, free
carriers, and self-heating. The model is the 5-state ODE system from Giron
Castro (2024) plus the de Cea (2019) voltage modulation, integrated with RK4 in
normalised time (tau = gamma_0 * t).

The experiments look at how fabrication variation (resonance offset and Q
spread) hits the modulator, and whether tuning the offset back out ("corrected"
vs "as-fabricated") actually helps. Metrics are OMA, ER, IL, and eye diagrams,
all at the through port.

## Requirements

MATLAB with the Parallel Computing Toolbox. Every experiment script uses
`parpool`/`parfor`. If you don't have the toolbox, swap the `parfor` loops for
plain `for` and they'll still run, just slower.

## Files

Core (you don't run these directly, the experiment scripts call them):

- `mod_rk4_step.m`  one RK4 step of the 5-state model, plus the coupled ODE
  right-hand side. This is the physics engine.
- `sim_mrr_modulator.m`  runs one full PRBS simulation and returns OMA/ER/IL at
  both ports. Does the adiabatic warmup so the cavity lands on the upper
  bistability branch.
- `gen_prbs7.m`  127-bit PRBS-7 generator (used by the two sweep scripts).
- `sim_mrr_eye.m`  runs the same simulation as `sim_mrr_modulator.m` but returns
  the per-step through/drop port power trace instead of the scalar metrics. The
  two eye scripts use it to build the waveforms.

Experiment scripts (these are the ones you run):

- `three_mW.m`
- `sweep_Q.m`
- `eye_diagram.m`
- `eye_diagram_Q_sweep.m`

## three_mW.m

12 configs (3 powers x 4 data rates) at the nominal Q. For each config it does a
nominal detuning sweep, then runs 2000 rings with random offset + Q spread, both
as-fabricated and corrected (offset tuned to zero).

Edit near the top: `powers`, `datarates`, `N_rings`, `sigma_lam` (offset spread)
and `sigma_Q_rel` (Q spread). `rng(42)` fixes the ring draws so runs repeat.

Output: OMA-vs-ring and IL-vs-ring figures plus the nominal sweep figures in
`paper/`, a console summary table, and `multi_config_results_full.mat`.

## sweep_Q.m

The bigger sweep. Outer loop over `Q_values = 3000:500:15000`, inner loop over 18
configs (3 powers x 6 data rates). For each target Q it derives `tau_c` from the
nominal loss, then runs the same nominal-sweep + 2000-ring experiment.

Edit near the top: `Q_values`, `powers`, `datarates`, `N_rings`, the two spread
parameters.

Output: one folder per Q under `Q_sweep_results/Q_XXXXX/` with the detuning-sweep
and ring-index figures, a `summary_table_QXXXXX.txt`, and a `.mat`. At the end it
writes a grand summary table across all Q values.

## eye_diagram.m

Through-port eye diagrams for 12 configs (3 powers x 4 rates) at one Q. Four
cases per config: {as-fabricated, corrected} x {zero detuning, OMA-optimal
detuning}. The optimal detuning is found by sweeping the laser off the corrected
ring, then both rings are run at that setting.

Edit near the top: `powers`, `datarates`, `df_grid` (detuning search),
`res_alpha_scale` / `res_tau_c_scale` (the residual Q both rings carry), and
`fab_offset_GHz` (the uncorrected offset on the as-fab ring).

Output: four eye-grid PNGs in `paper/`, a signed OMA / eye-opening table, and
`paper/eye_results.mat`. Uses `sim_mrr_eye.m` for the waveforms.

## eye_diagram_Q_sweep.m

Same as `eye_diagram.m` but swept over `Q_values = [3000 6000 9000 12000 15000]`.
Each Q gets its own subfolder under `Q sweep eye diagram/` with the four eye
grids, a table, and a `.mat`. Uses `sim_mrr_eye.m` for the waveforms.

## Notes

- All four experiment scripts share the same parameter block (Giron Castro 2024
  Table 1 for the ring, de Cea 2019 for the modulator). If you change a physical
  constant, change it everywhere.
- "Corrected" means the resonance offset is tuned to zero; the residual Q spread
  stays. "As-fabricated" keeps the offset. Same ring, they only differ by the
  offset.
- OMA is reported normalised to input power. In the eye scripts the OMA and eye
  opening are signed: a negative value means the eye inverted or closed.
