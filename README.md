# SpikingStableRecurrentNonlinearNets

A spiking **Hodgkin–Huxley** recurrent network with spike-frequency adaptation
(SFA), short-term synaptic depression (STD), and short-term facilitation (STF),
built to reproduce the dynamical mechanisms of the continuous rate model
`SRNNModelCellTypes` (in `../FractionalResevoir`) in a biophysical spiking model —
so that **edge-of-chaos** behaviour (largest Lyapunov exponent via Benettin's
method) can be compared between the two. Connectivity and short-term-plasticity
parameters come from Campagnola et al. 2022.

## Quick start (MATLAB)

```matlab
cd SpikingStableRecurrentNonlinearNetsB
run scripts/setup_paths.m          % add src/ and scripts/ to the path

model = SRNNModelHH( ...
    'n', 80, 'T_range', [0 1000], ...      % ms
    'level_of_chaos', 1.5, ...             % edge-of-chaos scan knob
    'n_a', 1, 'n_b', 1, 'n_u', 1, ...      % SFA + STD + STF
    'lya_method', 'benettin', 'lya_dt', 20, ...
    'store_full_state', true);
model.build(); model.run(); model.plot();
fprintf('LLE = %.4g /ms\n', model.lya_results.LLE);
```

Or run the example: `run scripts/examples/demo_hh_network.m`.

## Layout

- `src/SRNNModelHH.m` — the spiking model (HH + SFA/STD/STF, Campagnola connectivity).
- `src/SRNNModelBase.m` — shared lifecycle + Benettin/QR Lyapunov (adapted from `../FractionalResevoir`).
- `src/integrators/integrate_hh_hybrid.m` — event-aware RK4 (spike-triggered jumps).
- `src/channels/` — Traub–Miles HH gating kinetics.
- `src/connectivity/` — Campagnola 2022 parameter loader + CSV data.
- `scripts/tests/` — sentinel tests; run all with `run scripts/tests/run_all_tests.m`.
- `docs/EquationsParametersDocs/hh_system_equations.md` — governing equations.

See `CLAUDE.md` for architecture, conventions, and the MATLAB path gotcha.
