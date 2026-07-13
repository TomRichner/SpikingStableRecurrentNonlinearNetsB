# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this repo is

A spiking Hodgkin–Huxley (HH) recurrent network that reproduces the dynamical
mechanisms of the continuous rate model `SRNNModelCellTypes` (in the sibling repo
`../FractionalResevoir`) — spike-frequency adaptation (SFA), short-term synaptic
depression (STD), and short-term facilitation (STF) — so that **edge-of-chaos**
results (largest Lyapunov exponent via Benettin's method) can be compared between
the rate model and a biophysical spiking model. Connectivity and plasticity
parameters come from Campagnola et al. 2022 (Allen synaptic-physiology dataset).

Related read-only reference repos:
- `../FractionalResevoir` — the rate model, base-class lifecycle, Benettin/QR
  Lyapunov code, and the Campagnola loader this repo adapted/copied.
- `../CorticalSpreadDepolarizationModel/model_int.m` — source of the Traub–Miles
  HH gating kinetics and the singularity-handling idiom.

## Architecture

Handle-class lifecycle mirrored from FractionalResevoir:
`model = SRNNModelHH(...); model.build(); model.run(); model.plot();`

- `src/SRNNModelBase.m` — abstract base (adapted from FractionalResevoir). Owns
  the reflective name-value constructor, the `build → run` template, and the
  Lyapunov numerics (`benettin_algorithm_internal`, QR spectrum, Kaplan–Yorke)
  reused nearly verbatim so results are numerically comparable to the rate model.
  Concrete subclasses implement the abstract hooks: `set_defaults`,
  `build_network`, `build_stimulus`, `validate`, `get_params`,
  `decimate_and_unpack`, `eval_dynamics`, `eval_jacobian`.
- `src/SRNNModelHH.m` — the concrete spiking model (HH + SFA/STD/STF, Campagnola
  connectivity, jump-parameter packing, plotting).
- `src/integrators/integrate_hh_hybrid.m` — event-aware fixed-step RK4. Standalone
  and unit-testable; the model exposes it via the `run_integrator` seam
  (`obj.ode_solver`) so both `run()` and the Benettin reshoot apply identical
  spike jumps.
- `src/channels/hh_gating_rates.m`, `hh_gating_inf.m` — Traub–Miles α/β with
  removable-singularity patches, and steady-state gating for initialisation.
- `src/connectivity/load_campagnola_matrices.m` + `campagnola/` — copied verbatim
  from FractionalResevoir; the CSVs under `campagnola/` are the version-controlled
  source of truth.

Governing equations: `docs/EquationsParametersDocs/hh_system_equations.md`.

## Key conventions (specific to this codebase)

- **Units are milliseconds** (mV, µF/cm², mS/cm², µA/cm²). Campagnola time
  constants (s) are converted ×1000 to ms in `load_parameter_tables`. When
  setting `fs`, `T_range`, `lya_dt`, keep them all in ms.
- **Hybrid system.** HH neurons are a smooth ODE (no hard reset); the only
  discontinuities are synaptic/SFA jumps at presynaptic spikes. `eval_dynamics`
  returns ONLY the smooth between-spike RHS; all spike-triggered increments live
  in `integrate_hh_hybrid`. Do not add jumps to `dynamics_hh`.
- **Weight orientation:** `obj.W` is `(pre × post)`, `W(j,i)` = |conductance| from
  pre `j` to post `i`. This differs from `SRNNModelCellTypes` (post × pre). Signs
  are realized by synaptic reversal potential, not signed weights.
- **State layout** (cursor-walked, guarded blocks):
  `S = [V; m; h; n; a; b; p; g]`. Disabling a mechanism (`n_a/n_b/n_u = 0`)
  removes its block from `S` and from `N_sys_eqs`. `g` is always present.
- **Lyapunov:** only `lya_method='benettin'` is supported. `eval_jacobian` errors
  by design — the event-based (hybrid) QR spectrum would need saltation matrices
  (deferred). Set `obj.lya_dt` explicitly (ms); the base default (0.02) is only
  meaningful in seconds.
- No `arguments` blocks in constructors (matches FractionalResevoir); public
  property defaults + imperative `validate()` with namespaced error IDs + clamps.
- `snake_case` properties/vars, `PascalCase` classes, `test_<subject>.m` scripts
  printing a grep-able `<NAME>_PASS`/`<NAME>_FAIL` sentinel.

## Running / testing

There is no formal unit-test framework — tests are plain scripts under
`scripts/tests/` that print a sentinel. Run all of them with
`scripts/tests/run_all_tests.m` (prints `ALL_TESTS_PASS`/`ALL_TESTS_FAIL`).

Every entry-point script calls `setup_paths()` first. **Path gotcha:** the sibling
repo `FractionalResevoir` also has a `scripts/setup_paths.m`; if it is on the
MATLAB path, `setup_paths` may resolve to the wrong repo. When verifying via the
MATLAB MCP, first isolate the path:
```matlab
repoB = '.../SpikingStableRecurrentNonlinearNetsB';
restoredefaultpath;
addpath(genpath(fullfile(repoB,'src'))); addpath(genpath(fullfile(repoB,'scripts')));
cd(repoB);
```
Then `run('scripts/tests/run_all_tests.m')`. A single HH neuron rests quietly with
no drive and spikes under supra-threshold constant current; verify with
`scripts/tests/test_single_neuron.m` or the demo `scripts/examples/demo_hh_network.m`.

## Status / deferred

Core-first pass complete and verified: HH class + SFA/STD/STF + Campagnola
connectivity + Benettin LLE + tests. Deferred: QR full spectrum (needs jump
matrices), richer per-type plotting, and a parameter-sweep driver
(`ParamSpaceAnalysis2` analog) to map edge-of-chaos across `level_of_chaos` and
mechanism toggles. At the current default biophysical scaling the sampled LLEs are
negative (stable); locating the edge-of-chaos regime is the intended sweep work.

## Commit style

Match FractionalResevoir: plain commit messages, **no `Co-Authored-By: Claude`
trailer**. One commit per logical change.
