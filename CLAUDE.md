# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

Pure MATLAB (R2025b); no build/compile step, no package manager. The preferred
way to run code here is the **matlab MCP** tools (`run_matlab_file`,
`evaluate_matlab_code`, `check_matlab_code`) against a live session launched via
the `/launch-matlab` skill; `matlab -batch` works headless as a fallback.

```matlab
run scripts/tests/run_all_tests.m        % run every test_*.m; prints ALL_TESTS_PASS/FAIL
run scripts/tests/test_benettin.m        % run ONE test (prints its own <NAME>_PASS/FAIL)
run scripts/examples/demo_hh_network.m   % build+run+plot a network, print the LLE
```
Headless single test:
`matlab -batch "cd('<repo>'); addpath(genpath('src'),genpath('scripts')); run('scripts/tests/test_benettin.m')"`

There is **no unit-test framework** — tests are plain scripts that print a
grep-able `<NAME>_PASS`/`<NAME>_FAIL` sentinel; success is verified by grepping
stdout for `_FAIL`. Static-analyze a file with the matlab MCP `check_matlab_code`.

**Path gotcha (important):** the sibling repo `../FractionalResevoir` also has a
`scripts/setup_paths.m`; if it is on the MATLAB path, `setup_paths` can resolve to
the WRONG repo and code won't be found. Isolate first when the session is shared:
```matlab
repoB = '<abs path to this repo>';
restoredefaultpath;
addpath(genpath(fullfile(repoB,'src'))); addpath(genpath(fullfile(repoB,'scripts')));
cd(repoB);
```

## What this repo is

Spiking Hodgkin–Huxley (HH) recurrent networks that reproduce the dynamical
mechanisms of the continuous rate models in the sibling repo `../FractionalResevoir`
— spike-frequency adaptation (SFA), short-term synaptic depression (STD), and
short-term facilitation (STF) — so that **edge-of-chaos** results (largest Lyapunov
exponent via Benettin's method) can be compared between rate and biophysical spiking
models. There are **two concrete spiking models**, paralleling the two rate models:
- `SRNNModelHH` ↔ `SRNNModelCellTypes`: 4 Campagnola cell types, data-driven
  connectivity, per-(pre,post-type) single-timescale STD/STF.
- `SRNNModelHHEI` ↔ `SRNNModel2`: 2 populations (E/I), RMT/Harris connectivity,
  per-population multi-timescale SFA (K-current) and STD, no STF.
Connectivity/plasticity for the cell-type model come from Campagnola et al. 2022.

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
- `src/SRNNModelHH.m` — Campagnola cell-type model (HH + SFA/STD/STF, per-(pre,post-type)
  synaptic resources, Bernoulli+magnitude connectivity, sign via reversal potential).
- `src/SRNNModelHHEI.m` — E/I RMT model comparable to `SRNNModel2`: `RMTMatrix`
  connectivity (wrong-sign entries **zeroed**, magnitude + Dale sign via reversal;
  `indegree`/α is the sparsity knob), **per-population timescale counts**
  (`n_a_vec`, `n_b_vec` are K-vectors → ragged per-type state blocks), SFA as a
  K-current with **DC-gain-balanced pools** (increment `Δₗ = a_incr0/τₗ` so logspaced
  timescales contribute equally, matching the rate model's `aₗ*=r`; `a_incr0=1000/target_rate`),
  and multi-timescale STD with **product efficacy** (`Πₘ bₘ`), per-spike depletion
  `bₘ-=p0·bₘ`; STF absent. `p0` maps to the rate model's `τ_rel` at the target rate.
- `src/integrators/integrate_hh_events.m` — generic event-aware fixed-step RK4
  (threshold detection + hysteresis); takes a `jump_fn(y,spiked)` closure.
  `integrate_hh_hybrid.m` (→ `apply_hh_jumps`) and `SRNNModelHHEI` (→ `apply_hhei_jumps`)
  both drive it, so `run()` and the Benettin reshoot apply identical spike jumps.
- `src/channels/hh_gating_rates.m`, `hh_gating_inf.m` — Traub–Miles α/β with
  removable-singularity patches, and steady-state gating for initialisation.
- `src/connectivity/` — `load_campagnola_matrices.m` + `campagnola/` CSVs (cell-type
  model) and `RMTMatrix.m` (E/I model), copied from FractionalResevoir.

Governing equations: `docs/EquationsParametersDocs/hh_system_equations.md` (cell-type).
`SRNNModelHHEI` is documented in its class header + the plan/discussion.

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
  meaningful in seconds. The LLE is in 1/ms; `lya_results.LLE_per_s` reports it in
  1/s (via `time_units_per_second=1000`) for direct comparison with the
  seconds-based rate model.
- No `arguments` blocks in constructors (matches FractionalResevoir); public
  property defaults + imperative `validate()` with namespaced error IDs + clamps.
- `snake_case` properties/vars, `PascalCase` classes, `test_<subject>.m` scripts
  printing a grep-able `<NAME>_PASS`/`<NAME>_FAIL` sentinel.

## Status / deferred

Both spiking models complete and verified (`ALL_TESTS_PASS`): `SRNNModelHH`
(Campagnola) and `SRNNModelHHEI` (E/I RMT) + Benettin LLE + tests + demos. The E/I
RMT model gives positive LLEs at moderate `level_of_chaos` (chaotic), the Campagnola
model negative (stable) at default scaling. Deferred: QR full spectrum (needs jump
matrices); a parameter-sweep driver (`ParamSpaceAnalysis2` analog) to map
edge-of-chaos across `level_of_chaos`; hand-tuning the E/I operating point
(`bias`, `c_type`, `p0_type`, `g_syn_scale` at ~5 Hz — formulas in the class header,
not auto-solved); long-run Benettin over second-scale SFA timescales (the dominant
cost — the multi-scale ms-spikes / second-adaptation stiffness needs runs of tens of
seconds; `SRNNModelHHEI` demo/tests use shortened τ for tractability).

## Commit style

Match FractionalResevoir: plain commit messages, **no `Co-Authored-By: Claude`
trailer**. One commit per logical change.
