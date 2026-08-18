# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Master's thesis (University of Stuttgart, co-supervised by Kyoto University) on **mission-centric requirements analysis for Model Predictive Control in Spacecraft Rendezvous and Proximity Operations**. The code serves two purposes:

1. **Algorithm implementation and comparison** — LQR/LQI as baselines, standard MPC (regulation + tracking), and eMPC (regulation + NMC orbit acquisition/maintenance) are all implemented. eMPC is the primary focus and main use case of the thesis.
2. **Performance analysis and validation** — simulation outputs (JSON exports, figures) validate the requirements framework developed in the thesis. A **model-validity envelope** study additionally quantifies where the CW-based controller model breaks under eccentricity, separation, and perturbations, using an Aerospace-Toolbox high-fidelity truth plant.

The dynamics are modeled using Clohessy-Wiltshire-Hill (CWH) equations and more general relative motion formulations. The controllers always use a CW internal model; fidelity is stressed by swapping the *plant*, not the controller model.

## Getting Started

```matlab
% Once per MATLAB/VS Code session — opens the Project and loads all paths
setup.m
```

`setup.m` opens `2026-docking-empc-collab.prj` and `addpath(genpath(...))`s the `simulations/` and `tests/` trees (their subfolders are NOT on the path until this runs). No compile or build step.

**Toolboxes:** Control System Toolbox (`dlqr`, `c2d`), Optimization Toolbox (`quadprog`), and — for the model-validity truth plant — **Aerospace Toolbox** (`gravityzonal`, `atmosnrlmsise00`). Figure export requires `tools/matlab2tikz-master/`.

## Running Simulations

Simulations are self-contained scripts under `simulations/`, organized by topic:

```
simulations/
  foundations/     % baselines & building blocks
    sim_lqr.m              % LQR baseline
    sim_lqi.m              % LQI (LQR + integral action on position)
    sim_mpc_regulation.m   % MPC driving state to origin
    sim_mpc_tracking.m     % MPC tracking a reference trajectory
    sim_rpo.m              % Relative periodic orbit analysis
    sim_two_body.m         % 2-body propagator validation
  empc/            % economic MPC analysis
    sim_empc_regulation.m  % eMPC regulation — min-fuel, vs standard MPC
    sim_empc_nmc.m         % eMPC NMC orbit acquisition/maintenance
    sim_nmc_free_drift.m   % NMC natural (drift-free) motion check
  adrios/          % ADRIOS mission scenarios
    sim_scenario1_docking.m
    sim_scenario2_flyaround.m
    sim_scenario4_closing.m
    sim_adrios_combinations.m
    sim_adrios_full_mission.m
  model_validity/  % "where the CW model breaks" envelope study
    sim_model_validity_sweep.m       % Layer 1a: eccentricity × separation (open-loop)
    sim_perturbation_sweep.m         % Layer 1b: altitude × BC-ratio, J2+drag
    sim_closed_loop_robustness.m     % docking × eccentricity
    sim_closing_robustness.m         % closing × eccentricity
    sim_closing_robustness_drag.m    % closing × BC-ratio
    sim_nmc_robustness_empc.m        % faithful economic NMC × eccentricity
    sim_nmc_robustness_empc_drag.m   % faithful economic NMC × BC-ratio
```

Scripts under `simulations/*/` and moved tests reference project folders via `fullfile(script_dir, '..', '..', ...)` (two levels up from the subfolder).

Output: MATLAB figures, JSON exports to `exports/scenarios/`, TikZ exports to `results/figures/`.

## Architecture

```
src/
  dynamics/     % Orbital dynamics models (relative + absolute/ECI)
  control/      % Controller implementations + closed-loop run harnesses
  utils/        % Constants, discretization, integrators, frame transforms, references
simulations/    % Runnable scripts, split into foundations/ empc/ adrios/ model_validity/
tests/          % verification/ (correctness certificates) + analysis/ (turnpike studies)
visualization/  % Plotting and 3D viz utilities
exports/        % JSON scenario exports (results/ holds LaTeX/TikZ figures)
results/        % LaTeX figure exports (results/figures/, results/figures/rpo/)
tools/          % matlab2tikz for LaTeX figure generation
```

### Dynamics (`src/dynamics/`)

| File | Description |
|---|---|
| `clohessy_wiltshire.m` | CWH linearized relative equations (circular orbits) |
| `relative_motion_general.m` | Linearized relative motion for elliptic orbits (Schaub 14.20, LTV) |
| `relative_motion_nonlinear.m` | Exact nonlinear relative motion (Schaub 14.13) |
| `two_body.m` | Keplerian 2-body ODE for absolute orbit propagation |
| `truth_propagator_at.m` | High-fidelity single-craft ECI RHS (Aerospace Toolbox): two-body + optional J2–J4 zonal gravity + optional NRLMSISE-00 drag. Collapses exactly to `two_body.m` with perturbations off. |

### Control (`src/control/`)

| File | Description |
|---|---|
| `lqr_controller.m` | Discrete LQR via `dlqr` |
| `mpc_regulation.m` | Finite-horizon QP-based MPC (regulation to origin) |
| `mpc_tracking.m` | MPC with reference trajectory tracking |
| `constraints.m` | Builds QP constraint matrices (input bounds, state constraints) shared across MPC variants |
| `run_docking.m` / `run_closing.m` / `run_flyaround.m` | Single closed-loop runs on the CW plant (docking to port, far-range closing to a handover ball, economic NMC fly-around) |
| `run_docking_truth.m` / `run_closing_truth.m` / `run_flyaround_truth.m` | Same controllers, but the plant is the high-fidelity ECI truth propagator (`ode113`). The controller keeps its CW internal model; per step the CW-optimal Hill input is rotated to ECI and applied to the chaser. Used by the model-validity closed-loop sweeps. |

### Key Conventions

- **State vector**: `[x, y, z, vx, vy, vz]` in Hill/LVLH frame — x radial (R-bar), y along-track (V-bar), z cross-track
- **Units**: SI throughout (meters, kg, seconds); angles in radians
- **Orbit baseline**: ISS-like circular orbit at 400 km altitude (`src/utils/constants.m`)
- **Discretization**: Zero-order hold (ZOH) via `src/utils/discretize.m`
- **Integration**: RK4 via `src/utils/rk4_integrator.m` for the project models; the truth plant uses `ode113` at tight tolerances (deliberately a different integrator so model and truth do not share integration error)
- **Frame transforms**: `src/utils/eci2hill.m` / `hill2eci.m` (target-centred Hill/LVLH ↔ ECI), `src/utils/oe2eci.m` (classical orbital elements → ECI state, Vallado COE2RV), `src/utils/cw_stm.m` (closed-form CW state-transition matrix)
- **Reference trajectories**: V-bar, R-bar, R-bar+V-bar approach profiles from `src/utils/reference_trajectory.m`
- **NMC orbits**: closed-form Natural Motion Circumnavigation (2:1 ellipse) via `src/utils/nmc_trajectory.m`
- **Target kinematics**: ECI-state-based orbital kinematics (true anomaly rates) for general relative motion via `src/utils/target_kinematics.m`

### MPC / eMPC Formulation

All MPC variants solve a finite-horizon QP via `quadprog` with constraint matrices built by `constraints.m`. The regulation variant penalizes deviation from the origin; the tracking variant penalizes deviation from a pre-computed reference trajectory. Constraints (thrust limits, state bounds, terminal ball/equality) are enforced as inequality/equality constraints in the QP.

eMPC uses a fuel-minimizing stage cost `ℓ(u) = u'R_eco u` (no state tracking term) with a terminal cost `P` from a DARE with stabilizing weights (Amrit, Rawlings & Angeli 2011). The NMC fly-around (Scenario 2) selects the orbit with a **periodic terminal equality** `x(N) = Π*(k mod P)` rather than a tracking cost. Two eMPC scenarios: regulation to origin (`sim_empc_regulation.m`) and NMC acquisition/maintenance (`sim_empc_nmc.m`).

### Model-Validity Envelope Study

Purpose: quantify where the CW controller model loses fidelity, to bound which operations the thesis can defensibly tackle with CW + MPC. Method: propagate target and chaser absolutely in ECI (`truth_propagator_at.m`) and difference into LVLH (`eci2hill.m`) — the GMAT/STK "differential" truth standard. Open-loop layers compare CW / general / nonlinear models against the nonlinear (or differenced) reference; closed-loop sweeps run the real controllers against the truth plant and report Δv inflation, handover success, and orbit integrity. Headline finding: feedback guarantees stability, not efficiency — docking tolerates eccentricity cheaply (~5% Δv at e=0.3), closing pays heavily (~210%), and the economic NMC hits a hard QP-infeasibility cliff near e≈0.10.

### Tests (`tests/`)

| File | Description |
|---|---|
| `verification/verify_terminal_cost.m` | Checks the terminal cost matrix P satisfies the DARE and stability conditions |
| `verification/validate_truth_vs_nonlinear.m` | Regression anchor: truth propagator (perturbations off) vs exact nonlinear relative motion agree to sub-µm over 3 orbits |
| `verification/verify_dissipativity_lmi.m`, `verify_periodic_dissipativity_lmi.m` | Strict-dissipativity LMI certificates (eMPC/NMC stability, technical-report proofs) |
| `analysis/quantify_turnpike_tradeoff.m`, `plot_turnpike.m`, `synthesize_min_regularizer.m` | Turnpike behaviour and minimal-regularizer studies for the economic MPC |

### JSON Export Schema

Simulation scripts export results to `exports/scenarios/` as JSON with: metadata (date, sample time, horizon), performance metrics (peak thrust, per-axis delta-v, settling time, eigenvalues), and full state/input trajectories.
