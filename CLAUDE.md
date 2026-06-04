# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Master's thesis (University of Stuttgart, co-supervised by Kyoto University) on **mission-centric requirements analysis for Model Predictive Control in Spacecraft Rendezvous and Proximity Operations**. The code serves two purposes:

1. **Algorithm implementation and comparison** — LQR/LQI as baselines, standard MPC (regulation + tracking), and eMPC (regulation + NMC orbit acquisition/maintenance) are all implemented. eMPC is the primary focus and main use case of the thesis.
2. **Performance analysis and validation** — simulation outputs (JSON exports, figures) are used to validate the requirements framework developed in the thesis.

The dynamics are modeled using Clohessy-Wiltshire-Hill (CWH) equations and more general relative motion formulations.

## Getting Started

```matlab
% Once per MATLAB session — loads all project paths
setup.m
```

No compile or build step. Requires MATLAB with **Control System Toolbox** (`dlqr`, `c2d`) and **Optimization Toolbox** (`quadprog`). Figure export requires `tools/matlab2tikz-master/`.

## Running Simulations

All simulations are self-contained scripts in `simulations/`:

```matlab
sim_lqr.m               % LQR baseline
sim_lqi.m               % LQI (LQR + integral action on position)
sim_mpc_regulation.m    % MPC driving state to origin
sim_mpc_tracking.m      % MPC tracking a reference trajectory
sim_empc_regulation.m   % eMPC regulation — min-fuel, compared against standard MPC
sim_empc_nmc.m          % eMPC NMC orbit acquisition and maintenance
sim_rpo.m               % Relative periodic orbit analysis
sim_two_body.m          % 2-body problem validation
```

Output: MATLAB figures, JSON exports to `exports/scenarios/`, TikZ exports to `results/figures/`.

## Architecture

```
src/
  dynamics/     % Orbital dynamics models
  control/      % Controller implementations
  utils/        % Constants, discretization, integrators, reference trajectories
simulations/    % Top-level runnable scripts
tests/          % Validation and verification scripts
visualization/  % Plotting and 3D viz utilities
exports/        % JSON scenario exports and LaTeX figure exports
tools/          % matlab2tikz for LaTeX figure generation
```

### Dynamics (`src/dynamics/`)

| File | Description |
|---|---|
| `clohessy_wiltshire.m` | CWH linearized equations (circular orbits) |
| `relative_motion_general.m` | Linearized relative motion for elliptic orbits (Schaub 14.20) |
| `relative_motion_nonlinear.m` | Exact nonlinear relative motion (Schaub 14.13) |
| `two_body.m` | Keplerian 2-body ODE for absolute orbit propagation |

### Control (`src/control/`)

| File | Description |
|---|---|
| `lqr_controller.m` | Discrete LQR via `dlqr` |
| `mpc_regulation.m` | Finite-horizon QP-based MPC (regulation to origin) |
| `mpc_tracking.m` | MPC with reference trajectory tracking |
| `constraints.m` | Builds QP constraint matrices (input bounds, state constraints) shared across MPC variants |

### Key Conventions

- **State vector**: `[x, y, z, vx, vy, vz]` in Hill/LVLH frame — x radial, y along-track, z cross-track
- **Units**: SI throughout (meters, kg, seconds)
- **Orbit baseline**: ISS-like circular orbit at 400 km altitude (`src/utils/constants.m`)
- **Discretization**: Zero-order hold (ZOH) via `src/utils/discretize.m`
- **Integration**: RK4 via `src/utils/rk4_integrator.m`
- **Reference trajectories**: V-bar, R-bar, R-bar+V-bar approach profiles from `src/utils/reference_trajectory.m`
- **NMC orbits**: closed-form Natural Motion Circumnavigation (2:1 ellipse) via `src/utils/nmc_trajectory.m`
- **Target kinematics**: ECI-state-based orbital kinematics (true anomaly rates) for general relative motion via `src/utils/target_kinematics.m`

### MPC / eMPC Formulation

All MPC variants solve a finite-horizon QP via `quadprog` with constraint matrices built by `constraints.m`. The regulation variant penalizes deviation from the origin; the tracking variant penalizes deviation from a pre-computed reference trajectory. Constraints (thrust limits) are enforced via inequality constraints in the QP.

eMPC uses a fuel-minimizing stage cost `ℓ(u) = u'R_eco u` (no state tracking term) with a terminal cost `P` computed from a DARE with stabilizing weights (Amrit, Rawlings & Angeli 2011). Two scenarios are implemented: regulation to origin (`sim_empc_regulation.m`) and NMC orbit acquisition/maintenance (`sim_empc_nmc.m`).

### Tests (`tests/`)

| File | Description |
|---|---|
| `test_two_body.m` | Validates the two-body propagator |
| `verify_terminal_cost.m` | Checks that the terminal cost matrix P satisfies the DARE and stability conditions |

### JSON Export Schema

Simulation scripts export results to `exports/scenarios/` as JSON with: metadata (date, sample time, horizon), performance metrics (peak thrust, per-axis delta-v, settling time, eigenvalues), and full state/input trajectories.
