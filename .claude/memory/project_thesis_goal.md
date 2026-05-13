---
name: project-thesis-goal
description: Thesis topic and code purpose — mission-centric requirements analysis for MPC in spacecraft RPO
metadata:
  type: project
---

Thesis title: "Mission-Centric Requirements Analysis for Model Predictive Control in Spacecraft Rendezvous and Proximity Operations"

The code has two roles:
1. **Algorithm implementation** — LQR/LQI done, standard MPC (regulation + tracking) done, **eMPC is the primary next implementation** and the main thesis use case.
2. **Validation tool** — simulation outputs (JSON exports, figures) validate the requirements framework developed in the thesis.

**Why:** The thesis argues that MPC controller design for spacecraft RPO should be driven by mission requirements (approach corridors, delta-v budgets, safety constraints). eMPC explicitly optimizes an economic objective (fuel, time) rather than a quadratic surrogate, making it the most mission-aligned formulation.

**How to apply:** When suggesting new features or controller extensions, prioritize what serves eMPC implementation or requirements validation. Keep the JSON export schema extensible for new metrics.
