# Docking Control via Economic Model Predictive Control

Research software for spacecraft rendezvous and docking using Economic Model Predictive Control (EMPC).

---

## Overview

This repository contains software developed for the simulation and evaluation of spacecraft rendezvous and proximity operations using Economic Model Predictive Control.

Current components include:

- spacecraft rendezvous dynamics
- Economic MPC implementations
- simulation environment
- numerical experiments
- supporting documentation

The repository accompanies ongoing academic research in mission-centric spacecraft guidance and control.

---

## Associated Publication

The accompanying technical report is available in `report/`.

The citation information can be found in `CITATION.cff`.

---

## Citation

If this software contributes to your research, please cite the associated publication.

Citation metadata is provided in `CITATION.cff`.

---

## Collaboration

This repository was developed as part of a research collaboration between

- Kyoto University
- University of Stuttgart

---

## Repository Structure

```
src/            Core library: dynamics models, controllers, utilities
simulations/    Runnable scenario scripts (foundations/, empc/, adrios/, model_validity/)
tests/          Correctness certificates (verification/) and turnpike studies (analysis/)
visualization/  Plotting and 3D visualization utilities
exports/        JSON scenario exports
results/        LaTeX/TikZ figure exports
tools/          matlab2tikz (vendored) for figure generation
```

See `CLAUDE.md` for a detailed description of each module.

---

## Where to find what (thesis figure/table index)

Not every script here made it into the final thesis — this repo also holds
exploratory work, superseded formulations, and prototypes that informed the
thesis without being cited directly. The table below maps each thesis
figure/table to the script that generates it, for anyone trying to reproduce
a specific result.

| Thesis item | Script | Output |
|---|---|---|
| Fig. 2.4, 2.5 (NMC vs. free drift) | `simulations/empc/sim_nmc_free_drift.m` | `results/figures/nmc_free_drift/` |
| Fig. 3.1 (LQR baseline) | `simulations/foundations/sim_lqr.m` | `results/figures/lqr/` |
| Fig. 5.1, Table 5.2 (Op. 1 — far-range closing) | `simulations/adrios/sim_scenario4_closing.m` | `results/figures/scenario4/B_hill_frame.tikz`, `exports/scenarios/sim_scenario4_closing.json` |
| Fig. 5.2, Tables 5.3–5.4 (Op. 2 — NMC fly-around) | `simulations/adrios/sim_scenario2_flyaround.m` | `results/figures/s2_flyaround/`, `exports/scenarios/sim_scenario2_flyaround.json` |
| Fig. 5.3, Table 5.5 (Op. 3 — final docking) | `simulations/adrios/sim_scenario1_docking.m` | `results/figures/s1_docking/`, `exports/scenarios/sim_scenario1_docking.json` |
| Table 5.6 (end-to-end mission cost) | `simulations/adrios/sim_adrios_combinations.m` | `results/figures/s_full_mission/combinations_table.tex` |
| Table 6.1 ($M_1$ eigenvalues) | `tests/verification/verify_terminal_cost.m` | console output |
| Fig. 7.1 (Layer 1a envelope) | `simulations/model_validity/sim_model_validity_sweep.m` | `results/figures/validity_sweep/model_error_e_sep.tikz` |
| Fig. 7.2–7.3 (Layer 1b envelope) | `simulations/model_validity/sim_perturbation_sweep.m` | `results/figures/validity_sweep/model_error_alt_bc_{heatmap,corner}.tikz` |
| Fig. 7.4, Table 7.1 row (Op. 3 robustness) | `simulations/model_validity/sim_closed_loop_robustness.m` | `results/figures/validity_sweep/closed_loop_robustness.tikz` |
| Fig. 7.5 (Op. 3 Hill-frame vs. $e$) | `simulations/model_validity/sim_hillframe_docking.m` | `results/figures/validity_sweep/hillframe_docking.tikz` |
| Fig. 7.6 (Op. 1 robustness) | `simulations/model_validity/sim_closing_robustness.m` | `results/figures/validity_sweep/closing_robustness.tikz` |
| Fig. 7.7 (Op. 1 Hill-frame vs. $e$) | `simulations/model_validity/sim_hillframe_closing.m` | `results/figures/validity_sweep/hillframe_closing.tikz` |
| Fig. 7.8 (Op. 2 robustness) | `simulations/model_validity/sim_nmc_robustness_empc.m` | `results/figures/validity_sweep/nmc_robustness_empc.tikz` |
| Fig. 7.9 (Op. 2 Hill-frame vs. $e$) | `simulations/model_validity/sim_hillframe_nmc.m` | `results/figures/validity_sweep/hillframe_nmc.tikz` |

Scripts not listed above (e.g. `simulations/foundations/sim_lqi.m`,
`sim_rpo.m`, `sim_two_body.m`; `simulations/empc/sim_empc_regulation.m`,
`sim_empc_nmc.m`; `tests/analysis/*`; and the extra
`simulations/model_validity/sim_{closing_regularizer_comparison,closing_robustness_drag,nmc_robustness_empc_drag,nmc_terminal_comparison}.m`
sweeps) are exploratory or superseded work that is not cited by figure/table
number in the thesis, kept here for context and reproducibility rather than
deletion.

---

## Notes for Contributors

- Follow existing structure and conventions  
- Keep changes modular and documented  
- Coordinate major changes within the collaboration  

---

## License

This project is released under the BSD 3-Clause License.

See `LICENSE.txt`.

---

## Disclaimer

This software is intended for research and educational purposes.

It is experimental software and is provided without any guarantee of correctness, completeness, or suitability for safety-critical applications.
