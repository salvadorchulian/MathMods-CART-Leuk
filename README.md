# Mathematical Models for CAR-T Immunotherapy and CD19 Dynamics in Leukemia: A Comparative Analysis

This repository contains the **MATLAB** implementation of the mathematical frameworks, interactive interfaces, and sensitivity analysis tools developed for the research paper:
> **Mathematical Models for CAR-T immunotherapy and CD19 Dynamics in Leukemia: a comparative analysis** > *Salvador Chulián, Ana Niño-López, Rocío Picón-González, María Rosa* > Department of Mathematics, University of Cádiz & Biomedical Research and Innovation Institute of Cádiz (INiBICA), Spain.

---

## Project Overview

This project provides code for two distinct mathematical modeling frameworks—Ordinary Differential Equations (ODE) and Partial Differential Equations (PDE)—designed to analyze the long-term interaction kinetics between CAR-T cells, healthy B cells, and heterogeneous leukemic cell clones in B-cell Acute Lymphoblastic Leukemia (B-ALL). 

Unlike classical mathematical representations that treat antigen loss as an irreversible absorbing state, these frameworks integrate immune-regulated, bidirectional phenotypic transitions governed by Michaelis-Menten dynamics. This enables a detailed, mechanistic evaluation of tumor escape selection (antigen escape) and tumor resensitization under targeted immunotherapy pressure.

---

## Modeling Frameworks

**Compartmental Ordinary Differential Equation (ODE) Model** A time-dependent compartmental system mapping discrete, binary expression states of the CD19 antigen on tumor surfaces alongside functional CAR-T subpopulations:
* `LN(t)`: CD19-negative leukemic cells that evade direct CAR-T cytotoxicity.
* `LP(t)`: CD19-positive leukemic cells, the primary target of the immune response.
* `B(t)`: Healthy mature B cells, which constitutively express CD19 and face collateral depletion.
* `CA(t)` and `CM(t)`: Activated and memory CAR-T cell pools regulating expansion and durable immunological memory.

**Antigen-Structured Partial Differential Equation (PDE) Model** A continuous mathematical framework that models healthy and leukemic populations along a continuous spectrum of CD19 expression levels denoted by the spatial-like coordinate $x \in [0, 1]$. This formulation simulates fine-grained clonal evolution, continuous antigen downregulation, and spatial-like drift across the phenotype space induced by immune selection pressures.

---

## Repository Structure

The root directory contains the following verified MATLAB scripts and documentation:

**FOLDER Scenario_Simulations**
* **`interactive_ODE_model.m`**: A full-featured Graphical User Interface (GUI) containing interactive sliders (*UI sliders*). It allows users to manipulate biological parameters in real time and instantly visualize the resulting Lotka-Volterra type predator-prey oscillations.
* **`interactive_PDE_model.m`**: An interactive visual dashboard tailored for the antigen-structured PDE model. It simulates and animates density profiles across the $x$ antigen domain dynamically as parameters are modified.
* **`Simulations_PDE.pdf`**: A supplementary data file illustrating temporal trajectories, density distributions, and simulation panels across different initial clone configurations (e.g., single clone centered at 0.4 or 0.6, and mixed multi-clonal setups).

**FOLDER Antigen_Proportion_Analyses**
* **`CAR_T_ODE_Antigen_Proportion.m`**: Simulates the compartmental ODE model over a multi-year horizon (up to 10 years). It performs bivariate parameter sweeps across critical thresholds ($h$ and $k$) to construct multi-panel heatmaps of long-term tumor composition across different initial conditions.
* **`CAR_T_PDE_Antigen_Proportion.m`**: Executes parameter sweeps and spatial integration over the continuous antigen-structured PDE domain to evaluate long-term clonal composition, mean antigen profiles, and clearance rates.
* **`ODE_Ant_***.png files`**: Simulations from the ODE model using different scenarios: 100percent Antigen Negative, 100percent Antigen Positive and 50% and 50%  antigen positive and negative (mixed clones).
* **`PDE_Ant_***.png files`**: Simulations from the PDE model using different scenarios: a single clone centered at x=0.4, a single clone centered at x=0.6, and two clones centered at x=0.4 and x=0.6.

**FOLDER Sensitivity_Analyses**
* **`SobolODEfinal.m`**: Conducts a Global Sensitivity Analysis (GSA) using the Saltelli sampling scheme to compute first-order and total Sobol indices, identifying which biochemical parameters drive variance in the ODE populations.
* **`SobolPDEfinal.m`**: Implements the sensitivity analysis framework optimized for the continuous PDE model to track parameter importance.
* **`Sensitivity_ODE.pdf`** Heatmap summary of the sensitivity analyses from the ODE model.
* 
* **`Sensitivity_PDE.pdf`** Heatmap summary of the sensitivity analyses from the PDE model.

---

## Numerical Implementation & Safeguards

The computational engine integrates specialized design strategies to preserve structural and biological stability:
* **Adaptive Integration & Stiffness:** Time marching leverages adaptive solvers (`ode45` and `ode15s`) configured with tight absolute and relative error tolerances to stably resolve the stiff, coupled tracking of rapid immune expansion and slow phenotypic shifts.
* **Spatial Discretization:** The PDE continuum is discretized uniformly over $N_x = 100$ spatial mesh nodes. Global antigen burdens are calculated at each timestep using numerical trapezoidal quadrature.
* **Biological Truncation Threshold:** To eliminate non-physical numerical artifacts (such as sub-cellular fractional values) and correctly model the clinical biology of absolute disease eradication, any cellular population falling below **one cell unit** ($< 1$) during integration loops is automatically truncated to zero.
* **Minimal Residual Disease (MRD):** Evaluation scripts calculate the clinical manifestation of MRD, using a threshold of $0.1\%$ total leukemic fraction for clear graphic visibility in simulations.

---

## Citation

If you use these scripts, mathematical models, or simulation setups in an academic or industrial publication, please cite the core paper:

```bibtex
@article{chulian2026mathematical,
  title={Mathematical Models for CAR-T immunotherapy and CD19 Dynamics in Leukemia: a comparative analysis},
  author={Chuli{\'a}n, Salvador and Ni{\~n}o-L{\'o}pez, Ana and Pic{\'o}n-Gonz{\'a}lez, Roc{\'\i}o and Rosa, Mar{\'\i}a},
  journal={Arxiv},
  year={2026}
}
