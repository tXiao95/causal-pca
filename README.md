# Causal Sufficient Dimension Reduction for Multiple Continuous Exposures

Environmental health studies increasingly seek to estimate the joint effects of multiple continuous exposures, such as air pollution components, PFAS, and chemical mixtures. This problem has been identified as a major priority area by the National Institute of Environmental Health Sciences (NIEHS). However, existing approaches often enforce tradeoffs between flexibility, interpretability, computational feasibility, and causal interpretability. This repository contains code to reproduce the analyses from our study on causal sufficient dimension reduction (causal SDR), a framework for estimating and interpreting the causal effects of multiple continuous exposures through low-dimensional representations that preserve the causal exposure-response surface.

## Motivation

Many existing methods for environmental mixtures possess important limitations:

| Method | Flexible Estimation | Interpretable | Causally Valid | Scalable |
|---|---|---|---|---|
| BKMR | ✓ | ✗ | ✗ | ✗ |
| WQS | ✗ | ✓ | ✗ | ✓ |
| Quantile g-computation | ✗ | ✓ | ✓ | ✓ |
| Principal components pursuit | Partial | Partial | ✗ | ✓ |
| Proposed causal SDR | ✓ | ✓ | ✓ | ✓ |

To address these limitations, we develop **causal sufficient dimension reduction (causal SDR)**, a framework for estimating low-dimensional representations of multiple continuous exposures that preserve the causal exposure-response surface (ERS). The output is essentially a causal PCA, where the principal components instead of chasing variance respect the underlying causal structure. 

Building on classical sufficient dimension reduction (SDR) methodology, we define a **causal central mean subspace (CCMS)** as the smallest subspace that preserves the causal mean exposure-response function under standard causal identification assumptions.

## Background

Previous work by [Nabi et al. (2022)](https://proceedings.mlr.press/v180/nabi22a.html) proposed a marginal structural modeling (MSM) framework for causal SDR based on locally efficient score equations. While theoretically appealing, this approach requires iterative numerical optimization and substantial nuisance function estimation, including nuisance components that depend on the dimension reduction parameter itself.

To circumvent these computational challenges, we propose a two-stage pseudo-outcome implementation for causal SDR. This strategy:

- substantially simplifies implementation,
- reduces nuisance estimation burden,
- leverages existing tools from causal inference for continuous exposures and SDR, and
- enables flexible downstream estimation of the causal ERS using the reduced exposure representation.

## Main Contributions

Our theoretical and simulation results establish that:

- causal SDR consistently recovers the causal central mean subspace,
- the estimated structural dimension remains consistent even when pseudo-outcomes are estimated,
- the convergence rate depends on the quality of pseudo-outcome estimation, and
- downstream inference based on the reduced exposure is asymptotically equivalent to the oracle setting where the reduction is known.

![](results/jasa-initial-submission/main_paper_final_results_nnet/interaction/Simulation_Coverage.png)


## Real Data Application

We apply the proposed causal SDR framework to a four-component PFAS mixture to study the causal effect of maternal PFAS exposure on infant birthweight in a cohort of 305 African-American mothers from the Atlanta African-American Maternal-Child Cohort study. The analysis demonstrates how causal SDR can improve both interpretability of exposure mixtures and statistical efficiency in environmental mixture analyses. We additionally compare results with those of BKMR and qgcomp. 

![](results/jasa-initial-submission/ATL-AA/ERS.png)

---

# Repository Structure

| Directory | Description |
|---|---|
| `R/` | Functions implementing causal SDR and benchmark methods |
| `code/` | Analysis scripts for simulation studies and ATL-AA real data analysis |
| `manuscript/` | Posters, slides, and dissemination materials |
| `results/` | Main simulation and real data analysis outputs |

### Notes on Results

Within `results/`:

- `final_results_nnet/` contains the primary results reported in the manuscript.
- `main_paper_final_results_SL_xgboost/` and `main_paper_final_results_nonsparse_nnet_128x64x32/` correspond to supplementary analyses reported in the supplement.

---

# References

- Nabi, R., McNutt, T., & Shpitser, I. (2022). *Semiparametric causal sufficient dimension reduction of multidimensional treatments*. Proceedings of Machine Learning Research. https://proceedings.mlr.press/v180/nabi22a.html
