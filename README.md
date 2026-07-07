# TERGM.jl


[![Network Analysis](https://img.shields.io/badge/Network-Analysis-orange.svg)](https://github.com/statistical-network-analysis-with-Julia/TERGM.jl)
[![Build Status](https://github.com/statistical-network-analysis-with-Julia/TERGM.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/statistical-network-analysis-with-Julia/TERGM.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://statistical-network-analysis-with-Julia.github.io/TERGM.jl/stable/)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://statistical-network-analysis-with-Julia.github.io/TERGM.jl/dev/)
[![Julia](https://img.shields.io/badge/Julia-1.12+-purple.svg)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<p align="center">
  <img src="docs/src/assets/logo.svg" alt="TERGM.jl icon" width="160">
</p>

Separable Temporal ERGMs (STERGM) in Julia — a port of the R `tergm`
package (Krivitsky & Handcock 2014).

## The separable model

Each transition Y_{t−1} → Y_t factors into:

- a **formation** model on the formation network **Y⁺ = Y_{t−1} ∪ Y_t**
  (free dyads: the non-edges of Y_{t−1});
- a **dissolution** model on the dissolution network **Y⁻ = Y_{t−1} ∩ Y_t**
  (free dyads: the edges of Y_{t−1}), parameterized as **persistence** —
  positive coefficients mean ties last longer.

`formation_network` / `dissolution_network` expose the auxiliary
construction directly. Formation and dissolution models take any mix of
standard ERGM.jl terms (evaluated on Y⁺/Y⁻) and temporal terms
(`EdgeStability`, `Delrecip`, `PersistentEdge`, `NewEdge`) that also
condition on Y_{t−1}.

## Quick Start

```julia
using TERGM, ERGM, Network

networks = [net_t0, net_t1, net_t2]   # panel of same-sized Networks

result = stergm(networks,
                [Edges(), Mutual()],       # formation model
                [Edges()])                 # dissolution (persistence) model
println(result)

# Simulate forward from the last panel
future = simulate_stergm(result, 10)

# Transition-level goodness of fit
stergm_gof(result; n_sim = 100)
```

## Estimation

The default (and honest) estimator is **CMPLE**: pooled logistic
pseudo-likelihood over the free dyads of the auxiliary networks. For
dyad-independent terms this *is* the conditional MLE; for dyad-dependent
terms (`Triangle`, `Mutual`, ...) it is an approximation — `method =
:cmle` warns and falls back to CMPLE (MCMC-based CMLE is not
implemented), and `method = :egmme` raises an error rather than
fabricating estimates.

An edges-only model reproduces the analytic formation/persistence
log-odds exactly, and simulation→estimation round trips recover both
coefficient vectors (tested).

## Temporal descriptives

`edge_ages(networks)` / `mean_edge_age(networks)` report how long the
final panel's ties have been in place.

## References

1. Krivitsky, P.N. & Handcock, M.S. (2014). A separable model for dynamic
   networks. *JRSS-B*, 76(1), 29-46.

2. Krivitsky, P.N. & Handcock, M.S. tergm: Fit, Simulate and Diagnose
   Models for Network Evolution Based on Exponential-Family Random Graph
   Models. R package.
   [https://cran.r-project.org/package=tergm](https://cran.r-project.org/package=tergm)

## License

MIT License - see [LICENSE](LICENSE) for details.
