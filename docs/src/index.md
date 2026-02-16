# TERGM.jl

*Temporal Exponential Random Graph Models for Julia*

A Julia package for statistical modeling of dynamic networks using Separable Temporal ERGMs (STERGMs), with support for temporal terms, conditional maximum likelihood estimation, and network simulation.

## Overview

Temporal Exponential Random Graph Models (TERGMs) extend the ERGM framework to dynamic networks observed at multiple time points. TERGM.jl focuses on the Separable Temporal ERGM (STERGM), which decomposes network change into two separate processes: **edge formation** and **edge dissolution** (persistence).

TERGM.jl is a port of the R [tergm](https://cran.r-project.org/package=tergm) package from the [StatNet](https://statnet.org/) collection.

### What is a STERGM?

A Separable Temporal ERGM models the transition from network $Y_{t-1}$ to $Y_t$ as:

$$P(Y_t \mid Y_{t-1}) \propto P^+(Y_t \mid Y_{t-1}; \theta^+) \times P^-(Y_t \mid Y_{t-1}; \theta^-)$$

Where:

- $P^+$ governs **edge formation** (for non-edges at $t-1$): which new ties form?
- $P^-$ governs **edge persistence** (for edges at $t-1$): which existing ties persist?
- $\theta^+$ and $\theta^-$ are separate parameter vectors

This separability assumption makes estimation tractable and provides clear interpretations: formation parameters describe why ties form, and dissolution parameters describe why ties persist.

### Key Concepts

| Concept | Description |
|---------|-------------|
| **STERGM** | Separable Temporal ERGM with formation and dissolution models |
| **Formation** | Process governing new edge creation |
| **Dissolution** | Process governing edge persistence (vs. dissolution) |
| **Temporal Terms** | ERGM terms that reference the previous network state |
| **CMLE** | Conditional Maximum Likelihood Estimation |
| **CMPLE** | Conditional Maximum Pseudo-Likelihood Estimation |

### Applications

TERGMs are widely used in:

- **Social network evolution**: How do friendship networks change over time?
- **Organizational dynamics**: Modeling tie formation and dissolution in firms
- **Policy networks**: How do collaborations among organizations evolve?
- **Epidemiology**: Understanding contact network dynamics for disease modeling
- **International relations**: Modeling the formation and dissolution of alliances

## Features

- **STERGM framework**: Separate formation and dissolution models with distinct parameters
- **Temporal terms**: Edge stability, memory, delayed reciprocity, edge age, time-lagged terms
- **Multiple estimation methods**: CMLE, CMPLE, and experimental EGMME
- **Network simulation**: Simulate future network sequences from fitted models
- **Goodness-of-fit**: Diagnostic tools for model assessment
- **Standard ERGM terms**: Use any ERGM.jl term in formation or dissolution models

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Statistical-network-analysis-with-Julia/TERGM.jl")
```

Or for development:

```julia
using Pkg
Pkg.develop(path="/path/to/TERGM.jl")
```

## Quick Start

```julia
using Network
using ERGM
using TERGM

# Network sequence (panel data at 4 time points)
networks = [net_t1, net_t2, net_t3, net_t4]

# Define formation and dissolution terms
form_terms = [Edges(), Triangle()]
diss_terms = [Edges()]

# Fit STERGM via Conditional Maximum Likelihood
result = stergm(networks, form_terms, diss_terms)
println(result)

# Simulate 10 future time steps
future = simulate_stergm(result, 10)
```

## Choosing Terms

| Use Case | Formation Terms | Dissolution Terms |
|----------|----------------|-------------------|
| Basic density | [`Edges()`] | [`Edges()`] |
| Reciprocity effects | [`Edges()`, `Mutual()`] | [`Edges()`, `Mutual()`] |
| Clustering | [`Edges()`, `Triangle()`] | [`Edges()`] |
| Delayed reciprocity | [`Edges()`, `Delrecip()`] | [`Edges()`] |
| Edge persistence | [`Edges()`] | [`Edges()`, `EdgeStability()`] |

## Documentation

```@contents
Pages = [
    "getting_started.md",
    "guide/terms.md",
    "guide/estimation.md",
    "guide/simulation.md",
    "api/types.md",
    "api/terms.md",
    "api/estimation.md",
]
Depth = 2
```

## Theoretical Background

### The Separability Assumption

The key insight of STERGMs is that network transitions can be decomposed into two independent processes operating on disjoint sets of dyads:

1. **Formation**: Among dyads where no edge exists at $t-1$, which edges form at $t$?
2. **Persistence**: Among dyads where an edge exists at $t-1$, which edges persist to $t$?

This decomposition is exact when formation and dissolution probabilities depend only on the current network state and the edge's own status, not on what other edges are simultaneously forming or dissolving.

### Conditional Maximum Likelihood

CMLE conditions on the previous network state $Y_{t-1}$ and estimates parameters by maximizing the conditional likelihood of observing $Y_t$. For each transition:

- Formation parameters are estimated from the set of dyads that were non-edges at $t-1$
- Dissolution parameters are estimated from the set of dyads that were edges at $t-1$

Each sub-problem reduces to a logistic regression on change statistics.

## References

1. Krivitsky, P. N., & Handcock, M. S. (2014). A separable model for dynamic networks. *Journal of the Royal Statistical Society: Series B*, 76(1), 29-46.

2. Krivitsky, P. N., & Handcock, M. S. (2023). Modeling of dynamic networks based on egocentric data with durational information. *Sociological Methodology*, 53(2), 250-286.

3. Robins, G., & Pattison, P. (2001). Random graph models for temporal processes in social networks. *Journal of Mathematical Sociology*, 25(1), 5-41.

4. Hanneke, S., Fu, W., & Xing, E. P. (2010). Discrete temporal models of social networks. *Electronic Journal of Statistics*, 4, 585-605.
