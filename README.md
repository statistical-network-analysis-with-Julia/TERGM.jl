# TERGM.jl

Temporal Exponential Random Graph Models for Julia.

## Overview

TERGM.jl provides models for dynamic networks, including Separable Temporal ERGMs (STERGMs) that model edge formation and dissolution as separate processes. It enables statistical inference on network dynamics from panel data.

This package is a Julia port of the R `tergm` package from the StatNet collection.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Statistical-network-analysis-with-Julia/TERGM.jl")
```

## Features

- **STERGM**: Separable Temporal ERGMs with formation/dissolution models
- **Temporal terms**: Edge stability, memory, delayed reciprocity
- **Estimation**: CMLE, CMPLE methods
- **Simulation**: Generate network sequences from fitted models

## Quick Start

```julia
using Network
using TERGM

# Network sequence (panel data)
networks = [net_t1, net_t2, net_t3, net_t4]

# Define formation and dissolution terms
form_terms = [Edges(), Triangle()]
diss_terms = [Edges()]

# Fit STERGM
result = stergm(networks, form_terms, diss_terms)
println(result)

# Simulate future networks
future = simulate_stergm(result, 10)
```

## STERGM Model

The Separable Temporal ERGM models network change as:

```
P(Y_t | Y_{t-1}) ∝ P+(Y_t | Y_{t-1}; θ+) × P-(Y_t | Y_{t-1}; θ-)
```

Where:
- `P+` governs edge formation (for non-edges at t-1)
- `P-` governs edge persistence (for edges at t-1)

```julia
# Create STERGM specification
formula = STERGM(formation_terms, dissolution_terms)

# Create model with data
model = STERGMModel(formula, networks)
```

## Temporal Terms

### Edge Dynamics
```julia
EdgeStability()      # Edges persisting from previous time
PersistentEdge()     # Same as EdgeStability
NewEdge()            # Edges not present at previous time
```

### Memory Effects
```julia
Memory(theta)        # Weighted memory of previous state
EdgeAge()            # Age of edges (requires sequence)
```

### Reciprocity
```julia
Delrecip()           # Delayed reciprocity (reciprocate previous edges)
```

### Lagged Terms
```julia
TimeLag(term, lag)   # Apply term to lagged network
```

### Wrapping Standard Terms
```julia
FormationTerm(Edges())      # Apply only to formation
DissolutionTerm(Edges())    # Apply only to dissolution
```

## Estimation

```julia
# Conditional Maximum Likelihood
result = stergm(networks, form_terms, diss_terms; method=:cmle)

# Conditional Maximum Pseudo-Likelihood
result = stergm(networks, form_terms, diss_terms; method=:cmple)

# Access coefficients
result.formation_coef    # Formation coefficients
result.dissolution_coef  # Dissolution (persistence) coefficients
result.formation_se      # Standard errors
result.dissolution_se
```

## Simulation

```julia
# Simulate from fitted model
future_nets = simulate_stergm(result, n_steps; burnin=100)

# Simulate from parameters
nets = simulate_network_sequence(
    formula, init_net, n_steps;
    form_coef=θ_plus,
    diss_coef=θ_minus
)
```

## Goodness-of-Fit

```julia
gof_result = stergm_gof(result; n_sim=100)
```

## Example: Email Communication

```julia
# Weekly email networks
weeks = [email_week1, email_week2, email_week3, email_week4]

# Model: edges form/dissolve based on density and reciprocity
form = [Edges(), Mutual()]
diss = [Edges(), Mutual()]

result = stergm(weeks, form, diss)

# Positive formation mutual coefficient → reciprocity increases formation
# Positive dissolution mutual coefficient → reciprocity increases persistence
```

## References

- Krivitsky, P. N., & Handcock, M. S. (2014). A separable model for dynamic networks. Journal of the Royal Statistical Society: Series B, 76(1), 29-46.

## License

MIT License
