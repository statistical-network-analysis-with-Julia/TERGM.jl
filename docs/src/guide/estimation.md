# STERGM Estimation

This guide covers the estimation methods available in TERGM.jl for fitting Separable Temporal ERGMs.

## Overview

STERGM estimation leverages the separability assumption to decompose the problem into two independent logistic regressions:

1. **Formation model**: Logistic regression on non-edges at $t-1$, predicting which edges form at $t$
2. **Dissolution model**: Logistic regression on edges at $t-1$, predicting which edges persist at $t$

Each sub-problem is solved by maximizing a conditional likelihood function using Newton-Raphson optimization.

## Estimation Methods

### Conditional Maximum Likelihood (CMLE)

CMLE is the default and recommended estimation method. It conditions on the previous network $Y_{t-1}$ and maximizes the conditional likelihood of observing the transition to $Y_t$.

```julia
result = stergm(networks, form_terms, diss_terms; method=:cmle)
```

**How it works:**

For each transition from $Y_{t-1}$ to $Y_t$:

1. Identify the **non-edge set** at $t-1$ (candidates for formation)
2. Identify the **edge set** at $t-1$ (candidates for dissolution)
3. For each dyad in the non-edge set, compute formation change statistics and whether the edge formed
4. For each dyad in the edge set, compute dissolution change statistics and whether the edge persisted
5. Fit logistic regression via Newton-Raphson on each subset

All transitions are pooled to obtain a single set of formation and dissolution coefficients.

**Parameters:**

```julia
result = stergm(networks, form_terms, diss_terms;
    method=:cmle,
    maxiter=100,    # Maximum Newton-Raphson iterations
    tol=1e-6        # Convergence tolerance for gradient norm
)
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `maxiter` | Maximum Newton-Raphson iterations | 100 |
| `tol` | Convergence tolerance | 1e-6 |

### Conditional Maximum Pseudo-Likelihood (CMPLE)

CMPLE is computationally equivalent to CMLE for STERGMs because the separability assumption makes the conditional likelihood decompose into independent dyad-level contributions.

```julia
result = stergm(networks, form_terms, diss_terms; method=:cmple)
```

In practice, CMPLE and CMLE produce identical results for STERGMs. CMPLE is provided for API compatibility with the R tergm package.

### Equilibrium Generalized Method of Moments (EGMME)

EGMME is an experimental method for estimating STERGMs when the network is assumed to be at equilibrium. It matches observed network statistics to their expected values under the equilibrium distribution.

```julia
result = stergm(networks, form_terms, diss_terms; method=:egmme)
```

**Note**: EGMME is experimental and may not converge reliably. Use CMLE for production analyses.

## The STERGMModel Type

The estimation pipeline centers on the `STERGMModel` type:

```julia
# Create model specification
formula = STERGM(formation_terms, dissolution_terms)

# Bind to data
model = STERGMModel(formula, networks)

# Or with custom time labels
model = STERGMModel(formula, networks; times=[1.0, 2.0, 3.0, 4.0])
```

The `STERGMModel` validates that:
- At least 2 networks are provided
- All networks have the same vertex count
- Time labels match network count

## Understanding Results

### The STERGMResult Type

```julia
result = stergm(networks, form_terms, diss_terms)

# Formation model
result.formation_coef    # Coefficient vector
result.formation_se      # Standard error vector

# Dissolution model
result.dissolution_coef  # Coefficient vector
result.dissolution_se    # Standard error vector

# Model metadata
result.model             # The STERGMModel object
result.method            # :cmle, :cmple, or :egmme
result.converged         # Bool
result.loglik            # Log-likelihood (NaN for CMLE)
```

### Printing Results

```julia
println(result)
```

Output:

```text
STERGM Results
==============
Method: cmle
Converged: true
Log-likelihood: NaN

Formation Model:
----------------------------------------
  edges                   -3.2456 (SE: 0.4523)
  mutual                   1.2345 (SE: 0.3891)

Dissolution (Persistence) Model:
----------------------------------------
  edges                    2.1234 (SE: 0.5012)
  mutual                   0.8765 (SE: 0.4123)
```

### Interpreting Formation Coefficients

Formation coefficients are **log-odds** for edge formation among non-edges at $t-1$:

$$\text{logit}[P(\text{form}_{ij})] = \theta^+_1 \cdot \Delta x_1(i,j) + \theta^+_2 \cdot \Delta x_2(i,j) + \ldots$$

| Term | Coefficient | Probability Change |
|------|-------------|-------------------|
| `edges` = -3.0 | Baseline formation prob = 0.047 | -- |
| `mutual` = 1.5 | If $(j,i)$ exists: prob increases | exp(1.5) = 4.5x odds |
| `triangle` = 0.8 | Per shared partner | exp(0.8) = 2.2x odds |

### Interpreting Dissolution Coefficients

Dissolution coefficients are **log-odds** for edge **persistence** among edges at $t-1$:

$$\text{logit}[P(\text{persist}_{ij})] = \theta^-_1 \cdot \Delta x_1(i,j) + \theta^-_2 \cdot \Delta x_2(i,j) + \ldots$$

**Important**: Positive coefficients increase **persistence** (decrease dissolution):

| Term | Coefficient | Interpretation |
|------|-------------|----------------|
| `edges` = 2.0 | Baseline persistence = 0.88 | 88% of edges persist |
| `mutual` = 1.0 | Mutual edges persist more | exp(1) = 2.7x odds |
| `edges` = -1.0 | Baseline persistence = 0.27 | Only 27% persist (high turnover) |

### Computing Confidence Intervals

```julia
using Distributions

z = quantile(Normal(), 0.975)  # 1.96 for 95% CI

# Formation CIs
form_lower = result.formation_coef .- z .* result.formation_se
form_upper = result.formation_coef .+ z .* result.formation_se

# Dissolution CIs
diss_lower = result.dissolution_coef .- z .* result.dissolution_se
diss_upper = result.dissolution_coef .+ z .* result.dissolution_se

# Print
for (i, term) in enumerate(result.model.formula.formation_terms)
    println("Formation $(name(term)): $(round(result.formation_coef[i], digits=3)) ",
            "[$(round(form_lower[i], digits=3)), $(round(form_upper[i], digits=3))]")
end
```

### Computing Hazard Ratios

```julia
# Hazard ratio for formation
for (i, term) in enumerate(result.model.formula.formation_terms)
    hr = exp(result.formation_coef[i])
    println("Formation $(name(term)): HR = $(round(hr, digits=3))")
end

# Hazard ratio for persistence
for (i, term) in enumerate(result.model.formula.dissolution_terms)
    hr = exp(result.dissolution_coef[i])
    println("Persistence $(name(term)): HR = $(round(hr, digits=3))")
end
```

## Model Comparison

### Comparing Nested Models

```julia
# Model 1: edges only
form1 = [Edges()]
diss1 = [Edges()]
r1 = stergm(networks, form1, diss1)

# Model 2: edges + mutual
form2 = [Edges(), Mutual()]
diss2 = [Edges(), Mutual()]
r2 = stergm(networks, form2, diss2)

# Compare (informal: check coefficient significance and model fit)
println("Model 1:")
println(r1)
println("\nModel 2:")
println(r2)
```

### Goodness-of-Fit

```julia
gof = stergm_gof(result; n_sim=100)

println("Observed statistics: ", gof.observed)
println("Simulated mean: ", gof.simulated_mean)
println("Simulated SD: ", gof.simulated_sd)
println("Z-scores: ", gof.z_scores)

# Z-scores close to 0 indicate good fit
# |Z| > 2 suggests model misfit
```

## Convergence Issues

### Checking Convergence

```julia
if !result.converged
    @warn "Model did not converge"
    println("Formation coefficients may be unreliable: ", result.formation_coef)
    println("Dissolution coefficients may be unreliable: ", result.dissolution_coef)
end
```

### Common Causes and Solutions

| Issue | Symptom | Solution |
|-------|---------|----------|
| Perfect separation | Very large coefficients, NaN SEs | Remove problematic term or simplify model |
| Too few transitions | Non-convergence | Use more time points |
| Sparse edges | Formation model issues | Reduce formation model complexity |
| Dense edges | Dissolution model issues | Reduce dissolution model complexity |
| Multicollinearity | Large standard errors | Remove correlated terms |

### Increasing Iterations

```julia
# Allow more iterations
result = stergm(networks, form_terms, diss_terms;
    maxiter=500,
    tol=1e-8
)
```

### Simplifying the Model

```julia
# Start with the simplest possible model
result_simple = stergm(networks, [Edges()], [Edges()])

# If that converges, add one term at a time
result_mutual = stergm(networks, [Edges(), Mutual()], [Edges()])
result_full = stergm(networks, [Edges(), Mutual()], [Edges(), Mutual()])
```

## The Alias: fit_stergm

`fit_stergm` is an alias for `stergm` provided for naming consistency:

```julia
# These are equivalent
result = stergm(networks, form_terms, diss_terms)
result = fit_stergm(networks, form_terms, diss_terms)
```

## Advanced Usage

### Custom Time Labels

Associate time labels with network observations:

```julia
# Networks observed at specific dates
times = [2020.0, 2021.0, 2022.0, 2023.0]
result = stergm(networks, form_terms, diss_terms; times=times)
```

### Accessing the Model Object

The `STERGMResult` contains the full model specification:

```julia
model = result.model

# Access networks
model.networks       # Vector of observed networks
model.directed       # Whether directed

# Access formula
model.formula.formation_terms    # Formation term list
model.formula.dissolution_terms  # Dissolution term list
```

## Best Practices

1. **Start with `Edges()` only**: Establish baseline formation and dissolution rates
2. **Add terms one at a time**: Check convergence after each addition
3. **Check coefficient magnitudes**: Very large coefficients (|coef| > 10) suggest problems
4. **Verify with simulation**: Simulated networks should match observed data characteristics
5. **Use sufficient time points**: Aim for at least 3-4 transitions (4-5 time points)
6. **Consider network size**: CMLE scales with $O(n^2 \times T)$ where $n$ is vertex count and $T$ is transitions
7. **Report both models**: Always present both formation and dissolution results together
8. **Interpret cautiously**: Different mechanisms can produce similar coefficient patterns
