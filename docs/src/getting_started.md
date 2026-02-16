# Getting Started

This tutorial walks through common use cases for TERGM.jl, from preparing panel network data to fitting and interpreting STERGMs.

## Installation

Install TERGM.jl from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/Statistical-network-analysis-with-Julia/TERGM.jl")
```

TERGM.jl depends on ERGM.jl, Network.jl, and NetworkDynamic.jl, which will be installed automatically.

## Basic Workflow

The typical TERGM.jl workflow consists of five steps:

1. **Prepare panel data** - Organize your network observations across time
2. **Define formation terms** - Choose what drives new tie creation
3. **Define dissolution terms** - Choose what drives tie persistence
4. **Fit the STERGM** - Estimate separate formation and dissolution parameters
5. **Interpret and simulate** - Analyze results and simulate future networks

## Step 1: Prepare Panel Data

STERGMs require a sequence of static networks observed at discrete time points:

```julia
using Network
using ERGM
using TERGM

# Create networks at 4 time points
n = 10  # 10 actors

net1 = Network{Int}(; n=n, directed=true)
add_edge!(net1, 1, 2)
add_edge!(net1, 2, 3)
add_edge!(net1, 3, 1)

net2 = Network{Int}(; n=n, directed=true)
add_edge!(net2, 1, 2)   # persists from t1
add_edge!(net2, 2, 3)   # persists from t1
add_edge!(net2, 3, 1)   # persists from t1
add_edge!(net2, 4, 1)   # new edge
add_edge!(net2, 2, 4)   # new edge

net3 = Network{Int}(; n=n, directed=true)
add_edge!(net3, 1, 2)   # persists
add_edge!(net3, 4, 1)   # persists
add_edge!(net3, 2, 4)   # persists
add_edge!(net3, 3, 4)   # new edge
# edges (2,3) and (3,1) dissolved

net4 = Network{Int}(; n=n, directed=true)
add_edge!(net4, 1, 2)   # persists
add_edge!(net4, 4, 1)   # persists
add_edge!(net4, 3, 4)   # persists
add_edge!(net4, 5, 1)   # new edge
# edge (2,4) dissolved

networks = [net1, net2, net3, net4]
```

### Requirements

- All networks must have the **same number of vertices** (the actor set is fixed)
- At least **two time points** are required (one transition)
- Networks should be either all directed or all undirected

## Step 2: Define Formation Terms

Formation terms describe what drives new edges to form. These are standard ERGM terms applied to the subset of dyads that do not have edges at time $t-1$:

```julia
# Simple density-based formation
form_terms = [Edges()]

# Formation driven by reciprocity
form_terms = [Edges(), Mutual()]

# Formation driven by triadic closure
form_terms = [Edges(), Triangle()]
```

The `Edges()` term acts as an intercept for the formation model and controls the baseline formation rate.

## Step 3: Define Dissolution Terms

Dissolution terms describe what drives existing edges to persist (not dissolve). These are applied to the subset of dyads that have edges at time $t-1$:

```julia
# Simple persistence
diss_terms = [Edges()]

# Persistence driven by mutual ties
diss_terms = [Edges(), Mutual()]
```

A positive coefficient in the dissolution model means the statistic **increases persistence** (decreases dissolution).

## Step 4: Fit the STERGM

Use the `stergm` function to estimate the model:

```julia
result = stergm(networks, form_terms, diss_terms)
println(result)
```

### Output

```text
STERGM Results
==============
Method: cmle
Converged: true
Log-likelihood: NaN

Formation Model:
----------------------------------------
  edges                   -3.2456 (SE: 0.4523)
  triangle                 0.8912 (SE: 0.3201)

Dissolution (Persistence) Model:
----------------------------------------
  edges                    1.5678 (SE: 0.5012)
```

### Estimation Methods

```julia
# Conditional Maximum Likelihood (default)
result = stergm(networks, form_terms, diss_terms; method=:cmle)

# Conditional Maximum Pseudo-Likelihood (faster for large networks)
result = stergm(networks, form_terms, diss_terms; method=:cmple)

# Equilibrium Generalized Method of Moments (experimental)
result = stergm(networks, form_terms, diss_terms; method=:egmme)
```

| Method | Description | Speed | Accuracy |
|--------|-------------|-------|----------|
| `:cmle` | Conditional Maximum Likelihood | Moderate | High |
| `:cmple` | Conditional Maximum Pseudo-Likelihood | Fast | Good |
| `:egmme` | Equilibrium Generalized Method of Moments | Slow | Experimental |

## Step 5: Interpret Results

### Accessing Coefficients

```julia
# Formation coefficients and standard errors
result.formation_coef      # Vector of formation coefficients
result.formation_se        # Vector of formation standard errors

# Dissolution coefficients and standard errors
result.dissolution_coef    # Vector of dissolution coefficients
result.dissolution_se      # Vector of dissolution standard errors

# Model information
result.converged           # Whether estimation converged
result.method              # Estimation method used
result.loglik              # Log-likelihood
```

### Interpreting Formation Coefficients

Coefficients are **log-odds ratios** for edge formation:

| Coefficient | Interpretation |
|-------------|----------------|
| Formation `edges` = -3.0 | Baseline formation probability = exp(-3)/(1+exp(-3)) = 0.047 |
| Formation `mutual` = 1.5 | Reciprocal edges are exp(1.5) = 4.5x more likely to form |
| Formation `triangle` = 0.8 | Each shared partner increases formation odds by exp(0.8) = 2.2x |

### Interpreting Dissolution Coefficients

Coefficients are **log-odds ratios** for edge persistence:

| Coefficient | Interpretation |
|-------------|----------------|
| Dissolution `edges` = 2.0 | Baseline persistence probability = exp(2)/(1+exp(2)) = 0.88 |
| Dissolution `mutual` = 1.0 | Mutual edges are exp(1) = 2.7x more likely to persist |

### Computing Probabilities

```julia
# Baseline formation probability (from edges coefficient)
form_prob = 1.0 / (1.0 + exp(-result.formation_coef[1]))
println("Baseline formation prob: ", round(form_prob, digits=3))

# Baseline persistence probability (from edges coefficient)
pers_prob = 1.0 / (1.0 + exp(-result.dissolution_coef[1]))
println("Baseline persistence prob: ", round(pers_prob, digits=3))

# Dissolution probability = 1 - persistence probability
diss_prob = 1.0 - pers_prob
println("Baseline dissolution prob: ", round(diss_prob, digits=3))
```

## Complete Example

```julia
using Network
using ERGM
using TERGM
using Random

Random.seed!(42)

# Create a sequence of networks
n = 15
networks = Network{Int}[]

# Start with a sparse network
net = Network{Int}(; n=n, directed=true)
for i in 1:n
    for j in 1:n
        i == j && continue
        if rand() < 0.1  # 10% density
            add_edge!(net, i, j)
        end
    end
end
push!(networks, net)

# Evolve the network through several time points
for t in 2:5
    prev = networks[end]
    curr = Network{Int}(; n=n, directed=true)

    for i in 1:n, j in 1:n
        i == j && continue

        if has_edge(prev, i, j)
            # Existing edge: persist with probability 0.8
            if rand() < 0.8
                add_edge!(curr, i, j)
            end
        else
            # Non-edge: form with probability 0.05
            if rand() < 0.05
                add_edge!(curr, i, j)
            end
        end
    end

    push!(networks, curr)
end

# Fit STERGM
form_terms = [Edges()]
diss_terms = [Edges()]

result = stergm(networks, form_terms, diss_terms)
println(result)

# Check convergence
if result.converged
    println("\nModel converged successfully")

    # Formation rate
    form_prob = 1.0 / (1.0 + exp(-result.formation_coef[1]))
    println("Estimated formation prob: ", round(form_prob, digits=3))
    println("True formation prob: 0.05")

    # Persistence rate
    pers_prob = 1.0 / (1.0 + exp(-result.dissolution_coef[1]))
    println("Estimated persistence prob: ", round(pers_prob, digits=3))
    println("True persistence prob: 0.80")
else
    println("\nWarning: Model did not converge")
end
```

## Simulating Future Networks

Generate future network trajectories from a fitted model:

```julia
# Simulate 10 future time steps
future_nets = simulate_stergm(result, 10)

# Track network density over time
for (t, net) in enumerate(future_nets)
    density = ne(net) / (nv(net) * (nv(net) - 1))
    println("t+$t: density = $(round(density, digits=3))")
end
```

## Using Temporal Terms

TERGM.jl provides specialized temporal terms:

```julia
# Edge stability: edges persisting from previous time
form_terms = [Edges()]
diss_terms = [Edges(), EdgeStability()]

# Delayed reciprocity: reciprocate edges from previous time
form_terms = [Edges(), Delrecip()]
diss_terms = [Edges()]

# Wrap standard ERGM terms with formation/dissolution markers
form_terms = [FormationTerm(Edges()), FormationTerm(Triangle())]
diss_terms = [DissolutionTerm(Edges())]
```

See [Temporal Terms](guide/terms.md) for a complete guide.

## Goodness-of-Fit

Assess model fit by comparing observed and simulated network statistics:

```julia
gof = stergm_gof(result; n_sim=100)
println("Observed: ", gof.observed)
println("Simulated mean: ", gof.simulated_mean)
println("Z-scores: ", gof.z_scores)
```

## Best Practices

1. **Start simple**: Begin with `Edges()` in both formation and dissolution models
2. **Add terms gradually**: Add one term at a time and check convergence
3. **Check convergence**: Always verify `result.converged == true`
4. **Multiple time points**: More transitions provide better estimates (aim for 3+ time points)
5. **Consistent vertex set**: All networks must have the same vertices
6. **Validate with simulation**: Compare simulated and observed network statistics
7. **Consider model complexity**: Avoid over-parameterized models relative to data

## Next Steps

- Learn about [Temporal Terms](guide/terms.md) available for STERGM modeling
- Understand [STERGM Estimation](guide/estimation.md) methods in detail
- Master [Simulation](guide/simulation.md) from fitted models
