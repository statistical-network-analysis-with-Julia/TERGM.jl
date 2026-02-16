# Simulation

This guide covers how to simulate dynamic network sequences from fitted STERGM models or from specified parameters.

## Overview

Network simulation from STERGMs generates plausible future network trajectories based on estimated formation and dissolution dynamics. Simulation serves three main purposes:

1. **Prediction**: Generate possible future network states
2. **Model assessment**: Compare simulated networks with observed data (goodness-of-fit)
3. **Scenario analysis**: Explore "what-if" scenarios by modifying parameters

## Simulation from Fitted Models

### Basic Simulation

After fitting a STERGM, simulate future networks:

```julia
using Network
using ERGM
using TERGM

# Fit a model
result = stergm(networks, form_terms, diss_terms)

# Simulate 10 future time steps
future_nets = simulate_stergm(result, 10)

# future_nets is a Vector{Network{Int}} of length 10
println("Simulated $(length(future_nets)) networks")
for (t, net) in enumerate(future_nets)
    println("t+$t: $(ne(net)) edges, density=$(round(ne(net)/(nv(net)*(nv(net)-1)), digits=3))")
end
```

### Starting Network

By default, simulation starts from the last observed network. You can specify a different starting point:

```julia
# Start from a specific network
start_net = networks[1]  # Start from the first observed network
future = simulate_stergm(result, 10; start_net=start_net)
```

### Burn-in Period

Discard initial simulation steps to reduce dependence on the starting network:

```julia
# Discard first 100 steps, then keep 10
future = simulate_stergm(result, 10; burnin=100)
```

The burn-in period is useful when:
- The starting network may not be representative of the equilibrium
- You want the simulated networks to reflect the estimated dynamics rather than initial conditions

## Simulation from Parameters

### Specifying Formation and Dissolution Parameters

Simulate without a fitted model by providing parameters directly:

```julia
# Define model specification
formula = STERGM([Edges(), Mutual()], [Edges()])

# Initial network
init_net = Network{Int}(; n=20, directed=true)
for i in 1:20, j in 1:20
    if i != j && rand() < 0.1
        add_edge!(init_net, i, j)
    end
end

# Specify coefficients
form_coef = [-3.0, 1.0]   # Edges intercept, Mutual effect
diss_coef = [2.0]          # Edges intercept (high persistence)

# Simulate
nets = simulate_network_sequence(formula, init_net, 50;
    form_coef=form_coef,
    diss_coef=diss_coef,
    burnin=200
)
```

### Understanding Parameters

The formation and dissolution coefficients control the dynamics:

| Parameter | Effect |
|-----------|--------|
| Formation `edges` = -3.0 | Low baseline formation rate (~5%) |
| Formation `edges` = -1.0 | Moderate formation rate (~27%) |
| Formation `edges` = 0.0 | 50% formation rate (very dense) |
| Dissolution `edges` = 2.0 | High persistence (~88%) |
| Dissolution `edges` = 0.0 | 50% persistence (high turnover) |
| Dissolution `edges` = -2.0 | Low persistence (~12%) |

## The Simulation Algorithm

### One-Step Process

Each simulation step (`simulate_one_step`) proceeds as follows:

1. **Formation step**: For each non-edge $(i,j)$ in the current network:
   - Compute formation change statistics $\Delta x^+(i,j)$
   - Calculate formation probability: $p^+_{ij} = \text{logistic}(\theta^+ \cdot \Delta x^+_{ij})$
   - Form edge $(i,j)$ with probability $p^+_{ij}$

2. **Dissolution step**: For each edge $(i,j)$ in the current network (at the start of the step):
   - Compute dissolution change statistics $\Delta x^-(i,j)$
   - Calculate persistence probability: $p^-_{ij} = \text{logistic}(\theta^- \cdot \Delta x^-_{ij})$
   - Dissolve edge $(i,j)$ with probability $1 - p^-_{ij}$

```julia
# The formation probability for dyad (i,j)
delta = [change_stat(term, net, i, j) for term in formation_terms]
eta = dot(form_coef, delta)
form_prob = 1.0 / (1.0 + exp(-eta))

# The persistence probability for dyad (i,j)
delta = [change_stat(term, net, i, j) for term in dissolution_terms]
eta = dot(diss_coef, delta)
persist_prob = 1.0 / (1.0 + exp(-eta))
dissolve_prob = 1.0 - persist_prob
```

### Equilibrium Density

The equilibrium density of a STERGM with `Edges()` only is approximately:

$$d^* \approx \frac{p^+}{p^+ + (1 - p^-)}$$

Where:
- $p^+ = \text{logistic}(\theta^+_{\text{edges}})$ is the formation probability
- $p^- = \text{logistic}(\theta^-_{\text{edges}})$ is the persistence probability

```julia
# Compute expected equilibrium density
form_prob = 1.0 / (1.0 + exp(-form_coef[1]))
pers_prob = 1.0 / (1.0 + exp(-diss_coef[1]))
eq_density = form_prob / (form_prob + (1 - pers_prob))
println("Expected equilibrium density: $(round(eq_density, digits=3))")
```

## Analyzing Simulated Networks

### Tracking Network Statistics Over Time

```julia
# Simulate a long sequence
nets = simulate_stergm(result, 100; burnin=50)

# Track density, reciprocity, and edge count
for (t, net) in enumerate(nets)
    n = nv(net)
    density = ne(net) / (n * (n - 1))

    mutual = 0
    for e in edges(net)
        if has_edge(net, dst(e), src(e))
            mutual += 1
        end
    end
    recip = ne(net) > 0 ? mutual / ne(net) : 0.0

    if t % 10 == 0
        println("t=$t: density=$(round(density, digits=3)), ",
                "reciprocity=$(round(recip, digits=3)), ",
                "edges=$(ne(net))")
    end
end
```

### Comparing with Observed Data

```julia
# Observed statistics from last network
obs_net = networks[end]
obs_density = ne(obs_net) / (nv(obs_net) * (nv(obs_net) - 1))

# Simulated statistics
sim_densities = [ne(net) / (nv(net) * (nv(net) - 1)) for net in nets]
sim_mean = mean(sim_densities)
sim_sd = std(sim_densities)

println("Observed density: $(round(obs_density, digits=3))")
println("Simulated mean: $(round(sim_mean, digits=3)) (SD: $(round(sim_sd, digits=3)))")
println("Z-score: $(round((obs_density - sim_mean) / sim_sd, digits=2))")
```

## Goodness-of-Fit

### Using stergm_gof

The `stergm_gof` function automates goodness-of-fit assessment:

```julia
gof = stergm_gof(result; n_sim=100)

println("Observed statistics: ", gof.observed)
println("Simulated means: ", gof.simulated_mean)
println("Simulated SDs: ", gof.simulated_sd)
println("Z-scores: ", gof.z_scores)
```

### Custom GOF Statistics

Provide custom functions to compute additional statistics:

```julia
using Graphs

gof = stergm_gof(result;
    n_sim=100,
    statistics=[
        ne,                          # Edge count
        x -> mean(degree(x)),        # Mean degree
        x -> nv(x) > 2 ? global_clustering_coefficient(x) : 0.0,  # Clustering
    ]
)
```

### Interpreting GOF Results

| Z-score | Interpretation |
|---------|---------------|
| |Z| < 1 | Good fit |
| 1 < |Z| < 2 | Marginal fit |
| |Z| > 2 | Poor fit for this statistic |

## Scenario Analysis

### Modifying Formation Rates

Explore how changes in formation dynamics affect network evolution:

```julia
# Original model
nets_original = simulate_stergm(result, 50)

# Increase formation rate by modifying coefficients
modified_result = deepcopy(result)
modified_result.formation_coef[1] += 0.5  # Increase baseline formation

nets_faster = simulate_stergm(modified_result, 50)

# Compare densities
d_orig = mean([ne(n) / (nv(n) * (nv(n) - 1)) for n in nets_original])
d_fast = mean([ne(n) / (nv(n) * (nv(n) - 1)) for n in nets_faster])

println("Original density: $(round(d_orig, digits=3))")
println("Faster formation density: $(round(d_fast, digits=3))")
```

### Modifying Dissolution Rates

```julia
# Increase persistence (decrease dissolution)
formula = STERGM([Edges()], [Edges()])

# Low persistence
nets_low = simulate_network_sequence(formula, init_net, 50;
    form_coef=[-3.0], diss_coef=[0.5], burnin=100)

# High persistence
nets_high = simulate_network_sequence(formula, init_net, 50;
    form_coef=[-3.0], diss_coef=[3.0], burnin=100)

# Compare edge stability
println("Low persistence - mean edges: ",
    round(mean([ne(n) for n in nets_low]), digits=1))
println("High persistence - mean edges: ",
    round(mean([ne(n) for n in nets_high]), digits=1))
```

## Reproducibility

### Setting Random Seeds

For reproducible simulations, set the random seed:

```julia
using Random

Random.seed!(42)
nets1 = simulate_stergm(result, 10)

Random.seed!(42)
nets2 = simulate_stergm(result, 10)

# nets1 and nets2 are identical
```

### Multiple Replications

Run multiple independent simulations:

```julia
using Random

replications = 20
all_densities = Vector{Float64}[]

for rep in 1:replications
    Random.seed!(rep)
    nets = simulate_stergm(result, 50; burnin=100)
    densities = [ne(n) / (nv(n) * (nv(n) - 1)) for n in nets]
    push!(all_densities, densities)
end

# Summary statistics across replications
final_densities = [d[end] for d in all_densities]
println("Mean final density: $(round(mean(final_densities), digits=3))")
println("SD final density: $(round(std(final_densities), digits=3))")
```

## Performance Considerations

### Computational Cost

Each simulation step requires iterating over all possible dyads:
- Formation: $O(n^2 - m)$ where $n$ is vertex count and $m$ is edge count
- Dissolution: $O(m)$

For a simulation of $S$ steps: total cost is $O(S \times n^2)$.

### Tips for Large Networks

1. **Shorter simulations**: Use fewer steps when possible
2. **Larger burn-in**: Ensure equilibrium is reached before collecting results
3. **Fewer replications**: Balance statistical precision with computation
4. **Simpler terms**: Complex terms (triangles) are slower to compute change statistics for

## Best Practices

1. **Always use burn-in**: At least 100 steps for burn-in before collecting results
2. **Set random seeds**: For reproducibility and debugging
3. **Check equilibrium**: Verify that density and other statistics stabilize after burn-in
4. **Compare with data**: Always validate simulated networks against observed data
5. **Multiple replications**: Run at least 10-20 replications for reliable GOF assessment
6. **Monitor convergence**: Track statistics over simulation steps to verify stationarity
