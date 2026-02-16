# Temporal Terms

This guide covers all temporal terms available in TERGM.jl, including their interpretation and usage in STERGM models.

## Term Interface

All temporal terms in TERGM.jl extend the ERGM.jl term interface. They subtype `TemporalTerm <: AbstractERGMTerm` and implement:

```julia
compute(term, net, prev_net) -> Float64      # Full statistic value
change_stat(term, net, i, j, prev_net) -> Float64  # Change when toggling edge (i,j)
name(term) -> String                          # Human-readable term name
```

The key difference from standard ERGM terms is that temporal terms receive the **previous network** (`prev_net`) in addition to the current network, allowing them to capture dynamics.

## Standard ERGM Terms in STERGM

Any term from ERGM.jl can be used directly in formation or dissolution models. Common choices:

```julia
using ERGM

# Structural terms
Edges()       # Intercept (baseline rate)
Triangle()    # Transitivity (triadic closure)
GWESP(0.5)    # Geometrically weighted edgewise shared partners
GWDegree(0.5) # Geometrically weighted degree distribution

# Nodal terms
Mutual()        # Reciprocity
NodeMatch(:attr) # Homophily
NodeCov(:attr)   # Covariate main effect

# These terms ignore prev_net and compute only on the current network
```

When used in a formation model, the term is computed on the current network state restricted to non-edges at $t-1$. In a dissolution model, it is computed on the current network state restricted to edges at $t-1$.

## Edge Dynamics Terms

### EdgeStability

Counts edges that persist from the previous time step. This is the most fundamental temporal term.

```julia
EdgeStability()
```

**Statistic**: Number of edges present in both $Y_{t-1}$ and $Y_t$.

**Interpretation**: In the dissolution model, a positive coefficient indicates that edges are more stable when many other edges also persist (network inertia). In the formation model, this term has no effect (new edges cannot be "stable").

```julia
# Usage in dissolution model
diss_terms = [Edges(), EdgeStability()]
```

### PersistentEdge

Functionally equivalent to `EdgeStability`. Counts edges present in both current and previous networks.

```julia
PersistentEdge()
```

**Change statistic**: For edge $(i,j)$:
- If $(i,j)$ was in $Y_{t-1}$: toggling it changes the count by $\pm 1$
- If $(i,j)$ was not in $Y_{t-1}$: change is $0$ (only existing edges can persist)

### NewEdge

Counts edges present in the current network but not in the previous network. The complement of `PersistentEdge`.

```julia
NewEdge()
```

**Statistic**: Number of edges in $Y_t$ that were not in $Y_{t-1}$.

**Interpretation**: In the formation model, this acts as an alternative intercept counting only newly formed edges. In practice, `Edges()` is usually preferred as the formation intercept.

```julia
# Track new edge formation explicitly
form_terms = [NewEdge(), Triangle()]
```

## Memory Terms

### Memory

A weighted stability term where each persisting edge contributes a fixed weight $\theta$:

```julia
Memory(1.0)     # Default weight
Memory(0.5)     # Half weight per persisting edge
Memory(2.0)     # Double weight per persisting edge
```

**Statistic**: $\theta \times$ (number of edges in both $Y_t$ and $Y_{t-1}$).

**Interpretation**: Controls the strength of the memory effect. Higher $\theta$ values make past network structure more influential on the current state.

```julia
# Strong memory effect
diss_terms = [Edges(), Memory(2.0)]
```

### EdgeAge

Captures the age of edges -- how many consecutive time steps they have persisted. This requires the full network sequence, not just the previous time step.

```julia
EdgeAge()           # Default: no maximum age
EdgeAge(max_age=5)  # Cap age at 5
```

**Interpretation**: In the dissolution model, a positive coefficient indicates that older edges are more likely to persist ("the longer a tie exists, the stronger it becomes").

## Reciprocity Terms

### Delrecip (Delayed Reciprocity)

Counts edges $(i,j)$ in the current network where the reverse edge $(j,i)$ existed in the previous network. This captures delayed reciprocity -- the tendency to reciprocate edges from the previous time step.

```julia
Delrecip()
```

**Statistic**: $\sum_{(i,j) \in Y_t} \mathbb{1}[(j,i) \in Y_{t-1}]$

**Interpretation**: In the formation model, a positive coefficient means that actors tend to reciprocate edges that existed in the previous period.

```julia
# Test for delayed reciprocity in formation
form_terms = [Edges(), Delrecip()]
diss_terms = [Edges()]

result = stergm(networks, form_terms, diss_terms)
# Positive Delrecip coefficient → reciprocity drives formation
```

**Note**: This term is only meaningful for directed networks. For undirected networks, it always returns 0.

## Lagged Terms

### TimeLag

Apply any standard ERGM term to the **previous** network rather than the current one. This creates a lagged effect.

```julia
# Lagged edge count
TimeLag(Edges())

# Lagged triangle count
TimeLag(Triangle())

# With custom lag (default is 1)
TimeLag(Edges(), 1)
```

**Statistic**: `compute(wrapped_term, prev_net)` -- the term is computed on $Y_{t-1}$.

**Interpretation**: Tests whether the previous network's structural properties predict current network changes. For example, `TimeLag(Triangle())` in the formation model tests whether networks with more triangles at $t-1$ tend to form more edges at $t$.

```julia
# Does past transitivity affect current formation?
form_terms = [Edges(), Triangle(), TimeLag(Triangle())]
diss_terms = [Edges()]
```

## Formation and Dissolution Wrappers

### FormationTerm

Wraps a standard ERGM term to only apply when the dyad is eligible for formation (non-edge at $t-1$):

```julia
FormationTerm(Edges())
FormationTerm(Triangle())
FormationTerm(Mutual())
```

**Change statistic**: Returns `change_stat(wrapped_term, net, i, j)` if $(i,j)$ was not in $Y_{t-1}$, and $0$ otherwise.

**Use case**: Explicitly mark terms as belonging to the formation model. This is typically unnecessary since `stergm` already separates formation and dissolution, but it can be useful for combined models.

### DissolutionTerm

Wraps a standard ERGM term to only apply when the dyad is eligible for dissolution (edge at $t-1$):

```julia
DissolutionTerm(Edges())
DissolutionTerm(Mutual())
```

**Change statistic**: Returns `change_stat(wrapped_term, net, i, j)` if $(i,j)$ was in $Y_{t-1}$, and $0$ otherwise.

## Combining Terms in Models

### Basic Formation-Dissolution Model

```julia
# Formation: density + reciprocity
form_terms = [Edges(), Mutual()]

# Dissolution: density only (simple persistence)
diss_terms = [Edges()]

result = stergm(networks, form_terms, diss_terms)
```

### Rich Model with Temporal Effects

```julia
# Formation: new edges driven by reciprocity, transitivity, and delayed reciprocity
form_terms = [
    Edges(),           # Baseline formation rate
    Mutual(),          # Reciprocity
    Triangle(),        # Triadic closure
    Delrecip(),        # Delayed reciprocity from previous time
]

# Dissolution: persistence driven by reciprocity and stability
diss_terms = [
    Edges(),           # Baseline persistence
    Mutual(),          # Mutual ties persist longer
    EdgeStability(),   # Network inertia
]

result = stergm(networks, form_terms, diss_terms)
```

### Model with Node Attributes

```julia
using ERGM

# Define node attribute
gender = NodeAttribute(:gender, Dict(1=>"M", 2=>"F"), "Unknown")

form_terms = [
    Edges(),
    Mutual(),
    NodeMatch(gender),  # Homophily in formation
]

diss_terms = [
    Edges(),
    NodeMatch(gender),  # Homophily in persistence
]

result = stergm(networks, form_terms, diss_terms)
```

## Term Summary Table

| Term | Type | Description | Meaningful In |
|------|------|-------------|---------------|
| `Edges()` | Standard | Edge count (intercept) | Both |
| `Mutual()` | Standard | Mutual edges | Both |
| `Triangle()` | Standard | Triangle count | Both |
| `EdgeStability()` | Temporal | Persisting edge count | Dissolution |
| `PersistentEdge()` | Temporal | Same as EdgeStability | Dissolution |
| `NewEdge()` | Temporal | Newly formed edge count | Formation |
| `Memory(theta)` | Temporal | Weighted persistence | Dissolution |
| `EdgeAge()` | Temporal | Age of persisting edges | Dissolution |
| `Delrecip()` | Temporal | Delayed reciprocity | Formation |
| `TimeLag(term)` | Temporal | Lagged term on previous network | Both |
| `FormationTerm(term)` | Wrapper | Formation-only wrapper | Formation |
| `DissolutionTerm(term)` | Wrapper | Dissolution-only wrapper | Dissolution |

## Choosing Terms

### By Research Question

| Question | Recommended Terms |
|----------|-------------------|
| What drives tie formation? | `Edges`, `Mutual`, `Triangle`, `NodeMatch` |
| Why do ties persist? | `Edges`, `Mutual`, `EdgeStability`, `Memory` |
| Is there delayed reciprocity? | `Delrecip` in formation model |
| Does past structure matter? | `TimeLag(Triangle)`, `TimeLag(Mutual)` |
| Do similar actors form ties? | `NodeMatch` in formation model |
| Do similar actors maintain ties? | `NodeMatch` in dissolution model |

### Best Practices

1. **Always include `Edges()`**: This is the intercept for both formation and dissolution
2. **Start simple**: Begin with `Edges()` only, then add terms one at a time
3. **Check convergence**: Complex models may fail to converge with limited data
4. **Consider separation**: Formation and dissolution models operate on different subsets of dyads; some terms may have no variation in one subset
5. **Balance complexity**: More terms require more transitions (time points) for reliable estimation
6. **Use temporal terms sparingly**: They add complexity; only include them when substantively motivated
