# Temporal Terms

```@meta
CurrentModule = TERGM
```

Formation and dissolution models accept two kinds of terms:

1. **Standard ERGM.jl terms** (`Edges`, `Mutual`, `Triangle`,
   `NodeMatch`, ...) — evaluated on the auxiliary network (``Y^+`` for
   formation, ``Y^-`` for dissolution) with their usual add-direction
   change statistics.
2. **Temporal terms** ([`TemporalTerm`](@ref)s) — statistics of the pair
   (current network, previous network), implementing
   `compute(term, net, prev_net)` and
   `change_stat(term, net, i, j, prev_net)` with the same
   **state-independent add-direction** convention (verified brute-force
   in the tests).

## Available temporal terms

- [`EdgeStability`](@ref): dyads agreeing with the previous panel;
  ``\Delta = +1`` if the dyad had an edge at t−1, else ``-1``.
- [`Delrecip`](@ref): delayed reciprocity — `i→j` now with `j→i` before.
- [`PersistentEdge`](@ref): edges present in both panels.
- [`NewEdge`](@ref): edges present now but not before.

Calling a temporal term without `prev_net` raises an informative error.

## Descriptives

[`edge_ages`](@ref) and [`mean_edge_age`](@ref) summarize tie duration in
a panel sequence (ages of the final panel's edges in consecutive panels).
