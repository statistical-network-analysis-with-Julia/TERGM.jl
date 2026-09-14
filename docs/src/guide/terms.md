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

### Usable in a formula

- [`Delrecip`](@ref): delayed reciprocity, ``\sum_{i \ne j} y_{ij}\,
  y^{t-1}_{ji}`` — `i→j` now with `j→i` before; change statistic
  ``\mathbf{1}[y^{t-1}_{ji}]``. Directed networks only:
  `requires_directed(Delrecip()) == true`, so model construction refuses it
  on undirected panels, and `compute`/`change_stat` throw an
  `ArgumentError` on an undirected network. **Provenance:** this is
  btergm's `delrecip` (Leifeld, Cranmer & Desmarais; `?btergm::"tergm-terms"`,
  lag 1) applied to the auxiliary networks of a separable model. R tergm
  ≥ 4 — the version the golden fixture is generated with — has **no**
  `delrecip` term (`InitErgmTerm.delrecip not found`), so its coefficient
  cannot be pinned against tergm; its definition and change statistic are
  pinned by hand in the test suite instead.

### Descriptives only — refused in a formula

- [`EdgeStability`](@ref): dyads agreeing with the previous panel;
  ``\Delta = +1`` if the dyad had an edge at t−1, else ``-1``.
- [`PersistentEdge`](@ref): edges present in both panels (``\Delta = 1``
  on a prior tie, ``0`` otherwise).
- [`NewEdge`](@ref): edges present now but not before (``\Delta = 1`` on a
  prior non-tie, ``0`` otherwise).

In a **separable** model these three carry no information a CMPLE could
estimate. The formation model's free dyads are exactly the prior non-ties
and the persistence model's exactly the prior ties, and each of these
statistics depends on the dyad's previous state alone — so on the rows of
either model its change statistic is one constant:

| term | formation rows (prior non-ties) | persistence rows (prior ties) |
|---|---|---|
| `edge.stability` | ``-1`` everywhere: the column is ``-``edges | ``+1`` everywhere: the column is edges |
| `persistent.edges` | ``0`` everywhere | ``+1`` everywhere: edges under another name |
| `new.edges` | ``+1`` everywhere: edges under another name | ``0`` everywhere |

With `Edges()` in the formula the design is rank-deficient; without it the
term *is* `edges`. They are informative only in a non-separable
(btergm `memory`-style) TERGM, which TERGM.jl does not fit — so
[`STERGMModel`](@ref) refuses them in either formula:

```julia
using TERGM, ERGM, Networks
t0 = network(6; directed=true); add_edge!(t0, 1, 2); add_edge!(t0, 3, 4)
t1 = copy(t0); add_edge!(t1, 2, 3); rem_edge!(t1, 3, 4)
try
    stergm([t0, t1], [Edges(), EdgeStability()], [Edges()])
catch e
    println(e.msg)   # "formation model: term 'edge.stability' is constant on the free dyads of the formation model of a separable STERGM (its change statistic is -1 on every prior non-tie, i.e. the column is -edges), ..."
end
compute(PersistentEdge(), t1, t0), compute(NewEdge(), t1, t0)   # (1.0, 1.0): fine as descriptives
```

(`gof` reports the `persisted`/`formed` counts through these two terms.)
None of the three has a counterpart in R tergm 4.

Calling a temporal term without `prev_net` raises an informative error.

Temporal terms take part in ERGM.jl's term-trait protocol
(`requires_directed`, `requires_undirected`, `required_vertex_attributes`,
`is_dyad_dependent`), so ERGM.jl's formula validation and
`has_dyad_dependent(model)` see them exactly as they see the built-in
terms. The four built-in temporal terms are dyad-independent: they
condition only on the exogenous previous network.

## Descriptives

[`edge_ages`](@ref) and [`mean_edge_age`](@ref) summarize tie duration in
a panel sequence (ages of the final panel's edges in consecutive panels).
