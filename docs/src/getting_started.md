# Getting Started

Fit separate formation and persistence effects to the bundled s50 friendship panel, inspect coefficient labels, and simulate model checks. The examples also show how to assemble a panel by hand while keeping the actor set consistent.

!!! note "Before you begin"

    The implemented estimator is conditional MPLE (CMPLE) for loop-free, one-mode panels with aligned actors. Missing dyads are refused. Conditional MLE, equilibrium moment estimation, and arbitrary nonseparable temporal models are not implemented. Dependent terms retain the limitations of pseudo-likelihood inference.

## Installation

```@raw html
<p>Use Julia <strong>1.12 or newer</strong> and the <a href="/getting-started/">shared workspace installation guide</a>. These development packages are not yet registered; the guide prepares the required sibling checkouts and a Julia environment for the examples.</p>
```

Run the blocks below in order in that environment. They build on variables from earlier steps; stochastic examples use seeded random number generators where shown.

## A real panel: s50

RSiena's `s50` excerpt — 50 girls, three yearly waves of directed friendship
nominations — ships with Networks.jl as `load_dataset(:s50)`. Wave-1 smoking
becomes a vertex attribute so `NodeMatch(:smoke)` can ask whether ties
between girls with the same smoking status form and persist differently:

```julia
using TERGM, ERGM, Networks, Random

s50 = load_dataset(:s50)                       # (friendship, alcohol, smoke)
networks = [copy(w) for w in s50.friendship]   # three Network{Int,true} panels
for w in networks, v in 1:nv(w)
    set_vertex_attribute!(w, :smoke, v, s50.smoke[v, 1])
end

# R tergm: tergm(nets ~ Form(~edges + mutual + nodematch("smoke")) +
#                       Persist(~edges + nodematch("smoke")), estimate = "CMPLE")
result = stergm(networks,
                [Edges(), Mutual(), NodeMatch(:smoke)],   # formation      (Form)
                [Edges(), NodeMatch(:smoke)])             # dissolution (persistence)  (Persist)
println(result)                # header with AIC/BIC, then Formation: / Persistence: blocks, one legend
result.model                   # STERGMModel{Int64,true}: 3 panels of 50 vertices (directed), 4900 free dyads; formation: ...

formation_coef(result)         # log-odds of a prior non-tie forming
persistence_coef(result)       # log-odds of a prior tie persisting (tergm Persist())

tbl = coeftable(result)        # Networks.CoefficientTable, tergm's row labels
tbl.names                      # ["Form(1)~edges", "Form(1)~mutual", ..., "Persist(1)~nodematch.smoke"]
tbl["Persist(1)~nodematch.smoke"].estimate == persistence_coef(result)[2]   # true
confint(result)                # Wald 95% limits, one row per coefficient
confint(result; level = 0.9)   # narrower

# The usual StatsAPI verbs work on the stacked (formation, persistence) vector
coef(result); stderror(result); vcov(result); loglikelihood(result); aic(result); bic(result)

# Transition-level goodness of fit — gof(fit) is the family's one verb
gof(result; n_sim = 100, rng = Xoshiro(1))
```

`Mutual` is dyad-dependent, so the printed fit carries a caveat: a CMPLE's
inverse-Hessian standard errors are anticonservative there. Refit with
`se = :bootstrap` for per-transition block-bootstrap ones (see
[STERGM Estimation](@ref)).

Coming from R tergm: `Form(~…)` is the second argument, `Persist(~…)` the
third, `estimate = "CMPLE"` is the only estimator, `summary(fit)` is
`coeftable(fit)`, and `Diss()`'s coefficients are `-persistence_coef(fit)`
(the README has the full correspondence table). The term lists are anything
`fit_ergm` accepts, and the common slips are said in words:

```julia
using TERGM, ERGM, Networks
s50 = load_dataset(:s50)
networks = [copy(w) for w in s50.friendship]
terms = []                                   # built programmatically: a Vector{Any} is fine
push!(terms, Edges()); push!(terms, Mutual())
stergm(networks, terms, Edges())             # a bare term needs no brackets
try
    stergm(terms, [Edges()], networks)       # the panel goes first
catch e
    println(e.msg)   # "arguments are swapped: call stergm(networks, formation, dissolution) — the panel ... comes first ..."
end
```

## Building a panel by hand

```julia
using TERGM, ERGM, Networks, Random

# A panel of same-sized directed networks (synthetic for the demo)
rng = Xoshiro(5)
net_t0 = network(20; directed=true)
for i in 1:20, j in 1:20
    i != j && rand(rng) < 0.1 && add_edge!(net_t0, i, j)
end
net_t1 = copy(net_t0); net_t2 = copy(net_t1)
for w in (net_t1, net_t2), _ in 1:15
    i, j = rand(rng, 1:20), rand(rng, 1:20)
    i == j && continue
    has_edge(w, i, j) ? rem_edge!(w, i, j) : add_edge!(w, i, j)
end
networks = [net_t0, net_t1, net_t2]

result = stergm(networks,
                [Edges(), Mutual()],   # formation
                [Edges()])             # dissolution (persistence)

formation_coef(result)     # log-odds of a prior non-tie forming
persistence_coef(result)   # log-odds of a prior tie persisting (tergm Persist())

# The auxiliary networks are available directly
yplus  = formation_network(networks[1], networks[2])
yminus = dissolution_network(networks[1], networks[2])

# Simulate 10 steps forward and check fit
future = simulate_stergm(result, 10; rng = rng)
gof(result; n_sim = 100, rng = rng)   # Networks.gof generic (stergm_gof is deprecated)
```

Formation/dissolution models mix standard ERGM.jl terms with temporal
terms conditioned on the previous panel:

```julia
stergm(networks, [Edges(), Delrecip()], [Edges(), Mutual()])
```

A multi-level `NodeFactor(:attr)`, a `NodeMix(:attr)` and a degree range
`Degree(0:2)` expand to one coefficient per level / cell / degree with R's
labels (`nodefactor.attr.b`, `mix.attr.a.b`, `degree1`), exactly as
`fit_ergm` expands them.

Both formulas are validated against every panel with ERGM.jl's own
formula checks before anything is fit, so the common mistakes are loud
`ArgumentError`s naming the side, the panel and the term: a vertex
attribute missing from a panel (or not set on every vertex), a
directed-only term (`Mutual`, `Delrecip`) on undirected panels, an
undirected-only term (`Kstar`, `GWDegree`, `Degree`) on directed panels —
use `OStar`/`IStar`, `GWODegree`/`GWIDegree`, `ODegree`/`IDegree` — or an
`EdgeCov` matrix of the wrong size:

```julia
try
    stergm(networks, [Edges(), Kstar(2)], [Edges()])
catch e
    println(e.msg)   # "formation model, panel 1: term 'kstar2' is only defined for undirected networks ... Use `OStar(2)` / `IStar(2)` ..."
end
```
