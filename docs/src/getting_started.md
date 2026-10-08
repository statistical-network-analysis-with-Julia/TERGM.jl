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
nominations — ships with NetworkCore.jl as `load_dataset(:s50)`. Wave-1 smoking
becomes a vertex attribute so `NodeMatch(:smoke)` can ask whether ties
between girls with the same smoking status form and persist differently:

```julia
using TERGM, ERGM, NetworkCore, Random

s50 = load_dataset(:s50)                       # (friendship, alcohol, smoke)
networks = [copy(w) for w in s50.friendship]   # three Network{Int,true} panels
for w in networks, v in 1:nv(w)
    set_vertex_attribute!(w, :smoke, v, s50.smoke[v, 1])
end

# R tergm: tergm(nets ~ Form(~edges + mutual + nodematch("smoke")) +
#                       Persist(~edges + nodematch("smoke")), estimate = "CMPLE")
result = stergm(networks,
                [Edges(), Mutual(), NodeMatch(:smoke)],   # formation      (Form)
                [Edges(), NodeMatch(:smoke)];             # dissolution (persistence)  (Persist)
                method = :cmple)
println(result)                # header with AIC/BIC, then Formation: / Persistence: blocks, one legend
result.model                   # STERGMModel{Int64,true}: 3 panels of 50 vertices (directed), 4900 free dyads; formation: ...

formation_coef(result)         # log-odds of a prior non-tie forming
persistence_coef(result)       # log-odds of a prior tie persisting (tergm Persist())

tbl = coeftable(result)        # NetworkCore.CoefficientTable, tergm's row labels
tbl.names                      # ["Form(1)~edges", "Form(1)~mutual", ..., "Persist(1)~nodematch.smoke"]
tbl["Persist(1)~nodematch.smoke"].estimate == persistence_coef(result)[2]   # true
result.inference_withheld      # true: `Mutual` is dyad-dependent — no z, p or interval by default

# Inference for a dyad-dependent formula: bootstrap the CMPLE ...
boot = stergm(networks,
              [Edges(), Mutual(), NodeMatch(:smoke)],
              [Edges(), NodeMatch(:smoke)];
              method = :cmple, se = :bootstrap, n_boot = 200, rng = Xoshiro(2))
confint(boot)                  # Wald 95% limits, one row per coefficient
confint(boot; level = 0.9)     # narrower

# ... or fit the conditional MLE by MCMC — tergm's estimate = "CMLE", and the
# default here (method = :auto) because `Mutual` is dyad-dependent
mle = stergm(networks,
             [Edges(), Mutual(), NodeMatch(:smoke)],
             [Edges(), NodeMatch(:smoke)];
             rng = Xoshiro(3))
mle.method                     # :cmle
println(mle)                   # log-likelihood, MCMC diagnostics, both tables
mle.mcmc.persistence.exact     # true: that side is dyad-independent, fit without MCMC

# The usual StatsAPI verbs work on the stacked (formation, persistence) vector
coef(mle); stderror(mle); vcov(mle); loglikelihood(mle); aic(mle); bic(mle)

# Transition-level goodness of fit — gof(fit) is the family's one verb
gof(mle; n_sim = 100, rng = Xoshiro(1))
```

`Mutual` is dyad-dependent, so the CMPLE (`method = :cmple`, tergm's
`estimate = "CMPLE"`) is an approximation whose naive standard errors are
not calibrated: it prints estimates and standard errors, but no z or p,
and `confint` refuses. The parametric bootstrap (`se = :bootstrap`) gives
calibrated inference for the CMPLE; the default, `method = :auto`, fits
the conditional MLE itself here, as tergm does (see [STERGM Estimation](@ref)).
For a formula with no dyad-dependent term `:auto` fits the CMPLE, which is
then the exact conditional MLE.

Coming from R tergm: `Form(~…)` is the second argument, `Persist(~…)` the
third, `estimate = "CMLE"` (which R requires you to write: `tergm()` has
no default estimator) is TERGM.jl's default `method = :auto` for a
dyad-dependent formula (or explicitly `method = :cmle`) and `estimate = "CMPLE"`
is `method = :cmple`, `summary(fit)` is
`coeftable(fit)`, and `Diss()`'s coefficients are `-persistence_coef(fit)`
(the README has the full correspondence table). The term lists are anything
`fit_ergm` accepts, and the common slips are said in words:

```julia
using TERGM, ERGM, NetworkCore
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
using TERGM, ERGM, NetworkCore, Random

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
                [Edges()];             # dissolution (persistence)
                method = :cmple)       # the fast pseudo-likelihood fit

formation_coef(result)     # log-odds of a prior non-tie forming
persistence_coef(result)   # log-odds of a prior tie persisting (tergm Persist())

# The auxiliary networks are available directly
yplus  = formation_network(networks[1], networks[2])
yminus = dissolution_network(networks[1], networks[2])

# Simulate 10 steps forward and check fit
future = simulate_stergm(result, 10; rng = rng)
gof(result; n_sim = 100, rng = rng)   # the NetworkCore.gof generic
```

Formation/dissolution models mix standard ERGM.jl terms with temporal
terms conditioned on the previous panel:

```julia
stergm(networks, [Edges(), Delrecip()], [Edges(), Mutual()]; method = :cmple)
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
