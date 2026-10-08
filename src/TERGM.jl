"""
    TERGM.jl - Temporal Exponential Random Graph Models

Separable Temporal ERGMs (STERGM) for networks observed at discrete time
points (Krivitsky & Handcock 2014), following R `tergm`:

- the transition Y_{t-1} → Y_t factors into a **formation** model on the
  formation network Y⁺ = Y_{t-1} ∪ Y_t and a **dissolution** model on the
  dissolution network Y⁻ = Y_{t-1} ∩ Y_t;
- formation statistics are evaluated on Y⁺ (free dyads: the non-edges of
  Y_{t-1}), dissolution statistics on Y⁻ (free dyads: the edges of
  Y_{t-1});
- estimation is by CMPLE — conditional maximum *pseudo*-likelihood over
  the free dyads of the auxiliary networks, pooled across transitions (it
  equals the CMLE exactly for dyad-independent terms and is an
  approximation otherwise; `method = :cmple`) — or by CMLE, the
  Monte-Carlo conditional maximum likelihood on the same constrained
  sample spaces (`method = :cmle`). The default `method = :auto` is R
  tergm's choice: the CMLE for a dyad-dependent formula, the (exact) CMPLE
  for a dyad-independent one;
- simulation draws Y⁺ and Y⁻ with Metropolis samplers constrained to the
  free dyads and combines them into Y_t.

The dissolution model is parameterized in terms of **persistence**
(tergm's `Persist()`): positive coefficients mean ties are more likely to
persist. The fitted coefficients are `persistence_coef`; tergm's `Diss()`
parameterisation is their negation.
"""
module TERGM

using Distributions
using ERGM
using Graphs
using LinearAlgebra
using NetworkCore
using PrecompileTools: @setup_workload, @compile_workload
using Random
using Statistics

# ----------------------------------------------------------------------------
# Shared contracts, imported BY NAME from the package that owns them.
#
# Every verb or rule that more than one package needs is defined once in
# NetworkCore.jl (or, for the term protocol, in ERGM.jl) and extended here with
# methods for TERGM's own types — never re-implemented locally. The "No
# private cross-package reach-ins" testset reads this file as text and
# asserts that no `_`-prefixed name is taken from ERGM or NetworkCore, by
# `import` or by an `ERGM._x`/`NetworkCore._x` dotted reach-in: ERGM.jl's
# building blocks come from its extension API, `ERGM.Extension`.
# ----------------------------------------------------------------------------

# NetworkCore.jl: the shared numerics (the ONE bootstrap loop, the ONE z → p
# helper, the ONE `se=` validator; `newton_fit`/`logistic_derivatives`, the
# kernel behind ERGM's design fitter, are reached only through
# `ERGM.Extension.mple_fit_design`), the ONE `gof` generic (every model
# package adds methods for its own result types, so `gof(fit)` works uniformly and
# loading several model packages never collides on the name), and the seven
# accessors of the result-metadata protocol (`fit_metadata(fit)` collects
# them).
import NetworkCore: z_pvalues, check_se, bootstrap_cov
# The ecosystem's one threaded-loop helper: waits for every task and rethrows
# the first failure as the task's own exception, never a TaskFailedException
import NetworkCore: spawn_all
import NetworkCore: gof
import NetworkCore: estimand, objective, is_exact, se_method, missing_method,
                 approximations

# ERGM.jl: the term protocol (`name`, `compute`, `change_stat` and the
# traits TERGM extends for its temporal terms and its model type) and the
# exported Metropolis toggle kernel.
import ERGM: name, compute, change_stat, is_dyad_dependent,
             has_dyad_dependent, requires_directed, required_vertex_attributes,
             mh_toggle!, TermSet
# ERGM.jl's extension API (`ERGM.Extension`, the semver-covered surface for
# packages built on ERGM.jl), imported by name: formula validation (so
# TERGM's panels are checked by exactly the code that checks an `ERGMModel`),
# term materialization (so a multi-level
# `NodeFactor`/`NodeMix` or a `Degree(0:2)` expands into one statistic per
# level/cell/degree with R's labels, and the nodal terms snapshot their
# attribute into typed vectors, exactly as an `ERGMModel` does) and the
# dyad-scaled sampler defaults.
import ERGM.Extension: validate_formula, materialize, mcmc_defaults
# ERGM.jl's compressed-row pseudo-likelihood fitter: R ergm's
# boundary-statistic drop, the separation verdict (`converged=false`) and
# the two R-sentence warnings, prefixed with the caller's name through
# `context="cmple"`. CMPLE is the same logistic fit on the free dyads of
# Y⁺/Y⁻, so it fails in exactly the same way as `ERGM.mple`.
import ERGM.Extension: mple_fit_design
# The ONE term-list normaliser (`fit_ergm` runs it): a `Vector{Any}` built
# with `push!`, a bare term, a nested vector of terms are all accepted, and a
# non-term element is refused with its position, type and the "did you mean
# `Edges()`?" hint — so `stergm` fails exactly as `fit_ergm` does. And the
# ONE term expansion as a *specification* (`expand_terms`: levels, cells and
# degrees resolved against a network, as plain terms without a snapshot),
# which a model that scores T−1 auxiliary networks needs (see below).
import ERGM.Extension: collect_terms, expand_terms
# The Monte-Carlo MLE of ERGM.jl's extension API, imported by name: the
# iteration itself (`mcmle_solve`, written against a sampling callback —
# Hummel-stepped update, R ergm's `confidence` stopping rule with its
# sample-size boost, Fisher + Monte-Carlo covariance) and the path-sampling
# bridge integral (`bridge_integrate`, Simpson's rule). The CMLE supplies
# only what is TERGM's: the constrained sampler pooled over transitions and
# the exact dyad-independent reference of the bridge.
import ERGM.Extension: mcmle_solve, bridge_integrate
import StatsAPI
# The full StatsAPI surface (`check_statsapi(fit; strict=true)` pins it):
# `coeftable` and `confint` are the SAME bindings NetworkCore.jl and ERGM.jl
# re-export, so `using ERGM, TERGM` leaves every verb single-owner.
import StatsAPI: coef, stderror, vcov, confint, loglikelihood, aic, bic, nobs,
                 dof, coeftable, coefnames

# Model types
export STERGM, STERGMResult, STERGMModel

# Temporal terms (evaluated with reference to the previous network)
export EdgeStability, Delrecip, PersistentEdge, NewEdge

# Auxiliary-network construction
export formation_network, dissolution_network

# Estimation
# `egmme` is deliberately NOT exported: it is unimplemented and can only
# throw, and exporting a name advertises a capability. It remains reachable
# as `TERGM.egmme` so the error is informative rather than an UndefVarError.
export stergm, fit_stergm, cmple, cmle

# Coefficient accessors. The dissolution model is fit in tergm's `Persist()`
# parameterisation, so its coefficients are `persistence_*`.
export formation_coef, formation_se, persistence_coef, persistence_se

# Simulation
export simulate_stergm, simulate_network_sequence

# Diagnostics: `gof` is THE verb — NetworkCore.jl's shared generic, extended
# with a method for STERGMResult.
export gof

# Temporal descriptives
export edge_ages, mean_edge_age

# StatsAPI methods (re-exported so `coef(fit)` etc. work with just `using TERGM`;
# `coefnames` is the StatsAPI binding NetworkCore and ERGM export too)
export coef, stderror, vcov, confint, loglikelihood, aic, bic, nobs, dof, coeftable,
       coefnames

# =============================================================================
# Temporal Terms
# =============================================================================
#
# Temporal terms are statistics of (Y, Y_{t-1}) pairs. They implement
#   compute(term, net, prev_net)::Float64
#   change_stat(term, net, i, j, prev_net)::Float64
# where change_stat is the ADD-DIRECTION change g(y⁺ij) − g(y⁻ij) in `net`
# holding prev_net fixed — state-independent, exactly like ERGM.jl's
# convention. Standard ERGM terms may be mixed freely with temporal terms
# in formation/dissolution models; they are evaluated on the auxiliary
# network without reference to prev_net.

"""
    TemporalTerm <: AbstractERGMTerm

Base type for statistics of a network *transition*: they depend on the
current network and the previous one.
"""
abstract type TemporalTerm <: AbstractERGMTerm end

# The two-argument forms a temporal term cannot answer. Typed on the abstract
# term ONLY (no `::Int` on `i, j`, as ERGM.jl's fallbacks): a user term's
# untyped `change_stat(::MyTerm, net, i, j)` is then strictly more specific
# and never ambiguous with this method.
compute(t::TemporalTerm, net) = throw(ArgumentError(
    "$(name(t)) is a temporal term: its statistic depends on the previous " *
    "network, so call compute(term, net, prev_net) (or define " *
    "`compute(::$(nameof(typeof(t))), net)` if the term does not need it)"))
change_stat(t::TemporalTerm, net, i, j) = throw(ArgumentError(
    "$(name(t)) is a temporal term: its change statistic depends on the " *
    "previous network, so call change_stat(term, net, i, j, prev_net) (or " *
    "define `change_stat(::$(nameof(typeof(t))), net, i, j)` if the term does " *
    "not need it)"))

# Evaluate any term's change statistic / statistic on a transition: a
# temporal term sees the previous panel, a standard ERGM term does not.
# `_tcompute` is what `gof` uses to evaluate the formation and persistence
# model statistics on the observed and simulated auxiliary networks.
_tchange(t::TemporalTerm, net, i, j, prev) = change_stat(t, net, i, j, prev)
_tchange(t::AbstractERGMTerm, net, i, j, prev) = change_stat(t, net, i, j)

_tcompute(t::TemporalTerm, net, prev) = compute(t, net, prev)
_tcompute(t::AbstractERGMTerm, net, prev) = compute(t, net)

"""
    EdgeStability <: TemporalTerm

The number of dyads whose state agrees with the previous network:
`Σ_{ij} 1[y_ij = y^{t-1}_ij]`. Its change statistic is `+1` for a dyad
that had a tie at t−1 and `-1` for one that did not, whatever the
current state.

**A descriptive, not a model term.** In a separable STERGM the change
statistic is constant on the free dyads of either model — `-1` on every
prior non-tie (formation: the column is `-edges`), `+1` on every prior tie
(persistence: the column is `+edges`) — so the term carries no information
a CMPLE could estimate; [`STERGMModel`](@ref) refuses it in a formula with
an `ArgumentError` saying so. It is informative only in a non-separable
(btergm `memory`-style) TERGM, which TERGM.jl does not fit. Use
`compute(EdgeStability(), curr, prev)` for reporting.

# Example
```julia
using TERGM, ERGM
t0 = network(5; directed=true)
for (i, j) in [(1, 2), (2, 1), (3, 4), (4, 5)]; add_edge!(t0, i, j); end
t1 = network(5; directed=true)
for (i, j) in [(1, 2), (3, 4), (2, 3), (5, 1)]; add_edge!(t1, i, j); end
compute(EdgeStability(), t1, t0)               # 16.0: 20 ordered dyads, 4 changed
change_stat(EdgeStability(), t1, 2, 1, t0)     # 1.0: 2→1 was a tie at t0
change_stat(EdgeStability(), t1, 2, 3, t0)     # -1.0: 2→3 was not
```
"""
struct EdgeStability <: TemporalTerm end

name(::EdgeStability) = "edge.stability"

function compute(::EdgeStability, net, prev_net)
    n = nv(net)
    directed = is_directed(net)
    total = 0.0
    for i in 1:n
        for j in (directed ? (1:n) : (i+1:n))
            i == j && continue
            if has_edge(net, i, j) == has_edge(prev_net, i, j)
                total += 1.0
            end
        end
    end
    return total
end

# Adding (i,j): agreement becomes 1 if prev had the edge, 0 otherwise;
# without the edge it is the reverse → Δ = ±1, independent of current state
change_stat(::EdgeStability, net, i::Int, j::Int, prev_net) =
    has_edge(prev_net, i, j) ? 1.0 : -1.0

"""
    Delrecip <: TemporalTerm

Delayed reciprocity: `Σ_{i≠j} y_ij · y^{t-1}_ji`, the number of ordered
dyads with `i→j` in the current network and `j→i` in the previous one
(label `delrecip`; change statistic `1[y^{t-1}_ji]`, independent of the
current state). **Directed networks only**: `compute` and `change_stat`
throw an `ArgumentError` on an undirected network (a silent `0.0` would be
a wrong, all-zero design column), and model construction refuses the term
on undirected panels through `requires_directed`.

**Provenance.** This is btergm's `delrecip` (Leifeld, Cranmer & Desmarais;
`?btergm::"tergm-terms"`, lag 1), applied here to the auxiliary networks of
a separable model. R **tergm** (≥ 4, the version the golden fixture is
generated with) has **no `delrecip` term** — `tergm(... Form(~edges +
delrecip) ...)` fails with `InitErgmTerm.delrecip not found` — so no R
fixture can pin its coefficient; its definition and change statistic are
pinned by hand in the test suite instead.

# Example
```julia
using TERGM, ERGM
t0 = network(5; directed=true)
for (i, j) in [(1, 2), (2, 1), (3, 4), (4, 5)]; add_edge!(t0, i, j); end
t1 = network(5; directed=true)
for (i, j) in [(1, 2), (3, 4), (2, 3), (5, 1)]; add_edge!(t1, i, j); end
compute(Delrecip(), t1, t0)                    # 1.0: 1→2 at t1 answers 2→1 at t0
change_stat(Delrecip(), t1, 2, 1, t0)          # 1.0: t0 has 1→2
change_stat(Delrecip(), t1, 3, 4, t0)          # 0.0: t0 has no 4→3
requires_directed(Delrecip())                  # true
u0 = network(3; directed=false); u1 = copy(u0)
try; compute(Delrecip(), u1, u0); catch e; e isa ArgumentError; end   # true
```
"""
struct Delrecip <: TemporalTerm end

name(::Delrecip) = "delrecip"

# The ONE message for a directed-only temporal term on an undirected
# network — the same shape as ERGM.jl's `_undirected_only_message`.
_directed_only_message(t) =
    "term '$(name(t))' is only defined for directed networks, but the " *
    "network is undirected. Remove the term or use directed panels."

# Guard called by the directed-only terms' own compute/change_stat. The
# branch is on the network's type parameter, so it costs nothing in the
# hot loops of a (directed) model.
@inline function _refuse_undirected(term, net)
    is_directed(net) || _throw_directed_only(term)
    return nothing
end
@noinline _throw_directed_only(term) = throw(ArgumentError(_directed_only_message(term)))

function compute(t::Delrecip, net, prev_net)
    _refuse_undirected(t, net)
    total = 0.0
    for e in edges(net)
        has_edge(prev_net, dst(e), src(e)) && (total += 1.0)
    end
    return total
end

function change_stat(t::Delrecip, net, i::Int, j::Int, prev_net)
    _refuse_undirected(t, net)
    return has_edge(prev_net, j, i) ? 1.0 : 0.0
end

"""
    PersistentEdge <: TemporalTerm

The number of edges present in both the current and the previous network
(label `persistent.edges`). Its change statistic is `1` for a dyad that was
a tie at t−1 and `0` otherwise.

**A descriptive, not a model term** (`gof` reports it as the `persisted`
count). On the free dyads of a separable STERGM its change statistic is
constant — `0` on every prior non-tie (formation: an all-zero column), `1`
on every prior tie (persistence: the `edges` column again) — so
[`STERGMModel`](@ref) refuses it in a formula with an `ArgumentError`; see
[`EdgeStability`](@ref).

# Example
```julia
using TERGM, ERGM
t0 = network(5; directed=true)
for (i, j) in [(1, 2), (2, 1), (3, 4), (4, 5)]; add_edge!(t0, i, j); end
t1 = network(5; directed=true)
for (i, j) in [(1, 2), (3, 4), (2, 3), (5, 1)]; add_edge!(t1, i, j); end
compute(PersistentEdge(), t1, t0)              # 2.0: 1→2 and 3→4 persist
change_stat(PersistentEdge(), t1, 1, 2, t0)    # 1.0
change_stat(PersistentEdge(), t1, 2, 3, t0)    # 0.0
```
"""
struct PersistentEdge <: TemporalTerm end

name(::PersistentEdge) = "persistent.edges"

function compute(::PersistentEdge, net, prev_net)
    total = 0.0
    for e in edges(net)
        has_edge(prev_net, src(e), dst(e)) && (total += 1.0)
    end
    return total
end

change_stat(::PersistentEdge, net, i::Int, j::Int, prev_net) =
    has_edge(prev_net, i, j) ? 1.0 : 0.0

"""
    NewEdge <: TemporalTerm

The number of edges present now but absent in the previous network
(label `new.edges`). Its change statistic is `1` for a dyad that was not a
tie at t−1 and `0` otherwise — the complement of [`PersistentEdge`](@ref).

**A descriptive, not a model term** (`gof` reports it as the `formed`
count). On the free dyads of a separable STERGM its change statistic is
constant — `1` on every prior non-tie (formation: the `edges` column
again), `0` on every prior tie (persistence: an all-zero column) — so
[`STERGMModel`](@ref) refuses it in a formula with an `ArgumentError`; see
[`EdgeStability`](@ref).

# Example
```julia
using TERGM, ERGM
t0 = network(5; directed=true)
for (i, j) in [(1, 2), (2, 1), (3, 4), (4, 5)]; add_edge!(t0, i, j); end
t1 = network(5; directed=true)
for (i, j) in [(1, 2), (3, 4), (2, 3), (5, 1)]; add_edge!(t1, i, j); end
compute(NewEdge(), t1, t0)                     # 2.0: 2→3 and 5→1 are new
change_stat(NewEdge(), t1, 2, 3, t0)           # 1.0
change_stat(NewEdge(), t1, 1, 2, t0)           # 0.0
compute(NewEdge(), t1, t0) + compute(PersistentEdge(), t1, t0) == ne(t1)   # true
```
"""
struct NewEdge <: TemporalTerm end

name(::NewEdge) = "new.edges"

function compute(::NewEdge, net, prev_net)
    total = 0.0
    for e in edges(net)
        has_edge(prev_net, src(e), dst(e)) || (total += 1.0)
    end
    return total
end

change_stat(::NewEdge, net, i::Int, j::Int, prev_net) =
    has_edge(prev_net, i, j) ? 0.0 : 1.0

# Dependence classification (extends ERGM.is_dyad_dependent, whose fallback
# is the conservative `true`): the built-in temporal terms condition only on
# the *previous* network, which is exogenous within a transition, so their
# change statistics do not depend on other dyads of the current network.
is_dyad_dependent(::EdgeStability) = false
is_dyad_dependent(::Delrecip) = false
is_dyad_dependent(::PersistentEdge) = false
is_dyad_dependent(::NewEdge) = false

# Constancy on the free dyads of a separable model. The formation model
# sees only prior NON-ties and the persistence model only prior ties, so a
# change statistic that depends on the dyad's previous state alone is the
# same number on every row of the model — ±edges or identically zero — and
# a CMPLE cannot estimate it (rank-deficient beside `Edges`, `edges` under
# another name without it). `nothing` for an informative term; otherwise the
# sentence `_validate_terms` refuses the formula with. Terms carrying
# information about the *previous* configuration of OTHER dyads (`Delrecip`)
# are fine: their column varies across rows.
_separable_constant(::AbstractERGMTerm, ::Symbol) = nothing
_separable_constant(::EdgeStability, side::Symbol) = side === :formation ?
    "-1 on every prior non-tie, i.e. the column is -edges" :
    "+1 on every prior tie, i.e. the column is +edges"
_separable_constant(::PersistentEdge, side::Symbol) = side === :formation ?
    "0 on every prior non-tie, i.e. an all-zero column" :
    "+1 on every prior tie, i.e. the column is edges under another name"
_separable_constant(::NewEdge, side::Symbol) = side === :formation ?
    "+1 on every prior non-tie, i.e. the column is edges under another name" :
    "0 on every prior tie, i.e. an all-zero column"

# Delayed reciprocity is only meaningful on directed networks: extend
# ERGM.jl's public direction trait so formula validation (ERGM's
# `validate_formula`, which TERGM runs on every panel) rejects it on
# undirected panels exactly as it rejects `Mutual`.
requires_directed(::Delrecip) = true

# =============================================================================
# Temporal descriptives
# =============================================================================

"""
    edge_ages(networks::Vector{<:Network}) -> Dict{Tuple{Int,Int}, Int}

For the last network of the sequence, the age of each current edge: the
number of consecutive panels (ending at the last one) in which it has
been present. Keys are `(src, dst)` pairs; an empty sequence is an
`ArgumentError`.

# Example
```julia
using TERGM, ERGM
t0 = network(5; directed=true)
for (i, j) in [(1, 2), (2, 1), (3, 4), (4, 5)]; add_edge!(t0, i, j); end
t1 = network(5; directed=true)
for (i, j) in [(1, 2), (3, 4), (2, 3), (5, 1)]; add_edge!(t1, i, j); end
ages = edge_ages([t0, t1])
ages[(1, 2)]     # 2: present at t0 and t1
ages[(2, 3)]     # 1: formed at t1
length(ages)     # 4: one entry per edge of the last panel
```
"""
function edge_ages(networks::Vector{<:Network})
    isempty(networks) && throw(ArgumentError("empty network sequence"))
    last_net = networks[end]
    ages = Dict{Tuple{Int, Int}, Int}()
    for e in edges(last_net)
        i, j = Int(src(e)), Int(dst(e))
        age = 1
        for t in (length(networks)-1):-1:1
            has_edge(networks[t], i, j) || break
            age += 1
        end
        ages[(i, j)] = age
    end
    return ages
end

"""
    mean_edge_age(networks::Vector{<:Network}) -> Float64

Mean [`edge_ages`](@ref) of the final panel's edges (NaN if edgeless).

# Example
```julia
using TERGM, ERGM
t0 = network(5; directed=true)
for (i, j) in [(1, 2), (2, 1), (3, 4), (4, 5)]; add_edge!(t0, i, j); end
t1 = network(5; directed=true)
for (i, j) in [(1, 2), (3, 4), (2, 3), (5, 1)]; add_edge!(t1, i, j); end
mean_edge_age([t0, t1])                        # 1.5: ages 2, 2, 1, 1
mean_edge_age([t0, network(5; directed=true)]) # NaN: the last panel is edgeless
```
"""
function mean_edge_age(networks::Vector{<:Network})
    ages = edge_ages(networks)
    return isempty(ages) ? NaN : mean(values(ages))
end

# =============================================================================
# Auxiliary networks (Krivitsky & Handcock construction)
# =============================================================================

# Attribute-preserving copy: delegates to `Base.copy(::Network)`, which
# duplicates the graph and all vertex/edge/network attributes, so nodal
# terms (NodeMatch, NodeCov, ...) in formation/dissolution formulas keep
# seeing covariates on the auxiliary networks.
_copy_net(net::Network) = copy(net)

"""
    formation_network(prev_net, curr_net) -> Network

The formation network Y⁺ = Y_{t-1} ∪ Y_t. Formation statistics are
evaluated on Y⁺; its free dyads are the non-edges of Y_{t-1}. Vertex,
edge and network attributes of `prev_net` are preserved (nodal terms in a
formation formula see their covariates).

# Example
```julia
using TERGM, ERGM
t0 = network(5; directed=true)
for (i, j) in [(1, 2), (2, 1), (3, 4), (4, 5)]; add_edge!(t0, i, j); end
t1 = network(5; directed=true)
for (i, j) in [(1, 2), (3, 4), (2, 3), (5, 1)]; add_edge!(t1, i, j); end
yplus = formation_network(t0, t1)
ne(yplus)                                      # 6: the 4 ties of t0 plus the 2 that formed
has_edge(yplus, 2, 1) && has_edge(yplus, 5, 1) # true: dissolved AND new ties are both in
```
"""
function formation_network(prev_net::Network, curr_net::Network)
    yplus = _copy_net(prev_net)
    for e in edges(curr_net)
        add_edge!(yplus, src(e), dst(e))
    end
    return yplus
end

"""
    dissolution_network(prev_net, curr_net) -> Network

The dissolution network Y⁻ = Y_{t-1} ∩ Y_t. Dissolution (persistence)
statistics are evaluated on Y⁻; its free dyads are the edges of Y_{t-1}.
Attributes of `prev_net` are preserved, as for [`formation_network`](@ref).

# Example
```julia
using TERGM, ERGM
t0 = network(5; directed=true)
for (i, j) in [(1, 2), (2, 1), (3, 4), (4, 5)]; add_edge!(t0, i, j); end
t1 = network(5; directed=true)
for (i, j) in [(1, 2), (3, 4), (2, 3), (5, 1)]; add_edge!(t1, i, j); end
yminus = dissolution_network(t0, t1)
ne(yminus)                                     # 2: only 1→2 and 3→4 persisted
has_edge(yminus, 2, 1)                         # false: 2→1 dissolved
has_edge(yminus, 2, 3)                         # false: a new tie is not in Y⁻
```
"""
function dissolution_network(prev_net::Network, curr_net::Network)
    yminus = _copy_net(prev_net)
    for e in edges(prev_net)
        if !has_edge(curr_net, src(e), dst(e))
            rem_edge!(yminus, src(e), dst(e))
        end
    end
    return yminus
end

# =============================================================================
# Model Types
# =============================================================================

"""
    STERGM

A separable temporal ERGM formula: formation terms and dissolution
(persistence) terms. Terms may be standard ERGM.jl terms (evaluated on
the auxiliary networks) or [`TemporalTerm`](@ref)s (also conditioned on
Y_{t-1}). Each side is any collection of terms `fit_ergm` would accept — a
`Vector{<:AbstractERGMTerm}`, a `push!`-built `Vector{Any}`, a bare term, a
vector nesting vectors of terms — normalised by ERGM.jl's term-list
collector, so a non-term element is refused with its position and type
(and `Edges` instead of `Edges()` gets the "did you mean" hint). Both
sides must be non-empty (`ArgumentError` naming the side otherwise); the
fields are `formation` and `dissolution`. Terms are validated against the
panels when the formula is bound to data by [`STERGMModel`](@ref).

# Example
```julia
using TERGM, ERGM
f = STERGM([Edges(), Delrecip()], [Edges(), Mutual()])
f.formation                                    # AbstractERGMTerm[Edges(), Delrecip()]
length(f.dissolution)                          # 2
f                                              # STERGM(formation: edges + delrecip; persistence: edges + mutual)
terms = []; push!(terms, Edges()); push!(terms, Mutual())
STERGM(terms, Edges()).formation == f.dissolution   # true: Vector{Any} and a bare term are fine
try; STERGM(AbstractERGMTerm[], [Edges()]); catch e; e isa ArgumentError; end   # true
```
"""
struct STERGM
    formation::Vector{AbstractERGMTerm}
    dissolution::Vector{AbstractERGMTerm}

    function STERGM(formation, dissolution)
        new(_side_terms(formation, "formation"),
            _side_terms(dissolution, "dissolution"))
    end
end

# One side's term list through ERGM's collector; an empty side is refused
# in TERGM's own words (naming the side) before the collector sees it.
function _side_terms(terms, side::String)
    terms isa AbstractVector && isempty(terms) &&
        throw(ArgumentError("the $side model must contain at least one term " *
                            "(e.g. `[Edges()]`); formation and dissolution " *
                            "models must both be non-empty"))
    terms isa Union{AbstractERGMTerm, AbstractVector} ||
        throw(ArgumentError("the $side model must be a term or a vector of " *
                            "terms (e.g. `[Edges(), Mutual()]`); got " *
                            "$(repr(terms)) of type $(typeof(terms))"))
    try
        return collect_terms(terms)
    catch e
        e isa ArgumentError || rethrow()
        throw(ArgumentError("$side model: " * e.msg))
    end
end

Base.show(io::IO, f::STERGM) =
    print(io, "STERGM(formation: ", join((name(t) for t in f.formation), " + "),
          "; persistence: ", join((name(t) for t in f.dissolution), " + "), ")")

# Formula validation at model construction is ERGM.jl's, not ours: the
# `ERGM.Extension.validate_formula` — the very function the `ERGMModel`
# constructor runs — is applied to EVERY panel, so a formation/dissolution
# formula is held to exactly the rules an ERGM formula is held to (a
# declared vertex/edge attribute must exist AND be set on every vertex —
# statnet refuses NA; `requires_directed` terms such as `Mutual`/`Delrecip`
# are refused on undirected panels; `requires_undirected` terms such as
# `Kstar`/`GWDegree`/`Degree` are refused on directed panels with the
# directed variant named; an `EdgeCov` matrix of the wrong size is refused).
# Only the message prefix is TERGM's: it names the side and the panel.
#
# One check IS TERGM's, because only a separable model has it: a term whose
# change statistic is constant on the free dyads of the side it sits in
# (`_separable_constant`) is refused, naming the mechanism.
# ERGM.jl's `Offset(term, coef)` fixes a coefficient; the CMPLE, the CMLE
# and the samplers here would treat it as a free one, so it is refused. R
# tergm refuses the same model in a conditional fit ("The ERGM `formula=`
# argument of an `N()` operator may not have offsets"): its offsets belong
# to EGMME, which is not implemented here.
function _refuse_offsets(terms, side::String)
    any(t -> t isa Offset, terms) && throw(ArgumentError(
        "$side model: offset terms (`Offset(term, coef)`, statnet's " *
        "`offset()`) are not implemented in TERGM.jl — the coefficient " *
        "would be estimated instead of held fixed (R tergm refuses offsets " *
        "inside Form()/Persist() in a CMLE/CMPLE fit too; it uses them with " *
        "EGMME). Remove the offset (see README 'Not implemented')."))
    return nothing
end

function _validate_terms(terms::Vector{AbstractERGMTerm}, networks, side::String)
    _refuse_offsets(terms, side)
    for term in terms
        why = _separable_constant(term, side == "formation" ? :formation : :dissolution)
        why === nothing && continue
        throw(ArgumentError(
            "$side model: term '$(name(term))' is constant on the free dyads of " *
            "the $side model of a separable STERGM (its change statistic is $why), " *
            "so the CMPLE design would be rank-deficient with `Edges()` and the " *
            "term would be `edges` under another name without it. The memory terms " *
            "EdgeStability, PersistentEdge and NewEdge are only informative in a " *
            "non-separable (btergm `memory`-style) TERGM, which TERGM.jl does not " *
            "fit; drop the term from the formula and use it as a descriptive " *
            "(`compute(term, curr, prev)`) instead."))
    end
    ts = TermSet(terms)
    for (t, net) in enumerate(networks)
        try
            validate_formula(ts, net)
        catch e
            e isa ArgumentError || rethrow()
            throw(ArgumentError("$side model, panel $t: " * e.msg))
        end
    end
    return nothing
end

# ----------------------------------------------------------------------------
# Term expansion and materialization — ERGM.jl's, not ours.
#
# `ERGMModel` runs `ERGM.Extension.materialize` on its formula: a
# multi-level `NodeFactor`/`NodeMix` and a multi-degree `Degree`/`IDegree`/
# `ODegree` expand into one statistic per level/cell/degree (R's
# `nodefactor.grp.b`, `mix.grp.a.b`, `degree1`), and every attribute term is
# replaced by a typed twin holding a dense snapshot of the attribute, so the
# hot loops never read the attribute Dicts. A STERGM needs exactly the same
# two things, with one twist: its terms are evaluated on T−1 auxiliary
# networks (each carrying the previous panel's attributes), so ONE snapshot
# taken at construction would be the first panel's attribute frozen onto
# every transition. The model therefore stores the *specification* — the
# expanded formula as plain single-statistic terms (`NodeFactor(:grp;
# level="b")`, `NodeMix(:grp, "a", "b")`, `Degree(1)`; every other term
# unchanged) — and snapshots per transition, against the auxiliary network
# the statistics are evaluated on (`_materialized_tuple`).
#
# The specification is what ERGM's `Extension.expand_terms(terms, net)`
# returns: one plain term per statistic, in ERGM's (R's) order, levels
# resolved against `net` — the same expansion `materialize` performs,
# without the snapshot. (It replaced TERGM's own `_specification`, which
# unwrapped the twin's `base` field from outside ERGM.)
# ----------------------------------------------------------------------------

# The materialized twins of an expanded side, as a tuple, snapshotted from
# `net` — the network the change statistics are evaluated on. Every term is
# single-statistic here, so the tuple has one entry per coefficient.
function _materialized_tuple(terms::AbstractVector{<:AbstractERGMTerm}, net)
    mts = materialize(TermSet(terms), net).terms
    length(mts) == length(terms) ||
        error("TERGM: materializing an expanded formula changed its length " *
              "($(length(terms)) terms → $(length(mts)) statistics); the model " *
              "must hold single-statistic terms only")
    return mts
end

# Expand one side against the panel. Levels/cells/degrees are resolved from
# the FIRST panel (as `ERGMModel` resolves them from its one network). R
# tergm resolves them from the union of all panels — a level absent from
# one wave simply has an all-zero column on that transition — and the two
# rules coincide exactly when no later panel carries a value the first one
# lacks, so that is what is checked: for every term that expands by
# attribute level (a multi-level NodeFactor, an unresolved NodeMix), the
# values of its attribute on every later panel must occur on panel 1.
# A level missing from a later panel is fine (its single-level term
# materializes there as an all-zero column, as in R); a NEW level would be
# a column R fits and this model does not have, so it is refused, naming
# the panel and the value.
function _expand_side(terms::Vector{AbstractERGMTerm}, networks, side::String)
    expanded = try
        expand_terms(terms, networks[1])
    catch e
        e isa ArgumentError || rethrow()
        throw(ArgumentError("$side model, panel 1: " * e.msg))
    end
    for term in terms
        materialize(term, networks[1]) isa AbstractVector || continue   # expands by level
        for attr in required_vertex_attributes(term)
            first_levels = Set(Any[v for v in values(get_vertex_attribute(networks[1], attr))])
            for t in 2:length(networks)
                new_levels = sort!(unique(Any[v for v in values(get_vertex_attribute(networks[t], attr))
                                            if !(v in first_levels)]); by=string)
                isempty(new_levels) || throw(ArgumentError(
                    "$side model, panel $t: term '$(name(term))' resolves its levels " *
                    "from the attribute's values on panel 1 " *
                    "($(join(sort!(collect(first_levels); by=string), ", "))), but " *
                    ":$attr takes the value$(length(new_levels) == 1 ? "" : "s") " *
                    "$(join(repr.(new_levels), ", ")) on panel $t, which panel 1 " *
                    "lacks — R tergm would give that level its own coefficient. " *
                    "Make the first panel carry every level the attribute takes " *
                    "on any panel (a level absent from a LATER panel is fine: its " *
                    "column is zero on that transition), or name the levels " *
                    "explicitly (`NodeFactor(:attr; levels=[...])`)."))
            end
        end
    end
    return expanded
end

# The formula with both sides expanded against `net` (single-statistic
# terms; idempotent on an already-expanded formula).
_expand_formula(formula::STERGM, net) =
    STERGM(expand_terms(formula.formation, net), expand_terms(formula.dissolution, net))

"""
    STERGMModel{T,D}

A STERGM formula bound to an observed network panel — a vector of
`Network{T,D}` (vertex type `T`, directedness `D`, both taken from the
panels, so `is_directed(model) == D` is a type-level fact and the
`networks` field is a concretely typed vector).

Construction (`STERGMModel(formula, networks)`) requires at least two
panels that share size, directedness and vertex type; refuses a **two-mode
(bipartite)** panel and a panel containing a **self-loop** (the one-mode
CMPLE, the free-dyad lists and the sampler range over the off-diagonal
one-mode dyads only, exactly as ERGM.jl refuses the same networks); refuses
any panel with masked (missing) dyads — CMPLE would enumerate them at face
value; refuses the memory terms `EdgeStability`/`PersistentEdge`/`NewEdge`
in either formula (constant on the free dyads of a separable model) and
ERGM.jl's `Offset` terms (not implemented here: the coefficient would be
estimated, not held fixed); and
validates both formulas against **every** panel with ERGM.jl's own formula
validation (`ERGM.Extension.validate_formula`): a declared vertex attribute must
exist and be set on every vertex, `Mutual`/`Delrecip` are refused on
undirected panels, `Kstar`/`GWDegree`/`Degree` on directed ones (use
`OStar`/`IStar`, `GWODegree`/`GWIDegree`, `ODegree`/`IDegree`), and an
`EdgeCov` matrix must be `n×n`. The `ArgumentError` names the side
(`formation`/`dissolution`), the panel, the term and the fix.

Both formulas are then **expanded** exactly as an `ERGMModel` expands
its terms (ERGM.jl's `materialize`): a multi-level `NodeFactor(:grp)`
becomes one `nodefactor.grp.<level>` statistic per non-base level, a
multi-cell `NodeMix(:grp)` one `mix.grp.<l1>.<l2>` per selected cell, and
`Degree(0:2)`/`IDegree`/`ODegree` one `degree<d>` per degree — one
coefficient per statistic, labelled as R `tergm` labels them. The stored
`formula` holds the expanded, single-statistic terms (`model.formula.formation`
is what `formation_coef` is indexed by). Levels, cells and degrees are
resolved from the first panel, as R `tergm` resolves them from the union
of the panels: a level absent from a later panel has an all-zero column
on that transition, and a later panel carrying a level the first one
lacks is refused, naming the panel and the value. The attribute values
themselves are read per transition from the panel the transition starts
at, so a time-varying attribute is honoured.

`show` prints one line — panels, size, directedness, free dyads and the
two formulas' coefficient labels — never the panels themselves.

# Example
```julia
using TERGM, ERGM
t0 = network(5; directed=true); add_edge!(t0, 1, 2)
t1 = network(5; directed=true); add_edge!(t1, 1, 2); add_edge!(t1, 2, 3)
model = STERGMModel(STERGM([Edges()], [Edges()]), [t0, t1])
is_directed(model)          # true
typeof(model)               # STERGMModel{Int64, true}
model                       # STERGMModel{Int64,true}: 2 panels of 5 vertices (directed), 20 free dyads; formation: edges; persistence: edges
bip = [network(6; bipartite=3), network(6; bipartite=3)]
try; STERGMModel(STERGM([Edges()], [Edges()]), bip); catch e; e isa ArgumentError; end   # true
```
"""
struct STERGMModel{T,D}
    formula::STERGM
    networks::Vector{Network{T,D}}

    function STERGMModel{T,D}(formula::STERGM,
                              networks::Vector{Network{T,D}}) where {T,D}
        length(networks) >= 2 ||
            throw(ArgumentError("need at least two panels (one transition)"))
        n = nv(networks[1])
        all(net -> nv(net) == n, networks) ||
            throw(ArgumentError("all panels must share size, directedness and " *
                                "vertex type (sizes differ)"))
        for (t, net) in enumerate(networks)
            # Two-mode panels: the one-mode CMPLE would enumerate every i≠j
            # dyad, including the within-mode dyads that can never form, as
            # free formation rows (ERGM.jl refuses the same network).
            is_two_mode(net) && throw(ArgumentError(
                "TERGM (panel $t): two-mode (bipartite) panels are not supported " *
                "— the one-mode CMPLE would enumerate the impossible within-mode " *
                "dyads as free dyads, biasing the formation model; bipartite " *
                "terms and free-dyad sets are not implemented (see README 'Not " *
                "implemented'; ERGM.jl refuses the same network)."))
            # Self-loops: the statistics would count the loop, but the free
            # dyads, `nobs` and the sampler range over off-diagonal dyads only,
            # so the two halves of the package would disagree about the data.
            loops = [v for v in 1:Int(n) if has_edge(net, v, v)]
            isempty(loops) || throw(ArgumentError(
                "TERGM (panel $t): the panel contains " *
                "$(length(loops)) self-loop$(length(loops) == 1 ? "" : "s") (at " *
                "vertex $(join(loops, ", "))); STERGM statistics, free dyads and " *
                "the sampler range over off-diagonal dyads only, so a loop could " *
                "neither form nor dissolve yet would be counted by every " *
                "statistic (ERGM.jl refuses the same network). Remove it with " *
                "`rem_edge!(net, v, v)` or build the panels with `loops=false`."))
            # Panel missingness is not implemented: CMPLE enumerates every
            # free dyad of Y⁺/Y⁻ as observed, so a masked dyad would silently
            # enter the design matrix at its face value. Reject rather than
            # invent data.
            require_observed(net; context="TERGM (panel $t)", face_ok=false)
        end
        _validate_terms(formula.formation, networks, "formation")
        _validate_terms(formula.dissolution, networks, "dissolution")
        # The expanded specification: one plain term per statistic (R's
        # per-level / per-cell / per-degree columns), snapshotted per
        # transition by the estimator, sampler and gof — see `_expand_side`.
        expanded = STERGM(_expand_side(formula.formation, networks, "formation"),
                          _expand_side(formula.dissolution, networks, "dissolution"))
        new{T,D}(expanded, networks)
    end
end

# One line, like `ERGMModel`'s: never the panels themselves (a
# fifty-actor eight-wave panel is eight networks). Labels are the
# direction-aware R labels the fit will print.
function Base.show(io::IO, m::STERGMModel{T,D}) where {T,D}
    net = m.networks[1]
    print(io, "STERGMModel{$T,$D}: $(length(m.networks)) panels of $(nv(net)) ",
          "vertices ", D ? "(directed)" : "(undirected)", ", ",
          _n_free_dyads(m), " free dyads; formation: ",
          join(_labels(m.formula.formation, net), " + "),
          "; persistence: ", join(_labels(m.formula.dissolution, net), " + "))
    return nothing
end

# `T` and `D` come from the first panel; every other panel must be the same
# `Network{T,D}` — the directedness check is a type check now.
function STERGMModel(formula::STERGM, networks::AbstractVector{<:Network})
    length(networks) >= 2 ||
        throw(ArgumentError("need at least two panels (one transition)"))
    N = typeof(networks[1])
    all(net -> net isa N, networks) ||
        throw(ArgumentError("all panels must share size, directedness and " *
                            "vertex type (got $(join(unique(map(typeof, networks)), ", ")))"))
    return _stergm_model(formula, collect(N, networks))
end
_stergm_model(formula::STERGM, networks::Vector{Network{T,D}}) where {T,D} =
    STERGMModel{T,D}(formula, networks)

Graphs.is_directed(::STERGMModel{T,D}) where {T,D} = D
Graphs.is_directed(::Type{<:STERGMModel{T,D}}) where {T,D} = D

"""
    STERGMResult{T,D}

Fitted STERGM. The dissolution model is fit in tergm's **persistence**
(`Persist()`) parameterisation: `persistence_coef` are persistence log-odds
(positive = ties last longer; tergm's `Diss()` coefficients are their
negation).

Fields beyond the coefficients and standard errors:

- `method` — `:cmple` or `:cmle`.
- `loglik_formation`, `loglik_dissolution` — the maximized conditional
  pseudo-log-likelihoods (`:cmple`), or the path-sampled conditional
  log-likelihoods (`:cmle`; `NaN` when `bridge_rungs = 0`).
- `vcov` — the joint covariance of the stacked (formation, then persistence)
  coefficients: block-diagonal by separability for `se_type == :hessian`
  and `:fisher`, the empirical covariance of the refits for the two
  bootstraps.
- `se_type` — how the standard errors were obtained: `:hessian` (inverse
  observed pseudo-likelihood information), `:bootstrap` (parametric
  bootstrap), `:block_bootstrap` (per-transition block bootstrap) or `:fisher`
  (a CMLE fit: inverse Fisher information from the MCMC sample plus the
  Monte-Carlo component).
- `boot_replicates` — the `n_boot × p` refitted coefficients of a bootstrap
  (rows with a non-finite refit are kept as `NaN` and excluded from the
  covariance; empty otherwise).
- `n_kept` — the number of free dyads the *finite* coefficients were
  estimated on: every free dyad of every transition, or, after a statistic
  at the boundary of its attainable range fixed a coefficient at `±Inf`
  (R ergm's `drop`), the dyads the dropped statistics do not touch. It is
  the BIC sample size, as R's `logLik` `nobs` attribute is (`nobs(r)`
  itself stays every free dyad).
- `converged` — `false` when Newton exhausted `maxiter`, the
  pseudo-likelihood has no finite maximum (perfect separation: the
  coefficients are then the last iterate, not an estimate), or the CMLE's
  MCMC convergence tests did not pass; all are warned about at fit time and
  listed by `approximations(r)`.
- `inference_withheld` — `true` for a CMPLE fit (default `se`) of a
  dyad-dependent formula: the naive pseudo-likelihood standard errors are
  reported, but no z value, p-value or confidence interval is built on them
  (see [`cmple`](@ref)).
- `mcmc` — `nothing` for a CMPLE fit; for a CMLE fit a NamedTuple
  `(formation, persistence)` of per-side reports (`iterations`,
  `step_length`, `termination`, `termination_p`, `t_ratios`,
  `hotelling_p`, `n_eff`, `n_samples`, `mcmc_se`, `exact` — `true` for a
  dyad-independent side, which is fit exactly without MCMC).

Accessors: [`formation_coef`](@ref), [`formation_se`](@ref),
[`persistence_coef`](@ref), [`persistence_se`](@ref), plus the
full StatsAPI surface on the stacked (formation, then persistence)
coefficients: `coef`, `stderror`, `vcov`, [`confint`](@ref), `loglikelihood`,
`aic`, `bic`, `nobs`, `dof` and [`coeftable`](@ref) (labels `Form(1)~…` /
`Persist(1)~…`, as `summary(tergm)` prints them).

# Example
```julia
using TERGM, ERGM
t0 = network(5; directed=true)
for (i, j) in [(1, 2), (2, 1), (3, 4), (4, 5)]; add_edge!(t0, i, j); end
t1 = network(5; directed=true)
for (i, j) in [(1, 2), (3, 4), (2, 3), (5, 1)]; add_edge!(t1, i, j); end
fit = stergm([t0, t1], [Edges()], [Edges()])   # a STERGMResult{Int64, true}
fit.method, fit.se_type, fit.converged         # (:cmple, :hessian, true)
fit.n_kept == nobs(fit) == 20                  # true: 16 formation + 4 dissolution rows
size(fit.boot_replicates)                      # (0, 2): no bootstrap was run
fit.inference_withheld                         # false: an edges-only CMPLE is exact
approximations(fit)                            # String[]
```
"""
struct STERGMResult{T,D}
    model::STERGMModel{T,D}
    formation_coef::Vector{Float64}
    formation_se::Vector{Float64}
    persistence_coef::Vector{Float64}
    persistence_se::Vector{Float64}
    vcov::Matrix{Float64}
    loglik_formation::Float64
    loglik_dissolution::Float64
    converged::Bool
    method::Symbol
    se_type::Symbol
    boot_replicates::Matrix{Float64}
    n_kept::Int
    inference_withheld::Bool
    mcmc::Union{Nothing,NamedTuple}
end

"""
    formation_coef(r::STERGMResult) -> Vector{Float64}

The formation model's coefficients: log-odds contributions to a prior
non-tie forming, one per formation term (`r.model.formula.formation`).

# Example
```julia
using TERGM, ERGM
t0 = network(5; directed=true); add_edge!(t0, 1, 2); add_edge!(t0, 3, 4)
t1 = network(5; directed=true); add_edge!(t1, 1, 2); add_edge!(t1, 2, 3)
fit = stergm([t0, t1], [Edges()], [Edges()])
formation_coef(fit)[1] ≈ log(1 / 17)   # 1 of 18 prior non-ties formed: logit(1/18)
approximations(fit) == String[]        # true: both sides have a finite CMPLE
```
"""
formation_coef(r::STERGMResult) = r.formation_coef

"""
    formation_se(r::STERGMResult) -> Vector{Float64}

Standard errors of [`formation_coef`](@ref); `se_method(r)` says how they
were obtained (`:hessian`, `:bootstrap`, `:block_bootstrap`, or `:fisher` for a
CMLE fit).

# Example
```julia
using TERGM, ERGM
t0 = network(5; directed=true); add_edge!(t0, 1, 2); add_edge!(t0, 3, 4)
t1 = network(5; directed=true); add_edge!(t1, 1, 2); add_edge!(t1, 2, 3)
fit = stergm([t0, t1], [Edges()], [Edges()])
formation_se(fit) == stderror(fit)[1:1]   # true
formation_se(fit)[1] ≈ sqrt(1 / (18 * (1/18) * (17/18)))   # true: Bernoulli information of 18 rows
```
"""
formation_se(r::STERGMResult) = r.formation_se

"""
    persistence_coef(r::STERGMResult) -> Vector{Float64}

The dissolution model's coefficients in tergm's **persistence**
(`Persist()`) parameterisation: log-odds contributions to a prior tie
persisting, one per dissolution term (`r.model.formula.dissolution`).
Positive values mean ties last longer; tergm's `Diss()` coefficients are
`-persistence_coef(r)`.

# Example
```julia
using TERGM, ERGM
t0 = network(5; directed=true); add_edge!(t0, 1, 2); add_edge!(t0, 3, 4)
t1 = network(5; directed=true); add_edge!(t1, 1, 2); add_edge!(t1, 2, 3)
fit = stergm([t0, t1], [Edges()], [Edges()])
persistence_coef(fit)[1] ≈ 0.0   # 1 of 2 prior ties persisted: logit(1/2)
approximations(fit) == String[]  # true: both sides have a finite CMPLE
```
"""
persistence_coef(r::STERGMResult) = r.persistence_coef

"""
    persistence_se(r::STERGMResult) -> Vector{Float64}

Standard errors of [`persistence_coef`](@ref); `se_method(r)` says how they
were obtained (`:hessian`, `:bootstrap`, `:block_bootstrap`, or `:fisher` for a
CMLE fit).

# Example
```julia
using TERGM, ERGM
t0 = network(5; directed=true); add_edge!(t0, 1, 2); add_edge!(t0, 3, 4)
t1 = network(5; directed=true); add_edge!(t1, 1, 2); add_edge!(t1, 2, 3)
fit = stergm([t0, t1], [Edges()], [Edges()])
persistence_se(fit) == stderror(fit)[2:2]   # true
persistence_se(fit)[1] ≈ sqrt(2)             # true: 2 Bernoulli rows at p = 1/2
```
"""
persistence_se(r::STERGMResult) = r.persistence_se

"""
    has_dyad_dependent(model::STERGMModel) -> Bool

Method on ERGM.jl's exported predicate: whether **either** formula
(formation or dissolution) contains a dyad-dependent term
(`is_dyad_dependent`; the temporal terms are dyad-independent — they
condition only on the exogenous previous network).

This is THE predicate that decides whether CMPLE is an approximation: when
both formulas are dyad-independent the conditional pseudo-likelihood *is* the
conditional likelihood, so CMPLE is the exact CMLE and the inverse-Hessian
standard errors are the exact ML ones. Used by `show(::STERGMResult)`
(the prose caveat), by `is_exact(::STERGMResult)` (the machine-readable
answer) and by the rule that withholds naive inference, so the three cannot
drift apart.

# Example
```julia
using TERGM, ERGM
t0 = network(5; directed=true); add_edge!(t0, 1, 2)
t1 = network(5; directed=true); add_edge!(t1, 1, 2); add_edge!(t1, 2, 3)
has_dyad_dependent(STERGMModel(STERGM([Edges(), Delrecip()], [Edges()]), [t0, t1]))  # false
has_dyad_dependent(STERGMModel(STERGM([Edges()], [Edges(), Mutual()]), [t0, t1]))    # true
```
"""
has_dyad_dependent(model::STERGMModel) =
    any(is_dyad_dependent(t)
        for t in Iterators.flatten((model.formula.formation,
                                    model.formula.dissolution)))

# Coefficient labels resolved against a panel with the direction-aware
# two-argument `name(term, net)` (ERGM 0.2): a directed panel's `GWESP(0.5)`
# is `gwesp.OTP.fixed.0.5`, as R prints it, an undirected one's
# `gwesp.fixed.0.5`. Used by `show` and by `coeftable`.
_labels(terms, net) = [name(t, net) for t in terms]

# The two coefficient tables, built ONCE from the fit's own vectors —
# labels resolved against the panel, z = θ/se, p through `_pvalues` (a
# coefficient fixed at ∓Inf by a boundary statistic prints p = 0, as R
# does). `show` prints them as the `Formation:` / `Persistence:` blocks and
# `coeftable` stacks them under tergm's `Form(1)~` / `Persist(1)~` prefixes,
# so what is shown and what is inspected are the same numbers by
# construction, not by discipline. When the naive inference is withheld
# (`inference_withheld`), z and p of every finite coefficient are NaN.
function _side_tables(r::STERGMResult)
    net = r.model.networks[1]
    fz = r.formation_coef ./ r.formation_se
    dz = r.persistence_coef ./ r.persistence_se
    fp = _pvalues(r.formation_coef, fz)
    dp = _pvalues(r.persistence_coef, dz)
    # the ecosystem's separation policy: an unconverged CMPLE (separated, or
    # Newton cap exhausted) has no estimate to test — z and p are withheld
    if r.inference_withheld || _cmple_unconverged(r)
        for (θ, z, p) in ((r.formation_coef, fz, fp), (r.persistence_coef, dz, dp))
            for k in eachindex(θ)
                isfinite(θ[k]) || continue
                z[k] = NaN
                p[k] = NaN
            end
        end
    end
    form = CoefficientTable(_labels(r.model.formula.formation, net),
                            r.formation_coef, r.formation_se;
                            z_values=fz, p_values=fp)
    pers = CoefficientTable(_labels(r.model.formula.dissolution, net),
                            r.persistence_coef, r.persistence_se;
                            z_values=dz, p_values=dp)
    return form, pers
end

# One R-style block through the shared presentation layer
# (NetworkCore.print_coeftable), from a CoefficientTable's own vectors.
_print_block(io::IO, tbl::CoefficientTable; legend::Bool) =
    print_coeftable(io, tbl.names, tbl.estimates, tbl.std_errors, tbl.p_values;
                    z_values=tbl.z_values, legend=legend)

# How the standard errors were obtained, in words (the `Standard errors:`
# header line of `show`)
function _se_description(r::STERGMResult)
    B = size(r.boot_replicates, 1)
    r.se_type === :bootstrap && return "parametric bootstrap ($B replicates)"
    r.se_type === :block_bootstrap && return "per-transition block bootstrap ($B replicates)"
    r.se_type === :fisher && return _used_mcmc(r) ?
        "inverse Fisher information from the MCMC sample, plus the Monte-Carlo error" :
        "inverse Fisher information (exact: the formula is dyad-independent)"
    return "inverse Hessian of the pseudo-likelihood"
end

# The stopping rule's verdict on one side of a CMLE fit, in words
function _termination_detail(rep)
    p = isnan(rep.termination_p) ? "NaN" : string(round(rep.termination_p, sigdigits=3))
    return rep.termination === :confidence ?
        "$(round(Int, 100 * rep.conv_confidence))% equivalence test p $p (needs < " *
        "$(round(1 - rep.conv_confidence, sigdigits=3)))" :
        "Hotelling p $p, max t-ratio $(round(maximum(rep.t_ratios), digits=3))"
end

# Whether a CMLE fit drew an MCMC sample on at least one side (a
# dyad-independent side is fit exactly, without one)
_used_mcmc(r::STERGMResult) =
    r.mcmc !== nothing && !(r.mcmc.formation.exact && r.mcmc.persistence.exact)

function Base.show(io::IO, r::STERGMResult)
    n_panels = length(r.model.networks)
    println(io, "STERGM Results ($(r.method))")
    println(io, "="^40)
    println(io, "Panels: $n_panels ($(n_panels - 1) " *
                "$(n_panels == 2 ? "transition" : "transitions"), " *
                "$(nobs(r)) free dyads)")
    if r.method === :cmle
        if isnan(loglikelihood(r))
            println(io, "Log-likelihood: not estimated (bridge_rungs = 0)")
        else
            println(io, "Log-likelihood", _used_mcmc(r) ? " (path-sampled)" : "",
                    ": formation $(round(r.loglik_formation, digits=3)), " *
                    "persistence $(round(r.loglik_dissolution, digits=3))")
        end
    else
        println(io, "Pseudo-log-likelihood: formation $(round(r.loglik_formation, digits=3)), " *
                    "persistence $(round(r.loglik_dissolution, digits=3))")
    end
    # The same AIC/BIC line as `ERGMResult`'s header (pseudo-likelihood
    # based for CMPLE: heuristics under dyadic dependence, exact for
    # dyad-independent formulas where CMPLE is the CMLE)
    println(io, "AIC: $(round(aic(r), digits=2)), BIC: $(round(bic(r), digits=2))")
    println(io, "Standard errors: ", _se_description(r))
    # The verdict sits in the header, and an unconverged fit says so right
    # there — it must never look like a fit with a footnote (the same
    # sentence `approximations` reports, so the two cannot drift apart).
    println(io, "Converged: $(r.converged)")
    if !r.converged
        println(io, "  WARNING: ", _nonconvergence_caveat(r))
    end
    if _used_mcmc(r)
        for (side, rep) in (("formation", r.mcmc.formation), ("persistence", r.mcmc.persistence))
            rep.exact && continue
            println(io, "MCMC ($side): $(rep.iterations) iteration",
                    rep.iterations == 1 ? "" : "s", ", $(rep.n_samples) draws ",
                    "(effective sample size $(round(Int, rep.n_eff))); ",
                    _termination_detail(rep))
        end
    end
    println(io)

    # Shared ecosystem presentation layer: the two blocks ARE the two
    # halves of `coeftable(r)` (a NetworkCore.CoefficientTable each), rendered
    # through `print_coeftable` with the significance-code legend printed
    # once, below the second table.
    form, pers = _side_tables(r)
    println(io, "Formation:")
    _print_block(io, form; legend=false)
    println(io)
    println(io, "Persistence:")
    _print_block(io, pers; legend=true)

    # A coefficient fixed at ∓Inf by a boundary statistic is said under the
    # tables, from the same predicate `approximations` reports it by.
    fixed = _fixed_coefficient_note(r)
    if fixed !== nothing
        println(io)
        println(io, "Note: ", fixed)
    end

    # Honest-uncertainty caveat (the pattern of ERGM.jl's `show`): CMPLE
    # fits of dyad-dependent formulas have naive inverse-Hessian standard
    # errors. Dyad-independent formulas need no caveat — there CMPLE is the
    # CMLE.
    if r.method === :cmple && has_dyad_dependent(r.model)
        println(io)
        if r.se_type === :bootstrap
            println(io, "Note: this model contains dyad-dependent terms and was fit by")
            println(io, "conditional maximum pseudolikelihood (CMPLE). Standard errors are")
            println(io, "parametric-bootstrap estimates; the CMPLE point estimates may still")
            println(io, "be biased. Consider refitting with method=:cmle.")
        elseif r.se_type === :block_bootstrap
            println(io, "Note: this model contains dyad-dependent terms and was fit by")
            println(io, "conditional maximum pseudolikelihood (CMPLE). Standard errors are")
            println(io, "per-transition block-bootstrap estimates, which under-cover on short")
            println(io, "panels; the CMPLE point estimates may still be biased. Prefer")
            println(io, "se=:bootstrap or method=:cmle.")
        elseif r.inference_withheld
            println(io, "Note: z values and p-values are not reported (NaN). This model contains")
            println(io, "dyad-dependent terms and was fit by conditional maximum pseudolikelihood")
            println(io, "(CMPLE); the standard errors shown are the naive pseudolikelihood ones,")
            println(io, "which treat dependent dyads as independent (95% Wald intervals covered")
            println(io, "0.70-0.82 for GWESP and 0.89-0.96 for Mutual in simulation), so no test")
            println(io, "or interval is built on them. For inference refit with se=:bootstrap")
            println(io, "(parametric bootstrap) or method=:cmle; se=:hessian requests the naive")
            println(io, "Wald table explicitly.")
        else
            println(io, "Warning: this model contains dyad-dependent terms and was fit by")
            println(io, "conditional maximum pseudolikelihood (CMPLE). The standard errors are")
            println(io, "based on the naive pseudolikelihood and can be anticonservative; the")
            println(io, "p-values should not be trusted. Refit with method=:cmle, or use")
            println(io, "se=:bootstrap for parametric-bootstrap standard errors.")
        end
    end
    n_dropped = _n_dropped_replicates(r)
    if n_dropped > 0
        println(io)
        println(io, "Note: $n_dropped of $(size(r.boot_replicates, 1)) bootstrap " *
                    "refits had no finite CMPLE and were excluded from the standard errors. " *
                    _BOOT_EXCLUSION_BIAS)
    end
end

# THE sentence every bootstrap caller in the ERGM family uses to disclose what
# excluding failed refits does to the standard errors (the warning, `show`
# and `approximations` all carry it)
const _BOOT_EXCLUSION_BIAS =
    "The standard errors are conditional on a finite refit: the excluded " *
    "replicates are the extreme ones, so the standard errors are biased downward."

# Bootstrap replicates (rows of `boot_replicates`) that were excluded from
# the covariance because a refit had no finite coefficient vector.
_n_dropped_replicates(r::STERGMResult) =
    count(b -> !all(isfinite, view(r.boot_replicates, b, :)),
          1:size(r.boot_replicates, 1))

# p-values from z, with a coefficient fixed at ∓Inf by a boundary statistic
# (SE 0, z = ∓Inf) printed as p = 0, as R does.
function _pvalues(θ::AbstractVector, z::AbstractVector)
    p = z_pvalues(z)
    for k in eachindex(θ)
        isinf(θ[k]) && (p[k] = 0.0)
    end
    return p
end

# Coefficient labels of the stacked (formation, persistence) vector as tergm
# spells them (`Form(1)~edges`, `Persist(1)~nodematch.grp`, minus the lag
# index): what the fit-time warnings and the notes name a statistic by.
function _stacked_labels(model::STERGMModel)
    net = model.networks[1]
    return vcat("Form~" .* _labels(model.formula.formation, net),
                "Persist~" .* _labels(model.formula.dissolution, net))
end

# The one sentence for `converged == false`, shared by `show` and
# `approximations` so the printed and the machine-readable caveat agree.
# A CMPLE that did not converge — a separated design (the shared
# `NetworkCore.logistic_separation` verdict, through ERGM's design fitter) or
# an exhausted Newton cap — is a point on an asymptote, not an estimate: as
# across the ecosystem, z values, p-values and intervals are withheld (NaN)
_cmple_unconverged(r::STERGMResult) = r.method === :cmple && !r.converged

_nonconvergence_caveat(r::STERGMResult) = r.method === :cmle ?
    "the Monte-Carlo conditional MLE did not pass its convergence test " *
    "within maxiter iterations (see `fit.mcmc` for the stopping rule's " *
    "verdict and the t-ratios of the final MCMC sample): the coefficients " *
    "are the last iterate and the standard " *
    "errors are unreliable; raise maxiter or n_samples, or continue from " *
    "these coefficients with `init = coef(fit)`." :
    "the conditional maximum pseudo-likelihood estimate does not exist or was " *
    "not reached: the Newton iteration did not converge — either maxiter was " *
    "exhausted or the pseudo-likelihood has no finite maximum (perfect " *
    "separation; R tergm warns \"The MPLE does not exist!\" for the same " *
    "design). The coefficients are the last iterate, not an estimate, and " *
    "the standard errors are meaningless; z values, p-values and confidence " *
    "intervals are withheld (NaN). Remove or coarsen the term whose " *
    "statistic perfectly predicts the ties, or raise maxiter."

# A coefficient fixed at ∓Inf by a boundary statistic (R ergm's `drop`),
# read off the coefficients themselves; `nothing` when every coefficient is
# finite.
function _fixed_coefficient_note(r::STERGMResult)
    labels = _stacked_labels(r.model)
    θ = coef(r)
    lo = [labels[k] for k in eachindex(labels) if θ[k] == -Inf]
    hi = [labels[k] for k in eachindex(labels) if θ[k] == Inf]
    # A NaN coefficient beside finite ones on its side: a statistic that does
    # not vary on the free dyads, or a linear combination of the ones before
    # it (an all-NaN side is the empty-side case, said by the caveat)
    pf = length(r.formation_coef)
    side_nan(k) = isnan(θ[k]) &&
        !all(isnan, k <= pf ? r.formation_coef : r.persistence_coef)
    na = [labels[k] for k in eachindex(labels) if side_nan(k)]
    isempty(lo) && isempty(hi) && isempty(na) && return nothing
    notes = String[]
    if !(isempty(lo) && isempty(hi))
        parts = String[]
        isempty(lo) || push!(parts, "$(join(lo, ", ")) fixed at -Inf (observed " *
                                    "statistic at its smallest attainable value)")
        isempty(hi) || push!(parts, "$(join(hi, ", ")) fixed at +Inf (observed " *
                                    "statistic at its largest attainable value)")
        push!(notes, "coefficient(s) " * join(parts, "; ") * ": no finite estimate " *
                     "exists; the other coefficients are estimated on the free dyads " *
                     "these statistics do not touch, as R ergm does (drop=TRUE), with " *
                     "standard error 0 and p-value 0 recorded for the fixed ones")
    end
    isempty(na) || push!(notes, "coefficient(s) $(join(na, ", ")) not identifiable " *
                                "(the statistic does not vary on the free dyads, or is a " *
                                "linear combination of the statistics before it): " *
                                "reported as NaN, as R tergm's CMPLE reports NA")
    return join(notes, "; ")
end

# ============================================================================
# The shared result-metadata protocol (NetworkCore.jl `src/results.jl`)
# ============================================================================
#
# `fit_metadata(fit)` collects these accessors. They are derived from the SAME
# `has_dyad_dependent` predicate the prose caveat in `show` uses, so the
# machine-readable answer and the printed one can never disagree.

estimand(::STERGMResult) = :stergm

"""
    objective(r::STERGMResult) -> Symbol

`:conditional_pseudolikelihood` — the pooled logistic pseudo-likelihood over the
free dyads of the formation (Y⁺) and dissolution (Y⁻) auxiliary networks — for
a CMPLE fit and for a CMLE fit of a dyad-independent formula (where the two
coincide and no MCMC is run); `:mc_likelihood` for a CMLE fit that maximized
a Monte-Carlo approximation of the conditional likelihood.
"""
objective(r::STERGMResult) =
    _used_mcmc(r) ? :mc_likelihood : :conditional_pseudolikelihood

"""
    is_exact(r::STERGMResult) -> Bool

`true` iff the fit **converged**, every coefficient is **finite**, and
**every** term in **both** formulas is dyad-independent: there the
conditional pseudo-likelihood is the conditional likelihood, so CMPLE is the
exact CMLE. Add one dyad-dependent term to either formula and CMPLE reports
`false` — the approximation is a property of the fit, not of the method
name — and so does CMLE, which then maximizes a Monte-Carlo approximation
of the likelihood. A fit whose pseudo-likelihood has no finite maximum
(`converged == false`, or a coefficient fixed at `±Inf` by a boundary
statistic) is not an exact estimate of anything and reports `false` too.
"""
is_exact(r::STERGMResult) =
    r.method in (:cmple, :cmle) && r.converged && all(isfinite, coef(r)) &&
    !has_dyad_dependent(r.model)

se_method(r::STERGMResult) = r.se_type

# `STERGMModel` calls `require_observed` on every panel with the default
# `:error` policy: CMPLE would enumerate a masked dyad at its face value, so
# masked panels are refused outright rather than silently used.
missing_method(::STERGMResult) = :rejected

function approximations(r::STERGMResult)
    out = String[]
    if r.method === :cmple && has_dyad_dependent(r.model)
        push!(out, "conditional maximum pseudo-likelihood of a dyad-dependent " *
                   "formula: the dyad conditionals of each transition are " *
                   "multiplied as if independent, so the point estimates are " *
                   "biased in finite samples")
        r.se_type === :hessian &&
            push!(out, "inverse-Hessian standard errors of the naive " *
                       "conditional pseudo-likelihood: expected anticonservative " *
                       "under dyadic dependence")
        r.inference_withheld &&
            push!(out, "z values, p-values and confidence intervals withheld: " *
                       "the naive pseudo-likelihood standard errors are not " *
                       "calibrated under dyadic dependence (refit with " *
                       "se=:bootstrap or method=:cmle; se=:hessian opts in)")
    end
    r.se_type === :block_bootstrap &&
        push!(out, "per-transition block-bootstrap standard errors from " *
                   "$(length(r.model.networks) - 1) transitions: the resampling " *
                   "units are the transitions, so the standard errors under-cover " *
                   "on short panels (95% intervals covered " *
                   "$(_block_coverage_note(length(r.model.networks) - 1)) in simulation)")
    if _used_mcmc(r)
        push!(out, "CMLE: the conditional likelihood is approximated by an MCMC " *
                   "sample, so the estimates carry Monte-Carlo error (included in " *
                   "the standard errors; see `fit.mcmc`)")
        push!(out, isnan(loglikelihood(r)) ?
                   "log-likelihood not estimated (bridge_rungs=0): AIC/BIC are NaN" :
                   "the reported log-likelihood (and AIC/BIC) is a path-sampling " *
                   "bridge estimate from a dyad-independent reference model")
    end
    # Non-convergence and a coefficient fixed at ∓Inf are part of what the
    # fit actually did, so they are reported here as well as warned about at
    # fit time (never only in a log line).
    r.converged || push!(out, _nonconvergence_caveat(r))
    fixed = _fixed_coefficient_note(r)
    fixed === nothing || push!(out, fixed)
    n_dropped = _n_dropped_replicates(r)
    n_dropped > 0 &&
        push!(out, "bootstrap standard errors from $(size(r.boot_replicates, 1) - n_dropped) " *
                   "of $(size(r.boot_replicates, 1)) replicates: $n_dropped refits had no " *
                   "finite CMPLE and were excluded. " * _BOOT_EXCLUSION_BIAS)
    return out
end

# StatsAPI interface: methods on the shared statistics generics (mirroring
# ERGM.jl), so fitted STERGMs interoperate with StatsBase/GLM-style tooling.
# Coefficient-shaped accessors stack formation first, then persistence.

# Total number of free dyads pooled over transitions: every off-diagonal
# dyad of each transition is free in exactly one of the two models
function _n_free_dyads(model::STERGMModel)
    n = Int(nv(model.networks[1]))
    per = is_directed(model) ? n * (n - 1) : n * (n - 1) ÷ 2
    return (length(model.networks) - 1) * per
end

StatsAPI.coef(r::STERGMResult) = vcat(r.formation_coef, r.persistence_coef)
StatsAPI.stderror(r::STERGMResult) = vcat(r.formation_se, r.persistence_se)
StatsAPI.vcov(r::STERGMResult) = r.vcov
StatsAPI.loglikelihood(r::STERGMResult) = r.loglik_formation + r.loglik_dissolution
StatsAPI.aic(r::STERGMResult) = -2 * loglikelihood(r) + 2 * dof(r)
# The BIC sample size is the number of free dyads the finite coefficients
# were estimated on (`n_kept`: every free dyad, or, after R's drop of a
# boundary statistic, the dyads it does not touch — R's `logLik` nobs
# attribute); `nobs` itself stays every free dyad of every transition.
StatsAPI.bic(r::STERGMResult) = -2 * loglikelihood(r) + dof(r) * log(r.n_kept)
StatsAPI.nobs(r::STERGMResult) = _n_free_dyads(r.model)
# R's `logLik.ergm` df: a coefficient fixed at ∓Inf by a boundary statistic
# (R's drop) is not an estimated parameter
StatsAPI.dof(r::STERGMResult) = count(isfinite, coef(r))

"""
    coeftable(r::STERGMResult) -> NetworkCore.CoefficientTable

The R-style coefficient table (`Estimate`, `Std.Error`, `z value`,
`Pr(>|z|)`) of the stacked coefficients — formation first, then persistence
— as an inspectable `NetworkCore.CoefficientTable` (a method of
`StatsAPI.coeftable`, the ecosystem's ONE `coeftable` binding). Rows are
labelled as `summary(tergm)` labels them: `Form(1)~<term>` and
`Persist(1)~<term>`, with `<term>` ERGM.jl's direction-aware R label
(`gwesp.OTP.fixed.0.5` on a directed panel). The numbers are exactly the
ones `show(r)` prints in its `Formation:` / `Persistence:` blocks — both are
built from the same two tables — and a coefficient fixed at `±Inf` by a
boundary statistic carries `p = 0`, as R prints it. Rows can be read by
index or by name.

# Example
```julia
using TERGM, ERGM, NetworkCore
nets = [network(6; directed=true) for _ in 1:3]
for (t, (i, j)) in enumerate([(1, 2), (2, 3), (3, 4)]); add_edge!(nets[t], 1, 2); add_edge!(nets[t], i, j); end
fit = stergm(nets, [Edges()], [Edges()])
tbl = coeftable(fit)
tbl.names                                   # ["Form(1)~edges", "Persist(1)~edges"]
tbl["Persist(1)~edges"].estimate == persistence_coef(fit)[1]   # true
tbl[1].std_error == formation_se(fit)[1]    # true
```
"""
function StatsAPI.coeftable(r::STERGMResult)
    form, pers = _side_tables(r)
    return CoefficientTable(vcat("Form(1)~" .* form.names, "Persist(1)~" .* pers.names),
                            vcat(form.estimates, pers.estimates),
                            vcat(form.std_errors, pers.std_errors);
                            z_values=vcat(form.z_values, pers.z_values),
                            p_values=vcat(form.p_values, pers.p_values))
end

"""
    coefnames(r::STERGMResult) -> Vector{String}

The coefficient labels in `coef(r)` order — formation first, then
persistence — as `summary(tergm)` and R's `names(coef(fit))` give them:
`Form(1)~<term>` and `Persist(1)~<term>`. Identical to
`coeftable(r).names` (a method of `StatsAPI.coefnames`); a fresh vector.

# Example
```julia
using TERGM, ERGM, NetworkCore
nets = [network(6; directed=true) for _ in 1:3]
for (t, (i, j)) in enumerate([(1, 2), (2, 3), (3, 4)]); add_edge!(nets[t], 1, 2); add_edge!(nets[t], i, j); end
fit = stergm(nets, [Edges()], [Edges()])
coefnames(fit)                              # ["Form(1)~edges", "Persist(1)~edges"]
coefnames(fit) == coeftable(fit).names      # true
```
"""
function StatsAPI.coefnames(r::STERGMResult)
    net = r.model.networks[1]
    return vcat("Form(1)~" .* _labels(r.model.formula.formation, net),
                "Persist(1)~" .* _labels(r.model.formula.dissolution, net))
end

"""
    confint(r::STERGMResult; level=0.95) -> Matrix{Float64}

Normal-theory (Wald) confidence limits `θ̂ ± z_{(1+level)/2} · se`, one row
per stacked coefficient (formation first, then persistence; the order of
`coef(r)`), lower limit in column 1 and upper in column 2 (a method of
`StatsAPI.confint`). The standard errors are the ones the fit reports
(`r.se_type` says which). A CMPLE fit (default `se`) of a **dyad-dependent**
formula reports no interval — `confint` throws an `ArgumentError` naming
the ways to get one (`se = :bootstrap`, `method = :cmle`, or the explicit
opt-in `se = :hessian`, whose intervals inherit the naive pseudo-likelihood
standard errors; see [`cmple`](@ref)). A coefficient fixed at `±Inf` has a
degenerate interval `[±Inf, ±Inf]`, as it has standard error 0. `level`
must lie in `(0, 1)` (`ArgumentError` otherwise).

# Example
```julia
using TERGM, ERGM, NetworkCore
nets = [network(6; directed=true) for _ in 1:3]
for (t, (i, j)) in enumerate([(1, 2), (2, 3), (3, 4)]); add_edge!(nets[t], 1, 2); add_edge!(nets[t], i, j); end
fit = stergm(nets, [Edges()], [Edges()])
ci = confint(fit)                          # 2×2
all(ci[:, 1] .< coef(fit) .< ci[:, 2])     # true
confint(fit; level=0.9)                    # narrower
```
"""
function StatsAPI.confint(r::STERGMResult; level::Real=0.95)
    0 < level < 1 ||
        throw(ArgumentError("confint: level must be in (0, 1) (got $level)"))
    r.inference_withheld && throw(ArgumentError(
        "confint: no interval is reported for a CMPLE (default se) of a " *
        "dyad-dependent formula — its naive pseudo-likelihood standard errors " *
        "are not calibrated (95% Wald intervals covered 0.70-0.82 for GWESP in " *
        "simulation). Refit with se=:bootstrap (parametric bootstrap) or " *
        "method=:cmle, or pass se=:hessian explicitly to accept the naive " *
        "Wald intervals."))
    # an unconverged CMPLE (separation, Newton cap): withheld, as its z and p
    _cmple_unconverged(r) && return fill(NaN, length(coef(r)), 2)
    q = quantile(Normal(), 1 - (1 - level) / 2)
    θ, se = coef(r), stderror(r)
    return hcat(θ .- q .* se, θ .+ q .* se)
end

# =============================================================================
# Estimation
# =============================================================================

"""
    stergm(networks, formation, dissolution; method=:auto, kwargs...) -> STERGMResult
    fit_stergm(networks, formation, dissolution; method=:auto, kwargs...) -> STERGMResult

Fit a STERGM to a panel of networks (`fit_stergm === stergm`).

`method`:
- `:auto` (default) — resolved by ERGM.jl's rule (`ERGM.resolve_method`):
  the CMPLE when neither formula has a dyad-dependent term (there the
  conditional pseudo-likelihood IS the conditional likelihood, so the CMPLE
  is the exact CMLE and no MCMC is run), the CMLE otherwise. R's `tergm()`
  has no default estimator — it requires `estimate=` — so this default is
  TERGM.jl's choice: the conditional MLE whenever it differs from the
  CMPLE.
- `:cmple` — conditional maximum pseudo-likelihood on the
  formation/dissolution networks ([`cmple`](@ref)). Exact CMLE for
  dyad-independent terms; an approximation for dyad-dependent terms
  (`Mutual`, `GWESP`, ...), for which the default fit reports estimates and
  naive standard errors but no z, p or interval — pass `se = :bootstrap`
  for parametric-bootstrap inference. tergm's `estimate = "CMPLE"`.
- `:cmle` — Monte-Carlo conditional maximum likelihood ([`cmle`](@ref)):
  tergm's `estimate = "CMLE"`. For a dyad-independent formula it returns
  the CMPLE (the two coincide) without running MCMC.
- `:egmme` — **throws an `ArgumentError`** (EGMME is not implemented; see
  `TERGM.egmme`).

Remaining keyword arguments are forwarded to the estimator. A keyword that
estimator does not take is an `ArgumentError` naming the estimator that
does — e.g. `se = :bootstrap` (a CMPLE keyword) on a dyad-dependent formula,
which `method = :auto` fits by CMLE, asks for `method = :cmple`.

`networks` is the panel: a vector of at least two `Network`s observed at
successive times, first. `formation` and `dissolution` are each a term or
any collection of terms `fit_ergm` accepts (a `push!`-built `Vector{Any}`,
a bare `Edges()`, nested vectors), normalised by [`STERGM`](@ref). The
common slips are said in words rather than as a `MethodError`: the panel
in the wrong position (`stergm(terms, terms, networks)`) is
`ArgumentError: arguments are swapped …`, a single network instead of a
panel is `ArgumentError: stergm needs a panel …`.

**Coming from R tergm:** `tergm(nets ~ Form(~edges + mutual) +
Persist(~edges), estimate = "CMLE")` (R requires `estimate=`) is
`stergm(nets, [Edges(), Mutual()], [Edges()])` — `Form()` is the second
argument, `Persist()` the third; `estimate = "CMPLE"` is `method = :cmple`.

# Example
```julia
using TERGM, ERGM
t0 = network(5; directed=true); add_edge!(t0, 1, 2); add_edge!(t0, 3, 4)
t1 = network(5; directed=true); add_edge!(t1, 1, 2); add_edge!(t1, 2, 3)
fit = stergm([t0, t1], [Edges()], [Edges()])
fit.converged                     # true
formation_coef(fit)[1] ≈ log(1 / 17)   # 1 of 18 prior non-ties formed
persistence_coef(fit)[1] ≈ 0.0         # 1 of 2 prior ties persisted
coef(stergm([t0, t1], Edges(), Edges())) == coef(fit)   # true: a bare term needs no brackets
try; stergm([Edges()], [Edges()], [t0, t1]); catch e; occursin("swapped", e.msg); end   # true
```
"""
function stergm(networks::AbstractVector{<:Network}, formation, dissolution;
                method::Symbol=:auto, kwargs...)
    # an unknown method is refused before the model is built
    method in (:auto, :egmme) || _resolve_stergm_method(method, false)
    model = STERGMModel(STERGM(formation, dissolution), networks)
    method === :egmme && return egmme(model; kwargs...)
    est = _resolve_stergm_method(method, has_dyad_dependent(model))
    _check_estimator_keywords(est, method, keys(kwargs))
    return est === :cmple ? cmple(model; kwargs...) : cmle(model; kwargs...)
end

# ERGM.jl's default-estimator rule with TERGM's estimator names: `:auto` is
# the CMPLE (exact) for a dyad-independent formula and the CMLE otherwise
_resolve_stergm_method(method::Symbol, dependent::Bool) =
    ERGM.resolve_method(method, dependent; exact=:cmple, mcmc=:cmle,
                        methods=(:cmple, :cmle), context="stergm")

# The keywords each estimator accepts, read from its own signature so the
# list cannot drift from the code
_estimator_keywords(::Val{:cmple}) = Base.kwarg_decl(which(cmple, Tuple{STERGMModel}))
_estimator_keywords(::Val{:cmle}) = Base.kwarg_decl(which(cmle, Tuple{STERGMModel}))

# A keyword the chosen estimator does not take is refused in words — most
# often a CMPLE keyword (`se=:bootstrap`) on a dyad-dependent formula, which
# `method=:auto` sends to the CMLE
function _check_estimator_keywords(est::Symbol, method::Symbol, keys)
    accepted = _estimator_keywords(Val(est))
    bad = [k for k in keys if !(k in accepted)]
    isempty(bad) && return nothing
    other = est === :cmple ? :cmle : :cmple
    listed = join(("`$k`" for k in bad), ", ")
    why = method === :auto ?
        " method=:auto chose $(repr(est)) because the formula is " *
        (est === :cmle ? "dyad-dependent (the conditional MLE, R tergm's estimate = \"CMLE\")." :
                         "dyad-independent (the CMPLE is the exact CMLE there).") : ""
    hint = all(in(_estimator_keywords(Val(other))), bad) ?
        " $(length(bad) == 1 ? "It is a keyword" : "They are keywords") of " *
        "method=$(repr(other)); pass method=$(repr(other)) explicitly to use " *
        "$(length(bad) == 1 ? "it" : "them")." :
        " See `?cmple` and `?cmle` for the keywords of each estimator."
    throw(ArgumentError("stergm: keyword $listed is not accepted by " *
                        "method=$(repr(est)).$why$hint"))
end

# The common mis-orderings, said in words instead of a MethodError
# (mirroring `fit_ergm`'s swapped-argument methods): the panel in the third
# or the second position, or a term in the first; and a single network where
# a panel is needed. A panel in two slots is never a valid call either. The
# methods whose first argument is `AbstractVector{<:Network}` exist only to
# make these unambiguous against the real entry point
# (`Test.detect_ambiguities(TERGM)` is asserted empty).
const _SWAPPED_ARGS_MSG =
    "arguments are swapped: call stergm(networks, formation, dissolution) — " *
    "the panel (a vector of ≥ 2 networks observed at successive times) comes " *
    "first, then the formation terms, then the dissolution (persistence) " *
    "terms (R tergm: nets ~ Form(~…) + Persist(~…))"
stergm(::AbstractVector, ::Any, ::AbstractVector{<:Network}; kwargs...) =
    throw(ArgumentError(_SWAPPED_ARGS_MSG))
stergm(::AbstractVector{<:Network}, ::Any, ::AbstractVector{<:Network}; kwargs...) =
    throw(ArgumentError(_SWAPPED_ARGS_MSG))
stergm(::AbstractVector, ::AbstractVector{<:Network}, ::Any; kwargs...) =
    throw(ArgumentError(_SWAPPED_ARGS_MSG))
stergm(::AbstractVector{<:Network}, ::AbstractVector{<:Network}, ::Any; kwargs...) =
    throw(ArgumentError(_SWAPPED_ARGS_MSG))
stergm(::AbstractVector, ::AbstractVector{<:Network}, ::AbstractVector{<:Network}; kwargs...) =
    throw(ArgumentError(_SWAPPED_ARGS_MSG))
stergm(::AbstractVector{<:Network}, ::AbstractVector{<:Network}, ::AbstractVector{<:Network}; kwargs...) =
    throw(ArgumentError(_SWAPPED_ARGS_MSG))
stergm(::AbstractERGMTerm, ::Any, ::Any; kwargs...) =
    throw(ArgumentError(_SWAPPED_ARGS_MSG))
stergm(net::Network, ::Any, ::Any; kwargs...) =
    throw(ArgumentError(
        "stergm needs a panel — a vector of at least two networks observed at " *
        "successive times, `stergm([t0, t1, …], formation, dissolution)` — but " *
        "got a single $(typeof(net)). A STERGM is a model of the transitions " *
        "between panels; for one cross-section fit an ERGM (`ERGM.fit_ergm`)."))

const fit_stergm = stergm

# One side's logistic pseudo-likelihood fit over its free dyads. Returns a
# NamedTuple `(coef, se, vcov, loglik, converged, separated, n_kept)`.
#
# The fit is ERGM.jl's `mple_fit_design` — the very function behind
# `ERGM.mple` — on the Bernoulli-row design (`n_tot = 1`, `n_one = y`): the
# CMPLE over the free dyads of Y⁺/Y⁻ IS a logistic pseudo-likelihood, so it
# has to fail in exactly the ways ERGM's does, and R's:
#
# - a statistic at the boundary of its attainable range (a `NodeMatch` no
#   prior same-group tie of which persists, a `Triangle` on a triangle-free
#   auxiliary network, …) has no finite maximizer. As R ergm does under its
#   default `drop=TRUE`, the coefficient is fixed at ∓Inf (standard error 0)
#   and the remaining coefficients are the CMPLE on the free dyads the
#   dropped statistics do not touch — the exact limit of the
#   pseudo-likelihood — with `n_kept` counting those dyads (R's `logLik`
#   nobs attribute, the BIC sample size);
# - a design separated by a COMBINATION of statistics is returned with
#   `converged = false` (R: "The MPLE does not exist!") instead of the point
#   where Newton met its tolerance on the flat asymptote;
# - a Newton iteration that exhausts `maxiter` is `converged = false`.
#
# The derivatives are the shared `NetworkCore.logistic_derivatives` driven by
# `NetworkCore.newton_fit` inside `mple_fit_design`:
# workspaces allocated once, η = Xβ by gemv and −H = X'WX by gemm, never a
# per-row `x * x'` outer product. Never paste a logistic loop back in here.
#
# `warn=true` emits R's sentences through `mple_fit_design`'s own warnings
# (the drop, aliasing, NetworkCore's separation warning) under the `cmple`
# context; the bootstraps' refits pass
# `warn=false` (a boundary in a resampled or simulated panel is not a fact
# about the observed one, and the refits run on every thread) and are
# reported once, in aggregate.
function _logistic_fit(X::Matrix{Float64}, y::Vector{Bool},
                       names::Vector{String}=fill("", size(X, 2));
                       maxiter::Int=100, tol::Float64=1e-8, warn::Bool=false)
    n, p = size(X)
    n > 0 || return (coef=fill(NaN, p), se=fill(NaN, p), vcov=fill(NaN, p, p),
                     loglik=NaN, converged=false, separated=false, n_kept=0)

    n_tot, n_one = ones(n), Float64.(y)
    fit = mple_fit_design(X, n_tot, n_one, names;
                           maxiter=maxiter, tol=tol, warn=warn, context="cmple",
                           estimate="CMPLE")
    return (coef=fit.coefficients, se=fit.std_errors, vcov=fit.var_cov,
            loglik=fit.loglik, converged=fit.converged,
            separated=fit.separated, n_kept=Int(round(fit.n_kept)))
end

# Per-transition CMPLE design blocks: for transition t (2:T), the
# formation rows (prior non-edges, stats on Y⁺) and dissolution rows
# (prior edges, stats on Y⁻). Kept per transition so the bootstraps can
# resample or re-simulate whole transitions.
#
# Every off-diagonal dyad of a transition is free in exactly one of the two
# models, so the block sizes are known before a change statistic is
# evaluated: `ne(prev)` dissolution rows, the rest formation rows. The four
# arrays of a transition are allocated once at those sizes and filled in
# place through `_fill_blocks!`, a function barrier over the two term
# *tuples* — the row fill is the statically unrolled `_fill_row!`, so the
# design build allocates the matrices and vectors it returns (plus the two
# auxiliary networks) and nothing per row (pinned by `@allocated`).
function _cmple_blocks(model::STERGMModel{T,D}) where {T,D}
    n_trans = length(model.networks) - 1

    Xf_blocks = Vector{Matrix{Float64}}(undef, n_trans)
    yf_blocks = Vector{Vector{Bool}}(undef, n_trans)
    Xd_blocks = Vector{Matrix{Float64}}(undef, n_trans)
    yd_blocks = Vector{Vector{Bool}}(undef, n_trans)

    for t in 2:length(model.networks)
        Xf_blocks[t-1], yf_blocks[t-1], Xd_blocks[t-1], yd_blocks[t-1] =
            _transition_blocks(model.formula, model.networks[t-1], model.networks[t])
    end

    return Xf_blocks, yf_blocks, Xd_blocks, yd_blocks
end

# The design block of ONE transition `prev → curr` under an expanded
# formula: `(Xf, yf, Xd, yd)`. `curr` is the observed next panel in the
# CMPLE, or a panel simulated from `prev` in the parametric bootstrap.
function _transition_blocks(formula::STERGM, prev::Network{T,D},
                            curr::Network{T,D}) where {T,D}
    pf, pd = length(formula.formation), length(formula.dissolution)
    n = Int(nv(prev))
    n_dyads = D ? n * (n - 1) : n * (n - 1) ÷ 2
    yplus = formation_network(prev, curr)
    yminus = dissolution_network(prev, curr)

    # The attribute terms' typed snapshots, taken from the auxiliary
    # network of THIS transition (its attributes are Y_{t−1}'s), so the
    # row fill below reads dense vectors, never the attribute Dicts
    fterms = _materialized_tuple(formula.formation, yplus)
    dterms = _materialized_tuple(formula.dissolution, yminus)

    nd = Int(ne(prev))            # prior ties: the dissolution-free dyads
    nf = n_dyads - nd             # prior non-ties: the formation-free dyads
    Xf = Matrix{Float64}(undef, nf, pf)
    yf = Vector{Bool}(undef, nf)
    Xd = Matrix{Float64}(undef, nd, pd)
    yd = Vector{Bool}(undef, nd)
    _fill_blocks!(Xf, yf, Xd, yd, fterms, dterms, prev, curr, yplus, yminus,
                  n, Val(D))
    return Xf, yf, Xd, yd
end

# The row fill of one transition, typed on the term tuples: formation rows
# (prior non-edges, statistics on Y⁺) and dissolution rows (prior edges,
# statistics on Y⁻) in dyad order.
function _fill_blocks!(Xf::Matrix{Float64}, yf::Vector{Bool},
                       Xd::Matrix{Float64}, yd::Vector{Bool},
                       fterms::Tuple, dterms::Tuple, prev, curr, yplus, yminus,
                       n::Int, ::Val{D}) where {D}
    kf = 0
    kd = 0
    for i in 1:n
        for j in (D ? (1:n) : (i+1:n))
            i == j && continue
            if !has_edge(prev, i, j)
                kf += 1
                _fill_row!(Xf, kf, fterms, 1, yplus, i, j, prev)
                @inbounds yf[kf] = has_edge(curr, i, j)
            else
                kd += 1
                _fill_row!(Xd, kd, dterms, 1, yminus, i, j, prev)
                @inbounds yd[kd] = has_edge(curr, i, j)
            end
        end
    end
    (kf == size(Xf, 1) && kd == size(Xd, 1)) ||
        error("TERGM._fill_blocks!: free-dyad count mismatch ($kf formation rows " *
              "for $(size(Xf, 1)) allocated, $kd dissolution rows for " *
              "$(size(Xd, 1))); `ne(prev)` disagrees with the dyad enumeration")
    return nothing
end

# One design row, term by term, statically unrolled over the tuple
@inline function _fill_row!(X::Matrix{Float64}, r::Int, terms::Tuple, k::Int,
                            net, i::Int, j::Int, prev)
    isempty(terms) && return nothing
    @inbounds X[r, k] = _tchange(first(terms), net, i, j, prev)
    return _fill_row!(X, r, Base.tail(terms), k + 1, net, i, j, prev)
end

_stack(blocks::AbstractVector{Matrix{Float64}}, p::Int) =
    isempty(blocks) ? Matrix{Float64}(undef, 0, p) : reduce(vcat, blocks)
_stack(blocks::AbstractVector{Vector{Bool}}) =
    isempty(blocks) ? Bool[] : reduce(vcat, blocks)

# Refit both sides on a set of transition blocks (silently: a boundary in a
# resampled or simulated panel is reported in aggregate). A NaN row when
# either side has no finite CMPLE — a statistic at the boundary of its
# attainable range, a separated design, an empty side — so that a ∓Inf never
# poisons a bootstrap covariance.
function _refit_blocks(Xf_blocks, yf_blocks, Xd_blocks, yd_blocks, pf::Int, pd::Int;
                       maxiter::Int, tol::Float64)
    f = _logistic_fit(_stack(Xf_blocks, pf), _stack(yf_blocks);
                      maxiter=maxiter, tol=tol, warn=false)
    d = _logistic_fit(_stack(Xd_blocks, pd), _stack(yd_blocks);
                      maxiter=maxiter, tol=tol, warn=false)
    finite = f.converged && d.converged &&
             all(isfinite, f.coef) && all(isfinite, d.coef)
    return finite ? vcat(f.coef, d.coef) : fill(NaN, pf + pd)
end

# What both bootstraps do with their replicates: the covariance of the
# finite refits; a replicate without a finite refit is excluded, counted
# once in a warning and kept as a NaN row of the returned replicates.
function _bootstrap_result(boot, n_boot::Int, kind::String, se::Symbol)
    replicates = boot.replicates
    ok = [all(isfinite, view(replicates, b, :)) for b in 1:n_boot]
    n_ok = count(ok)
    n_ok == n_boot && return boot.vcov, replicates

    n_ok >= 2 || throw(ArgumentError(
        "cmple: se=:$se — only $n_ok of the $n_boot $kind refits produced " *
        "finite coefficients (a replicate separated a statistic or left one " *
        "side without free dyads); a covariance needs at least 2. Increase " *
        "n_boot, simplify the formula, or use more panels."))
    @warn "cmple: se=:$se — $(n_boot - n_ok) of the $n_boot $kind " *
          "refits had no finite coefficients and were excluded; the standard " *
          "errors are the empirical covariance of the $n_ok finite refits. " *
          "$_BOOT_EXCLUSION_BIAS " *
          "`fit.boot_replicates` holds every refit (NaN rows excluded); " *
          "`approximations(fit)` records the exclusion."
    return Matrix{Float64}(cov(replicates[ok, :])), replicates
end

# ----------------------------------------------------------------------------
# The per-transition block bootstrap (the btergm approach): resample whole
# time-transitions with replacement, refit both logistic pseudo-likelihoods
# on the resampled rows, and return the empirical covariance of the stacked
# (formation, persistence) coefficients.
#
# Its resampling units are the T transitions, so it needs many of them. In
# simulation (30 actors, `edges + nodematch` on both sides, 300 panels per
# length; 95 % Wald intervals) it covered 0.57–0.65 at 2 transitions,
# 0.74–0.76 at 3, 0.77–0.88 at 4–5, 0.87–0.90 at 8, 0.90–0.92 at 10,
# 0.91–0.95 at 15–20 and 0.93–0.94 at 30, while the inverse-Hessian errors
# covered 0.92–0.97 at every length. (Multiplying the standard error by
# √(T/(T−1)) and using t(T−1) quantiles restored 0.93–0.97 at every length
# in the same simulation; the estimator is kept as btergm defines it.) So
# `se=:block_bootstrap` is refused below `_BLOCK_BOOTSTRAP_MIN_TRANSITIONS` and
# warned about below `_BLOCK_BOOTSTRAP_WARN_TRANSITIONS`.
#
# The loop itself is `NetworkCore.bootstrap_cov` — the ONE shared bootstrap of
# the ecosystem. This function supplies only the two callbacks that are
# TERGM's: `simulate` draws the `n_boot` transition resamples (all of them
# in one call, from `rng`), and `refit` refits the CMPLE on one resample.
# ----------------------------------------------------------------------------
const _BLOCK_BOOTSTRAP_MIN_TRANSITIONS = 10
const _BLOCK_BOOTSTRAP_WARN_TRANSITIONS = 20

# The simulated 95 % coverage of the block bootstrap near `n_trans`
# transitions, for the messages
_block_coverage_note(n_trans::Int) =
    n_trans < 3 ? "0.57-0.65 at 2 transitions" :
    n_trans < 4 ? "0.74-0.76 at 3 transitions" :
    n_trans < 8 ? "0.77-0.88 at 4-5 transitions" :
    n_trans < 15 ? "0.87-0.92 at 8-10 transitions" :
    n_trans < 30 ? "0.91-0.95 at 15-20 transitions" : "0.93-0.94 at 30 transitions"

# The panel-length gate of `se=:block_bootstrap` (in `cmple`; the resampling loop
# below is the scheme itself and runs at any length ≥ 2, which is what the
# golden fixture pins against R on 7 transitions)
function _check_block_bootstrap(n_trans::Int)
    n_trans >= _BLOCK_BOOTSTRAP_MIN_TRANSITIONS || throw(ArgumentError(
        "cmple: se=:block_bootstrap resamples whole time-transitions, and the panel " *
        "has only $n_trans transition$(n_trans == 1 ? "" : "s"); with fewer than " *
        "$(_BLOCK_BOOTSTRAP_MIN_TRANSITIONS) the block-bootstrap standard errors " *
        "are too small (95% intervals covered " *
        "$(_block_coverage_note(n_trans)) in simulation, against 0.94-0.97 for " *
        "the alternatives). Use se=:bootstrap (parametric bootstrap from the " *
        "fitted model, calibrated at any panel length), the default " *
        "inverse-Hessian standard errors (exact for a dyad-independent " *
        "formula), or method=:cmle."))
    n_trans >= _BLOCK_BOOTSTRAP_WARN_TRANSITIONS ||
        @warn "cmple: se=:block_bootstrap with $n_trans transitions — the block " *
              "bootstrap resamples whole transitions and still under-covers " *
              "at this panel length (95% intervals covered " *
              "$(_block_coverage_note(n_trans)) in simulation). Prefer " *
              "se=:bootstrap or method=:cmle; `approximations(fit)` records this."
    return nothing
end

function _cmple_block_bootstrap(Xf_blocks, yf_blocks, Xd_blocks, yd_blocks,
                                θ̂::Vector{Float64}, pf::Int, pd::Int;
                                n_boot::Int, maxiter::Int, tol::Float64,
                                rng::Random.AbstractRNG)
    n_trans = length(Xf_blocks)
    n_trans >= 2 || throw(ArgumentError(
        "the block bootstrap resamples time-transitions and needs at least 2 " *
        "(got $n_trans)"))
    n_boot >= 2 || throw(ArgumentError("n_boot must be at least 2"))

    simulate(rng, B) = [rand(rng, 1:n_trans, n_trans) for _ in 1:B]
    refit(idx) = _refit_blocks(Xf_blocks[idx], yf_blocks[idx], Xd_blocks[idx],
                               yd_blocks[idx], pf, pd; maxiter=maxiter, tol=tol)
    boot = bootstrap_cov(refit, simulate, θ̂; n_boot=n_boot, rng=rng)
    return _bootstrap_result(boot, n_boot, "block-bootstrap", :block_bootstrap)
end

# ----------------------------------------------------------------------------
# The parametric bootstrap of the CMPLE (the analogue of ERGM.jl's
# `mple(se=:bootstrap)`): simulate each transition from the fitted model,
# conditional on the OBSERVED previous panel — Y*_t ~ STERGM(θ̂ | Y_{t−1}),
# exactly the conditional distribution the CMPLE's likelihood is about —
# refit the CMPLE on the simulated transitions (Y_{t−1} → Y*_t), and take
# the empirical covariance of the refits. Unlike the block bootstrap its
# replicates are not limited by the number of transitions, so it is
# calibrated on a two-wave panel too.
#
# The loop is `NetworkCore.bootstrap_cov`. `simulate(rng, B)` draws one seed per
# (replicate, transition) from `rng` up front; `refit(seeds)` simulates and
# refits one replicate on its own `Xoshiro(seed)` streams — deterministic
# given its seeds, so the threaded loop is reproducible from `rng` alone and
# thread-count independent.
# ----------------------------------------------------------------------------
function _cmple_parametric_bootstrap(model::STERGMModel, θf::Vector{Float64},
                                     θd::Vector{Float64};
                                     n_boot::Int, burnin::Union{Nothing,Int},
                                     maxiter::Int, tol::Float64,
                                     rng::Random.AbstractRNG)
    n_boot >= 2 || throw(ArgumentError("n_boot must be at least 2"))
    nets, formula = model.networks, model.formula
    n_trans = length(nets) - 1
    pf, pd = length(θf), length(θd)

    simulate(rng, B) = [rand(rng, UInt64, n_trans) for _ in 1:B]
    function refit(seeds)
        Xf = Vector{Matrix{Float64}}(undef, n_trans)
        yf = Vector{Vector{Bool}}(undef, n_trans)
        Xd = Vector{Matrix{Float64}}(undef, n_trans)
        yd = Vector{Vector{Bool}}(undef, n_trans)
        for t in 1:n_trans
            sim = _simulate_transition(nets[t], formula, θf, θd, burnin,
                                       Random.Xoshiro(seeds[t]))
            Xf[t], yf[t], Xd[t], yd[t] = _transition_blocks(formula, nets[t], sim)
        end
        return _refit_blocks(Xf, yf, Xd, yd, pf, pd; maxiter=maxiter, tol=tol)
    end
    boot = bootstrap_cov(refit, simulate, vcat(θf, θd); n_boot=n_boot, rng=rng)
    return _bootstrap_result(boot, n_boot, "parametric-bootstrap", :bootstrap)
end

"""
    cmple(model::STERGMModel; maxiter=100, tol=1e-8, se=nothing,
          n_boot=200, burnin=nothing, rng=Random.default_rng()) -> STERGMResult

Conditional maximum pseudo-likelihood: for every transition t,

- **formation** rows are the non-edges of Y_{t-1}; the response is edge
  presence in Y_t and change statistics are evaluated on the formation
  network Y⁺ = Y_{t-1} ∪ Y_t;
- **dissolution** rows are the edges of Y_{t-1}; the response is
  persistence into Y_t and change statistics are evaluated on the
  dissolution network Y⁻ = Y_{t-1} ∩ Y_t.

Rows pool across transitions; the two logistic likelihoods are maximized
independently (separability) by ERGM.jl's MPLE design fitter on the shared
`NetworkCore.newton_fit`/`logistic_derivatives` kernel (`maxiter` Newton
iterations, tolerance `tol`).

**Nothing about a bad fit is silent.** A fit that exhausts `maxiter` is
returned with `converged == false`, a warning
(`cmple: the Newton iteration did not converge in maxiter=…`), the caveat
under `show`'s tables and an entry in `approximations(fit)`.

**A statistic at the boundary of its attainable range** — a `NodeMatch`
no prior same-group tie of which persists, a `Triangle` on a triangle-free
auxiliary network, … — has no finite CMPLE. As R ergm does (its default
`drop=TRUE`), `cmple` warns "observed statistic(s) … are at their smallest
attainable values. Their coefficients will be fixed at -Inf", returns that
coefficient as `-Inf` (`+Inf` at the largest value) with standard error 0
and p-value 0, and fits the remaining coefficients on the free dyads the
dropped statistic does not touch — the exact limit of the
pseudo-likelihood. `dof(fit)` counts the finite coefficients and `bic`
uses the dyads they were estimated on (`fit.n_kept`), as R's `logLik`
does; `approximations(fit)` and `show` name the fixed coefficients.
(R `tergm` 4.2 itself does not drop: its operator terms bypass ergm's
attainable-range check, so it warns "The MPLE does not exist!" and returns
a large finite value on the asymptote, e.g. −19.6 with a standard error of
1600 — the same *finite* coefficients to 1e-6, but an arbitrary number
where TERGM.jl reports `-Inf`.)

**Perfect separation by a combination of statistics** — which the boundary
test cannot see — leaves the pseudo-likelihood without a maximum. R warns
"The MPLE does not exist!"; `cmple` decides it with the ecosystem's shared
verdict (`NetworkCore.logistic_separation`, an exact linear programme, via
ERGM.jl's design fitter), warns "cmple: the CMPLE does not exist
(separation)" naming the separating terms, and returns the fit with
`converged == false`; z values, p-values and `confint` are then `NaN`, as
for every unconverged CMPLE.

# Standard errors

`se` selects the standard errors (validated by `NetworkCore.check_se`):

- `:hessian` — inverse observed pseudo-likelihood information. Exact for a
  dyad-independent formula (there the CMPLE is the CMLE). The default
  (`se = nothing`) is `:hessian` with the inference it supports decided by
  the formula — see "Inference under dyadic dependence".
- `:bootstrap` — **parametric bootstrap**: each of `n_boot` replicates
  simulates every transition from the fitted model, conditional on the
  observed previous panel ([`simulate_stergm`](@ref), `burnin` toggles per
  side), refits the CMPLE on the simulated transitions, and the standard
  errors are the empirical covariance of the refits (through the shared
  `NetworkCore.bootstrap_cov`; reproducible from `rng`, identical at any
  thread count). It is the same keyword, and the same procedure, as ERGM.jl's `mple(se = :bootstrap)`
  and is calibrated under dyadic dependence at any panel length (coverage
  0.93–0.98 in the simulations below). The point estimates are unchanged;
  `vcov` becomes the joint bootstrap covariance and `boot_replicates` holds
  the refits. A replicate without a finite refit is excluded with one
  warning and recorded in `approximations(fit)`. Refused (`ArgumentError`)
  for a fit that did not converge or carries a `±Inf` coefficient — a
  panel cannot be simulated at it.
- `:block_bootstrap` — per-transition **block bootstrap** (the `btergm` scheme):
  resample the time-transitions with replacement `n_boot` times and refit.
  Its resampling units are the transitions, so it **under-covers on short
  panels** and is not a remedy for anything there: it is **refused below
  10 transitions** (`ArgumentError`) and warned about below 20. Refused
  also when the point estimate carries a `±Inf` coefficient.

# Inference under dyadic dependence

For a formula with a dyad-dependent term (`has_dyad_dependent(model)`) the
naive inverse-Hessian standard errors are not calibrated in general. In
simulation (30 actors, 200 panels per configuration, 1–10 transitions, 95 %
Wald intervals) they covered 0.89–0.96 for `Form(edges + mutual) +
Persist(edges + mutual)` but only 0.70–0.82 on the formation side of
`edges + GWESP(0.5)`; the parametric bootstrap covered 0.91–0.98 in both.
So, as ERGM.jl's `mple` does, **a CMPLE fit (default `se`) of a dyad-dependent
formula reports its point estimates and naive standard errors but no
inference built on them**: the z and p columns are `NaN` (with a note
saying why), `confint` throws an `ArgumentError`, `fit.inference_withheld`
is `true` and `approximations(fit)` records it. For inference use
`se = :bootstrap` or `method = :cmle`. Passing `se = :hessian`
**explicitly** is the written opt-in to the naive Wald table (what R
`tergm`'s CMPLE summary prints), shown with the pseudo-likelihood caveat.
Dyad-independent formulas are unaffected.

# Example
```julia
using TERGM, ERGM, Random
rng = Xoshiro(1)
panels = [network(12; directed=true) for _ in 1:4]
for i in 1:12, j in 1:12
    i != j && rand(rng) < 0.2 && add_edge!(panels[1], i, j)
end
for t in 2:4, i in 1:12, j in 1:12
    i == j && continue
    keep = has_edge(panels[t-1], i, j) ? rand(rng) < 0.7 : rand(rng) < 0.1
    keep && add_edge!(panels[t], i, j)
end
model = STERGMModel(STERGM([Edges(), Mutual()], [Edges()]), panels)
fit = cmple(model)
fit.inference_withheld          # true: a dyad-dependent CMPLE, default se
boot = cmple(model; se=:bootstrap, n_boot=50, rng=Xoshiro(2))
se_method(boot)                 # :bootstrap
size(boot.boot_replicates)      # (50, 3)
coef(boot) == coef(fit)         # true: only the standard errors change
size(confint(boot))             # (3, 2)
```
"""
function cmple(model::STERGMModel; maxiter::Int=100, tol::Float64=1e-8,
               se::Union{Nothing,Symbol}=nothing, n_boot::Int=200,
               burnin::Union{Nothing,Int}=nothing,
               rng::Random.AbstractRNG=Random.default_rng())
    # `se=nothing` (the default) is `:hessian` whose inference is withheld
    # under dyadic dependence; an explicit `se=:hessian` opts in to the
    # naive Wald table (see "Inference under dyadic dependence")
    naive_opt_in = se === :hessian
    se = something(se, :hessian)
    check_se(se, (:hessian, :bootstrap, :block_bootstrap); context="cmple")
    burnin === nothing || burnin >= 0 ||
        throw(ArgumentError("cmple: burnin must be ≥ 0 (got $burnin)"))
    # Refused before anything is fitted: the answer depends only on the panel
    se === :block_bootstrap && _check_block_bootstrap(length(model.networks) - 1)

    fterms = model.formula.formation
    dterms = model.formula.dissolution
    pf, pd = length(fterms), length(dterms)

    Xf_blocks, yf_blocks, Xd_blocks, yd_blocks = _cmple_blocks(model)
    labels = _stacked_labels(model)

    f = _logistic_fit(_stack(Xf_blocks, pf), _stack(yf_blocks), labels[1:pf];
                      maxiter=maxiter, tol=tol, warn=true)
    d = _logistic_fit(_stack(Xd_blocks, pd), _stack(yd_blocks), labels[pf+1:end];
                      maxiter=maxiter, tol=tol, warn=true)
    fcoef, fse, dcoef, dse = f.coef, f.se, d.coef, d.se
    converged = f.converged && d.converged

    # An unconverged fit is a loud result, never a silent one: it is still
    # returned (`converged == false` is recorded and `approximations` lists
    # it), but said. A separated design has already been warned about by
    # `_logistic_fit` (R's "The MPLE does not exist!"); the Newton cap is the
    # other way to get here. An empty side (no free dyad on any transition)
    # is a NaN fit, said in its own words.
    for (side, fit) in (("formation", f), ("dissolution", d))
        fit.converged && continue
        fit.separated && continue
        if all(isnan, fit.coef)
            @warn "cmple: the $side model has no free dyads on any transition " *
                  "(no prior $(side == "formation" ? "non-tie" : "tie") to " *
                  "form or persist), so its coefficients are NaN and " *
                  "`converged == false`."
        else
            @warn "cmple: the Newton iteration did not converge in " *
                  "maxiter=$maxiter iterations on the $side model (the " *
                  "pseudo-likelihood may be unbounded — perfect separation or " *
                  "a degenerate statistic). The returned coefficients are the " *
                  "last iterate and `converged == false`; raise maxiter, or " *
                  "check the model for a statistic that is constant or " *
                  "perfectly predicts the ties."
        end
    end

    θ = vcat(fcoef, dcoef)
    if se in (:block_bootstrap, :bootstrap) && any(isinf, θ)
        # A coefficient fixed at ∓Inf by a boundary statistic cannot be
        # bootstrapped: the statistic is at its boundary on every resample
        # of the transitions, and a panel cannot be simulated at an infinite
        # coefficient. Say so instead of excluding every replicate.
        fixed = [labels[k] for k in eachindex(θ) if isinf(θ[k])]
        throw(ArgumentError(
            "cmple: se=:$se is not available when a coefficient is " *
            "fixed at ±Inf by a statistic at the boundary of its attainable " *
            "range ($(join(fixed, ", "))): the statistic is at that boundary " *
            "on every resample of the transitions, and a panel cannot be " *
            "simulated at an infinite coefficient, so no bootstrap " *
            "refit can have a finite CMPLE. Remove the term (as R's " *
            "drop=TRUE does) or keep the default se=:hessian, which reports " *
            "standard error 0 for the fixed coefficient and the " *
            "inverse-Hessian errors of the rest."))
    end
    if se in (:block_bootstrap, :bootstrap) && converged && any(isnan, θ)
        # A statistic with no identifiable coefficient (NaN): it does not
        # vary, or is a linear combination of others, so no replicate has a
        # coefficient for it either
        bad = [labels[k] for k in eachindex(θ) if isnan(θ[k])]
        throw(ArgumentError(
            "cmple: se=:$se is not available when a coefficient is not " *
            "identifiable ($(join(bad, ", ")): reported as NaN): remove the " *
            "term(s) and refit."))
    end
    if se === :block_bootstrap
        vcov_joint, replicates =
            _cmple_block_bootstrap(Xf_blocks, yf_blocks, Xd_blocks, yd_blocks,
                                   θ, pf, pd; n_boot=n_boot, maxiter=maxiter,
                                   tol=tol, rng=rng)
        ses = sqrt.(max.(diag(vcov_joint), 0.0))
        fse, dse = ses[1:pf], ses[pf+1:end]
    elseif se === :bootstrap
        converged || throw(ArgumentError(
            "cmple: se=:bootstrap simulates panels from the fitted model, but " *
            "the fit did not converge (see the warning above): its coefficients " *
            "are not an estimate to simulate from."))
        vcov_joint, replicates =
            _cmple_parametric_bootstrap(model, fcoef, dcoef; n_boot=n_boot,
                                        burnin=burnin, maxiter=maxiter, tol=tol,
                                        rng=rng)
        ses = sqrt.(max.(diag(vcov_joint), 0.0))
        fse, dse = ses[1:pf], ses[pf+1:end]
    else
        # Joint covariance of the stacked coefficients: block-diagonal,
        # since the formation and dissolution likelihoods are maximized
        # independently
        vcov_joint = [f.vcov zeros(pf, pd); zeros(pd, pf) d.vcov]
        replicates = Matrix{Float64}(undef, 0, pf + pd)
    end

    # Under dyadic dependence the naive inverse-Hessian errors are not
    # calibrated: unless the caller asked for them in writing, no z, p or
    # interval is built on them
    withheld = se === :hessian && !naive_opt_in && has_dyad_dependent(model)

    return STERGMResult(model, fcoef, fse, dcoef, dse, vcov_joint, f.loglik,
                        d.loglik, converged, :cmple, se, replicates,
                        f.n_kept + d.n_kept, withheld, nothing)
end

# =============================================================================
# Conditional maximum likelihood (Monte-Carlo MLE on the constrained spaces)
# =============================================================================
#
# Separability makes the conditional likelihood of a panel a product of two
# exponential families, each a product over transitions:
#
#   formation:    Y⁺_t ∈ {Y ⊇ Y_{t−1}},  P(Y⁺_t) ∝ exp(θ⁺·g(Y⁺_t, Y_{t−1}))
#   persistence:  Y⁻_t ∈ {Y ⊆ Y_{t−1}},  P(Y⁻_t) ∝ exp(θ⁻·g(Y⁻_t, Y_{t−1}))
#
# so the sufficient statistic of a side is the SUM over transitions of the
# model statistics of its auxiliary network, and the MLE solves
# Σ_t E_θ[g(Y_t)] = Σ_t g(y_t^obs). The expectation is over T−1 independent
# constrained networks: each is sampled by the Metropolis kernel of
# `simulate_stergm` (a uniformly drawn free dyad, toggled through
# `ERGM.mh_toggle!`), and a draw of the pooled statistic is the sum of one
# draw per transition. Statistics are tracked RELATIVE to the observed
# auxiliary network, s(Y) = g(Y) − g(y^obs), so the target of the estimating
# equation is the zero vector and the log-likelihood is −log Σ_Y exp(θ·s(Y)).

# One side of the conditional likelihood: its expanded terms, and per
# transition the previous panel, the observed auxiliary network (the chains'
# starting state) and the free dyads.
struct _CMLESide{N<:Network}
    side::Symbol
    terms::Vector{AbstractERGMTerm}
    prevs::Vector{N}
    starts::Vector{N}
    frees::Vector{Vector{NTuple{2,Int}}}
end

function _CMLESide(model::STERGMModel{T,D}, side::Symbol) where {T,D}
    nets = model.networks
    terms = side === :formation ? model.formula.formation : model.formula.dissolution
    prevs = nets[1:end-1]
    starts = [side === :formation ? formation_network(nets[t-1], nets[t]) :
                                    dissolution_network(nets[t-1], nets[t])
              for t in 2:length(nets)]
    frees = [_free_dyads(prev, side) for prev in prevs]
    return _CMLESide{Network{T,D}}(side, terms, prevs, starts, frees)
end

# The free dyads of one side given the previous panel: its non-edges for the
# formation model, its edges for the dissolution model (dyad order).
function _free_dyads(prev::Network{T,D}, constrain::Symbol) where {T,D}
    n = Int(nv(prev))
    free = NTuple{2, Int}[]
    for i in 1:n
        for j in (D ? (1:n) : (i+1:n))
            i == j && continue
            if constrain === :formation
                has_edge(prev, i, j) || push!(free, (i, j))
            else
                has_edge(prev, i, j) && push!(free, (i, j))
            end
        end
    end
    return free
end

# The thinning interval of the CMLE's chains: two proposals per free dyad
# between recorded draws. A uniformly drawn free dyad is re-proposed after
# `n_free` steps on average, and one Metropolis proposal nearly redraws a
# dyad from its conditional, so draws `2·n_free` steps apart are close to
# independent (the autocorrelation left is measured by the Geyer effective
# sample size behind the convergence tests and the Monte-Carlo errors).
_cmle_interval(n_free::Int) = max(2 * n_free, 1)

# `n × p` draws of the pooled relative statistic Σ_t [g(Y_t) − g(y_t^obs)] at
# coefficients θ: one chain per transition, started at the observed auxiliary
# network, each on its own `Xoshiro(seed)` with the seeds drawn from `rng` up
# front — reproducible from `rng` alone and independent of the thread count.
function _cmle_sample(ctx::_CMLESide, θ::Vector{Float64}, n::Int,
                      burnin::Union{Nothing,Int}, interval::Union{Nothing,Int},
                      rng::Random.AbstractRNG; default_interval::F=_cmle_interval) where {F}
    n_trans = length(ctx.prevs)
    p = length(θ)
    seeds = rand(rng, UInt64, n_trans)
    parts = Vector{Matrix{Float64}}(undef, n_trans)
    spawn_all(n_trans) do t
        out = zeros(n, p)
        free = ctx.frees[t]
        if !isempty(free)
            b = burnin === nothing ? mcmc_defaults(length(free)).burnin : burnin
            iv = interval === nothing ? default_interval(length(free)) : interval
            _mh_constrained_stats!(Random.Xoshiro(seeds[t]), _copy_net(ctx.starts[t]),
                                   ctx.prevs[t], free,
                                   _materialized_tuple(ctx.terms, ctx.prevs[t]),
                                   θ, b, iv, out)
        end
        parts[t] = out
        return nothing
    end
    total = parts[1]
    for t in 2:n_trans
        total .+= parts[t]
    end
    return total
end

# The constrained Metropolis chain of `_mh_constrained!`, recording the
# running statistics (relative to the starting state) of every `interval`-th
# state into the rows of `out`. The statistics are kept current by adding or
# subtracting the accepted move's change statistics.
function _mh_constrained_stats!(rng::Random.AbstractRNG, net::Network, prev::Network,
                                free::Vector{NTuple{2,Int}}, terms::Tuple,
                                θ::Vector{Float64}, burnin::Int, interval::Int,
                                out::Matrix{Float64})
    p = length(terms)
    delta = Vector{Float64}(undef, p)
    stats = zeros(p)
    propose = rng -> @inbounds free[rand(rng, 1:length(free))]
    change! = function (delta, move)
        i, j = move
        _fill_delta!(delta, terms, 1, net, i, j, prev)
        return has_edge(net, i, j)
    end
    apply! = function (move, removal)
        i, j = move
        if removal
            rem_edge!(net, i, j)
            @inbounds for k in 1:p
                stats[k] -= delta[k]
            end
        else
            add_edge!(net, i, j)
            @inbounds for k in 1:p
                stats[k] += delta[k]
            end
        end
        return nothing
    end
    on_sample = function (k)
        @inbounds for c in 1:p
            out[k, c] = stats[c]
        end
        return nothing
    end
    mh_toggle!(rng, θ, delta, propose, change!, apply!, on_sample;
               burnin=burnin, interval=interval, n_samples=size(out, 1))
    return out
end

# The thinning interval of the bridge's chains: ERGM.jl's dyad-scaled rule.
# A rung needs only the MEAN of the statistics, for which closely spaced
# draws waste nothing, so the bridge thins far less than the estimation
# chains (whose draws feed a covariance and the convergence tests).
_bridge_interval(n_free::Int) = mcmc_defaults(n_free).interval

# The logistic log-likelihood Σ [y·η − log(1 + e^η)], η = Xθ, of one side's
# CMPLE design — the exact conditional log-likelihood at a θ whose
# dyad-dependent coordinates are zero.
function _logistic_loglik(X::Matrix{Float64}, y::Vector{Bool}, θ::Vector{Float64})
    ll = 0.0
    η = X * θ
    @inbounds for r in eachindex(y)
        ll += (y[r] ? η[r] : 0.0) - log1p(exp(-abs(η[r]))) - max(η[r], 0.0)
    end
    return ll
end

# Path-sampled conditional log-likelihood of one side at θ̂. With the
# statistics measured relative to the observed network, ℓ(θ) = −log Z(θ)
# and d/dθ log Z = E_θ[s], so along θ(u) = θ_ref + u·(θ̂ − θ_ref)
#   ℓ(θ̂) = ℓ(θ_ref) − ∫₀¹ E_{θ(u)}[s]·(θ̂ − θ_ref) du,
# integrated by Simpson's rule over `rungs` intervals, each rung an MCMC
# sample. The reference is dyad-independent — the CMPLE of the side's
# dyad-independent terms alone (the dyad-dependent coefficients at zero) —
# where the likelihood is the logistic one and `ℓ(θ_ref)` is exact.
function _cmle_bridge_loglik(ctx::_CMLESide, X::Matrix{Float64}, y::Vector{Bool},
                             θ̂::Vector{Float64}; rungs::Int, n_samples::Int,
                             burnin, interval, maxiter::Int, tol::Float64,
                             rng::Random.AbstractRNG)
    p = length(θ̂)
    indep = findall(t -> !is_dyad_dependent(t), ctx.terms)
    θ_ref = zeros(p)
    if !isempty(indep)
        sub = _logistic_fit(X[:, indep], y; maxiter=maxiter, tol=tol, warn=false)
        if sub.converged && all(isfinite, sub.coef)
            θ_ref[indep] .= sub.coef
        else
            θ_ref[indep] .= θ̂[indep]
        end
    end
    ℓ_ref = _logistic_loglik(X, y, θ_ref)
    # ERGM.jl's bridge integral log Z(θ̂) − log Z(θ_ref); each rung's mean
    # comes from the CMLE's own pooled sampler (rungs run in order, drawing
    # their seeds from `rng`)
    mean_stats = (θu, k) -> vec(mean(_cmle_sample(ctx, Vector{Float64}(θu), n_samples,
                                                  burnin, interval, rng;
                                                  default_interval=_bridge_interval);
                                     dims=1))
    return ℓ_ref - bridge_integrate(mean_stats, θ_ref, θ̂; rungs=rungs)
end

"""
    cmle(model::STERGMModel; n_samples=1000, burnin=nothing, interval=nothing,
         maxiter=60, termination=:confidence, conv_precision=0.1,
         conv_confidence=0.99, conv_threshold=0.1, hotelling_alpha=0.05,
         gamma0=0.1, max_step_norm=5.0, max_n_samples=16 * n_samples,
         init=nothing, bridge_rungs=16, bridge_samples=n_samples,
         newton_maxiter=100, tol=1e-8, rng=Random.default_rng(),
         verbose=false) -> STERGMResult

Conditional maximum likelihood by MCMC — R `tergm`'s `estimate = "CMLE"`.

By separability the conditional likelihood is a product of a formation and
a persistence likelihood, each an exponential family on a constrained
sample space: Y⁺_t ranges over the networks containing Y_{t−1} with
statistics evaluated on Y⁺_t, Y⁻_t over the sub-networks of Y_{t−1}, and
the sufficient statistic of a side is its model statistics summed over
transitions. Each side is fit separately:

- A side with **no dyad-dependent term** is fit exactly: its CMPLE is its
  CMLE, so no MCMC is run. For a dyad-independent formula `cmle` therefore
  returns the coefficients, standard errors and log-likelihood of
  [`cmple`](@ref) (`fit.mcmc.formation.exact == true`).
- A side with a dyad-dependent term is fit by **Monte-Carlo MLE**, started
  at its CMPLE (or `init`). Every iteration draws `n_samples` pooled
  statistics — one Metropolis chain per transition on the free dyads of
  that side (the kernel of [`simulate_stergm`](@ref), `ERGM.mh_toggle!`
  with a uniformly drawn free dyad), started at the observed auxiliary
  network, `burnin` toggles (default: ERGM.jl's dyad-scaled rule, 20 per
  free dyad), one draw every `interval` toggles (default: 2 per free
  dyad) — and runs ERGM.jl's MCMLE iteration (`ERGM.Extension.mcmle_solve`): the update
  `θ += γ·Σ̂⁻¹(g_obs − ḡ)` with Hummel's step length `γ` (from `gamma0`,
  at most doubling per iteration, norm-capped at `max_step_norm`).
  The stopping rule is ERGM.jl's, hence R ergm 4's: at full step length,
  `termination = :confidence` (default) declares convergence when the
  estimating equation at the *updated* coefficients lies, with confidence
  `conv_confidence`, inside the tolerance region of precision
  `conv_precision` (`ERGM.Extension.confidence_test`); when the test fails near the
  solution the next sample is enlarged as R does, up to `max_n_samples`.
  `termination = :hotelling` is the older rule (every t-ratio
  `|g_obs − ḡ|/sd` below `conv_threshold` and a Hotelling T² test
  non-significant at `hotelling_alpha`).

**Standard errors** are the inverse Fisher information estimated from the
final MCMC sample plus the Monte-Carlo component of the estimate
(`ERGM.Extension.mcmle_covariance`); `se_method(fit) == :fisher`, and
`fit.mcmc.<side>.mcmc_se` is the Monte-Carlo standard error alone.

**The log-likelihood** (hence AIC/BIC) is path-sampled: Simpson's rule over
`bridge_rungs` intervals (`bridge_samples` draws per rung, thinned by
ERGM.jl's dyad-scaled interval unless `interval` is given) from a
dyad-independent reference — the CMPLE of the side's dyad-independent
terms — where the likelihood is exact. `bridge_rungs = 0` skips it
(`loglikelihood`, `aic`, `bic` are `NaN`).

`init` is a stacked start (formation, then persistence — `coef(fit)` of an
earlier fit, to continue an unconverged one). `newton_maxiter` and `tol`
control the CMPLE starts and the exact sides. All randomness comes from
`rng`: one seed per chain is drawn from it up front, so a fit is
reproducible and identical at any thread count.

**Refusals and loud results.** A dyad-dependent side whose CMPLE start
does not exist (a statistic at the boundary of its attainable range — no
finite MLE exists either — a separated design, or an empty side) is
refused with an `ArgumentError`. A fit that exhausts `maxiter` is returned
with `converged == false`, a warning quoting the stopping rule's verdict,
the caveat in `show` and an entry in `approximations(fit)`; `fit.mcmc`
holds the per-side diagnostics (`iterations`, `step_length`,
`termination`, `termination_p`, and — recomputed on the final sample —
`t_ratios`, `hotelling_p`, `n_eff`, `n_samples`, `mcmc_se`).

**Validation.** With `Mutual` as the only dyad-dependent term the
conditional likelihood factorises over pairs of actors and the exact CMLE
is computable by enumeration; the test suite compares `cmle` with it and
with `tergm`'s CMLE on two panels (fixture `cmle_stergm.toml`).

# Example
```julia
using TERGM, ERGM, Random
rng = Xoshiro(1)
panels = [network(12; directed=true) for _ in 1:4]
for i in 1:12, j in 1:12
    i != j && rand(rng) < 0.2 && add_edge!(panels[1], i, j)
end
for t in 2:4, i in 1:12, j in 1:12
    i == j && continue
    keep = has_edge(panels[t-1], i, j) ? rand(rng) < 0.7 : rand(rng) < 0.1
    keep && add_edge!(panels[t], i, j)
end
model = STERGMModel(STERGM([Edges(), Mutual()], [Edges()]), panels)
fit = cmle(model; n_samples=500, rng=Xoshiro(2))
fit.method, se_method(fit)            # (:cmle, :fisher)
fit.mcmc.formation.exact              # false: Mutual needs MCMC
fit.mcmc.persistence.exact            # true: edges alone is fit exactly
persistence_coef(fit) == persistence_coef(cmple(model))   # true
size(confint(fit))                    # (3, 2)

exact = cmle(STERGMModel(STERGM([Edges()], [Edges()]), panels))
coef(exact) == coef(cmple(STERGMModel(STERGM([Edges()], [Edges()]), panels)))   # true
```
"""
function cmle(model::STERGMModel; n_samples::Int=1000,
              burnin::Union{Nothing,Int}=nothing,
              interval::Union{Nothing,Int}=nothing,
              maxiter::Int=60, termination::Symbol=:confidence,
              conv_precision::Float64=0.1, conv_confidence::Float64=0.99,
              conv_threshold::Float64=0.1,
              hotelling_alpha::Float64=0.05, gamma0::Float64=0.1,
              max_step_norm::Float64=5.0,
              max_n_samples::Union{Nothing,Int}=nothing,
              init::Union{Nothing,AbstractVector{<:Real}}=nothing,
              bridge_rungs::Int=16, bridge_samples::Union{Nothing,Int}=nothing,
              newton_maxiter::Int=100, tol::Float64=1e-8,
              rng::Random.AbstractRNG=Random.default_rng(),
              verbose::Bool=false)
    n_samples >= 2 || throw(ArgumentError("cmle: n_samples must be ≥ 2 (got $n_samples)"))
    maxiter >= 1 || throw(ArgumentError("cmle: maxiter must be ≥ 1 (got $maxiter)"))
    burnin === nothing || burnin >= 0 ||
        throw(ArgumentError("cmle: burnin must be ≥ 0 (got $burnin)"))
    interval === nothing || interval >= 1 ||
        throw(ArgumentError("cmle: interval must be ≥ 1 (got $interval)"))
    bridge_rungs >= 0 || throw(ArgumentError(
        "cmle: bridge_rungs must be ≥ 0 (got $bridge_rungs); 0 skips the " *
        "log-likelihood estimate (loglik/AIC/BIC are then NaN)"))
    n_max = something(max_n_samples, 16 * n_samples)
    n_max >= n_samples || throw(ArgumentError(
        "cmle: max_n_samples ($n_max) must be ≥ n_samples ($n_samples)"))
    0 < gamma0 <= 1 || throw(ArgumentError("cmle: gamma0 must be in (0, 1] (got $gamma0)"))
    termination in (:confidence, :hotelling) || throw(ArgumentError(
        "cmle: termination must be :confidence (R ergm's equivalence test) or " *
        ":hotelling (t-ratios + Hotelling test); got $(repr(termination))"))
    conv_precision > 0 || throw(ArgumentError(
        "cmle: conv_precision must be positive (got $conv_precision)"))
    0 < conv_confidence < 1 || throw(ArgumentError(
        "cmle: conv_confidence must be in (0, 1) (got $conv_confidence)"))

    pf, pd = length(model.formula.formation), length(model.formula.dissolution)
    init === nothing || length(init) == pf + pd || throw(ArgumentError(
        "cmle: init has length $(length(init)) but the model has $(pf + pd) " *
        "coefficients ($pf formation, then $pd persistence)"))
    Xf_blocks, yf_blocks, Xd_blocks, yd_blocks = _cmple_blocks(model)
    labels = _stacked_labels(model)
    nb = something(bridge_samples, n_samples)

    function fit_side(side::Symbol, X, y, names, start)
        ctx = _CMLESide(model, side)
        word = side === :formation ? "formation" : "persistence"
        dependent = any(is_dyad_dependent, ctx.terms)
        # The CMPLE of the side: the answer for a dyad-independent side, the
        # starting point otherwise
        pl = _logistic_fit(X, y, names; maxiter=newton_maxiter, tol=tol,
                           warn=!dependent)
        if !dependent
            pl.converged || pl.separated || @warn "cmle: the $word model is " *
                "dyad-independent and was fit exactly by its CMPLE, which did " *
                "not converge (an empty side or an exhausted Newton cap); " *
                "`converged == false`."
            p = length(pl.coef)
            return (coef=pl.coef, se=pl.se, vcov=pl.vcov, loglik=pl.loglik,
                    converged=pl.converged, n_kept=pl.n_kept,
                    report=(exact=true, iterations=0, step_length=1.0,
                            termination=:none, termination_p=NaN,
                            conv_confidence=conv_confidence,
                            t_ratios=zeros(p), hotelling_p=1.0, n_eff=NaN,
                            n_samples=0, mcmc_se=zeros(p)))
        end
        if start === nothing
            (pl.converged && all(isfinite, pl.coef)) || throw(ArgumentError(
                "cmle: the CMPLE of the $word model, used as the starting point, " *
                "does not exist — " *
                (any(isinf, pl.coef) ?
                    "a statistic is at the boundary of its attainable range " *
                    "($(join([names[k] for k in eachindex(names) if isinf(pl.coef[k])], ", "))), " *
                    "so its coefficient has no finite MLE either" :
                 all(isnan, pl.coef) ?
                    "the side has no free dyads on any transition" :
                 any(isnan, pl.coef) ?
                    "statistic(s) " *
                    join([names[k] for k in eachindex(names) if isnan(pl.coef[k])], ", ") *
                    " do not vary on the free dyads, or are linear combinations of " *
                    "the statistics before them, so the CMPLE has no coefficient for " *
                    "them (NaN; R reports NA, and R's CMLE a finite value the data do " *
                    "not identify)" :
                    "the design is separated or Newton did not converge (R: " *
                    "\"The MPLE does not exist!\")") *
                ". Remove or coarsen the term, fit by `cmple` (which reports " *
                "the ±Inf coefficient), or supply a start with `init=`."))
            θ0 = copy(pl.coef)
        else
            all(isfinite, start) || throw(ArgumentError(
                "cmle: init must be finite (got $(start) for the $word model)"))
            θ0 = Vector{Float64}(start)
        end
        sample = (θ, n) -> _cmle_sample(ctx, θ, n, burnin, interval, rng)
        sol = mcmle_solve(sample, θ0; labels=names, n_samples=n_samples,
                           max_n_samples=n_max, maxiter=maxiter,
                           termination=termination,
                           conv_precision=conv_precision,
                           conv_confidence=conv_confidence,
                           conv_threshold=conv_threshold,
                           hotelling_alpha=hotelling_alpha, gamma0=gamma0,
                           max_step_norm=max_step_norm, verbose=verbose,
                           context="cmle ($word)")
        report = (exact=false, iterations=sol.iterations, step_length=sol.step_length,
                  termination=termination, termination_p=sol.termination_p,
                  conv_confidence=conv_confidence, t_ratios=sol.tests.t_ratios,
                  hotelling_p=sol.tests.hotelling_p, n_eff=sol.tests.n_eff,
                  n_samples=size(sol.final.samples, 1), mcmc_se=sol.mcmc_se)
        sol.converged || @warn "cmle: the $word model did not converge in " *
            "maxiter=$maxiter iterations ($(_termination_detail(report)); " *
            "last max t-ratio " *
            "$(round(maximum(report.t_ratios), digits=3)), step length γ " *
            "$(round(report.step_length, digits=3))): the estimates are the " *
            "last iterate and the standard errors are unreliable; raise " *
            "maxiter/n_samples, or continue with `init = coef(fit)`."
        return (coef=sol.coef, se=sol.se, vcov=sol.vcov, loglik=NaN,
                converged=sol.converged, n_kept=length(y), report=report,
                ctx=ctx, X=X, y=y)
    end

    Xf, yf = _stack(Xf_blocks, pf), _stack(yf_blocks)
    Xd, yd = _stack(Xd_blocks, pd), _stack(yd_blocks)
    f = fit_side(:formation, Xf, yf, labels[1:pf],
                 init === nothing ? nothing : init[1:pf])
    d = fit_side(:dissolution, Xd, yd, labels[pf+1:end],
                 init === nothing ? nothing : init[pf+1:end])

    # The path-sampled log-likelihoods, run only now so that they consume
    # randomness after the estimates: coefficients and standard errors are
    # identical with and without the bridge
    function side_loglik(s)
        s.report.exact && return s.loglik
        bridge_rungs == 0 && return NaN
        return _cmle_bridge_loglik(s.ctx, s.X, s.y, s.coef; rungs=bridge_rungs,
                                   n_samples=nb, burnin=burnin, interval=interval,
                                   maxiter=newton_maxiter, tol=tol, rng=rng)
    end
    ll_f = side_loglik(f)
    ll_d = side_loglik(d)

    vcov_joint = [f.vcov zeros(pf, pd); zeros(pd, pf) d.vcov]
    return STERGMResult(model, f.coef, f.se, d.coef, d.se, vcov_joint, ll_f, ll_d,
                        f.converged && d.converged, :cmle, :fisher,
                        Matrix{Float64}(undef, 0, pf + pd), f.n_kept + d.n_kept,
                        false, (formation=f.report, persistence=d.report))
end

"""
    TERGM.egmme(model::STERGMModel; kwargs...)

Equilibrium Generalized Method of Moments estimation is **not implemented**
and this function is **not exported** — exporting it would advertise an
estimator that does not exist. It is retained, and reachable as
`TERGM.egmme`, only so that calling it raises an explanatory `ArgumentError`
rather than an `UndefVarError`; `stergm(...; method = :egmme)` throws the
same error.

Use [`cmple`](@ref) or [`cmle`](@ref) on panel data.
"""
function egmme(model::STERGMModel; kwargs...)
    throw(ArgumentError(
        "EGMME (tergm's estimate = \"EGMME\", for a single cross-section " *
        "observed at equilibrium) is not implemented in TERGM.jl; fit panel " *
        "data with method = :cmle (the default for a dyad-dependent formula) or " *
        "method = :cmple instead"))
end

# =============================================================================
# Simulation
# =============================================================================

# Fill `delta` with the add-direction change statistics of every term of a
# tuple — a statically unrolled recursion, so the sampler's inner step does
# not dispatch dynamically through an abstract term vector.
@inline function _fill_delta!(delta::Vector{Float64}, terms::Tuple, k::Int,
                              net, i::Int, j::Int, prev)
    isempty(terms) && return nothing
    @inbounds delta[k] = _tchange(first(terms), net, i, j, prev)
    return _fill_delta!(delta, Base.tail(terms), k + 1, net, i, j, prev)
end

# Metropolis sampling of Y⁺ (constrain=:formation, free dyads = non-edges
# of prev) or Y⁻ (constrain=:dissolution, free dyads = edges of prev),
# starting from Y_{t-1}. `steps === nothing` resolves to ERGM.jl's
# dyad-scaled default (`mcmc_defaults`: 20 toggles per free dyad).
function _sample_constrained(prev::Network{T,D}, terms, θ::Vector{Float64},
                             constrain::Symbol, steps::Union{Nothing,Int},
                             rng::Random.AbstractRNG) where {T,D}
    net = _copy_net(prev)
    free = _free_dyads(prev, constrain)
    isempty(free) && return net

    burnin = steps === nothing ? mcmc_defaults(length(free)).burnin : steps
    # Typed snapshots of the attribute terms from `prev` (the network the
    # change statistics condition on; the draw starts as its copy)
    _mh_constrained!(rng, net, prev, free, _materialized_tuple(terms, prev), θ, burnin)
    return net
end

# The Metropolis loop is ERGM.jl's exported `mh_toggle!` kernel: one toggle
# loop for the family. A move is a dyad drawn uniformly from `free` — the
# symmetric proposal on the constrained sample space (every free dyad is
# equally likely to be proposed from either of its states), so no Hastings
# correction is needed and the chain's stationary distribution is exactly
# the constrained ERGM (pinned against enumeration in the tests). ERGM.jl's
# tie/no-tie proposal is not used here: on the free dyads of one side a
# single uniform proposal already nearly redraws the dyad from its
# conditional. `change!` fills the temporal change statistics and
# reports whether the dyad is currently a tie (a removal); `apply!` toggles
# it. No sample is recorded (`n_samples=0`): the state after `burnin` toggles
# is the draw. The kernel consumes `rng` in the order the pre-0.2 hand-written
# loop did — the proposal index first, then one uniform for the acceptance
# test — so the sampled sequence is bit-identical. A function barrier: the
# term tuple is concretely typed here, so the closures and the kernel are
# fully typed and allocation-free per step.
function _mh_constrained!(rng::Random.AbstractRNG, net::Network, prev::Network,
                          free::Vector{NTuple{2,Int}}, terms::Tuple,
                          θ::Vector{Float64}, burnin::Int)
    delta = Vector{Float64}(undef, length(terms))
    propose = rng -> @inbounds free[rand(rng, 1:length(free))]
    change! = function (delta, move)
        i, j = move
        _fill_delta!(delta, terms, 1, net, i, j, prev)
        return has_edge(net, i, j)
    end
    apply! = function (move, removal)
        i, j = move
        removal ? rem_edge!(net, i, j) : add_edge!(net, i, j)
        return nothing
    end
    mh_toggle!(rng, θ, delta, propose, change!, apply!, _ -> nothing;
               burnin=burnin, interval=1, n_samples=0)
    return net
end

"""
    simulate_stergm(prev_net, formula::STERGM, θ_form, θ_diss;
                    burnin=nothing, rng=Random.default_rng()) -> Network

Simulate one STERGM transition from `prev_net`: draw the formation
network Y⁺ (Metropolis over the non-edges of `prev_net`, under the
formation coefficients `θ_form`) and the dissolution network Y⁻ (over its
edges, under the **persistence** coefficients `θ_diss`), then combine
`Y_t = (Y⁺ \\ Y_{t−1}) ∪ Y⁻`.

`burnin` is the number of Metropolis toggles per side; `nothing` (the
default) resolves to ERGM.jl's dyad-scaled rule (`ERGM.Extension.mcmc_defaults`:
20 toggles per free dyad of that side). Both samplers run on the shared
`ERGM.mh_toggle!` kernel and draw only from `rng`. A move toggles one dyad
drawn uniformly from the free dyads of the side — the symmetric proposal on
the constrained sample space, whose stationary distribution is exactly the
constrained model (verified against enumeration in the test suite).

`formula` is validated against `prev_net` with ERGM.jl's formula rules (a
nodal term's attribute must exist there; `ArgumentError` naming the side
and the term otherwise) and expanded exactly as [`STERGMModel`](@ref)
expands it, so a raw multi-level `NodeFactor(:grp)` takes one coefficient
per level and a fitted model's expanded formula passes through unchanged.
`θ_form`/`θ_diss` must have one entry per expanded term; the
`ArgumentError` otherwise names the side's terms, both counts and the
per-side accessors (`formation_coef(fit)` / `persistence_coef(fit)`, not
the stacked `coef(fit)`).

# Example
```julia
using TERGM, ERGM, Random
prev = network(10; directed=true)
for i in 1:10, j in 1:10
    i != j && (i + j) % 3 == 0 && add_edge!(prev, i, j)
end
formula = STERGM([Edges()], [Edges()])
yt = simulate_stergm(prev, formula, [-30.0], [30.0]; rng=Xoshiro(1))
ne(yt) == ne(prev)     # true: no formation, full persistence
```
"""
function simulate_stergm(prev_net::Network, formula::STERGM,
                         θ_form::Vector{Float64}, θ_diss::Vector{Float64};
                         burnin::Union{Nothing,Int}=nothing,
                         rng::Random.AbstractRNG=Random.default_rng())
    # The formula is held to ERGM's formula rules on the starting network
    # (a nodal term's attribute must exist there) and expanded to one term
    # per statistic, exactly as `STERGMModel` does — so a raw multi-level
    # `NodeFactor(:grp)` takes one coefficient per level here too, and the
    # expanded formula of a fitted model passes through unchanged.
    for (side, terms) in (("formation", formula.formation),
                          ("dissolution", formula.dissolution))
        _refuse_offsets(terms, "simulate_stergm: $side")
        try
            validate_formula(TermSet(terms), prev_net)
        catch e
            e isa ArgumentError || rethrow()
            throw(ArgumentError("simulate_stergm: $side model: " * e.msg))
        end
    end
    formula = _expand_formula(formula, prev_net)
    net = prev_net
    for (side, terms, θ, accessor) in (("formation", formula.formation, θ_form,
                                        "formation_coef(fit)"),
                                       ("dissolution (persistence)", formula.dissolution,
                                        θ_diss, "persistence_coef(fit)"))
        length(θ) == length(terms) || throw(ArgumentError(
            "simulate_stergm: the $side model has $(length(terms)) term(s) " *
            "($(join(_labels(terms, net), ", "))) but $(length(θ)) $side " *
            "coefficient(s) were given ($(θ)); pass one per term — " *
            "$accessor for a fitted model, not the stacked coef(fit)."))
    end
    burnin === nothing || burnin >= 0 ||
        throw(ArgumentError("burnin must be ≥ 0 (got $burnin)"))
    return _simulate_transition(prev_net, formula, θ_form, θ_diss, burnin, rng)
end

# One simulated transition under an already validated and expanded formula
# (what `simulate_stergm` runs after its checks, and the parametric
# bootstrap runs per replicate).
function _simulate_transition(prev_net::Network, formula::STERGM,
                              θ_form::Vector{Float64}, θ_diss::Vector{Float64},
                              burnin::Union{Nothing,Int}, rng::Random.AbstractRNG)
    yplus = _sample_constrained(prev_net, formula.formation, θ_form,
                                :formation, burnin, rng)
    yminus = _sample_constrained(prev_net, formula.dissolution, θ_diss,
                                 :dissolution, burnin, rng)

    # Y_t = new formations ∪ persisted ties (attribute-preserving: start
    # from a copy of Y_{t-1}, drop dissolved ties, add formed ones)
    yt = _copy_net(prev_net)
    for e in edges(prev_net)
        has_edge(yminus, src(e), dst(e)) || rem_edge!(yt, src(e), dst(e))
    end
    for e in edges(yplus)
        has_edge(prev_net, src(e), dst(e)) || add_edge!(yt, src(e), dst(e))
    end
    return yt
end

"""
    simulate_stergm(result::STERGMResult, n_steps; kwargs...) -> Vector{Network}

Simulate `n_steps` transitions forward from the last observed panel of a
fitted STERGM, at its fitted formation and persistence coefficients
(keywords as for [`simulate_network_sequence`](@ref)).

# Example
```julia
using TERGM, ERGM, Random
t0 = network(6; directed=true); add_edge!(t0, 1, 2); add_edge!(t0, 3, 4)
t1 = network(6; directed=true); add_edge!(t1, 1, 2); add_edge!(t1, 2, 3)
fit = stergm([t0, t1], [Edges()], [Edges()])
future = simulate_stergm(fit, 3; rng=Xoshiro(1))
length(future)          # 4: the last observed panel plus 3 simulated ones
```
"""
function simulate_stergm(result::STERGMResult, n_steps::Int; kwargs...)
    return simulate_network_sequence(result.model.formula,
                                     result.model.networks[end], n_steps,
                                     result.formation_coef,
                                     result.persistence_coef; kwargs...)
end

"""
    simulate_network_sequence(formula, init_net, n_steps, θ_form, θ_diss;
                              burnin=nothing, rng=Random.default_rng()) -> Vector{Network}

Simulate a sequence of `n_steps` STERGM transitions starting from
`init_net` (the returned sequence includes the initial network as its first
element and has the concrete element type of `init_net`). `θ_diss` are
**persistence** coefficients; `burnin`/`rng` as for
[`simulate_stergm`](@ref).

# Example
```julia
using TERGM, ERGM, Random
init = network(8; directed=true)
for i in 1:8, j in 1:8
    i != j && (i + j) % 4 == 0 && add_edge!(init, i, j)
end
seq = simulate_network_sequence(STERGM([Edges()], [Edges()]), init, 5,
                                [-2.0], [1.0]; rng=Xoshiro(3))
length(seq)                      # 6
eltype(seq)                      # Network{Int64, true}
```
"""
function simulate_network_sequence(formula::STERGM, init_net::Network,
                                   n_steps::Int,
                                   θ_form::Vector{Float64},
                                   θ_diss::Vector{Float64};
                                   burnin::Union{Nothing,Int}=nothing,
                                   rng::Random.AbstractRNG=Random.default_rng())
    seq = typeof(init_net)[_copy_net(init_net)]
    for _ in 1:n_steps
        push!(seq, simulate_stergm(seq[end], formula, θ_form, θ_diss;
                                   burnin=burnin, rng=rng))
    end
    return seq
end

# =============================================================================
# Goodness of fit
# =============================================================================

"""
    gof(result::STERGMResult; n_sim=50, burnin=nothing, rng=...) -> GOFResult

Transition-level goodness of fit, after `gof.tergm`: for each observed
transition Y_{t−1} → Y_t, simulate `n_sim` STERGM transitions from the same
starting panel Y_{t−1} at the fitted coefficients, and compare — pooled
over transitions — what was observed to what the model reproduces.

Extends NetworkCore.jl's shared `gof` generic and returns the shared
`NetworkCore.GOFResult` container with these `GOFStatistic` panels, in order:

1. `"tie changes"` (levels `formed`, `persisted`): the pooled counts of
   ties that formed and ties that persisted. **Not a diagnostic for a
   formula containing `Edges()` on that side** — those counts are the
   sufficient statistics of `Form~edges` / `Persist~edges`, which the CMPLE
   reproduces in expectation by construction (p-values near 1 are expected);
   it is a sanity check of the simulator, and the one panel that speaks for
   an edges-free formula.
2. `"formation statistics"` and `"persistence statistics"`: every term of
   the formation model evaluated on the observed formation network Y⁺ =
   Y_{t−1} ∪ Y_t against its value on the simulated Y⁺, and the same for
   the persistence model on Y⁻ = Y_{t−1} ∩ Y_t (labels are the fitted
   coefficients'). The `edges` level is again reproduced by construction;
   the others are the model-statistic check of `gof.tergm` — the first place
   a **dyad-dependent CMPLE** (a `Mutual`, `Triangle`, `GWESP`, … whose
   pseudo-likelihood point estimate may be biased) shows a misfit.
3. The degree distribution of the simulated Y_t against the observed Y_t,
   pooled over transitions — `"idegree"` and `"odegree"` (levels `0` … `n−1`)
   for a directed panel, `"degree"` for an undirected one: the
   out-of-model check, informative for every formula (tergm's default
   `gof` statistics are the degree and duration distributions).

Two-sided Monte-Carlo p-values use the `(1 + k)/(N + 1)` estimator (never
exactly zero). One seed per (transition, simulation) is drawn from `rng` up
front and each simulation runs on its own `Xoshiro(seed)` on every thread,
so the result is reproducible from `rng` alone and thread-count
independent. `burnin` is forwarded to [`simulate_stergm`](@ref). `gof(fit)`
is the one goodness-of-fit verb of the model family.

# Example
```julia
using TERGM, ERGM, Random
t0 = network(6; directed=true); add_edge!(t0, 1, 2); add_edge!(t0, 3, 4)
t1 = network(6; directed=true); add_edge!(t1, 1, 2); add_edge!(t1, 2, 3)
fit = stergm([t0, t1], [Edges()], [Edges()])
g = gof(fit; n_sim=20, rng=Xoshiro(1))
[s.name for s in g.statistics]   # ["tie changes", "formation statistics", "persistence statistics", "idegree", "odegree"]
g.statistics[1].labels           # ["formed", "persisted"]
g.statistics[2].labels           # ["edges"]
```
"""
function gof(result::STERGMResult; n_sim::Int=50,
             burnin::Union{Nothing,Int}=nothing,
             rng::Random.AbstractRNG=Random.default_rng())
    model = result.model
    nets = model.networks
    formula = model.formula
    fterms, dterms = formula.formation, formula.dissolution
    pf, pd = length(fterms), length(dterms)

    n_sim >= 1 || throw(ArgumentError("gof: n_sim must be at least 1 (got $n_sim)"))
    n_trans = length(nets) - 1
    n = Int(nv(nets[1]))
    directed = is_directed(model)
    modes = directed ? (:in, :out) : (:total,)

    # Observed: pooled tie changes, formation/persistence model statistics on
    # the observed auxiliary networks, and the degree distribution of Y_t
    obs_changes = zeros(2)
    obs_f = zeros(pf)
    obs_d = zeros(pd)
    obs_deg = zeros(n, length(modes))
    counts = zeros(n)
    for t in 2:length(nets)
        prev, curr = nets[t-1], nets[t]
        obs_changes[1] += compute(NewEdge(), curr, prev)
        obs_changes[2] += compute(PersistentEdge(), curr, prev)
        _transition_stats!(obs_f, obs_d, fterms, dterms, prev, curr)
        for (m, mode) in enumerate(modes)
            _degree_counts!(counts, curr, mode)
            obs_deg[:, m] .+= counts
        end
    end

    # One seed per (transition, simulation) drawn from `rng` up front, in a
    # fixed order, and each simulation run on its own `Xoshiro(seed)`: the
    # draws are reproducible from `rng` alone and independent of how many
    # threads execute them, so `gof(fit; n_sim=100)` is as fast as the
    # machine allows and gives the same answer on every machine.
    seeds = rand(rng, UInt64, n_trans * n_sim)
    changes = zeros(n_trans, n_sim, 2)
    sim_f = zeros(n_trans, n_sim, pf)
    sim_d = zeros(n_trans, n_sim, pd)
    sim_deg = zeros(n_trans, n_sim, n, length(modes))
    spawn_all(n_trans * n_sim) do k
        t, s = divrem(k - 1, n_sim) .+ 1
        prev = nets[t]
        sim = simulate_stergm(prev, formula, result.formation_coef,
                              result.persistence_coef; burnin=burnin,
                              rng=Random.Xoshiro(seeds[k]))
        changes[t, s, 1] = compute(NewEdge(), sim, prev)
        changes[t, s, 2] = compute(PersistentEdge(), sim, prev)
        _transition_stats!(view(sim_f, t, s, :), view(sim_d, t, s, :),
                           fterms, dterms, prev, sim)
        for (m, mode) in enumerate(modes)
            _degree_counts!(view(sim_deg, t, s, :, m), sim, mode)
        end
        return nothing
    end
    pooled(a) = dropdims(sum(a; dims=1); dims=1)   # (n_sim × levels)

    net = nets[1]
    panels = GOFStatistic[
        GOFStatistic("tie changes", ["formed", "persisted"], obs_changes, pooled(changes)),
        GOFStatistic("formation statistics", _labels(fterms, net), obs_f, pooled(sim_f)),
        GOFStatistic("persistence statistics", _labels(dterms, net), obs_d, pooled(sim_d)),
    ]
    deg_pooled = pooled(sim_deg)                    # (n_sim × n × modes)
    for (m, mode) in enumerate(modes)
        push!(panels, GOFStatistic(mode === :in ? "idegree" : mode === :out ? "odegree" : "degree",
                                   string.(0:n-1), obs_deg[:, m], deg_pooled[:, :, m]))
    end
    return GOFResult(panels; model="STERGM")
end

# The formation and persistence model statistics of one transition
# Y_{t−1} → Y_t, accumulated into `acc_f`/`acc_d`: every formation term on
# Y⁺ = Y_{t−1} ∪ Y_t and every persistence term on Y⁻ = Y_{t−1} ∩ Y_t (a
# temporal term also sees Y_{t−1}) — the statistics the CMPLE's formation
# and persistence models are fit to.
function _transition_stats!(acc_f, acc_d, fterms, dterms, prev, curr)
    yplus = formation_network(prev, curr)
    yminus = dissolution_network(prev, curr)
    for (k, term) in enumerate(_materialized_tuple(fterms, yplus))
        acc_f[k] += _tcompute(term, yplus, prev)
    end
    for (k, term) in enumerate(_materialized_tuple(dterms, yminus))
        acc_d[k] += _tcompute(term, yminus, prev)
    end
    return nothing
end

# Degree histogram of `net` over the bins 0 … n−1 (`counts` has length n),
# for in-, out- or total degree
function _degree_counts!(counts, net, mode::Symbol)
    fill!(counts, 0.0)
    for v in 1:Int(nv(net))
        d = mode === :in ? length(inneighbors(net, v)) :
            mode === :out ? length(outneighbors(net, v)) :
            length(neighbors(net, v))
        counts[d + 1] += 1.0
    end
    return nothing
end

# ----------------------------------------------------------------------------
# Precompile workload: the fit → show → bootstrap → simulate → gof path a
# first session takes, on a 6-actor 3-wave panel, so the method instances
# are cached in the package image instead of compiled at the user's first
# call. Warnings
# the tiny fits emit are routed to a devnull logger: a precompile-time warning
# is never about the user's data.
# ----------------------------------------------------------------------------
@setup_workload begin
    _pc_nets = Network{Int,true}[]
    for edges_t in (((1, 2), (3, 4), (2, 5), (5, 1)), ((1, 2), (2, 3), (3, 4), (5, 1)),
                    ((1, 2), (2, 3), (5, 6), (6, 1)))
        _pc_net = network(6; directed=true)
        set_vertex_attribute!(_pc_net, :grp,
            Dict(1 => "A", 2 => "A", 3 => "B", 4 => "B", 5 => "A", 6 => "B"))
        for (i, j) in edges_t
            add_edge!(_pc_net, i, j)
        end
        push!(_pc_nets, _pc_net)
    end
    _pc_null = Base.CoreLogging.ConsoleLogger(devnull, Base.CoreLogging.Warn)
    @compile_workload begin
        Base.CoreLogging.with_logger(_pc_null) do
            _pc_rng = Random.Xoshiro(20260912)
            _pc_fit = stergm(_pc_nets, [Edges(), NodeMatch(:grp)], [Edges(), Delrecip()])
            sprint(show, _pc_fit); coef(_pc_fit); stderror(_pc_fit); bic(_pc_fit)
            coeftable(_pc_fit); confint(_pc_fit); confint(_pc_fit; level=0.9)
            approximations(_pc_fit); is_exact(_pc_fit); fit_metadata(_pc_fit)
            stergm(_pc_nets, [Edges()], [Edges()]; se=:bootstrap, n_boot=2, rng=_pc_rng)
            _pc_dep = stergm(_pc_nets, [Edges(), Mutual()], [Edges()]; method=:cmple)
            sprint(show, _pc_dep); approximations(_pc_dep)
            # the default (`method=:auto`) of a dyad-dependent formula is the
            # CMLE: one short iteration compiles its path (a refusal of the
            # tiny panel's start is not an error of the workload)
            try
                sprint(show, stergm(_pc_nets, [Edges(), Mutual()], [Edges()];
                                    n_samples=20, maxiter=1, bridge_rungs=1,
                                    bridge_samples=10, rng=_pc_rng))
            catch err
                err isa ArgumentError || rethrow()
            end
            _pc_cmle = stergm(_pc_nets, [Edges(), Delrecip()], [Edges()]; method=:cmle)
            sprint(show, _pc_cmle)
            simulate_stergm(_pc_fit, 1; burnin=5, rng=_pc_rng)
            gof(_pc_fit; n_sim=2, burnin=5, rng=_pc_rng)
        end
    end
end

end # module
