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
  the free dyads of the auxiliary networks, pooled across transitions.
  This equals the CMLE exactly for dyad-independent terms and is an
  approximation otherwise (MCMC-based CMLE is not implemented);
- simulation draws Y⁺ and Y⁻ with Metropolis samplers constrained to the
  free dyads and combines them into Y_t.

The dissolution model is parameterized in terms of **persistence**:
positive coefficients mean ties are more likely to persist.
"""
module TERGM

using Distributions
using ERGM
using Graphs
using LinearAlgebra
using Networks
using Random
using Statistics

import ERGM: name, compute, change_stat, is_dyad_dependent, newton_fit,
              logistic_derivatives
import StatsAPI
import StatsAPI: coef, stderror, vcov, loglikelihood, aic, bic, nobs, dof

# `gof` extends the ONE shared Networks.jl generic (every model package adds
# methods for its own result types), so `gof(fit)` works uniformly across the
# ecosystem and loading several model packages never collides on the name.
import Networks: gof

# The shared result-metadata protocol (Networks.jl `src/results.jl`): the seven
# generic accessors that say what a fit actually did. Imported by name because
# TERGM adds methods for `STERGMResult`; `fit_metadata(fit)` collects them.
import Networks: estimand, objective, is_exact, se_method, missing_method,
                 approximations

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

# Simulation
export simulate_stergm, simulate_network_sequence

# Diagnostics (`gof` is Networks.jl's shared generic, extended with a method
# for STERGMResult; `stergm_gof` is kept as a backward-compatible alias)
export gof, stergm_gof

# Temporal descriptives
export edge_ages, mean_edge_age

# StatsAPI methods (re-exported so `coef(fit)` etc. work with just `using TERGM`)
export coef, stderror, vcov, loglikelihood, aic, bic, nobs, dof

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

compute(t::TemporalTerm, net) =
    error("$(name(t)) is a temporal term; call compute(term, net, prev_net)")
change_stat(t::TemporalTerm, net, i::Int, j::Int) =
    error("$(name(t)) is a temporal term; call change_stat(term, net, i, j, prev_net)")

# Evaluate any term's change stat on a transition
_tchange(t::TemporalTerm, net, i, j, prev) = change_stat(t, net, i, j, prev)
_tchange(t::AbstractERGMTerm, net, i, j, prev) = change_stat(t, net, i, j)

_tcompute(t::TemporalTerm, net, prev) = compute(t, net, prev)
_tcompute(t::AbstractERGMTerm, net, prev) = compute(t, net)

"""
    EdgeStability <: TemporalTerm

The number of dyads whose state agrees with the previous network:
`Σ_{ij} 1[y_ij = y^{t-1}_ij]`.
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

Delayed reciprocity: the number of ordered dyads (i, j) with `i→j` now and
`j→i` in the previous network. Directed networks only.
"""
struct Delrecip <: TemporalTerm end

name(::Delrecip) = "delrecip"

function compute(::Delrecip, net, prev_net)
    is_directed(net) || return 0.0
    total = 0.0
    for e in edges(net)
        has_edge(prev_net, dst(e), src(e)) && (total += 1.0)
    end
    return total
end

function change_stat(::Delrecip, net, i::Int, j::Int, prev_net)
    is_directed(net) || return 0.0
    return has_edge(prev_net, j, i) ? 1.0 : 0.0
end

"""
    PersistentEdge <: TemporalTerm

The number of edges present in both the current and the previous network.
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

The number of edges present now but absent in the previous network.
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

# Delayed reciprocity is only meaningful on directed networks; extend
# ERGM.jl's direction trait so model construction rejects it on undirected
# panels (like `Mutual`).
ERGM._requires_directed(::Delrecip) = true

# =============================================================================
# Temporal descriptives
# =============================================================================

"""
    edge_ages(networks::Vector{<:Network}) -> Dict{Tuple{Int,Int}, Int}

For the last network of the sequence, the age of each current edge: the
number of consecutive panels (ending at the last one) in which it has
been present.
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
evaluated on Y⁺; its free dyads are the non-edges of Y_{t-1}.
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
Y_{t-1}).
"""
struct STERGM
    formation::Vector{AbstractERGMTerm}
    dissolution::Vector{AbstractERGMTerm}

    function STERGM(formation::Vector{<:AbstractERGMTerm},
                    dissolution::Vector{<:AbstractERGMTerm})
        (isempty(formation) || isempty(dissolution)) &&
            throw(ArgumentError("formation and dissolution models must be non-empty"))
        new(collect(AbstractERGMTerm, formation),
            collect(AbstractERGMTerm, dissolution))
    end
end

# Formula validation at model construction, mirroring ERGM.jl's
# `ERGMModel` checks (its `_validate_formula` helper is internal, so the
# two checks are replicated here on top of the term traits ERGM.jl defines
# as extension points): attribute-based terms must reference a vertex
# attribute that exists on every panel, and direction-incompatible terms
# (`Mutual`, `Delrecip`, ...) are rejected on undirected panels. Failing
# loudly here prevents silently wrong fits from all-zero design columns.
function _validate_terms(terms::Vector{AbstractERGMTerm}, networks, side::String)
    for term in terms
        attr = ERGM._vertex_attribute(term)
        if attr !== nothing
            for (t, net) in enumerate(networks)
                available = list_vertex_attributes(net)
                if !(attr in available)
                    avail_str = isempty(available) ? "(none)" :
                        join(sort!([":" * String(a) for a in available]), ", ")
                    throw(ArgumentError(
                        "$side term '$(name(term))' refers to vertex attribute " *
                        ":$attr, which does not exist on panel $t. Available " *
                        "vertex attributes: $avail_str."))
                end
            end
        end
        if ERGM._requires_directed(term) && !is_directed(networks[1])
            throw(ArgumentError(
                "$side term '$(name(term))' is only defined for directed " *
                "networks, but the panels are undirected. Remove the term or " *
                "use directed panels."))
        end
    end
    return nothing
end

"""
    STERGMModel

A STERGM formula bound to an observed network panel.

Construction validates the formula against the panels — every
attribute-based term's vertex attribute must exist on each panel, and
direction-incompatible terms (`Mutual`, `Delrecip`, ...) are rejected on
undirected panels — throwing an `ArgumentError` otherwise (mirroring
`ERGM.ERGMModel`).
"""
struct STERGMModel{T}
    formula::STERGM
    networks::Vector{Network{T}}

    function STERGMModel(formula::STERGM, networks::Vector{<:Network{T}}) where T
        length(networks) >= 2 ||
            throw(ArgumentError("need at least two panels (one transition)"))
        n = nv(networks[1])
        dir = is_directed(networks[1])
        for net in networks
            (nv(net) == n && is_directed(net) == dir) ||
                throw(ArgumentError("all panels must share size and directedness"))
        end
        # Panel missingness is not implemented: CMPLE enumerates every free
        # dyad of Y⁺/Y⁻ as observed, so a masked dyad would silently enter the
        # design matrix at its face value. Reject rather than invent data.
        for (t, net) in enumerate(networks)
            require_observed(net; context="TERGM (panel $t)", face_ok=false)
        end
        _validate_terms(formula.formation, networks, "formation")
        _validate_terms(formula.dissolution, networks, "dissolution")
        new{T}(formula, networks)
    end
end

"""
    STERGMResult

Fitted STERGM. Dissolution coefficients are **persistence** log-odds.
`loglik_*` are the maximized conditional pseudo-log-likelihoods. `vcov` is
the joint covariance of the stacked (formation, then dissolution)
coefficients — block-diagonal by separability for `se_type == :hessian`,
the empirical bootstrap covariance for `se_type == :bootstrap`. `se_type`
records how the standard errors were obtained (`:hessian` for the inverse
observed pseudo-likelihood information, `:bootstrap` for the
per-transition block bootstrap).
"""
struct STERGMResult{T}
    model::STERGMModel{T}
    formation_coef::Vector{Float64}
    formation_se::Vector{Float64}
    dissolution_coef::Vector{Float64}
    dissolution_se::Vector{Float64}
    vcov::Matrix{Float64}
    loglik_formation::Float64
    loglik_dissolution::Float64
    converged::Bool
    method::Symbol
    se_type::Symbol
end

# Backward-compatible constructor from before `se_type` existed
STERGMResult(model, formation_coef, formation_se, dissolution_coef,
             dissolution_se, vcov, loglik_formation, loglik_dissolution,
             converged, method::Symbol) =
    STERGMResult(model, formation_coef, formation_se, dissolution_coef,
                 dissolution_se, vcov, loglik_formation, loglik_dissolution,
                 converged, method, :hessian)

"""
    _has_dyad_dependent(model::STERGMModel) -> Bool

Whether **either** formula (formation or dissolution) contains a dyad-dependent
term (`ERGM.is_dyad_dependent`; the temporal terms are dyad-independent — they
condition only on the exogenous previous network).

This is THE predicate that decides whether CMPLE is an approximation: when both
formulas are dyad-independent the conditional pseudo-likelihood *is* the
conditional likelihood, so CMPLE is the exact CMLE and the inverse-Hessian
standard errors are the exact ML ones. Defined once and used by both
`show(::STERGMResult)` (the prose caveat) and `is_exact(::STERGMResult)` (the
machine-readable answer), so the two cannot drift apart.
"""
_has_dyad_dependent(model::STERGMModel) =
    any(is_dyad_dependent(t)
        for t in Iterators.flatten((model.formula.formation,
                                    model.formula.dissolution)))

function Base.show(io::IO, r::STERGMResult)
    println(io, "STERGM Results ($(r.method))")
    println(io, "="^40)
    println(io, "Panels: $(length(r.model.networks)); converged: $(r.converged)")
    println(io, "Pseudo-log-likelihood: formation $(round(r.loglik_formation, digits=3)), " *
                "dissolution $(round(r.loglik_dissolution, digits=3))")
    println(io)

    # Shared ecosystem presentation layer (Networks.print_coeftable): one
    # R-style table per model, with the significance-code legend printed
    # once below the second table
    fz = r.formation_coef ./ r.formation_se
    dz = r.dissolution_coef ./ r.dissolution_se
    println(io, "Formation:")
    print_coeftable(io, name.(r.model.formula.formation), r.formation_coef,
                    r.formation_se, ERGM._z_pvalues(fz);
                    z_values=fz, legend=false)
    println(io)
    println(io, "Dissolution (persistence):")
    print_coeftable(io, name.(r.model.formula.dissolution), r.dissolution_coef,
                    r.dissolution_se, ERGM._z_pvalues(dz); z_values=dz)

    # Honest-uncertainty caveat (mirroring ERGM.jl's show): CMPLE fits of
    # dyad-dependent formulas have suspect inverse-Hessian standard errors.
    # Dyad-independent formulas need no caveat — there CMPLE is the CMLE.
    if r.method in (:cmple, :cmle) && _has_dyad_dependent(r.model)
        println(io)
        if r.se_type == :bootstrap
            println(io, "Note: this model contains dyad-dependent terms and was fit by")
            println(io, "conditional maximum pseudolikelihood (CMPLE). Standard errors are")
            println(io, "per-transition block-bootstrap estimates; the CMPLE point estimates")
            println(io, "may still be biased.")
        else
            println(io, "Warning: this model contains dyad-dependent terms and was fit by")
            println(io, "conditional maximum pseudolikelihood (CMPLE). The standard errors")
            println(io, "are based on the naive pseudolikelihood and are suspect (typically")
            println(io, "anticonservative); the p-values should not be trusted. Use")
            println(io, "se=:bootstrap for per-transition block-bootstrap standard errors.")
        end
    end
end

# ============================================================================
# The shared result-metadata protocol (Networks.jl `src/results.jl`)
# ============================================================================
#
# `fit_metadata(fit)` collects these accessors. They are derived from the SAME
# `_has_dyad_dependent` predicate the prose caveat in `show` uses, so the
# machine-readable answer and the printed one can never disagree.

estimand(::STERGMResult) = :stergm

"""
    objective(r::STERGMResult) -> Symbol

`:conditional_pseudolikelihood` — the pooled logistic pseudo-likelihood over the
free dyads of the formation (Y⁺) and dissolution (Y⁻) auxiliary networks.
"""
objective(::STERGMResult) = :conditional_pseudolikelihood

"""
    is_exact(r::STERGMResult) -> Bool

`true` iff **every** term in **both** formulas is dyad-independent: there the
conditional pseudo-likelihood is the conditional likelihood, so CMPLE is the
exact CMLE. Add one dyad-dependent term to either formula and the same estimator
reports `false` — the approximation is a property of the fit, not of the method
name.
"""
is_exact(r::STERGMResult) =
    r.method in (:cmple, :cmle) && !_has_dyad_dependent(r.model)

se_method(r::STERGMResult) = r.se_type

# `STERGMModel` calls `require_observed` on every panel with the default
# `:error` policy: CMPLE would enumerate a masked dyad at its face value, so
# masked panels are refused outright rather than silently used.
missing_method(::STERGMResult) = :rejected

function approximations(r::STERGMResult)
    out = String[]
    if _has_dyad_dependent(r.model)
        push!(out, "conditional maximum pseudo-likelihood of a dyad-dependent " *
                   "formula: the dyad conditionals of each transition are " *
                   "multiplied as if independent, so the point estimates are " *
                   "biased in finite samples")
        r.se_type === :hessian &&
            push!(out, "inverse-Hessian standard errors of the naive " *
                       "conditional pseudo-likelihood: expected anticonservative " *
                       "under dyadic dependence")
    end
    return out
end

# StatsAPI interface: methods on the shared statistics generics (mirroring
# ERGM.jl), so fitted STERGMs interoperate with StatsBase/GLM-style tooling.
# Coefficient-shaped accessors stack formation first, then dissolution.

# Total number of free dyads pooled over transitions: every off-diagonal
# dyad of each transition is free in exactly one of the two models
function _n_free_dyads(model::STERGMModel)
    n = Int(nv(model.networks[1]))
    per = is_directed(model.networks[1]) ? n * (n - 1) : n * (n - 1) ÷ 2
    return (length(model.networks) - 1) * per
end

StatsAPI.coef(r::STERGMResult) = vcat(r.formation_coef, r.dissolution_coef)
StatsAPI.stderror(r::STERGMResult) = vcat(r.formation_se, r.dissolution_se)
StatsAPI.vcov(r::STERGMResult) = r.vcov
StatsAPI.loglikelihood(r::STERGMResult) = r.loglik_formation + r.loglik_dissolution
StatsAPI.aic(r::STERGMResult) = -2 * loglikelihood(r) + 2 * dof(r)
StatsAPI.bic(r::STERGMResult) = -2 * loglikelihood(r) + dof(r) * log(nobs(r))
StatsAPI.nobs(r::STERGMResult) = _n_free_dyads(r.model)
StatsAPI.dof(r::STERGMResult) =
    length(r.formation_coef) + length(r.dissolution_coef)

# =============================================================================
# Estimation
# =============================================================================

"""
    stergm(networks, formation, dissolution; method=:cmple, kwargs...) -> STERGMResult

Fit a STERGM to a panel of networks.

`method`:
- `:cmple` (default) — conditional maximum pseudo-likelihood on the
  formation/dissolution networks. Exact CMLE for dyad-independent terms;
  an approximation for dyad-dependent terms (`Triangle`, ...).
- `:cmle` — alias for `:cmple` with a one-time warning; MCMC-based exact
  CMLE is not implemented.

Remaining keyword arguments are forwarded to [`cmple`](@ref) — in
particular `se=:bootstrap` for per-transition block-bootstrap standard
errors (with `n_boot` and `rng`).
"""
function stergm(networks::Vector{<:Network},
                formation::Vector{<:AbstractERGMTerm},
                dissolution::Vector{<:AbstractERGMTerm};
                method::Symbol=:cmple, kwargs...)
    model = STERGMModel(STERGM(formation, dissolution), networks)

    if method == :cmple
        return cmple(model; kwargs...)
    elseif method == :cmle
        return cmle(model; kwargs...)
    elseif method == :egmme
        return egmme(model; kwargs...)
    else
        throw(ArgumentError("Unknown method: $method"))
    end
end

const fit_stergm = stergm

# Logistic pseudo-likelihood fit via the shared ERGM.newton_fit optimizer
# (Newton–Raphson with step halving); returns (coef, se, vcov, loglik,
# converged)
function _logistic_fit(X::Matrix{Float64}, y::Vector{Bool};
                       maxiter::Int=100, tol::Float64=1e-8)
    n, p = size(X)
    n > 0 || return (fill(NaN, p), fill(NaN, p), fill(NaN, p, p), NaN, false)

    # The derivatives come from the shared `ERGM.logistic_derivatives` (review
    # finding 15): the CMPLE over the free dyads IS a logistic likelihood, and
    # its derivatives are gemv/gemm over the whole design — η = Xβ, ∇ = X'r,
    # −H = X'WX — rather than a per-row `x * x'` outer product allocating a p×p
    # matrix on every one of the n rows of every Newton evaluation. Never paste
    # the loop back in; ERGMMulti and ERGMRank run on the same one.
    fit = newton_fit(logistic_derivatives(X, y), zeros(p); maxiter=maxiter,
                     tol=tol)
    return (fit.θ, fit.se, fit.vcov, fit.loglik, fit.converged)
end

# Per-transition CMPLE design blocks: for transition t (2:T), the
# formation rows (prior non-edges, stats on Y⁺) and dissolution rows
# (prior edges, stats on Y⁻). Kept per transition so the block bootstrap
# can resample whole transitions.
function _cmple_blocks(model::STERGMModel)
    fterms = model.formula.formation
    dterms = model.formula.dissolution
    pf, pd = length(fterms), length(dterms)
    directed = is_directed(model.networks[1])
    n = Int(nv(model.networks[1]))

    Xf_blocks = Matrix{Float64}[]
    yf_blocks = Vector{Bool}[]
    Xd_blocks = Matrix{Float64}[]
    yd_blocks = Vector{Bool}[]

    for t in 2:length(model.networks)
        prev = model.networks[t-1]
        curr = model.networks[t]
        yplus = formation_network(prev, curr)
        yminus = dissolution_network(prev, curr)

        Xf_rows = Vector{Vector{Float64}}()
        yf = Bool[]
        Xd_rows = Vector{Vector{Float64}}()
        yd = Bool[]
        for i in 1:n
            for j in (directed ? (1:n) : (i+1:n))
                i == j && continue
                if !has_edge(prev, i, j)
                    # Formation-free dyad: stats on Y⁺
                    push!(Xf_rows, [_tchange(term, yplus, i, j, prev)
                                    for term in fterms])
                    push!(yf, has_edge(curr, i, j))
                else
                    # Dissolution-free dyad: stats on Y⁻
                    push!(Xd_rows, [_tchange(term, yminus, i, j, prev)
                                    for term in dterms])
                    push!(yd, has_edge(curr, i, j))
                end
            end
        end

        push!(Xf_blocks, isempty(Xf_rows) ? Matrix{Float64}(undef, 0, pf) :
                         permutedims(reduce(hcat, Xf_rows)))
        push!(yf_blocks, yf)
        push!(Xd_blocks, isempty(Xd_rows) ? Matrix{Float64}(undef, 0, pd) :
                         permutedims(reduce(hcat, Xd_rows)))
        push!(yd_blocks, yd)
    end

    return Xf_blocks, yf_blocks, Xd_blocks, yd_blocks
end

_stack(blocks::Vector{Matrix{Float64}}, p::Int) =
    isempty(blocks) ? Matrix{Float64}(undef, 0, p) : reduce(vcat, blocks)
_stack(blocks::Vector{Vector{Bool}}) =
    isempty(blocks) ? Bool[] : reduce(vcat, blocks)

# Per-transition block bootstrap of the CMPLE (the btergm approach):
# resample whole time-transitions with replacement, refit both logistic
# pseudo-likelihoods on the resampled rows, and return the empirical
# covariance of the stacked (formation, dissolution) coefficients.
# Replicates with non-finite coefficients (e.g. separation in a resample)
# are dropped with a warning.
function _cmple_block_bootstrap(Xf_blocks, yf_blocks, Xd_blocks, yd_blocks,
                                pf::Int, pd::Int;
                                n_boot::Int, maxiter::Int, tol::Float64,
                                rng::Random.AbstractRNG)
    n_trans = length(Xf_blocks)
    n_trans >= 2 ||
        throw(ArgumentError("se=:bootstrap resamples time-transitions and " *
                            "needs at least 3 panels (2 transitions); got " *
                            "$(n_trans) transition(s)"))
    n_boot >= 2 || throw(ArgumentError("n_boot must be at least 2"))

    boot = Matrix{Float64}(undef, n_boot, pf + pd)
    n_ok = 0
    for _ in 1:n_boot
        idx = rand(rng, 1:n_trans, n_trans)
        fcoef, _, _, _, _ = _logistic_fit(_stack(Xf_blocks[idx], pf),
                                          _stack(yf_blocks[idx]);
                                          maxiter=maxiter, tol=tol)
        dcoef, _, _, _, _ = _logistic_fit(_stack(Xd_blocks[idx], pd),
                                          _stack(yd_blocks[idx]);
                                          maxiter=maxiter, tol=tol)
        θb = vcat(fcoef, dcoef)
        all(isfinite, θb) || continue
        n_ok += 1
        boot[n_ok, :] = θb
    end

    n_ok >= 2 ||
        error("block bootstrap failed: fewer than 2 of $n_boot replicates " *
              "produced finite coefficients")
    n_ok < n_boot &&
        @warn "block bootstrap: dropped $(n_boot - n_ok) of $n_boot " *
              "replicates with non-finite coefficients"

    return Matrix(cov(@view boot[1:n_ok, :]))
end

"""
    cmple(model::STERGMModel; maxiter=100, tol=1e-8, se=:hessian,
          n_boot=200, rng=Random.default_rng()) -> STERGMResult

Conditional maximum pseudo-likelihood: for every transition t,

- **formation** rows are the non-edges of Y_{t-1}; the response is edge
  presence in Y_t and change statistics are evaluated on the formation
  network Y⁺ = Y_{t-1} ∪ Y_t;
- **dissolution** rows are the edges of Y_{t-1}; the response is
  persistence into Y_t and change statistics are evaluated on the
  dissolution network Y⁻ = Y_{t-1} ∩ Y_t.

Rows pool across transitions; the two logistic likelihoods are maximized
independently (separability).

`se` selects the standard errors:
- `:hessian` (default) — inverse observed pseudo-likelihood information.
  **Caution:** anticonservative for dyad-dependent terms.
- `:bootstrap` — per-transition block bootstrap (the `btergm` approach):
  resample the time-transitions with replacement `n_boot` times, refit the
  CMPLE on each resample, and use the empirical covariance of the refitted
  stacked coefficients (needs at least 2 transitions, i.e. 3 panels).
  `rng` seeds the resampling. Point estimates are unchanged; `vcov` becomes
  the full joint bootstrap covariance (no longer block-diagonal).
"""
function cmple(model::STERGMModel{T}; maxiter::Int=100, tol::Float64=1e-8,
               se::Symbol=:hessian, n_boot::Int=200,
               rng::Random.AbstractRNG=Random.default_rng()) where T
    se in (:hessian, :bootstrap) ||
        throw(ArgumentError("se must be :hessian or :bootstrap, got :$se"))

    fterms = model.formula.formation
    dterms = model.formula.dissolution
    pf, pd = length(fterms), length(dterms)

    Xf_blocks, yf_blocks, Xd_blocks, yd_blocks = _cmple_blocks(model)

    fcoef, fse, fvcov, fll, fconv = _logistic_fit(_stack(Xf_blocks, pf),
                                                  _stack(yf_blocks);
                                                  maxiter=maxiter, tol=tol)
    dcoef, dse, dvcov, dll, dconv = _logistic_fit(_stack(Xd_blocks, pd),
                                                  _stack(yd_blocks);
                                                  maxiter=maxiter, tol=tol)

    if se == :bootstrap
        vcov_joint = _cmple_block_bootstrap(Xf_blocks, yf_blocks,
                                            Xd_blocks, yd_blocks, pf, pd;
                                            n_boot=n_boot, maxiter=maxiter,
                                            tol=tol, rng=rng)
        ses = sqrt.(max.(diag(vcov_joint), 0.0))
        fse, dse = ses[1:pf], ses[pf+1:end]
    else
        # Joint covariance of the stacked coefficients: block-diagonal,
        # since the formation and dissolution likelihoods are maximized
        # independently
        vcov_joint = [fvcov zeros(pf, pd); zeros(pd, pf) dvcov]
    end

    return STERGMResult(model, fcoef, fse, dcoef, dse, vcov_joint, fll, dll,
                        fconv && dconv, :cmple, se)
end

"""
    cmle(model::STERGMModel; kwargs...)

Conditional maximum likelihood (MCMC-based) is **not implemented** and this
function **throws**.

It previously warned and silently returned a [`cmple`](@ref) fit labelled
`:cmle`. That is unsafe for research use: a warning is easy to miss in a script
that continues, and the returned object then misreports its own estimator.
CMPLE coincides with CMLE only for dyad-independent formulas; for dependent
terms it is a different estimator with anticonservative standard errors.

Call [`cmple`](@ref) explicitly (`method = :cmple`, the default) to fit by
conditional maximum *pseudo*-likelihood.
"""
function cmle(model::STERGMModel; kwargs...)
    throw(ArgumentError(
        "cmle: MCMC-based conditional maximum likelihood is not implemented in " *
        "TERGM.jl. It is not the same estimator as CMPLE except for " *
        "dyad-independent formulas, so this no longer falls back silently.\n" *
        "  • use `cmple(model)` (or `method = :cmple`, the default) to fit by " *
        "conditional maximum pseudo-likelihood; pass `se = :bootstrap` for " *
        "per-transition block-bootstrap standard errors."))
end

"""
    TERGM.egmme(model::STERGMModel; kwargs...)

Equilibrium Generalized Method of Moments estimation is **not implemented**
and this function is **not exported** — exporting it would advertise an
estimator that does not exist. It is retained, and reachable as
`TERGM.egmme`, only so that calling it raises an informative error rather
than an `UndefVarError`; `stergm(...; method = :egmme)` throws the same
error.

Use [`cmple`](@ref) (`method = :cmple`, the default) on panel data.
"""
function egmme(model::STERGMModel; kwargs...)
    error("EGMME is not implemented in TERGM.jl; fit panel data with " *
          "method = :cmple instead")
end

# =============================================================================
# Simulation
# =============================================================================

# Metropolis sampling of Y⁺ (constrain=:formation, free dyads = non-edges
# of prev) or Y⁻ (constrain=:dissolution, free dyads = edges of prev),
# starting from Y_{t-1}
function _sample_constrained(prev::Network, terms, θ::Vector{Float64},
                             constrain::Symbol, steps::Int,
                             rng::Random.AbstractRNG)
    net = _copy_net(prev)
    n = Int(nv(net))
    directed = is_directed(net)

    free = NTuple{2, Int}[]
    for i in 1:n
        for j in (directed ? (1:n) : (i+1:n))
            i == j && continue
            if constrain == :formation
                has_edge(prev, i, j) || push!(free, (i, j))
            else
                has_edge(prev, i, j) && push!(free, (i, j))
            end
        end
    end
    isempty(free) && return net

    delta = Vector{Float64}(undef, length(terms))
    for _ in 1:steps
        i, j = free[rand(rng, 1:length(free))]
        for (k, t) in enumerate(terms)
            delta[k] = _tchange(t, net, i, j, prev)
        end
        log_accept = dot(θ, delta)
        if has_edge(net, i, j)
            log_accept = -log_accept
        end
        if log(rand(rng)) < log_accept
            if has_edge(net, i, j)
                rem_edge!(net, i, j)
            else
                add_edge!(net, i, j)
            end
        end
    end

    return net
end

"""
    simulate_stergm(prev_net, formula::STERGM, θ_form, θ_diss;
                    burnin=..., rng=...) -> Network

Simulate one STERGM transition from `prev_net`: draw the formation
network Y⁺ (Metropolis over the non-edges of `prev_net`) and the
dissolution network Y⁻ (over its edges), then combine
`Y_t = (Y⁺ \\ Y_{t-1}) ∪ Y⁻`.
"""
function simulate_stergm(prev_net::Network, formula::STERGM,
                         θ_form::Vector{Float64}, θ_diss::Vector{Float64};
                         burnin::Int=3000,
                         rng::Random.AbstractRNG=Random.default_rng())
    length(θ_form) == length(formula.formation) ||
        throw(ArgumentError("need one formation coefficient per term"))
    length(θ_diss) == length(formula.dissolution) ||
        throw(ArgumentError("need one dissolution coefficient per term"))

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
fitted STERGM.
"""
function simulate_stergm(result::STERGMResult, n_steps::Int; kwargs...)
    return simulate_network_sequence(result.model.formula,
                                     result.model.networks[end], n_steps,
                                     result.formation_coef,
                                     result.dissolution_coef; kwargs...)
end

"""
    simulate_network_sequence(formula, init_net, n_steps, θ_form, θ_diss;
                              burnin=..., rng=...) -> Vector{Network}

Simulate a sequence of `n_steps` STERGM transitions starting from
`init_net` (returned sequence includes the initial network as its first
element).
"""
function simulate_network_sequence(formula::STERGM, init_net::Network,
                                   n_steps::Int,
                                   θ_form::Vector{Float64},
                                   θ_diss::Vector{Float64};
                                   burnin::Int=3000,
                                   rng::Random.AbstractRNG=Random.default_rng())
    seq = Network{Int}[_copy_net(init_net)]
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
    gof(result::STERGMResult; n_sim=50, rng=...) -> GOFResult

Transition-level goodness of fit: for each observed transition, simulate
`n_sim` STERGM transitions from the same starting panel and compare the
observed number of formed and persisted ties (pooled over transitions) to
their simulated distributions.

Extends Networks.jl's shared `gof` generic and returns the shared
`Networks.GOFResult` container: one `GOFStatistic` panel named
`"tie changes"` with levels `"formed"` and `"persisted"`, and two-sided
Monte-Carlo p-values computed with the `(1 + k)/(N + 1)` estimator (never
exactly zero). [`stergm_gof`](@ref) is a backward-compatible alias.
"""
function gof(result::STERGMResult; n_sim::Int=50,
             rng::Random.AbstractRNG=Random.default_rng())
    model = result.model
    nets = model.networks
    formula = model.formula

    obs = zeros(2)              # pooled (formed, persisted) tie counts
    sims = zeros(n_sim, 2)

    for t in 2:length(nets)
        prev, curr = nets[t-1], nets[t]
        obs[1] += compute(NewEdge(), curr, prev)
        obs[2] += compute(PersistentEdge(), curr, prev)

        for s in 1:n_sim
            sim = simulate_stergm(prev, formula, result.formation_coef,
                                  result.dissolution_coef; rng=rng)
            sims[s, 1] += compute(NewEdge(), sim, prev)
            sims[s, 2] += compute(PersistentEdge(), sim, prev)
        end
    end

    stat = GOFStatistic("tie changes", ["formed", "persisted"], obs, sims)
    return GOFResult(stat; model="STERGM")
end

"""
    stergm_gof(result::STERGMResult; n_sim=50, rng=...) -> GOFResult

Backward-compatible alias for [`gof`](@ref) (Networks.jl's shared generic).
"""
const stergm_gof = gof

end # module
