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
using Network
using Random
using Statistics

import ERGM: name, compute, change_stat

# Model types
export STERGM, STERGMResult, STERGMModel

# Temporal terms (evaluated with reference to the previous network)
export EdgeStability, Delrecip, PersistentEdge, NewEdge

# Auxiliary-network construction
export formation_network, dissolution_network

# Estimation
export stergm, fit_stergm, cmple, cmle, egmme

# Simulation
export simulate_stergm, simulate_network_sequence

# Diagnostics
export stergm_gof

# Temporal descriptives
export edge_ages, mean_edge_age

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

function _copy_net(net::Network{T}) where T
    c = Network{T}(; n=Int(nv(net)), directed=is_directed(net))
    for e in edges(net)
        add_edge!(c, src(e), dst(e))
    end
    return c
end

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
    yminus = Network{Int}(; n=Int(nv(prev_net)), directed=is_directed(prev_net))
    for e in edges(prev_net)
        if has_edge(curr_net, src(e), dst(e))
            add_edge!(yminus, src(e), dst(e))
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

"""
    STERGMModel

A STERGM formula bound to an observed network panel.
"""
struct STERGMModel{T}
    formula::STERGM
    networks::Vector{Network{T}}

    function STERGMModel(formula::STERGM, networks::Vector{Network{T}}) where T
        length(networks) >= 2 ||
            throw(ArgumentError("need at least two panels (one transition)"))
        n = nv(networks[1])
        dir = is_directed(networks[1])
        for net in networks
            (nv(net) == n && is_directed(net) == dir) ||
                throw(ArgumentError("all panels must share size and directedness"))
        end
        new{T}(formula, networks)
    end
end

"""
    STERGMResult

Fitted STERGM. Dissolution coefficients are **persistence** log-odds.
`loglik_*` are the maximized conditional pseudo-log-likelihoods.
"""
struct STERGMResult{T}
    model::STERGMModel{T}
    formation_coef::Vector{Float64}
    formation_se::Vector{Float64}
    dissolution_coef::Vector{Float64}
    dissolution_se::Vector{Float64}
    loglik_formation::Float64
    loglik_dissolution::Float64
    converged::Bool
    method::Symbol
end

function Base.show(io::IO, r::STERGMResult)
    println(io, "STERGM Results ($(r.method))")
    println(io, "="^40)
    println(io, "Panels: $(length(r.model.networks)); converged: $(r.converged)")
    println(io, "Pseudo-log-likelihood: formation $(round(r.loglik_formation, digits=3)), " *
                "dissolution $(round(r.loglik_dissolution, digits=3))")
    println(io)
    println(io, "Formation:")
    for (k, t) in enumerate(r.model.formula.formation)
        println(io, "  $(rpad(name(t), 22)) $(lpad(round(r.formation_coef[k], digits=4), 10)) " *
                    "(SE: $(round(r.formation_se[k], digits=4)))")
    end
    println(io, "Dissolution (persistence):")
    for (k, t) in enumerate(r.model.formula.dissolution)
        println(io, "  $(rpad(name(t), 22)) $(lpad(round(r.dissolution_coef[k], digits=4), 10)) " *
                    "(SE: $(round(r.dissolution_se[k], digits=4)))")
    end
end

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

# Weighted logistic Newton-Raphson with step-halving; returns
# (coef, se, loglik, converged)
function _logistic_fit(X::Matrix{Float64}, y::Vector{Bool};
                       maxiter::Int=100, tol::Float64=1e-8)
    n, p = size(X)
    n > 0 || return (fill(NaN, p), fill(NaN, p), NaN, false)

    function derivatives(β)
        ll = 0.0
        grad = zeros(p)
        hess = zeros(p, p)
        for r in 1:n
            η = dot(β, @view X[r, :])
            pr = 1.0 / (1.0 + exp(-η))
            ll += y[r] ? (η < 0 ? η - log1p(exp(η)) : -log1p(exp(-η))) :
                         (η < 0 ? -log1p(exp(η)) : -η - log1p(exp(-η)))
            resid = (y[r] ? 1.0 : 0.0) - pr
            x = @view X[r, :]
            grad .+= resid .* x
            hess .-= (pr * (1 - pr)) .* (x * x')
        end
        return ll, grad, hess
    end

    β = zeros(p)
    ll, grad, hess = derivatives(β)
    converged = false

    for _ in 1:maxiter
        step = try
            -hess \ grad
        catch
            break
        end

        stepsize = 1.0
        ll_new, grad_new, hess_new = ll, grad, hess
        for _ in 1:10
            ll_new, grad_new, hess_new = derivatives(β .+ stepsize .* step)
            ll_new >= ll && break
            stepsize /= 2
        end

        β .+= stepsize .* step
        ll_change = abs(ll_new - ll)
        ll, grad, hess = ll_new, grad_new, hess_new

        if ll_change < tol && norm(grad) < sqrt(tol)
            converged = true
            break
        end
    end

    se = try
        sqrt.(abs.(diag(pinv(-hess))))
    catch
        fill(NaN, p)
    end

    return (β, se, ll, converged)
end

"""
    cmple(model::STERGMModel; maxiter=100, tol=1e-8) -> STERGMResult

Conditional maximum pseudo-likelihood: for every transition t,

- **formation** rows are the non-edges of Y_{t-1}; the response is edge
  presence in Y_t and change statistics are evaluated on the formation
  network Y⁺ = Y_{t-1} ∪ Y_t;
- **dissolution** rows are the edges of Y_{t-1}; the response is
  persistence into Y_t and change statistics are evaluated on the
  dissolution network Y⁻ = Y_{t-1} ∩ Y_t.

Rows pool across transitions; the two logistic likelihoods are maximized
independently (separability). Standard errors are pseudo-likelihood SEs —
anticonservative for dyad-dependent terms.
"""
function cmple(model::STERGMModel{T}; maxiter::Int=100, tol::Float64=1e-8) where T
    fterms = model.formula.formation
    dterms = model.formula.dissolution
    pf, pd = length(fterms), length(dterms)
    directed = is_directed(model.networks[1])
    n = Int(nv(model.networks[1]))

    Xf_rows = Vector{Vector{Float64}}()
    yf = Bool[]
    Xd_rows = Vector{Vector{Float64}}()
    yd = Bool[]

    for t in 2:length(model.networks)
        prev = model.networks[t-1]
        curr = model.networks[t]
        yplus = formation_network(prev, curr)
        yminus = dissolution_network(prev, curr)

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
    end

    Xf = isempty(Xf_rows) ? Matrix{Float64}(undef, 0, pf) :
         permutedims(reduce(hcat, Xf_rows))
    Xd = isempty(Xd_rows) ? Matrix{Float64}(undef, 0, pd) :
         permutedims(reduce(hcat, Xd_rows))

    fcoef, fse, fll, fconv = _logistic_fit(Xf, yf; maxiter=maxiter, tol=tol)
    dcoef, dse, dll, dconv = _logistic_fit(Xd, yd; maxiter=maxiter, tol=tol)

    return STERGMResult(model, fcoef, fse, dcoef, dse, fll, dll,
                        fconv && dconv, :cmple)
end

"""
    cmle(model::STERGMModel; kwargs...) -> STERGMResult

Alias for [`cmple`](@ref). True MCMC-based CMLE is **not implemented**;
CMPLE coincides with it only for dyad-independent terms, so a one-time
warning is emitted.
"""
function cmle(model::STERGMModel; kwargs...)
    @warn "cmle: MCMC-based CMLE is not implemented; falling back to CMPLE " *
          "(exact only for dyad-independent terms)" maxlog = 1
    result = cmple(model; kwargs...)
    return STERGMResult(result.model, result.formation_coef, result.formation_se,
                        result.dissolution_coef, result.dissolution_se,
                        result.loglik_formation, result.loglik_dissolution,
                        result.converged, :cmle)
end

"""
    egmme(model::STERGMModel; kwargs...)

Equilibrium Generalized Method of Moments estimation is **not
implemented**; this raises an error rather than returning placeholder
coefficients. Use CMPLE on panel data.
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

    # Y_t = new formations ∪ persisted ties
    yt = Network{Int}(; n=Int(nv(prev_net)), directed=is_directed(prev_net))
    for e in edges(yplus)
        has_edge(prev_net, src(e), dst(e)) || add_edge!(yt, src(e), dst(e))
    end
    for e in edges(yminus)
        add_edge!(yt, src(e), dst(e))
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
    stergm_gof(result::STERGMResult; n_sim=50, rng=...) -> NamedTuple

Transition-level goodness of fit: for each observed transition, simulate
`n_sim` STERGM transitions from the same starting panel and compare the
observed number of formed and persisted ties to their simulated
distributions (two-sided Monte Carlo p-values, pooled over transitions).
"""
function stergm_gof(result::STERGMResult; n_sim::Int=50,
                    rng::Random.AbstractRNG=Random.default_rng())
    model = result.model
    nets = model.networks
    formula = model.formula

    obs_formed = 0.0
    obs_persisted = 0.0
    sim_formed = zeros(n_sim)
    sim_persisted = zeros(n_sim)

    for t in 2:length(nets)
        prev, curr = nets[t-1], nets[t]
        obs_formed += compute(NewEdge(), curr, prev)
        obs_persisted += compute(PersistentEdge(), curr, prev)

        for s in 1:n_sim
            sim = simulate_stergm(prev, formula, result.formation_coef,
                                  result.dissolution_coef; rng=rng)
            sim_formed[s] += compute(NewEdge(), sim, prev)
            sim_persisted[s] += compute(PersistentEdge(), sim, prev)
        end
    end

    mc_p = (sims, obs) -> min(1.0, 2.0 * min(mean(sims .>= obs), mean(sims .<= obs)))

    return (
        observed = (formed = obs_formed, persisted = obs_persisted),
        simulated_mean = (formed = mean(sim_formed), persisted = mean(sim_persisted)),
        simulated_sd = (formed = std(sim_formed), persisted = std(sim_persisted)),
        p_values = (formed = mc_p(sim_formed, obs_formed),
                    persisted = mc_p(sim_persisted, obs_persisted)),
        n_sim = n_sim
    )
end

end # module
