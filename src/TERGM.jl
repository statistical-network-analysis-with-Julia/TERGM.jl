"""
    TERGM.jl - Temporal Exponential Random Graph Models

Provides models for dynamic networks including:
- Separable Temporal ERGM (STERGM) for discrete-time networks
- Formation and dissolution models
- Temporal terms (stability, memory, etc.)

Port of the R tergm package from the StatNet collection.
"""
module TERGM

using Distributions
using ERGM
using Graphs
using LinearAlgebra
using Network
using NetworkDynamic
using Optim
using Random
using Statistics
using StatsBase

# Model types
export STERGM, STERGMResult, STERGMModel
export FormationModel, DissolutionModel

# Temporal terms
export EdgeAge, EdgeStability, Memory, Delrecip
export FormationTerm, DissolutionTerm
export TimeLag, PersistentEdge, NewEdge

# Estimation
export stergm, fit_stergm
export egmme, cmle, cmple

# Simulation
export simulate_stergm, simulate_network_sequence

# Diagnostics
export stergm_gof

# =============================================================================
# Temporal Term Types
# =============================================================================

"""
    TemporalTerm <: AbstractERGMTerm

Base type for terms that reference previous network state.
"""
abstract type TemporalTerm <: AbstractERGMTerm end

"""
    FormationTerm{T<:AbstractERGMTerm}

Wrapper that marks a term as part of the formation model.
"""
struct FormationTerm{T<:AbstractERGMTerm} <: TemporalTerm
    term::T
end

name(ft::FormationTerm) = "Form(" * name(ft.term) * ")"

function compute(ft::FormationTerm, net, prev_net=nothing)
    return compute(ft.term, net)
end

function change_stat(ft::FormationTerm, net, i::Int, j::Int, prev_net=nothing)
    # Formation term only counts for edges not present in previous network
    if !isnothing(prev_net) && has_edge(prev_net, i, j)
        return 0.0
    end
    return change_stat(ft.term, net, i, j)
end

"""
    DissolutionTerm{T<:AbstractERGMTerm}

Wrapper that marks a term as part of the dissolution model.
"""
struct DissolutionTerm{T<:AbstractERGMTerm} <: TemporalTerm
    term::T
end

name(dt::DissolutionTerm) = "Diss(" * name(dt.term) * ")"

function compute(dt::DissolutionTerm, net, prev_net=nothing)
    return compute(dt.term, net)
end

function change_stat(dt::DissolutionTerm, net, i::Int, j::Int, prev_net=nothing)
    # Dissolution term only counts for edges present in previous network
    if !isnothing(prev_net) && !has_edge(prev_net, i, j)
        return 0.0
    end
    return change_stat(dt.term, net, i, j)
end

"""
    EdgeStability <: TemporalTerm

Term counting edges that persist from previous time step.
"""
struct EdgeStability <: TemporalTerm end

name(::EdgeStability) = "edgestability"

function compute(::EdgeStability, net, prev_net)
    isnothing(prev_net) && return 0.0
    count = 0
    for e in edges(net)
        if has_edge(prev_net, src(e), dst(e))
            count += 1
        end
    end
    return Float64(count)
end

function change_stat(::EdgeStability, net, i::Int, j::Int, prev_net)
    isnothing(prev_net) && return 0.0
    if !has_edge(prev_net, i, j)
        return 0.0  # New edge, doesn't contribute to stability
    end
    return has_edge(net, i, j) ? -1.0 : 1.0
end

"""
    Memory{T} <: TemporalTerm

Memory term capturing edge persistence with a decay parameter.
"""
struct Memory{T<:Real} <: TemporalTerm
    theta::T  # Decay parameter

    Memory(theta::T=1.0) where T<:Real = new{T}(theta)
end

name(m::Memory) = "memory.$(m.theta)"

function compute(m::Memory, net, prev_net)
    isnothing(prev_net) && return 0.0
    count = 0.0
    for e in edges(net)
        if has_edge(prev_net, src(e), dst(e))
            count += m.theta
        end
    end
    return count
end

"""
    EdgeAge <: TemporalTerm

Term based on age of edges (how many time steps they've persisted).
Requires network sequence.
"""
struct EdgeAge <: TemporalTerm
    max_age::Int
    EdgeAge(max_age::Int=Inf) = new(max_age)
end

name(ea::EdgeAge) = "edgeage"

"""
    Delrecip <: TemporalTerm

Delayed reciprocity: reciprocated edges based on previous time step.
For directed networks only.
"""
struct Delrecip <: TemporalTerm end

name(::Delrecip) = "delrecip"

function compute(::Delrecip, net, prev_net)
    isnothing(prev_net) && return 0.0
    !is_directed(net) && return 0.0

    count = 0
    for e in edges(net)
        # Current edge (i,j) reciprocates previous edge (j,i)
        if has_edge(prev_net, dst(e), src(e))
            count += 1
        end
    end
    return Float64(count)
end

function change_stat(::Delrecip, net, i::Int, j::Int, prev_net)
    isnothing(prev_net) && return 0.0
    if has_edge(prev_net, j, i)
        return has_edge(net, i, j) ? -1.0 : 1.0
    end
    return 0.0
end

"""
    TimeLag{T<:AbstractERGMTerm} <: TemporalTerm

Lagged version of a standard ERGM term computed on previous network.
"""
struct TimeLag{T<:AbstractERGMTerm} <: TemporalTerm
    term::T
    lag::Int

    TimeLag(term::T, lag::Int=1) where T<:AbstractERGMTerm = new{T}(term, lag)
end

name(tl::TimeLag) = "lag$(tl.lag)." * name(tl.term)

function compute(tl::TimeLag, net, prev_net)
    isnothing(prev_net) && return 0.0
    return compute(tl.term, prev_net)
end

"""
    PersistentEdge <: TemporalTerm

Counts edges that exist in both current and previous network.
"""
struct PersistentEdge <: TemporalTerm end

name(::PersistentEdge) = "persistent.edges"

function compute(::PersistentEdge, net, prev_net)
    isnothing(prev_net) && return 0.0
    count = 0
    for e in edges(net)
        if has_edge(prev_net, src(e), dst(e))
            count += 1
        end
    end
    return Float64(count)
end

function change_stat(::PersistentEdge, net, i::Int, j::Int, prev_net)
    isnothing(prev_net) && return 0.0
    if !has_edge(prev_net, i, j)
        return 0.0
    end
    return has_edge(net, i, j) ? -1.0 : 1.0
end

"""
    NewEdge <: TemporalTerm

Counts edges that exist in current but not previous network.
"""
struct NewEdge <: TemporalTerm end

name(::NewEdge) = "new.edges"

function compute(::NewEdge, net, prev_net)
    isnothing(prev_net) && return Float64(ne(net))
    count = 0
    for e in edges(net)
        if !has_edge(prev_net, src(e), dst(e))
            count += 1
        end
    end
    return Float64(count)
end

function change_stat(::NewEdge, net, i::Int, j::Int, prev_net)
    isnothing(prev_net) && return has_edge(net, i, j) ? -1.0 : 1.0
    if has_edge(prev_net, i, j)
        return 0.0
    end
    return has_edge(net, i, j) ? -1.0 : 1.0
end

# =============================================================================
# STERGM Model Types
# =============================================================================

"""
    STERGM

Separable Temporal ERGM specification.

Models formation and dissolution as separate processes:
- P(Y_t | Y_{t-1}) ∝ P+(Y_t | Y_{t-1}; θ+) × P-(Y_t | Y_{t-1}; θ-)

Where P+ governs edge formation and P- governs edge persistence (1 - dissolution).

# Fields
- `formation_terms::Vector{AbstractERGMTerm}`: Terms for formation model
- `dissolution_terms::Vector{AbstractERGMTerm}`: Terms for dissolution (persistence) model
- `constraints::Vector`: Optional constraints
"""
struct STERGM
    formation_terms::Vector{AbstractERGMTerm}
    dissolution_terms::Vector{AbstractERGMTerm}
    constraints::Vector{Any}

    function STERGM(formation::Vector{<:AbstractERGMTerm},
                    dissolution::Vector{<:AbstractERGMTerm};
                    constraints::Vector=Any[])
        new(formation, dissolution, constraints)
    end
end

"""
    STERGMModel{T}

STERGM model with observed network sequence.
"""
struct STERGMModel{T}
    formula::STERGM
    networks::Vector{Network{T}}
    directed::Bool
    times::Vector{Float64}

    function STERGMModel(formula::STERGM, networks::Vector{Network{T}};
                         times::Vector{Float64}=Float64.(1:length(networks))) where T
        length(networks) >= 2 || throw(ArgumentError("Need at least 2 time points"))
        all(nv(n) == nv(networks[1]) for n in networks) ||
            throw(ArgumentError("All networks must have same vertex count"))
        new{T}(formula, networks, is_directed(networks[1]), times)
    end
end

"""
    STERGMResult

Results from fitting a STERGM.
"""
struct STERGMResult{T}
    model::STERGMModel{T}
    formation_coef::Vector{Float64}
    dissolution_coef::Vector{Float64}
    formation_se::Vector{Float64}
    dissolution_se::Vector{Float64}
    loglik::Float64
    method::Symbol
    converged::Bool
end

function Base.show(io::IO, result::STERGMResult)
    println(io, "STERGM Results")
    println(io, "==============")
    println(io, "Method: $(result.method)")
    println(io, "Converged: $(result.converged)")
    println(io, "Log-likelihood: $(round(result.loglik, digits=4))")
    println(io)

    println(io, "Formation Model:")
    println(io, "-"^40)
    for (i, term) in enumerate(result.model.formula.formation_terms)
        println(io, "  $(rpad(name(term), 20)) $(lpad(round(result.formation_coef[i], digits=4), 10)) " *
                    "(SE: $(round(result.formation_se[i], digits=4)))")
    end
    println(io)

    println(io, "Dissolution (Persistence) Model:")
    println(io, "-"^40)
    for (i, term) in enumerate(result.model.formula.dissolution_terms)
        println(io, "  $(rpad(name(term), 20)) $(lpad(round(result.dissolution_coef[i], digits=4), 10)) " *
                    "(SE: $(round(result.dissolution_se[i], digits=4)))")
    end
end

# =============================================================================
# Estimation Functions
# =============================================================================

"""
    stergm(networks, formation_terms, dissolution_terms; kwargs...) -> STERGMResult

Fit a Separable Temporal ERGM.

# Arguments
- `networks`: Vector of networks across time
- `formation_terms`: ERGM terms for edge formation
- `dissolution_terms`: ERGM terms for edge persistence

# Keyword Arguments
- `method::Symbol=:cmle`: Estimation method (:cmle, :cmple, :egmme)
- `times::Vector`: Time points corresponding to networks
"""
function stergm(networks::Vector{<:Network},
                formation_terms::Vector{<:AbstractERGMTerm},
                dissolution_terms::Vector{<:AbstractERGMTerm};
                method::Symbol=:cmle,
                times::Vector{Float64}=Float64.(1:length(networks)),
                kwargs...)

    formula = STERGM(formation_terms, dissolution_terms)
    model = STERGMModel(formula, networks; times=times)

    if method == :cmle
        return cmle(model; kwargs...)
    elseif method == :cmple
        return cmple(model; kwargs...)
    elseif method == :egmme
        return egmme(model; kwargs...)
    else
        throw(ArgumentError("Unknown method: $method. Use :cmle, :cmple, or :egmme"))
    end
end

fit_stergm = stergm

"""
    cmle(model::STERGMModel; kwargs...) -> STERGMResult

Conditional Maximum Likelihood Estimation for STERGM.
"""
function cmle(model::STERGMModel{T};
              maxiter::Int=100,
              tol::Float64=1e-6) where T

    n_form = length(model.formula.formation_terms)
    n_diss = length(model.formula.dissolution_terms)
    n_transitions = length(model.networks) - 1

    # Initialize coefficients
    form_coef = zeros(n_form)
    diss_coef = zeros(n_diss)

    # Build transition data for formation and dissolution
    # Formation: edges that could form (not in Y_{t-1}) and whether they did
    # Dissolution: edges that existed (in Y_{t-1}) and whether they persisted

    for iter in 1:maxiter
        form_grad = zeros(n_form)
        diss_grad = zeros(n_diss)
        form_hess = zeros(n_form, n_form)
        diss_hess = zeros(n_diss, n_diss)

        for t in 1:n_transitions
            prev_net = model.networks[t]
            curr_net = model.networks[t + 1]
            n = nv(prev_net)

            # Formation model: condition on non-edges at t-1
            for i in 1:n, j in (model.directed ? (1:n) : (i+1:n))
                i == j && continue
                has_edge(prev_net, i, j) && continue  # Only non-edges can form

                # Compute change statistics
                delta = [change_stat(term, curr_net, i, j) for term in model.formula.formation_terms]
                formed = has_edge(curr_net, i, j) ? 1.0 : 0.0

                # Logistic contribution
                eta = dot(form_coef, delta)
                prob = 1.0 / (1.0 + exp(-eta))

                form_grad .+= delta .* (formed - prob)
                form_hess .-= (prob * (1 - prob)) .* (delta * delta')
            end

            # Dissolution model: condition on edges at t-1
            for i in 1:n, j in (model.directed ? (1:n) : (i+1:n))
                i == j && continue
                !has_edge(prev_net, i, j) && continue  # Only edges can dissolve

                # Compute change statistics for persistence
                delta = [change_stat(term, curr_net, i, j) for term in model.formula.dissolution_terms]
                persisted = has_edge(curr_net, i, j) ? 1.0 : 0.0

                eta = dot(diss_coef, delta)
                prob = 1.0 / (1.0 + exp(-eta))

                diss_grad .+= delta .* (persisted - prob)
                diss_hess .-= (prob * (1 - prob)) .* (delta * delta')
            end
        end

        # Newton-Raphson update
        if det(form_hess) != 0 && !any(isnan, form_hess)
            form_step = -form_hess \ form_grad
            form_coef .+= form_step
        end

        if det(diss_hess) != 0 && !any(isnan, diss_hess)
            diss_step = -diss_hess \ diss_grad
            diss_coef .+= diss_step
        end

        # Check convergence
        if maximum(abs.(form_grad)) < tol && maximum(abs.(diss_grad)) < tol
            # Compute standard errors from Hessian
            form_se = sqrt.(abs.(diag(pinv(-form_hess))))
            diss_se = sqrt.(abs.(diag(pinv(-diss_hess))))

            return STERGMResult{T}(
                model, form_coef, diss_coef, form_se, diss_se,
                NaN, :cmle, true
            )
        end
    end

    # Did not converge
    form_se = fill(NaN, n_form)
    diss_se = fill(NaN, n_diss)

    return STERGMResult{T}(
        model, form_coef, diss_coef, form_se, diss_se,
        NaN, :cmle, false
    )
end

"""
    cmple(model::STERGMModel; kwargs...) -> STERGMResult

Conditional Maximum Pseudo-Likelihood Estimation for STERGM.
"""
function cmple(model::STERGMModel{T}; kwargs...) where T
    # CMPLE is essentially the same as CMLE for STERGM since the
    # model is already separable
    return cmle(model; kwargs...)
end

"""
    egmme(model::STERGMModel; kwargs...) -> STERGMResult

Equilibrium Generalized Method of Moments Estimation for STERGM.
Used when the network sequence is a single cross-sectional observation
assumed to be from equilibrium.
"""
function egmme(model::STERGMModel{T};
               target_stats::Union{Nothing, Vector{Float64}}=nothing,
               maxiter::Int=100) where T

    @warn "EGMME estimation is experimental"

    n_form = length(model.formula.formation_terms)
    n_diss = length(model.formula.dissolution_terms)

    # For EGMME, we match the observed network statistics
    # assuming the process is at equilibrium

    # Compute observed statistics from the final network
    obs_net = model.networks[end]

    form_stats = [compute(term, obs_net) for term in model.formula.formation_terms]
    diss_stats = [compute(term, obs_net) for term in model.formula.dissolution_terms]

    # Simple initialization
    form_coef = zeros(n_form)
    diss_coef = zeros(n_diss)

    # For now, return placeholder
    form_se = fill(NaN, n_form)
    diss_se = fill(NaN, n_diss)

    return STERGMResult{T}(
        model, form_coef, diss_coef, form_se, diss_se,
        NaN, :egmme, false
    )
end

# =============================================================================
# Simulation Functions
# =============================================================================

"""
    simulate_stergm(result::STERGMResult, n_steps::Int; kwargs...) -> Vector{Network}

Simulate a network sequence from a fitted STERGM.
"""
function simulate_stergm(result::STERGMResult{T}, n_steps::Int;
                         start_net::Union{Nothing, Network{T}}=nothing,
                         burnin::Int=0) where T

    # Start from provided network or last observed
    current = isnothing(start_net) ? deepcopy(result.model.networks[end]) : deepcopy(start_net)
    networks = Network{T}[]

    for step in 1:(burnin + n_steps)
        current = simulate_one_step(current, result)
        if step > burnin
            push!(networks, deepcopy(current))
        end
    end

    return networks
end

"""
    simulate_one_step(net::Network, result::STERGMResult) -> Network

Simulate one time step of STERGM dynamics.
"""
function simulate_one_step(net::Network{T}, result::STERGMResult{T}) where T
    new_net = deepcopy(net)
    n = nv(net)
    directed = is_directed(net)

    # Formation step: for each non-edge, consider forming
    for i in 1:n, j in (directed ? (1:n) : (i+1:n))
        i == j && continue
        has_edge(net, i, j) && continue

        # Compute formation probability
        delta = [change_stat(term, new_net, i, j) for term in result.model.formula.formation_terms]
        eta = dot(result.formation_coef, delta)
        prob = 1.0 / (1.0 + exp(-eta))

        if rand() < prob
            add_edge!(new_net, i, j)
        end
    end

    # Dissolution step: for each edge, consider dissolving
    for i in 1:n, j in (directed ? (1:n) : (i+1:n))
        i == j && continue
        !has_edge(net, i, j) && continue  # Only consider edges that existed

        # Compute persistence probability
        delta = [change_stat(term, new_net, i, j) for term in result.model.formula.dissolution_terms]
        eta = dot(result.dissolution_coef, delta)
        persist_prob = 1.0 / (1.0 + exp(-eta))

        if rand() > persist_prob  # Dissolve with prob 1 - persist_prob
            rem_edge!(new_net, i, j)
        end
    end

    return new_net
end

"""
    simulate_network_sequence(formula::STERGM, init_net::Network, n_steps::Int;
                              form_coef, diss_coef, kwargs...) -> Vector{Network}

Simulate a network sequence from STERGM parameters.
"""
function simulate_network_sequence(formula::STERGM, init_net::Network{T}, n_steps::Int;
                                   form_coef::Vector{Float64},
                                   diss_coef::Vector{Float64},
                                   burnin::Int=100) where T
    current = deepcopy(init_net)
    networks = Network{T}[]
    n = nv(init_net)
    directed = is_directed(init_net)

    for step in 1:(burnin + n_steps)
        new_net = deepcopy(current)

        # Formation
        for i in 1:n, j in (directed ? (1:n) : (i+1:n))
            i == j && continue
            has_edge(current, i, j) && continue

            delta = [change_stat(term, new_net, i, j) for term in formula.formation_terms]
            eta = dot(form_coef, delta)
            prob = 1.0 / (1.0 + exp(-eta))

            if rand() < prob
                add_edge!(new_net, i, j)
            end
        end

        # Dissolution
        for i in 1:n, j in (directed ? (1:n) : (i+1:n))
            i == j && continue
            !has_edge(current, i, j) && continue

            delta = [change_stat(term, new_net, i, j) for term in formula.dissolution_terms]
            eta = dot(diss_coef, delta)
            persist_prob = 1.0 / (1.0 + exp(-eta))

            if rand() > persist_prob
                rem_edge!(new_net, i, j)
            end
        end

        current = new_net
        if step > burnin
            push!(networks, deepcopy(current))
        end
    end

    return networks
end

# =============================================================================
# Diagnostics
# =============================================================================

"""
    stergm_gof(result::STERGMResult; kwargs...) -> NamedTuple

Goodness-of-fit diagnostics for STERGM.
"""
function stergm_gof(result::STERGMResult{T};
                    n_sim::Int=100,
                    statistics::Vector{Function}=[ne, x -> mean(degree(x))]) where T

    # Simulate networks
    sim_nets = simulate_stergm(result, n_sim)

    # Compute statistics on observed final network
    obs_net = result.model.networks[end]
    obs_stats = [f(obs_net) for f in statistics]

    # Compute statistics on simulated networks
    sim_stats = [mean([f(net) for net in sim_nets]) for f in statistics]
    sim_sds = [std([f(net) for net in sim_nets]) for f in statistics]

    return (
        observed=obs_stats,
        simulated_mean=sim_stats,
        simulated_sd=sim_sds,
        z_scores=[(obs - sim) / sd for (obs, sim, sd) in zip(obs_stats, sim_stats, sim_sds)]
    )
end

end # module
