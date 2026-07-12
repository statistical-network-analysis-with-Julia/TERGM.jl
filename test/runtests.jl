using TERGM
using ERGM
using Network
using Graphs: src, dst
using Random
using Statistics
using Test

# Two-panel directed fixture on 5 actors
# t0: 1→2, 2→1, 3→4, 4→5
# t1: 1→2 (persists), 2→1 dissolves, 3→4 (persists), 4→5 dissolves,
#     2→3 forms, 5→1 forms
function fixture_panels()
    t0 = network(5)
    for (i, j) in [(1, 2), (2, 1), (3, 4), (4, 5)]
        add_edge!(t0, i, j)
    end
    t1 = network(5)
    for (i, j) in [(1, 2), (3, 4), (2, 3), (5, 1)]
        add_edge!(t1, i, j)
    end
    return [t0, t1]
end

# T-panel directed sequence on n actors with planted persistence
function random_panels(rng; n::Int=10, T::Int=5, p0::Float64=0.15,
                       p_diss::Float64=0.3, p_form::Float64=0.08)
    panels = Network.Network{Int}[]
    net = network(n)
    for i in 1:n, j in 1:n
        i != j && rand(rng) < p0 && add_edge!(net, i, j)
    end
    push!(panels, net)
    for _ in 2:T
        nxt = TERGM._copy_net(panels[end])
        for i in 1:n, j in 1:n
            i == j && continue
            if has_edge(panels[end], i, j)
                rand(rng) < p_diss && rem_edge!(nxt, i, j)
            elseif rand(rng) < p_form
                add_edge!(nxt, i, j)
            end
        end
        push!(panels, nxt)
    end
    return panels
end

@testset "TERGM.jl" begin
    @testset "Auxiliary networks (Krivitsky-Handcock)" begin
        nets = fixture_panels()
        yplus = formation_network(nets[1], nets[2])
        yminus = dissolution_network(nets[1], nets[2])

        # Y⁺ = union: 6 edges
        @test ne(yplus) == 6
        @test has_edge(yplus, 2, 1)  # from t0
        @test has_edge(yplus, 2, 3)  # from t1

        # Y⁻ = intersection: 1→2 and 3→4
        @test ne(yminus) == 2
        @test has_edge(yminus, 1, 2)
        @test has_edge(yminus, 3, 4)
        @test !has_edge(yminus, 2, 1)
    end

    @testset "Temporal term values" begin
        nets = fixture_panels()
        prev, curr = nets[1], nets[2]

        # Agreement: 20 dyads − (2 dissolved + 2 formed) = 16
        @test compute(EdgeStability(), curr, prev) == 16.0
        @test compute(PersistentEdge(), curr, prev) == 2.0
        @test compute(NewEdge(), curr, prev) == 2.0
        # Delayed reciprocity: 2→3? prev 3→2 no. 5→1? prev 1→5 no.
        # 1→2? prev 2→1 yes ✓. 3→4? prev 4→3 no. → 1
        @test compute(Delrecip(), curr, prev) == 1.0

        # Temporal terms demand prev_net
        @test_throws ErrorException compute(EdgeStability(), curr)
        @test_throws ErrorException change_stat(EdgeStability(), curr, 1, 2)
    end

    @testset "Temporal change stats are state-independent add-direction" begin
        nets = fixture_panels()
        prev, curr = nets[1], nets[2]

        for term in [EdgeStability(), Delrecip(), PersistentEdge(), NewEdge()]
            for i in 1:5, j in 1:5
                i == j && continue
                # Brute force on a copy
                work = TERGM._copy_net(curr)
                had = has_edge(work, i, j)
                had && rem_edge!(work, i, j)
                s0 = compute(term, work, prev)
                add_edge!(work, i, j)
                s1 = compute(term, work, prev)
                expected = s1 - s0

                @test change_stat(term, curr, i, j, prev) ≈ expected atol = 1e-12
            end
        end
    end

    @testset "Edge ages" begin
        nets = fixture_panels()
        ages = edge_ages(nets)
        @test ages[(1, 2)] == 2   # present in both panels
        @test ages[(2, 3)] == 1   # formed at t1
        @test mean_edge_age(nets) == 1.5  # ages 2,2,1,1

        @test isnan(mean_edge_age([network(3)]))
        @test_throws ArgumentError edge_ages(Network.Network{Int}[])
    end

    @testset "CMPLE analytic (edges-only)" begin
        nets = fixture_panels()
        result = stergm(nets, [Edges()], [Edges()])

        @test result.converged
        @test result.method == :cmple

        # Formation: 16 prev non-edges, 2 formed → logit(2/16)
        @test result.formation_coef[1] ≈ log((2 / 16) / (1 - 2 / 16)) atol = 1e-5
        # Persistence: 4 prev edges, 2 persisted → logit(1/2) = 0
        @test result.dissolution_coef[1] ≈ 0.0 atol = 1e-5

        @test isfinite(result.loglik_formation)
        @test isfinite(result.loglik_dissolution)
        @test all(result.formation_se .> 0)
    end

    @testset "Auxiliary networks preserve attributes; nodal terms estimable" begin
        # Regression: _copy_net/dissolution_network used to drop vertex
        # attributes, so nodal terms (NodeMatch, ...) in formation and
        # dissolution formulas produced all-zero design-matrix columns
        rng = Random.Xoshiro(7)
        n = 10
        t0 = network(n)
        set_vertex_attribute!(t0, :group,
                              Dict(v => (v <= 5 ? "a" : "b") for v in 1:n))
        for i in 1:n, j in 1:n
            i != j && rand(rng) < 0.1 && add_edge!(t0, i, j)
        end
        # Planted homophily: same-group ties form readily, cross-group rarely
        t1 = copy(t0)
        for i in 1:n, j in 1:n
            (i == j || has_edge(t0, i, j)) && continue
            same = (i <= 5) == (j <= 5)
            rand(rng) < (same ? 0.7 : 0.05) && add_edge!(t1, i, j)
        end

        # Y⁺ and Y⁻ carry the vertex attributes
        yplus = formation_network(t0, t1)
        yminus = dissolution_network(t0, t1)
        @test get_vertex_attribute(yplus, :group, 1) == "a"
        @test get_vertex_attribute(yminus, :group, 6) == "b"

        # The nodal term's design-matrix column — NodeMatch change stats on
        # Y⁺ over the formation-free dyads — must not be identically zero
        col = [change_stat(NodeMatch(:group), yplus, i, j)
               for i in 1:n for j in 1:n if i != j && !has_edge(t0, i, j)]
        @test any(!iszero, col)

        # CMPLE recovers the planted homophily sign
        result = stergm([t0, t1], [Edges(), NodeMatch(:group)], [Edges()])
        @test result.converged
        @test result.formation_coef[2] > 0

        # Simulated transitions carry attributes forward too
        yt = simulate_stergm(t1, STERGM([Edges()], [Edges()]), [-1.0], [0.5];
                             burnin=500, rng=rng)
        @test get_vertex_attribute(yt, :group, 3) == "a"
    end

    @testset "Temporal and dyad-dependent terms are estimable" begin
        rng = Random.Xoshiro(31)
        # Build a 3-panel sequence with some persistence
        panels = Network.Network{Int}[]
        net = network(10)
        for i in 1:10, j in 1:10
            i != j && rand(rng) < 0.15 && add_edge!(net, i, j)
        end
        push!(panels, net)
        for t in 2:3
            nxt = TERGM._copy_net(panels[end])
            for i in 1:10, j in 1:10
                i == j && continue
                if has_edge(panels[end], i, j)
                    rand(rng) < 0.3 && rem_edge!(nxt, i, j)
                elseif rand(rng) < 0.08
                    add_edge!(nxt, i, j)
                end
            end
            push!(panels, nxt)
        end

        # Formerly a MethodError: bare temporal term + dyad-dependent term
        result = stergm(panels, [Edges(), Delrecip()], [Edges(), Mutual()])
        @test result isa STERGMResult
        @test all(isfinite, result.formation_coef)
        @test all(isfinite, result.dissolution_coef)

        # cmle warns and falls back to CMPLE
        r2 = @test_logs (:warn, r"CMPLE") match_mode = :any begin
            stergm(panels, [Edges()], [Edges()]; method=:cmle)
        end
        @test r2.method == :cmle

        # EGMME raises instead of returning zeros
        @test_throws ErrorException stergm(panels, [Edges()], [Edges()];
                                           method=:egmme)
    end

    @testset "Simulation respects the separable structure" begin
        rng = Random.Xoshiro(17)
        prev = network(12)
        for i in 1:12, j in 1:12
            i != j && rand(rng) < 0.2 && add_edge!(prev, i, j)
        end
        formula = STERGM([Edges()], [Edges()])

        # Extreme coefficients: no formation, full persistence
        yt = simulate_stergm(prev, formula, [-30.0], [30.0]; rng=rng)
        @test ne(yt) == ne(prev)
        @test all(has_edge(yt, src(e), dst(e)) for e in edges(prev))

        # No formation, no persistence → empty
        yt2 = simulate_stergm(prev, formula, [-30.0], [-30.0]; rng=rng)
        @test ne(yt2) == 0

        # Moderate model: formation rate among free dyads ≈ σ(θ⁺)
        θf, θd = -2.0, 1.0
        n_free = 12 * 11 - ne(prev)
        formed = Float64[]
        persisted = Float64[]
        for _ in 1:30
            y = simulate_stergm(prev, formula, [θf], [θd]; burnin=2000, rng=rng)
            push!(formed, compute(NewEdge(), y, prev))
            push!(persisted, compute(PersistentEdge(), y, prev))
        end
        @test mean(formed) ≈ n_free / (1 + exp(-θf)) rtol = 0.25
        @test mean(persisted) ≈ ne(prev) / (1 + exp(-θd)) rtol = 0.25
    end

    @testset "Simulation-estimation round trip" begin
        rng = Random.Xoshiro(23)
        init = network(14)
        for i in 1:14, j in 1:14
            i != j && rand(rng) < 0.15 && add_edge!(init, i, j)
        end
        formula = STERGM([Edges()], [Edges()])
        θf_true, θd_true = -1.8, 0.8

        seq = simulate_network_sequence(formula, init, 6, [θf_true], [θd_true];
                                        burnin=4000, rng=rng)
        result = stergm(seq, [Edges()], [Edges()])

        @test result.converged
        @test result.formation_coef[1] ≈ θf_true atol = 0.35
        @test result.dissolution_coef[1] ≈ θd_true atol = 0.45
    end

    @testset "Goodness of fit" begin
        rng = Random.Xoshiro(41)
        nets = fixture_panels()
        result = stergm(nets, [Edges()], [Edges()])

        g = gof(result; n_sim=30, rng=rng)

        # gof extends Network.jl's shared generic and returns the shared
        # GOFResult container; stergm_gof remains as an alias
        @test stergm_gof === gof
        @test g isa GOFResult
        @test n_simulations(g) == 30
        stat = only(g.statistics)
        @test stat.labels == ["formed", "persisted"]
        @test stat.observed == [2.0, 2.0]
        @test all(0.0 .< stat.p_values .<= 1.0)
        # The saturated edges-only model should fit its own data
        @test stat.p_values[1] > 0.01
        @test stat.p_values[2] > 0.01

        # ... and renders through the shared formatted display
        out = sprint(show, g)
        @test occursin("Goodness-of-fit assessment: STERGM", out)
        @test occursin("MC p-value", out)
    end

    @testset "StatsAPI accessors" begin
        nets = fixture_panels()
        r = stergm(nets, [Edges()], [Edges()])

        @test coef(r) == vcat(r.formation_coef, r.dissolution_coef)
        @test stderror(r) == vcat(r.formation_se, r.dissolution_se)
        V = vcov(r)
        @test size(V) == (2, 2)
        @test V[1, 2] == 0.0  # separable → block-diagonal
        @test sqrt(V[1, 1]) ≈ r.formation_se[1]
        @test sqrt(V[2, 2]) ≈ r.dissolution_se[1]
        @test loglikelihood(r) ≈ r.loglik_formation + r.loglik_dissolution
        @test dof(r) == 2
        @test nobs(r) == 20  # 5·4 ordered dyads × 1 transition
        @test aic(r) ≈ -2 * loglikelihood(r) + 4
        @test bic(r) ≈ -2 * loglikelihood(r) + 2 * log(20)
    end

    @testset "Formula validation at model construction" begin
        nets = fixture_panels()

        # Attribute-based term whose attribute is missing → ArgumentError
        # listing the available attributes
        err = try
            STERGMModel(STERGM([Edges(), NodeMatch(:group)], [Edges()]), nets)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin(":group", err.msg)
        @test occursin("Available vertex attributes", err.msg)
        @test occursin("formation", err.msg)

        # ... and the same through the stergm() entry point, on the
        # dissolution side
        err2 = try
            stergm(nets, [Edges()], [Edges(), NodeCov(:wealth)])
            nothing
        catch e
            e
        end
        @test err2 isa ArgumentError
        @test occursin("dissolution", err2.msg)
        @test occursin(":wealth", err2.msg)

        # Attribute present on every panel → constructs (and lists it when
        # another attribute is missing)
        for net in nets
            set_vertex_attribute!(net, :grp,
                                  Dict(v => (v <= 2 ? "a" : "b") for v in 1:5))
        end
        @test STERGMModel(STERGM([Edges(), NodeMatch(:grp)], [Edges()]),
                          nets) isa STERGMModel
        err3 = try
            STERGMModel(STERGM([Edges(), NodeMatch(:other)], [Edges()]), nets)
            nothing
        catch e
            e
        end
        @test err3 isa ArgumentError
        @test occursin(":grp", err3.msg)

        # Direction-incompatible terms on undirected panels are rejected
        un = [network(4; directed=false), network(4; directed=false)]
        @test_throws ArgumentError STERGMModel(
            STERGM([Edges()], [Edges(), Mutual()]), un)
        @test_throws ArgumentError STERGMModel(
            STERGM([Edges(), Delrecip()], [Edges()]), un)
        # ... but accepted on directed panels
        @test STERGMModel(STERGM([Edges(), Delrecip()], [Edges(), Mutual()]),
                          fixture_panels()) isa STERGMModel
    end

    @testset "Per-transition block-bootstrap SEs (btergm-style)" begin
        panels = random_panels(Random.Xoshiro(99))

        r_h = stergm(panels, [Edges()], [Edges()])
        r_b = stergm(panels, [Edges()], [Edges()];
                     se=:bootstrap, n_boot=60, rng=Random.Xoshiro(1))

        # Bootstrap changes only the uncertainty, not the point estimates
        @test r_b.formation_coef == r_h.formation_coef
        @test r_b.dissolution_coef == r_h.dissolution_coef
        @test r_h.se_type == :hessian
        @test r_b.se_type == :bootstrap

        @test all(isfinite, stderror(r_b))
        @test all(stderror(r_b) .> 0)
        V = vcov(r_b)
        @test size(V) == (2, 2)
        @test V ≈ V'
        @test sqrt(V[1, 1]) ≈ r_b.formation_se[1]
        @test sqrt(V[2, 2]) ≈ r_b.dissolution_se[1]

        # Reproducible under the same rng seed
        r_b2 = stergm(panels, [Edges()], [Edges()];
                      se=:bootstrap, n_boot=60, rng=Random.Xoshiro(1))
        @test stderror(r_b2) == stderror(r_b)

        # Resampling transitions needs at least two of them
        @test_throws ArgumentError stergm(fixture_panels(), [Edges()], [Edges()];
                                          se=:bootstrap)
        # Unknown se choice fails loudly
        @test_throws ArgumentError stergm(panels, [Edges()], [Edges()];
                                          se=:jackknife)
        @test_throws ArgumentError stergm(panels, [Edges()], [Edges()];
                                          se=:bootstrap, n_boot=1)
    end

    @testset "Dyad-dependence caveat in show()" begin
        # Temporal terms condition only on the exogenous previous network
        @test !is_dyad_dependent(EdgeStability())
        @test !is_dyad_dependent(Delrecip())
        @test !is_dyad_dependent(PersistentEdge())
        @test !is_dyad_dependent(NewEdge())
        @test is_dyad_dependent(Mutual())

        panels = random_panels(Random.Xoshiro(31))

        # Dyad-dependent formula + naive SEs → warning with pointer to
        # se=:bootstrap
        r_dep = stergm(panels, [Edges()], [Edges(), Mutual()])
        out = sprint(show, r_dep)
        @test occursin("dyad-dependent", out)
        @test occursin("anticonservative", out)
        @test occursin("se=:bootstrap", out)

        # Both blocks render through the shared Network.jl coefficient
        # printer (R-style columns, one significance-code legend)
        @test occursin("Formation:", out)
        @test occursin("Dissolution (persistence):", out)
        @test count("Pr(>|z|)", out) == 2
        @test count("Signif. codes:", out) == 1

        # Dyad-independent formula (incl. temporal terms) → no caveat
        r_ind = stergm(panels, [Edges(), Delrecip()], [Edges()])
        @test !occursin("dyad-dependent", sprint(show, r_ind))

        # Bootstrap SEs → milder note, no anticonservative warning
        r_boot = stergm(panels, [Edges()], [Edges(), Mutual()];
                        se=:bootstrap, n_boot=30, rng=Random.Xoshiro(7))
        ob = sprint(show, r_boot)
        @test occursin("block-bootstrap", ob)
        @test !occursin("anticonservative", ob)
    end

    @testset "Input validation" begin
        @test_throws ArgumentError STERGM(AbstractERGMTerm[], [Edges()])
        @test_throws ArgumentError STERGMModel(STERGM([Edges()], [Edges()]),
                                               [network(3)])
        nets = [network(3), network(4)]
        @test_throws ArgumentError STERGMModel(STERGM([Edges()], [Edges()]), nets)
        @test fit_stergm === stergm
    end
end
