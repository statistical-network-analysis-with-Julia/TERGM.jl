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

        g = stergm_gof(result; n_sim=30, rng=rng)
        @test g.n_sim == 30
        @test g.observed.formed == 2.0
        @test g.observed.persisted == 2.0
        @test 0.0 <= g.p_values.formed <= 1.0
        # The saturated edges-only model should fit its own data
        @test g.p_values.formed > 0.01
        @test g.p_values.persisted > 0.01
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
