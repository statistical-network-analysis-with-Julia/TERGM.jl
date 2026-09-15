using TERGM
using ERGM
using Networks
using Graphs: src, dst, inneighbors
using Random
using Statistics
using LinearAlgebra: dot, diag
using TOML
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
    panels = Network{Int}[]
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
        @test_throws ArgumentError edge_ages(Network{Int}[])
    end

    @testset "CMPLE analytic (edges-only)" begin
        nets = fixture_panels()
        result = stergm(nets, [Edges()], [Edges()])

        @test result.converged
        @test result.method == :cmple

        # Formation: 16 prev non-edges, 2 formed → logit(2/16)
        @test result.formation_coef[1] ≈ log((2 / 16) / (1 - 2 / 16)) atol = 1e-5
        # Persistence: 4 prev edges, 2 persisted → logit(1/2) = 0
        @test result.persistence_coef[1] ≈ 0.0 atol = 1e-5

        @test isfinite(result.loglik_formation)
        @test isfinite(result.loglik_dissolution)
        @test all(result.formation_se .> 0)
    end

    # ------------------------------------------------------------------
    # Allocation regression on the CMPLE derivative loop (review finding 15)
    #
    # `_logistic_fit` used to carry its own copy of the logistic loop with a
    # per-row `(pr*(1-pr)) .* (x * x')` inside it: a fresh p×p matrix on every
    # one of the design rows of every Newton evaluation (470 KB per evaluation on
    # a 25-actor, 8-wave panel). It now runs on the shared, workspace-backed
    # `Networks.logistic_derivatives` (hosted in Networks.jl since the 2026-09
    # panel; ERGM.jl re-exports the same binding) — the same builder ERGM's,
    # ERGMMulti's and ERGMRank's MPLEs use, in the binomial-row form
    # `(X, n_tot, n_one)` that `_logistic_fit` hands to ERGM's design fitter.
    # This test is what stops the outer product coming back, and it measures the
    # closure the fitter builds from the package's OWN CMPLE design rows.
    # ------------------------------------------------------------------
    @testset "CMPLE derivative evaluations allocate O(p²), not O(rows · p²)" begin
        function evaluation_allocs(n_actors, n_waves)
            rng = Random.Xoshiro(6)
            nets = Network{Int}[]
            for _ in 1:n_waves
                net = network(n_actors; directed=true)
                for i in 1:n_actors, j in 1:n_actors
                    i != j && rand(rng) < 0.15 && add_edge!(net, i, j)
                end
                push!(nets, net)
            end
            model = STERGMModel(STERGM([Edges(), Mutual()], [Edges()]), nets)
            Xf_blocks, yf_blocks, _, _ = TERGM._cmple_blocks(model)
            X = reduce(vcat, Xf_blocks)
            y = reduce(vcat, yf_blocks)
            @test Networks.logistic_derivatives === ERGM.logistic_derivatives
            d = Networks.logistic_derivatives(X, ones(length(y)), Float64.(y))
            β = [0.1, 0.05]
            d(β)                    # warm up: @allocated on a first call
            return size(X, 1), @allocated d(β)   # would measure compilation
        end

        rows_small, a_small = evaluation_allocs(6, 2)
        rows_big, a_big = evaluation_allocs(25, 8)
        @test rows_big > 50 * rows_small
        # 50x the design rows, the same allocations: only the (p) gradient and
        # (p×p) Hessian handed back to `newton_fit` are allocated per evaluation.
        @test a_small <= 512
        @test a_big <= 512
        @test a_big <= a_small + 64
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
        panels = Network{Int}[]
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

        # Formerly a MethodError: bare temporal term + dyad-dependent term.
        # On this panel every prior mutual tie persists, so `Persist~mutual`
        # sits at its largest attainable value: no finite CMPLE exists for
        # it, and the fit says so (R ergm's drop: +Inf, SE 0) instead of
        # "converging" on the asymptote as it did before 0.2
        result = @test_logs (:warn, r"Persist~mutual are at their largest attainable values") match_mode=:any stergm(
            panels, [Edges(), Delrecip()], [Edges(), Mutual()])
        @test result isa STERGMResult
        @test result.converged
        @test all(isfinite, result.formation_coef)
        @test isfinite(result.persistence_coef[1])
        @test result.persistence_coef[2] == Inf
        @test result.persistence_se[2] == 0.0
        @test dof(result) == 3

        # cmle raises rather than silently falling back to CMPLE (it used to
        # warn, return a CMPLE fit, and stamp it `:cmle` — so the result
        # misreported its own estimator). See the "cmle does not masquerade"
        # testset for the full contract.
        @test_throws ArgumentError stergm(panels, [Edges()], [Edges()];
                                          method=:cmle)

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
        @test result.persistence_coef[1] ≈ θd_true atol = 0.45
    end

    @testset "Goodness of fit" begin
        rng = Random.Xoshiro(41)
        nets = fixture_panels()
        result = stergm(nets, [Edges()], [Edges()])

        g = gof(result; n_sim=30, rng=rng)

        # gof extends Networks.jl's shared generic and returns the shared
        # GOFResult container. `gof(fit)` is the ONE verb of the model family
        # (panel item 16): `stergm_gof` is a deprecated forwarding method,
        # not an alias — it warns, forwards every keyword and is gone next
        # release.
        @test TERGM.gof === Networks.gof === ERGM.gof
        @test stergm_gof !== gof
        g_old = @test_deprecated stergm_gof(result; n_sim=3, rng=Random.Xoshiro(5))
        @test g_old isa GOFResult
        @test n_simulations(g_old) == 3
        # ... and forwards to exactly gof: same seed, same draws
        g_new = gof(result; n_sim=3, rng=Random.Xoshiro(5))
        @test [s.simulated for s in g_old.statistics] == [s.simulated for s in g_new.statistics]
        @test g isa GOFResult
        @test n_simulations(g) == 30
        stat = g.statistics[1]
        @test stat.name == "tie changes"
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

    # ------------------------------------------------------------------
    # gof panels (panel 2026-09 round 2): the pooled formed/persisted counts
    # are the sufficient statistics of Form~edges / Persist~edges, which a
    # CMPLE reproduces by construction, so `gof` also compares the formation
    # and persistence MODEL statistics on Y⁺/Y⁻ and the degree distribution
    # of Y_t — after gof.tergm.
    # ------------------------------------------------------------------
    @testset "gof compares model statistics and degree distributions" begin
        panels = random_panels(Random.Xoshiro(20); n=10, T=4)
        fit = stergm(panels, [Edges(), Delrecip()], [Edges(), Mutual()])
        g = gof(fit; n_sim=40, burnin=300, rng=Random.Xoshiro(3))
        @test [s.name for s in g.statistics] ==
              ["tie changes", "formation statistics", "persistence statistics",
               "idegree", "odegree"]
        form, pers, ideg, odeg = g.statistics[2:5]
        # labels are the fitted coefficients' (R's labels, resolved on the panel)
        @test form.labels == ["edges", "delrecip"]
        @test pers.labels == ["edges", "mutual"]
        @test ideg.labels == odeg.labels == string.(0:9)
        # observed model statistics: every term on the observed Y⁺ / Y⁻,
        # pooled over transitions
        exp_f = zeros(2); exp_d = zeros(2)
        for t in 2:length(panels)
            prev, curr = panels[t-1], panels[t]
            yplus, yminus = formation_network(prev, curr), dissolution_network(prev, curr)
            exp_f .+= [compute(Edges(), yplus), compute(Delrecip(), yplus, prev)]
            exp_d .+= [compute(Edges(), yminus), compute(Mutual(), yminus)]
        end
        @test form.observed == exp_f
        @test pers.observed == exp_d
        # the `edges` levels are the tie-change counts plus the fixed prior ties
        # (Y⁺ ⊇ Y_{t−1}) / exactly the persisted ties (Y⁻ ⊆ Y_{t−1}) — the same
        # numbers as the "tie changes" panel, so the CMPLE reproduces them
        n_prior = sum(ne(panels[t]) for t in 1:length(panels)-1)
        @test form.observed[1] == g.statistics[1].observed[1] + n_prior
        @test pers.observed[1] == g.statistics[1].observed[2]
        @test form.simulated[:, 1] == g.statistics[1].simulated[:, 1] .+ n_prior
        @test pers.simulated[:, 1] == g.statistics[1].simulated[:, 2]
        # the degree histograms count every vertex of every transition's Y_t
        @test sum(ideg.observed) == sum(odeg.observed) == 10 * (length(panels) - 1)
        @test all(sum(ideg.simulated; dims=2) .== 10 * (length(panels) - 1))
        @test size(ideg.simulated) == (40, 10)
        obs_in = zeros(10)
        for t in 2:length(panels), v in 1:10
            obs_in[length(inneighbors(panels[t], v)) + 1] += 1
        end
        @test ideg.observed == obs_in
        @test all(0 .< s.p_values[k] <= 1 for s in g.statistics for k in eachindex(s.p_values))

        # An undirected panel has one "degree" panel
        un = Network{Int,false}[]
        rng = Random.Xoshiro(8)
        u = network(8; directed=false)
        for i in 1:8, j in i+1:8
            rand(rng) < 0.3 && add_edge!(u, i, j)
        end
        push!(un, u)
        for _ in 1:2
            w = copy(un[end])
            for i in 1:8, j in i+1:8
                if has_edge(un[end], i, j)
                    rand(rng) < 0.3 && rem_edge!(w, i, j)
                elseif rand(rng) < 0.1
                    add_edge!(w, i, j)
                end
            end
            push!(un, w)
        end
        gu = gof(stergm(un, [Edges()], [Edges()]); n_sim=5, burnin=50, rng=Random.Xoshiro(1))
        @test [s.name for s in gu.statistics] ==
              ["tie changes", "formation statistics", "persistence statistics", "degree"]
        @test gu.statistics[4].labels == string.(0:7)
        @test occursin("Goodness-of-fit for degree", sprint(show, gu))
    end

    @testset "StatsAPI accessors" begin
        nets = fixture_panels()
        r = stergm(nets, [Edges()], [Edges()])

        @test coef(r) == vcat(r.formation_coef, r.persistence_coef)
        @test stderror(r) == vcat(r.formation_se, r.persistence_se)
        V = vcov(r)
        @test size(V) == (2, 2)
        @test V[1, 2] == 0.0  # separable → block-diagonal
        @test sqrt(V[1, 1]) ≈ r.formation_se[1]
        @test sqrt(V[2, 2]) ≈ r.persistence_se[1]
        @test loglikelihood(r) ≈ r.loglik_formation + r.loglik_dissolution
        @test dof(r) == 2
        @test nobs(r) == 20  # 5·4 ordered dyads × 1 transition
        @test aic(r) ≈ -2 * loglikelihood(r) + 4
        @test bic(r) ≈ -2 * loglikelihood(r) + 2 * log(20)

        # The FULL StatsAPI surface (panel 2026-09, item 15), pinned by
        # Networks' checker for a Hessian fit and a bootstrap fit alike; the
        # verbs are the ONE StatsAPI binding, re-exported (never a local
        # same-named function)
        @test check_statsapi(r; strict=true) !== nothing
        @test all(check_statsapi(r))
        panels = random_panels(Random.Xoshiro(3); n=8, T=4)
        rb = stergm(panels, [Edges(), Delrecip()], [Edges()];
                    se=:bootstrap, n_boot=20, rng=Random.Xoshiro(1))
        @test check_statsapi(rb; strict=true) !== nothing
        @test TERGM.coeftable === TERGM.StatsAPI.coeftable === Networks.coeftable === ERGM.coeftable
        @test TERGM.confint === TERGM.StatsAPI.confint === ERGM.confint
        @test TERGM.coef === TERGM.StatsAPI.coef

        # confint: Wald limits from the SEs the fit reports, in coef(r)'s order
        ci = confint(r)
        @test size(ci) == (2, 2)
        q975 = 1.959963984540054
        @test maximum(abs.((ci[:, 2] .- ci[:, 1]) .- 2 * q975 .* stderror(r))) < 1e-12
        @test ci[:, 1] ≈ coef(r) .- q975 .* stderror(r)
        @test all(ci[:, 1] .< coef(r) .< ci[:, 2])
        ci90 = confint(r; level=0.9)
        @test all(ci90[:, 2] .- ci90[:, 1] .< ci[:, 2] .- ci[:, 1])
        @test_throws ArgumentError confint(r; level=1.0)
        @test_throws ArgumentError confint(r; level=0.0)
        # a bootstrap fit's intervals use the bootstrap SEs
        cib = confint(rb)
        @test cib[:, 2] .- cib[:, 1] ≈ 2 * q975 .* stderror(rb)

        # coeftable: tergm's `Form(1)~`/`Persist(1)~` labels on ERGM's R term
        # labels, one row per stacked coefficient, the same numbers show prints
        tbl = coeftable(r)
        @test tbl isa CoefficientTable
        @test tbl.names == ["Form(1)~edges", "Persist(1)~edges"]
        @test tbl.estimates == coef(r)
        @test tbl.std_errors == stderror(r)
        @test tbl["Persist(1)~edges"].estimate == persistence_coef(r)[1]
        @test tbl[1].std_error == formation_se(r)[1]
        @test tbl.p_values == z_pvalues(tbl.z_values)
        tb = coeftable(rb)
        @test tb.names == ["Form(1)~edges", "Form(1)~delrecip", "Persist(1)~edges"]
        @test tb.std_errors == stderror(rb)
        # what is shown IS what is inspected: every printed estimate/SE of
        # both blocks appears in the table, and the table renders through the
        # same presentation layer (four R-style columns, one legend)
        out = sprint(show, rb)
        for row in tb
            @test occursin(string(round(row.estimate, digits=4)), out)
            @test occursin(string(round(row.std_error, digits=4)), out)
        end
        tout = sprint(show, tb)
        @test occursin("Form(1)~delrecip", tout)
        @test count("Pr(>|z|)", tout) == 1
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

        # ---- Validation is ERGM.jl's `_validate_formula`, run on EVERY panel;
        # only the "side, panel" prefix is TERGM's (panel 2026-09, item 12).
        msg_of(f) = try; f(); "no error"; catch e; e isa ArgumentError ? e.msg : rethrow(); end

        # (a) NA-completeness: an attribute set on vertices 1..4 only of panel 2
        # is refused, naming the panel, the term, the attribute and vertex 5
        nets_na = fixture_panels()
        set_vertex_attribute!(nets_na[1], :grp,
                              Dict(v => (v <= 2 ? "a" : "b") for v in 1:5))
        set_vertex_attribute!(nets_na[2], :grp,
                              Dict(v => (v <= 2 ? "a" : "b") for v in 1:4))
        m_na = msg_of(() -> STERGMModel(STERGM([Edges(), NodeMatch(:grp)], [Edges()]),
                                        nets_na))
        @test occursin("formation model, panel 2", m_na)
        @test occursin("nodematch.grp", m_na)
        @test occursin(":grp", m_na)
        @test occursin("every vertex", m_na)
        @test occursin("vertices 5", m_na)
        # ... and the same attribute complete on both panels is accepted
        set_vertex_attribute!(nets_na[2], :grp, 5, "b")
        @test STERGMModel(STERGM([Edges(), NodeMatch(:grp)], [Edges()]),
                          nets_na) isa STERGMModel

        # (b) Undirected-only terms on directed panels are refused with the
        # directed variant named (before: a directed panel silently got out-stars)
        for (term, hint) in ((Kstar(2), "OStar"), (GWDegree(0.5), "GWODegree"),
                             (Degree(1), "ODegree"))
            m = msg_of(() -> STERGMModel(STERGM([Edges(), term], [Edges()]),
                                         fixture_panels()))
            @test occursin("formation model, panel 1", m)
            @test occursin("undirected networks", m)
            @test occursin(hint, m)
            # ... on either side
            m2 = msg_of(() -> STERGMModel(STERGM([Edges()], [Edges(), term]),
                                          fixture_panels()))
            @test occursin("dissolution model, panel 1", m2)
        end

        # (c) Their directed counterparts, and GWESP at decay 0, are accepted
        @test STERGMModel(STERGM([Edges(), OStar(2), GWODegree(0.5)],
                                 [Edges(), GWESP(0.0)]),
                          fixture_panels()) isa STERGMModel

        # (d) Directed-only terms on undirected panels: message names the side,
        # the panel and the term
        un2 = [network(4; directed=false), network(4; directed=false)]
        m_un = msg_of(() -> STERGMModel(STERGM([Edges()], [Edges(), Delrecip()]), un2))
        @test occursin("dissolution model, panel 1", m_un)
        @test occursin("delrecip", m_un)
        @test occursin("directed networks", m_un)

        # (e) EdgeCov with a wrong-sized covariate matrix is refused (ERGM's
        # covariate-size check comes for free)
        m_cov = msg_of(() -> STERGMModel(STERGM([Edges(), EdgeCov(zeros(4, 4); name="cov")],
                                                [Edges()]), fixture_panels()))
        @test occursin("formation model, panel 1", m_cov)
        @test occursin("4×4", m_cov)
        @test occursin("5 vertices", m_cov)
        @test STERGMModel(STERGM([Edges(), EdgeCov(zeros(5, 5); name="cov")], [Edges()]),
                          fixture_panels()) isa STERGMModel

        # A non-ArgumentError raised inside validation is not re-wrapped
        # (the wrapper only prefixes ERGM's ArgumentErrors)
        @test_throws ArgumentError STERGMModel(STERGM([Edges(), NodeMatch(:none)],
                                                      [Edges()]), fixture_panels())
    end

    @testset "Delrecip refuses undirected networks" begin
        @test requires_directed(Delrecip())
        @test !requires_directed(EdgeStability())
        un = network(4; directed=false)
        add_edge!(un, 1, 2)
        @test_throws ArgumentError compute(Delrecip(), un, un)
        @test_throws ArgumentError change_stat(Delrecip(), un, 1, 2, un)
        msg = try; compute(Delrecip(), un, un); catch e; e.msg; end
        @test occursin("delrecip", msg)
        @test occursin("directed networks", msg)
        # The directed fixture values are unchanged
        nets = fixture_panels()
        @test compute(Delrecip(), nets[2], nets[1]) == 1.0
    end

    @testset "STERGMModel{T,D} is directedness-typed (panel item 19)" begin
        m = STERGMModel(STERGM([Edges()], [Edges()]), fixture_panels())
        @test m isa STERGMModel{Int, true}
        @test isconcretetype(fieldtype(typeof(m), :networks))
        @test is_directed(m) == true
        @test is_directed(typeof(m)) == true

        un = [network(4; directed=false), network(4; directed=false)]
        add_edge!(un[1], 1, 2); add_edge!(un[2], 1, 2); add_edge!(un[2], 2, 3)
        mu = STERGMModel(STERGM([Edges()], [Edges()]), un)
        @test mu isa STERGMModel{Int, false}
        @test is_directed(mu) == false
        @test nobs(cmple(mu)) == 6     # 4·3/2 unordered dyads × 1 transition

        # The result carries the parameters too
        r = cmple(m)
        @test r isa STERGMResult{Int, true}
        @test is_directed(r.model)

        # Mixed directedness is a type mismatch, refused at construction
        @test_throws ArgumentError STERGMModel(STERGM([Edges()], [Edges()]),
                                               [fixture_panels()[1], un[1]])
        # An abstractly typed input vector is collected to the concrete panel type
        mixed_vec = Network[fixture_panels()...]
        @test STERGMModel(STERGM([Edges()], [Edges()]), mixed_vec) isa STERGMModel{Int, true}

        # Simulated sequences are concretely typed as well
        seq = simulate_network_sequence(STERGM([Edges()], [Edges()]),
                                        fixture_panels()[1], 2, [-1.0], [0.5];
                                        burnin=50, rng=Random.Xoshiro(1))
        @test isconcretetype(eltype(seq))
        @test eltype(seq) === Network{Int, true}
        @test length(seq) == 3
    end

    @testset "has_dyad_dependent is ERGM's exported generic" begin
        @test hasmethod(has_dyad_dependent, Tuple{STERGMModel})
        @test !isdefined(TERGM, :_has_dyad_dependent)
        panels = random_panels(Random.Xoshiro(20); n=10, T=4)
        indep = STERGMModel(STERGM([Edges(), Delrecip()], [Edges()]), panels)
        dep = STERGMModel(STERGM([Edges(), Triangle()], [Edges()]), panels)
        dep_diss = STERGMModel(STERGM([Edges()], [Edges(), Triangle()]), panels)
        @test !has_dyad_dependent(indep)
        @test has_dyad_dependent(dep)
        @test has_dyad_dependent(dep_diss)
        # The same binding ERGM exports and the deprecated private alias resolve to
        @test TERGM.has_dyad_dependent === ERGM.has_dyad_dependent
        @test ERGM._has_dyad_dependent === ERGM.has_dyad_dependent
    end

    @testset "Coefficient labels are R's, resolved against the panel" begin
        # (on the 5-actor fixture both GW statistics are at their smallest
        # attainable value — no shared partner on Y⁺/Y⁻ — so the fit drops
        # them, naming them by R's directed labels)
        r = @test_logs (:warn, r"Form~gwesp.OTP.fixed.0.5 are at their smallest") (:warn, r"Persist~gwdsp.OTP.fixed.1 are at their smallest") match_mode=:any stergm(
            fixture_panels(), [Edges(), GWESP(0.5)], [Edges(), GWDSP(1)])
        out = sprint(show, r)
        @test occursin("gwesp.OTP.fixed.0.5", out)   # directed panel: R's label
        @test occursin("gwdsp.OTP.fixed.1", out)     # integer decay prints as 1
        @test TERGM._labels(r.model.formula.formation, r.model.networks[1]) ==
              ["edges", "gwesp.OTP.fixed.0.5"]

        un = [network(4; directed=false), network(4; directed=false)]
        add_edge!(un[1], 1, 2); add_edge!(un[2], 1, 2); add_edge!(un[2], 2, 3)
        @test TERGM._labels([Edges(), GWESP(0.5), GWDegree(0.5)], un[1]) ==
              ["edges", "gwesp.fixed.0.5", "gwdeg.fixed.0.5"]

        # p-values are Networks.jl's ONE z → p helper (floored, NaN-aware); a
        # coefficient fixed at ∓Inf (SE 0) prints p = 0, as R does
        fz = r.formation_coef ./ r.formation_se
        @test TERGM.z_pvalues === Networks.z_pvalues
        p = TERGM._pvalues(r.formation_coef, fz)
        @test all(0 .<= p .<= 1)
        @test 0 < p[1] <= 1
        @test p[2] == 0.0
    end

    # ------------------------------------------------------------------
    # Expanding terms (panel 2026-09 round 3): a multi-level NodeFactor /
    # NodeMix and a Degree range expand into one statistic per level / cell
    # / degree with R's labels — through ERGM's public `_materialize`, the
    # step `ERGMModel` runs — instead of the raw term's pooled single column
    # under a label R never prints. Attribute terms are snapshotted per
    # transition, against the auxiliary network the statistics are
    # evaluated on, so a time-varying attribute is honoured.
    # ------------------------------------------------------------------
    @testset "Expanding terms materialize per statistic with R's labels" begin
        rng = Random.Xoshiro(23)
        n = 12
        nets = Network{Int,true}[]
        for t in 1:4
            net = network(n; directed=true)
            for v in 1:n
                set_vertex_attribute!(net, :grp, v, ("a", "b", "c")[mod1(v, 3)])
            end
            for i in 1:n, j in 1:n
                i == j && continue
                p = t == 1 ? 0.25 : (has_edge(nets[t-1], i, j) ? 0.7 : 0.12)
                rand(rng) < p && add_edge!(net, i, j)
            end
            push!(nets, net)
        end

        # (a) NodeFactor on a 3-level attribute: K−1 columns under R's names,
        # on either side, exactly what ERGM.jl's fit_ergm produces
        fit = stergm(nets, [Edges()], [Edges(), NodeFactor(:grp)])
        @test coeftable(fit).names == ["Form(1)~edges", "Persist(1)~edges",
                                       "Persist(1)~nodefactor.grp.b",
                                       "Persist(1)~nodefactor.grp.c"]
        @test length(persistence_coef(fit)) == 3
        @test fit.converged && isempty(approximations(fit))
        @test coeftable(fit).names[3:4] ==
              "Persist(1)~" .* ERGMModel(ERGMFormula([NodeFactor(:grp)]), nets[1]).formula.terms.names
        # the stored formula IS the expansion (as `ERGMModel.formula` is), in
        # plain single-statistic terms; the user's raw formula prints the
        # unexpanded label
        @test typeof.(fit.model.formula.dissolution) == [Edges, NodeFactor, NodeFactor]
        @test TERGM._labels(fit.model.formula.dissolution, nets[1]) ==
              ["edges", "nodefactor.grp.b", "nodefactor.grp.c"]
        @test sprint(show, STERGM([Edges()], [Edges(), NodeFactor(:grp)])) ==
              "STERGM(formation: edges; persistence: edges + nodefactor.grp)"
        @test occursin("persistence: edges + nodefactor.grp.b + nodefactor.grp.c",
                       sprint(show, fit.model))
        # ... and the fit is exactly the fit of the hand-expanded formula
        by_hand = stergm(nets, [Edges()],
                         [Edges(), NodeFactor(:grp; level="b"), NodeFactor(:grp; level="c")])
        @test coef(by_hand) == coef(fit)
        @test stderror(by_hand) == stderror(fit)
        @test loglikelihood(by_hand) == loglikelihood(fit)
        # the raw term's single column was the SUM of the two (the pooled
        # model the package used to fit silently)
        yminus = dissolution_network(nets[1], nets[2])
        @test compute(NodeFactor(:grp), yminus) ==
              compute(NodeFactor(:grp; level="b"), yminus) +
              compute(NodeFactor(:grp; level="c"), yminus)
        @test !has_dyad_dependent(fit.model)
        @test is_exact(fit)

        # (b) NodeMix: every selected cell, in statnet's order (first cell
        # dropped), directed 3 levels → 8 cells
        mix = stergm(nets, [Edges(), NodeMix(:grp)], [Edges()])
        @test coeftable(mix).names[2:9] == "Form(1)~" .* [
            "mix.grp.b.a", "mix.grp.c.a", "mix.grp.a.b", "mix.grp.b.b",
            "mix.grp.c.b", "mix.grp.a.c", "mix.grp.b.c", "mix.grp.c.c"]
        @test coeftable(mix).names[2:9] ==
              "Form(1)~" .* ERGMModel(ERGMFormula([NodeMix(:grp)]), nets[1]).formula.terms.names
        @test length(formation_coef(mix)) == 9
        # a resolved single cell stays a single column
        @test coeftable(stergm(nets, [Edges(), NodeMix(:grp, "a", "b")], [Edges()])).names ==
              ["Form(1)~edges", "Form(1)~mix.grp.a.b", "Persist(1)~edges"]

        # (c) degree ranges: Degree(0:2) on an undirected panel is three
        # `degree<d>` rows (it used to construct and then crash in the row
        # fill with ERGM's "Expand via ERGMModel" message); IDegree/ODegree
        # ranges on a directed one
        un = Network{Int,false}[]
        for t in 1:3
            net = network(10; directed=false)
            for i in 1:10, j in (i+1):10
                p = t == 1 ? 0.3 : (has_edge(un[t-1], i, j) ? 0.7 : 0.15)
                rand(rng) < p && add_edge!(net, i, j)
            end
            push!(un, net)
        end
        deg = stergm(un, [Edges(), Degree(1:2)], [Edges(), Degree(0:1)])
        @test coeftable(deg).names == ["Form(1)~edges", "Form(1)~degree1", "Form(1)~degree2",
                                       "Persist(1)~edges", "Persist(1)~degree0",
                                       "Persist(1)~degree1"]
        @test deg.model.formula.formation == AbstractERGMTerm[Edges(), Degree(1), Degree(2)]
        model_d = STERGMModel(STERGM([Edges(), Degree(0:2)], [Edges()]), un)
        @test occursin("formation: edges + degree0 + degree1 + degree2", sprint(show, model_d))
        @test length(TERGM._cmple_blocks(model_d)[1][1][1, :]) == 4
        odeg = STERGMModel(STERGM([Edges(), ODegree(0:1)], [Edges(), IDegree([1, 2])]), nets)
        @test TERGM._labels(odeg.formula.formation, nets[1]) == ["edges", "odegree0", "odegree1"]
        @test TERGM._labels(odeg.formula.dissolution, nets[1]) == ["edges", "idegree1", "idegree2"]

        # (d) levels are resolved from panel 1, as R resolves them from the
        # union of the panels: a level ABSENT from panel 3 is fine — its
        # column is zero on the transition that starts there — ...
        nets_c = [copy(w) for w in nets]
        for v in 1:n
            get_vertex_attribute(nets_c[3], :grp, v) == "c" &&
                set_vertex_attribute!(nets_c[3], :grp, v, "a")
        end
        model_c = STERGMModel(STERGM([Edges(), NodeFactor(:grp)], [Edges()]), nets_c)
        @test TERGM._labels(model_c.formula.formation, nets_c[1]) ==
              ["edges", "nodefactor.grp.b", "nodefactor.grp.c"]
        Xfc, _, _, _ = TERGM._cmple_blocks(model_c)
        @test all(iszero, Xfc[3][:, 3])          # transition 3→4 starts at panel 3: no "c"
        @test any(!iszero, Xfc[1][:, 3])
        # ... whereas a level panel 1 LACKS would be a column R fits and this
        # model does not have: refused, naming the side, the panel and the value
        nets_d = [copy(w) for w in nets]
        set_vertex_attribute!(nets_d[3], :grp, 1, "d")
        e = try
            STERGMModel(STERGM([Edges(), NodeFactor(:grp)], [Edges()]), nets_d)
            nothing
        catch err
            err
        end
        @test e isa ArgumentError
        @test startswith(e.msg, "formation model, panel 3: term 'nodefactor.grp'")
        @test occursin("(a, b, c)", e.msg)
        @test occursin("takes the value \"d\" on panel 3", e.msg)
        @test occursin("levels=", e.msg)
        e2 = try; STERGMModel(STERGM([Edges()], [Edges(), NodeMix(:grp)]), nets_d); nothing; catch err; err; end
        @test e2 isa ArgumentError && startswith(e2.msg, "dissolution model, panel 3: term 'mix.grp'")
        # a resolved single level, a NodeMatch or a NodeCov does not expand by
        # level, so a new value on a later panel is just data there
        @test STERGMModel(STERGM([Edges(), NodeFactor(:grp; level="b"), NodeMatch(:grp)],
                                 [Edges()]), nets_d) isa STERGMModel
        # a single-level attribute cannot expand at all (ERGM's message, prefixed)
        one = [copy(w) for w in nets]
        for w in one, v in 1:n
            set_vertex_attribute!(w, :grp, v, "a")
        end
        e1 = try; STERGMModel(STERGM([Edges()], [Edges(), NodeFactor(:grp)]), one); nothing; catch err; err; end
        @test e1 isa ArgumentError
        @test startswith(e1.msg, "dissolution model, panel 1")
        @test occursin("no levels remain", e1.msg)

        # (e) the snapshot is per transition: transition t's rows read the
        # attribute of panel t−1 (the auxiliary networks carry it), so a
        # time-varying attribute enters the design as the raw term would
        # have read it
        tv = [copy(w) for w in nets]
        for v in 1:n
            set_vertex_attribute!(tv[2], :grp, v, ("c", "a", "b")[mod1(v, 3)])   # relabelled at t=2
        end
        model_tv = STERGMModel(STERGM([Edges(), NodeMatch(:grp)], [Edges(), NodeFactor(:grp)]), tv)
        Xf, _, Xd, _ = TERGM._cmple_blocks(model_tv)
        for t in 2:4
            yplus = formation_network(tv[t-1], tv[t])
            yminus = dissolution_network(tv[t-1], tv[t])
            @test Xf[t-1][:, 2] == [change_stat(NodeMatch(:grp), yplus, i, j)
                                    for i in 1:n for j in 1:n if i != j && !has_edge(tv[t-1], i, j)]
            @test Xd[t-1][:, 2] == [change_stat(NodeFactor(:grp; level="b"), yminus, i, j)
                                    for i in 1:n for j in 1:n if i != j && has_edge(tv[t-1], i, j)]
        end

        # (f) simulation and gof see the expansion: a raw multi-level formula
        # takes one coefficient per level, the fitted model's expanded formula
        # passes through unchanged, and gof labels the expanded statistics
        raw = STERGM([Edges()], [Edges(), NodeFactor(:grp)])
        y = simulate_stergm(nets[end], raw, [-2.0], [1.0, 0.3, -0.3]; burnin=200, rng=rng)
        @test y isa Network{Int,true}
        @test simulate_stergm(fit, 2; burnin=200, rng=rng) isa Vector{Network{Int,true}}
        g = gof(fit; n_sim=3, burnin=100, rng=rng)
        @test g.statistics[3].labels == ["edges", "nodefactor.grp.b", "nodefactor.grp.c"]
        obs = zeros(1); obs_d = zeros(3)
        TERGM._transition_stats!(obs, obs_d, fit.model.formula.formation,
                                 fit.model.formula.dissolution, nets[1], nets[2])
        @test obs_d == [compute(t, dissolution_network(nets[1], nets[2]))
                        for t in (Edges(), NodeFactor(:grp; level="b"), NodeFactor(:grp; level="c"))]
    end

    @testset "simulate_stergm names a coefficient-count mismatch and an absent attribute" begin
        nets = fixture_panels()
        f = STERGM([Edges()], [Edges()])
        msg_of(g) = try; g(); "no error"; catch e; e isa ArgumentError ? e.msg : rethrow(); end
        m = msg_of(() -> simulate_stergm(nets[1], f, [1.0, 2.0], [1.0]))
        @test occursin("formation model has 1 term(s) (edges)", m)
        @test occursin("2 formation coefficient(s) were given ([1.0, 2.0])", m)
        @test occursin("formation_coef(fit)", m)
        @test occursin("not the stacked coef(fit)", m)
        m = msg_of(() -> simulate_stergm(nets[1], f, [1.0], Float64[]))
        @test occursin("dissolution (persistence) model has 1 term(s) (edges)", m)
        @test occursin("0 dissolution (persistence) coefficient(s) were given", m)
        @test occursin("persistence_coef(fit)", m)
        # the expanded count is what is asked for
        for w in nets
            set_vertex_attribute!(w, :grp, Dict(v => (v <= 2 ? "a" : "b") for v in 1:5))
        end
        m = msg_of(() -> simulate_stergm(nets[1], STERGM([Edges(), NodeMix(:grp)], [Edges()]),
                                         [1.0], [1.0]))
        @test occursin("formation model has 4 term(s) (edges, mix.grp.b.a, mix.grp.a.b, mix.grp.b.b)", m)
        # a nodal term whose attribute the starting network lacks is refused
        # in ERGM's words (not a KeyError deep in the sampler)
        m = msg_of(() -> simulate_stergm(network(5; directed=true),
                                         STERGM([Edges(), NodeMatch(:grp)], [Edges()]),
                                         [1.0, 1.0], [1.0]))
        @test startswith(m, "simulate_stergm: formation model: term 'nodematch.grp'")
        @test occursin("Available vertex attributes: (none)", m)
    end

    @testset "Persistence accessors; dissolution_* deprecated (item 31)" begin
        nets = fixture_panels()
        r = stergm(nets, [Edges(), Delrecip()], [Edges()])
        @test formation_coef(r) === r.formation_coef
        @test formation_se(r) === r.formation_se
        @test persistence_coef(r) === r.persistence_coef
        @test persistence_se(r) === r.persistence_se
        @test length(persistence_coef(r)) == 1
        @test coef(r) == vcat(formation_coef(r), persistence_coef(r))
        @test stderror(r) == vcat(formation_se(r), persistence_se(r))
        @test :persistence_coef in propertynames(r)
        @test !(:dissolution_coef in propertynames(r))

        # The deprecated spellings return the very same vectors — and WARN
        # (Pkg.test runs with --depwarn=yes, so the warnings are observable;
        # `@test_deprecated` is a no-op under --depwarn=no and the equalities
        # still hold there)
        @test (@test_deprecated dissolution_coef(r)) === persistence_coef(r)
        @test (@test_deprecated dissolution_se(r)) === persistence_se(r)
        @test (@test_logs (:warn, r"deprecated") match_mode=:any r.dissolution_coef) ===
              persistence_coef(r)
        @test (@test_logs (:warn, r"deprecated") match_mode=:any r.dissolution_se) ===
              persistence_se(r)
        # the shim leaves normal field access alone: no warning on the real
        # fields, and the property names are the real ones
        @test_logs r.persistence_coef
        @test_logs r.formation_coef
        @test propertynames(r) == fieldnames(STERGMResult)

        # The persistence log-odds are what an edges-only model reproduces:
        # 4 prior ties, 2 persisted → logit(1/2) = 0
        @test persistence_coef(stergm(nets, [Edges()], [Edges()]))[1] ≈ 0.0 atol = 1e-6
    end

    @testset "Block bootstrap runs on Networks.bootstrap_cov" begin
        panels = random_panels(Random.Xoshiro(99))
        r_b = stergm(panels, [Edges()], [Edges()];
                     se=:bootstrap, n_boot=40, rng=Random.Xoshiro(1))
        @test size(r_b.boot_replicates) == (40, 2)
        @test all(isfinite, r_b.boot_replicates)
        @test TERGM._n_dropped_replicates(r_b) == 0
        # The reported covariance IS the empirical covariance of the replicates
        @test vcov(r_b) ≈ cov(r_b.boot_replicates)
        @test stderror(r_b) ≈ sqrt.(diag(cov(r_b.boot_replicates)))
        @test isempty(approximations(r_b))
        # A Hessian fit carries no replicates
        r_h = stergm(panels, [Edges()], [Edges()])
        @test size(r_h.boot_replicates) == (0, 2)

        # Thread-count independence: the resamples are drawn from `rng` in
        # one call and the refits are deterministic, so a serial run of the
        # same shared loop gives the same replicates
        Xf, yf, Xd, yd = TERGM._cmple_blocks(r_h.model)
        serial = Networks.bootstrap_cov(
            idx -> vcat(TERGM._logistic_fit(TERGM._stack(Xf[idx], 1), TERGM._stack(yf[idx])).coef,
                        TERGM._logistic_fit(TERGM._stack(Xd[idx], 1), TERGM._stack(yd[idx])).coef),
            (rng, B) -> [rand(rng, 1:length(Xf), length(Xf)) for _ in 1:B],
            coef(r_h); n_boot=40, rng=Random.Xoshiro(1), threaded=false)
        @test serial.replicates == r_b.boot_replicates

        # A replicate without a finite refit is excluded, warned about ONCE
        # and recorded: panels 2 and 3 are empty, so a resample made only of
        # transitions 2 and 3 has no dissolution rows at all
        t0 = network(6); add_edge!(t0, 1, 2); add_edge!(t0, 2, 3); add_edge!(t0, 3, 4)
        t1 = network(6)
        t2 = network(6)
        t3 = network(6); add_edge!(t3, 1, 2); add_edge!(t3, 4, 5); add_edge!(t3, 5, 6)
        t4 = network(6); add_edge!(t4, 1, 2); add_edge!(t4, 4, 5); add_edge!(t4, 2, 6)
        logs, r_drop = Test.collect_test_logs() do
            stergm([t0, t1, t2, t3, t4], [Edges()], [Edges()];
                   se=:bootstrap, n_boot=60, rng=Random.Xoshiro(3))
        end
        excluded = [l for l in logs if l.level == Base.CoreLogging.Warn &&
                                       occursin("excluded", string(l.message))]
        @test length(excluded) == 1
        n_drop = TERGM._n_dropped_replicates(r_drop)
        @test n_drop > 0
        @test all(isfinite, stderror(r_drop))
        @test any(occursin("$n_drop refits", a) for a in approximations(r_drop))
        @test occursin("$n_drop of 60 block-bootstrap refits", sprint(show, r_drop))
        ok = [all(isfinite, r_drop.boot_replicates[b, :]) for b in 1:60]
        @test vcov(r_drop) ≈ cov(r_drop.boot_replicates[ok, :])

        # `se=` is validated by Networks.check_se, so the message has the
        # ecosystem's one shape
        msg = try; cmple(r_h.model; se=:sandwich); catch e; e.msg; end
        @test occursin("cmple: se must be one of (:hessian, :bootstrap)", msg)
        @test occursin(":sandwich", msg)
    end

    @testset "Sampler defaults are dyad-scaled; kernel is ERGM.mh_toggle!" begin
        # Default burnin resolves to ERGM's `_mcmc_defaults` for the FREE dyads
        # of each side: 20 toggles per free dyad
        @test ERGM._mcmc_defaults(16).burnin == 320
        rng = Random.Xoshiro(17)
        prev = network(12)
        for i in 1:12, j in 1:12
            i != j && rand(rng) < 0.2 && add_edge!(prev, i, j)
        end
        set_vertex_attribute!(prev, :grp, Dict(v => (isodd(v) ? "a" : "b") for v in 1:12))
        n_free_form = 12 * 11 - ne(prev)
        formula = STERGM([Edges()], [Edges()])
        # An explicit burnin equal to the resolved default reproduces the
        # default draw exactly (same rng, same number of toggles per side)
        y_default = simulate_stergm(prev, formula, [-2.0], [1.0]; rng=Random.Xoshiro(5))
        y_explicit_form = TERGM._sample_constrained(prev, formula.formation, [-2.0],
                                                    :formation,
                                                    ERGM._mcmc_defaults(n_free_form).burnin,
                                                    Random.Xoshiro(5))
        y_default_form = TERGM._sample_constrained(prev, formula.formation, [-2.0],
                                                   :formation, nothing, Random.Xoshiro(5))
        edgeset(g) = sort([(Int(src(e)), Int(dst(e))) for e in edges(g)])
        @test edgeset(y_explicit_form) == edgeset(y_default_form)
        @test y_default isa Network{Int, true}
        @test_throws ArgumentError simulate_stergm(prev, formula, [-2.0], [1.0]; burnin=-1)

        # The constrained sampler is bit-identical to the pre-0.2 hand loop: a
        # reference implementation of that loop, run on the same rng
        function reference_loop(prev, terms, θ, constrain, steps, rng)
            net = TERGM._copy_net(prev)
            n = Int(nv(net))
            free = NTuple{2,Int}[]
            for i in 1:n, j in 1:n
                i == j && continue
                if constrain == :formation
                    has_edge(prev, i, j) || push!(free, (i, j))
                else
                    has_edge(prev, i, j) && push!(free, (i, j))
                end
            end
            delta = zeros(length(terms))
            for _ in 1:steps
                i, j = free[rand(rng, 1:length(free))]
                for (k, t) in enumerate(terms)
                    delta[k] = TERGM._tchange(t, net, i, j, prev)
                end
                la = dot(θ, delta)
                has_edge(net, i, j) && (la = -la)
                if log(rand(rng)) < la
                    has_edge(net, i, j) ? rem_edge!(net, i, j) : add_edge!(net, i, j)
                end
            end
            return net
        end
        # ... with an attribute term in the formula: the reference loop reads the
        # raw `NodeMatch` through the attribute Dict, the sampler its
        # materialized twin — same draws, same statistics
        terms = [Edges(), Delrecip(), NodeMatch(:grp)]
        θ = [-1.5, 0.8, 0.6]
        for side in (:formation, :dissolution)
            a = TERGM._sample_constrained(prev, terms, θ, side, 500, Random.Xoshiro(11))
            b = reference_loop(prev, terms, θ, side, 500, Random.Xoshiro(11))
            @test edgeset(a) == edgeset(b)
        end

        # ... and allocation-free per step (the kernel and the typed closures),
        # attribute term included (round 3: the raw `NodeMatch` cost 112 B per
        # toggle through the attribute Dict; the materialized twin costs 0).
        # Measured at a coefficient no proposal is accepted at, so the only
        # per-step work is the kernel's: an accepted toggle grows the network's
        # adjacency lists, which is the data structure's cost, not the sampler's.
        θ_never = [-Inf, 0.0, 0.0]
        net = copy(prev)
        free = [(i, j) for i in 1:Int(nv(prev)) for j in 1:Int(nv(prev))
                if i != j && !has_edge(prev, i, j)]
        snapshots = TERGM._materialized_tuple(terms, prev)
        rng = Random.Xoshiro(1)
        # Measure the kernel on prepared inputs. Copying attribute Dicts in
        # `_sample_constrained` has a small, hash-layout-dependent setup cost.
        free_steps(steps) = (TERGM._mh_constrained!(rng, net, prev, free, snapshots,
                                                   θ_never, steps); nothing)
        free_steps(10)
        a1 = @allocated free_steps(1_000)
        a2 = @allocated free_steps(2_000)
        @test a2 == a1
        @test edgeset(net) == edgeset(prev)
    end

    # ------------------------------------------------------------------
    # Cross-package hygiene (panel 2026-09, item 13): every `_`-prefixed name
    # this package takes from ERGM.jl or Networks.jl must be declared `public`
    # there, and no `ERGM._x` / `Networks._x` dotted reach-in may remain in
    # the source — shared helpers are imported by name.
    # ------------------------------------------------------------------
    @testset "No private cross-package reach-ins" begin
        src_text = read(joinpath(@__DIR__, "..", "src", "TERGM.jl"), String)
        mods = Dict("ERGM" => ERGM, "Networks" => Networks)

        # Code only: docstrings and comments may *mention* `ERGM._validate_formula`
        code = replace(src_text, r"\"\"\"[\s\S]*?\"\"\"" => "")
        code = join((replace(l, r"#.*$" => "") for l in split(code, '\n')), '\n')
        dotted = [(m.captures[1], "_" * m.captures[2])
                  for m in eachmatch(r"\b(ERGM|Networks)\._([A-Za-z0-9_!]+)", code)]
        @test isempty(dotted)

        # Names imported by name: `import ERGM: a, b,` continued over lines
        imported = Tuple{String,String}[]
        lines = split(src_text, '\n')
        k = 1
        while k <= length(lines)
            m = match(r"^import (ERGM|Networks):\s*(.*)$", lines[k])
            if m !== nothing
                rest = replace(String(m.captures[2]), r"#.*$" => "")
                while endswith(strip(rest), ",") && k < length(lines)
                    k += 1
                    rest *= " " * replace(strip(lines[k]), r"#.*$" => "")
                end
                for nm in split(rest, ',')
                    s = String(strip(nm))
                    isempty(s) && continue
                    push!(imported, (String(m.captures[1]), s))
                end
            end
            k += 1
        end
        @test !isempty(imported)
        # Every imported `_`-name must be `public` (or exported) where it
        # lives. The allow-list of names imported while a cross-repo `public`
        # request was pending (`_mple_fit_design`, `_warn_boundary`,
        # `_warn_separated`, `_collect_terms`) is empty since the 2026-09
        # reconciliation declared them; the marker comment went with it.
        for (modname, nm) in imported
            mod = mods[modname]
            sym = Symbol(nm)
            @test isdefined(mod, sym)
            startswith(nm, "_") && @test Base.ispublic(mod, sym)
        end
        @test !occursin("cross-repo request pending", src_text)
        for nm in ("_validate_formula", "_materialize", "_mcmc_defaults",
                   "_mple_fit_design", "_collect_terms", "_expand_terms")
            @test ("ERGM", nm) in imported
        end
        # The two R-sentence warnings are `_mple_fit_design`'s own (through its
        # `context="cmple"`); TERGM neither imports nor re-emits them
        @test !(("ERGM", "_warn_boundary") in imported)
        @test !(("ERGM", "_warn_separated") in imported)
        # ... and no local twin of ERGM's term expansion remains
        @test !isdefined(TERGM, :_specification) && !isdefined(TERGM, :_polish_newton)
        @test TERGM._expand_terms === ERGM._expand_terms
        # The deprecated private aliases are gone from this package
        # (as whole identifiers: `set_vertex_attribute!` is not `_vertex_attribute`)
        for stale in ("_z_pvalues", "_requires_directed", "_vertex_attribute",
                      "_has_dyad_dependent")
            @test !occursin(Regex("(?<![A-Za-z0-9_!])" * stale * "(?![A-Za-z0-9_!])"),
                            src_text)
        end
    end

    # ------------------------------------------------------------------
    # Nothing about a bad fit is silent (panel 2026-09, item 24 applied to
    # CMPLE): an exhausted Newton cap, a separated design and an empty side
    # are warned about at fit time, recorded in `approximations`, printed by
    # `show`, and make `is_exact` false.
    # ------------------------------------------------------------------
    @testset "Unconverged and degenerate CMPLE fits are loud" begin
        panels = random_panels(Random.Xoshiro(99))

        # (a) the Newton cap: warned (naming maxiter), recorded, printed
        unconv = @test_logs (:warn, r"did not converge in maxiter=1") match_mode=:any stergm(
            panels, [Edges()], [Edges()]; maxiter=1)
        @test !unconv.converged
        @test any(occursin("did not converge", a) for a in approximations(unconv))
        @test !is_exact(unconv)
        @test occursin("did not converge", sprint(show, unconv))
        @test !fit_metadata(unconv).is_exact
        # ... and a converged fit says nothing of the sort
        ok = stergm(panels, [Edges()], [Edges()])
        @test ok.converged
        @test !any(occursin("did not converge", a) for a in approximations(ok))
        @test !occursin("Warning:", sprint(show, ok))
        @test is_exact(ok)

        # (b) perfect separation by a COMBINATION of statistics (which the
        # boundary test cannot see): on this panel no cross-group tie persists,
        # so persistence edges → -Inf and nodematch → +Inf with their sum
        # finite. R tergm warns "The MPLE does not exist!"; so does cmple, and
        # the fit comes back `converged == false` rather than as the point where
        # Newton met its tolerance on the asymptote.
        grp = Dict(v => (isodd(v) ? "a" : "b") for v in 1:6)
        mk() = (net = network(6); set_vertex_attribute!(net, :grp, grp); net)
        u0 = mk()
        for (i, j) in ((1, 3), (3, 5), (2, 4), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6))
            add_edge!(u0, i, j)
        end
        u1 = mk()          # cross-group prior ties all dissolve; (1,3),(2,4) persist
        for (i, j) in ((1, 3), (2, 4), (3, 1), (6, 4), (2, 6))
            add_edge!(u1, i, j)
        end
        u2 = mk()          # (2,6) is cross-group and dissolves; (1,3),(6,4) persist
        for (i, j) in ((1, 3), (6, 4), (1, 4), (5, 2))
            add_edge!(u2, i, j)
        end
        sep = @test_logs (:warn, r"cmple: the MPLE does not exist \(perfect separation\)") match_mode=:any stergm(
            [u0, u1, u2], [Edges()], [Edges(), NodeMatch(:grp)])
        @test !sep.converged
        @test !is_exact(sep)
        @test any(occursin("does not exist", a) for a in approximations(sep))
        @test occursin("does not exist", sprint(show, sep))
        # Only the persistence side is separated: the formation side is a
        # finite, converged fit (7 of 22 + 6 of 25 prior non-ties formed)
        @test isfinite(formation_coef(sep)[1])
        @test all(isfinite, persistence_coef(sep))    # the last iterate, not ±Inf
        @test size(sep.boot_replicates) == (0, 3)

        # (c) an unidentified design is refused at construction now (see
        # "Memory terms are refused in a separable formula"); an identical
        # column made by hand — `EdgeCov` of the previous panel's complement
        # is −1 on every prior non-tie, i.e. −edges — still comes back as a
        # loud non-convergence rather than `converged == true` with standard
        # errors of 1e7 (which is what it did before 0.2)
        n_c = nv(panels[1])
        cov = [Float64(i != j && !has_edge(panels[1], i, j)) for i in 1:n_c, j in 1:n_c]
        two = panels[1:2]
        coll = @test_logs (:warn, r"did not converge") match_mode=:any stergm(
            two, [Edges(), EdgeCov(cov; name="notie")], [Edges()])
        @test !coll.converged
        @test !is_exact(coll)
        @test any(occursin("did not converge", a) for a in approximations(coll))

        # (d) a side with no free dyad on any transition: NaN coefficients,
        # said in so many words
        e0 = network(4)
        e1 = network(4); add_edge!(e1, 1, 2)
        emp = @test_logs (:warn, r"dissolution model has no free dyads") match_mode=:any stergm(
            [e0, e1], [Edges()], [Edges()])
        @test !emp.converged
        @test all(isnan, persistence_coef(emp))
        @test isfinite(formation_coef(emp)[1])
    end

    # ------------------------------------------------------------------
    # Frozen pre-refactor sampler output. Before `_sample_constrained` was
    # ported to ERGM's `mh_toggle!` kernel (2026-09), this exact draw was
    # taken with the committed hand-written loop (git 42dd0ee, 2026-07-17):
    # the kernel draws the proposal index, then one uniform, in the same
    # order, so the output must stay bit-identical. A literal, so the pin
    # does not depend on a reference loop written in the same sprint.
    # ------------------------------------------------------------------
    @testset "Frozen pre-refactor sampler output" begin
        rng = Random.Xoshiro(17)
        prev = network(12)
        for i in 1:12, j in 1:12
            i != j && rand(rng) < 0.2 && add_edge!(prev, i, j)
        end
        @test ne(prev) == 31
        y = simulate_stergm(prev, STERGM([Edges(), Mutual()], [Edges()]),
                            [-2.0, 0.5], [1.0]; burnin=500, rng=Random.Xoshiro(17))
        es = sort([(Int(src(e)), Int(dst(e))) for e in edges(y)])
        frozen = [(2, 6), (2, 7), (2, 12), (3, 4), (3, 6), (3, 9), (4, 1), (4, 5),
                  (4, 9), (4, 11), (5, 1), (5, 3), (5, 6), (5, 8), (6, 1), (6, 2),
                  (6, 7), (6, 10), (7, 2), (7, 5), (7, 8), (7, 9), (7, 10), (8, 3),
                  (8, 7), (9, 1), (9, 3), (9, 4), (9, 7), (9, 12), (10, 3), (10, 4),
                  (10, 12), (11, 6), (12, 1), (12, 7), (12, 9)]
        @test es == frozen
    end

    # ------------------------------------------------------------------
    # Keyword vocabulary (panel 2026-09, item 16): `maxiter`, `n_sim`,
    # `n_boot`, `rng` — never `max_iter`, `n_sims`, `seed`, ...
    # ------------------------------------------------------------------
    @testset "Keyword vocabulary follows the root CLAUDE.md" begin
        banned = (:max_iter, :max_iterations, :n_sims, :n_simulations, :seed,
                  :nboot, :n_bootstrap)
        own(f) = [m for m in methods(f) if m.module === TERGM]
        for f in (cmple, stergm, gof, simulate_stergm, simulate_network_sequence)
            ms = own(f)
            @test !isempty(ms)
            for m in ms
                kw = Base.kwarg_decl(m)
                @test isempty(intersect(kw, banned))
            end
        end
        @test :maxiter in Base.kwarg_decl(only(own(cmple)))
        @test :rng in Base.kwarg_decl(only(own(cmple)))
        @test :n_boot in Base.kwarg_decl(only(own(cmple)))
        gof_kw = Base.kwarg_decl(only(own(gof)))
        @test :n_sim in gof_kw && :rng in gof_kw && :burnin in gof_kw
        for m in own(simulate_stergm)
            kw = Base.kwarg_decl(m)
            # the (result, n_steps; kwargs...) forwarder declares only the splat
            any(endswith(string(k), "...") for k in kw) && continue
            @test :burnin in kw && :rng in kw
        end
    end

    # ------------------------------------------------------------------
    # `gof` draws one seed per (transition, simulation) from `rng` up front
    # and runs the simulations on every thread: reproducible from `rng`
    # alone and thread-count independent.
    # ------------------------------------------------------------------
    @testset "gof is seeded per simulation and thread-count independent" begin
        panels = random_panels(Random.Xoshiro(20); n=8, T=4)
        fit = stergm(panels, [Edges()], [Edges()])
        g1 = gof(fit; n_sim=12, burnin=50, rng=Random.Xoshiro(41))
        g2 = gof(fit; n_sim=12, burnin=50, rng=Random.Xoshiro(41))
        @test [s.simulated for s in g1.statistics] == [s.simulated for s in g2.statistics]
        @test [s.p_values for s in g1.statistics] == [s.p_values for s in g2.statistics]
        # A serial reproduction of the same seeding scheme (the only place
        # randomness enters) gives the very same simulated statistics
        n_trans = length(panels) - 1
        seeds = rand(Random.Xoshiro(41), UInt64, n_trans * 12)
        serial = zeros(12, 2)
        for k in 1:(n_trans * 12)
            t, s = divrem(k - 1, 12) .+ 1
            sim = simulate_stergm(panels[t], fit.model.formula, formation_coef(fit),
                                  persistence_coef(fit); burnin=50,
                                  rng=Random.Xoshiro(seeds[k]))
            serial[s, 1] += compute(NewEdge(), sim, panels[t])
            serial[s, 2] += compute(PersistentEdge(), sim, panels[t])
        end
        @test g1.statistics[1].simulated == serial
        @test_throws ArgumentError gof(fit; n_sim=0)
    end

    @testset "Per-transition block-bootstrap SEs (btergm-style)" begin
        panels = random_panels(Random.Xoshiro(99))

        r_h = stergm(panels, [Edges()], [Edges()])
        r_b = stergm(panels, [Edges()], [Edges()];
                     se=:bootstrap, n_boot=60, rng=Random.Xoshiro(1))

        # Bootstrap changes only the uncertainty, not the point estimates
        @test r_b.formation_coef == r_h.formation_coef
        @test r_b.persistence_coef == r_h.persistence_coef
        @test r_h.se_type == :hessian
        @test r_b.se_type == :bootstrap

        @test all(isfinite, stderror(r_b))
        @test all(stderror(r_b) .> 0)
        V = vcov(r_b)
        @test size(V) == (2, 2)
        @test V ≈ V'
        @test sqrt(V[1, 1]) ≈ r_b.formation_se[1]
        @test sqrt(V[2, 2]) ≈ r_b.persistence_se[1]

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

        # Both blocks render through the shared Networks.jl coefficient
        # printer (R-style columns, one significance-code legend), under the
        # two model names — `Persistence:`, tergm's Persist() parameterisation
        # (the pre-0.2 header said "Dissolution"), and the header says how the
        # standard errors were obtained
        @test occursin("Formation:", out)
        @test occursin("Persistence:", out)
        @test !occursin("Dissolution", out)
        @test count("Pr(>|z|)", out) == 2
        @test count("Signif. codes:", out) == 1
        @test occursin("Converged: true", out)
        @test occursin("inverse Hessian", out)
        @test occursin("Panels: $(length(panels)) ($(length(panels) - 1) transitions", out)
        # the AIC/BIC line ERGMResult's header carries, with the same numbers
        # the StatsAPI verbs report
        @test occursin("AIC: $(round(aic(r_dep), digits=2)), BIC: $(round(bic(r_dep), digits=2))", out)
        # the two printed blocks are the two halves of coeftable(r_dep)
        tbl = coeftable(r_dep)
        for row in tbl
            @test occursin(string(round(row.estimate, digits=4)), out)
        end

        # Dyad-independent formula (incl. temporal terms) → no caveat
        r_ind = stergm(panels, [Edges(), Delrecip()], [Edges()])
        @test !occursin("dyad-dependent", sprint(show, r_ind))

        # Bootstrap SEs → milder note, no anticonservative warning
        r_boot = stergm(panels, [Edges()], [Edges(), Mutual()];
                        se=:bootstrap, n_boot=30, rng=Random.Xoshiro(7))
        ob = sprint(show, r_boot)
        @test occursin("block-bootstrap", ob)
        @test occursin("block bootstrap (30 replicates)", ob)
        @test !occursin("anticonservative", ob)

        # An unconverged fit is loud in the HEADER, not in a footnote: the
        # verdict line reads false and the WP2 caveat sentence follows it
        # before any table is printed
        unconv = @test_logs (:warn, r"did not converge") match_mode=:any stergm(
            panels, [Edges()], [Edges()]; maxiter=1)
        ou = sprint(show, unconv)
        @test occursin("Converged: false", ou)
        @test occursin("WARNING:", ou)
        @test findfirst("WARNING:", ou).start < findfirst("Formation:", ou).start
        @test occursin("did not converge", ou)
        @test !occursin("Converged: false", out)
        @test !occursin("WARNING:", out)
    end

    @testset "Input validation" begin
        @test_throws ArgumentError STERGM(AbstractERGMTerm[], [Edges()])
        @test_throws ArgumentError STERGMModel(STERGM([Edges()], [Edges()]),
                                               [network(3)])
        nets = [network(3), network(4)]
        @test_throws ArgumentError STERGMModel(STERGM([Edges()], [Edges()]), nets)
        @test fit_stergm === stergm
    end

    @testset "cmle does not masquerade as CMPLE" begin
        # `cmle` used to warn and return a CMPLE fit stamped `:cmle`, so the
        # result object misreported its own estimator. It now throws.
        nets = [network(5), network(5)]
        add_edge!(nets[1], 1, 2)
        add_edge!(nets[2], 1, 2); add_edge!(nets[2], 2, 3)
        model = STERGMModel(STERGM([Edges()], [Edges()]), nets)

        @test_throws ArgumentError cmle(model)
        @test_throws ArgumentError stergm(nets, [Edges()], [Edges()]; method=:cmle)

        # The message must point at the estimator that does exist
        msg = try
            cmle(model)
        catch e
            sprint(showerror, e)
        end
        @test occursin("cmple", msg)

        # EGMME is unexported (it can only throw, so it must not be
        # advertised) but stays reachable — and still raises rather than
        # returning placeholders.
        @test !(:egmme in names(TERGM))
        @test_throws ErrorException TERGM.egmme(model)

        # CMPLE itself still works and labels itself honestly
        res = cmple(model)
        @test res.method == :cmple
    end

    @testset "Missing dyads are rejected" begin
        # CMPLE enumerates every free dyad of Y⁺/Y⁻ as observed, so a masked
        # dyad would enter the design matrix at its face value. Reject it.
        nets = [network(5), network(5)]
        add_edge!(nets[1], 1, 2)
        add_edge!(nets[2], 1, 2)

        # Absent-face masked dyad
        set_missing_dyad!(nets[2], 3, 4)
        @test_throws ArgumentError STERGMModel(STERGM([Edges()], [Edges()]), nets)
        @test_throws ArgumentError stergm(nets, [Edges()], [Edges()])

        # The error names the offending panel
        msg = try
            STERGMModel(STERGM([Edges()], [Edges()]), nets)
        catch e
            sprint(showerror, e)
        end
        @test occursin("panel 2", msg)

        # Present-face masked dyad is rejected just the same
        clear_missing_dyads!(nets[2])
        set_missing_dyad!(nets[2], 1, 2)
        @test_throws ArgumentError STERGMModel(STERGM([Edges()], [Edges()]), nets)

        # Clearing the mask restores a fittable model
        clear_missing_dyads!(nets[2])
        @test STERGMModel(STERGM([Edges()], [Edges()]), nets) isa STERGMModel

        # The missing-data VOCABULARY (panel item 5): TERGM exposes no
        # `missing=` keyword, so its policy set is Networks' default
        # `(:error,)`, it does not claim `supports_missing`, and the refusal
        # must not advertise a `:face` policy the caller cannot pass
        # (`face_ok=false`)
        @test missing_policies(stergm) == (:error,)
        @test missing_policies(cmple) == (:error,)
        @test !supports_missing(stergm)
        @test !supports_missing(cmple)
        @test all(!(:missing in Base.kwarg_decl(m)) for m in methods(stergm))
        @test !(:missing in Base.kwarg_decl(only(methods(cmple))))
        set_missing_dyad!(nets[2], 3, 4)
        msg_face = try
            stergm(nets, [Edges()], [Edges()])
        catch e
            sprint(showerror, e)
        end
        @test !occursin(":face", msg_face)
        @test occursin("clear_missing_dyads!", msg_face)
        @test occursin("TERGM (panel 2)", msg_face)
        clear_missing_dyads!(nets[2])
    end

    @testset "Shared verbs are single-owner; no method ambiguities" begin
        # Every verb TERGM extends is THE binding of the package that owns it,
        # so `using ERGM, TERGM` (and SNA/REM beside them) leaves each of
        # them defined exactly once
        @test TERGM.gof === Networks.gof
        @test TERGM.coeftable === TERGM.StatsAPI.coeftable
        @test TERGM.confint === TERGM.StatsAPI.confint
        @test TERGM.name === ERGM.name
        @test TERGM.compute === ERGM.compute === Networks.compute
        @test TERGM.change_stat === ERGM.change_stat
        @test TERGM.has_dyad_dependent === ERGM.has_dyad_dependent
        @test TERGM.z_pvalues === Networks.z_pvalues
        @test TERGM.bootstrap_cov === Networks.bootstrap_cov
        @test TERGM.newton_fit === Networks.newton_fit === ERGM.newton_fit
        @test TERGM.mh_toggle! === ERGM.mh_toggle!
        @test isempty(Test.detect_ambiguities(TERGM))
        @test isempty(Test.detect_ambiguities(TERGM, ERGM, Networks))
        # both spellings of the entry point, one definition
        @test fit_stergm === stergm
    end

    # ------------------------------------------------------------------
    # Release engineering (panel 2026-09, item 29): the CI/Documentation
    # workflows must never regain a hand-maintained sibling list. The clone
    # step parses the `[sources]` tables with TOML, so the list is exactly
    # what the path dependencies point at.
    # ------------------------------------------------------------------
    @testset "Workflows derive their clone lists from [sources]" begin
        pkgdir = normpath(joinpath(@__DIR__, ".."))
        # Every [sources] entry is a sibling checkout `../<Pkg>.jl`, or the
        # package itself (`..`) from docs/
        expected = Set{String}()
        for (f, self_ok) in ((joinpath(pkgdir, "Project.toml"), false),
                             (joinpath(pkgdir, "docs", "Project.toml"), true))
            srcs = get(TOML.parsefile(f), "sources", Dict{String,Any}())
            @test !isempty(srcs)
            for (pkg, spec) in srcs
                path = spec["path"]
                if pkg == "TERGM"
                    @test self_ok && path == ".."
                else
                    @test occursin(Regex("^(\\.\\./)+" * pkg * "\\.jl\$"), path)
                    push!(expected, pkg)
                end
            end
        end
        @test expected == Set(["Networks", "ERGM"])   # what TERGM depends on
        @test !("NetworkDynamic" in expected)          # the over-clone item 29 named
        for wf in ("CI.yml", "Documentation.yml")
            yaml = read(joinpath(pkgdir, ".github", "workflows", wf), String)
            # no literal package list survives: the only `for pkg in` iterates
            # the shell variable filled from the TOML parse
            @test !occursin(r"for pkg in [A-Z]", yaml)
            @test occursin(r"for pkg in \$pkgs", yaml)
            @test occursin("TOML.parsefile", yaml)
            @test occursin("TERGM.jl/Project.toml", yaml)
            @test occursin("TERGM.jl/docs/Project.toml", yaml)
            # `julia` must be on the path before the clone step runs it
            @test findfirst("setup-julia", yaml) !== nothing
            @test first(findfirst("setup-julia", yaml)) <
                  first(findfirst("Reconstruct the ecosystem layout", yaml))
            # the clone is exactly one-per-sibling; no sibling name is spelled out
            for stale in ("NetworkDynamic", "SNA.jl", "REM.jl")
                @test !occursin(Regex("clone[^\\n]*" * stale), yaml)
            end
        end
        # The multithreaded cell (bootstrap/gof thread-independence testsets)
        ci = read(joinpath(pkgdir, ".github", "workflows", "CI.yml"), String)
        @test occursin("JULIA_NUM_THREADS", ci)
    end

    # ------------------------------------------------------------------
    # Usability (grade-A criterion 5): every export TERGM documents itself
    # carries a runnable `# Example` block. The StatsAPI verbs are documented
    # by StatsAPI; TERGM's own methods of `coeftable`/`confint`/`gof`/
    # `has_dyad_dependent` are checked through the module's own doc store.
    # ------------------------------------------------------------------
    @testset "Every export is documented with an example" begin
        meta = Base.Docs.meta(TERGM)
        own = Dict{Symbol,String}()
        for (binding, multidoc) in meta
            txt = join((string(get(d.data, :fields, ""), "\n",
                               join(d.text, "\n")) for d in values(multidoc.docs)), "\n")
            own[binding.var] = txt
        end
        statsapi = (:coef, :stderror, :vcov, :loglikelihood, :aic, :bic, :nobs, :dof)
        for n in names(TERGM)
            n in (:TERGM, :fit_stergm) && continue   # module; const alias of stergm
            n in statsapi && continue                 # documented by StatsAPI
            @test haskey(own, n)
            haskey(own, n) || continue
            if n === :stergm_gof
                @test occursin(r"deprecated"i, own[n])   # a deprecation, not a feature
            else
                @test occursin("# Example", own[n])
                @test occursin("```julia", own[n])
            end
        end
    end

    @testset "Result metadata protocol" begin
        panels = random_panels(Random.Xoshiro(20); n=10, T=4)

        # THE assertion this protocol exists for: ONE estimator (CMPLE), two
        # formulas. CMPLE of a dyad-independent formula IS the exact CMLE; the
        # same estimator on a formula with one dyad-dependent term is not.
        # (`Delrecip`, not `EdgeStability`: on the formation side every prior
        # non-tie's `EdgeStability` change is −1, i.e. the column is −Edges,
        # so that formula is unidentified — and since 0.2 the fit says so;
        # see "Unconverged and degenerate CMPLE fits are loud")
        indep = stergm(panels, [Edges(), Delrecip()], [Edges()])
        dep = stergm(panels, [Edges(), Triangle()], [Edges()])

        md_indep = fit_metadata(indep)
        @test md_indep.estimand == :stergm
        @test md_indep.objective == :conditional_pseudolikelihood
        @test md_indep.is_exact            # both formulas dyad-independent
        @test md_indep.se_method == :hessian
        @test md_indep.missing_method == :rejected
        @test isempty(md_indep.approximations)

        md_dep = fit_metadata(dep)
        @test md_dep.objective == :conditional_pseudolikelihood   # same estimator
        @test !md_dep.is_exact                                    # different formula
        @test any(occursin("anticonservative", a) for a in md_dep.approximations)

        # A dyad-dependent term in the DISSOLUTION formula alone is enough:
        # exactness needs every term of BOTH formulas to be dyad-independent
        dep_diss = stergm(panels, [Edges()], [Edges(), Triangle()])
        @test !is_exact(dep_diss)

        # `show`'s prose caveat and the protocol are driven by the same
        # predicate, so they agree on every fit
        for fit in (indep, dep, dep_diss)
            printed = sprint(show, fit)
            @test occursin("pseudolikelihood", printed) == !is_exact(fit)
        end

        # Bootstrap standard errors are reported as such, and do not change
        # whether the objective is exact
        boot = stergm(panels, [Edges(), Triangle()], [Edges()];
                      se=:bootstrap, n_boot=5, rng=Random.Xoshiro(21))
        @test se_method(boot) == :bootstrap
        @test !is_exact(boot)
        @test !any(occursin("anticonservative", a) for a in approximations(boot))
    end

    # ------------------------------------------------------------------
    # Golden fixture: a REAL statnet `tergm` CMPLE fit, with provenance
    # (issue #8). test/fixtures/r/panel_stergm.R regenerates it.
    #
    # The formulas are dyad-independent on purpose. Under dyad independence the
    # conditional pseudo-likelihood IS the conditional likelihood, so CMPLE is
    # the exact conditional MLE and both packages are maximizing the same convex
    # logistic likelihood on the same rows. That is what lets this testset assert
    # at 1e-6 instead of hand-waving about "close enough": a failure means the
    # auxiliary networks (Y+/Y-), the free-dyad selection, or a change statistic
    # is wrong.
    # ------------------------------------------------------------------
    @testset "Golden fixture: statnet tergm CMPLE on a simulated panel" begin
        g = load_golden(joinpath(@__DIR__, "fixtures", "panel_stergm.toml"))
        @test g.provenance["tergm_version"] == "4.2.2"

        # Rebuild R's exact panel from the frozen edge lists. Nothing is
        # re-simulated in Julia — the two packages see identical networks.
        n = Int(g.values["n_actors"])
        grp = String.(g.values["grp"])
        nets = Network{Int}[]
        for w in 1:Int(g.values["n_waves"])
            net = Network{Int}(; n=n, directed=true)
            for v in 1:n
                set_vertex_attribute!(net, :grp, v, grp[v])
            end
            s = Int.(g.values["wave_src"][w])
            d = Int.(g.values["wave_dst"][w])
            for k in eachindex(s)
                add_edge!(net, s[k], d[k])
            end
            push!(nets, net)
        end
        @test [ne(x) for x in nets] == Int.(g.values["edge_counts"])

        terms() = [Edges(), NodeMatch(:grp)]
        fit = stergm(nets, terms(), terms())
        @test fit.converged
        coefs = vcat(fit.formation_coef, fit.persistence_coef)
        ses = vcat(fit.formation_se, fit.persistence_se)

        # (0) The accessors BY NAME against tergm's own Form()/Persist()
        # blocks: `persistence_coef` is tergm ≥ 4's `Persist()` sign (positive
        # = ties last longer), never `Diss()`'s negation — pinned at the same
        # tolerances as the stacked vectors. And `coeftable` labels the rows
        # exactly as `summary(tergm)` does.
        for (key, val) in (("formation_coefficients", formation_coef(fit)),
                           ("formation_std_errors", formation_se(fit)),
                           ("persistence_coefficients", persistence_coef(fit)),
                           ("persistence_std_errors", persistence_se(fit)))
            @test check_golden(g, key, val) || error(golden_report(g, key, val))
        end
        @test all(persistence_coef(fit) .> 0)     # Persist() sign, not Diss()
        @test coeftable(fit).names == String.(g.values["term_names"])
        @test coeftable(fit).estimates == coefs
        @test check_statsapi(fit; strict=true) !== nothing

        # (1) Against tergm AS SHIPPED. Its CMPLE runs R's glm at the default
        # epsilon = 1e-8 and stops ~1.5e-8 (coefficients) / ~6e-6 (standard
        # errors) short of the exact optimum, so `std_errors` is compared at
        # 1e-4 — a tolerance that measures R's slack, not Julia's error.
        @test check_golden(g, "coefficients", coefs) ||
              error(golden_report(g, "coefficients", coefs))
        @test check_golden(g, "std_errors", ses) ||
              error(golden_report(g, "std_errors", ses))

        # (2) Against the SAME estimator taken to convergence (R glm at
        # epsilon = 1e-14). This is the assertion with teeth: TERGM.jl lands on
        # the exact conditional MLE, reproducing R to ~1e-10 in both the
        # coefficients and the standard errors.
        @test check_golden(g, "exact_coefficients", coefs) ||
              error(golden_report(g, "exact_coefficients", coefs))
        @test check_golden(g, "exact_std_errors", ses) ||
              error(golden_report(g, "exact_std_errors", ses))
        # ...and it is CLOSER to the exact optimum than tergm itself is.
        exact_c = Float64.(g.values["exact_coefficients"])
        exact_se = Float64.(g.values["exact_std_errors"])
        @test maximum(abs.(coefs .- exact_c)) <
              Float64(g.values["cmple_vs_glm_max_abs_diff"])
        @test maximum(abs.(ses .- exact_se)) <
              Float64(g.values["cmple_se_vs_glm_max_abs_diff"])

        # (3) Block-bootstrap standard errors. Both sides resample the 7
        # transitions with replacement (btergm's scheme), but they draw
        # DIFFERENT resamples, so this can only be compared at the resolution of
        # the bootstrap's own noise — which the fixture measures by rerunning the
        # R bootstrap under five further seeds (`boot_se_seed_sd`). The Julia
        # side averages five seeded bootstraps for the same reason. Seeded on
        # both sides, so the comparison is deterministic.
        boot_ses = [begin
                        b = stergm(nets, terms(), terms(); se=:bootstrap,
                                   n_boot=400, rng=Random.Xoshiro(s))
                        vcat(b.formation_se, b.persistence_se)
                    end for s in (101, 202, 303, 404, 505)]
        boot_se = mean(boot_ses)
        @test check_golden(g, "boot_std_errors", boot_se) ||
              error(golden_report(g, "boot_std_errors", boot_se))
        # The gap to R is no larger than the bootstrap's own seed-to-seed noise,
        # times a small factor — the two bootstraps differ by about as much as
        # either differs from itself.
        r_sd = Float64.(g.values["boot_se_seed_sd"])
        @test all(abs.(boot_se .- Float64.(g.values["boot_std_errors"])) .< 3 .* r_sd)
    end

    # ------------------------------------------------------------------
    # Golden fixture, second panel: a statistic at the boundary of its
    # attainable range. No same-group tie ever persists, so the persistence
    # nodematch has no finite CMPLE. TERGM.jl applies R ergm's `drop`
    # semantics (coefficient -Inf, SE 0, the rest fit on the untouched rows);
    # tergm 4.2.2 warns "The MPLE does not exist!" and returns a point on the
    # asymptote (-19.6, SE 1621 — recorded in the fixture, NOT asserted). The
    # finite coefficients are the same estimator on the same rows in both
    # packages, and are pinned at the main panel's tolerances.
    # ------------------------------------------------------------------
    @testset "Golden fixture: boundary statistic follows R ergm's drop" begin
        g = load_golden(joinpath(@__DIR__, "fixtures", "panel_stergm.toml"))
        n = Int(g.values["boundary_n_actors"])
        grp = String.(g.values["boundary_grp"])
        nets = Network{Int}[]
        for w in 1:Int(g.values["boundary_n_waves"])
            net = Network{Int}(; n=n, directed=true)
            for v in 1:n
                set_vertex_attribute!(net, :grp, v, grp[v])
            end
            s = Int.(g.values["boundary_wave_src"][w])
            d = Int.(g.values["boundary_wave_dst"][w])
            for k in eachindex(s)
                add_edge!(net, s[k], d[k])
            end
            push!(nets, net)
        end
        @test [ne(x) for x in nets] == Int.(g.values["boundary_edge_counts"])
        # The boundary really is there: no same-group prior tie persists
        @test all(!has_edge(nets[t], src(e), dst(e))
                  for t in 2:length(nets) for e in edges(nets[t-1])
                  if grp[src(e)] == grp[dst(e)])

        terms() = [Edges(), NodeMatch(:grp)]
        # R ergm's sentence, under the cmple context, naming the statistic
        fit = @test_logs (:warn, r"cmple: observed statistic\(s\) Persist~nodematch.grp are at their smallest attainable values") match_mode=:any stergm(
            nets, terms(), terms())
        θ = coef(fit)
        ses = stderror(fit)
        fixed = Int(g.values["boundary_fixed_index"])
        fin = Int.(g.values["boundary_finite_index"])
        @test θ[fixed] == (Int(g.values["boundary_fixed_sign"]) < 0 ? -Inf : Inf)
        @test ses[fixed] == 0.0
        @test fit.converged                       # the reduced fit converged
        @test all(isfinite, θ[fin])
        # tergm's row labels, and the fixed row carries p = 0 in the table as
        # in the printed block; the Wald interval of a fixed coefficient is
        # degenerate rather than a finite number invented from SE 0
        tbl = coeftable(fit)
        @test tbl.names == String.(g.values["boundary_term_names"])
        @test tbl.p_values[fixed] == 0.0
        @test tbl.estimates[fixed] == θ[fixed]
        @test all(confint(fit)[fixed, :] .== θ[fixed])
        @test all(isfinite, confint(fit)[fin, :])
        @test check_statsapi(fit; strict=true) !== nothing

        # The finite coefficients: against tergm as shipped and against the
        # exact limit of the pseudo-likelihood, at the main panel's tolerances
        for key in ("boundary_finite_coefficients", "boundary_exact_finite_coefficients")
            @test check_golden(g, key, θ[fin]) || error(golden_report(g, key, θ[fin]))
        end
        for key in ("boundary_finite_std_errors", "boundary_exact_finite_std_errors")
            @test check_golden(g, key, ses[fin]) || error(golden_report(g, key, ses[fin]))
        end
        # The log-likelihood: R's differs from the exact dropped limit by the
        # asymptote's tail only
        ll = loglikelihood(fit)
        @test check_golden(g, "boundary_loglik", ll) ||
              error(golden_report(g, "boundary_loglik", ll))
        @test ll ≈ Float64(g.values["boundary_dropped_loglik"]) atol = 1e-8
        tail = Float64(g.values["boundary_asymptote_tail"])
        @test abs(ll - Float64(g.values["boundary_loglik"])) < 2 * abs(tail)
        # R's warning is the fixture's record of what tergm does instead
        @test "The MPLE does not exist!" in String.(g.values["boundary_warnings"])
        @test sign(Float64(g.values["boundary_coefficients"][fixed])) ==
              Int(g.values["boundary_fixed_sign"])

        # StatsAPI after a drop: nobs is every free dyad (R's logLik nobs),
        # dof counts the finite coefficients (R ergm's drop convention; tergm
        # itself reports 4 because it did not drop) and BIC's sample size is
        # the untouched free dyads
        @test nobs(fit) == Int(g.values["boundary_loglik_nobs"])
        @test dof(fit) == length(fin) == Int(g.values["boundary_loglik_df"]) - 1
        @test fit.n_kept == Int(g.values["boundary_n_kept"])
        @test bic(fit) ≈ -2 * ll + dof(fit) * log(fit.n_kept)
        @test aic(fit) ≈ -2 * ll + 2 * dof(fit)

        # Said everywhere: approximations, show, p-value 0, is_exact false
        @test !is_exact(fit)
        @test any(occursin("Persist~nodematch.grp fixed at -Inf", a)
                  for a in approximations(fit))
        out = sprint(show, fit)
        @test occursin("fixed at -Inf", out)
        @test occursin("drop=TRUE", out)
        @test TERGM._pvalues(θ, θ ./ ses)[fixed] == 0.0

        # The block bootstrap is refused: the statistic is at its boundary on
        # every resample of the transitions
        msg = try
            stergm(nets, terms(), terms(); se=:bootstrap, n_boot=10)
            "no error"
        catch e
            e isa ArgumentError ? e.msg : rethrow()
        end
        @test occursin("se=:bootstrap is not available", msg)
        @test occursin("Persist~nodematch.grp", msg)

        # Dropping the term gives the same finite coefficients on the same
        # rows — but the persistence edges coefficient is then fit on ALL
        # prior ties, not the untouched ones, so it differs: the drop is not
        # "the model without the term"
        alt = stergm(nets, terms(), [Edges()])
        @test formation_coef(alt) ≈ θ[1:2] atol = 1e-10
        @test !(persistence_coef(alt)[1] ≈ θ[3])
    end

    # Rebuild a frozen panel of the golden fixture from its edge lists (the
    # main panel has no key prefix; the others are `boundary_`, `undirected_`)
    function golden_panel(g, prefix::String; directed::Bool)
        n = Int(g.values[prefix * "n_actors"])
        grp = String.(g.values[prefix * "grp"])
        nets = directed ? Network{Int,true}[] : Network{Int,false}[]
        for w in 1:Int(g.values[prefix * "n_waves"])
            net = Network{Int}(; n=n, directed=directed)
            for v in 1:n
                set_vertex_attribute!(net, :grp, v, grp[v])
            end
            s = Int.(g.values[prefix * "wave_src"][w])
            d = Int.(g.values[prefix * "wave_dst"][w])
            for k in eachindex(s)
                add_edge!(net, s[k], d[k])
            end
            push!(nets, net)
        end
        @test [ne(x) for x in nets] == Int.(g.values[prefix * "edge_counts"])
        return nets
    end

    # Assert one frozen tergm fit (`<prefix>_coefficients`, `_std_errors`,
    # `_loglik`, `_loglik_nobs`, `_term_names`, `_warnings`) against a
    # TERGM.jl fit of the same panel and formula
    function check_golden_fit(g, prefix::String, fit)
        @test fit.converged
        @test isempty(String.(g.values[prefix * "_warnings"]))   # R converged silently too
        @test coeftable(fit).names == String.(g.values[prefix * "_term_names"])
        for (key, val) in ((prefix * "_coefficients", coef(fit)),
                           (prefix * "_std_errors", stderror(fit)),
                           (prefix * "_loglik", loglikelihood(fit)))
            @test check_golden(g, key, val) || error(golden_report(g, key, val))
        end
        @test nobs(fit) == Int(g.values[prefix * "_loglik_nobs"])
        @test dof(fit) == Int(g.values[prefix * "_loglik_df"])
        @test check_statsapi(fit; strict=true) !== nothing
        return nothing
    end

    # ------------------------------------------------------------------
    # Golden fixture, dyad-dependent formulas on the main panel (panel
    # 2026-09 round 2). tergm's CMPLE for `mutual`/`triangle` is the same
    # estimator on the same rows: what this pins is the evaluation of a
    # dyad-dependent change statistic ON Y⁺/Y⁻ (not on Y_t), which nothing
    # in the repo asserted against R before. The triangle fit is also the
    # regression test for the convergence verdict: `Networks.newton_fit`
    # returned it `converged == false` one Newton step short of the optimum
    # (log-likelihood flat at the noise floor, gradient norm 2e-4 in the
    # triangle column's units) where R's glm calls the same point converged.
    # ------------------------------------------------------------------
    @testset "Golden fixture: dyad-dependent CMPLE (mutual, triangle) matches tergm" begin
        g = load_golden(joinpath(@__DIR__, "fixtures", "panel_stergm.toml"))
        nets = golden_panel(g, ""; directed=true)

        mut = stergm(nets, [Edges(), Mutual()], [Edges(), Mutual()])
        check_golden_fit(g, "dyaddep_mutual", mut)
        @test has_dyad_dependent(mut.model)
        @test !is_exact(mut)                       # a CMPLE, and it says so

        tri = @test_logs stergm(nets, [Edges(), Triangle()], [Edges()])   # no warning at all
        check_golden_fit(g, "dyaddep_triangle", tri)
        # ... and it is the exact optimum, not the point one step short of it:
        # the gradient of the pooled formation pseudo-likelihood vanishes there
        Xf, yf, _, _ = TERGM._cmple_blocks(tri.model)
        X = reduce(vcat, Xf); y = reduce(vcat, yf)
        _, grad, _ = Networks.logistic_derivatives(X, ones(length(y)), Float64.(y))(formation_coef(tri))
        @test maximum(abs, grad) < 1e-8
        # (`Networks.newton_fit` itself declares convergence on the Newton
        # decrement — the polish TERGM used to carry as `_polish_newton` — and
        # a fit stopped by the Newton cap stays unconverged, see "Unconverged
        # and degenerate CMPLE fits are loud")
        @test Networks.newton_fit(Networks.logistic_derivatives(X, ones(length(y)), Float64.(y)),
                                  zeros(2); maxiter=1).converged == false
    end

    # ------------------------------------------------------------------
    # Golden fixture, undirected panel: the n(n−1)/2 unordered free-dyad
    # enumeration (nobs 364 = R's logLik nobs), `nodematch`/`triangle`/
    # `kstar(2)` on undirected Y⁺/Y⁻ and R's undirected labels.
    # ------------------------------------------------------------------
    @testset "Golden fixture: undirected panel matches tergm" begin
        g = load_golden(joinpath(@__DIR__, "fixtures", "panel_stergm.toml"))
        nets = golden_panel(g, "undirected_"; directed=false)
        @test all(!is_directed(x) for x in nets)

        fit = stergm(nets, [Edges(), NodeMatch(:grp)], [Edges(), Triangle()])
        @test fit isa STERGMResult{Int,false}
        check_golden_fit(g, "undirected", fit)
        @test nobs(fit) == (length(nets) - 1) * 14 * 13 ÷ 2

        ks = stergm(nets, [Edges(), Kstar(2)], [Edges()])
        check_golden_fit(g, "undirected_kstar", ks)
        @test coeftable(ks).names[2] == "Form(1)~kstar2"
    end

    # ------------------------------------------------------------------
    # Golden fixture, expanding terms (panel 2026-09 round 3): statnet
    # resolves nodefactor / nodemix / degree(a:b) into one column per level
    # / cell / degree at model construction, and tergm's operators inherit
    # it. `check_golden_fit` asserts the term names verbatim first, so the
    # pooled single column the package used to fit fails on the labels
    # before any number is compared.
    # ------------------------------------------------------------------
    @testset "Golden fixture: expanding terms (nodefactor, nodemix, degree) match tergm" begin
        g = load_golden(joinpath(@__DIR__, "fixtures", "panel_stergm.toml"))
        @test occursin("nodefactor", g.provenance["expand_models"])

        un = golden_panel(g, "undirected_"; directed=false)
        grp3 = String.(g.values["undirected_grp3"])
        @test sort(unique(grp3)) == ["a", "b", "c"]
        for w in un, v in 1:nv(w)
            set_vertex_attribute!(w, :grp3, v, grp3[v])
        end
        fit = stergm(un, [Edges(), NodeFactor(:grp3)], [Edges(), NodeMix(:grp3)])
        @test length(coef(fit)) == 9                       # 1 + 2 + 1 + 5
        check_golden_fit(g, "expand_undirected", fit)
        deg = stergm(un, [Edges(), Degree(1:2)], [Edges(), Degree(0:1)])
        check_golden_fit(g, "expand_undirected_degree", deg)

        nets = golden_panel(g, ""; directed=true)
        dfit = stergm(nets, [Edges(), NodeFactor(:grp)], [Edges(), NodeMix(:grp)])
        @test length(coef(dfit)) == 6                      # 1 + 1 + 1 + 3
        check_golden_fit(g, "expand_directed", dfit)
        ddeg = stergm(nets, [Edges(), NodeMix(:grp)], [Edges(), ODegree(0:1)])
        check_golden_fit(g, "expand_directed_degree", ddeg)
    end

    # ------------------------------------------------------------------
    # The memory terms are refused in a CMPLE formula (panel 2026-09 round 2):
    # on the free dyads of a separable model their change statistic is
    # ±edges or identically zero, so the design is rank-deficient beside
    # `Edges()` and `edges` under another name without it. They remain
    # descriptives (`compute`, and `gof`'s formed/persisted counts).
    # ------------------------------------------------------------------
    @testset "Memory terms are refused in a separable formula" begin
        nets = fixture_panels()
        msg_of(f) = try; f(); "no error"; catch e; e isa ArgumentError ? e.msg : rethrow(); end
        for (term, label) in ((EdgeStability(), "edge.stability"),
                              (PersistentEdge(), "persistent.edges"),
                              (NewEdge(), "new.edges"))
            for (side, f, d) in (("formation", [Edges(), term], [Edges()]),
                                 ("dissolution", [Edges()], [Edges(), term]))
                m = msg_of(() -> STERGMModel(STERGM(f, d), nets))
                @test startswith(m, "$side model: term '$label' is constant on the free dyads")
                @test occursin("separable", m)
                @test occursin("btergm", m)
                @test occursin("compute(term, curr, prev)", m)
                # ... through the entry point too, and alone (without Edges)
                @test_throws ArgumentError stergm(nets, f, d)
                @test_throws ArgumentError stergm(nets, side == "formation" ? [term] : [Edges()],
                                                  side == "formation" ? [Edges()] : [term])
            end
            # the term itself keeps working as a descriptive
            @test compute(term, nets[2], nets[1]) isa Float64
            @test TERGM._separable_constant(term, :formation) isa String
            @test TERGM._separable_constant(term, :dissolution) isa String
        end
        # the mechanism, stated: on the free dyads the change statistic IS
        # constant — −1 / +1 / 0 exactly as the messages say
        prev, curr = nets[1], nets[2]
        yplus = formation_network(prev, curr); yminus = dissolution_network(prev, curr)
        free_f = [(i, j) for i in 1:5 for j in 1:5 if i != j && !has_edge(prev, i, j)]
        free_d = [(i, j) for i in 1:5 for j in 1:5 if i != j && has_edge(prev, i, j)]
        @test all(change_stat(EdgeStability(), yplus, i, j, prev) == -1.0 for (i, j) in free_f)
        @test all(change_stat(EdgeStability(), yminus, i, j, prev) == 1.0 for (i, j) in free_d)
        @test all(change_stat(PersistentEdge(), yplus, i, j, prev) == 0.0 for (i, j) in free_f)
        @test all(change_stat(PersistentEdge(), yminus, i, j, prev) == 1.0 for (i, j) in free_d)
        @test all(change_stat(NewEdge(), yplus, i, j, prev) == 1.0 for (i, j) in free_f)
        @test all(change_stat(NewEdge(), yminus, i, j, prev) == 0.0 for (i, j) in free_d)
        # informative temporal terms are untouched
        @test TERGM._separable_constant(Delrecip(), :formation) === nothing
        @test TERGM._separable_constant(Edges(), :dissolution) === nothing
        @test STERGMModel(STERGM([Edges(), Delrecip()], [Edges()]), nets) isa STERGMModel
    end

    # ------------------------------------------------------------------
    # Two-mode panels and self-loops are refused at construction (panel
    # 2026-09 round 2): the one-mode CMPLE, the free-dyad lists and the
    # sampler range over the off-diagonal one-mode dyads only — ERGM.jl
    # refuses the same networks.
    # ------------------------------------------------------------------
    @testset "Two-mode and self-loop panels are refused" begin
        msg_of(f) = try; f(); "no error"; catch e; e isa ArgumentError ? e.msg : rethrow(); end
        bip = [network(6; bipartite=3) for _ in 1:3]
        for (t, (i, j)) in enumerate(((1, 4), (2, 5), (3, 6)))
            add_edge!(bip[t], 1, 4); add_edge!(bip[t], i, j)
        end
        @test all(is_two_mode, bip)
        m = msg_of(() -> STERGMModel(STERGM([Edges()], [Edges()]), bip))
        @test occursin("TERGM (panel 1)", m)
        @test occursin("two-mode (bipartite)", m)
        @test occursin("within-mode", m)
        @test_throws ArgumentError stergm(bip, [Edges()], [Edges()])
        # a two-mode panel later in the sequence is named by its index
        mixed = [network(6), network(6), network(6; bipartite=3)]
        @test occursin("TERGM (panel 3)", msg_of(() -> STERGMModel(STERGM([Edges()], [Edges()]), mixed)))

        lp = [network(4; directed=true, loops=true) for _ in 1:2]
        for w in lp
            add_edge!(w, 1, 1); add_edge!(w, 1, 2)
        end
        add_edge!(lp[2], 3, 3)
        m2 = msg_of(() -> STERGMModel(STERGM([Edges()], [Edges()]), lp))
        @test occursin("TERGM (panel 1)", m2)
        @test occursin("1 self-loop (at vertex 1)", m2)
        @test occursin("rem_edge!(net, v, v)", m2)
        rem_edge!(lp[1], 1, 1)
        m3 = msg_of(() -> STERGMModel(STERGM([Edges()], [Edges()]), lp))
        @test occursin("TERGM (panel 2)", m3)
        @test occursin("2 self-loops (at vertex 1, 3)", m3)
        # loops *allowed* but absent is fine
        rem_edge!(lp[2], 1, 1); rem_edge!(lp[2], 3, 3)
        @test nobs(stergm(lp, [Edges()], [Edges()])) == 12
    end

    # ------------------------------------------------------------------
    # Model-specification ergonomics (panel 2026-09 round 2; ERGM.jl item
    # 31): programmatic term lists, bare terms, swapped arguments and a
    # single network are said in words, never a MethodError.
    # ------------------------------------------------------------------
    @testset "stergm accepts what fit_ergm accepts and names the common slips" begin
        nets = random_panels(Random.Xoshiro(3); n=8, T=4)
        ref = stergm(nets, [Edges(), Delrecip()], [Edges()])
        terms = []
        push!(terms, Edges()); push!(terms, Delrecip())
        @test terms isa Vector{Any}
        @test coef(stergm(nets, terms, Any[Edges()])) == coef(ref)
        @test coef(stergm(nets, terms, Edges())) == coef(ref)          # a bare term
        @test coef(fit_stergm(nets, [Edges(), [Delrecip()]], [Edges()])) == coef(ref)   # nested
        @test coef(stergm(nets, Edges(), Edges())) == coef(stergm(nets, [Edges()], [Edges()]))
        f = STERGM(terms, Edges())
        @test f.formation isa Vector{AbstractERGMTerm}
        @test f.formation == AbstractERGMTerm[Edges(), Delrecip()]

        err(f) = try; f(); nothing; catch e; e; end
        for bad in (() -> stergm([Edges()], [Edges()], nets),
                    () -> stergm(Edges(), Edges(), nets),
                    () -> stergm(Any[Edges()], [Edges()], nets),
                    () -> stergm(nets, [Edges()], nets),
                    () -> fit_stergm([Edges()], [Edges()], nets),
                    # the panel in the second position ("model first"), alone
                    # or beside a panel in another slot (round 3)
                    () -> stergm([Edges()], nets, [Edges()]),
                    () -> stergm(Any[Edges()], nets, Edges()),
                    () -> stergm(nets, nets, [Edges()]),
                    () -> stergm([Edges()], nets, nets),
                    () -> stergm(nets, nets, nets))
            e = err(bad)
            @test e isa ArgumentError
            @test occursin("arguments are swapped", e.msg)
            @test occursin("stergm(networks, formation, dissolution)", e.msg)
            @test occursin("Form(~", e.msg)
        end
        e = err(() -> stergm(nets[1], [Edges()], [Edges()]))
        @test e isa ArgumentError
        @test occursin("needs a panel", e.msg)
        @test occursin("fit_ergm", e.msg)
        # non-terms are named with their position and side, and a term TYPE
        # gets ERGM's "did you mean" hint
        e = err(() -> stergm(nets, [Edges(), 3], [Edges()]))
        @test e isa ArgumentError
        @test startswith(e.msg, "formation model: element 2")
        @test occursin("Int64", e.msg)
        e = err(() -> stergm(nets, [Edges()], [Edges]))
        @test e isa ArgumentError
        @test startswith(e.msg, "dissolution model:")
        @test occursin("did you mean `Edges()`", e.msg)
        e = err(() -> STERGM([Edges()], AbstractERGMTerm[]))
        @test e isa ArgumentError
        @test occursin("dissolution model must contain at least one term", e.msg)
        e = err(() -> STERGM(:edges, [Edges()]))
        @test e isa ArgumentError
        @test occursin("formation model must be a term or a vector of terms", e.msg)
        @test isempty(Test.detect_ambiguities(TERGM))
    end

    @testset "STERGMModel and STERGM print compactly; the result header has AIC/BIC" begin
        nets = fixture_panels()
        for w in nets
            set_vertex_attribute!(w, :grp, Dict(v => (v <= 2 ? "a" : "b") for v in 1:5))
        end
        f = STERGM([Edges(), Mutual(), NodeMatch(:grp)], [Edges(), GWESP(0.5)])
        @test sprint(show, f) ==
              "STERGM(formation: edges + mutual + nodematch.grp; persistence: edges + gwesp.fixed.0.5)"
        m = STERGMModel(f, nets)
        out = sprint(show, m)
        @test out == "STERGMModel{Int64,true}: 2 panels of 5 vertices (directed), 20 free " *
                     "dyads; formation: edges + mutual + nodematch.grp; persistence: " *
                     "edges + gwesp.OTP.fixed.0.5"
        @test !occursin("Network{", out)            # never the panels themselves
        @test length(out) < 200
        un = [network(4; directed=false), network(4; directed=false)]
        add_edge!(un[1], 1, 2); add_edge!(un[2], 1, 2)
        @test sprint(show, STERGMModel(STERGM([Edges()], [Edges()]), un)) ==
              "STERGMModel{Int64,false}: 2 panels of 4 vertices (undirected), 6 free " *
              "dyads; formation: edges; persistence: edges"
        # `result.model` is the first thing inspected after a fit
        r = stergm(nets, [Edges()], [Edges()])
        @test sprint(show, r.model) == sprint(show, STERGMModel(STERGM([Edges()], [Edges()]), nets))
        ro = sprint(show, r)
        @test occursin("AIC: $(round(aic(r), digits=2)), BIC: $(round(bic(r), digits=2))", ro)
        @test findfirst("AIC:", ro).start < findfirst("Standard errors:", ro).start
    end

    # ------------------------------------------------------------------
    # The CMPLE design build allocates the blocks it returns and nothing per
    # row (panel 2026-09 round 2): the row fill is a typed barrier over the
    # term tuples, not a `Vector{Float64}` per free dyad through an abstract
    # term vector followed by `hcat`.
    # ------------------------------------------------------------------
    @testset "CMPLE design build allocates O(blocks), not O(rows · p)" begin
        # The model carries an attribute term (round 3): the raw `NodeMatch`
        # used to cost 160 B per row through the attribute Dict — 548 KB on
        # the README's 4 900-row s50 design, which the structural-terms-only
        # pin never saw. The design build now snapshots it per transition
        # (`_materialized_tuple`: one O(n) codes vector and a small Dict per
        # attribute term per side per transition) and fills rows at 0 B.
        function build(n_actors, n_waves)
            rng = Random.Xoshiro(6)
            nets = Network{Int}[]
            for _ in 1:n_waves
                net = network(n_actors; directed=true)
                set_vertex_attribute!(net, :grp,
                                      Dict(v => (isodd(v) ? "a" : "b") for v in 1:n_actors))
                for i in 1:n_actors, j in 1:n_actors
                    i != j && rand(rng) < 0.15 && add_edge!(net, i, j)
                end
                push!(nets, net)
            end
            return STERGMModel(STERGM([Edges(), Mutual(), NodeMatch(:grp)], [Edges()]), nets)
        end
        function overhead(model)
            TERGM._cmple_blocks(model)                       # warm up
            total = @allocated TERGM._cmple_blocks(model)
            Xf, yf, Xd, yd = TERGM._cmple_blocks(model)
            arrays = sum(sizeof, Xf) + sum(sizeof, yf) + sum(sizeof, Xd) + sum(sizeof, yd)
            # the two auxiliary networks of every transition are copies of
            # Y_{t−1}: the data structure's cost, not the design build's
            copies = sum(@allocated((formation_network(model.networks[t-1], model.networks[t]),
                                     dissolution_network(model.networks[t-1], model.networks[t])))
                         for t in 2:length(model.networks))
            rows = sum(size.(Xf, 1)) + sum(size.(Xd, 1))
            return rows, total - arrays - copies
        end
        n_trans = 7
        rows_small, over_small = overhead(build(25, n_trans + 1))
        rows_big, over_big = overhead(build(50, n_trans + 1))
        @test rows_big > 4 * rows_small
        # 4x the rows, the same per-transition overhead: the only part that
        # may grow with the actors is the attribute snapshot — one Int per
        # vertex per attribute term (here: one term, 25 more vertices) — never
        # anything per row (which would add ≥ 8 × 3 × 4 000 rows ≈ 96 KB here)
        @test over_big <= over_small + n_trans * (8 * 25) + 1024
        @test over_big < n_trans * 4 * 1024
        # ... because the fill itself allocates nothing, attribute term included
        model = build(25, 8)
        prev, curr = model.networks[1], model.networks[2]
        yplus, yminus = formation_network(prev, curr), dissolution_network(prev, curr)
        nd = ne(prev); nf = 25 * 24 - nd
        Xf = Matrix{Float64}(undef, nf, 3); yf = Vector{Bool}(undef, nf)
        Xd = Matrix{Float64}(undef, nd, 1); yd = Vector{Bool}(undef, nd)
        fterms = TERGM._materialized_tuple(model.formula.formation, yplus)
        dterms = TERGM._materialized_tuple(model.formula.dissolution, yminus)
        TERGM._fill_blocks!(Xf, yf, Xd, yd, fterms, dterms, prev, curr, yplus, yminus, 25, Val(true))
        @test (@allocated TERGM._fill_blocks!(Xf, yf, Xd, yd, fterms, dterms, prev, curr,
                                              yplus, yminus, 25, Val(true))) == 0
        # and the rows are the ones the old row-by-row build produced (the
        # materialized NodeMatch column equals the raw term's change statistic)
        expect = [Float64(!has_edge(prev, i, j) ? 1.0 : 0.0) for i in 1:25 for j in 1:25 if i != j && !has_edge(prev, i, j)]
        @test Xf[:, 1] == expect
        @test Xf[:, 2] == [change_stat(Mutual(), yplus, i, j) for i in 1:25 for j in 1:25 if i != j && !has_edge(prev, i, j)]
        @test Xf[:, 3] == [change_stat(NodeMatch(:grp), yplus, i, j) for i in 1:25 for j in 1:25 if i != j && !has_edge(prev, i, j)]
        @test yf == [has_edge(curr, i, j) for i in 1:25 for j in 1:25 if i != j && !has_edge(prev, i, j)]
        @test yd == [has_edge(curr, i, j) for i in 1:25 for j in 1:25 if i != j && has_edge(prev, i, j)]
    end
end
