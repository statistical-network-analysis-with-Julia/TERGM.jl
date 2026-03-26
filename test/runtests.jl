using TERGM
using Test

@testset "TERGM.jl" begin
    @testset "Module loading" begin
        @test @isdefined(TERGM)
    end

    @testset "Temporal term types" begin
        @test EdgeStability() isa EdgeStability
        @test PersistentEdge() isa PersistentEdge
        @test NewEdge() isa NewEdge
        @test Delrecip() isa Delrecip
        @test TERGM.Memory() isa TERGM.Memory
        @test TERGM.Memory(0.5) isa TERGM.Memory
        @test EdgeAge(10) isa EdgeAge
    end

    @testset "Formation and Dissolution wrappers" begin
        # FormationTerm and DissolutionTerm need an AbstractERGMTerm
        # EdgeStability is a TemporalTerm <: AbstractERGMTerm
        ft = FormationTerm(EdgeStability())
        @test ft isa FormationTerm
        dt = DissolutionTerm(EdgeStability())
        @test dt isa DissolutionTerm
    end

    @testset "TimeLag" begin
        tl = TimeLag(EdgeStability())
        @test tl isa TimeLag
        @test tl.lag == 1

        tl2 = TimeLag(NewEdge(), 2)
        @test tl2.lag == 2
    end

    @testset "STERGM construction" begin
        form_terms = [EdgeStability()]
        diss_terms = [PersistentEdge()]
        s = STERGM(form_terms, diss_terms)
        @test s isa STERGM
        @test length(s.formation_terms) == 1
        @test length(s.dissolution_terms) == 1
        @test isempty(s.constraints)
    end

    @testset "Estimation API" begin
        @test stergm === fit_stergm
        @test isdefined(TERGM, :cmle)
        @test isdefined(TERGM, :cmple)
        @test isdefined(TERGM, :egmme)
    end

    @testset "Simulation API" begin
        @test isdefined(TERGM, :simulate_stergm)
        @test isdefined(TERGM, :simulate_network_sequence)
    end

    @testset "Diagnostics" begin
        @test isdefined(TERGM, :stergm_gof)
    end
end
