# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

TERGM.jl is a Julia port of the R [tergm](https://cran.r-project.org/package=tergm) package from the StatNet suite, implementing Separable Temporal Exponential Random Graph Models (STERGMs) for statistical modeling of dynamic networks observed at discrete time points.

## Development Commands

- **Run tests:** `julia --project=. -e 'include("test/runtests.jl")'` (or `include("test/runtests.jl")` from the Julia REPL)
- **Build docs:** Uses Documenter.jl; workflow in `.github/workflows/Documentation.yml`
- **Add package for development:** `using Pkg; Pkg.develop(path=".")`

## Architecture

The entire package lives in a single file: `src/TERGM.jl`. It is organized into four sections:

1. **Temporal Term Types** -- `TemporalTerm <: AbstractERGMTerm` subtypes: `FormationTerm`, `DissolutionTerm`, `EdgeStability`, `PersistentEdge`, `NewEdge`, `Memory`, `EdgeAge`, `Delrecip`, `TimeLag`. Each implements `name()`, `compute()`, and optionally `change_stat()`.
2. **STERGM Model Types** -- `STERGM` (formula holding formation/dissolution terms), `STERGMModel` (formula + observed network sequence), `STERGMResult` (fitted coefficients, SEs, convergence info).
3. **Estimation Functions** -- `stergm()` dispatches to `cmle()`, `cmple()`, or `egmme()`. CMLE uses Newton-Raphson on logistic change-stat likelihoods. CMPLE currently delegates to CMLE. EGMME is experimental/placeholder.
4. **Simulation & Diagnostics** -- `simulate_stergm()`, `simulate_one_step()`, `simulate_network_sequence()`, and `stergm_gof()`.

The main entry point is `stergm(networks, formation_terms, dissolution_terms; method=:cmle)`.

## Key Dependencies

- **ERGM.jl** -- provides `AbstractERGMTerm`, standard ERGM terms (`Edges`, `Triangle`, `Mutual`), and `change_stat`/`compute` interface
- **Network.jl** -- network data structure (`Network{T}`), `has_edge`, `add_edge!`, `rem_edge!`, `nv`, `ne`, `edges`, `is_directed`
- **NetworkDynamic.jl** -- dynamic network support
- **Optim.jl**, **Distributions.jl**, **StatsBase.jl** -- numerical optimization and statistics

Requires Julia 1.9+.

## Conventions

- All temporal terms subtype `TemporalTerm <: AbstractERGMTerm` and must implement `name()` and `compute(term, net, prev_net)`.
- Change statistics follow the ERGM.jl convention: `change_stat(term, net, i, j)` returns the change in the statistic when toggling edge (i,j). Temporal variants accept an additional `prev_net` argument.
- Formation terms operate on dyads that were non-edges at t-1; dissolution terms operate on dyads that were edges at t-1.
- The dissolution model estimates **persistence** (not dissolution) coefficients -- positive values mean higher persistence.
- The package uses `Float64` for all coefficients and statistics.
- Exports are declared at the top of the module; no re-exports from dependencies.
