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

1. **Temporal Term Types** -- `TemporalTerm <: AbstractERGMTerm` subtypes: `EdgeStability`, `Delrecip`, `PersistentEdge`, `NewEdge`. Each implements `compute(term, net, prev_net)` and `change_stat(term, net, i, j, prev_net)` (ADD-DIRECTION, state-independent; calling without prev_net raises an informative error). Standard ERGM.jl terms mix freely (dispatch via `_tchange`/`_tcompute`). Tie-duration descriptives: `edge_ages`, `mean_edge_age`. The old `FormationTerm`/`DissolutionTerm` wrappers and Memory/EdgeAge/TimeLag stubs were removed — the formation/dissolution split is expressed by which term list a term is in.
2. **STERGM Model Types** -- `STERGM` (formula holding formation/dissolution terms), `STERGMModel` (formula + observed network sequence; construction validates the formulas — attribute-based terms must reference vertex attributes present on every panel, and direction-incompatible terms (`Mutual`, `Delrecip`) are rejected on undirected panels, mirroring `ERGM.ERGMModel`), `STERGMResult` (fitted coefficients, SEs, convergence info, and `se_type` — `:hessian` or `:bootstrap`; `show()` prints a statnet-style caveat for CMPLE fits of dyad-dependent formulas). The temporal terms extend `ERGM.is_dyad_dependent` with `false` (they condition only on the exogenous previous network).
3. **Estimation** -- `stergm()` defaults to `cmple()`: pooled logistic pseudo-likelihood over the FREE DYADS of the auxiliary networks (formation rows = prior non-edges with stats on Y⁺ = Y_{t−1} ∪ Y_t; dissolution rows = prior edges with stats on Y⁻ = Y_{t−1} ∩ Y_t; `formation_network`/`dissolution_network` expose the construction), maximized with the shared `ERGM.newton_fit` optimizer. Design rows are built per transition (`_cmple_blocks`) so `se=:bootstrap` can run a per-transition block bootstrap (btergm-style: resample transitions with replacement, refit, empirical covariance of the stacked coefficients; needs ≥ 2 transitions, seeded by `rng`; point estimates unchanged). `cmle` warns and falls back to CMPLE (exact only for dyad-independent terms; MCMC CMLE not implemented). `egmme` raises an error instead of returning placeholders.
4. **Simulation & Diagnostics** -- `simulate_stergm()` samples Y⁺ (Metropolis over prior non-edges) and Y⁻ (over prior edges) and combines Y_t = (Y⁺ \\ Y_{t−1}) ∪ Y⁻; `simulate_network_sequence()` iterates; `stergm_gof()` compares formed/persisted tie counts per transition with MC p-values.

The main entry point is `stergm(networks, formation_terms, dissolution_terms; method=:cmple)`.

## Key Dependencies

- **ERGM.jl** -- provides `AbstractERGMTerm`, standard ERGM terms (`Edges`, `Triangle`, `Mutual`), and `change_stat`/`compute` interface
- **Network.jl** -- network data structure (`Network{T}`), `has_edge`, `add_edge!`, `rem_edge!`, `nv`, `ne`, `edges`, `is_directed`
- **Distributions.jl** -- distribution utilities

Requires Julia 1.12+.

## Conventions

- All temporal terms subtype `TemporalTerm <: AbstractERGMTerm` and must implement `name()` and `compute(term, net, prev_net)`.
- Change statistics follow the ERGM.jl ADD-DIRECTION convention: `change_stat(term, net, i, j)` = g(y⁺ij) − g(y⁻ij), state-independent. Temporal variants take a trailing `prev_net` argument and obey the same convention (brute-force verified in tests).
- Formation terms operate on dyads that were non-edges at t-1; dissolution terms operate on dyads that were edges at t-1.
- The dissolution model estimates **persistence** (not dissolution) coefficients -- positive values mean higher persistence.
- The package uses `Float64` for all coefficients and statistics.
- Exports are declared at the top of the module; no re-exports from dependencies.
