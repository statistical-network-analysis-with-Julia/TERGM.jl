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
3. **Estimation** -- `stergm()` defaults to `cmple()`: pooled logistic pseudo-likelihood over the FREE DYADS of the auxiliary networks (formation rows = prior non-edges with stats on Y⁺ = Y_{t−1} ∪ Y_t; dissolution rows = prior edges with stats on Y⁻ = Y_{t−1} ∩ Y_t; `formation_network`/`dissolution_network` expose the construction), maximized with the shared `ERGM.newton_fit` optimizer. Design rows are built per transition (`_cmple_blocks`) so `se=:bootstrap` can run a per-transition block bootstrap (btergm-style: resample transitions with replacement, refit, empirical covariance of the stacked coefficients; needs ≥ 2 transitions, seeded by `rng`; point estimates unchanged). `cmle` throws an `ArgumentError` (MCMC CMLE is not implemented, and CMPLE is only the same estimator for dyad-independent terms, so it does not silently stand in). `egmme` raises an error instead of returning placeholders, and is **not exported** — an estimator that can only throw must not advertise itself in the public API; it stays reachable as `TERGM.egmme` so the error is informative rather than an `UndefVarError`, and `method=:egmme` still raises.
4. **Simulation & Diagnostics** -- `simulate_stergm()` samples Y⁺ (Metropolis over prior non-edges) and Y⁻ (over prior edges) and combines Y_t = (Y⁺ \\ Y_{t−1}) ∪ Y⁻; `simulate_network_sequence()` iterates; `stergm_gof()` compares formed/persisted tie counts per transition with MC p-values.

The main entry point is `stergm(networks, formation_terms, dissolution_terms; method=:cmple)`.

## The CMPLE derivative loop is ERGM.jl's, not ours (review finding 15)

`_logistic_fit` used to carry its own copy of the logistic loop with `hess .-= (pr*(1-pr)) .* (x * x')` inside it: a fresh `p×p` outer product on every one of the design rows of every Newton evaluation — **471 KB and 0.46 ms per evaluation** on a 25-actor, 8-wave panel. It now runs on the shared **`ERGM.logistic_derivatives(X, y)`** (the same builder ERGMMulti's MPLE and ERGMRank's swap MPLE use): workspaces allocated once, `η = Xβ` by gemv and `−H = X'WX` by gemm, **192 bytes and 0.11 ms per evaluation** — 4.0x faster, with allocations independent of the number of rows. **Never paste the loop back in**; if it needs a feature, add it to `ERGM.logistic_derivatives`. The `@allocated` regression test measures the closure built from TERGM's own `_cmple_blocks` design rows, so it guards this package's path and not just ERGM's function. The summation order moves from row-wise accumulation to BLAS, so the arithmetic is not bit-identical — but the *fit* is: measured against the old loop on the same design, **max|Δθ| = 0.0**. (Newton's last step is quadratically convergent; a last-ulp difference in the gradient and Hessian does not move its fixed point.) The golden fixture is unmoved at its 1e-6 exact-CMLE tolerance, and TERGM.jl remains closer to the exact optimum than `tergm` itself.

## Golden fixtures (statnet tergm)

`test/fixtures/panel_stergm.toml`, generated by `test/fixtures/r/panel_stergm.R` and loaded with Networks.jl's `load_golden` (which throws without `[provenance]`). An 8-wave, 25-actor directed panel is frozen as edge lists, so both packages fit identical networks; the model is `Form(~edges + nodematch("grp")) + Persist(~edges + nodematch("grp"))` under `estimate="CMPLE"`.

The formulas are **dyad-independent on purpose**: the conditional pseudo-likelihood is then the conditional likelihood, CMPLE is the exact conditional MLE, and agreement can be asserted at 1e-6 instead of hand-waved. A dyad-dependent formula would have bought a prettier model and destroyed the ability to say anything exact.

**What it found — and it is R's slack, not ours.** `tergm`'s CMPLE runs R's `glm` at the default `epsilon = 1e-8` and stops **1.5e-8** (coefficients) / **6.1e-6** (standard errors) short of the exact optimum. TERGM.jl's Newton–Raphson lands *on* the optimum, reproducing a tightly-converged (`epsilon = 1e-14`) R `glm` to ~1e-10. So the fixture freezes **both**: `coefficients`/`std_errors` are tergm as shipped (SEs compared at 1e-4 — a tolerance that measures R), and `exact_coefficients`/`exact_std_errors` are the same estimator taken to convergence (compared at 1e-6 — the tolerance that measures Julia). The testset also asserts TERGM.jl is *closer to the exact optimum than tergm itself is*. Do not "fix" a red test here by loosening the exact tolerance.

Block-bootstrap SEs are compared at the resolution the bootstrap actually has: the R script reruns btergm's transition-resampling scheme under five further seeds and freezes the seed-to-seed sd of every bootstrap SE, and the tolerance is a multiple of that. Two bootstraps drawing different resamples cannot be compared any other way.

## Key Dependencies

- **ERGM.jl** -- provides `AbstractERGMTerm`, standard ERGM terms (`Edges`, `Triangle`, `Mutual`), and `change_stat`/`compute` interface
- **Networks.jl** -- network data structure (`Network{T}`), `has_edge`, `add_edge!`, `rem_edge!`, `nv`, `ne`, `edges`, `is_directed`
- **Distributions.jl** -- distribution utilities

Requires Julia 1.12+.

## Conventions

- All temporal terms subtype `TemporalTerm <: AbstractERGMTerm` and must implement `name()` and `compute(term, net, prev_net)`.
- Change statistics follow the ERGM.jl ADD-DIRECTION convention: `change_stat(term, net, i, j)` = g(y⁺ij) − g(y⁻ij), state-independent. Temporal variants take a trailing `prev_net` argument and obey the same convention (brute-force verified in tests).
- Formation terms operate on dyads that were non-edges at t-1; dissolution terms operate on dyads that were edges at t-1.
- The dissolution model estimates **persistence** (not dissolution) coefficients -- positive values mean higher persistence.
- The package uses `Float64` for all coefficients and statistics.
- Exports are declared at the top of the module; no re-exports from dependencies.
