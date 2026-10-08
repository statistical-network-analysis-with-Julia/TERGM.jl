# Changelog

All notable changes to TERGM.jl are documented in this file. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
package adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - Unreleased

First public release of TERGM.jl, a Julia port of R's `tergm` (statnet):
separable temporal ERGMs (Krivitsky & Handcock 2014) for panels of
networks, fitted by conditional maximum pseudo-likelihood (CMPLE) or
Monte-Carlo conditional maximum likelihood (CMLE), with simulation and
goodness of fit. Estimates, standard errors, labels and log-likelihoods are
checked against `tergm` 4.2.2 by provenanced fixtures. The notes below are
relative to 0.1.0, which was not publicly released.

**Dependency renamed:** the foundation package is now `NetworkCore` (developed as `Networks`); write `using NetworkCore` where code said `using Networks`. Types and functions keep their names.

### Highlights

- **A default estimator**: `stergm` takes `method = :auto` — the CMLE
  (tergm's `estimate = "CMLE"`) when a term is dyad-dependent, and the
  exact CMPLE otherwise. R's `tergm()` has no default (it requires
  `estimate=`), so this is TERGM.jl's choice.
- **CMPLE on the Krivitsky–Handcock auxiliary networks**
  (`stergm(networks, formation, dissolution)`): formation statistics on
  Y⁺ = Y_{t−1} ∪ Y_t over prior non-ties, persistence statistics on
  Y⁻ = Y_{t−1} ∩ Y_t over prior ties, pooled over transitions. For a
  dyad-independent formula it is the exact conditional MLE and agrees with
  a tightly converged R `glm` to 1e-6; dyad-dependent, undirected,
  expanding-term and boundary panels are pinned against `tergm` as shipped
  (`test/fixtures/panel_stergm.toml`).
- **CMLE** (`method = :cmle`, tergm's `estimate = "CMLE"`): Monte-Carlo
  maximum likelihood on the constrained sample spaces, with ERGM.jl's
  Hummel-stepped update, R ergm's `confidence` stopping rule and sample-size
  boost, Fisher-plus-Monte-Carlo standard errors and a path-sampled
  log-likelihood. A dyad-independent side is fit exactly, without MCMC.
  Pinned in `test/fixtures/cmle_stergm.toml` against the **exact**
  conditional MLE of `edges + mutual` models (computed by enumeration) and
  against `tergm`'s CMLE of those and of the tergm tutorial's
  `edges + mutual + cyclicalties + transitiveties` model.
- **Calibrated inference under dyadic dependence.** A CMPLE fit
  (`method = :cmple`) of a dyad-dependent formula reports estimates and naive standard errors
  but no z, p or confidence interval; `se = :bootstrap` (parametric
  bootstrap: every transition simulated from the fitted model given the
  observed previous panel, CMPLE refitted) and `method = :cmle` give
  calibrated ones; `se = :hessian` is the explicit opt-in to the naive
  Wald table.
- **Every ERGM.jl term is a STERGM term**, validated and expanded exactly
  as `fit_ergm` does it (R's labels and columns for
  `nodefactor`/`nodemix`/`degree`), including `Concurrent`, `DegRange`,
  `GWNSP`, `CyclicalTies` and `TransitiveTies`; `Delrecip` (delayed
  reciprocity) is the one temporal term.
- **`simulate_stergm`, `simulate_network_sequence` and `gof`**: an exact
  constrained Metropolis sampler (checked against enumeration, and against
  `tergm`'s one-step `simulate` in `test/fixtures/simulate_stergm.toml`), and
  transition-level goodness of fit after `gof.tergm` (tie changes, model
  statistics on Y⁺/Y⁻, degree distributions) with Monte-Carlo p-values.
- **The full StatsAPI surface** on `STERGMResult` (`coef`, `stderror`,
  `vcov`, `confint`, `loglikelihood`, `nobs`, `dof`, `aic`, `bic`,
  `coeftable` with tergm's `Form(1)~…` / `Persist(1)~…` labels).
- **Reproducible and thread-count independent**: all randomness comes from
  an `rng` keyword; bootstraps, `gof` and the CMLE draw one seed per task
  up front.
- **Nothing about a bad fit is silent**: non-convergence, separation, a
  coefficient fixed at `±Inf`, excluded bootstrap replicates and withheld
  inference are warned about, printed by `show` and listed by
  `approximations(fit)`.

### Added

- `coefnames(fit)` (StatsAPI) returns the coefficient labels
  (`Form(1)~edges`, `Persist(1)~edges`, …), the same as `coeftable(fit).names`.
- `cmle(model)` / `stergm(...; method = :cmle)`: the Monte-Carlo
  conditional MLE described above. Keywords `n_samples`, `burnin`,
  `interval`, `maxiter`, `termination` (`:confidence` or `:hotelling`),
  `conv_precision`, `conv_confidence`, `init`, `bridge_rungs`, `rng`, ….
  `fit.mcmc` holds the per-side diagnostics.
- `cmple(model; se = :bootstrap, n_boot, burnin, rng)`: parametric
  bootstrap standard errors. In simulation (30 actors, 1–10 transitions)
  its 95 % intervals covered 0.91–0.98 for models with `Mutual` or
  `GWESP`.
- `STERGMResult` fields `inference_withheld`, `mcmc`, `boot_replicates`,
  `n_kept`, `se_type` (`:hessian`, `:bootstrap`, `:block_bootstrap`, `:fisher`).
- `coeftable(fit)`, `confint(fit; level)`, and the remaining StatsAPI
  verbs; the result-metadata protocol (`objective`, `is_exact`,
  `se_method`, `approximations`, `fit_metadata`).
- `gof(fit)` panels: tie changes, formation/persistence model statistics,
  and the degree distribution(s) of Y_t.
- `formation_network`, `dissolution_network`, `edge_ages`,
  `mean_edge_age`; the descriptives `EdgeStability`, `PersistentEdge`,
  `NewEdge`.
- `STERGMModel{T,D}`: a formula bound to a concretely typed panel, with
  `is_directed(model)`; one-line `show` for `STERGM` and `STERGMModel`.
- Golden fixtures with provenance and checked-in R scripts:
  `panel_stergm.toml` (CMPLE: dyad-independent, `mutual`, `triangle`,
  undirected, expanding terms, `concurrent`/`degrange`/`gwnsp`, a boundary
  panel, block-bootstrap scheme) and `cmle_stergm.toml` (CMLE).
- A precompile workload: the first `stergm` fit takes about 0.04 s instead
  of 1.6 s.

### Changed

- **ERGM.jl's building blocks come from its extension API.** The formula
  pipeline (`validate_formula`, `materialize`, `collect_terms`,
  `expand_terms`), the CMPLE fit (`mple_fit_design`), the CMLE iteration
  and bridge (`mcmle_solve`, `bridge_integrate`) and the sampler budget
  (`mcmc_defaults`) are imported from `ERGM.Extension` by name, not from
  `public` underscore names of ERGM, which no longer exist. No behaviour
  change.
- **`stergm`/`fit_stergm` default to `method = :auto`**: the CMLE
  (tergm's `estimate = "CMLE"`; R's `tergm()` has no default and requires
  `estimate=`) when either formula has a dyad-dependent term, the
  CMPLE when none has — there the CMPLE is the exact conditional MLE and no
  MCMC is run. The rule is ERGM.jl's `ERGM.resolve_method`, shared with
  `fit_ergm`. Before, every formula was fitted by CMPLE unless
  `method = :cmle` was passed. Migration: pass `method = :cmple` for the
  pseudo-likelihood fit (its `se = :bootstrap` / `:block_bootstrap` /
  `:hessian` options are CMPLE keywords: with `:auto` on a dyad-dependent
  formula they are an `ArgumentError` naming `method = :cmple`).
- **EGMME** (`method = :egmme`, `TERGM.egmme`) throws an explanatory
  `ArgumentError`, not an `ErrorException`; so do `compute(term, net)` and
  `change_stat(term, net, i, j)` of a temporal term called without the
  previous network. That `change_stat` fallback no longer types `i, j` as
  `Int`, so a user temporal term's untyped four-argument method is no
  longer reported as ambiguous.
- A separated CMPLE design is decided by the ecosystem's shared verdict
  (`NetworkCore.logistic_separation`, through ERGM.jl's design fitter) and
  warned about as "cmple: the CMPLE does not exist (separation)", naming
  the separating terms. An unconverged CMPLE (separated, or Newton cap
  exhausted) now withholds its z values, p-values and confidence intervals
  (`NaN`), as across the ecosystem.
- An error raised inside `gof`'s simulations or the CMLE's per-transition
  chains reaches the caller as the original exception, not wrapped in a
  `TaskFailedException` (the loops run on NetworkCore's `spawn_all`).
- The warning, `show` and `approximations` of a bootstrap fit that
  excluded refits without a finite estimate say that the standard errors
  are conditional on a finite refit and therefore biased downward.
- **`se = :bootstrap` is the parametric bootstrap**, as in ERGM.jl and
  the rest of the ecosystem. The per-transition block bootstrap
  (`btergm`'s scheme) is `se = :block_bootstrap`.
- **The block bootstrap (`se = :block_bootstrap`) is refused below 10
  transitions and warned about below 20**, and nothing recommends it. It resamples whole transitions, and in simulation its 95 %
  intervals covered 0.57–0.65 at 2 transitions, 0.74–0.88 at 3–5,
  0.87–0.92 at 8–10 and 0.91–0.95 at 15–30, against 0.92–0.97 for the
  inverse-Hessian errors of the same dyad-independent fits.
- **The CMPLE of a dyad-dependent formula withholds z, p and `confint`
  by default** (see Highlights). In simulation the naive standard errors
  covered 0.89–0.96 with `Mutual` but 0.70–0.82 on the formation side of
  `edges + GWESP(0.5)`.
- `method = :cmle` used to throw; it now fits the conditional MLE.
- The dissolution model is reported in tergm's **persistence**
  (`Persist()`) parameterisation only: `persistence_coef`,
  `persistence_se`. The never-released `dissolution_coef`/`dissolution_se`
  spellings and `stergm_gof` are removed (use `persistence_*` and `gof`).
- `stergm` takes the panel first, then the formation and the persistence
  terms; any term collection `fit_ergm` accepts is accepted, and swapped
  arguments or a single network are `ArgumentError`s in words.
  `simulate_stergm` / `simulate_network_sequence` take the coefficients
  positionally.
- Formula validation is ERGM.jl's, on every panel: missing or incomplete
  vertex attributes, wrong-direction terms and wrong-size `EdgeCov`
  matrices are refused, naming the side, the panel and the term.
- A statistic at the boundary of its attainable range follows R ergm's
  `drop`: the coefficient is `±Inf` with standard error 0 and the others
  are fitted on the untouched dyads. A separated design or an exhausted
  Newton cap is `converged == false` with a warning.
- Sampler `burnin` defaults to ERGM.jl's dyad-scaled rule (20 toggles per
  free dyad of the side).
- Requires Julia 1.12.

### Fixed

- The workflow testset failed outside the sibling layout (a lone checkout
  or a registry install); it now skips the CI layout step there, with a
  message.
- A statistic whose change statistics are all zero on a side's free dyads
  (`Form(~triangle)` when the union network has no two-path) made the CMPLE
  return 0 for every coefficient of that side, unconverged. Its coefficient
  is now `NaN`, with R's "not varying" warning, and the side's other
  coefficients are estimated, as tergm's CMPLE reports `NA`; `show` and
  `approximations` say so, the bootstraps refuse it, and the CMLE, which
  has no CMPLE start for that side, says why.

- Multi-level `NodeFactor`, `NodeMix` and degree ranges were fitted as one
  pooled, mislabelled column; they now expand to R's columns.
- Auxiliary networks dropped vertex attributes, so nodal terms had
  all-zero columns.
- A dyad-dependent CMPLE could be reported unconverged one Newton step
  short of the optimum (the convergence verdict is now scale-free).
- `Delrecip` returned a silent 0 on undirected networks; it now throws.

### Performance

- The CMPLE derivative evaluation allocates O(p²) instead of O(rows · p²)
  (the shared `NetworkCore.logistic_derivatives` kernel), and the design
  build allocates per transition, not per row.
- The constrained Metropolis step and the CMLE's statistic sampler
  allocate nothing per toggle (ERGM.jl's `mh_toggle!` kernel).
- Bootstrap refits, `gof` simulations and the CMLE's per-transition chains
  run on all threads.

### Differences from R tergm

- `Diss()` is not offered; `persistence_coef` has `Persist()`'s sign.
- The CMPLE drops a boundary statistic (`-Inf`, as R `ergm` does) where
  `tergm` 4.2 returns a large finite value with a warning.
- TERGM.jl's CMPLE is converged to the exact optimum; `tergm`'s stops at
  `glm`'s default tolerance (differences of 1e-8 in coefficients, 1e-5 in
  standard errors).
- The summary of a dyad-dependent CMPLE (`method = :cmple`) prints no z
  or p by default; `tergm` prints the naive Wald table (`se = :hessian`
  here).
- Multi-category term levels are resolved from the first panel; a later
  panel carrying a new level is refused.
- `se = :bootstrap` is a parametric bootstrap, not `btergm`'s resampling of
  time steps (that is `se = :block_bootstrap`).

### Known limitations

- **EGMME** is not implemented (`method = :egmme` throws an
  `ArgumentError`): no equilibrium
  estimation from a cross-section plus durations.
- **CMLE refinements of `tergm`**: a fixed number of MCMC samples per
  iteration at a fixed thinning interval (enlarged only by the stopping
  rule's boost); no effective-sample-size-adaptive sampling, no
  missing-data CMLE. A side whose CMPLE start does not exist (boundary
  statistic, separation) is refused.
- **Offsets and constraints**: ERGM.jl's `Offset` is refused in a STERGM
  formula (as `tergm` refuses offsets in a CMLE/CMPLE fit); there is no
  `constraints=`.
- **Terms ERGM.jl does not have** are not available; a STERGM formula
  takes ERGM.jl's terms plus `Delrecip`.
- **The block bootstrap on short panels**: `se = :block_bootstrap` is refused
  below 10 transitions and warned about below 20.
- **Missing (masked) dyads** in a panel are refused; there is no
  `missing=` keyword.
- **Two-mode (bipartite) panels** are refused.
- **Self-loops**: a panel containing a loop is refused.
- **Non-separable (btergm-style) TERGMs**: the memory terms
  `EdgeStability`, `PersistentEdge`, `NewEdge` are refused in a formula
  and kept as descriptives.
- **Panels of changing composition and `networkDynamic` input** are not
  supported: every panel must have the same vertex set.

## [0.1.0] - 2026-02-09

Initial development version (not publicly released): STERGM formulation
with formation/dissolution formulas, pseudo-likelihood estimation and
sequence simulation.
