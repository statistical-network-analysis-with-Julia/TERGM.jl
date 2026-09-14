# Changelog

All notable changes to TERGM.jl are documented in this file. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
package adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - Unreleased

Release driven by the 2026-07 expert-panel review: CMPLE now really fits the
Krivitsky–Handcock separable formation/dissolution construction on
attribute-preserving auxiliary networks (fixing silently-zero nodal-term
columns), with block-bootstrap SEs, StatsAPI accessors, and a shared
`GOFResult` for goodness of fit.

### Breaking

- **`stergm_gof` returns a `Networks.GOFResult`** instead of a NamedTuple
  `(observed, simulated_mean, simulated_sd, z_scores)` of user-supplied
  network statistics. The result has a "tie changes" panel (levels
  `formed`/`persisted`) with Monte-Carlo p-values; the `statistics=` keyword
  is gone (use `n_sim`/`rng`). *Migration:* stop indexing
  `.observed`/`.z_scores` — use the `GOFResult` accessors or `show`; compute
  custom statistics on `simulate_network_sequence` output instead.
- **`STERGM` formula fields renamed** `formation_terms`/`dissolution_terms`
  to `formation`/`dissolution`; the unused `constraints` field/keyword was
  dropped, and empty formation or dissolution models now throw.
  *Migration:* access `formula.formation`; remove `constraints=`.
- **`STERGMResult` layout changed:** new `vcov`, `loglik_formation`,
  `loglik_dissolution`, `se_type` fields replace the scalar `loglik`.
  *Migration:* replace `result.loglik` with `loglikelihood(result)` (or the
  per-model fields).
- **`simulate_stergm` / `simulate_network_sequence` take coefficients
  positionally:** `simulate_network_sequence(formula, init_net, n_steps,
  θ_form, θ_diss; burnin, rng)` (was `form_coef=`/`diss_coef=` keywords);
  `simulate_stergm(prev_net, formula, θ_form, θ_diss; ...)` is the
  single-transition primitive, with a `(result, n_steps; ...)` convenience.
  The literal default `burnin` was raised 100 → 3000 and then, later in
  this release, replaced by the dyad-scaled default described below.
  *Migration:* move coefficients out of keywords.
- **Default estimation method is `:cmple`** (the honest name for what was
  fit all along). `cmle`/`method=:cmle` now throws an `ArgumentError` instead
  of silently returning a CMPLE fit stamped `:cmle`: MCMC-based CMLE is not
  implemented, and CMPLE is only the same estimator for dyad-independent
  formulas. `egmme` throws instead of returning placeholder zero
  coefficients, and is **no longer exported** — an unimplemented estimator
  should not advertise itself in the public API. It is still reachable as
  `TERGM.egmme` (so the error is informative, not an `UndefVarError`), and
  `stergm(...; method=:egmme)` still raises. *Migration:* use `method=:cmple`
  (numerically identical to the old default path); replace bare `egmme` with
  `TERGM.egmme` if you were relying on the export.
- **Removed term types:** `EdgeAge`, `Memory`, `TimeLag`, and the
  `FormationTerm`/`DissolutionTerm` wrappers (plus the phantom
  `FormationModel`/`DissolutionModel` exports). Standard ERGM terms now go
  directly into the formation/dissolution vectors and are evaluated on the
  auxiliary networks. *Migration:* drop the wrappers and pass terms
  directly; for edge-age descriptives use the new
  `edge_ages`/`mean_edge_age`.
- **`EdgeStability` redefined** as the count of dyads agreeing with the
  previous network (label `edge.stability`; was a persisting-edge count
  labeled `edgestability`). *Migration:* use `PersistentEdge` for the old
  persisting-edges meaning.
- **`STERGMModel` dropped the `times=` keyword** and the
  `.directed`/`.times` fields. *Migration:* remove `times=`.
- **`STERGMModel{T,D}` / `STERGMResult{T,D}` are parameterised on the
  panels' directedness** (panel 2026-09, item 19), mirroring
  `ERGM.ERGMModel{T,D}`: `networks::Vector{Network{T,D}}` is a concretely
  typed field, `is_directed(model)` (a `Graphs.is_directed` method) is a
  type-level fact, and every panel must be the same `Network{T,D}` — a
  mixed-directedness or mixed-vertex-type panel is refused at construction.
  `STERGMModel(formula, networks)` accepts any `AbstractVector{<:Network}`
  and collects it to the concrete panel type; `simulate_network_sequence`
  returns a `Vector{typeof(init_net)}` (was an abstract-element
  `Vector{Network{Int}}`). *Migration:* a `where T` signature written
  against `STERGMModel{T}` must become `where {T,D}`.
- **Formula validation is ERGM.jl's** (panel 2026-09, item 12): both
  formulas are now checked against **every** panel by the `public`
  `ERGM._validate_formula` — the function the `ERGMModel` constructor runs
  — instead of TERGM's own two re-implemented checks. New refusals, each an
  `ArgumentError` prefixed `formation model, panel t:` /
  `dissolution model, panel t:` and naming the term and the fix:
  a declared vertex attribute that is not set on **every** vertex of a
  panel (statnet refuses NA attribute values; before, the vertex was
  silently zero-filled); the undirected-only `Kstar`/`GWDegree`/`Degree` on
  directed panels, with `OStar`/`IStar`, `GWODegree`/`GWIDegree`,
  `ODegree`/`IDegree` named (before, a directed panel silently got
  out-stars); an `EdgeCov` covariate matrix that is not `n×n`. The
  missing-attribute and directed-only-term messages are ERGM's wording
  (`term 'nodematch.grp' refers to vertex attribute :grp, which does not
  exist on the network. Available vertex attributes: ...`). *Migration:*
  switch `Kstar`/`GWDegree`/`Degree` to their directed variants on directed
  panels; set attributes on every vertex.
- **`Delrecip` throws on undirected networks.** `compute(Delrecip(), net,
  prev)` and `change_stat(Delrecip(), net, i, j, prev)` raise an
  `ArgumentError` ("only defined for directed networks") on an undirected
  `net` instead of silently returning `0.0` (the July P1-11 silent-zero
  class); `requires_directed(Delrecip()) == true` on ERGM.jl's public trait.
- **Persistence is the one parameterisation, and the one name** (panel
  2026-09, item 31). `STERGMResult` fields `dissolution_coef`/`dissolution_se`
  are renamed `persistence_coef`/`persistence_se`, with exported accessor
  functions `formation_coef`, `formation_se`, `persistence_coef`,
  `persistence_se`. The `dissolution_*` spellings are **deprecated**, not
  removed: `result.dissolution_coef` (a `getproperty` shim) and the exported
  `dissolution_coef(result)`/`dissolution_se(result)` return the very same
  persistence vectors with a deprecation warning. tergm's `Diss()` sign
  convention (the negation) is deliberately not offered; the docs state the
  relation. `show` labels the second table `Persistence:` (was
  `Dissolution (persistence):`). *Migration:* rename `dissolution_coef` →
  `persistence_coef` (same numbers).
- **`stergm_gof` is no longer the same binding as `gof`** (`stergm_gof ===
  gof` was true in 0.1): it is a deprecated forwarding method that warns
  and calls `gof(result; kwargs...)`. `gof(fit)` is the one goodness-of-fit
  verb of the model family (panel 2026-09, item 16). *Migration:* rename
  `stergm_gof(fit; ...)` → `gof(fit; ...)`; the keywords are unchanged.
- **`show(::STERGMResult)` layout.** The header now reads `Panels: T (T−1
  transitions, n free dyads)`, `Pseudo-log-likelihood: formation …,
  persistence …` (was `dissolution …`), a `Standard errors:` line naming
  how they were obtained (`inverse Hessian of the pseudo-likelihood` or
  `per-transition block bootstrap (n replicates)`) and `Converged:
  true/false` on its own line; an unconverged fit prints `WARNING:` with
  the caveat sentence **in the header**, before any table (it used to be a
  footnote under the tables). The two blocks are `Formation:` and
  `Persistence:`. *Migration:* anything parsing the printed header; use
  `coeftable(fit)` / `confint(fit)` / `fit.converged` instead of scraping
  `show`.
- **Coefficient labels are R's, resolved against the panel.** `show` names
  coefficients with ERGM.jl 0.2's direction-aware `name(term, net)`: a
  directed panel's `GWESP(0.5)` is `gwesp.OTP.fixed.0.5`, an undirected
  one's `gwesp.fixed.0.5`; integer decays print as `gwesp.fixed.1` and
  `GWDegree` as `gwdeg.fixed.<d>` (ERGM.jl's label changes).
- **`STERGMResult` gains a `boot_replicates::Matrix{Float64}` field** (the
  `n_boot × p` block-bootstrap refits; `0 × p` for `se=:hessian`). A
  replicate whose refit has no finite coefficient is kept as a `NaN` row,
  excluded from the covariance, warned about **once**, listed by
  `approximations(fit)` and noted by `show` (before: dropped silently
  inside the loop with a warning and no record). Fewer than two finite
  replicates is now an `ArgumentError` (was an `ErrorException`).
  `STERGMResult` has one (positional, 13-field) constructor — the layout is
  new in this release, so there is no released layout to stay compatible
  with. *Migration:* construct results through `stergm`/`cmple`.
- **The memory terms `EdgeStability`, `PersistentEdge` and `NewEdge` are
  refused in a formula** (panel 2026-09 round 2). On the free dyads of a
  separable model each one's change statistic is a constant — `-1`/`+1`
  (edge.stability), `0`/`+1` (persistent.edges), `+1`/`0` (new.edges) on
  the formation/persistence rows — so beside `Edges()` the CMPLE design
  was rank-deficient (returned `converged == false` with the generic
  "perfect separation or a degenerate statistic" warning and the cause
  never named) and without it the term was `edges` under another name.
  `STERGMModel` now throws `ArgumentError: formation model: term
  'edge.stability' is constant on the free dyads of the formation model of
  a separable STERGM (its change statistic is -1 on every prior non-tie,
  i.e. the column is -edges) …`, naming the mechanism and that the terms
  are only informative in a non-separable (btergm `memory`-style) TERGM,
  which TERGM.jl does not fit. The three stay exported as descriptives
  (`compute(term, curr, prev)`; `gof` reports formed/persisted through
  them); none has a counterpart in R tergm 4. *Migration:* drop them from
  formulas; `Delrecip` is the one temporal term a formula can carry.
- **Two-mode (bipartite) panels and panels with self-loops are refused**
  at construction (panel 2026-09 round 2), as ERGM.jl refuses the same
  networks. Before, `network(6; bipartite=3)` panels fitted with every
  `i≠j` dyad as a free dyad — the 6 impossible within-mode dyads of each
  transition entered the formation model as rows that could never form
  (`nobs` 30 where 18 cross-mode free dyads exist; `Form~edges` biased
  downward) — and a `(v, v)` loop was counted by `compute(Edges())`,
  `PersistentEdge`/`NewEdge` and `gof`'s observed counts while the design,
  `nobs` and the sampler skipped it. Now `ArgumentError: TERGM (panel t):
  two-mode (bipartite) panels are not supported — the one-mode CMPLE would
  enumerate the impossible within-mode dyads as free dyads …` and
  `ArgumentError: TERGM (panel t): the panel contains 1 self-loop (at
  vertex v) …` (naming `rem_edge!(net, v, v)`). Documented under "Not
  implemented" in the README and the estimation guide.
- **`gof` returns five panels, not one** (panel 2026-09 round 2, after
  `gof.tergm`): `"tie changes"` (as before, and first), then
  `"formation statistics"` / `"persistence statistics"` — every formula
  term evaluated on the observed vs simulated Y⁺ / Y⁻, pooled over
  transitions, labelled by the fitted coefficients — and the degree
  distribution of Y_t (`"idegree"` + `"odegree"` on a directed panel,
  `"degree"` on an undirected one; levels `0` … `n−1`). The formed/persisted
  counts are the sufficient statistics of `Form~edges` / `Persist~edges`,
  which a CMPLE reproduces by construction, so the old single panel had
  no diagnostic power for any formula containing `Edges()`; the model
  statistics are where a dyad-dependent CMPLE's misfit shows, the degree
  panels the out-of-model check. The `gof` docstring and the simulation
  guide say which panel is informative for which formula. *Migration:*
  `only(g.statistics)` → `g.statistics[1]`; the seeded simulations of a
  given `rng` are unchanged (same seeds, same draws), only more is
  computed from each.
- **`stergm(networks, formation, dissolution)` takes any term collection
  and names the common slips** (panel 2026-09 round 2; ERGM.jl item 31).
  `formation`/`dissolution` are each a bare term or any collection
  `fit_ergm` accepts — a `push!`-built `Vector{Any}`, nested vectors —
  normalised by ERGM's term-list collector (`STERGM` likewise), so a
  non-term element is refused with its side, position and type and a term
  *type* gets the "did you mean `Edges()`?" hint. `stergm(terms, terms,
  networks)` and `stergm(Edges(), …, networks)` throw `ArgumentError:
  arguments are swapped: call stergm(networks, formation, dissolution) —
  the panel … comes first … (R tergm: nets ~ Form(~…) + Persist(~…))`; a
  single `Network` throws `ArgumentError: stergm needs a panel …` pointing
  at `ERGM.fit_ergm`. Before, all four were raw `MethodError`s with a wall
  of `Closest candidates`. An empty side is
  `ArgumentError: the formation model must contain at least one term …`.
- **`show(::STERGMResult)` prints `AIC: …, BIC: …`** after the
  pseudo-log-likelihood line, as `ERGMResult`'s header does; `STERGMModel`
  and `STERGM` gained compact one-line `show` methods (`STERGMModel{Int64,true}:
  3 panels of 50 vertices (directed), 4900 free dyads; formation: edges +
  mutual + nodematch.smoke; persistence: edges + nodematch.smoke`;
  `STERGM(formation: edges + mutual; persistence: edges)`) instead of the
  default struct dump that printed every panel and the raw term
  constructors. *Migration:* anything parsing those printouts.
- **Sampler `burnin` defaults to `nothing`**, resolved per side by ERGM.jl's
  dyad-scaled rule `ERGM._mcmc_defaults` — 20 Metropolis toggles per free
  dyad of that side (the non-edges of Y_{t−1} for formation, its edges for
  dissolution) — instead of the literal `3000` for every network size
  (panel 2026-09, item 24e). `simulate_stergm`, `simulate_network_sequence`
  and `gof` (which forwards `burnin`) all resolve it the same way; an
  explicit integer is honoured as given. *Migration:* pass `burnin=3000` to
  reproduce the old draws.
- **A statistic at the boundary of its attainable range follows R ergm's
  `drop` semantics** (panel 2026-09, items 24/28 applied to CMPLE). A
  `NodeMatch` no prior same-group tie of which persists, a `Triangle` on a
  triangle-free auxiliary network, … has no finite CMPLE; `cmple` used to
  "converge" on the flat asymptote and return an arbitrary large
  coefficient with a huge standard error, silently. It now warns with R
  ergm's sentence (`cmple: observed statistic(s) Persist~nodematch.grp are
  at their smallest attainable values. Their coefficients will be fixed at
  -Inf …`), returns the coefficient as `∓Inf` with standard error 0 and
  p-value 0, and fits the remaining coefficients on the free dyads the
  dropped statistic does not touch — the exact limit of the
  pseudo-likelihood. `dof(fit)` now counts the *finite* coefficients and
  `bic(fit)` uses the dyads they were estimated on (new
  `STERGMResult.n_kept` field; `nobs` stays every free dyad), as R's
  `logLik` does; `is_exact(fit)` is `false`; `approximations(fit)` and
  `show` name the fixed coefficient; `se=:bootstrap` throws an
  `ArgumentError` naming the statistic (it is at its boundary on every
  resample of the transitions, so no refit could be finite). *R tergm 4.2.2
  itself does not drop* — its `Form()`/`Persist()` operator terms bypass
  ergm's attainable-range check, so it warns `The MPLE does not exist!` and
  returns a point on the asymptote (`-19.57`, SE `1621` on the fixture's
  boundary panel). The finite coefficients agree to 1e-6 (golden fixture,
  `boundary_*` keys); where tergm prints an artefact of `glm`'s stopping
  rule, TERGM.jl prints `-Inf`.
- **A perfectly separated design is `converged == false`, and an
  unconverged fit is loud.** Separation by a *combination* of statistics
  (no cross-group tie persists: edges `→ -Inf`, nodematch `→ +Inf`) used to
  come back `converged == true` at the point where Newton met its
  tolerance on the asymptote. `cmple` now detects the asymptote (ERGM.jl's
  `_separated` test), warns with R's sentence (`cmple: the MPLE does not
  exist (perfect separation) …`) and returns `converged == false`. A fit
  that exhausts `maxiter` is warned about (`cmple: the Newton iteration did
  not converge in maxiter=… iterations on the formation model …`), as is a
  side with no free dyad on any transition (a `NaN` fit); every
  `converged == false` result carries the caveat in `approximations(fit)`
  and under `show`'s tables, and `is_exact(fit)` is `false`.
- **`gof` seeds every simulation from `rng` up front and runs them on all
  threads.** One `UInt64` seed per (transition, simulation) is drawn from
  `rng` in a fixed order and each simulation runs on its own
  `Xoshiro(seed)`, so the result is reproducible from `rng` alone and
  identical on every thread count; the simulated statistics of a given
  `rng` differ from 0.1's (which threaded `rng` through the simulations
  serially). `n_sim < 1` is an `ArgumentError`.
- **Minimum Julia raised to 1.12**; package UUID regenerated;
  `NetworkDynamic`/`Optim`/`StatsBase` dependencies dropped. *Migration:*
  upgrade Julia and re-resolve environments pinning the old UUID.

### Added

- **Golden rows for expanding terms** (round 3): `panel_stergm.toml` gains
  `undirected_grp3` (a 3-level attribute on the frozen undirected panel)
  and four tergm 4.2.2 fits — `expand_undirected` (`Form(~edges +
  nodefactor("grp3")) + Persist(~edges + nodemix("grp3"))`),
  `expand_undirected_degree` (`degree(1:2)` / `degree(0:1)`),
  `expand_directed` (`nodefactor("grp")` / `nodemix("grp")`) and
  `expand_directed_degree` (`nodemix("grp")` / `odegree(0:1)`) — term
  names, coefficients, standard errors, log-likelihood, df and nobs, all
  warning-free in R (`stopifnot` in the script). Read by the "Golden
  fixture: expanding terms (nodefactor, nodemix, degree) match tergm"
  testset.
- **Testsets** (round 3): "Expanding terms materialize per statistic with
  R's labels", "simulate_stergm names a coefficient-count mismatch and an
  absent attribute"; the slips testset gains the panel-in-second-position
  calls; the two `@allocated` pins carry an attribute term.
- **Provenanced golden fixture against a real statnet `tergm` CMPLE fit**
  (issue #8). `test/fixtures/panel_stergm.toml` freezes a tergm 4.2.2 fit of a
  simulated 8-wave, 25-actor directed panel (`Form(~edges + nodematch("grp")) +
  Persist(~edges + nodematch("grp"))`), regenerable with
  `Rscript test/fixtures/r/panel_stergm.R > test/fixtures/panel_stergm.toml`. The
  eight waves are frozen as edge lists, so both packages fit identical networks.

  The formulas are **dyad-independent on purpose**: the conditional
  pseudo-likelihood is then the conditional likelihood, CMPLE is the exact
  conditional MLE, and agreement can be asserted at 1e-6 rather than hand-waved.

  **Finding, and it is R's, not ours:** `tergm`'s CMPLE runs R's `glm` at the
  default `epsilon = 1e-8` and stops there — **1.5e-8** short of the exact optimum
  in the coefficients and **6.1e-6** short in the standard errors. TERGM.jl's
  Newton–Raphson lands on the exact optimum, reproducing a tightly-converged
  (`epsilon = 1e-14`) R `glm` to **~1e-10** in both. The fixture therefore freezes
  *both*: `coefficients`/`std_errors` are tergm as shipped (SEs compared at 1e-4,
  a tolerance that measures tergm's slack), and `exact_coefficients`/
  `exact_std_errors` are the same estimator taken to convergence (compared at
  1e-6, where TERGM.jl passes with four orders of magnitude to spare). The testset
  also asserts TERGM.jl is *closer to the exact optimum than tergm itself is*.

  Block-bootstrap SEs are compared too, at the resolution the bootstrap actually
  has: the R script reruns its btergm-style transition resampling under five
  further seeds and freezes the seed-to-seed sd of every bootstrap SE
  (0.0011–0.0042). TERGM.jl's five-seed mean lands 0.0006–0.0071 from R's — the
  two bootstraps differ by about as much as either differs from itself.

- `formation_network(prev, curr)` and `dissolution_network(prev, curr)`
  exported constructors for the Y⁺/Y⁻ auxiliary networks.
- Per-transition block-bootstrap standard errors (the btergm approach):
  `stergm(...; se=:bootstrap, n_boot, rng)`, producing a full joint `vcov`;
  `se_type` records `:hessian` vs `:bootstrap`, and `show` prints an
  honest-uncertainty caveat for dyad-dependent CMPLE fits.
- StatsAPI accessors on `STERGMResult`: `coef`, `stderror`, `vcov`,
  `loglikelihood`, `aic`, `bic`, `nobs`, `dof`.
- `gof` method on the ecosystem-wide `Networks.gof` generic (`stergm_gof`
  deprecated; see Deprecated below).
- Formula validation at model construction: attribute-based terms must
  reference a vertex attribute present on every panel, and
  direction-incompatible terms (`Mutual`, `Delrecip`) are rejected on
  undirected panels — errors instead of silently wrong fits.
- Temporal descriptives `edge_ages(networks)` and
  `mean_edge_age(networks)`.
- `has_dyad_dependent(model::STERGMModel)` — a method on ERGM.jl's
  **exported** generic (replacing TERGM's private same-named
  `_has_dyad_dependent`), `true` when either formula contains a
  dyad-dependent term; `show`, `is_exact` and `approximations` all dispatch
  on it.
- Accessor functions `formation_coef`, `formation_se`, `persistence_coef`,
  `persistence_se` (and the deprecated `dissolution_coef`/`dissolution_se`),
  each with a runnable docstring example; `simulate_stergm`, `stergm`,
  `cmple`, `cmle`, `gof` docstrings gained runnable examples.
- CI runs one cell with `JULIA_NUM_THREADS=4` (the block bootstrap refits
  on every thread through `Networks.bootstrap_cov`, and the suite asserts
  the replicates are thread-count independent).
- Testsets: "No private cross-package reach-ins" (reads `src/TERGM.jl` as
  text: no `ERGM._x`/`Networks._x` dotted reach-in survives, and every
  `_`-prefixed name imported by name from ERGM/Networks is declared `public`
  there), the NA-completeness / undirected-only / `EdgeCov`-size refusals,
  `STERGMModel{T,D}` concreteness, the persistence accessors and their
  deprecated aliases, the `bootstrap_cov`-backed bootstrap (replicates
  recorded, dropped replicates warned once and listed, serial == threaded),
  the `mh_toggle!`-backed sampler (bit-identical to a reference of the
  pre-0.2 hand loop; 0 bytes per step), and R's directed `gwesp` labels.

- **Golden fixture, second panel (`boundary_*` keys):** an 8-actor, 3-wave
  directed panel in which no same-group tie ever persists, fitted by
  `tergm(..., estimate="CMPLE")` with `Form(~edges + nodematch("grp")) +
  Persist(~edges + nodematch("grp"))` (R 4.6.1, tergm 4.2.2). Frozen: R's
  warning (`The MPLE does not exist!`), the finite coefficients and their
  standard errors both as shipped and at `glm` convergence (compared at
  1e-6, the same exact-CMLE argument as the main panel — measured slack
  8e-16), the log-likelihood (1e-5: R's differs from the exact dropped
  limit by the asymptote's tail, frozen as `boundary_asymptote_tail` =
  −1.4e-7), R's `logLik` df/nobs, and — as a record, not an assertion —
  tergm's asymptote value. The testset asserts the fixed coefficient by
  sign of `Inf`, the finite ones at the fixture's tolerances, and every
  loud-failure contract above (`@test_logs`, `approximations`, `show`,
  `dof`/`bic`/`nobs`, the bootstrap refusal).
- **PrecompileTools workload** (panel 2026-09, item 18): the fit → `show` →
  bootstrap → simulate → `gof` path on a 6-actor, 3-wave panel is compiled
  into the package image. Measured (Julia 1.12.6, one thread): first
  `stergm` fit **1.6 s → 0.04 s**, first `show` 0.57 s → 0.01 s, first
  `se=:bootstrap` 0.39 s → 0.01 s, first `simulate_stergm` 0.12 s → 0.03 s,
  first `gof` 0.26 s → 0.01 s; `using TERGM` 0.79 s → 0.71 s; package
  precompilation 2.3 s → 4.9 s. New dependency `PrecompileTools` (compat
  `1`). The workload also covers `coeftable` and `confint` (both levels):
  first `coeftable(fit)` **0.047 s → 0.000 s**, first `confint(fit)`
  **0.10 s → 0.007 s** (measured with the two lines removed from the
  workload vs. present; `using TERGM` unchanged at ~0.9 s on this machine).
- **The full StatsAPI surface on `STERGMResult`** (panel 2026-09, item 15):
  `coeftable(fit)` returns a `Networks.CoefficientTable` — `Estimate`,
  `Std.Error`, `z value`, `Pr(>|z|)` — with one row per stacked coefficient
  labelled as `summary(tergm)` labels them (`Form(1)~edges`,
  `Persist(1)~nodematch.grp`; the term part is ERGM.jl's direction-aware R
  label), built from the very same two tables `show` prints
  (`_side_tables`), with `p = 0` on a coefficient fixed at `±Inf`; the
  golden fixture pins `coeftable(fit).names == term_names` for both panels.
  `confint(fit; level=0.95)` gives Wald limits from the standard errors
  the fit reports (`level ∉ (0, 1)` is an `ArgumentError`; a fixed `±Inf`
  coefficient has the degenerate interval `[±Inf, ±Inf]`). Both are the
  ONE `StatsAPI` binding Networks.jl and ERGM.jl re-export (so `using
  ERGM, TERGM` keeps every verb single-owner) and are exported alongside
  `coef`, `stderror`, `vcov`, `loglikelihood`, `aic`, `bic`, `nobs`, `dof`.
  The test suite pins the whole surface with
  `Networks.check_statsapi(fit; strict=true)` for a Hessian fit, a
  bootstrap fit and the boundary (dropped-statistic) fit, and asserts
  `coeftable`'s numbers appear verbatim in `show`'s output.
- **The golden fixture pins the persistence accessors by name**: the R
  script's `[tolerance]` block gains `formation_coefficients`,
  `formation_std_errors`, `persistence_coefficients`,
  `persistence_std_errors` (the same numbers as the stacked vectors split
  by side, at the same 1e-6 / 1e-4 tolerances), so
  `persistence_coef(fit)`/`persistence_se(fit)` are asserted against
  tergm's `Persist()` block directly — the sign convention is tergm ≥ 4's
  `Persist()`, never `Diss()`'s negation. Fixture regenerated with Rscript
  (R 4.6.1, tergm 4.2.2); every frozen value is unchanged.
- **The missing-data vocabulary is pinned** (panel 2026-09, item 5):
  `missing_policies(stergm) == missing_policies(cmple) == (:error,)`
  (Networks' default — TERGM exposes no `missing=` keyword, asserted via
  `Base.kwarg_decl`), `supports_missing` is `false`, and the refusal a
  masked panel produces names `clear_missing_dyads!` and the panel but
  never a `:face` policy the caller could not pass (`face_ok=false`).
- **Real teaching data.** The README Quick Start, Getting Started and the
  estimation guide fit RSiena's `s50` panel (50 girls, three waves of
  directed friendship nominations, `Networks.load_dataset(:s50)`) with a
  wave-1 `:smoke` vertex attribute, instead of a synthetic toggled panel
  (panel 2026-09, item 22; the synthetic construction stays as "Building a
  panel by hand"). Every block is executed by the site's snippet checker.
- **Testsets** "Shared verbs are single-owner; no method ambiguities"
  (`TERGM.gof === Networks.gof`, `TERGM.coeftable === StatsAPI.coeftable`,
  …, `Test.detect_ambiguities(TERGM, ERGM, Networks)` empty,
  `fit_stergm === stergm`), the StatsAPI surface (`check_statsapi`,
  `confint` width `2·1.96·se` at 1e-12 and level monotonicity, `coeftable`
  labels/rows/name lookup, show ⊇ coeftable), the deprecations
  (`@test_deprecated stergm_gof`/`dissolution_coef`/`dissolution_se`,
  `@test_logs (:warn, r"deprecated")` on the property shim — `Pkg.test`
  runs with `--depwarn=yes`, so the warnings are observed, not assumed;
  `propertynames(r) == fieldnames(STERGMResult)`), the missing-data
  vocabulary, and the unconverged `show` header. 556 tests (was 461).
- **Testsets** "Unconverged and degenerate CMPLE fits are loud" (Newton
  cap, separation, empty side), "Frozen pre-refactor sampler output" (the
  sorted edge list of a seeded `simulate_stergm` draw taken with the
  committed pre-`mh_toggle!` loop, as a literal), "Keyword vocabulary
  follows the root CLAUDE.md" (`Base.kwarg_decl` of every entry point
  contains no `max_iter`/`n_sims`/`seed`), "gof is seeded per simulation
  and thread-count independent" (serial reproduction of the seeding scheme
  equals the threaded run), and the boundary golden testset. The "No
  private cross-package reach-ins" testset gains an allow-list for the
  three ERGM.jl helpers imported while a `public` request is pending
  (`_mple_fit_design`, `_warn_boundary`, `_warn_separated`), each asserted
  *not* public so the list goes red — and is emptied — the day ERGM
  declares them.

- **Golden fixture, dyad-dependent and undirected panels** (panel 2026-09
  round 2; `dyaddep_*` and `undirected_*` keys, R 4.6.1, tergm 4.2.2). The
  only pinned fit had been the directed dyad-independent panel; tergm's
  CMPLE for dyad-dependent terms and for undirected panels — the same
  estimator written by other people — was unpinned. Frozen now, as tergm
  ships them and compared at the main panel's as-shipped tolerances
  (1e-6 coefficients, 1e-4 standard errors, 1e-5 log-likelihood, with the
  justification in the fixture's `[tolerance]` block): `Form(~edges +
  mutual) + Persist(~edges + mutual)` and `Form(~edges + triangle) +
  Persist(~edges)` on the existing 25-actor directed panel, and — on a new
  frozen 14-actor, 5-wave **undirected** panel — `Form(~edges +
  nodematch("grp")) + Persist(~edges + triangle)` and `Form(~edges +
  kstar(2)) + Persist(~edges)`, each with R's `logLik`, `df`, `nobs` (364 =
  the `n(n−1)/2` unordered free-dyad enumeration), term names and the
  absence of warnings. Measured agreement: mutual 5e-10, triangle 2.5e-10,
  undirected 1e-9 on the coefficients. The testsets also assert
  `fit.converged`, which the triangle fit failed before the fix below.
- **Testsets** "gof compares model statistics and degree distributions",
  "Memory terms are refused in a separable formula" (all three, both
  sides, and the constancy itself), "Two-mode and self-loop panels are
  refused", "stergm accepts what fit_ergm accepts and names the common
  slips", "STERGMModel and STERGM print compactly; the result header has
  AIC/BIC", "CMPLE design build allocates O(blocks), not O(rows · p)"
  (`_fill_blocks!` at 0 bytes; the build's overhead beyond the returned
  arrays and the auxiliary networks independent of the row count), and the
  two golden testsets above. 860 tests (was 663).

### Changed

- **`_polish_newton` and `_specification` are gone** (reconciliation
  round). `Networks.newton_fit` adopted the scale-free convergence verdict
  TERGM requested — the Newton decrement `½·gᵀ(−H)⁻¹g < tol` beside the
  gradient norm, rounding-noise steps accepted while the gradient shrinks,
  the last full step taken at a stall whose predicted gain is below `tol` —
  so `_logistic_fit` is one call to ERGM's now-`public` `_mple_fit_design`
  with `context="cmple"` (both R sentences carry TERGM's name from there;
  the `_warn_boundary`/`_warn_separated` imports and TERGM's re-emission of
  them are deleted). ERGM's `public` `_expand_terms(terms, net)` returns the
  specification TERGM's `_specification` used to extract by unwrapping the
  twin's `base` field from outside ERGM; `_collect_terms` is `public` too,
  and the "No private cross-package reach-ins" allow-list is empty. The
  golden panel is unchanged at its tolerances (`dyaddep_triangle` still
  `converged`, gradient < 1e-8).
- **`gof` runs the formation/persistence model statistics through
  `_tcompute`** (a temporal term sees Y_{t−1}, a standard term does not) —
  the helper was unused before.
- **The CMPLE design build fills preallocated blocks through a typed
  barrier** (`_fill_blocks!` over the two term tuples, statically unrolled
  per row by `_fill_row!`), replacing one `Vector{Float64}` per free dyad
  through an abstract term vector followed by `hcat`/`permutedims` —
  ~520 B per row before (9.7 MiB for a 100-actor, 8-wave panel's design),
  now the four arrays of each transition (plus its two auxiliary networks)
  and a per-transition constant, whatever the row count. **Round-3
  correction:** the "~3 KB constant" was measured on structural terms
  only; with the raw `NodeMatch` the README teaches, every row still went
  through the attribute Dict (160 B per call; 548 KB of overhead on the
  s50 Quick-Start's 4 900 rows, 86.6 KB per 30-actor transition in
  `_fill_blocks!`). Attribute terms are now snapshotted per transition
  (`_materialized_tuple`: one `Int` per vertex plus a small Dict per
  attribute term per side per transition — ~2.2 KB per transition on the
  12-actor smoke test), and `_fill_blocks!` is **0 bytes** with
  `[Edges(), Mutual(), NodeMatch(:grp)]` (measured; the pin now uses that
  formula and bounds the build's overhead per transition plus the O(n)
  snapshot).
- **Convergence is judged by the Newton decrement when the shared kernel
  stops at the noise floor** (round 2: `_polish_newton`, see Fixed; since
  the reconciliation round the criterion lives in `Networks.newton_fit`
  itself and the local polish is deleted — see the entry above).
- `Delrecip`'s directed-only message no longer claims "R tergm's delrecip is
  likewise directed-only": tergm ≥ 4 has no `delrecip` term. Its docstring
  and the terms guide now state the definition (`Σ y_ij · y^{t−1}_ji`) and
  the provenance (btergm's `delrecip`, lag 1) explicitly, and say that no
  tergm fixture can pin it.
- `ERGM._collect_terms` joins `_mple_fit_design`/`_warn_boundary`/
  `_warn_separated` on the pending-`public` allow-list of the "No private
  cross-package reach-ins" testset (cross-repo request: `public
  _collect_terms` in ERGM.jl).
- Built-in temporal terms (`Delrecip`, `PersistentEdge`, `NewEdge`) use the
  ecosystem-wide state-independent add-direction `change_stat` convention
  and R-style labels (`delrecip`, `persistent.edges`, `new.edges`);
  `Delrecip` is directed-only.
- The dissolution model is explicitly documented as the persistence
  parameterization (positive coefficient = ties persist), following
  Krivitsky & Handcock.
- **Every shared contract is imported by name from the package that owns
  it; no private reach-in remains** (panel 2026-09, items 13 and 28).
  `Networks.z_pvalues` (the ONE floored, NaN-aware z → p helper) replaces
  the deprecated `ERGM._z_pvalues`; `Networks.check_se` validates `se=`
  (`cmple: se must be one of (:hessian, :bootstrap) (got :sandwich)`)
  instead of a hand-typed check; `Networks.newton_fit`/`logistic_derivatives`
  are imported from Networks.jl (ERGM.jl re-exports the same bindings);
  `ERGM.requires_directed(::Delrecip)` replaces `ERGM._requires_directed`;
  the two remaining `_`-names taken from ERGM.jl (`_validate_formula`,
  `_mcmc_defaults`) are declared `public` there and imported by name.
- **The block bootstrap is `Networks.bootstrap_cov`** — the ONE shared
  resampling loop — with TERGM supplying only its two callbacks (draw the
  `n_boot` transition resamples from `rng` in one call; refit the CMPLE on
  one resample). The resamples are drawn in the order the old serial loop
  drew them, so the bootstrap standard errors are unchanged; the refits run
  on every thread and are thread-count independent (all randomness lives in
  the resample draw). The golden tergm/btergm fixture is untouched.
- **The constrained Metropolis samplers run on ERGM.jl's exported
  `mh_toggle!` kernel** — one toggle loop for the whole family. The kernel
  consumes `rng` in the order the pre-0.2 hand-written loop did (proposal
  index, then one uniform for the acceptance test), so the sampled sequence
  is bit-identical (pinned against a reference implementation of the old
  loop); the term tuple is passed through a function barrier so the step is
  statically dispatched and allocation-free (pinned by `@allocated`).
- **`_logistic_fit` is ERGM.jl's MPLE design fitter** (`_mple_fit_design`
  on the Bernoulli-row design `n_tot = 1`, `n_one = y`; the same
  `Networks.newton_fit`/`logistic_derivatives` kernel as before): the
  CMPLE over the free dyads of Y⁺/Y⁻ is a logistic pseudo-likelihood and
  now fails in exactly the ways `ERGM.mple` does (boundary drop, separation
  detection, R's warning sentences under the `cmple` context). The golden
  fixture is unmoved at its 1e-6 exact-CMLE tolerance. The block bootstrap's
  refits run it with `warn=false` and turn an unconverged or `±Inf` refit
  into a `NaN` row (excluded, warned once, recorded), never a `∓Inf` that
  would poison the covariance.
- `stergm` accepts any `AbstractVector{<:Network}` of panels;
  `is_directed(model)` replaces `is_directed(model.networks[1])` throughout.
- **`show(::STERGMResult)` and `coeftable` are built from the same two
  `CoefficientTable`s** (`_side_tables(r)`: labels resolved against the
  panel, `z = θ/se`, `p` through `_pvalues`); `show` renders each through
  `Networks.print_coeftable` with the legend printed once. The
  `Converged:` verdict, the non-convergence caveat, the standard-error
  method and the replicate count moved into the header; the dyad-dependence
  caveat, the fixed-coefficient note and the excluded-replicate note stay
  under the tables.
- `Base.getproperty(::STERGMResult, ::Symbol)` (the `dissolution_*` shim)
  is `@inline`, so access to the real fields compiles to a plain
  `getfield`.
- `gof`'s docstring and the guides say `stergm_gof` is deprecated; the
  README and guides gain "Inspecting a fit" and "Deprecated" sections
  describing the header, the two blocks, `coeftable`, `confint` and the
  exact deprecation behaviour.
- **CI and Documentation workflows derive their sibling clone list from
  the `[sources]` tables** of `Project.toml`/`docs/Project.toml` (panel
  2026-09, item 29): `setup-julia` now runs first and the clone step parses
  both files with the `TOML` stdlib (skipping the package's own
  `{path = ".."}` entry), so exactly Networks and ERGM are cloned — the
  hand-maintained list also cloned NetworkDynamic, which TERGM does not
  depend on; `docs/Project.toml` drops its unused NetworkDynamic dependency
  and source (it was why the local docs environment failed to instantiate).
  The "Workflows derive their clone lists from [sources]" testset reads
  both YAML files and fails on any literal `for pkg in <names>` list or any
  `[sources]` path that is not a `../<Pkg>.jl` sibling.
- **README, Getting Started and the estimation guide say what a user gets
  today** (panel item 8): `method = :cmle` **throws** (the prose said
  "warns and falls back"); `using TERGM, ERGM, Networks, Random` (the module
  was renamed from `Network`; the three snippets failed on their first
  line); the install section describes the ordered `Pkg.add(url=…)` install
  and, for development, side-by-side clones with `Pkg.develop(path=
  "../TERGM.jl")` instead of pointing at a root workspace project that is
  not published; a **Missing data** section in the README and the
  estimation guide states the exact refusal a masked panel produces, that
  no `missing=` keyword exists, and why (CMPLE would enumerate an
  unobserved dyad as an observed row); formula validation, persistence
  naming, the block bootstrap's recorded replicates and the "Not
  implemented" surface (MCMC CMLE, EGMME, masked panels) are documented
  with the exact error a user sees.
- **Every export has a docstring with a runnable `# Example`:**
  `EdgeStability`, `Delrecip`, `PersistentEdge`, `NewEdge`, `edge_ages`,
  `mean_edge_age`, `formation_network`, `dissolution_network`, `STERGM`
  and `STERGMResult` gained examples on the 5-actor two-panel fixture of
  the test suite (each executed and its stated values verified), joining
  the accessors, `stergm`/`cmple`/`cmle`, `simulate_*`, `gof`,
  `coeftable`/`confint` and `STERGMModel` from earlier in this release; the
  "Every export is documented with an example" testset pins it. With the
  workflow testset, 663 tests (was 556).
- `CLAUDE.md` rewritten to the current behaviour (development commands
  including the strict docs build, the snippet checker and the clean-depot
  check; the four architecture sections; the shared kernels; the boundary
  fixture panel; cross-package hygiene; conventions). The stale "no
  re-exports from dependencies" convention is gone — the StatsAPI verbs
  *are* re-exported, as the one shared binding.

### Deprecated

Both spellings keep working for this release, emit a deprecation warning
(`Base.depwarn`; visible under `--depwarn=yes`, which `Pkg.test` sets) and
return exactly what the new name returns; they are removed in the release
after 0.2.0.

- `stergm_gof(result; kwargs...)` → `gof(result; kwargs...)`
  (`Base.@deprecate`; still exported so old code warns instead of failing).
- `dissolution_coef(result)` / `dissolution_se(result)` and the
  `result.dissolution_coef` / `result.dissolution_se` properties →
  `persistence_coef(result)` / `persistence_se(result)` (same vectors).

### Fixed

- **Multi-column terms were fit as ONE pooled, mislabelled column** (panel
  2026-09 round 3, critical). `STERGMModel` bound the raw terms and skipped
  the expansion step every `ERGMModel` runs (ERGM.jl's public
  `_materialize`), so `NodeFactor(:grp)` on a 3-level attribute became one
  column labelled `nodefactor.grp` holding the SUM of R's `nodefactor.grp.b`
  and `nodefactor.grp.c`, `NodeMix(:grp)` one `mix.grp` column instead of
  the mixing cells, and `Degree(0:2)` constructed a model that crashed in
  the row fill with ERGM's "Expand via ERGMModel" message. Against tergm
  4.2.2, `Persist(~edges + istar(2) + nodefactor("grp"))` returned two
  nodefactor coefficients where TERGM returned one, with a different
  `edges`/`istar2` estimate — a silently different model, with no error,
  warning or documentation. The constructor now expands both sides through
  `_materialize` and stores the expansion in `model.formula` as plain
  single-statistic terms (`NodeFactor(:grp; level="b")`, `NodeMix(:grp,
  "a", "b")`, `Degree(1)`; every other term unchanged), so `_cmple_blocks`,
  `_sample_constrained`, `_transition_stats!`, `_labels`, `show`,
  `coeftable` and `has_dyad_dependent` all see one term per statistic under
  R's labels. Levels are resolved from panel 1, which coincides with
  tergm's union-of-panels resolution exactly when no later panel carries a
  value panel 1 lacks: a level absent from a later panel has an all-zero
  column on that transition (as in R), a new one is refused
  (`ArgumentError` naming the side, the panel and the value); the attribute values are
  snapshotted **per transition** from the auxiliary network the statistics
  are evaluated on, so a time-varying attribute is honoured exactly as the
  raw term used to read it. `simulate_stergm(prev, formula, θf, θd)`
  validates and expands the formula against `prev` the same way, so a raw
  multi-level formula takes one coefficient per level. Pinned by the
  "Expanding terms materialize per statistic with R's labels" testset
  (labels equal to `ERGMModel`'s for NodeFactor/NodeMix/Degree ranges,
  the fit equal to the hand-expanded formula's, the per-transition
  snapshot, the level-set refusal) and by four new provenanced rows of
  `panel_stergm.toml` (`expand_*`: `nodefactor` + `nodemix` and
  `degree(1:2)` + `degree(0:1)` on the undirected panel given a 3-level
  `grp3`, `nodefactor` + `nodemix` and `nodemix` + `odegree(0:1)` on the
  directed panel), term names asserted verbatim and the numbers at the
  as-shipped 1e-6 / 1e-4 / 1e-5 tolerances.
- **The "Coming from R tergm" table equated tergm's `nsim` with
  `simulate_stergm(fit, k)`'s steps** (round 3). `nsim` is the number of
  replications (independent chains) and `time.slices` the steps per chain;
  `simulate_stergm(fit, k)` is one chain `k` steps forward, i.e.
  `simulate(fit, nsim = 1, time.slices = k)`. The row now says so and a
  second row gives the `k`-replicate idiom; the simulation guide states it
  too.
- **Six accessor docstring examples fit a boundary panel** (round 3):
  `formation_coef`/`formation_se`/`persistence_coef`/`persistence_se`/
  `dissolution_coef`/`dissolution_se` each built a two-panel fixture on
  which one side had no finite CMPLE, so the first thing a reader tried
  from `?formation_coef` printed an unexplained `±Inf` warning. All six now
  use the two-sided-finite fixture of `stergm`'s own example (t0: 1→2, 3→4;
  t1: 1→2, 2→3), with `approximations(fit) == String[]` as the visible
  proof, and the closed-form values as comments.
- **`stergm([Edges()], networks, [Edges()])` — the panel in the second
  position — was a raw `MethodError`** (round 3) where every other slip is
  an `ArgumentError` in words. Two more throwing methods (`(AbstractVector,
  AbstractVector{<:Network}, Any)` and the panel in both first slots) plus
  the disambiguators route it, and a panel in two or three slots, to the
  "arguments are swapped" message; `Test.detect_ambiguities(TERGM)` stays
  empty.
- **`simulate_stergm`'s coefficient-length error named neither the counts
  nor the terms** (round 3). It now reads `simulate_stergm: the formation
  model has 1 term(s) (edges) but 2 formation coefficient(s) were given
  ([1.0, 2.0]); pass one per term — formation_coef(fit) for a fitted
  model, not the stacked coef(fit).` (same for the dissolution side, with
  the expanded term list), and a nodal term whose attribute the starting
  network lacks is refused in ERGM's words (`simulate_stergm: formation
  model: term 'nodematch.grp' refers to vertex attribute :grp, which does
  not exist …`) instead of failing inside the sampler. Both pinned.
- **A dyad-dependent CMPLE was reported `converged == false` at the
  optimum** (panel 2026-09 round 2). `Form(~edges + triangle)` on the
  golden panel came back unconverged after `maxiter=100` with the generic
  "perfect separation or a degenerate statistic" warning, while R's `glm`
  converged on the same point: `Networks.newton_fit` rejects the last
  Newton step because the log-likelihood is flat at its floating-point
  noise floor (a 5e-12 "decrease"), then declares convergence only if
  `norm(grad) < sqrt(tol)` — a criterion in the gradient's own units,
  which for a triangle column running to 33 over 3 000 rows is 2e-4 one
  step short of the optimum. `_logistic_fit` now evaluates the scale-free
  Newton decrement `½·gᵀ(−H)⁻¹g` at such an iterate (via the same shared
  `logistic_derivatives`); when it is below `tol` the last full step is
  taken and `newton_fit` itself issues the verdict and covariance at the
  polished point (one linear solve, not a second optimizer; a fit stopped
  by the Newton cap or with a large decrement stays unconverged). The
  triangle fit now agrees with tergm to 2.5e-10 with `converged == true`
  (`dyaddep_triangle_*` fixture keys).
- Dead code removed: the two positional `STERGMResult` compatibility
  constructors (no released layout ever had 11 or 12 fields).
- **Attribute-preserving auxiliary networks (critical).** `_copy_net` (and
  the Y⁺/Y⁻ builders) previously copied only edges, so nodal terms
  (`NodeMatch`, `NodeCov`, ...) in formation/dissolution formulas — an
  advertised feature — produced all-zero design-matrix columns and garbage
  coefficients under the default CMPLE path. Copies now go through
  `Base.copy(::Network)` and covariates survive.
- CMPLE now fits formation on Y⁺ (free dyads = prior non-edges) and
  dissolution on Y⁻ (free dyads = prior edges) with change statistics
  evaluated on the auxiliary networks; the previous code evaluated change
  statistics on the current network directly.

### Performance

- **The CMPLE derivative loop no longer allocates (review finding 15).**
  `_logistic_fit` carried its own copy of the logistic loop with a per-row
  `(pr*(1-pr)) .* (x * x')` inside it — a fresh `p×p` matrix on every one of the
  design rows of every Newton evaluation, **471 KB per evaluation** on a
  25-actor, 8-wave panel. It now runs on the shared `ERGM.logistic_derivatives`
  (the same builder ERGMMulti and ERGMRank use): **192 bytes** per evaluation,
  independent of the number of rows, and **4.0x faster** (0.460 ms -> 0.114 ms).
  Pinned by an `@allocated` regression test. The summation order moves from
  row-wise accumulation to BLAS, so the arithmetic is not bit-identical — but
  the fitted coefficients are: measured against the old loop on the same
  design, **max|Δθ| = 0.0**. Newton's last step is quadratically convergent, so
  a last-ulp difference in the gradient and Hessian does not move the fixed
  point. The golden `tergm` CMPLE fixture is unmoved at its 1e-6 exact-CMLE
  tolerance.
- **The constrained Metropolis step allocates nothing** (panel 2026-09, §6:
  the `Vector{AbstractERGMTerm}` dispatch at the old `:853-855`). Measured
  on the 12-actor test network with `[Edges(), Delrecip()]` at a
  never-accept coefficient: **16 bytes per toggle before** (the committed
  hand loop; 79 B/step in the panel's three-term measurement) → **0 bytes
  per toggle** on `ERGM.mh_toggle!` through the typed term-tuple barrier.
  Pinned by `@allocated` (2 000 toggles allocate exactly what 1 000 do).
  Output bit-identical (frozen literal testset). **Round-3 correction:**
  that measurement, too, was structural-terms-only; a raw `NodeMatch` in
  the formula cost **112 bytes per toggle** through the attribute Dict
  (5.2 MB over the s50 Quick-Start's 46 740 default toggles per side).
  The sampler now snapshots the attribute terms once per draw
  (`_materialized_tuple` against Y_{t−1}) and is **0 bytes per toggle
  with `[Edges(), Delrecip(), NodeMatch(:grp)]`** — the pin uses that
  formula, and the bit-identity reference loop reads the raw term.
- **The block bootstrap is thread-parallel.** 400 replicates on the golden
  25-actor, 8-wave panel: **0.307 s** (committed serial loop) → **0.314 s**
  on 1 thread (same work, now `Networks.bootstrap_cov`) → **0.114 s on 4
  threads** (2.7x), with identical standard errors to 1e-15 (the resamples
  are drawn from `rng` up front in the old order). `gof(fit; n_sim=100)`
  is threaded the same way.
- Estimation builds per-transition design blocks once and fits via the
  shared `ERGM.newton_fit` (Newton with step halving), with numerically
  stable `log1p(exp(...))` log-likelihoods, replacing the hand-rolled
  per-iteration Newton loop.
- First `coeftable(fit)` 0.047 s → 0.000 s and first `confint(fit)`
  0.10 s → 0.007 s from the extended precompile workload (see Added).

## [0.1.0] - 2026-02-09

Initial release: STERGM formulation with formation/dissolution formulas,
CMLE-labeled pseudo-likelihood estimation, and sequence simulation.
