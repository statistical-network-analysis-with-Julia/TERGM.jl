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

- **`stergm_gof` returns a `Network.GOFResult`** instead of a NamedTuple
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
  Default `burnin` raised 100 → 3000. *Migration:* move coefficients out of
  keywords.
- **Default estimation method is `:cmple`** (the honest name for what was
  fit all along); `method=:cmle` still runs but warns once that MCMC-based
  CMLE is not implemented and falls back to CMPLE, and `egmme` now throws
  instead of returning placeholder zero coefficients. *Migration:* use
  `method=:cmple` (numerically identical to the old default path).
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
- **Minimum Julia raised to 1.12**; package UUID regenerated;
  `NetworkDynamic`/`Optim`/`StatsBase` dependencies dropped. *Migration:*
  upgrade Julia and re-resolve environments pinning the old UUID.

### Added

- `formation_network(prev, curr)` and `dissolution_network(prev, curr)`
  exported constructors for the Y⁺/Y⁻ auxiliary networks.
- Per-transition block-bootstrap standard errors (the btergm approach):
  `stergm(...; se=:bootstrap, n_boot, rng)`, producing a full joint `vcov`;
  `se_type` records `:hessian` vs `:bootstrap`, and `show` prints an
  honest-uncertainty caveat for dyad-dependent CMPLE fits.
- StatsAPI accessors on `STERGMResult`: `coef`, `stderror`, `vcov`,
  `loglikelihood`, `aic`, `bic`, `nobs`, `dof`.
- `gof` method on the ecosystem-wide `Network.gof` generic (`stergm_gof`
  kept as an alias).
- Formula validation at model construction: attribute-based terms must
  reference a vertex attribute present on every panel, and
  direction-incompatible terms (`Mutual`, `Delrecip`) are rejected on
  undirected panels — errors instead of silently wrong fits.
- Temporal descriptives `edge_ages(networks)` and
  `mean_edge_age(networks)`.

### Changed

- Built-in temporal terms (`Delrecip`, `PersistentEdge`, `NewEdge`) use the
  ecosystem-wide state-independent add-direction `change_stat` convention
  and R-style labels (`delrecip`, `persistent.edges`, `new.edges`);
  `Delrecip` is directed-only.
- The dissolution model is explicitly documented as the persistence
  parameterization (positive coefficient = ties persist), following
  Krivitsky & Handcock.

### Fixed

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

- Estimation builds per-transition design blocks once and fits via the
  shared `ERGM.newton_fit` (Newton with step halving), with numerically
  stable `log1p(exp(...))` log-likelihoods, replacing the hand-rolled
  per-iteration Newton loop.

## [0.1.0] - 2026-02-09

Initial release: STERGM formulation with formation/dissolution formulas,
CMLE-labeled pseudo-likelihood estimation, and sequence simulation.
