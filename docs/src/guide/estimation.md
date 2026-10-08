# STERGM Estimation

## CMPLE

[`stergm`](@ref) defaults to `method = :auto`: [`cmle`](@ref) (tergm's
`estimate = "CMLE"`) when either formula has a dyad-dependent term,
[`cmple`](@ref) otherwise. R's `tergm()` has no default estimator — it
requires `estimate=` — so this default is TERGM.jl's choice — for a
dyad-independent formula the two are the same estimator, and the CMPLE
needs no MCMC. The rule is ERGM.jl's (`ERGM.resolve_method`), so `:auto`
means the same for `fit_ergm`. A keyword that only the other estimator
takes (`se = :bootstrap` on a dyad-dependent formula) is an `ArgumentError`
naming it: pass `method = :cmple` explicitly.

[`cmple`](@ref) is conditional maximum
pseudo-likelihood on the Krivitsky–Handcock auxiliary networks:

- **formation rows**: for every transition, the non-edges of ``Y_{t-1}``;
  response = presence in ``Y_t``; change statistics evaluated on
  ``Y^+ = Y_{t-1} \cup Y_t``;
- **dissolution rows**: the edges of ``Y_{t-1}``; response = persistence
  into ``Y_t``; change statistics evaluated on
  ``Y^- = Y_{t-1} \cap Y_t``.

Rows pool across transitions; separability makes the two logistic
likelihoods independent, each maximized by ERGM.jl's MPLE design fitter on
the ecosystem's shared `NetworkCore.newton_fit` optimizer (Newton-Raphson
with step-halving; `maxiter`, `tol`) and `NetworkCore.logistic_derivatives`.
A CMPLE is a logistic pseudo-likelihood over the free dyads, so it fails
in exactly the ways ERGM.jl's MPLE does — and none of them is silent (see
[Degenerate fits](@ref) below).

Model construction validates both formulas against **every** panel with
ERGM.jl's own `ERGM.Extension.validate_formula` — the function the `ERGMModel`
constructor runs — so a STERGM formula is held to exactly the rules an
ERGM formula is held to. The `ArgumentError` is prefixed with the side and
the panel (`formation model, panel 2: ...`) and names the term and the fix:

- a vertex attribute a term declares must exist on the panel **and be set
  on every vertex** — statnet refuses NA attribute values, and a silently
  zero-filled vertex would be a wrong design column;
- directed-only terms (`Mutual`, `Delrecip`) are refused on undirected
  panels; `Delrecip` also throws from `compute`/`change_stat` on an
  undirected network instead of returning `0.0`;
- undirected-only terms (`Kstar`, `GWDegree`, `Degree`) are refused on
  directed panels with the directed variant named (`OStar`/`IStar`,
  `GWODegree`/`GWIDegree`, `ODegree`/`IDegree`), as R ergm does;
- an `EdgeCov` covariate matrix must be `n×n`.

Every ERGM.jl term is therefore a STERGM term, including the 0.2
additions `Concurrent`, `DegRange`/`IDegRange`/`ODegRange`, `GWNSP`,
`MeanDeg`, `Density`, `Sender`/`Receiver`, `CyclicalTies`/`TransitiveTies`
and `ERGM.TriadCensus`. One ERGM
term is refused: `Offset(term, coef)` (`formation model: offset terms …
are not implemented in TERGM.jl`), because the estimators here would
estimate its coefficient instead of holding it fixed.

Both formulas are then **expanded** as an `ERGMModel` expands its terms
(ERGM.jl's `materialize`): a multi-level `NodeFactor(:grp)` becomes one
`nodefactor.grp.<level>` statistic per non-base level, a multi-cell
`NodeMix(:grp)` one `mix.grp.<l1>.<l2>` per selected cell (statnet's cell
order, first cell dropped), and `Degree(0:2)` / `IDegree` / `ODegree` one
`degree<d>` per degree — one coefficient per statistic, under the labels
`summary(tergm)` prints. `result.model.formula` holds the expanded,
single-statistic terms, so `formation_coef(result)` is indexed by them.
Levels, cells and degrees are resolved from the first panel, exactly as
`tergm` resolves them from the union of the panels: a level absent from a
later panel has an all-zero column on the transition that starts there,
and a later panel carrying a level the first one lacks — a column R would
fit and this model does not have — is refused, naming the side, the panel
and the value. The attribute *values* are snapshotted
per transition from the panel it starts at (the auxiliary networks carry
Y_{t−1}'s attributes), so a time-varying attribute enters the design as
it would be read on that transition. The golden fixture pins the expansion
against tergm on both a directed and an undirected panel (`nodefactor`,
`nodemix`, `degree(1:2)`, `odegree(0:1)`).

Three refusals are TERGM's own, because only a separable panel model has
them (each an `ArgumentError` naming the panel or the side):

- a **two-mode (bipartite) panel** (`TERGM (panel t): two-mode (bipartite)
  panels are not supported …`) — the one-mode CMPLE would enumerate the
  impossible within-mode dyads as free formation rows and bias the fit;
- a panel containing a **self-loop** (`TERGM (panel t): the panel contains
  1 self-loop (at vertex v) …`) — the statistics would count it while the
  free dyads, `nobs` and the sampler never touch the diagonal (ERGM.jl
  refuses both of these networks too);
- the memory terms `EdgeStability`, `PersistentEdge`, `NewEdge` in either
  formula (`formation model: term 'edge.stability' is constant on the free
  dyads …`) — see [Temporal Terms](@ref).

The term lists themselves are anything `fit_ergm` accepts — a
`push!`-built `Vector{Any}`, a bare `Edges()`, nested vectors — and the
common slips are words, not `MethodError`s: `stergm(terms, terms,
networks)` is `arguments are swapped: call stergm(networks, formation,
dissolution) …`, a single network is `stergm needs a panel …`, a non-term
element is named with its position and side.

Coefficient labels are R's, resolved against the panel with ERGM.jl's
direction-aware `name(term, net)`: a directed panel's `GWESP(0.5)` is
`gwesp.OTP.fixed.0.5`, an undirected one's `gwesp.fixed.0.5`.

**Validated against tergm 4.2.2** (`test/fixtures/panel_stergm.toml`, with
provenance): the dyad-independent `edges + nodematch` panel at 1e-6 on the
exact optimum, and — since 0.2 — the dyad-dependent `Form(~edges + mutual)
+ Persist(~edges + mutual)` and `Form(~edges + triangle) + Persist(~edges)`
fits on the same directed panel, plus `Form(~edges + nodematch) +
Persist(~edges + triangle)` and `Form(~edges + kstar(2)) + Persist(~edges)`
on a 14-actor **undirected** 5-wave panel (`nobs` 364 = R's), all at 1e-6 /
1e-4 on tergm as shipped, and `concurrent`, `degrange` and `gwnsp` fits on
both panels (`terms_*` keys).

!!! warning "CMPLE vs CMLE"
    For **dyad-independent** terms CMPLE *is* the conditional MLE. For
    dyad-dependent terms (`Mutual`, `GWESP`, ...) it is an approximation:
    the point estimates are biased in finite samples and the naive
    pseudo-likelihood standard errors are not calibrated, so the default
    fit prints no z, p or interval (see
    [Standard errors and inference](@ref)). Use `se = :bootstrap` or
    `method = :cmle`.

## CMLE

[`cmle`](@ref) (`stergm(...; method = :cmle)`) is the conditional maximum
likelihood estimator, R `tergm`'s `estimate = "CMLE"`. By separability the
conditional likelihood is the product of a formation and a persistence
likelihood, each an exponential family on a constrained sample space —
``Y^+_t`` ranges over the networks containing ``Y_{t-1}``, ``Y^-_t`` over
its sub-networks — whose sufficient statistic is the model statistics
summed over transitions. Each side is fit on its own:

- a side **without a dyad-dependent term** is fit exactly — its CMPLE is
  its CMLE — and no MCMC is run. For a dyad-independent formula `cmle`
  returns exactly the numbers of `cmple`, as `tergm` does;
- a side **with a dyad-dependent term** is fit by Monte-Carlo MLE, started
  at its CMPLE. Every iteration runs one Metropolis chain per transition
  on the free dyads of that side (the sampler of
  [`simulate_stergm`](@ref)), started at the observed auxiliary network,
  and records `n_samples` draws of the pooled statistics, `interval`
  toggles apart (default: two per free dyad). The update is ERGM.jl's:
  ``\theta \leftarrow \theta + \gamma\,\hat\Sigma^{-1}(g_{obs} - \bar g)``
  with Hummel's step length ``\gamma``. The stopping rule is ERGM.jl's
  too, hence R ergm 4's `confidence` rule: at ``\gamma = 1`` the fit has
  converged when the estimating equation at the updated coefficients lies,
  with 99 % confidence, inside a tolerance region (`conv_precision`,
  `conv_confidence`); when the test fails near the solution, the next
  sample is enlarged as R does. `termination = :hotelling` selects the
  older rule (t-ratios below `conv_threshold` and a non-significant
  Hotelling T² test).

Standard errors are the inverse Fisher information estimated from the
final sample plus the Monte-Carlo error of the estimate
(`se_method(fit) == :fisher`; `fit.mcmc.formation.mcmc_se` is the
Monte-Carlo part alone). The log-likelihood is path-sampled by Simpson's
rule from a dyad-independent reference (`bridge_rungs = 0` skips it).
`fit.mcmc` holds the per-side diagnostics that `show` prints. All
randomness comes from `rng` — one seed per chain is drawn up front — so a
fit is reproducible and identical at any thread count.

```julia
using TERGM, ERGM, NetworkCore, Random

s50 = load_dataset(:s50)
networks = [copy(w) for w in s50.friendship]
fit = stergm(networks, [Edges(), Mutual()], [Edges(), Mutual()];
             method = :cmle, rng = Xoshiro(1))
fit.converged                        # true
se_method(fit)                       # :fisher
fit.mcmc.formation.iterations        # MCMLE iterations of the formation side
confint(fit)                         # 4×2 Wald limits
cmple_fit = stergm(networks, [Edges(), Mutual()], [Edges(), Mutual()]; method = :cmple)
coef(fit) .- coef(cmple_fit)         # what the pseudo-likelihood approximation costs
```

**Validated against the exact answer and against tergm.** With `mutual` as
the only dyad-dependent term the conditional likelihood factorises over
pairs of actors, so the exact CMLE, its standard errors and its maximum
can be computed by enumeration. `test/fixtures/cmle_stergm.toml` (with
provenance) freezes them for two panels, beside `tergm`'s CMLE under 11
seeds; `cmle` is compared with the exact maximiser at four times `tergm`'s
own seed-to-seed standard deviation. On the fixture's reciprocity panel
the CMPLE of `Persist~mutual` is 1.37 where the exact CMLE is 1.60. The
same fixture freezes `tergm`'s CMLE of the tergm tutorial's formula,
`Form(~edges + mutual + cyclicalties + transitiveties) + Persist(~the
same)`, which `cmle` reproduces within `tergm`'s seed-to-seed spread
(`[Edges(), Mutual(), CyclicalTies(), TransitiveTies()]` on both sides).

A side whose CMPLE start does not exist — a statistic at the boundary of
its attainable range (no finite MLE exists either), a separated design —
is refused with an `ArgumentError`; a fit that exhausts `maxiter` is
returned with `converged == false`, a warning and the caveat in `show`.
`method = :egmme` raises: EGMME is not implemented, and no placeholder
estimates are returned. `egmme` is not exported; it remains callable as
`TERGM.egmme` purely so that the error explains itself.

## Degenerate fits

Nothing about a bad fit is silent:

- **The Newton cap.** A fit that exhausts `maxiter` is returned with
  `converged == false`, a warning (`cmple: the Newton iteration did not
  converge in maxiter=… iterations on the formation model …`), the caveat
  under `show`'s tables and an entry in `approximations(result)`;
  `is_exact(result)` is `false`. The verdict is scale-free: when the shared
  Newton kernel stops because the log-likelihood is flat at its
  floating-point noise floor, the fit counts as converged only if the
  remaining Newton decrement ``\tfrac12 g^\top(-H)^{-1}g`` is below `tol`
  — in which case the last full step is taken and the optimum reported —
  so a `Triangle` formation model (a change statistic in the tens over
  thousands of rows) converges where R's `glm` does instead of being
  reported unconverged one step short of the same point.
- **A statistic at the boundary of its attainable range** — a `NodeMatch`
  no prior same-group tie of which persists, a `Triangle` on a
  triangle-free auxiliary network, … — has no finite CMPLE. As R ergm does
  under its default `drop=TRUE`, `cmple` warns

  ```
  cmple: observed statistic(s) Persist~nodematch.grp are at their smallest
  attainable values. Their coefficients will be fixed at -Inf (no finite
  maximum pseudo-likelihood estimate exists; R ergm reports the same). …
  ```

  returns that coefficient as `-Inf` (`+Inf` at the largest value) with
  standard error 0 and p-value 0, and fits the remaining coefficients on
  the free dyads the dropped statistic does not touch — the exact limit of
  the pseudo-likelihood. `dof(result)` counts the finite coefficients and
  `bic` uses the dyads they were estimated on (`result.n_kept`), as R's
  `logLik` does; `nobs` stays every free dyad. `approximations(result)`
  and `show` name the fixed coefficient (`coefficient(s)
  Persist~nodematch.grp fixed at -Inf …`), and both bootstraps are
  refused with an `ArgumentError` — a panel cannot be simulated at an
  infinite coefficient, and the statistic is at its boundary on every
  resample of the transitions.

  R `tergm` 4.2 itself does *not* drop: its `Form()`/`Persist()` operator
  terms bypass ergm's attainable-range check, so it warns `The MPLE does
  not exist!` and returns whatever point `glm`'s default stopping rule
  reached on the asymptote (about `-19.6` with a standard error of
  `1600`). The finite coefficients are identical in both packages to
  1e-6 — the golden fixture's `boundary_*` panel pins them — but where
  tergm prints an arbitrary large number TERGM.jl prints `-Inf`.
- **Perfect separation by a combination of statistics**, which the
  boundary test cannot see (no *cross*-group tie persists, say: edges
  `→ -Inf`, nodematch `→ +Inf`, their sum finite), leaves the
  pseudo-likelihood without a maximum. R warns `The MPLE does not exist!`;
  `cmple` decides it exactly (the ecosystem's shared linear-programming
  verdict, `NetworkCore.logistic_separation`, through ERGM.jl's design
  fitter), warns (`cmple: the CMPLE does not exist (separation) …`, naming
  the separating terms) and returns the fit with `converged == false` and
  its z values, p-values and confidence intervals withheld (NaN), rather
  than the point where Newton met its tolerance on the flat asymptote.
- **A side with no free dyad on any transition** (no prior tie to persist,
  or no prior non-tie to form) is a `NaN` fit, warned about in those
  words, with `converged == false`.

## Missing data

A panel with masked (unobserved) dyads is refused when the formula is
bound to it: `STERGMModel` (and so `stergm`/`cmple`/`cmle`) throws
`ArgumentError: TERGM (panel t) does not support missing (unobserved)
dyads, but the network has k masked dyad(s). …`, naming the panel and
`clear_missing_dyads!`. There is no `missing=`
keyword on any TERGM entry point (`missing_policies(stergm) == (:error,)`,
`supports_missing(stergm) == false`; a fitted result reports
`missing_method(fit) == :rejected`), and that is deliberate rather than a
gap in the keyword list: CMPLE enumerates every free dyad of ``Y^+``/``Y^-``
as one observed logistic row, so a masked dyad would enter the design at
its face value and be fitted as data. ERGM.jl's constrained missing-data
MCMLE conditions on the observed dyads of a single network; a
pseudo-likelihood has no such conditioning to offer. Clear or impute the
dyads before fitting.

```julia
using TERGM, ERGM, NetworkCore

t0 = network(6; directed=true); add_edge!(t0, 1, 2); add_edge!(t0, 3, 4)
t1 = copy(t0); add_edge!(t1, 2, 3)
set_missing_dyad!(t1, 4, 5)                  # 4→5 unobserved at wave 2
try
    stergm([t0, t1], [Edges()], [Edges()])
catch e
    println(e.msg)                           # "TERGM (panel 2) does not support missing (unobserved) dyads, but the network has 1 masked dyad. ... `clear_missing_dyads!(net)` ..."
end
missing_policies(stergm)                     # (:error,)
```

## Standard errors and inference

`cmple`'s `se` keyword selects the standard errors (validated by
`NetworkCore.check_se`):

| `se =` | Standard errors | z, p, `confint` |
|---|---|---|
| default | inverse Hessian of the pseudo-likelihood | reported for a dyad-independent formula (exact); **withheld** for a dyad-dependent one |
| `:bootstrap` | parametric bootstrap | reported |
| `:hessian` | the same inverse Hessian | reported — the naive Wald table, by explicit request |
| `:block_bootstrap` | block bootstrap over transitions | reported; refused below 10 transitions |

**Why the default withholds.** For a formula with a dyad-dependent term
the pseudo-likelihood treats dependent dyads as independent observations,
and its inverse Hessian is not a calibrated variance. How badly depends on
the term. In simulation (30 actors, 200–300 panels per configuration,
coverage of nominal 95 % Wald intervals):

| | 1–2 transitions | 3–5 | 8–10 | 15–30 |
|---|---|---|---|---|
| Hessian, dyad-independent formula | 0.91–0.97 | 0.92–0.98 | 0.94–0.97 | 0.93–0.97 |
| Hessian, `Mutual` on both sides | 0.89–0.96 | 0.91–0.96 | 0.90–0.93 | – |
| Hessian, `GWESP(0.5)` on both sides | 0.70–0.91 | 0.78–0.92 | 0.80–0.90 | – |
| parametric bootstrap, `Mutual` | 0.95–0.98 | 0.94–0.96 | 0.91–0.95 | – |
| parametric bootstrap, `GWESP(0.5)` | 0.94–0.98 | 0.95–0.96 | 0.94–0.98 | – |
| block bootstrap, dyad-independent | 0.57–0.65 | 0.74–0.88 | 0.87–0.92 | 0.91–0.95 |

The naive errors are close to nominal for `Mutual` and well short of it
for `GWESP` (the low end of each GWESP range is the formation side). The
package cannot tell which case a formula is in, so — as ERGM.jl's `mple`
does — the default fit of a dyad-dependent formula prints its estimates
and naive standard errors with `NaN` in the z and p columns and a note,
`confint` throws an `ArgumentError`, `fit.inference_withheld` is `true`
and `approximations(fit)` records it. `se = :hessian` is the written
opt-in to the naive table (what `summary()` of a `tergm` CMPLE fit shows).

**The parametric bootstrap** (`se = :bootstrap`) simulates, `n_boot`
times, every transition from the fitted model *conditional on the observed
previous panel* — the conditional distribution the CMPLE's likelihood is
about — refits the CMPLE on the simulated transitions and reports the
empirical covariance of the refits. It works on a two-wave panel, and it
is the same keyword, and the same procedure, as ERGM.jl's
`mple(se = :bootstrap)`: across the family `se = :bootstrap` is the
parametric bootstrap. The point estimates
are unchanged (they may still be biased: the bootstrap calibrates the
uncertainty of the CMPLE, it does not turn it into the CMLE);
`vcov(result)` becomes the joint bootstrap covariance and
`result.boot_replicates` holds the `n_boot × p` refits. `burnin` is the
number of Metropolis toggles per simulated side (default: 20 per free
dyad).

**The block bootstrap** (`se = :block_bootstrap`, `btergm`'s scheme) resamples
whole transitions with replacement. Its resampling units are the
transitions, so on a short panel there is little to resample — two
transitions have three distinct resamples — and its standard errors are
too small (last row of the table). It is therefore **refused below 10
transitions** and warned about below 20, and it is not a remedy for naive
standard errors. (In the same simulation, multiplying its standard errors
by ``\sqrt{T/(T-1)}`` and using ``t_{T-1}`` quantiles restored 0.93–0.97
coverage at every length; the estimator is kept as `btergm` defines it.)

Both bootstraps run through the ecosystem's one shared loop,
`NetworkCore.bootstrap_cov`: everything random is drawn from `rng` up front,
the refits run on every thread, and the result is thread-count
independent. A replicate whose refit has no finite coefficient is excluded
from the covariance, warned about once, and recorded in
`approximations(result)`.

```julia
using TERGM, ERGM, NetworkCore, Random

# A panel of four same-sized directed networks (synthetic for the demo):
# each wave is the previous one with 20 random dyads toggled
rng = Xoshiro(5)
net_t0 = network(20; directed=true)
for i in 1:20, j in 1:20
    i != j && rand(rng) < 0.12 && add_edge!(net_t0, i, j)
end
networks = [net_t0]
for t in 2:4
    w = copy(networks[end])
    for _ in 1:20
        i, j = rand(rng, 1:20), rand(rng, 1:20)
        i == j && continue
        has_edge(w, i, j) ? rem_edge!(w, i, j) : add_edge!(w, i, j)
    end
    push!(networks, w)
end

formation = [Edges(), Mutual()]
dissolution = [Edges()]

naive = stergm(networks, formation, dissolution; method = :cmple)
naive.inference_withheld                     # true: Mutual is dyad-dependent
result = stergm(networks, formation, dissolution;
                method = :cmple, se = :bootstrap, n_boot = 100, rng = Xoshiro(42))
se_method(result)                            # :bootstrap
confint(result)                              # 3×2
try
    stergm(networks, formation, dissolution; method = :cmple, se = :block_bootstrap)
catch e
    println(e.msg)                           # "cmple: se=:block_bootstrap resamples whole time-transitions, and the panel has only 3 transitions; ..."
end
```

## Interpretation

The dissolution model is fit in tergm's **persistence** (`Persist()`)
parameterisation: `persistence_coef(result)` are persistence log-odds, so
in an edges-only model `persistence_coef(result)[1] = logit(P(tie
persists))` and `formation_coef(result)[1] = logit(P(non-tie forms))` —
both reproduced exactly by the test suite, along with
simulation→estimation round trips. tergm's `Diss()` coefficients are the
negation, `-persistence_coef(result)`; TERGM.jl offers one sign convention
only.

## Inspecting a fit

`show(result)` prints a header — panels, transitions and free dyads, the
(pseudo-)log-likelihoods, how the standard errors were obtained (`inverse
Hessian of the pseudo-likelihood`, `parametric bootstrap (n replicates)`,
…), `Converged: true/false`, and for a CMLE fit one line of MCMC
diagnostics per side — then one R-style coefficient block per model,
`Formation:` and `Persistence:`, rendered through the ecosystem's shared
`NetworkCore.print_coeftable` with the significance-code legend printed once.
An unconverged fit says `WARNING:` in the header, before any table; a
coefficient fixed at `±Inf`, a dyad-dependent formula with naive standard
errors, and excluded bootstrap replicates are noted under the tables.

The printed blocks are the two halves of [`coeftable`](@ref)`(result)`: a
`NetworkCore.CoefficientTable` (`Estimate`, `Std.Error`, `z value`,
`Pr(>|z|)`) with one row per stacked coefficient, labelled exactly as
`summary(tergm)` labels them — `Form(1)~edges`, `Persist(1)~nodematch.grp`
— with the term part being ERGM.jl's direction-aware R label. Rows are read
by index or by name (`tbl["Persist(1)~edges"].p_value`).
[`confint`](@ref)`(result; level = 0.95)` gives Wald limits `θ̂ ± z·se`
from the standard errors the fit reports (`result.se_type` says which),
and throws for a fit whose inference is withheld. The full StatsAPI surface
is defined on
`STERGMResult` — `coef`, `stderror`, `vcov`, `confint`, `loglikelihood`,
`nobs`, `dof`, `aic`, `bic`, `coeftable`, on the stacked (formation, then
persistence) vector — and pinned by
`NetworkCore.check_statsapi(result; strict = true)`; every verb is the one
`StatsAPI` binding NetworkCore.jl and ERGM.jl re-export, so `using ERGM, TERGM`
leaves each of them defined once.

```julia
using TERGM, ERGM, NetworkCore

s50 = load_dataset(:s50)
networks = [copy(w) for w in s50.friendship]
for w in networks, v in 1:nv(w)
    set_vertex_attribute!(w, :smoke, v, s50.smoke[v, 1])
end
fit = stergm(networks, [Edges(), NodeMatch(:smoke)], [Edges(), NodeMatch(:smoke)])

tbl = coeftable(fit)
tbl.names                                   # ["Form(1)~edges", "Form(1)~nodematch.smoke", "Persist(1)~edges", "Persist(1)~nodematch.smoke"]
tbl["Persist(1)~nodematch.smoke"].estimate == persistence_coef(fit)[2]   # true
ci = confint(fit)                           # 4×2, rows in coef(fit)'s order
all(ci[:, 1] .< coef(fit) .< ci[:, 2])      # true
check_statsapi(fit; strict = true)          # every verb present and consistent
```

## Not implemented

- **EGMME** — `TERGM.egmme(model)` and `method = :egmme` throw an
  `ArgumentError`; `egmme` is not exported. There is no equilibrium
  estimation from a cross-section plus durations.
- **CMLE refinements of `tergm`** — `cmle` draws a fixed number of MCMC
  samples per iteration at a fixed thinning interval (enlarged only by the
  stopping rule's boost); `tergm`'s effective-sample-size-adaptive
  sampling is not implemented, and neither are missing-data CMLE or
  offsets. A side whose CMPLE start does not exist is refused.
- **Offsets and constraints** — ERGM.jl's `Offset(term, coef)` is refused
  in a STERGM formula (R `tergm` also refuses offsets inside
  `Form()`/`Persist()` in a CMLE/CMPLE fit; it uses them with EGMME); there
  is no `constraints=`.
- **Terms ERGM.jl does not have** (see its README): a STERGM formula takes
  exactly the terms ERGM.jl provides, plus `Delrecip`.
- **The block bootstrap on short panels** — `se = :block_bootstrap` is refused
  below 10 transitions and warned about below 20; use `se = :bootstrap`.
- **Missing (masked) dyads** — `STERGMModel` refuses a panel with masked
  dyads (`ArgumentError: TERGM (panel t): ...`) and exposes no `missing=`
  keyword (`missing_method(fit) == :rejected`); see [Missing data](@ref).
- **Two-mode (bipartite) panels** — refused at construction
  (`ArgumentError: TERGM (panel t): two-mode (bipartite) panels are not
  supported — the one-mode CMPLE would enumerate the impossible within-mode
  dyads as free dyads …`); bipartite terms and two-mode free-dyad sets are
  not implemented.
- **Self-loops** — a panel containing a loop is refused (`ArgumentError:
  TERGM (panel t): the panel contains 1 self-loop (at vertex v) …`); remove
  it with `rem_edge!(net, v, v)`.
- **Non-separable (btergm-style) TERGMs** — the memory terms
  `EdgeStability`, `PersistentEdge`, `NewEdge` are refused in a formula
  (`ArgumentError: formation model: term 'edge.stability' is constant on
  the free dyads …`) and kept as descriptives; see [Temporal Terms](@ref).
- **Panels of changing composition and `networkDynamic` input** — every
  panel must have the same vertex set; pass a vector of `Network`s.
