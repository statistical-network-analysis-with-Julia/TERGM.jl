# STERGM Estimation

## CMPLE

[`cmple`](@ref) (the default of [`stergm`](@ref)) is conditional maximum
pseudo-likelihood on the Krivitsky–Handcock auxiliary networks:

- **formation rows**: for every transition, the non-edges of ``Y_{t-1}``;
  response = presence in ``Y_t``; change statistics evaluated on
  ``Y^+ = Y_{t-1} \cup Y_t``;
- **dissolution rows**: the edges of ``Y_{t-1}``; response = persistence
  into ``Y_t``; change statistics evaluated on
  ``Y^- = Y_{t-1} \cap Y_t``.

Rows pool across transitions; separability makes the two logistic
likelihoods independent, each maximized by ERGM.jl's MPLE design fitter on
the ecosystem's shared `Networks.newton_fit` optimizer (Newton-Raphson
with step-halving; `maxiter`, `tol`) and `Networks.logistic_derivatives`.
A CMPLE is a logistic pseudo-likelihood over the free dyads, so it fails
in exactly the ways ERGM.jl's MPLE does — and none of them is silent (see
[Degenerate fits](@ref) below).

Model construction validates both formulas against **every** panel with
ERGM.jl's own `ERGM._validate_formula` — the function the `ERGMModel`
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

Both formulas are then **expanded** as an `ERGMModel` expands its terms
(ERGM.jl's `_materialize`): a multi-level `NodeFactor(:grp)` becomes one
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
1e-4 on tergm as shipped.

!!! warning "CMPLE vs CMLE"
    For **dyad-independent** terms CMPLE *is* the conditional MLE. For
    dyad-dependent terms (`Triangle`, `Mutual`, ...) it is an
    approximation, and the pseudo-likelihood standard errors are
    anticonservative (fitted results print a caveat). `method = :cmle`
    raises an `ArgumentError` — MCMC-based CMLE is not implemented, and
    CMPLE is not the same estimator except for dyad-independent formulas,
    so it will not silently stand in for it. `method = :egmme` likewise
    raises: EGMME is not implemented, and no placeholder estimates are
    returned. `egmme` is not exported (an unimplemented estimator should
    not appear in the public API); it remains callable as `TERGM.egmme`
    purely so that the error explains itself.

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
  Persist~nodematch.grp fixed at -Inf …`), and `se = :bootstrap` is
  refused with an `ArgumentError` — the statistic is at its boundary on
  every resample of the transitions, so no refit could be finite.

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
  `cmple` detects the asymptote, warns (`cmple: the MPLE does not exist
  (perfect separation) …`) and returns the fit with `converged == false`
  rather than the point where Newton met its tolerance on the flat
  asymptote.
- **A side with no free dyad on any transition** (no prior tie to persist,
  or no prior non-tie to form) is a `NaN` fit, warned about in those
  words, with `converged == false`.

## Missing data

A panel with masked (unobserved) dyads is refused when the formula is
bound to it: `STERGMModel` (and so `stergm`/`cmple`) throws
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
using TERGM, ERGM, Networks

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

## Block-bootstrap standard errors

`se = :bootstrap` replaces the naive pseudo-likelihood standard errors
with a per-transition block bootstrap (the `btergm` approach): the
time-transitions are resampled with replacement `n_boot` times, the
CMPLE is refit on each resample, and the empirical covariance of the
refitted coefficients is reported. The loop is the ecosystem's one shared
bootstrap, `Networks.bootstrap_cov` (the refits run on every thread; the
resamples are drawn from `rng` up front, so the result is thread-count
independent and reproducible). Point estimates are unchanged;
`vcov(result)` becomes the full joint bootstrap covariance and
`result.boot_replicates` holds the `n_boot × p` refits. A replicate whose
refit has no finite coefficient (a resample that separates a statistic or
leaves one side without free dyads) is excluded from the covariance, warned
about once, and recorded in `approximations(result)`. Needs at least two
transitions (three panels); `rng` seeds the resampling; `se` is validated
by `Networks.check_se`.

```julia
using TERGM, ERGM, Networks, Random

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

result = stergm(networks, formation, dissolution;
                se = :bootstrap, n_boot = 100, rng = Xoshiro(42))
```

## Interpretation

The dissolution model is fit in tergm's **persistence** (`Persist()`)
parameterisation: `persistence_coef(result)` are persistence log-odds, so
in an edges-only model `persistence_coef(result)[1] = logit(P(tie
persists))` and `formation_coef(result)[1] = logit(P(non-tie forms))` —
both reproduced exactly by the test suite, along with
simulation→estimation round trips. tergm's `Diss()` coefficients are the
negation, `-persistence_coef(result)`; TERGM.jl offers one sign convention
only. The pre-0.2 `dissolution_coef`/`dissolution_se` (field and function
spellings) are deprecated aliases of `persistence_coef`/`persistence_se`
and return the same vectors with a deprecation warning.

## Inspecting a fit

`show(result)` prints a header — panels, transitions and free dyads, how
the standard errors were obtained (`inverse Hessian of the
pseudo-likelihood` or `per-transition block bootstrap (n replicates)`),
`Converged: true/false` — then one R-style coefficient block per model,
`Formation:` and `Persistence:`, rendered through the ecosystem's shared
`Networks.print_coeftable` with the significance-code legend printed once.
An unconverged fit says `WARNING:` in the header, before any table; a
coefficient fixed at `±Inf`, a dyad-dependent formula with naive standard
errors, and excluded bootstrap replicates are noted under the tables.

The printed blocks are the two halves of [`coeftable`](@ref)`(result)`: a
`Networks.CoefficientTable` (`Estimate`, `Std.Error`, `z value`,
`Pr(>|z|)`) with one row per stacked coefficient, labelled exactly as
`summary(tergm)` labels them — `Form(1)~edges`, `Persist(1)~nodematch.grp`
— with the term part being ERGM.jl's direction-aware R label. Rows are read
by index or by name (`tbl["Persist(1)~edges"].p_value`).
[`confint`](@ref)`(result; level = 0.95)` gives Wald limits `θ̂ ± z·se`
from the standard errors the fit reports (inverse-Hessian or bootstrap, as
`result.se_type` says). The full StatsAPI surface is defined on
`STERGMResult` — `coef`, `stderror`, `vcov`, `confint`, `loglikelihood`,
`nobs`, `dof`, `aic`, `bic`, `coeftable`, on the stacked (formation, then
persistence) vector — and pinned by
`Networks.check_statsapi(result; strict = true)`; every verb is the one
`StatsAPI` binding Networks.jl and ERGM.jl re-export, so `using ERGM, TERGM`
leaves each of them defined once.

```julia
using TERGM, ERGM, Networks

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

## Deprecated

Both spellings below still work for one release, emit a deprecation
warning and return exactly what the new name returns; they are removed in
the release after 0.2.0.

- `stergm_gof(result; kwargs...)` → [`gof`](@ref)`(result; kwargs...)`: the
  model family has one goodness-of-fit verb (ERGM, ERGMCount, Siena, REM,
  … all answer to `gof`).
- `dissolution_coef(result)` / `dissolution_se(result)`, and the
  `result.dissolution_coef` / `result.dissolution_se` properties →
  [`persistence_coef`](@ref) / [`persistence_se`](@ref).

## Not implemented

- **MCMC-based CMLE** — `cmle(model)` and `method = :cmle` throw an
  `ArgumentError` whose message points at `cmple`.
- **EGMME** — `TERGM.egmme(model)` and `method = :egmme` throw an
  `ErrorException`; `egmme` is not exported.
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
