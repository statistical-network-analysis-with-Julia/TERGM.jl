# Simulation

[`simulate_stergm`](@ref) draws one transition:

1. sample ``Y^+`` by Metropolis toggles over the **non-edges** of
   ``Y_{t-1}`` under the formation coefficients (``Y^+`` always contains
   ``Y_{t-1}``);
2. sample ``Y^-`` by toggles over the **edges** of ``Y_{t-1}`` under the
   dissolution coefficients (``Y^-`` is always contained in ``Y_{t-1}``);
3. combine: ``Y_t = (Y^+ \setminus Y_{t-1}) \cup Y^-``.

Both samplers run on ERGM.jl's exported `mh_toggle!` Metropolis kernel
(one toggle loop for the whole family) and draw only from `rng`. The
`burnin` keyword (toggles per side) defaults to `nothing`, which resolves
to ERGM.jl's dyad-scaled rule `ERGM._mcmc_defaults`: 20 toggles per free
dyad of that side.

[`simulate_network_sequence`](@ref) iterates transitions;
`simulate_stergm(result, n_steps)` continues from the last observed
panel — **one chain, `n_steps` steps**: R's `simulate(fit, nsim = 1,
time.slices = n_steps)`. (tergm's `nsim` is the number of *replications*;
`k` independent one-step draws are `[simulate_stergm(result, 1;
rng = rng)[end] for _ in 1:k]`.)

The formula handed to `simulate_stergm(prev, formula, θ_form, θ_diss)` is
validated against `prev` (a nodal term's attribute must exist there) and
expanded exactly as [`STERGMModel`](@ref) expands it, so a raw
`NodeFactor(:grp)` takes one coefficient per level; the coefficient-count
check names the side's terms and the counts it saw, and points at
`formation_coef(fit)` / `persistence_coef(fit)` (not the stacked
`coef(fit)`).

## Goodness of fit

[`gof`](@ref) — TERGM.jl's method on the shared `Networks.gof` generic,
returning a `Networks.GOFResult` ([`stergm_gof`](@ref) is a deprecated
alias that warns and forwards) — simulates `n_sim` transitions from every
observed Y_{t−1} at the fitted coefficients and compares, pooled over
transitions, with two-sided Monte Carlo p-values. After `gof.tergm`, it
reports several panels, and it matters which one you read:

- **`tie changes`** (`formed`, `persisted`): the pooled counts of ties that
  formed and ties that persisted. For a formula with `Edges()` on that side
  these are the sufficient statistics of `Form~edges` / `Persist~edges`,
  which the CMPLE reproduces in expectation *by construction* — p-values
  near 1 here are expected and say nothing about fit. The panel is a
  simulator sanity check, and the only one that speaks for an edges-free
  formula.
- **`formation statistics`** / **`persistence statistics`**: every term of
  the formation model evaluated on the observed formation network
  ``Y^+ = Y_{t-1} \cup Y_t`` against its value on the simulated ``Y^+``,
  and the same for the persistence model on ``Y^- = Y_{t-1} \cap Y_t``
  (labels are the fitted coefficients'). `edges` is again reproduced by
  construction; the *other* rows are the model-statistic check — the first
  place a **dyad-dependent CMPLE** (a `Mutual`, `Triangle`, `GWESP`, …
  whose pseudo-likelihood point estimate may be biased) shows a misfit.
- **`idegree`** / **`odegree`** (directed) or **`degree`** (undirected):
  the degree distribution of the simulated ``Y_t`` against the observed
  ``Y_t``, levels `0` … `n−1`, pooled over transitions — the out-of-model
  check that is informative for every formula (tergm's default `gof`
  statistics are the degree and duration distributions).

`gof` draws one seed per (transition, simulation) from `rng` up front and
runs the simulations on every available thread, so `gof(result; n_sim =
100)` is as fast as the machine allows and gives the same answer whatever
the thread count; its `burnin` is forwarded to `simulate_stergm`.

```julia
using TERGM, ERGM, Networks, Random
s50 = load_dataset(:s50)
networks = [copy(w) for w in s50.friendship]
fit = stergm(networks, [Edges(), Mutual()], [Edges()])
g = gof(fit; n_sim = 100, rng = Xoshiro(1))
[s.name for s in g.statistics]      # ["tie changes", "formation statistics", "persistence statistics", "idegree", "odegree"]
g.statistics[2].labels              # ["edges", "mutual"] — the mutual row is the one to read
g.statistics[4]                     # in-degree distribution, observed vs simulated envelope
```
