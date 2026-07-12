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
likelihoods independent, each maximized with the shared
`ERGM.newton_fit` optimizer (Newton-Raphson with step-halving).

Model construction validates the formulas: attribute-based terms whose
vertex attribute is missing from a panel, and direction-incompatible
terms (`Mutual`, `Delrecip`) on undirected panels, throw an
`ArgumentError` listing the available attributes.

!!! warning "CMPLE vs CMLE"
    For **dyad-independent** terms CMPLE *is* the conditional MLE. For
    dyad-dependent terms (`Triangle`, `Mutual`, ...) it is an
    approximation, and the pseudo-likelihood standard errors are
    anticonservative (fitted results print a caveat). `method = :cmle`
    warns and falls back to CMPLE — MCMC-based CMLE is not implemented.
    `method = :egmme` raises an error: EGMME is not implemented, and no
    placeholder estimates are returned.

## Block-bootstrap standard errors

`se = :bootstrap` replaces the naive pseudo-likelihood standard errors
with a per-transition block bootstrap (the `btergm` approach): the
time-transitions are resampled with replacement `n_boot` times, the
CMPLE is refit on each resample, and the empirical covariance of the
refitted coefficients is reported. Point estimates are unchanged;
`vcov(result)` becomes the full joint bootstrap covariance. Needs at
least two transitions (three panels); `rng` seeds the resampling.

```julia
using TERGM, ERGM, Network, Random

# A panel of same-sized directed networks (synthetic for the demo)
rng = Xoshiro(5)
net_t0 = network(20; directed=true)
for i in 1:20, j in 1:20
    i != j && rand(rng) < 0.1 && add_edge!(net_t0, i, j)
end
net_t1 = copy(net_t0); net_t2 = copy(net_t1)
for w in (net_t1, net_t2), _ in 1:15
    i, j = rand(rng, 1:20), rand(rng, 1:20)
    i == j && continue
    has_edge(w, i, j) ? rem_edge!(w, i, j) : add_edge!(w, i, j)
end
networks = [net_t0, net_t1, net_t2]

formation = [Edges(), Mutual()]
dissolution = [Edges()]

result = stergm(networks, formation, dissolution;
                se = :bootstrap, n_boot = 100, rng = Xoshiro(42))
```

## Interpretation

Dissolution coefficients are **persistence** log-odds: in an edges-only
model, `dissolution_coef[1] = logit(P(tie persists))` and
`formation_coef[1] = logit(P(non-tie forms))` — both reproduced exactly
by the test suite, along with simulation→estimation round trips.
