# Estimation and Simulation

```@meta
CurrentModule = TERGM
```

`fit_stergm` is an exported alias for [`stergm`](@ref). The StatsAPI verbs
`coef`, `stderror`, `vcov`, `loglikelihood`, `nobs`, `dof`, `aic` and `bic`
are defined on [`STERGMResult`](@ref) (stacked formation-then-persistence
vector) alongside the two documented below; `Networks.check_statsapi(fit;
strict = true)` pins the full surface.

```@docs
stergm
cmple
cmle
egmme
formation_coef
formation_se
persistence_coef
persistence_se
dissolution_coef
dissolution_se
has_dyad_dependent(::STERGMModel)
coeftable(::STERGMResult)
confint(::STERGMResult)
objective(::STERGMResult)
is_exact(::STERGMResult)
simulate_stergm
simulate_network_sequence
gof(::STERGMResult)
stergm_gof
```
