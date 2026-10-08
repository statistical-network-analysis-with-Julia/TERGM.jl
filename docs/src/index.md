# TERGM.jl

Model tie formation and persistence between successive observations of the same actors. TERGM.jl fits separable temporal ERGMs to binary network panels, with distinct terms and coefficients for new ties and surviving ties.

| First analysis | Learn the model or data | Reference and detail |
|:--|:--|:--|
| [Fit the s50 friendship panel](getting_started.md) | [Choose temporal terms](guide/terms.md) | [Understand conditional estimation](guide/estimation.md) |

!!! note "Supported scope"

    The implemented estimators are conditional MPLE (CMPLE) and Monte-Carlo conditional MLE (CMLE, the default for a dyad-dependent formula — R's `tergm()` has no default and requires `estimate=`; a dyad-independent one is fitted by its CMPLE, which is then exact) for loop-free, one-mode panels with aligned actors. Missing dyads are refused. Equilibrium moment estimation, offsets and arbitrary nonseparable temporal models are not implemented. A CMPLE of dependent terms is an approximation: its default fit prints no z, p or interval; use the parametric bootstrap or the CMLE for inference.

## Installation

```@raw html
<p>Use Julia <strong>1.12 or newer</strong> and the <a href="/getting-started/">shared workspace installation guide</a>. These development packages are not yet registered; the guide prepares the required sibling checkouts and a Julia environment for the examples.</p>
```

## Quick Start

Separate formation and persistence baselines in the three-wave s50 friendship panel:

```julia
using NetworkCore, ERGM, TERGM

waves = load_dataset(:s50).friendship
fit = stergm(waves, [Edges()], [Edges()])
display(fit)
println((formation=formation_coef(fit), persistence=persistence_coef(fit)))
```

The formation coefficient concerns ties absent at the previous wave. The persistence coefficient concerns ties that were present: a larger value means greater persistence. These observations reveal between-wave transitions, not the unobserved order of changes within each survey interval.

## The separable factorization

```math
P(Y_t \mid Y_{t-1}) = P^+(Y^+ \mid Y_{t-1}) \times P^-(Y^- \mid Y_{t-1})
```

with the **formation network** ``Y^+ = Y_{t-1} \cup Y_t`` and the
**dissolution network** ``Y^- = Y_{t-1} \cap Y_t``. Formation statistics
are evaluated on ``Y^+`` (free dyads: prior non-edges), dissolution
statistics on ``Y^-`` (free dyads: prior edges). Dissolution coefficients
measure **persistence**.

## Contents

```@contents
Pages = ["getting_started.md", "guide/terms.md", "guide/estimation.md",
         "guide/simulation.md", "api/types.md", "api/terms.md",
         "api/estimation.md"]
Depth = 2
```

## Citation

If you use TERGM.jl in your work, please cite it using the entry in
[`CITATION.bib`](https://github.com/statistical-network-analysis-with-Julia/TERGM.jl/blob/main/CITATION.bib).
Please also cite R `tergm` and the methods paper the model comes from
(Krivitsky & Handcock 2014, *JRSS-B* 76(1), 29–46) — see
[how to cite the ecosystem](https://statistical-network-analysis-with-julia.github.io/citing/):

```biblatex
@misc{SNWJTERGMJL,
  author = {Santoni, Simone},
  title = {TERGM.jl: Separable Temporal Exponential Random Graph Models in Julia},
  year = {2026},
  url = {https://github.com/statistical-network-analysis-with-Julia/TERGM.jl},
  note = {Homepage: https://statistical-network-analysis-with-Julia.github.io/TERGM.jl; GitHub: https://github.com/statistical-network-analysis-with-Julia}
}
```
