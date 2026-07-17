# TERGM.jl

Separable Temporal ERGMs (STERGM) for discrete-time network panels,
following Krivitsky & Handcock (2014) and R `tergm`.

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
[`CITATION.bib`](https://github.com/statistical-network-analysis-with-Julia/TERGM.jl/blob/main/CITATION.bib):

```biblatex
@misc{SNWJTERGMJL,
  author = {{Statistical Network Analysis with Julia}},
  title = {TERGM.jl: Separable Temporal Exponential Random Graph Models in Julia},
  year = {2026},
  url = {https://github.com/statistical-network-analysis-with-Julia/TERGM.jl},
  note = {Homepage: https://statistical-network-analysis-with-Julia.github.io/TERGM.jl; GitHub: https://github.com/statistical-network-analysis-with-Julia}
}
```
