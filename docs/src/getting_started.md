# Getting Started

```julia
using TERGM, ERGM, Network

# A panel of same-sized directed networks
networks = [net_t0, net_t1, net_t2]

result = stergm(networks,
                [Edges(), Mutual()],   # formation
                [Edges()])             # dissolution (persistence)

result.formation_coef
result.dissolution_coef   # persistence log-odds

# The auxiliary networks are available directly
yplus  = formation_network(networks[1], networks[2])
yminus = dissolution_network(networks[1], networks[2])

# Simulate 10 steps forward and check fit
future = simulate_stergm(result, 10)
stergm_gof(result; n_sim = 100)
```

Formation/dissolution models mix standard ERGM.jl terms with temporal
terms conditioned on the previous panel:

```julia
stergm(networks, [Edges(), Delrecip()], [Edges(), Mutual()])
```
