# Getting Started

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
gof(result; n_sim = 100)   # Networks.gof generic; stergm_gof is an alias
```

Formation/dissolution models mix standard ERGM.jl terms with temporal
terms conditioned on the previous panel:

```julia
stergm(networks, [Edges(), Delrecip()], [Edges(), Mutual()])
```
