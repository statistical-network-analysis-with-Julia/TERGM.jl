# Simulation

[`simulate_stergm`](@ref) draws one transition:

1. sample ``Y^+`` by Metropolis toggles over the **non-edges** of
   ``Y_{t-1}`` under the formation coefficients (``Y^+`` always contains
   ``Y_{t-1}``);
2. sample ``Y^-`` by toggles over the **edges** of ``Y_{t-1}`` under the
   dissolution coefficients (``Y^-`` is always contained in ``Y_{t-1}``);
3. combine: ``Y_t = (Y^+ \setminus Y_{t-1}) \cup Y^-``.

[`simulate_network_sequence`](@ref) iterates transitions;
`simulate_stergm(result, n_steps)` continues from the last observed
panel. [`stergm_gof`](@ref) compares observed formed/persisted tie counts
per transition against their simulated distributions with two-sided Monte
Carlo p-values.
