# Golden fixture: one-step STERGM simulation, statnet `tergm` against
# TERGM.jl's `simulate_stergm` (and hence `gof`, which draws from it).
#
# Regenerate from the package root (~3 min):
#
#   Rscript test/fixtures/r/simulate_stergm.R > test/fixtures/simulate_stergm.toml
#
# WHAT IS BEING VALIDATED
#
# `simulate_stergm(prev, formula, theta_form, theta_persist)` draws Y_t from
# the separable model given a fixed Y_{t-1}. R's
# `simulate(nw ~ Form(~...) + Persist(~...), coef = ..., time.slices = 1,
# dynamic = TRUE)` draws from the same conditional distribution, with a
# different sampler written by other people. Both are Monte-Carlo, so what is
# compared is the MEAN of statistics of the simulated transition over many
# draws -- ties formed, ties persisted, the model statistics of Y+ and Y-, and
# the out-/in-degree distribution of Y_t (the panels `gof` reports) -- at a
# tolerance set by the Monte-Carlo error, which the fixture measures: every
# statistic's per-draw sd, and the spread of R's mean across 5 seeds.
#
# TWO MODELS, FROM ONE FROZEN STARTING NETWORK (20 actors, directed)
#
# (a) `mutual_*`: Form(~edges + mutual) + Persist(~edges + mutual). The
#     conditional distribution factorises over unordered pairs, so the EXACT
#     expectation of the pair-additive statistics (formed, persisted, edges
#     and mutual of Y+ and Y-) is computed here by enumeration and frozen as
#     `mutual_exact_means`. tergm's simulated means are asserted (stopifnot)
#     to agree with it; the Julia testset compares against the exact values.
# (b) `tri_*`: the tergm tutorial's Form(~edges + mutual + cyclicalties +
#     transitiveties) + Persist(~the same). No exact answer: the target is
#     tergm's mean over 5 x 1000 draws.
#
# Each time step runs 20 000 MCMC proposals (about 60 per free dyad), far past
# burn-in for both samplers, so the comparison is of stationary distributions.

suppressMessages({
  .libPaths(c(path.expand("~/R/library"), .libPaths()))
  library(network)
  library(ergm)
  library(tergm)
})

seed <- 20261003
set.seed(seed)
n <- 20
A <- matrix(0L, n, n)
A[row(A) != col(A)] <- rbinom(n * (n - 1), 1, 0.12)
nw <- network(A, directed = TRUE)

sim_seeds <- 2000 + 1:5
n_draws <- 1000
burnin <- 20000
max_deg <- 8          # degree bins 0..7, and ">= 8"

theta_mutual <- c(-2.6, 1.8, 0.2, 1.5)
theta_tri <- c(-2.25, 1.85, 0.05, -0.25, 0.0, 1.6, 0.15, 0.1)

deg_hist <- function(d) c(tabulate(pmin(d, max_deg) + 1, nbins = max_deg + 1))

# statistics of one simulated Y_t
draw_stats <- function(B, triadic) {
  yplus <- pmax(A, B); yminus <- A * B
  f <- if (triadic) ~edges + mutual + cyclicalties + transitiveties else ~edges + mutual
  sp <- summary(as.formula(paste("network(yplus, directed = TRUE)", deparse(f))))
  sm <- summary(as.formula(paste("network(yminus, directed = TRUE)", deparse(f))))
  c(sum(B == 1 & A == 0), sum(B == 1 & A == 1), sp, sm,
    deg_hist(rowSums(B)), deg_hist(colSums(B)))
}
stat_names <- function(triadic) {
  t <- if (triadic) c("edges", "mutual", "cyclicalties", "transitiveties") else c("edges", "mutual")
  c("formed", "persisted", paste0("form.", t), paste0("persist.", t),
    paste0("odegree", 0:max_deg), paste0("idegree", 0:max_deg))
}

simulate_model <- function(rhs, theta, triadic) {
  per_seed <- lapply(sim_seeds, function(s) {
    set.seed(s)
    sims <- NULL
    invisible(capture.output(suppressMessages(
      sims <- simulate(as.formula(paste("nw", rhs)), coef = theta, time.slices = 1,
                       dynamic = TRUE, nsim = n_draws, output = "final",
                       control = control.simulate.formula.tergm(
                         MCMC.burnin.min = burnin, MCMC.burnin.max = burnin))),
      type = "output"))
    t(sapply(sims, function(x) draw_stats(as.matrix(x), triadic)))
  })
  all <- do.call(rbind, per_seed)
  list(mean = colMeans(all), sd = apply(all, 2, sd),
       seed_means = t(sapply(per_seed, colMeans)),
       seed_sd = apply(t(sapply(per_seed, colMeans)), 2, sd))
}

sim_mutual <- simulate_model("~ Form(~edges + mutual) + Persist(~edges + mutual)",
                             theta_mutual, FALSE)
sim_tri <- simulate_model(
  "~ Form(~edges + mutual + cyclicalties + transitiveties) + Persist(~edges + mutual + cyclicalties + transitiveties)",
  theta_tri, TRUE)

# --- exact expectations for model (a), by pair enumeration --------------------
exact <- c(formed = 0, persisted = 0, fe = 0, fm = 0, pe = 0, pm = 0)
exact_var <- exact
for (i in 1:(n - 1)) for (j in (i + 1):n) {
  a <- A[i, j]; b <- A[j, i]
  side <- function(theta, formation) {
    us <- if (formation) a:1 else 0:a
    vs <- if (formation) b:1 else 0:b
    g <- expand.grid(u = us, v = vs)
    w <- exp(theta[1] * (g$u + g$v) + theta[2] * (g$u * g$v)); p <- w / sum(w)
    e <- g$u + g$v; m <- g$u * g$v
    chg <- if (formation) e - (a + b) else e          # ties formed / ties persisted
    mom <- function(x) c(sum(p * x), sum(p * x^2) - sum(p * x)^2)
    list(chg = mom(chg), e = mom(e), m = mom(m))
  }
  f <- side(theta_mutual[1:2], TRUE); d <- side(theta_mutual[3:4], FALSE)
  exact <- exact + c(f$chg[1], d$chg[1], f$e[1], f$m[1], d$e[1], d$m[1])
  exact_var <- exact_var + c(f$chg[2], d$chg[2], f$e[2], f$m[2], d$e[2], d$m[2])
}
# tergm's simulation has the exact conditional distribution as its target:
# its 5000-draw means sit within 4 standard errors of the exact expectations
stopifnot(all(abs(sim_mutual$mean[1:6] - exact) <= 4 * sqrt(exact_var / (length(sim_seeds) * n_draws))))
# ... and its seed-to-seed spread is what independent draws would give
stopifnot(all(sim_mutual$seed_sd[1:6] < 3 * sqrt(exact_var / n_draws)))

num <- function(x) paste(sprintf("%.17g", x), collapse = ", ")
strs <- function(x) paste(sprintf('"%s"', x), collapse = ", ")
e <- which(A == 1L, arr.ind = TRUE)

cat('name = "simulate_stergm"\n\n')
cat("[provenance]\n")
cat(sprintf('r_version = "%s"\n', as.character(getRversion())))
cat(sprintf('tergm_version = "%s"\n', as.character(packageVersion("tergm"))))
cat(sprintf('ergm_version = "%s"\n', as.character(packageVersion("ergm"))))
cat(sprintf('network_version = "%s"\n', as.character(packageVersion("network"))))
cat(sprintf("seed = %d\n", seed))
cat('script = "test/fixtures/r/simulate_stergm.R"\n')
cat(sprintf('date = "%s"\n', format(Sys.Date())))
cat('dataset = "simulated: one 20-actor directed starting network, Bernoulli(0.12), frozen below as an edge list"\n')
cat(sprintf('model = "simulate(nw ~ Form(~...) + Persist(~...), coef = theta, time.slices = 1, dynamic = TRUE, nsim = %d, output = \\"final\\") under each of the seeds %s, MCMC.burnin.min = MCMC.burnin.max = %d; (a) edges + mutual, (b) edges + mutual + cyclicalties + transitiveties"\n',
            n_draws, paste(sim_seeds, collapse = ","), burnin))
cat('exact = "model (a): exact expectations and variances of formed, persisted, and edges/mutual of Y+ and Y- by enumeration of the pair states"\n')
cat("\n")
cat("[tolerance]\n")
cat("# Means of simulated statistics are compared at their Monte-Carlo error:\n")
cat("# the Julia testset uses 4 x sd x sqrt(1/N_julia + 1/N_ref) per statistic,\n")
cat("# with `*_sds` the per-draw standard deviations frozen in [values] and\n")
cat("# N_ref = 5000 for tergm's mean (N_ref = Inf for the exact expectations).\n")
cat("# `*_seed_sd` records how much tergm's 1000-draw mean moves between its\n")
cat("# five seeds -- the same quantity, measured: sd / sqrt(1000). No fixed\n")
cat("# number is given here because the tolerance depends on the Julia sample\n")
cat("# size; `mc_sigmas` is the multiplier.\n")
cat("mc_sigmas = 4\n")
cat("\n")
cat("[values]\n")
cat(sprintf("n_actors = %d\n", n))
cat(sprintf("edge_count = %d\n", sum(A)))
cat(sprintf("src = [%s]\n", paste(e[, 1], collapse = ", ")))
cat(sprintf("dst = [%s]\n", paste(e[, 2], collapse = ", ")))
cat(sprintf("n_draws = %d\n", length(sim_seeds) * n_draws))
cat(sprintf("burnin = %d\n", burnin))
cat(sprintf("max_degree_bin = %d\n", max_deg))
emit <- function(prefix, theta, sim, triadic) {
  k <- length(theta) / 2
  cat(sprintf("%s_formation_coefficients = [%s]\n", prefix, num(theta[1:k])))
  cat(sprintf("%s_persistence_coefficients = [%s]\n", prefix, num(theta[(k + 1):(2 * k)])))
  cat(sprintf("%s_statistics = [%s]\n", prefix, strs(stat_names(triadic))))
  cat(sprintf("%s_means = [%s]\n", prefix, num(sim$mean)))
  cat(sprintf("%s_sds = [%s]\n", prefix, num(sim$sd)))
  cat(sprintf("%s_seed_sd = [%s]\n", prefix, num(sim$seed_sd)))
}
cat("# --- (a) edges + mutual ------------------------------------------------------\n")
emit("mutual", theta_mutual, sim_mutual, FALSE)
cat("# exact, by pair enumeration: formed, persisted, form.edges, form.mutual, persist.edges, persist.mutual\n")
cat(sprintf("mutual_exact_means = [%s]\n", num(exact)))
cat(sprintf("mutual_exact_sds = [%s]\n", num(sqrt(exact_var))))
cat("# --- (b) the tutorial's triadic formula -------------------------------------\n")
emit("tri", theta_tri, sim_tri, TRUE)
