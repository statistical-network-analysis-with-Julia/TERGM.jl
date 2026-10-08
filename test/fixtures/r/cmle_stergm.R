# Golden fixture: statnet `tergm` conditional MLE (estimate = "CMLE") of a
# dyad-dependent separable model, with the EXACT conditional MLE beside it.
#
# Regenerate from the package root (~3 min):
#
#   Rscript test/fixtures/r/cmle_stergm.R > test/fixtures/cmle_stergm.toml
#
# WHAT IS BEING VALIDATED
#
# TERGM.jl's `stergm(...; method = :cmle)` maximises the conditional likelihood
# of a separable model by MCMC (Monte-Carlo MLE on the free dyads of Y+ and Y-,
# pooled over transitions). R's `tergm(..., estimate = "CMLE")` is the same
# estimator, written by other people. Both are Monte-Carlo: two runs of either
# package differ by their seeds, so the comparison is made at the resolution
# the estimator has -- R is refitted under 11 seeds and the seed-to-seed
# standard deviation of every coefficient is frozen beside the mean.
#
# WHY THE MODEL IS edges + mutual: AN EXACT ANSWER EXISTS
#
# With `mutual` as the only dyad-dependent term, the conditional likelihood of
# each side factorises over the n(n-1)/2 unordered PAIRS {i, j}: the two ties
# i->j and j->i depend on each other and on nothing else. The normalising
# constant of a pair has at most four states, so the exact conditional
# likelihood, its exact maximiser, its exact Fisher information and its exact
# maximum are computed here by enumeration and Newton's method -- no Monte
# Carlo. They are frozen as `*_exact_*` and are what the Julia testset
# compares to. tergm's 11-seed mean is asserted (in this script, `stopifnot`)
# to agree with that exact maximiser within its own Monte-Carlo error, which
# is what makes "the exact conditional MLE" and "what tergm estimates" the
# same target.
#
# TWO PANELS
#
# (a) `main_*`: the 25-actor, 8-wave directed panel of `panel_stergm.toml`,
#     regenerated here from the same seed and code (its edge counts are
#     asserted). It was drawn from a dyad-INDEPENDENT model, so the mutual
#     coefficients are near zero and the CMLE is close to the CMPLE.
# (b) `recip_*`: a 20-actor, 5-wave directed panel drawn EXACTLY (pair by
#     pair, by enumeration) from Form(~edges + mutual) + Persist(~edges +
#     mutual) with strong reciprocity. Here the CMLE and the CMPLE differ by
#     more than the Monte-Carlo error, so a CMLE that silently returned the
#     CMPLE would fail. The waves are frozen below as edge lists.
#
# For a dyad-independent formula the CMLE is the CMPLE; `main_indep_*` records
# that tergm's CMLE of edges + nodematch returns its CMPLE.
#
# THE TRIADIC MODEL (`tri_*`)
#
# The formula of the statnet tergm tutorial (Goodreau et al.), Form(~edges +
# mutual + cyclicalties + transitiveties) + Persist(~the same), fitted to
# the reciprocity panel. No exact answer exists for a triadic model, so the
# target is tergm itself: its CMLE under the same 11 seeds (mean, seed-to-seed
# sd, mean standard errors), and its CMPLE (deterministic: it pins the change
# statistics of the triadic terms on Y+ and Y- at 1e-6).

suppressMessages({
  .libPaths(c(path.expand("~/R/library"), .libPaths()))
  library(network)
  library(tergm)
})

n_seeds <- 11
cmle_seeds <- 1000 + seq_len(n_seeds)

# --- (a) the main panel of panel_stergm.toml, regenerated ---------------------
seed <- 20260713
set.seed(seed)
n_actors <- 25
n_waves <- 8
grp <- rep(c("a", "b"), length.out = n_actors)
same <- outer(grp, grp, "==")
diag(same) <- FALSE
logistic <- function(x) 1 / (1 + exp(-x))
theta_form <- c(edges = -2.2, nodematch.grp = 0.8)
theta_persist <- c(edges = 0.5, nodematch.grp = 0.9)
P_form <- logistic(theta_form[1] + theta_form[2] * same)
P_persist <- logistic(theta_persist[1] + theta_persist[2] * same)
offdiag <- which(row(matrix(0, n_actors, n_actors)) !=
                 col(matrix(0, n_actors, n_actors)))
M <- matrix(0L, n_actors, n_actors)
M[offdiag] <- rbinom(length(offdiag), 1, 0.10)
mats <- list(M)
for (t in 2:n_waves) {
  U <- matrix(runif(n_actors * n_actors), n_actors, n_actors)
  N <- ifelse(M == 1L, as.integer(U < P_persist), as.integer(U < P_form))
  diag(N) <- 0L
  mats[[t]] <- N
  M <- N
}
# The panel frozen in panel_stergm.toml
stopifnot(identical(as.integer(sapply(mats, sum)),
                    c(68L, 121L, 165L, 183L, 209L, 213L, 218L, 212L)))
as_nw <- function(A, g = NULL) {
  nw <- network(A, directed = TRUE)
  if (!is.null(g)) nw %v% "grp" <- g
  nw
}
nets <- lapply(mats, as_nw, g = grp)

# --- the exact conditional likelihood of edges + mutual, by pair enumeration --
# For a transition prev -> cur and an unordered pair {i, j} with prior state
# (a, b) = (prev[i, j], prev[j, i]):
#   formation:   the pair's state (u, v) in Y+ ranges over u >= a, v >= b;
#   persistence: the pair's state (u, v) in Y- ranges over u <= a, v <= b;
# the statistics are edges = u + v and mutual = u * v. Returns, per side, the
# observed statistics, and functions giving the log-likelihood, the expected
# statistics and their covariance at theta.
pair_table <- function(mats, side) {
  rows <- list()
  for (t in 2:length(mats)) {
    prev <- mats[[t - 1]]; cur <- mats[[t]]
    n <- nrow(prev)
    for (i in 1:(n - 1)) for (j in (i + 1):n) {
      a <- prev[i, j]; b <- prev[j, i]
      if (side == "formation") {
        u <- max(a, cur[i, j]); v <- max(b, cur[j, i])
      } else {
        u <- min(a, cur[i, j]); v <- min(b, cur[j, i])
      }
      rows[[length(rows) + 1]] <- c(a, b, u, v)
    }
  }
  do.call(rbind, rows)
}
exact_side <- function(mats, side) {
  tab <- pair_table(mats, side)
  # collapse to the (a, b) prior types and the observed statistics
  obs <- c(sum(tab[, 3] + tab[, 4]), sum(tab[, 3] * tab[, 4]))
  types <- unique(tab[, 1:2, drop = FALSE])
  counts <- apply(types, 1, function(ab) sum(tab[, 1] == ab[1] & tab[, 2] == ab[2]))
  states <- function(ab) {
    us <- if (side == "formation") ab[1]:1 else 0:ab[1]
    vs <- if (side == "formation") ab[2]:1 else 0:ab[2]
    g <- expand.grid(u = us, v = vs)
    cbind(g$u + g$v, g$u * g$v)
  }
  moments <- function(theta) {
    logZ <- 0; mu <- c(0, 0); V <- matrix(0, 2, 2)
    for (k in seq_len(nrow(types))) {
      S <- states(types[k, ])
      w <- exp(S %*% theta); Z <- sum(w); p <- as.numeric(w / Z)
      m <- colSums(S * p)
      logZ <- logZ + counts[k] * log(Z)
      mu <- mu + counts[k] * m
      V <- V + counts[k] * (crossprod(S * sqrt(p)) - tcrossprod(m))
    }
    list(logZ = logZ, mu = mu, V = V)
  }
  theta <- c(0, 0)
  for (it in 1:200) {                       # Newton on a strictly concave likelihood
    mm <- moments(theta)
    step <- solve(mm$V, obs - mm$mu)
    theta <- theta + step
    if (max(abs(step)) < 1e-13) break
  }
  mm <- moments(theta)
  stopifnot(max(abs(obs - mm$mu)) < 1e-8)
  list(coef = as.numeric(theta), se = sqrt(diag(solve(mm$V))),
       loglik = sum(theta * obs) - mm$logZ, obs = obs, n_pairs = nrow(tab))
}

# --- tergm CMLE under `n_seeds` seeds ------------------------------------------
f_mutual <- ~ Form(~edges + mutual) + Persist(~edges + mutual)
f_triadic <- ~ Form(~edges + mutual + cyclicalties + transitiveties) +
  Persist(~edges + mutual + cyclicalties + transitiveties)
fit_cmle <- function(nets, s, rhs = f_mutual) {
  f <- NULL
  fml <- as.formula(paste("nets", paste(deparse(rhs), collapse = " ")))
  invisible(capture.output(suppressMessages(suppressWarnings(
    f <- tergm(fml, estimate = "CMLE", times = seq_along(nets),
               control = control.tergm(seed = s)))), type = "output"))
  ll <- NA_real_
  invisible(capture.output(suppressMessages(suppressWarnings(
    ll <- as.numeric(logLik(f)))), type = "output"))
  list(coef = as.numeric(coef(f)), se = sqrt(diag(vcov(f))), names = names(coef(f)),
       loglik = ll)
}
run_seeds <- function(nets, rhs = f_mutual) {
  fits <- lapply(cmle_seeds, function(s) fit_cmle(nets, s, rhs))
  coefs <- t(sapply(fits, `[[`, "coef"))
  ses <- t(sapply(fits, `[[`, "se"))
  lls <- sapply(fits, `[[`, "loglik")
  list(names = fits[[1]]$names, mean = colMeans(coefs), sd = apply(coefs, 2, sd),
       se_mean = colMeans(ses), se_sd = apply(ses, 2, sd),
       loglik_mean = mean(lls), loglik_sd = sd(lls))
}
fit_cmple <- function(nets, rhs = f_mutual) {
  f <- NULL
  fml <- as.formula(paste("nets", paste(deparse(rhs), collapse = " ")))
  invisible(capture.output(suppressMessages(
    f <- tergm(fml, estimate = "CMPLE", times = seq_along(nets))), type = "output"))
  list(coef = as.numeric(coef(f)), se = sqrt(diag(vcov(f))),
       loglik = as.numeric(logLik(f)))
}

main_exact_f <- exact_side(mats, "formation")
main_exact_p <- exact_side(mats, "persistence")
main_exact <- c(main_exact_f$coef, main_exact_p$coef)
main_exact_se <- c(main_exact_f$se, main_exact_p$se)
main_cmle <- run_seeds(nets)
main_cmple <- fit_cmple(nets)

# The dyad-independent formula: tergm's CMLE is its CMPLE
indep <- function(est) {
  f <- NULL
  invisible(capture.output(suppressMessages(suppressWarnings(
    f <- tergm(nets ~ Form(~edges + nodematch("grp")) + Persist(~edges + nodematch("grp")),
               estimate = est, times = 1:n_waves,
               control = control.tergm(seed = 1)))), type = "output"))
  as.numeric(coef(f))
}
indep_cmle <- indep("CMLE")
indep_cmple <- indep("CMPLE")

# --- (b) the reciprocity panel: an exact draw from edges + mutual --------------
recip_seed <- 20261002
set.seed(recip_seed)
r_n <- 20
r_waves <- 5
r_theta_form <- c(edges = -2.6, mutual = 1.8)
r_theta_persist <- c(edges = 0.2, mutual = 1.5)
draw_pair <- function(a, b, theta, side) {
  us <- if (side == "formation") a:1 else 0:a
  vs <- if (side == "formation") b:1 else 0:b
  g <- expand.grid(u = us, v = vs)
  w <- exp(theta[1] * (g$u + g$v) + theta[2] * (g$u * g$v))
  k <- sample.int(nrow(g), 1, prob = w)
  c(g$u[k], g$v[k])
}
R0 <- matrix(0L, r_n, r_n)
R0[row(R0) != col(R0)] <- rbinom(r_n * (r_n - 1), 1, 0.12)
r_mats <- list(R0)
for (t in 2:r_waves) {
  prev <- r_mats[[t - 1]]
  N <- matrix(0L, r_n, r_n)
  for (i in 1:(r_n - 1)) for (j in (i + 1):r_n) {
    a <- prev[i, j]; b <- prev[j, i]
    plus <- draw_pair(a, b, r_theta_form, "formation")       # Y+ on the pair
    minus <- draw_pair(a, b, r_theta_persist, "persistence") # Y- on the pair
    # Y_t = (Y+ \ Y_{t-1}) u Y-
    N[i, j] <- if (a == 1L) minus[1] else plus[1]
    N[j, i] <- if (b == 1L) minus[2] else plus[2]
  }
  r_mats[[t]] <- N
}
r_nets <- lapply(r_mats, as_nw)
recip_exact_f <- exact_side(r_mats, "formation")
recip_exact_p <- exact_side(r_mats, "persistence")
recip_exact <- c(recip_exact_f$coef, recip_exact_p$coef)
recip_exact_se <- c(recip_exact_f$se, recip_exact_p$se)
recip_cmle <- run_seeds(r_nets)
recip_cmple <- fit_cmple(r_nets)

# --- (c) the tutorial's triadic formula on the reciprocity panel ---------------
tri_cmle <- run_seeds(r_nets, f_triadic)
tri_cmple <- fit_cmple(r_nets, f_triadic)
stopifnot(length(tri_cmle$mean) == 8)

# tergm's CMLE estimates the exact conditional MLE: its 11-seed mean sits
# within 4 standard errors of the mean (sd / sqrt(11)) of the exact maximiser,
# with a floor of 2% of the exact standard error for the O(1/sample) bias of
# a Monte-Carlo MLE.
agree <- function(cm, ex, ex_se)
  all(abs(cm$mean - ex) <= 4 * cm$sd / sqrt(n_seeds) + 0.02 * ex_se)
stopifnot(agree(main_cmle, main_exact, main_exact_se))
stopifnot(agree(recip_cmle, recip_exact, recip_exact_se))
# ... and on the reciprocity panel the CMPLE is a different number: further
# from the exact CMLE than 4 seed standard deviations on some coefficient
stopifnot(any(abs(recip_cmple$coef - recip_exact) > 4 * recip_cmle$sd))

num <- function(x) paste(sprintf("%.17g", x), collapse = ", ")
strs <- function(x) paste(sprintf('"%s"', x), collapse = ", ")
edge_arrays <- function(ms, sel) paste(sapply(ms, function(A) {
  e <- which(A == 1L, arr.ind = TRUE)
  paste0("[", paste(e[, sel], collapse = ", "), "]")
}), collapse = ", ")

cat('name = "cmle_stergm"\n\n')
cat("[provenance]\n")
cat(sprintf('r_version = "%s"\n', as.character(getRversion())))
cat(sprintf('tergm_version = "%s"\n', as.character(packageVersion("tergm"))))
cat(sprintf('ergm_version = "%s"\n', as.character(packageVersion("ergm"))))
cat(sprintf('network_version = "%s"\n', as.character(packageVersion("network"))))
cat(sprintf("seed = %d\n", seed))
cat(sprintf("recip_seed = %d\n", recip_seed))
cat('script = "test/fixtures/r/cmle_stergm.R"\n')
cat(sprintf('date = "%s"\n', format(Sys.Date())))
cat('dataset = "(a) main: the 25-actor, 8-wave directed panel of panel_stergm.toml, regenerated from the same seed (edge counts asserted); (b) recip: 20 actors, 5 directed waves drawn exactly, pair by pair, from Form(~edges + mutual) + Persist(~edges + mutual) at [values].recip_generating_*, frozen below"\n')
cat('model = "Form(~edges + mutual) + Persist(~edges + mutual), estimate=\\"CMLE\\", control.tergm(seed = s) for each s in cmle_seeds; logLik() by tergm\'s bridge sampler; tri_*: Form(~edges + mutual + cyclicalties + transitiveties) + Persist(~the same) on the recip panel, CMLE under the same seeds and CMPLE"\n')
cat(sprintf('cmle_seeds = "%s"\n', paste(cmle_seeds, collapse = ",")))
cat('exact = "exact conditional MLE by enumeration of the <= 4 states of every unordered pair and Newton\'s method (no Monte Carlo); exact standard errors from the exact Fisher information"\n')
cat("\n")

cat("[tolerance]\n")
cat("# A Monte-Carlo MLE is a random variable: two runs differ by their seeds.\n")
cat("# `*_cmle_seed_sd` in [values] is tergm, refitted under 11 seeds,\n")
cat("# disagreeing with itself. The Julia testset compares ONE seeded TERGM.jl\n")
cat("# fit with the EXACT conditional MLE (`*_exact_coefficients`, computed here\n")
cat("# without Monte Carlo) at 4 x that seed sd, floored at 5% of the exact\n")
cat("# standard error: tergm's own seed sd on some coefficients is so small\n")
cat("# (it adapts its MCMC sample size to an effective-size target, which\n")
cat("# TERGM.jl does not) that 4 x sd would otherwise demand more precision\n")
cat("# than a fixed-sample Monte-Carlo MLE has. 5% of a standard error is far\n")
cat("# below anything that changes an inference, and on the reciprocity panel\n")
cat("# it is far below the CMPLE-to-CMLE gap (`recip_cmple_minus_exact`), so a\n")
cat("# fit that returned the CMPLE fails. Standard errors are compared with the\n")
cat("# exact ones at 10% relative (the Fisher information is estimated from\n")
cat("# the MCMC sample, plus the Monte-Carlo term). The log-likelihood is\n")
cat("# compared with the exact maximum at 4 x tergm's seed sd of its bridge\n")
cat("# estimate, floored at 0.5.\n")
tol_coef <- function(cm, ex_se) pmax(4 * cm$sd, 0.05 * ex_se)
cat(sprintf("main_exact_coefficients = %.4g\n", max(tol_coef(main_cmle, main_exact_se))))
cat(sprintf("recip_exact_coefficients = %.4g\n", max(tol_coef(recip_cmle, recip_exact_se))))
cat(sprintf("main_exact_loglik = %.4g\n", max(4 * main_cmle$loglik_sd, 0.5)))
cat(sprintf("recip_exact_loglik = %.4g\n", max(4 * recip_cmle$loglik_sd, 0.5)))
cat("#\n")
cat("# THE TRIADIC MODEL (`tri_*`). No exact answer exists, so one seeded\n")
cat("# TERGM.jl fit is compared with the MEAN of tergm's 11 fits. Both are\n")
cat("# Monte-Carlo estimates of the same maximiser; the tolerance per\n")
cat("# coefficient (`tri_coefficient_tolerances`) is 4 x tergm's seed sd x\n")
cat("# sqrt(1 + 1/11) (the sd of a single fit minus an 11-fit mean), floored at\n")
cat("# 10% of tergm's mean standard error -- the floor again because tergm's\n")
cat("# adaptive sampling makes its seed sd smaller than a fixed-sample MCMLE's.\n")
cat("# The CMPLE is deterministic and compared at 1e-6, as in panel_stergm.toml.\n")
cat("tri_cmple_coefficients = 1e-6\n")
cat("tri_cmple_std_errors = 1e-4\n")
cat("tri_cmple_loglik = 1e-5\n")
cat("\n")

cat("[values]\n")
emit <- function(prefix, cm, ex, ex_se, exf, exp_, cmple) {
  cat(sprintf("%s_term_names = [%s]\n", prefix, strs(cm$names)))
  cat("# the exact conditional MLE (pair enumeration + Newton): the target\n")
  cat(sprintf("%s_exact_coefficients = [%s]\n", prefix, num(ex)))
  cat(sprintf("%s_exact_std_errors = [%s]\n", prefix, num(ex_se)))
  cat(sprintf("%s_exact_loglik = %.17g\n", prefix, exf$loglik + exp_$loglik))
  cat(sprintf("%s_exact_loglik_formation = %.17g\n", prefix, exf$loglik))
  cat(sprintf("%s_exact_loglik_persistence = %.17g\n", prefix, exp_$loglik))
  cat("# per-coefficient tolerance: max(4 x tergm's seed sd, 5% of the exact SE)\n")
  cat(sprintf("%s_coefficient_tolerances = [%s]\n", prefix,
              num(pmax(4 * cm$sd, 0.05 * ex_se))))
  cat("# observed sufficient statistics (edges, mutual) of Y+ and of Y-, pooled\n")
  cat(sprintf("%s_observed_formation_statistics = [%s]\n", prefix, num(exf$obs)))
  cat(sprintf("%s_observed_persistence_statistics = [%s]\n", prefix, num(exp_$obs)))
  cat("# tergm estimate=\"CMLE\" under 11 seeds: mean, seed-to-seed sd\n")
  cat(sprintf("%s_cmle_coefficients = [%s]\n", prefix, num(cm$mean)))
  cat(sprintf("%s_cmle_seed_sd = [%s]\n", prefix, num(cm$sd)))
  cat(sprintf("%s_cmle_std_errors = [%s]\n", prefix, num(cm$se_mean)))
  cat(sprintf("%s_cmle_std_errors_seed_sd = [%s]\n", prefix, num(cm$se_sd)))
  cat(sprintf("%s_cmle_loglik = %.17g\n", prefix, cm$loglik_mean))
  cat(sprintf("%s_cmle_loglik_seed_sd = %.17g\n", prefix, cm$loglik_sd))
  cat("# tergm estimate=\"CMPLE\" of the same model, and its distance from the exact CMLE\n")
  cat(sprintf("%s_cmple_coefficients = [%s]\n", prefix, num(cmple$coef)))
  cat(sprintf("%s_cmple_minus_exact = [%s]\n", prefix, num(cmple$coef - ex)))
}
cat("# --- (a) the main panel (frozen in panel_stergm.toml) ---------------------\n")
cat(sprintf("main_edge_counts = [%s]\n", paste(sapply(mats, sum), collapse = ", ")))
emit("main", main_cmle, main_exact, main_exact_se, main_exact_f, main_exact_p, main_cmple)
cat("# dyad-independent formula (edges + nodematch): tergm's CMLE is its CMPLE\n")
cat(sprintf("main_indep_cmle_coefficients = [%s]\n", num(indep_cmle)))
cat(sprintf("main_indep_cmple_coefficients = [%s]\n", num(indep_cmple)))
cat(sprintf("main_indep_cmle_vs_cmple_max_abs_diff = %.17g\n", max(abs(indep_cmle - indep_cmple))))
cat("\n# --- (b) the reciprocity panel, frozen -------------------------------------\n")
cat(sprintf("recip_n_actors = %d\n", r_n))
cat(sprintf("recip_n_waves = %d\n", r_waves))
cat(sprintf("recip_edge_counts = [%s]\n", paste(sapply(r_mats, sum), collapse = ", ")))
cat(sprintf("recip_wave_src = [%s]\n", edge_arrays(r_mats, 1)))
cat(sprintf("recip_wave_dst = [%s]\n", edge_arrays(r_mats, 2)))
cat(sprintf("recip_generating_formation = [%s]\n", num(r_theta_form)))
cat(sprintf("recip_generating_persistence = [%s]\n", num(r_theta_persist)))
emit("recip", recip_cmle, recip_exact, recip_exact_se, recip_exact_f, recip_exact_p, recip_cmple)
cat("\n# --- (c) the tutorial's triadic formula on the reciprocity panel -----------\n")
cat(sprintf("tri_term_names = [%s]\n", strs(tri_cmle$names)))
cat("# tergm estimate=\"CMLE\" under 11 seeds: mean, seed-to-seed sd, mean SE\n")
cat(sprintf("tri_cmle_coefficients = [%s]\n", num(tri_cmle$mean)))
cat(sprintf("tri_cmle_seed_sd = [%s]\n", num(tri_cmle$sd)))
cat(sprintf("tri_cmle_std_errors = [%s]\n", num(tri_cmle$se_mean)))
cat(sprintf("tri_cmle_std_errors_seed_sd = [%s]\n", num(tri_cmle$se_sd)))
cat(sprintf("tri_cmle_loglik = %.17g\n", tri_cmle$loglik_mean))
cat(sprintf("tri_cmle_loglik_seed_sd = %.17g\n", tri_cmle$loglik_sd))
cat("# max(4 x seed sd x sqrt(1 + 1/11), 10% of the mean SE)\n")
cat(sprintf("tri_coefficient_tolerances = [%s]\n",
            num(pmax(4 * tri_cmle$sd * sqrt(1 + 1 / n_seeds), 0.10 * tri_cmle$se_mean))))
cat("# tergm estimate=\"CMPLE\" of the same model (deterministic)\n")
cat(sprintf("tri_cmple_coefficients = [%s]\n", num(tri_cmple$coef)))
cat(sprintf("tri_cmple_std_errors = [%s]\n", num(tri_cmple$se)))
cat(sprintf("tri_cmple_loglik = %.17g\n", tri_cmple$loglik))
