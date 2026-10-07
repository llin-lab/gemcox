## gemcox() must reproduce the engine it was ported from.
## The reference is a verbatim copy of the PATCHED functions
## (fixtures/patched_reference.R). Every argument is passed explicitly on
## both sides, so a change of default in either cannot hide a difference.

ref <- new.env(parent = asNamespace("gemcox"))
sys.source(test_path("fixtures", "patched_reference.R"), envir = ref)

max_abs <- function(a, b) max(abs(unname(as.matrix(a)) - unname(as.matrix(b))))

scenarios <- list(
  list(n = 300, p = 4, mu_sep = 1,   beta_sep = 2, rho = 0,   gamma = 1, normalize = FALSE, seed = 101),
  list(n = 400, p = 6, mu_sep = 0,   beta_sep = 3, rho = 0,   gamma = 1, normalize = FALSE, seed = 102),
  list(n = 250, p = 3, mu_sep = 2,   beta_sep = 1, rho = 0.3, gamma = 1, normalize = FALSE, seed = 103),
  list(n = 300, p = 5, mu_sep = 0.5, beta_sep = 2, rho = 0,   gamma = 0, normalize = FALSE, seed = 104),
  list(n = 300, p = 8, mu_sep = 1,   beta_sep = 2, rho = 0,   gamma = 1, normalize = TRUE,  seed = 105)
)

new_fit <- function(d, s, baseline, min_eff_events, n_starts = 1) {
  suppressWarnings(gemcox(
    X_gmm = d$X, X_cox = d$X, time = d$time, status = d$status, K = 2,
    gamma = s$gamma, lambda = 0.05, alpha = 0, baseline = baseline,
    covariance = "diagonal", normalize_gmm_by_dim = s$normalize, temp = 1,
    n_starts = n_starts, max_iter = 100, tol = 1e-4, init = "kmeans", seed = 1,
    verbose = FALSE,
    control = gemcox_control(gmm_ridge = 1e-6, pi_floor = 0.02,
                             min_eff_events = min_eff_events, max_jump = 10,
                             max_cumhaz = 500, denom_floor = 1e-8, eta_cap = 12,
                             eta_clamp = 30)))
}

expect_same_fit <- function(new, old, label) {
  expect_equal(length(new$loglik), length(old$loglik), label = paste(label, "trace length"))
  expect_lt(max_abs(new$loglik, old$loglik), 1e-6, label = paste(label, "loglik trace"))
  expect_lt(max_abs(new$beta, old$beta), 1e-6, label = paste(label, "beta"))
  expect_lt(max_abs(new$tau, old$tau), 1e-6, label = paste(label, "tau"))
  expect_lt(max_abs(new$pi, old$pi), 1e-6, label = paste(label, "pi"))
  expect_lt(max_abs(new$mu, old$mu_gmm), 1e-6, label = paste(label, "mu"))
}

test_that("shared baseline reproduces gemcox_full_shared_v4()", {
  for (s in scenarios) {
    d <- gemcox_simulate(s$n, s$p, s$mu_sep, s$beta_sep, rho = s$rho, seed = s$seed)
    old <- ref$gemcox_full_shared_v4(
      X_gmm = d$X, X_cox = d$X, time = d$time, status = d$status, K = 2,
      lambda = 0.05, alpha = 0, max_iter = 100, tol = 1e-4, verbose = FALSE,
      gmm_ridge = 1e-6, gmm_diag = TRUE, normalize_gmm_by_dim = s$normalize,
      surv_weight = s$gamma, temp = 1, pi_floor = 0.02, min_eff_events = 1,
      init_method = "kmeans", init_eta_weight = 0.75, init_seed = 1)
    expect_same_fit(new_fit(d, s, "shared", 1), old, paste("shared, seed", s$seed))
  }
})

test_that("cluster baseline reproduces gemcox_full()", {
  for (s in scenarios[1:3]) {
    d <- gemcox_simulate(s$n, s$p, s$mu_sep, s$beta_sep, rho = s$rho, seed = s$seed)
    old <- ref$gemcox_full(
      X_gmm = d$X, X_cox = d$X, time = d$time, status = d$status, K = 2,
      lambda = 0.05, alpha = 0, max_iter = 100, tol = 1e-4, verbose = FALSE,
      gmm_ridge = 1e-6, gmm_diag = TRUE, normalize_gmm_by_dim = s$normalize,
      surv_weight = s$gamma, temp = 1, pi_floor = 0.02, min_eff_events = 0.5,
      init_method = "kmeans", init_eta_weight = 0.75, init_seed = 1)
    expect_same_fit(new_fit(d, s, "cluster", 0.5), old, paste("cluster, seed", s$seed))
  }
})

test_that("n_starts > 1 reproduces gemcox_full_multistart_shared_v4()", {
  s <- scenarios[[1]]
  d <- gemcox_simulate(s$n, s$p, s$mu_sep, s$beta_sep, seed = s$seed)
  old <- ref$gemcox_full_multistart_shared_v4(
    X_gmm = d$X, X_cox = d$X, time = d$time, status = d$status, K = 2,
    lambda = 0.05, alpha = 0, max_iter = 100, tol = 1e-4,
    gmm_ridge = 1e-6, gmm_diag = TRUE, normalize_gmm_by_dim = FALSE,
    surv_weight = 1, temp = 1, pi_floor = 0.02, min_eff_events = 1,
    init_method = "kmeans", init_eta_weight = 0.75,
    n_starts = 3, init_seeds = 1:3, verbose = FALSE)
  new <- new_fit(d, s, "shared", 1, n_starts = 3)
  expect_equal(new$best_start, old$multistart$best_start)
  expect_lt(max_abs(new$starts$final_score, old$multistart$final_loglik), 1e-6)
  expect_same_fit(new, old, "multistart")
})
