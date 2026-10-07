## Version 0.2.0: EM runs to convergence by default; the old setting is
## available through options() (the simulation pipeline uses it). The suite
## sets the old setting in helper-em-setting.R, so each test here clears it.
clear_em_options <- function(env = parent.frame()) {
  old <- options(gemcox.tol = NULL, gemcox.max_iter = NULL)
  withr::defer(options(old), envir = env)
}

test_that("defaults are tol = 1e-8 and max_iter = 1000, overridable by options", {
  clear_em_options()
  d <- gemcox_simulate(n = 200, p = 4, mu_sep = 1, beta_sep = 2, seed = 11)
  f <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, K = 2, seed = 1))
  expect_equal(f$settings$tol, 1e-8)
  expect_equal(f$settings$max_iter, 1000L)
  options(gemcox.tol = 1e-4, gemcox.max_iter = 100L)
  g <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, K = 2, seed = 1))
  expect_equal(g$settings$tol, 1e-4)
  expect_equal(g$settings$max_iter, 100L)
  ## the old setting is what the explicit arguments give
  h <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, K = 2, seed = 1,
                               tol = 1e-4, max_iter = 100))
  expect_identical(g$beta, h$beta)
  expect_lte(g$iterations, f$iterations)
})

test_that("a fit stopped at max_iter warns with class gemcox_not_converged", {
  d <- gemcox_simulate(n = 200, p = 4, mu_sep = 1, beta_sep = 2, seed = 12)
  expect_warning(gemcox(d$X, time = d$time, status = d$status, K = 2, seed = 1, max_iter = 2),
                 class = "gemcox_not_converged")
  f <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, K = 2, seed = 1, max_iter = 2))
  expect_false(f$converged)
  expect_warning(gemcox:::lcc_fit(d$X, time = d$time, status = d$status, seed = 1, max_iter = 2),
                 class = "gemcox_not_converged")
})

test_that("the test and the CV criterion summarise non-convergence once", {
  d <- gemcox_simulate(n = 150, p = 3, mu_sep = 0, beta_sep = 2, seed = 13)
  w <- character(0)
  ht <- withCallingHandlers(
    gemcox_heterogeneity_test(X_gmm = d$X, time = d$time, status = d$status, B = 3, seed = 2,
                              max_iter = 2),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("^gemcox: EM stopped", w)))       # no per-fit warnings
  expect_true(any(grepl("stopped at max_iter without converging", w)))
  expect_gte(ht$n_not_converged, 4)
  w2 <- character(0)
  withCallingHandlers(
    gemcox_cv_loglik(d$X, time = d$time, status = d$status, nfolds = 3, seed = 1, max_iter = 2),
    warning = function(x) { w2 <<- c(w2, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("^gemcox: EM stopped", w2)))
  expect_true(any(grepl("fold fits stopped at max_iter", w2)))
})

## At K = 1 the fitted Cox model does not depend on gamma, but the stopping
## rule watches the log-likelihood trace, which does. With gamma = 0 the trace
## holds only the Gaussian-mixture part, so EM can stop before the Cox block
## updates have fully settled: coefficients differ by about 1e-6 and
## iteration counts may differ. At the old setting every gamma stopped at
## iteration 2 with identical coefficients (test-k1.R, which runs at that
## setting).
test_that("K = 1 is invariant to gamma at the converged default, to 1e-5", {
  clear_em_options()
  d <- gemcox_simulate(n = 400, p = 5, mu_sep = 0, beta_sep = 0, K_true = 1, seed = 31)
  b <- lapply(c(0, 0.5, 1, 2), function(g)
    coef(suppressWarnings(gemcox(d$X, time = d$time, status = d$status, K = 1, gamma = g))))
  for (x in b[-1]) expect_lt(max(abs(x - b[[1]])), 1e-5)
})
