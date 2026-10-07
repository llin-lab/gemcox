d <- gemcox_simulate(n = 200, p = 3, mu_sep = 1, beta_sep = 2, seed = 111)
fit <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, n_starts = 2))

test_that("coef returns the p x K coefficient matrix", {
  expect_identical(coef(fit), fit$beta)
  expect_equal(dim(coef(fit)), c(3L, 2L))
})

test_that("print, summary and plot run and return their input", {
  expect_output(expect_identical(print(fit), fit), "GeM-Cox fit: 2 clusters")
  s <- summary(fit)
  expect_s3_class(s, "summary.gemcox")
  expect_equal(s$contrast, fit$beta[, 1] - fit$beta[, 2])
  expect_output(print(s), "C1 - C2")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_identical(plot(fit), fit)
  expect_identical(plot(fit, which = "coef"), fit)
})

test_that("the fit records settings, starts and data", {
  expect_equal(fit$settings$gamma, 1)
  expect_equal(fit$settings$lambda, 0.05)
  expect_false(fit$settings$normalize_gmm_by_dim)
  expect_identical(fit$settings$init, "kmeans")
  expect_equal(nrow(fit$starts), 2)
  expect_equal(fit$starts$final_score[fit$best_start], max(fit$starts$final_score))
  expect_true(fit$x_cox_is_x_gmm)
  expect_equal(sum(fit$pi), 1)
  expect_equal(rowSums(fit$tau), rep(1, 200))
})

test_that("K = 1 fits print and summarise without a contrast", {
  f1 <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, K = 1))
  expect_output(print(f1), "1 cluster")
  expect_null(summary(f1)$contrast)
})

test_that("seed controls initialisation without touching the caller's RNG", {
  set.seed(3)
  before <- .Random.seed
  a <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, init = "random", seed = 4))
  expect_identical(.Random.seed, before)
  b <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, init = "random", seed = 4))
  expect_identical(coef(a), coef(b))
})
