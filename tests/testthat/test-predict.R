d <- gemcox_simulate(n = 300, p = 4, mu_sep = 2, beta_sep = 1, seed = 51)
fit0 <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, gamma = 0))
fit1 <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, gamma = 1))

test_that("tau recomputed on training X_gmm equals fit$tau when gamma = 0", {
  expect_lt(max(abs(predict(fit0, d$X, type = "tau") - fit0$tau)), 1e-6)
  for (s in 52:54) {
    dd <- gemcox_simulate(n = 250, p = 5, mu_sep = 1, beta_sep = 2, seed = s)
    f <- suppressWarnings(gemcox(dd$X, time = dd$time, status = dd$status, gamma = 0,
                                 n_starts = 2))
    expect_lt(max(abs(predict(f, dd$X) - f$tau)), 1e-6)
  }
})

test_that("with gamma = 1 the training tau uses the outcome and prediction does not", {
  expect_gt(max(abs(predict(fit1, d$X) - fit1$tau)), 1e-3)
})

test_that("predict matches columns by name", {
  X_shuffled <- d$X[, c(3, 1, 4, 2)]
  expect_equal(predict(fit1, X_shuffled), predict(fit1, d$X))
  expect_equal(predict(fit1, X_shuffled, type = "lp"), predict(fit1, d$X, type = "lp"))
  X_bad <- d$X
  colnames(X_bad)[1] <- "zzz"
  expect_error(predict(fit1, X_bad), "missing column")
})

test_that("lp is the centred mixture predictor; the uncentred one differs by subject", {
  tau <- predict(fit1, d$X)
  xc <- sweep(d$X, 2, fit1$cox_scaler$center)
  expect_equal(predict(fit1, d$X, type = "lp"), rowSums(tau * (xc %*% coef(fit1))),
               tolerance = 1e-10)
  uncentred <- rowSums(tau * (d$X %*% coef(fit1)))
  gap <- uncentred - predict(fit1, d$X, type = "lp")
  expect_gt(stats::sd(gap), 0)   # not a constant shift when tau varies
})

test_that("survival predictions are probabilities, non-increasing in time", {
  tt <- c(0.01, stats::quantile(d$time, c(.25, .5, .75)), 1e5)
  S <- predict(fit1, d$X[1:20, ], type = "survival", times = tt)
  expect_equal(dim(S), c(20L, 5L))
  expect_true(all(S >= 0 & S <= 1))
  expect_true(all(apply(S, 1, function(s) all(diff(s) <= 1e-12))))
  expect_equal(unname(S[, 1]), rep(1, 20), tolerance = 1e-12)   # before the first event
})

test_that("predict gives informative errors", {
  expect_error(predict(fit1), "newX_gmm is required")
  expect_error(predict(fit1, d$X, type = "survival"), "times is required")
  expect_error(predict(fit1, d$X[, 1:3]), "3 columns")
  f2 <- suppressWarnings(gemcox(d$X[, 1:2], d$X, time = d$time, status = d$status))
  expect_error(predict(f2, d$X[, 1:2], type = "lp"), "newX_cox is required")
  expect_error(predict(f2, d$X[1:5, 1:2], d$X[1:4, ], type = "lp"), "same number of rows")
})
