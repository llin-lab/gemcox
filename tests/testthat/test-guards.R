## Every numerical guard must warn when it binds (known defect 5).

test_that("shared Breslow warns when max_jump binds", {
  expect_warning(
    b <- shared_breslow(c(1, 2, 3), c(1, 1, 1), list(c(-5, -5, -5)), matrix(1, 3, 1)),
    class = "gemcox_guard_warning", regexp = "max_jump")
  expect_true(all(b$jump <= 10))
})

test_that("per-cluster Breslow warns when max_jump binds", {
  expect_warning(breslow_weighted(c(1, 2, 3), c(1, 1, 1), c(-5, -5, -5), c(1, 1, 1)),
                 class = "gemcox_guard_warning", regexp = "max_jump")
})

test_that("max_cumhaz, eta_clamp and denom_floor warn when they bind", {
  d <- gemcox_simulate(n = 100, p = 2, mu_sep = 0, beta_sep = 0, K_true = 1, seed = 71)
  expect_warning(shared_breslow(d$time, d$status, list(rep(0, 100)), matrix(1, 100, 1),
                                control = gemcox_control(max_cumhaz = 0.5)),
                 class = "gemcox_guard_warning", regexp = "max_cumhaz")
  expect_warning(breslow_weighted(c(1, 2), c(1, 1), c(40, 0), c(1, 1)),
                 class = "gemcox_guard_warning", regexp = "eta_clamp")
  expect_warning(breslow_weighted(c(1, 2), c(1, 1), c(0, 0), c(1e-9, 0)),
                 class = "gemcox_guard_warning", regexp = "denom_floor")
})

test_that("denom_floor does not warn when no weighted event is affected", {
  expect_silent(breslow_weighted(c(1, 2), c(1, 1), c(0, 0), c(1, 0)))
})

test_that("eta_cap warns when a linear predictor is capped", {
  X <- matrix(c(-3, 0, 3), ncol = 1)
  expect_warning(e <- eta_from_scaled(X, 20, make_x_scaler(X)),
                 class = "gemcox_guard_warning", regexp = "eta_cap")
  expect_equal(max(abs(e)), 12)
})

test_that("gemcox() reports guards once, records them, and they do not change the fit", {
  d <- gemcox_simulate(n = 200, p = 3, mu_sep = 1, beta_sep = 1, seed = 72)
  ctl <- gemcox_control(max_jump = 0.01)
  w <- NULL
  fit <- withCallingHandlers(
    gemcox(d$X, time = d$time, status = d$status, control = ctl),
    warning = function(cond) {
      w <<- c(w, list(cond))
      invokeRestart("muffleWarning")
    })
  expect_length(w, 1)
  expect_s3_class(w[[1]], "gemcox_guard_summary")
  expect_true("max_jump" %in% fit$guards$guard)
  quiet <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, control = ctl))
  expect_identical(coef(quiet), coef(fit))
})
