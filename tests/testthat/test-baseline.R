no_caps <- gemcox_control(max_jump = 1e12, max_cumhaz = 1e12)

test_that("get_cumhaz_at is 0 before the grid and the step value after", {
  got <- get_cumhaz_at(list(time = c(2, 4, 6), cumhaz = c(.1, .5, .9)), c(1, 3, 5, 7))
  expect_identical(got, c(0, .1, .5, .9))
  expect_length(got, 4)
})

test_that("get_cumhaz_at keeps order and length with pre-grid times anywhere", {
  b <- list(time = c(2, 4, 6), cumhaz = c(.1, .5, .9))
  tt <- c(7, 1, 0.5, 3, 1.5, 6, 2)
  expect_identical(get_cumhaz_at(b, tt), c(.9, 0, 0, .1, 0, .9, .1))
  expect_identical(get_cumhaz_at(list(time = numeric(0), cumhaz = numeric(0)), 1:3),
                   c(0, 0, 0))
  ## The legacy expression gets this wrong: it is the defect being guarded.
  legacy <- function(baseline, times) {
    idx <- findInterval(times, baseline$time)
    ifelse(idx > 0, baseline$cumhaz[idx], 0)
  }
  expect_false(identical(legacy(b, tt), get_cumhaz_at(b, tt)))
})

test_that("shared Breslow at K = 1 equals survival::basehaz at event times", {
  d <- gemcox_simulate(n = 300, p = 4, mu_sep = 0, beta_sep = 0, K_true = 1, seed = 21)
  cf <- coxph(Surv(d$time, d$status) ~ d$X, ties = "breslow")
  eta <- as.numeric(d$X %*% coef(cf))
  b <- shared_breslow(d$time, d$status, list(eta), matrix(1, 300, 1), control = no_caps)
  bh <- survival::basehaz(cf, centered = FALSE)
  ref <- bh$hazard[match(b$time, bh$time)]
  expect_false(anyNA(ref))
  expect_lt(max(abs(b$cumhaz - ref)), 1e-8)
})

test_that("per-cluster weighted Breslow with unit weights equals basehaz", {
  d <- gemcox_simulate(n = 300, p = 4, mu_sep = 0, beta_sep = 0, K_true = 1, seed = 22)
  cf <- coxph(Surv(d$time, d$status) ~ d$X, ties = "breslow")
  eta <- as.numeric(d$X %*% coef(cf))
  b <- breslow_weighted(d$time, d$status, eta, rep(1, 300), control = no_caps)
  bh <- survival::basehaz(cf, centered = FALSE)
  expect_lt(max(abs(b$cumhaz - bh$hazard[match(b$time, bh$time)])), 1e-8)
})
