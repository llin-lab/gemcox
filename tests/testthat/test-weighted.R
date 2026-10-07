## The tau-weighted Cox M-step must be a weighted Cox fit.

wdata <- function(seed) {
  d <- gemcox_simulate(n = 300, p = 4, mu_sep = 1, beta_sep = 1, seed = seed)
  d$w <- with_local_seed(seed, stats::runif(300, 0.05, 1))
  d$Xs <- apply_x_scaler(d$X, make_x_scaler(d$X))
  d
}
no_caps <- gemcox_control(max_jump = 1e12, max_cumhaz = 1e12)

## The engine calls glmnet at its default convergence threshold (as PATCHED
## did), which limits agreement with coxph to ~2e-4. At a tight threshold the
## same weighted objective agrees to < 1e-6, showing the objectives coincide.
test_that("per-cluster M-step (glmnet Cox, lambda = 0) equals coxph(weights = w)", {
  for (s in 41:43) {
    d <- wdata(s)
    cf <- cox_glmnet_weighted(d$time, d$status, d$Xs, d$w, lambda = 0, alpha = 0,
                              min_eff_events = 0, control = no_caps)
    ref <- as.numeric(coef(coxph(Surv(d$time, d$status) ~ d$Xs, weights = d$w,
                                 ties = "breslow")))
    expect_lt(max(abs(cf$beta_s - ref)), 5e-4)
    tight <- glmnet::glmnet(d$Xs, Surv(d$time, d$status), family = "cox", weights = d$w,
                            alpha = 0, lambda = 0, standardize = FALSE,
                            thresh = 1e-14, maxit = 1e6)
    expect_lt(max(abs(as.numeric(coef(tight, s = 0)) - ref)), 1e-6)
  }
})

test_that("weighted Breslow at the weighted Cox solution equals the direct formula", {
  d <- wdata(44)
  ref <- as.numeric(coef(coxph(Surv(d$time, d$status) ~ d$Xs, weights = d$w,
                               ties = "breslow")))
  eta <- as.numeric(d$Xs %*% ref)
  b <- breslow_weighted(d$time, d$status, eta, d$w, control = no_caps)
  direct <- vapply(b$time, function(t) {
    sum(d$w[d$time == t & d$status == 1]) / sum(d$w[d$time >= t] * exp(eta[d$time >= t]))
  }, 0)
  expect_lt(max(abs(b$jump - direct)), 1e-12)
})

test_that("shared M-step (weighted Poisson working model) has the weighted Cox solution as fixed point", {
  for (s in 45:47) {
    d <- wdata(s)
    ref <- as.numeric(coef(coxph(Surv(d$time, d$status) ~ d$Xs, weights = d$w,
                                 ties = "breslow")))
    base <- breslow_weighted(d$time, d$status, as.numeric(d$Xs %*% ref), d$w,
                             control = no_caps)
    pf <- cox_poisson_weighted(d$time, d$status, d$Xs, d$w, base, lambda = 0,
                               alpha = 0, min_eff_events = 0)
    expect_true(pf$ok)
    expect_lt(max(abs(pf$beta_s - ref)), 1e-4)
  }
})
