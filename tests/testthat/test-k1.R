## K = 1 must reduce to an ordinary Cox model. This is the assertion that
## caught the coefficient-attenuation bug. The comparison with coxph is at
## lambda = 0; at the default lambda = 0.05 the reference is glmnet's ridge
## Cox with the same penalty (the default shrinks by roughly 10%).

k1_sets <- list(
  gemcox_simulate(n = 400, p = 5, mu_sep = 0, beta_sep = 0, K_true = 1, seed = 31),
  gemcox_simulate(n = 600, p = 8, mu_sep = 0, beta_sep = 0, K_true = 1, seed = 32),
  gemcox_simulate(n = 300, p = 3, mu_sep = 0, beta_sep = 0, K_true = 1, rho = 0.4, seed = 33)
)

fit_k1 <- function(d, ...) {
  suppressWarnings(gemcox(d$X, time = d$time, status = d$status, K = 1, ...))
}

for (bl in c("shared", "cluster")) {
  test_that(sprintf("K = 1 reproduces coxph at lambda = 0 (%s baseline)", bl), {
    for (d in k1_sets) {
      ref <- as.numeric(coef(coxph(Surv(d$time, d$status) ~ d$X, ties = "breslow")))
      b <- as.numeric(coef(fit_k1(d, lambda = 0, baseline = bl)))
      expect_lt(abs(sqrt(sum(b^2)) / sqrt(sum(ref^2)) - 1), 0.02)
      expect_lt(max(abs(b - ref)), 0.02)
    }
  })

  test_that(sprintf("K = 1 is invariant to gamma (%s baseline)", bl), {
    d <- k1_sets[[1]]
    fits <- lapply(c(0, 0.5, 1, 2), function(g) fit_k1(d, gamma = g, baseline = bl))
    for (f in fits[-1]) {
      expect_equal(coef(f), coef(fits[[1]]), tolerance = 1e-10)
      expect_equal(f$iterations, fits[[1]]$iterations)
    }
  })
}

test_that("K = 1 at the default penalty equals glmnet ridge Cox with the same lambda", {
  for (d in k1_sets) {
    sc <- make_x_scaler(d$X)
    g <- glmnet::glmnet(apply_x_scaler(d$X, sc), Surv(d$time, d$status), family = "cox",
                        alpha = 0, lambda = 0.05, standardize = FALSE,
                        thresh = 1e-12, maxit = 1e6)
    ref <- as.numeric(coef(g, s = 0.05)) / sc$scale
    b <- as.numeric(coef(fit_k1(d)))
    expect_lt(abs(sqrt(sum(b^2)) / sqrt(sum(ref^2)) - 1), 0.02)
    expect_lt(max(abs(b - ref)), 1e-4)
  }
})
