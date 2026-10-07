## Latent-class Cox competitor (lcc_fit), added after sim-freeze-v2.

unit <- function(x) x / sqrt(sum(x^2))

test_that("K = 1 reproduces coxph at lambda = 0 (shared and cluster baselines)", {
  d <- gemcox_simulate(n = 400, p = 5, mu_sep = 0, beta_sep = 0, K_true = 1, seed = 31)
  ref <- as.numeric(coef(coxph(Surv(d$time, d$status) ~ d$X, ties = "breslow")))
  for (bl in c("shared", "cluster")) {
    b <- as.numeric(suppressWarnings(lcc_fit(d$X, time = d$time, status = d$status, K = 1,
                                             lambda = 0, baseline = bl))$beta[, 1])
    expect_lt(abs(sqrt(sum(b^2)) / sqrt(sum(ref^2)) - 1), 0.02)
    expect_lt(max(abs(b - ref)), 0.02)
  }
})

test_that("E-step tau matches a hand computation, incl. subjects before the first event", {
  x <- matrix(c(-1, 0.5, 2, -0.3), ncol = 1)
  time <- c(0.5, 1, 2, 3)            # subject 1 is censored before the first event time
  status <- c(0, 1, 0, 1)
  b <- list(time = c(1, 3), jump = c(0.2, 0.5), cumhaz = c(0.2, 0.7))
  gate <- list(a0 = c(0, 0.3), A = matrix(c(0, -0.8), nrow = 1))
  eta <- cbind(c(0.1, -0.2, 0.4, 0.3), c(-0.3, 0.5, 0, -0.1))
  ## by hand: gating softmax; H = (0, 0.2, 0.2, 0.7); dLambda at the events = (-, 0.2, -, 0.5)
  lin <- cbind(0 + 0 * x[, 1], 0.3 - 0.8 * x[, 1])
  lp <- lin - log(rowSums(exp(lin)))
  H <- c(0, 0.2, 0.2, 0.7)
  lh <- c(0, log(0.2), 0, log(0.5))
  lf <- status * (lh + eta) - H * exp(eta)
  e <- exp(lp + lf)
  hand <- e / rowSums(e)
  engine <- tau_from_evidence(lcc_log_evidence(x, time, status, eta, list(b, b), gate), 1)
  expect_equal(unname(engine), unname(hand), tolerance = 1e-10)
})

test_that("the gating M-step is the weighted multinomial MLE fitted by glmnet", {
  d <- gemcox_simulate(n = 300, p = 4, mu_sep = 1.5, beta_sep = 1, seed = 41)
  Xg <- apply_x_scaler(d$X, make_x_scaler(d$X))
  tau <- with_local_seed(3, { u <- stats::runif(300); cbind(C1 = u, C2 = 1 - u) })
  g <- gate_fit(Xg, tau, lambda_gate = 0.05)
  ## (a) the same glmnet call made directly
  direct <- glmnet::glmnet(Xg, tau, family = "multinomial", alpha = 0, lambda = 0.05,
                           standardize = FALSE)
  pd <- predict(direct, Xg, s = 0.05, type = "response")[, , 1]
  expect_lt(max(abs(exp(gate_logprob(g, Xg)) - pd)), 1e-12)
  ## (b) an independent formulation: each subject once per class, weighted by tau
  big <- glmnet::glmnet(rbind(Xg, Xg), factor(rep(1:2, each = 300)), family = "multinomial",
                        weights = c(tau[, 1], tau[, 2]), alpha = 0, lambda = 0.05,
                        standardize = FALSE, thresh = 1e-12, maxit = 1e6)
  pb <- predict(big, Xg, s = 0.05, type = "response")[, , 1]
  expect_lt(max(abs(exp(gate_logprob(g, Xg)) - pb)), 1e-4)
})

## Monotonicity. glmnet normalises the tau weights within each Cox update,
## so cluster k's ridge term is n_k(tau) * lambda / 2 * |beta_k|^2 and its
## weight changes with tau at every iteration. The EM inequality holds with
## the same weights on both sides, but the tracked trace changes weights
## between iterations, so it can fall slightly. Seeds 701-706 are the six
## datasets examined in development, including the two with small decreases
## (705: -7.5e-9 relative; 706: -3.0e-6 relative). A diagnostic before
## Phase C found the same decreases with glmnet's threshold at 1e-14, and
## none when both sides use the same weights. Tolerance: 1e-5 relative,
## set after the decreases were observed.
test_that("the penalised log-likelihood is non-decreasing across EM iterations", {
  for (s in 1:6) {
    d <- gemcox_simulate(n = 400, p = 6, mu_sep = 1, beta_sep = 2, seed = 700 + s)
    f <- suppressWarnings(lcc_fit(d$X, time = d$time, status = d$status, tol = 1e-8,
                                  max_iter = 200))
    expect_gte(min(diff(f$pen_loglik)) / abs(utils::tail(f$pen_loglik, 1)), -1e-5)
  }
})

test_that("recovers the coefficient contrast and gating direction on a large easy dataset", {
  d <- gemcox_simulate(n = 5000, p = 10, mu_sep = 3, beta_sep = 3, seed = 42)
  f <- suppressWarnings(lcc_fit(d$X, time = d$time, status = d$status))
  ctr <- unit(f$beta[, 1] - f$beta[, 2])
  tru <- unit(d$beta_list[[1]] - d$beta_list[[2]])
  expect_gt(abs(sum(ctr * tru)), 0.95)
  gdir <- unit(f$gate_slopes_orig[, 2] - f$gate_slopes_orig[, 1])
  prof <- unit(colMeans(d$X[d$Z == 2, ]) - colMeans(d$X[d$Z == 1, ]))
  expect_gt(abs(sum(gdir * prof)), 0.9)
})

test_that("intercept-only gating reduces the E-step to constant mixing proportions", {
  d <- gemcox_simulate(n = 300, p = 4, mu_sep = 1, beta_sep = 2, seed = 43)
  Xg <- apply_x_scaler(d$X, make_x_scaler(d$X))
  tau <- with_local_seed(4, { u <- stats::runif(300); cbind(C1 = u, C2 = 1 - u) })
  g <- gate_fit(Xg, tau, lambda_gate = 0.05, slopes = FALSE)
  P <- exp(gate_logprob(g, Xg))
  expect_equal(max(abs(sweep(P, 2, P[1, ]))), 0, tolerance = 1e-14)
  expect_equal(unname(P[1, ]), unname(colMeans(tau)), tolerance = 1e-12)
  f <- suppressWarnings(lcc_fit(d$X, time = d$time, status = d$status, gate_slopes = FALSE))
  Pf <- lcc_predict_tau(f, d$X)
  expect_equal(max(abs(sweep(Pf, 2, Pf[1, ]))), 0, tolerance = 1e-14)
})

test_that("settings record the gating penalty", {
  d <- gemcox_simulate(n = 200, p = 3, mu_sep = 1, beta_sep = 1, seed = 44)
  f <- suppressWarnings(lcc_fit(d$X, time = d$time, status = d$status))
  expect_equal(f$settings$lambda_gate, 0.05)
  expect_equal(f$settings$gate_alpha, 0)
  expect_equal(f$settings$lambda, 0.05)
})

## Harness identity: adding the competitor must not change any stored
## dataset. The hashes below were computed with the package at the
## sim-freeze-v2 tag, before lcc was added, for replicates whose stored
## event counts were checked against the frozen results files.
test_that("adding the competitor does not change any stored simulated dataset", {
  old_kind <- RNGkind()
  on.exit(RNGkind(old_kind[1], old_kind[2], old_kind[3]))
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  mk_seed <- function(e, c, r) 20260101L + e * 1e6 + c * 1e3 + r
  h <- function(d) {
    f <- tempfile(); on.exit(unlink(f))
    saveRDS(list(X = d$X, time = d$time, status = d$status, Z = d$Z), f, compress = FALSE)
    unname(tools::md5sum(f))
  }
  frozen <- list(
    list(1, 3, 17,  list(n = 400,  p = 10, mu_sep = 0,   beta_sep = 2), "3a1ee11f90c7da11b590a0f4f8d4ccdd"),
    list(2, 6, 101, list(n = 800,  p = 10, mu_sep = 0.5, beta_sep = 2), "8c6e580198741adfeddf5b29d468b1b8"),
    list(3, 1, 5,   list(n = 200,  p = 10, mu_sep = 0.5, beta_sep = 2), "2ec8dcf0e5fc30f85511920fc78568df"),
    list(31, 4, 60, list(n = 800,  p = 80, mu_sep = 0.5, beta_sep = 2), "f4bb1bb015ca5ccaa500ab533733d0ce"),
    list(4, 4, 3,   list(n = 6400, p = 10, mu_sep = 0,   beta_sep = 2), "331a8d76e3dc1172437d4a13758dd9fc"),
    list(5, 1, 400, list(n = 400,  p = 10, mu_sep = 0,   beta_sep = 0, K_true = 1), "1fbe2179c00fdc577b31eb7e050db9b5"),
    list(5, 2, 9,   list(n = 400,  p = 10, mu_sep = 0,   beta_sep = 0, K_true = 1, rho = 0.5), "a116f7bd292f5353a290aecd4c4dde69"),
    list(24, 4, 50, list(n = 400,  p = 10, mu_sep = 1,   beta_sep = 0, rho = 0.5), "86b01cb66785424df0e74fc3cf923bf1"),
    list(25, 12, 7, list(n = 400,  p = 10, mu_sep = 0,   beta_sep = 3, K_true = 3), "4f31051ce32306a72304ecb757947961"),
    list(23, 1, 2,  list(n = 117,  p = 9,  mu_sep = 0,   beta_sep = 0, K_true = 1, event_rate = 0.44),
         "507ffb06dfa775f5170bf688a0bff975"))
  for (x in frozen) {
    d <- do.call(gemcox_simulate, c(x[[4]], list(seed = mk_seed(x[[1]], x[[2]], x[[3]]))))
    expect_identical(h(d), x[[5]], label = sprintf("dataset exp %d cell %d rep %d", x[[1]], x[[2]], x[[3]]))
  }
})
