## E-step correctness. The K = 1 test cannot see E-step errors (tau is
## identically 1), and the log_surv_density() indexing defect survived in
## GeMCox_combined.R after get_cumhaz_at() was fixed. These tests would have
## caught it.

# Independent per-subject computation, written in the test on purpose.
direct_lsd <- function(time, status, eta, b) {
  vapply(seq_along(time), function(i) {
    past <- b$time <= time[i]
    H <- sum(b$jump[past])
    h <- if (status[i] == 1) b$jump[b$time == time[i]] else 1
    status[i] * (log(h) + eta[i]) - H * exp(eta[i])
  }, 0)
}

# The legacy (buggy) evaluation, kept to show the test detects it.
legacy_lsd <- function(time, status, eta, b) {
  idx <- findInterval(time, b$time)
  Lambda <- ifelse(idx > 0, b$cumhaz[idx], 0)
  haz0 <- rep(1, length(time))
  ev <- which(status == 1)
  haz0[ev] <- pmax(b$jump[pmax(idx[ev], 1L)], 1e-12)
  status * (log(haz0) + eta) - Lambda * exp(eta)
}

test_that("K = 2 survival log-density equals a direct computation for every subject", {
  d <- gemcox_simulate(n = 300, p = 4, mu_sep = 1, beta_sep = 2, seed = 61)
  first_ev <- min(d$time[d$status == 1])
  ## five censored subjects before the first event, at the start and middle
  pos <- c(1, 2, 100, 101, 250)
  X <- d$X
  time <- d$time
  status <- d$status
  X[pos, ] <- X[pos + 1, ]
  time[pos] <- first_ev * c(0.1, 0.3, 0.5, 0.7, 0.9)
  status[pos] <- 0L
  expect_gt(sum(time < min(time[status == 1])), 0)

  fit <- suppressWarnings(gemcox(X, time = time, status = status, K = 2))
  E <- predict_eta(fit, X)
  for (k in 1:2) {
    b <- fit$baseline[[k]]
    got <- log_surv_density(time, status, E[, k], b)
    expect_lt(max(abs(got - direct_lsd(time, status, E[, k], b))), 1e-10)

    ## row-order invariance: the recycling defect breaks exactly this
    o <- with_local_seed(k, sample(length(time)))
    expect_equal(log_surv_density(time[o], status[o], E[o, k], b), got[o],
                 tolerance = 1e-12)

    ## and the legacy indexing fails the same comparison
    expect_gt(max(abs(legacy_lsd(time, status, E[, k], b) - got)), 1e-3)
  }
})

test_that("E-step tau matches a hand computation on a tiny dataset", {
  x <- matrix(c(-1, 0.5, 2), ncol = 1)
  time <- c(1, 2, 3)
  status <- c(1, 0, 1)
  pi_k <- c(0.4, 0.6)
  mu <- matrix(c(-0.5, 1), nrow = 1)
  Sigma <- list(matrix(1), matrix(2))
  b <- list(time = c(1, 3), jump = c(0.2, 0.5), cumhaz = c(0.2, 0.7))
  eta <- cbind(c(0.1, -0.2, 0.4), c(-0.3, 0.5, 0))

  ## by hand: log pi_k + log N(x; mu_k, Sigma_k) + gamma * log f_k, where
  ## Lambda0 = (0.2, 0.2, 0.7) and dLambda0 at the events = (0.2, -, 0.5)
  lg <- cbind(log(0.4) - 0.5 * log(2 * pi * 1) - (x - (-0.5))^2 / 2,
              log(0.6) - 0.5 * log(2 * pi * 2) - (x - 1)^2 / 4)
  ls <- status * (c(log(0.2), 0, log(0.5)) + eta) - c(0.2, 0.2, 0.7) * exp(eta)
  hand <- function(gamma, temp) {
    e <- exp((lg + gamma * ls) / temp)
    e / rowSums(e)
  }
  engine <- function(gamma, temp) {
    L <- estep_log_evidence(x, time, status, eta, list(b, b), pi_k, mu, Sigma,
                            ridge = 0, gmm_scale = 1, gamma = gamma)
    tau_from_evidence(L, temp)
  }
  golden <- list(
    g1 = rbind(c(0.758273633060432, 0.241726366939568),
               c(0.418144902477907, 0.581855097522093),
               c(0.053243368673719, 0.946756631326281)),
    g0 = rbind(c(0.693409655058754, 0.306590344941246),
               c(0.378389049476261, 0.621610950523740),
               c(0.050503389011044, 0.949496610988956)),
    g2 = rbind(c(0.813113223087430, 0.186886776912570),
               c(0.458993483621453, 0.541006516378547),
               c(0.056123214851056, 0.943876785148944)),
    t2 = rbind(c(0.639136660351624, 0.360863339648376),
               c(0.458794488653990, 0.541205511346010),
               c(0.191687174032211, 0.808312825967789)))

  expect_equal(engine(1, 1), hand(1, 1), tolerance = 1e-12)
  expect_equal(engine(1, 1), golden$g1, tolerance = 1e-12)
  expect_equal(engine(0, 1), golden$g0, tolerance = 1e-12)
  expect_equal(engine(2, 1), golden$g2, tolerance = 1e-12)
  expect_equal(engine(1, 2), golden$t2, tolerance = 1e-12)
})
