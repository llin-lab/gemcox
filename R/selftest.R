## selftest.R -- gemcox_selftest(): quick engine checks for the simulation
## pipeline. The same assertions are in tests/testthat with more cases.

# Independent per-subject survival log-density, written without
# findInterval() or get_cumhaz_at(), for cross-checking log_surv_density().
direct_log_surv_density <- function(time, status, eta, baseline) {
  vapply(seq_along(time), function(i) {
    on_grid <- baseline$time <= time[i]
    H <- sum(baseline$jump[on_grid])
    h <- if (status[i] == 1) baseline$jump[max(which(on_grid))] else 1
    status[i] * (log(h) + eta[i]) - H * exp(eta[i])
  }, numeric(1))
}

#' Engine self-test
#'
#' Runs fast checks of the core engine on simulated data:
#' 1. `get_cumhaz_at()` on a grid with query times before the first event
#'    (the indexing defect that attenuated earlier fits);
#' 2. the E-step survival term at K = 2 equals an independent per-subject
#'    computation, including subjects observed before the first event;
#' 3. K = 1 reproduces `survival::coxph` (Breslow ties, `lambda = 0`) for the
#'    shared and the cluster-specific baseline;
#' 4. K = 1 does not depend on `gamma`;
#' 5. with `gamma = 0`, membership weights recomputed by [predict.gemcox()]
#'    on the training features equal the fitted weights.
#'
#' @param verbose Print each check.
#' @return `TRUE` if every check passes, otherwise `FALSE`; the per-check
#'   results are in `attr(, "checks")`.
#' @examples
#' \donttest{
#' gemcox_selftest(verbose = TRUE)
#' }
#' @export
gemcox_selftest <- function(verbose = FALSE) {
  res <- list()
  add <- function(name, ok, detail) {
    res[[length(res) + 1]] <<- data.frame(check = name, passed = isTRUE(ok),
                                          detail = detail, stringsAsFactors = FALSE)
    if (verbose) cat(sprintf("  [%s] %s: %s\n", if (isTRUE(ok)) "PASS" else "FAIL",
                             name, detail))
  }
  safely <- function(name, expr) {
    tryCatch(expr, error = function(e) add(name, FALSE, paste("error:", conditionMessage(e))))
  }

  ## 1. get_cumhaz_at indexing
  safely("get_cumhaz_at", {
    got <- get_cumhaz_at(list(time = c(2, 4, 6), cumhaz = c(.1, .5, .9)), c(1, 3, 5, 7))
    add("get_cumhaz_at", identical(got, c(0, .1, .5, .9)),
        paste(format(got), collapse = " "))
  })

  ## 2. E-step survival term at K = 2, with subjects before the first event
  safely("log_surv_density (K = 2)", {
    d <- gemcox_simulate(n = 300, p = 4, mu_sep = 1, beta_sep = 2, seed = 11)
    first_ev <- min(d$time[d$status == 1])
    d$X <- rbind(d$X, d$X[1:3, ])
    d$time <- c(d$time, first_ev * c(0.2, 0.5, 0.8))
    d$status <- c(d$status, 0L, 0L, 0L)
    fit <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, K = 2))
    E <- predict_eta(fit, d$X)
    dev <- max(vapply(1:2, function(k) {
      max(abs(log_surv_density(d$time, d$status, E[, k], fit$baseline[[k]]) -
                direct_log_surv_density(d$time, d$status, E[, k], fit$baseline[[k]])))
    }, 0))
    add("log_surv_density (K = 2)", dev < 1e-10,
        sprintf("max |difference| %.2e over %d subjects (3 before first event)",
                dev, length(d$time)))
  })

  ## 3-4. K = 1 reproduces coxph, for both baselines, independent of gamma
  d1 <- gemcox_simulate(n = 400, p = 5, mu_sep = 0, beta_sep = 0, K_true = 1, seed = 12)
  ref <- as.numeric(coef(coxph(Surv(d1$time, d1$status) ~ d1$X, ties = "breslow")))
  for (bl in c("shared", "cluster")) {
    nm <- sprintf("K = 1 matches coxph (%s)", bl)
    safely(nm, {
      b <- as.numeric(coef(gemcox(d1$X, time = d1$time, status = d1$status, K = 1,
                                  lambda = 0, baseline = bl)))
      ratio <- sqrt(sum(b^2)) / sqrt(sum(ref^2))
      mad <- max(abs(b - ref))
      add(nm, abs(ratio - 1) < 0.02 && mad < 0.02,
          sprintf("norm ratio %.4f, max |diff| %.4f", ratio, mad))
    })
  }
  safely("K = 1 invariant to gamma", {
    bs <- sapply(c(0, 0.5, 1, 2), function(g) {
      as.numeric(coef(gemcox(d1$X, time = d1$time, status = d1$status, K = 1,
                             lambda = 0, gamma = g)))
    })
    dev <- max(abs(bs - bs[, 1]))
    add("K = 1 invariant to gamma", dev < 1e-8, sprintf("max |diff| %.2e", dev))
  })

  ## 5. predict(type = "tau") reproduces the fitted weights at gamma = 0
  safely("predict tau = fit tau (gamma = 0)", {
    d <- gemcox_simulate(n = 300, p = 4, mu_sep = 2, beta_sep = 1, seed = 13)
    fit <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, gamma = 0))
    dev <- max(abs(predict(fit, d$X, type = "tau") - fit$tau))
    add("predict tau = fit tau (gamma = 0)", dev < 1e-6, sprintf("max |diff| %.2e", dev))
  })

  checks <- do.call(rbind, res)
  ok <- all(checks$passed)
  if (verbose) cat(if (ok) "gemcox_selftest: all checks passed.\n" else
    "gemcox_selftest: FAILED.\n")
  structure(ok, checks = checks)
}
