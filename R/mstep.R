## mstep.R -- M-step: closed-form GMM update (pi, mu, Sigma), tau-weighted
## ridge Cox updates, and the per-cluster weighted Breslow baseline.
##
## A "coxfit" is list(beta_s, baseline, lambda, ok, reason, eff_events):
## beta_s are coefficients on the standardized Cox features (the scaler is
## held once at fit level), baseline is the step function used for that
## cluster's survival density.

update_pi <- function(tau, pi_floor) {
  pi_k <- colMeans(tau)
  pi_k <- pmax(pi_k, pi_floor)
  pi_k / sum(pi_k)
}

# Weighted means and covariances of the scaled GMM features.
# initial = TRUE sets an empty cluster to N(0, I); otherwise it keeps its
# previous parameters.
update_gmm <- function(X_gmm_s, tau, mu, Sigma, diagonal, ridge, initial = FALSE) {
  q <- ncol(X_gmm_s)
  for (k in seq_len(ncol(tau))) {
    w <- tau[, k]
    wsum <- sum(w)
    if (!is.finite(wsum) || wsum < 1e-8) {
      if (initial) {
        mu[, k] <- rep(0, q)
        Sigma[[k]] <- diag(1, q)
      }
      next
    }
    mu[, k] <- colSums(X_gmm_s * w) / wsum
    xc <- sweep(X_gmm_s, 2, mu[, k], "-")
    S <- (t(xc * sqrt(w)) %*% (xc * sqrt(w))) / wsum
    S <- (S + t(S)) / 2
    if (diagonal) S <- diag(diag(S), q)
    Sigma[[k]] <- as.matrix(S) + diag(ridge, q)
  }
  list(mu = mu, Sigma = Sigma)
}

# Weighted penalised Cox fit (glmnet, family = "cox") with its own weighted
# Breslow baseline. Used for the global fit, for initialisation, and as the
# M-step when baseline = "cluster".
cox_glmnet_weighted <- function(time, status, Xs, w, lambda, alpha,
                                min_eff_events, control, fallback_fit = NULL) {
  status <- as.integer(status)
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  n_eff_event <- sum(w * status)

  if (!is.finite(n_eff_event) || n_eff_event < min_eff_events) {
    if (!is.null(fallback_fit)) {
      out <- fallback_fit
      out$ok <- FALSE
      out$reason <- "fallback_global_low_events"
      out$eff_events <- n_eff_event
      return(out)
    }
    return(list(
      beta_s = rep(0, ncol(Xs)),
      baseline = breslow_weighted(time, status, rep(0, length(time)), w, control),
      lambda = NA_real_, ok = FALSE, reason = "zero_effective_events",
      eff_events = n_eff_event))
  }

  fit <- tryCatch({
    gfit <- glmnet::glmnet(x = Xs, y = Surv(time, status), family = "cox",
                           weights = w, alpha = alpha, lambda = lambda,
                           standardize = FALSE, maxit = 1e5)
    list(beta_s = as.numeric(coef(gfit, s = lambda)), lambda = lambda,
         ok = TRUE, reason = "ok")
  }, error = function(e) NULL)

  if (is.null(fit)) {
    gfit <- tryCatch(glmnet::glmnet(x = Xs, y = Surv(time, status), family = "cox",
                                    weights = w, alpha = 0, lambda = lambda * 2,
                                    standardize = FALSE, maxit = 1e5),
                     error = function(e) NULL)
    if (!is.null(gfit)) {
      fit <- list(beta_s = as.numeric(coef(gfit, s = lambda * 2)),
                  lambda = lambda * 2, ok = TRUE, reason = "ridge_fallback")
    } else if (!is.null(fallback_fit)) {
      out <- fallback_fit
      out$ok <- FALSE
      out$reason <- "fallback_global_glmnet_failed"
      out$eff_events <- n_eff_event
      return(out)
    } else {
      fit <- list(beta_s = rep(0, ncol(Xs)), lambda = NA_real_, ok = FALSE,
                  reason = "glmnet_failed")
    }
  }

  beta_s <- fit$beta_s
  beta_s[!is.finite(beta_s)] <- 0
  eta <- clamp_eta(as.numeric(Xs %*% beta_s), control$eta_cap, "eta_cap")
  list(beta_s = beta_s,
       baseline = breslow_weighted(time, status, eta, w, control),
       lambda = fit$lambda, ok = isTRUE(fit$ok), reason = fit$reason,
       eff_events = n_eff_event)
}

# One cluster's coefficient update against a fixed shared baseline, via the
# Poisson working model: status_i ~ Poisson(Lambda0(t_i) exp(x_i' beta_k)),
# weights tau_ik, offset log Lambda0(t_i).
#
# intercept = FALSE: the shared baseline supplies the only intercept. A
# per-cluster intercept would let clusters differ in baseline risk, which
# contradicts the shared-baseline model (clusters differ only through beta).
cox_poisson_weighted <- function(time, status, Xs, w, baseline, lambda, alpha,
                                 min_eff_events, fallback_fit = NULL) {
  status <- as.integer(status)
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  eff_events <- sum(w * status)

  if (!is.finite(eff_events) || eff_events < min_eff_events) {
    if (!is.null(fallback_fit)) {
      out <- fallback_fit
      out$baseline <- baseline
      out$ok <- FALSE
      out$reason <- "fallback_low_events"
      out$eff_events <- eff_events
      return(out)
    }
    return(list(beta_s = rep(0, ncol(Xs)), baseline = baseline,
                lambda = NA_real_, ok = FALSE, reason = "low_events",
                eff_events = eff_events))
  }

  offset <- log(pmax(get_cumhaz_at(baseline, time), 1e-10))

  fit <- tryCatch({
    gfit <- glmnet::glmnet(x = Xs, y = status, family = "poisson",
                           intercept = FALSE, weights = w, offset = offset,
                           alpha = alpha, lambda = lambda,
                           standardize = FALSE, maxit = 1e5)
    list(beta_s = as.numeric(coef(gfit, s = lambda))[-1], ok = TRUE,
         reason = "ok", lambda = lambda)
  }, error = function(e) NULL)

  if (is.null(fit)) {
    if (!is.null(fallback_fit)) {
      out <- fallback_fit
      out$baseline <- baseline
      out$ok <- FALSE
      out$reason <- "fallback_poisson_failed"
      out$eff_events <- eff_events
      return(out)
    }
    fit <- list(beta_s = rep(0, ncol(Xs)), ok = FALSE, reason = "poisson_failed",
                lambda = lambda)
  }

  beta_s <- fit$beta_s
  beta_s[!is.finite(beta_s)] <- 0
  list(beta_s = beta_s, baseline = baseline, lambda = fit$lambda,
       ok = isTRUE(fit$ok), reason = fit$reason, eff_events = eff_events)
}
