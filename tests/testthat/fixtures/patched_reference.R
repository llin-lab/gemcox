## Reference implementation for the equivalence tests (test-equivalence.R).
##
## Verbatim copies of the functions gemcox() was ported from, extracted by
## source reference (not deparsed) from GeMCox_combined_PATCHED.R:
##   sha256 5e42cb0fac6ddc7aa706ba9e86b64453de6db5884d303447809bfc5373e9b5cb
##   git blob f6f767efc6605b494ab2f14396a90a6bb334b2d6
## Do not edit. If PATCHED changes, re-extract rather than hand-editing.
## Sourced into an environment whose parent is the gemcox namespace, so
## glmnet(), Surv() and coef() resolve through the package imports.

clamp <- function(x, lo = -30, hi = 30) pmin(pmax(x, lo), hi)

make_x_scaler <- function(X) {
  X <- as.matrix(X)
  center <- colMeans(X, na.rm = TRUE)
  scalev <- apply(X, 2, stats::sd, na.rm = TRUE)
  scalev[!is.finite(scalev) | scalev <= 0] <- 1
  list(center = center, scale = scalev)
}

apply_x_scaler <- function(X, scaler) {
  X <- as.matrix(X)
  sweep(sweep(X, 2, scaler$center, "-"), 2, scaler$scale, "/")
}

eta_from_scaled <- function(X, beta_s, scaler, eta_cap = 12) {
  Xs <- apply_x_scaler(X, scaler)
  eta <- as.numeric(Xs %*% beta_s)
  clamp(eta, -eta_cap, eta_cap)
}

rowLogSumExp <- function(A) {
  m <- apply(A, 1, max)
  m[!is.finite(m)] <- 0
  s <- rowSums(exp(A - m))
  s <- pmax(s, 1e-300)
  m + log(s)
}

dmvnorm_log <- function(X, mu, Sigma, ridge = 1e-6) {
  X <- as.matrix(X)
  mu <- as.numeric(mu)
  d <- length(mu)
  xc <- sweep(X, 2, mu, "-")

  S2 <- as.matrix(Sigma) + diag(ridge, d)
  cholS <- tryCatch(chol(S2), error = function(e) NULL)
  if (is.null(cholS)) {
    S2 <- as.matrix(Sigma) + diag(1e-3, d)
    cholS <- chol(S2)
  }

  Sinv <- chol2inv(cholS)
  logdet <- 2 * sum(log(diag(cholS)))
  quad <- rowSums((xc %*% Sinv) * xc)
  -0.5 * (d * log(2 * pi) + logdet + quad)
}

breslow_weighted <- function(time, status, eta, w,
                             denom_floor = 1e-8,
                             max_jump = 10,
                             max_cumhaz = 500) {
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  status <- as.integer(status)
  eta <- clamp(as.numeric(eta), -30, 30)
  exp_eta <- exp(eta)

  if (sum(w * status) < 1e-10) {
    return(list(time = numeric(0), jump = numeric(0), cumhaz = numeric(0)))
  }

  uniqT <- sort(unique(time[status == 1]))
  jump <- numeric(length(uniqT))

  for (m in seq_along(uniqT)) {
    t_m <- uniqT[m]
    num <- sum(w[time == t_m & status == 1])
    at_risk <- (time >= t_m)
    denom <- sum(w[at_risk] * exp_eta[at_risk])
    denom <- max(denom, denom_floor)
    jump[m] <- min(num / denom, max_jump)
  }

  cumhaz <- pmin(cumsum(jump), max_cumhaz)
  list(time = uniqT, jump = jump, cumhaz = cumhaz)
}

get_cumhaz_at <- function(baseline, times) {
  ## FIXED: x[0] silently drops elements, which made ifelse() recycle a
  ## short vector and assign wrong cumulative hazards to every subject
  ## after the first pre-grid time. Index only the valid positions.
  tt <- as.numeric(times)
  if (length(baseline$time) == 0) return(rep(0, length(tt)))
  idx <- findInterval(tt, baseline$time)
  out <- numeric(length(tt))
  ok  <- idx > 0
  out[ok] <- baseline$cumhaz[idx[ok]]
  out
}

cox_glmnet_weighted <- function(time, status, X, w,
                                lambda = 0.1,
                                alpha = 0,
                                min_eff_events = 0.5,
                                denom_floor = 1e-8,
                                max_jump = 10,
                                max_cumhaz = 500,
                                fallback_fit = NULL) {
  X <- as.matrix(X)
  status <- as.integer(status)
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0

  scaler <- make_x_scaler(X)
  Xs <- apply_x_scaler(X, scaler)
  n_eff_event <- sum(w * status)

  if (!is.finite(n_eff_event) || n_eff_event < min_eff_events) {
    if (!is.null(fallback_fit)) {
      out <- fallback_fit
      out$ok <- FALSE
      out$reason <- "fallback_global_low_events"
      out$eff_events <- n_eff_event
      return(out)
    }
    beta_s <- rep(0, ncol(X))
    baseline <- breslow_weighted(
      time, status, eta = rep(0, length(time)), w = w,
      denom_floor = denom_floor,
      max_jump = max_jump,
      max_cumhaz = max_cumhaz
    )
    return(list(
      beta_s = beta_s,
      beta_orig = beta_s / scaler$scale,
      scaler = scaler,
      baseline = baseline,
      lambda = NA_real_,
      ok = FALSE,
      reason = "zero_effective_events",
      eff_events = n_eff_event
    ))
  }

  if (is.character(lambda)) {
    warning("cox_glmnet_weighted: lambda should be a numeric scalar. Falling back to 0.1.")
    lambda <- 0.1
  }

  fit <- tryCatch({
    gfit <- glmnet(
      x = Xs,
      y = Surv(time, status),
      family = "cox",
      weights = w,
      alpha = alpha,
      lambda = lambda,
      standardize = FALSE,
      maxit = 1e5
    )
    beta_s <- as.numeric(coef(gfit, s = lambda))
    list(beta_s = beta_s, lambda = lambda, ok = TRUE, reason = "ok")
  }, error = function(e) NULL)

  if (is.null(fit)) {
    gfit <- tryCatch(glmnet(
      x = Xs,
      y = Surv(time, status),
      family = "cox",
      weights = w,
      alpha = 0,
      lambda = lambda * 2,
      standardize = FALSE,
      maxit = 1e5
    ), error = function(e) NULL)

    if (!is.null(gfit)) {
      fit <- list(
        beta_s = as.numeric(coef(gfit, s = lambda * 2)),
        lambda = lambda * 2,
        ok = TRUE,
        reason = "ridge_fallback"
      )
    } else if (!is.null(fallback_fit)) {
      out <- fallback_fit
      out$ok <- FALSE
      out$reason <- "fallback_global_glmnet_failed"
      out$eff_events <- n_eff_event
      return(out)
    } else {
      fit <- list(
        beta_s = rep(0, ncol(X)),
        lambda = NA_real_,
        ok = FALSE,
        reason = "glmnet_failed"
      )
    }
  }

  beta_s <- fit$beta_s
  beta_s[!is.finite(beta_s)] <- 0
  beta_orig <- beta_s / scaler$scale

  eta <- eta_from_scaled(X, beta_s, scaler, eta_cap = 12)
  baseline <- breslow_weighted(
    time, status, eta = eta, w = w,
    denom_floor = denom_floor,
    max_jump = max_jump,
    max_cumhaz = max_cumhaz
  )

  list(
    beta_s = beta_s,
    beta_orig = beta_orig,
    scaler = scaler,
    baseline = baseline,
    lambda = fit$lambda,
    ok = isTRUE(fit$ok),
    reason = fit$reason,
    eff_events = n_eff_event
  )
}

log_surv_density <- function(time, status, X, coxfit) {
  X <- as.matrix(X)
  status <- as.integer(status)
  n <- nrow(X)

  eta <- eta_from_scaled(X, coxfit$beta_s, coxfit$scaler, eta_cap = 12)
  exp_eta <- exp(eta)

  ev_t <- coxfit$baseline$time
  jump <- coxfit$baseline$jump
  cumhaz <- coxfit$baseline$cumhaz

  if (length(ev_t) == 0) {
    return(rep(0, n))
  }

  idx <- findInterval(time, ev_t)
  Lambda <- numeric(n)
  positive_idx <- idx > 0L
  Lambda[positive_idx] <- cumhaz[idx[positive_idx]]

  haz0 <- rep(1, n)
  ev_idx <- which(status == 1)
  if (length(ev_idx) > 0) {
    safe_idx <- pmax(idx[ev_idx], 1L)    # guard: idx=0 if test time < first training event
    haz0[ev_idx] <- pmax(jump[safe_idx], 1e-12)
  }

  status * (log(haz0) + eta) - Lambda * exp_eta
}

initialize_tau <- function(X_gmm_s, X_cox, global_coxfit,
                           K = 2,
                           init_method = c("kmeans", "supervised"),
                           init_eta_weight = 0.75,
                           init_seed = 1) {
  init_method <- match.arg(init_method)
  set.seed(init_seed)

  if (K == 1) {
    tau <- matrix(1, nrow(X_gmm_s), 1)
    colnames(tau) <- "C1"
    return(tau)
  }

  if (init_method == "supervised") {
    eta0 <- eta_from_scaled(X_cox, global_coxfit$beta_s, global_coxfit$scaler, eta_cap = 12)
    eta0 <- as.numeric(scale(eta0))
    eta0[!is.finite(eta0)] <- 0
    X_init <- cbind(X_gmm_s, init_eta_weight * eta0)
  } else {
    X_init <- X_gmm_s
  }

  km <- kmeans(X_init, centers = K, nstart = 50)
  Z0 <- km$cluster

  tau <- matrix(0, nrow(X_gmm_s), K)
  tau[cbind(seq_len(nrow(X_gmm_s)), Z0)] <- 1
  colnames(tau) <- paste0("C", seq_len(K))
  tau
}

gemcox_full <- function(
    X_gmm, X_cox, time, status,
    K = 2,
    lambda = NULL,
    alpha = 0,
    max_iter = 50,
    tol = 1e-4,
    verbose = TRUE,
    gmm_ridge = 1e-6,
    gmm_diag = TRUE,
    normalize_gmm_by_dim = TRUE,
    surv_weight = 1.0,
    temp = 1.0,
    pi_floor = 0.02,
    min_eff_events = 0.5,
    init_method = c("supervised", "kmeans"),
    init_eta_weight = 0.75,
    init_seed = 1) {

  init_method <- match.arg(init_method)

  X_gmm <- as.matrix(X_gmm)
  X_cox <- as.matrix(X_cox)
  time <- as.numeric(time)
  status <- as.integer(status)

  n <- nrow(X_gmm)
  stopifnot(nrow(X_cox) == n, length(time) == n, length(status) == n)
  stopifnot(all(status %in% c(0, 1)))

  q <- ncol(X_gmm)
  p <- ncol(X_cox)
  gmm_log_density_scale <- if (isTRUE(normalize_gmm_by_dim)) 1 / max(q, 1) else 1

  colnames(X_gmm) <- make.names(colnames(X_gmm), unique = TRUE)
  colnames(X_cox) <- make.names(colnames(X_cox), unique = TRUE)

  gmm_scaler <- make_x_scaler(X_gmm)
  X_gmm_s <- apply_x_scaler(X_gmm, gmm_scaler)

  if (is.null(lambda) || is.character(lambda)) {
    if (verbose) cat("Pre-selecting lambda via cv.glmnet (K=1, unweighted)...\n")
    lambda_use <- select_lambda(time, status, X_cox, alpha = alpha)
    if (verbose) cat(sprintf("  Selected lambda = %.5f\n", lambda_use))
  } else {
    lambda_use <- as.numeric(lambda)
  }

  global_coxfit <- cox_glmnet_weighted(
    time = time,
    status = status,
    X = X_cox,
    w = rep(1, n),
    lambda = lambda_use,
    alpha = alpha,
    min_eff_events = 0
  )

  tau <- initialize_tau(
    X_gmm_s = X_gmm_s,
    X_cox = X_cox,
    global_coxfit = global_coxfit,
    K = K,
    init_method = init_method,
    init_eta_weight = init_eta_weight,
    init_seed = init_seed
  )

  pi_c <- pmax(colMeans(tau), 1e-8)
  pi_c <- pi_c / sum(pi_c)

  mu_gmm <- matrix(0, q, K)
  Sigma_list <- vector("list", K)
  coxfit <- vector("list", K)
  beta <- matrix(0, p, K, dimnames = list(colnames(X_cox), paste0("C", seq_len(K))))

  ## Initial M-step
  for (k in seq_len(K)) {
    w <- tau[, k]
    wsum <- sum(w)
    if (!is.finite(wsum) || wsum < 1e-8) {
      mu_gmm[, k] <- rep(0, q)
      Sigma_list[[k]] <- diag(1, q)
    } else {
      mu_gmm[, k] <- colSums(X_gmm_s * w) / wsum
      xc <- sweep(X_gmm_s, 2, mu_gmm[, k], "-")
      S <- (t(xc * sqrt(w)) %*% (xc * sqrt(w))) / wsum
      S <- (S + t(S)) / 2
      if (gmm_diag) S <- diag(diag(S), q)
      Sigma_list[[k]] <- as.matrix(S) + diag(gmm_ridge, q)
    }

    coxfit[[k]] <- cox_glmnet_weighted(
      time = time,
      status = status,
      X = X_cox,
      w = w,
      lambda = lambda_use,
      alpha = alpha,
      min_eff_events = min_eff_events,
      fallback_fit = global_coxfit
    )
    beta[, k] <- coxfit[[k]]$beta_orig
  }

  loglik_trace <- numeric(max_iter)
  temp_use <- if (!is.finite(temp) || temp <= 0) 1.0 else temp

  for (iter in seq_len(max_iter)) {
    pi_c <- pmax(pi_c, 1e-8)
    pi_c <- pi_c / sum(pi_c)

    ## E-step
    log_resp <- matrix(NA_real_, n, K)
    for (k in seq_len(K)) {
      lgmm <- dmvnorm_log(X_gmm_s, mu_gmm[, k], Sigma_list[[k]], ridge = gmm_ridge)
      lsurv <- log_surv_density(time, status, X_cox, coxfit[[k]])
      log_resp[, k] <- log(pi_c[k]) + gmm_log_density_scale * lgmm + surv_weight * lsurv
    }
    log_resp[!is.finite(log_resp)] <- -700

    if (abs(temp_use - 1) < 1e-9) {
      loglik_i <- rowLogSumExp(log_resp)
    } else {
      loglik_i <- rowLogSumExp(log_resp / temp_use) * temp_use
    }
    loglik <- sum(loglik_i)
    loglik_trace[iter] <- loglik

    log_resp_t <- log_resp / temp_use
    lse_t <- rowLogSumExp(log_resp_t)
    tau <- exp(sweep(log_resp_t, 1, lse_t, "-"))
    tau <- pmax(tau, 1e-12)
    tau <- tau / rowSums(tau)

    if (verbose) cat(sprintf("Iter %d: loglik = %.3f\n", iter, loglik))

    if (iter > 1) {
      rel <- abs(loglik - loglik_trace[iter - 1]) / (abs(loglik_trace[iter - 1]) + 1e-8)
      if (rel < tol) {
        if (verbose) cat("Converged.\n")
        loglik_trace <- loglik_trace[seq_len(iter)]
        break
      }
    }

    ## M-step: mixing proportions
    pi_c <- colMeans(tau)
    pi_c <- pmax(pi_c, pi_floor)
    pi_c <- pi_c / sum(pi_c)

    ## M-step: GMM
    for (k in seq_len(K)) {
      w <- tau[, k]
      wsum <- sum(w)
      if (!is.finite(wsum) || wsum < 1e-8) next
      mu_gmm[, k] <- colSums(X_gmm_s * w) / wsum
      xc <- sweep(X_gmm_s, 2, mu_gmm[, k], "-")
      S <- (t(xc * sqrt(w)) %*% (xc * sqrt(w))) / wsum
      S <- (S + t(S)) / 2
      if (gmm_diag) S <- diag(diag(S), q)
      Sigma_list[[k]] <- as.matrix(S) + diag(gmm_ridge, q)
    }

    ## M-step: Cox
    for (k in seq_len(K)) {
      coxfit[[k]] <- cox_glmnet_weighted(
        time = time,
        status = status,
        X = X_cox,
        w = tau[, k],
        lambda = lambda_use,
        alpha = alpha,
        min_eff_events = min_eff_events,
        fallback_fit = global_coxfit
      )
      beta[, k] <- coxfit[[k]]$beta_orig
    }
  }

  clusterid <- max.col(tau, ties.method = "first")
  eff_events <- colSums(tau * status)

  out <- structure(list(
    tau = tau,
    clusterid = clusterid,
    pi = pi_c,
    gmm_scaler = gmm_scaler,
    mu_gmm = mu_gmm,
    Sigma_list = Sigma_list,
    beta = beta,
    beta_list = lapply(seq_len(K), function(k) beta[, k]),
    coxfit = coxfit,
    global_coxfit = global_coxfit,
    loglik = loglik_trace,
    eff_events = eff_events,
    lambda_used = lambda_use,
    K = K,
    gmm_diag = gmm_diag,
    gmm_ridge = gmm_ridge,
    gmm_log_density_scale = gmm_log_density_scale,
    normalize_gmm_by_dim = normalize_gmm_by_dim,
    surv_weight = surv_weight,
    temp = temp_use,
    min_eff_events = min_eff_events,
    init_method = init_method,
    init_eta_weight = init_eta_weight,
    call = match.call()
  ), class = "gemcox_fit")
  finalize_gemcox_fit(out, X_gmm, X_cox, time, status, iter, tol, max_iter)
}

gemcox_full_multistart <- function(...,
                                   n_starts = 10,
                                   init_seeds = NULL,
                                   verbose = FALSE) {
  if (is.null(init_seeds)) init_seeds <- seq_len(n_starts)
  if (length(init_seeds) < n_starts) {
    init_seeds <- rep(init_seeds, length.out = n_starts)
  }

  fits <- vector("list", n_starts)
  final_ll <- rep(NA_real_, n_starts)

  for (s in seq_len(n_starts)) {
    fit_s <- tryCatch(
      gemcox_full(..., init_seed = init_seeds[s], verbose = FALSE),
      error = function(e) NULL
    )
    fits[[s]] <- fit_s
    final_ll[s] <- if (is.null(fit_s)) NA_real_ else fit_s$final_estep_score
  }

  ok <- is.finite(final_ll)
  if (!any(ok)) stop("All multistart fits failed.")

  best_idx <- which.max(final_ll)
  best_fit <- fits[[best_idx]]
  best_fit$multistart <- list(
    n_starts = n_starts,
    init_seeds = init_seeds,
    final_loglik = final_ll,
    best_start = best_idx
  )

  if (verbose) {
    cat(sprintf("Selected multistart solution with final loglik = %.3f (start %d)\n",
                final_ll[best_idx], best_idx))
  }
  best_fit
}

predict_tau <- function(fit, X_gmm_new) {
  X_gmm_new <- as.matrix(X_gmm_new)
  Xg <- apply_x_scaler(X_gmm_new, fit$gmm_scaler)

  K <- length(fit$pi)
  log_resp <- matrix(NA_real_, nrow(Xg), K)
  gmm_scale <- if (!is.null(fit$gmm_log_density_scale)) fit$gmm_log_density_scale else 1

  for (k in seq_len(K)) {
    log_resp[, k] <- log(pmax(fit$pi[k], 1e-12)) +
      gmm_scale * dmvnorm_log(Xg, fit$mu_gmm[, k], fit$Sigma_list[[k]], ridge = fit$gmm_ridge)
  }

  log_resp[!is.finite(log_resp)] <- -700
  lse <- rowLogSumExp(log_resp)
  tau <- exp(sweep(log_resp, 1, lse, "-"))
  tau <- pmax(tau, 1e-12)
  tau / rowSums(tau)
}

shared_breslow_v4 <- function(time, status, eta_list, tau,
                              denom_floor = 1e-8,
                              max_jump = 10,
                              max_cumhaz = 500) {
  status <- as.integer(status)
  uniqT <- sort(unique(time[status == 1]))
  if (length(uniqT) == 0) {
    return(list(time = numeric(0), jump = numeric(0), cumhaz = numeric(0)))
  }
  
  C <- ncol(tau)
  exp_eta_list <- lapply(eta_list, function(eta) exp(clamp(as.numeric(eta), -30, 30)))
  
  jump <- numeric(length(uniqT))
  for (m in seq_along(uniqT)) {
    t_m <- uniqT[m]
    num <- sum(status[time == t_m])
    at_risk <- (time >= t_m)
    denom <- 0
    for (c in seq_len(C)) {
      denom <- denom + sum(tau[at_risk, c] * exp_eta_list[[c]][at_risk])
    }
    denom <- max(denom, denom_floor)
    jump[m] <- min(num / denom, max_jump)
  }
  
  list(time = uniqT, jump = jump, cumhaz = pmin(cumsum(jump), max_cumhaz))
}

cox_poisson_weighted_v4 <- function(time, status, X, w, baseline,
                                    lambda = 0.1,
                                    alpha = 0.5,
                                    min_eff_events = 0.5,
                                    fallback_fit = NULL) {
  X <- as.matrix(X)
  status <- as.integer(status)
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  
  scaler <- make_x_scaler(X)
  Xs <- apply_x_scaler(X, scaler)
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
    return(list(
      beta_s = rep(0, ncol(X)),
      beta_orig = rep(0, ncol(X)),
      scaler = scaler,
      baseline = baseline,
      lambda = NA_real_,
      ok = FALSE,
      reason = "low_events",
      eff_events = eff_events
    ))
  }
  
  Lambda_t <- get_cumhaz_at(baseline, time)
  Lambda_t <- pmax(Lambda_t, 1e-10)
  offset <- log(Lambda_t)
  
  fit <- tryCatch({
    gfit <- glmnet::glmnet(
      x = Xs,
      y = status,
      family = "poisson",
      intercept = FALSE, # shared baseline supplies the common intercept
      weights = w,
      offset = offset,
      alpha = alpha,
      lambda = lambda,
      standardize = FALSE,
      maxit = 1e5
    )
    beta_s <- as.numeric(stats::coef(gfit, s = lambda))[-1]
    list(beta_s = beta_s, ok = TRUE, reason = "ok", lambda = lambda)
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
    fit <- list(beta_s = rep(0, ncol(X)), ok = FALSE, reason = "poisson_failed", lambda = lambda)
  }
  
  beta_s <- fit$beta_s
  beta_s[!is.finite(beta_s)] <- 0
  
  list(
    beta_s = beta_s,
    beta_orig = beta_s / scaler$scale,
    scaler = scaler,
    baseline = baseline,
    lambda = fit$lambda,
    ok = isTRUE(fit$ok),
    reason = fit$reason,
    eff_events = eff_events
  )
}

gemcox_full_shared_v4 <- function(
    X_gmm, X_cox, time, status,
    K = 2,
    lambda = NULL,
    alpha = 0.5,
    max_iter = 60,
    tol = 1e-4,
    verbose = TRUE,
    gmm_ridge = 1e-6,
    gmm_diag = TRUE,
    normalize_gmm_by_dim = FALSE,
    surv_weight = 1.0,
    temp = 1.0,
    pi_floor = 0.02,
    min_eff_events = 1.0,
    init_method = c("supervised", "kmeans"),
    init_eta_weight = 0.75,
    init_seed = 1) {
  
  init_method <- match.arg(init_method)
  X_gmm <- as.matrix(X_gmm)
  X_cox <- as.matrix(X_cox)
  time <- as.numeric(time)
  status <- as.integer(status)
  
  n <- nrow(X_gmm)
  stopifnot(nrow(X_cox) == n, length(time) == n, length(status) == n)
  
  q <- ncol(X_gmm)
  p <- ncol(X_cox)
  gmm_log_density_scale <- if (isTRUE(normalize_gmm_by_dim)) 1 / max(q, 1) else 1
  
  colnames(X_gmm) <- make.names(colnames(X_gmm), unique = TRUE)
  colnames(X_cox) <- make.names(colnames(X_cox), unique = TRUE)
  
  gmm_scaler <- make_x_scaler(X_gmm)
  X_gmm_s <- apply_x_scaler(X_gmm, gmm_scaler)
  
  if (is.null(lambda) || is.character(lambda)) {
    lambda_use <- select_lambda(time, status, X_cox, alpha = alpha)
  } else {
    lambda_use <- as.numeric(lambda)
  }
  
  global_coxfit <- cox_glmnet_weighted(
    time = time,
    status = status,
    X = X_cox,
    w = rep(1, n),
    lambda = lambda_use,
    alpha = alpha,
    min_eff_events = 0
  )
  
  tau <- initialize_tau(
    X_gmm_s = X_gmm_s,
    X_cox = X_cox,
    global_coxfit = global_coxfit,
    K = K,
    init_method = init_method,
    init_eta_weight = init_eta_weight,
    init_seed = init_seed
  )
  
  pi_c <- pmax(colMeans(tau), 1e-8)
  pi_c <- pi_c / sum(pi_c)
  
  mu_gmm <- matrix(0, q, K)
  Sigma_list <- vector("list", K)
  coxfit <- vector("list", K)
  beta <- matrix(0, p, K, dimnames = list(colnames(X_cox), paste0("C", seq_len(K))))
  
  ## Initial GMM moments
  for (k in seq_len(K)) {
    w <- tau[, k]
    wsum <- sum(w)
    if (!is.finite(wsum) || wsum < 1e-8) {
      mu_gmm[, k] <- rep(0, q)
      Sigma_list[[k]] <- diag(1, q)
    } else {
      mu_gmm[, k] <- colSums(X_gmm_s * w) / wsum
      xc <- sweep(X_gmm_s, 2, mu_gmm[, k], "-")
      S <- (t(xc * sqrt(w)) %*% (xc * sqrt(w))) / wsum
      S <- (S + t(S)) / 2
      if (gmm_diag) S <- diag(diag(S), q)
      Sigma_list[[k]] <- as.matrix(S) + diag(gmm_ridge, q)
    }
  }
  
  ## Initial Cox betas from weighted Cox, then convert to shared baseline
  for (k in seq_len(K)) {
    coxfit[[k]] <- cox_glmnet_weighted(
      time = time,
      status = status,
      X = X_cox,
      w = tau[, k],
      lambda = lambda_use,
      alpha = alpha,
      min_eff_events = min_eff_events,
      fallback_fit = global_coxfit
    )
  }
  eta_list <- lapply(seq_len(K), function(k) {
    eta_from_scaled(X_cox, coxfit[[k]]$beta_s, coxfit[[k]]$scaler, eta_cap = 12)
  })
  shared_base <- shared_breslow_v4(time, status, eta_list, tau)
  for (k in seq_len(K)) {
    global_fb <- global_coxfit
    global_fb$baseline <- shared_base
    coxfit[[k]] <- cox_poisson_weighted_v4(
      time = time,
      status = status,
      X = X_cox,
      w = tau[, k],
      baseline = shared_base,
      lambda = lambda_use,
      alpha = alpha,
      min_eff_events = min_eff_events,
      fallback_fit = global_fb
    )
    beta[, k] <- coxfit[[k]]$beta_orig
  }
  
  loglik_trace <- numeric(max_iter)
  temp_use <- if (!is.finite(temp) || temp <= 0) 1.0 else temp
  
  for (iter in seq_len(max_iter)) {
    pi_c <- pmax(pi_c, 1e-8)
    pi_c <- pi_c / sum(pi_c)
    
    ## E-step
    log_resp <- matrix(NA_real_, n, K)
    for (k in seq_len(K)) {
      lgmm <- dmvnorm_log(X_gmm_s, mu_gmm[, k], Sigma_list[[k]], ridge = gmm_ridge)
      lsurv <- log_surv_density(time, status, X_cox, coxfit[[k]])
      log_resp[, k] <- log(pi_c[k]) + gmm_log_density_scale * lgmm + surv_weight * lsurv
    }
    log_resp[!is.finite(log_resp)] <- -700
    
    if (abs(temp_use - 1) < 1e-9) {
      loglik_i <- rowLogSumExp(log_resp)
    } else {
      loglik_i <- rowLogSumExp(log_resp / temp_use) * temp_use
    }
    loglik <- sum(loglik_i)
    loglik_trace[iter] <- loglik
    
    log_resp_t <- log_resp / temp_use
    lse_t <- rowLogSumExp(log_resp_t)
    tau <- exp(sweep(log_resp_t, 1, lse_t, "-"))
    tau <- pmax(tau, 1e-12)
    tau <- tau / rowSums(tau)
    
    if (verbose) cat(sprintf("Iter %d: loglik = %.3f\n", iter, loglik))
    
    if (iter > 1) {
      rel <- abs(loglik - loglik_trace[iter - 1]) / (abs(loglik_trace[iter - 1]) + 1e-8)
      if (rel < tol) {
        if (verbose) cat("Converged.\n")
        loglik_trace <- loglik_trace[seq_len(iter)]
        break
      }
    }
    
    ## M-step: pi
    pi_c <- colMeans(tau)
    pi_c <- pmax(pi_c, pi_floor)
    pi_c <- pi_c / sum(pi_c)
    
    ## M-step: GMM
    for (k in seq_len(K)) {
      w <- tau[, k]
      wsum <- sum(w)
      if (!is.finite(wsum) || wsum < 1e-8) next
      mu_gmm[, k] <- colSums(X_gmm_s * w) / wsum
      xc <- sweep(X_gmm_s, 2, mu_gmm[, k], "-")
      S <- (t(xc * sqrt(w)) %*% (xc * sqrt(w))) / wsum
      S <- (S + t(S)) / 2
      if (gmm_diag) S <- diag(diag(S), q)
      Sigma_list[[k]] <- as.matrix(S) + diag(gmm_ridge, q)
    }
    
    ## M-step: shared baseline then cluster-wise Poisson updates
    eta_list <- lapply(seq_len(K), function(k) {
      eta_from_scaled(X_cox, coxfit[[k]]$beta_s, coxfit[[k]]$scaler, eta_cap = 12)
    })
    shared_base <- shared_breslow_v4(time, status, eta_list, tau)
    
    for (k in seq_len(K)) {
      global_fb <- global_coxfit
      global_fb$baseline <- shared_base
      coxfit[[k]] <- cox_poisson_weighted_v4(
        time = time,
        status = status,
        X = X_cox,
        w = tau[, k],
        baseline = shared_base,
        lambda = lambda_use,
        alpha = alpha,
        min_eff_events = min_eff_events,
        fallback_fit = global_fb
      )
      beta[, k] <- coxfit[[k]]$beta_orig
    }
  }
  
  clusterid <- max.col(tau, ties.method = "first")
  eff_events <- colSums(tau * status)
  
  out <- structure(list(
    tau = tau,
    clusterid = clusterid,
    pi = pi_c,
    gmm_scaler = gmm_scaler,
    mu_gmm = mu_gmm,
    Sigma_list = Sigma_list,
    beta = beta,
    beta_list = lapply(seq_len(K), function(k) beta[, k]),
    coxfit = coxfit,
    global_coxfit = global_coxfit,
    loglik = loglik_trace,
    eff_events = eff_events,
    lambda_used = lambda_use,
    K = K,
    gmm_diag = gmm_diag,
    gmm_ridge = gmm_ridge,
    gmm_log_density_scale = gmm_log_density_scale,
    normalize_gmm_by_dim = normalize_gmm_by_dim,
    surv_weight = surv_weight,
    temp = temp_use,
    min_eff_events = min_eff_events,
    init_method = init_method,
    init_eta_weight = init_eta_weight,
    baseline_mode = "shared",
    call = match.call()
  ), class = "gemcox_fit")
  finalize_gemcox_fit(out, X_gmm, X_cox, time, status, iter, tol, max_iter)
}

gemcox_full_multistart_shared_v4 <- function(...,
                                             n_starts = 10,
                                             init_seeds = NULL,
                                             verbose = FALSE) {
  if (is.null(init_seeds)) init_seeds <- seq_len(n_starts)
  if (length(init_seeds) < n_starts) init_seeds <- rep(init_seeds, length.out = n_starts)
  
  fits <- vector("list", n_starts)
  final_ll <- rep(NA_real_, n_starts)
  
  for (s in seq_len(n_starts)) {
    fit_s <- tryCatch(
      gemcox_full_shared_v4(..., init_seed = init_seeds[s], verbose = FALSE),
      error = function(e) NULL
    )
    fits[[s]] <- fit_s
    final_ll[s] <- if (is.null(fit_s)) NA_real_ else fit_s$final_estep_score
  }
  
  ok <- is.finite(final_ll)
  if (!any(ok)) stop("All multistart fits failed.")
  best_idx <- which.max(final_ll)
  best_fit <- fits[[best_idx]]
  best_fit$multistart <- list(
    n_starts = n_starts,
    init_seeds = init_seeds,
    final_loglik = final_ll,
    best_start = best_idx
  )
  if (verbose) {
    cat(sprintf("Selected shared-baseline multistart solution with final loglik = %.3f (start %d)\n",
                final_ll[best_idx], best_idx))
  }
  best_fit
}

finalize_gemcox_fit <- function(fit, X_gmm, X_cox, time, status,
                                iterations, tol, max_iter) {
  Xg <- apply_x_scaler(as.matrix(X_gmm), fit$gmm_scaler)
  L <- matrix(0, nrow(Xg), fit$K)
  for (k in seq_len(fit$K)) {
    L[, k] <- log(pmax(fit$pi[k], 1e-12)) +
      fit$gmm_log_density_scale *
        dmvnorm_log(Xg, fit$mu_gmm[, k], fit$Sigma_list[[k]], fit$gmm_ridge) +
      fit$surv_weight * log_surv_density(time, status, X_cox, fit$coxfit[[k]])
  }
  if (any(!is.finite(L))) stop("Nonfinite final membership evidence.")
  Lt <- L / fit$temp
  fit$tau <- exp(sweep(Lt, 1, rowLogSumExp(Lt), "-"))
  fit$tau <- pmax(fit$tau, 1e-12)
  fit$tau <- fit$tau / rowSums(fit$tau)
  fit$clusterid <- max.col(fit$tau, ties.method = "first")
  fit$eff_events <- colSums(fit$tau * status)
  fit$final_estep_score <- sum(rowLogSumExp(Lt) * fit$temp)
  fit$iterations <- iterations
  tr <- fit$loglik
  fit$relative_score_change <- if (length(tr) > 1L)
    abs(tail(tr, 1) - tr[length(tr)-1L]) / (abs(tr[length(tr)-1L]) + 1e-8) else Inf
  fit$converged_score <- is.finite(fit$relative_score_change) &&
    fit$relative_score_change < tol
  fit$hit_iteration_limit <- iterations >= max_iter
  fit$engine_patch <- "2026-09-21 scoped audit"
  fit
}

