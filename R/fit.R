## fit.R -- gemcox(): the single user-facing fitter.
##
## Ported from GeMCox_combined_PATCHED.R:
##   gemcox_full_shared_v4()             -> baseline = "shared"
##   gemcox_full()                       -> baseline = "cluster"
##   gemcox_full_multistart_shared_v4()  -> n_starts > 1
## The EM iteration order, numerical floors and convergence rule are
## unchanged; tests/testthat/test-equivalence.R checks this against a
## verbatim copy of the PATCHED functions.

#' Numerical controls for gemcox()
#'
#' Low-level settings that rarely need changing. Each clipping bound emits a
#' warning when it binds, because silent clipping can attenuate coefficients.
#'
#' @param gmm_ridge Ridge added to each Gaussian covariance.
#' @param pi_floor Lower bound on mixing proportions after each M-step.
#' @param min_eff_events Minimum tau-weighted events for a cluster's Cox
#'   update; below it the cluster falls back to the pooled (K = 1) fit. `NULL`
#'   uses the legacy value for the chosen baseline (1 for `"shared"`, 0.5 for
#'   `"cluster"`).
#' @param max_jump,max_cumhaz Caps on a Breslow hazard jump and on the
#'   cumulative hazard.
#' @param denom_floor Floor on the Breslow risk-set denominator.
#' @param eta_cap Cap on |linear predictor| for survival terms and prediction.
#' @param eta_clamp Cap on |linear predictor| inside the Breslow estimators.
#' @return A named list of settings, for the `control` argument of [gemcox()].
#' @examples
#' gemcox_control(max_jump = 20)
#' @export
gemcox_control <- function(gmm_ridge = 1e-6, pi_floor = 0.02, min_eff_events = NULL,
                           max_jump = 10, max_cumhaz = 500, denom_floor = 1e-8,
                           eta_cap = 12, eta_clamp = 30) {
  check_scalar(gmm_ridge, "gmm_ridge", lower = 0)
  check_scalar(pi_floor, "pi_floor", lower = 0, upper = 0.5)
  if (!is.null(min_eff_events)) check_scalar(min_eff_events, "min_eff_events", lower = 0)
  check_scalar(max_jump, "max_jump", lower = 0, lower_open = TRUE)
  check_scalar(max_cumhaz, "max_cumhaz", lower = 0, lower_open = TRUE)
  check_scalar(denom_floor, "denom_floor", lower = 0)
  check_scalar(eta_cap, "eta_cap", lower = 0, lower_open = TRUE)
  check_scalar(eta_clamp, "eta_clamp", lower = 0, lower_open = TRUE)
  list(gmm_ridge = gmm_ridge, pi_floor = pi_floor, min_eff_events = min_eff_events,
       max_jump = max_jump, max_cumhaz = max_cumhaz, denom_floor = denom_floor,
       eta_cap = eta_cap, eta_clamp = eta_clamp)
}

#' Fit a Gaussian mixture of Cox models (GeM-Cox)
#'
#' Fits K latent clusters. Within cluster k the clustering features follow
#' `N(mu_k, Sigma_k)` and survival follows a Cox model with coefficients
#' `beta_k`. By default all clusters share one Breslow baseline hazard, so
#' clusters differ only in their survival *mechanism* (`beta_k`) and in
#' their feature profile (`mu_k`, `Sigma_k`).
#'
#' The EM algorithm alternates:
#' * **E-step.** `tau_ik` is proportional to
#'   `pi_k N(x_gmm_i | mu_k, Sigma_k) f_k(t_i, d_i | x_cox_i)^gamma`.
#' * **M-step.** Closed-form `pi`, `mu`, `Sigma`; the shared baseline by
#'   weighted Breslow; then each `beta_k` by a tau-weighted ridge fit of the
#'   Poisson working model with offset `log Lambda0(t)`.
#'
#' **What gamma does.** `gamma` sets how much a subject's own outcome
#' counts when deciding which cluster they belong to. With `gamma = 0`,
#' membership comes from the clustering features alone, so clusters that
#' differ only in mechanism (and not in feature profile) cannot be found.
#' With `gamma = 1` (the generative model, and the default) the outcome is
#' weighted as in the joint likelihood. `gamma` is fixed rather than tuned:
#' data-driven selectors either chose `gamma = 0` or inflated false
#' discovery in simulation.
#'
#' **Baseline and centring.** Cox features are centred at their training
#' means before fitting, so the baseline hazard is the hazard of a subject
#' whose `X_cox` equals the training means, and cluster k's hazard is
#' `lambda0(t) exp((x - xbar)' beta_k)`. With a shared baseline this matters:
#' because `xbar' beta_k` differs between clusters, the centred and uncentred
#' mixture predictors differ by a subject-specific amount, and only the
#' centred one pairs with the fitted baseline. The means are stored in
#' `cox_scaler$center`. The overall level of the linear predictor is not
#' identified separately from the baseline, so compare linear predictors
#' only up to a constant.
#'
#' **What is estimable.** When clusters differ only in mechanism, individual
#' labels are not consistently recoverable, but the coefficient contrast
#' `beta_1 - beta_2` is. Test K = 2 against K = 1 with
#' [gemcox_heterogeneity_test()], not by cross-validated C-index.
#'
#' @param X_gmm Numeric matrix or data frame of clustering features (n x q).
#'   Column names are kept exactly as given.
#' @param X_cox Numeric matrix or data frame of survival features (n x p).
#'   Defaults to `X_gmm`.
#' @param time Non-negative follow-up times.
#' @param status Event indicator, 1 = event, 0 = censored.
#' @param K Number of clusters.
#' @param gamma Survival weight in the E-step (see Details). Default 1.
#' @param lambda Ridge/elastic-net penalty for the Cox updates, on the glmnet
#'   scale with standardized features. Default 0.05.
#' @param alpha Elastic-net mixing (0 = ridge). Default 0.
#' @param baseline `"shared"` (one Breslow baseline, default) or `"cluster"`
#'   (a separate baseline per cluster).
#' @param covariance `"diagonal"` (default) or `"full"` Gaussian covariance.
#'   Full covariance is rarely estimable at small n.
#' @param normalize_gmm_by_dim If `TRUE`, the Gaussian log-density is divided
#'   by q in the E-step, which softens the membership weights. Its effect is
#'   regime-dependent; the default is `FALSE`. In simulation, `TRUE`
#'   recovered the coefficient contrast much better when feature profiles
#'   barely differed, and worse once they clearly differed. A rule choosing
#'   between them from the data did not meet the pre-specified criterion for
#'   becoming the default (see NEWS).
#' @param temp Temperature of the E-step softmax (1 = standard EM).
#' @param n_starts Number of EM starts; the start with the highest final
#'   E-step score is kept.
#' @param max_iter,tol Maximum EM iterations, and the convergence tolerance:
#'   EM stops when the relative change in the log-likelihood trace falls
#'   below `tol`. The defaults are `tol = 1e-8` and `max_iter = 1000`
#'   (since version 0.2.0; set with `options(gemcox.tol = , gemcox.max_iter = )`).
#'   - In a tolerance study, looser settings often stopped well before the
#'     optimum. At the earlier default (1e-4, 100 iterations), only 32% of
#'     fits gave a coefficient contrast within 0.01 of the 1e-8 fit.
#'   - The tight setting costs roughly 7 times more fitting time.
#'   - Results of versions before 0.2.0, including the frozen simulation
#'     study, used `tol = 1e-4, max_iter = 100`. To reproduce them, set
#'     `options(gemcox.tol = 1e-4, gemcox.max_iter = 100)`.
#'   - A fit that stops at `max_iter` without converging gives a warning.
#' @param init Initial clustering: `"kmeans"` on the scaled `X_gmm` (default),
#'   or `"random"` balanced random labels. Neither uses the outcome, so the
#'   `gamma = 0` fit stays a purely covariate-based comparator.
#' @param seed Seed for initialisation. Start s uses `seed + s - 1`; `NULL`
#'   uses seeds `1, ..., n_starts`. The caller's random number stream is left
#'   unchanged.
#' @param verbose Print the log-likelihood at each iteration.
#' @param control Numerical settings from [gemcox_control()].
#'
#' @return An object of class `"gemcox"`, a list including:
#' \describe{
#'   \item{beta}{p x K coefficients on the original feature scale.}
#'   \item{pi}{Mixing proportions.}
#'   \item{mu, Sigma}{Cluster means (q x K) and list of covariances, on the
#'     standardized `X_gmm` scale (see `gmm_scaler`).}
#'   \item{baseline}{List of K step functions `list(time, jump, cumhaz)`;
#'     identical across clusters when `baseline = "shared"`. They pair with
#'     linear predictors centred at the training means of `X_cox`.}
#'   \item{tau}{n x K membership weights, recomputed from the final
#'     parameters (they use the outcome when `gamma > 0`).}
#'   \item{loglik}{Trace of the E-step objective; with `gamma = 1` and
#'     `temp = 1` this is the observed-data log-likelihood.}
#'   \item{converged, iterations}{Convergence flag and iteration count.}
#'   \item{gmm_scaler, cox_scaler}{Centres and scales of the training
#'     features, needed to score new data.}
#'   \item{starts, guards}{Per-start results, and any numerical guards that
#'     bound during fitting.}
#'   \item{call, settings}{The call and all settings used.}
#' }
#' @seealso [predict.gemcox()], [gemcox_heterogeneity_test()],
#'   [gemcox_simulate()].
#' @examples
#' d <- gemcox_simulate(n = 300, p = 4, mu_sep = 1, beta_sep = 2, seed = 1)
#' fit <- gemcox(d$X, time = d$time, status = d$status, K = 2)
#' fit
#' coef(fit)
#' @export
gemcox <- function(X_gmm, X_cox = X_gmm, time, status, K = 2,
                   gamma = 1, lambda = 0.05, alpha = 0,
                   baseline = c("shared", "cluster"),
                   covariance = c("diagonal", "full"),
                   normalize_gmm_by_dim = FALSE, temp = 1,
                   n_starts = 1, max_iter = getOption("gemcox.max_iter", 1000L),
                   tol = getOption("gemcox.tol", 1e-8),
                   init = c("kmeans", "random"), seed = NULL, verbose = FALSE,
                   control = gemcox_control()) {
  cl <- match.call()
  baseline <- match.arg(baseline)
  covariance <- match.arg(covariance)
  init <- match.arg(init)

  X_gmm <- as_feature_matrix(X_gmm, "X_gmm")
  X_cox <- as_feature_matrix(X_cox, "X_cox")
  n <- nrow(X_gmm)
  if (nrow(X_cox) != n) {
    stop(sprintf("X_gmm has %d rows but X_cox has %d.", n, nrow(X_cox)), call. = FALSE)
  }
  y <- validate_outcome(time, status, n)
  check_scalar(K, "K", lower = 1, upper = n, integer = TRUE)
  check_scalar(gamma, "gamma", lower = 0)
  check_scalar(lambda, "lambda", lower = 0)
  check_scalar(alpha, "alpha", lower = 0, upper = 1)
  check_scalar(temp, "temp", lower = 0, lower_open = TRUE)
  check_scalar(n_starts, "n_starts", lower = 1, integer = TRUE)
  check_scalar(max_iter, "max_iter", lower = 1, integer = TRUE)
  check_scalar(tol, "tol", lower = 0)
  if (!is.null(seed)) check_scalar(seed, "seed", integer = TRUE)
  if (!is.logical(normalize_gmm_by_dim) || length(normalize_gmm_by_dim) != 1 ||
      is.na(normalize_gmm_by_dim)) {
    stop("normalize_gmm_by_dim must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(control$min_eff_events)) {
    control$min_eff_events <- if (baseline == "shared") 1 else 0.5
  }

  settings <- list(K = as.integer(K), gamma = gamma, lambda = lambda, alpha = alpha,
                   baseline = baseline, covariance = covariance,
                   normalize_gmm_by_dim = normalize_gmm_by_dim,
                   gmm_log_density_scale = if (normalize_gmm_by_dim) 1 / ncol(X_gmm) else 1,
                   temp = temp, n_starts = as.integer(n_starts), max_iter = as.integer(max_iter),
                   tol = tol, init = init, seed = seed, control = control)
  seeds <- if (is.null(seed)) seq_len(n_starts) else seed + seq_len(n_starts) - 1L

  run <- collect_guards({
    fits <- vector("list", n_starts)
    errs <- rep(NA_character_, n_starts)
    for (s in seq_len(n_starts)) {
      fits[[s]] <- tryCatch(
        em_fit(X_gmm, X_cox, y$time, y$status, settings, seeds[s],
               verbose = verbose && n_starts == 1),
        error = function(e) {
          if (n_starts == 1) stop(e)
          errs[s] <<- conditionMessage(e)
          NULL
        })
    }
    list(fits = fits, errs = errs)
  })
  fits <- run$value$fits
  score <- vapply(fits, function(f) if (is.null(f)) NA_real_ else f$final_score, 0)
  if (!any(is.finite(score))) {
    stop("All EM starts failed: ", paste(unique(stats::na.omit(run$value$errs)), collapse = "; "),
         call. = FALSE)
  }
  best <- which.max(score)
  fit <- fits[[best]]

  out <- c(fit, list(
    starts = data.frame(start = seq_len(n_starts), seed = seeds, final_score = score,
                        converged = vapply(fits, function(f) isTRUE(f$converged), NA),
                        error = run$value$errs, stringsAsFactors = FALSE),
    best_start = best,
    guards = run$guards,
    settings = settings,
    data = list(X_gmm = X_gmm, X_cox = X_cox, time = y$time, status = y$status),
    x_cox_is_x_gmm = identical(X_gmm, X_cox),
    call = cl))
  class(out) <- "gemcox"
  warn_guard_summary(run$guards)
  warn_not_converged(out$converged, out$iterations, tol)
  if (verbose && n_starts > 1) {
    message(sprintf("Kept start %d of %d (final score %.3f).", best, n_starts, score[best]))
  }
  out
}

# Initial membership weights (never outcome-informed).
init_tau <- function(X_gmm_s, K, init, seed) {
  n <- nrow(X_gmm_s)
  if (K == 1) return(matrix(1, n, 1, dimnames = list(NULL, "C1")))
  Z0 <- with_local_seed(seed, {
    if (init == "kmeans") {
      kmeans(X_gmm_s, centers = K, nstart = 50)$cluster
    } else {
      sample(rep(seq_len(K), length.out = n))
    }
  })
  tau <- matrix(0, n, K, dimnames = list(NULL, paste0("C", seq_len(K))))
  tau[cbind(seq_len(n), Z0)] <- 1
  tau
}

# Linear predictors (n x K), centred at the training means and capped.
eta_matrix <- function(Xs, coxfit, eta_cap) {
  B <- vapply(coxfit, function(cf) cf$beta_s, numeric(ncol(Xs)))
  B <- matrix(B, nrow = ncol(Xs))
  E <- Xs %*% B
  matrix(clamp_eta(as.numeric(E), eta_cap, "eta_cap"), nrow(Xs))
}

# Cox M-step for all clusters, for either baseline mode.
cox_mstep <- function(time, status, Xs, tau, coxfit, global, st) {
  K <- ncol(tau)
  ctl <- st$control
  if (st$baseline == "cluster") {
    for (k in seq_len(K)) {
      coxfit[[k]] <- cox_glmnet_weighted(time, status, Xs, tau[, k], st$lambda,
                                         st$alpha, ctl$min_eff_events, ctl,
                                         fallback_fit = global)
    }
    return(coxfit)
  }
  ## shared: baseline from the current betas, then Poisson updates per cluster
  E <- eta_matrix(Xs, coxfit, ctl$eta_cap)
  base <- shared_breslow(time, status, lapply(seq_len(K), function(k) E[, k]), tau, ctl)
  global_fb <- global
  global_fb$baseline <- base
  for (k in seq_len(K)) {
    coxfit[[k]] <- cox_poisson_weighted(time, status, Xs, tau[, k], base, st$lambda,
                                        st$alpha, ctl$min_eff_events,
                                        fallback_fit = global_fb)
  }
  coxfit
}

# One EM run from one initialisation.
em_fit <- function(X_gmm, X_cox, time, status, st, seed, verbose = FALSE) {
  K <- st$K
  ctl <- st$control
  n <- nrow(X_gmm)
  q <- ncol(X_gmm)
  diagonal <- st$covariance == "diagonal"

  gmm_scaler <- make_x_scaler(X_gmm)
  X_gmm_s <- apply_x_scaler(X_gmm, gmm_scaler)
  cox_scaler <- make_x_scaler(X_cox)
  Xs <- apply_x_scaler(X_cox, cox_scaler)

  global <- cox_glmnet_weighted(time, status, Xs, rep(1, n), st$lambda, st$alpha,
                                min_eff_events = 0, control = ctl)
  tau <- init_tau(X_gmm_s, K, st$init, seed)
  pi_k <- pmax(colMeans(tau), 1e-8)
  pi_k <- pi_k / sum(pi_k)

  g <- update_gmm(X_gmm_s, tau, matrix(0, q, K), vector("list", K), diagonal,
                  ctl$gmm_ridge, initial = TRUE)
  mu <- g$mu
  Sigma <- g$Sigma

  ## initial Cox fits: weighted glmnet Cox (then, if shared, one shared update)
  coxfit <- lapply(seq_len(K), function(k) {
    cox_glmnet_weighted(time, status, Xs, tau[, k], st$lambda, st$alpha,
                        ctl$min_eff_events, ctl, fallback_fit = global)
  })
  if (st$baseline == "shared") coxfit <- cox_mstep(time, status, Xs, tau, coxfit, global, st)

  trace <- numeric(st$max_iter)
  converged <- FALSE
  for (iter in seq_len(st$max_iter)) {
    pi_k <- pmax(pi_k, 1e-8)
    pi_k <- pi_k / sum(pi_k)

    ## E-step
    E <- eta_matrix(Xs, coxfit, ctl$eta_cap)
    L <- estep_log_evidence(X_gmm_s, time, status, E, lapply(coxfit, `[[`, "baseline"),
                            pi_k, mu, Sigma, ctl$gmm_ridge, st$gmm_log_density_scale,
                            st$gamma)
    L <- sanitize_evidence(L)
    trace[iter] <- evidence_score(L, st$temp)
    tau <- tau_from_evidence(L, st$temp)
    if (verbose) cat(sprintf("Iter %d: loglik = %.3f\n", iter, trace[iter]))

    if (iter > 1) {
      rel <- abs(trace[iter] - trace[iter - 1]) / (abs(trace[iter - 1]) + 1e-8)
      if (rel < st$tol) {
        converged <- TRUE
        break
      }
    }

    ## M-step
    pi_k <- update_pi(tau, ctl$pi_floor)
    g <- update_gmm(X_gmm_s, tau, mu, Sigma, diagonal, ctl$gmm_ridge)
    mu <- g$mu
    Sigma <- g$Sigma
    coxfit <- cox_mstep(time, status, Xs, tau, coxfit, global, st)
  }
  trace <- trace[seq_len(iter)]

  ## Final E-step at the returned parameters, so tau matches them.
  E <- eta_matrix(Xs, coxfit, ctl$eta_cap)
  L <- estep_log_evidence(X_gmm_s, time, status, E, lapply(coxfit, `[[`, "baseline"),
                          pi_k, mu, Sigma, ctl$gmm_ridge, st$gmm_log_density_scale,
                          st$gamma, pi_min = 1e-12)
  if (any(!is.finite(L))) stop("Non-finite membership evidence at the final parameters.")
  tau <- tau_from_evidence(L, st$temp)

  cn <- paste0("C", seq_len(K))
  colnames(tau) <- cn
  B_s <- matrix(vapply(coxfit, `[[`, numeric(ncol(Xs)), "beta_s"), ncol(Xs),
                dimnames = list(colnames(X_cox), cn))
  dimnames(mu) <- list(colnames(X_gmm), cn)
  names(pi_k) <- cn
  Sigma <- lapply(Sigma, function(S) {
    dimnames(S) <- list(colnames(X_gmm), colnames(X_gmm))
    S
  })
  names(Sigma) <- cn
  baselines <- lapply(coxfit, `[[`, "baseline")
  names(baselines) <- cn

  list(
    beta = B_s / cox_scaler$scale,
    beta_scaled = B_s,
    pi = pi_k,
    mu = mu,
    Sigma = Sigma,
    baseline = baselines,
    tau = tau,
    cluster = max.col(tau, ties.method = "first"),
    eff_events = stats::setNames(colSums(tau * status), cn),
    loglik = trace,
    final_score = sum(rowLogSumExp(L / st$temp) * st$temp),
    converged = converged,
    iterations = iter,
    cox_status = data.frame(cluster = cn,
                            ok = vapply(coxfit, function(cf) isTRUE(cf$ok), NA),
                            reason = vapply(coxfit, `[[`, "", "reason"),
                            stringsAsFactors = FALSE),
    global = list(beta = global$beta_s / cox_scaler$scale, baseline = global$baseline),
    gmm_scaler = gmm_scaler,
    cox_scaler = cox_scaler
  )
}
