## lcc.R -- latent-class Cox with logistic gating (a mixture of experts):
## the competitor for the simulation study, added after sim-freeze-v2.
## Internal, not exported; used by the simulation harness only.
##
##   P(Z = k | x)  = softmax_k(a_k0 + a_k' x_gate)          gating
##   h_k(t | x)    = lambda0(t) exp(beta_k' x_cox)          experts
##   E-step: tau_ik propto P(Z = k | x_i) f_k(t_i, d_i | x_i)
##
## It differs from GeM-Cox in one respect only: membership is modelled
## discriminatively (logistic in x) instead of generatively (a Gaussian
## mixture for x). Everything else is GeM-Cox's own machinery, reused: the
## tau-weighted Cox M-step (cox_mstep: shared Breslow baseline and
## weighted Poisson updates with the same ridge), log_surv_density() via
## get_cumhaz_at(), eta_matrix() with its caps, k-means initialisation,
## the convergence rule, multi-start handling and the guard warnings.
## Gating: glmnet multinomial with the membership weights as soft labels,
## unpenalised intercept, ridge (alpha = 0) on the slopes with penalty
## lambda_gate (default: the beta penalty). Gating features are X_gmm,
## standardised with the same scaler as GeM-Cox's clustering features.

# Gating M-step: weighted multinomial logistic regression of the soft
# labels tau on the scaled gating features. Returns intercepts a0 (K) and
# slopes A (q x K) on the scaled scale.
gate_fit <- function(Xg_s, tau, lambda_gate, slopes = TRUE, prev = NULL) {
  K <- ncol(tau)
  q <- ncol(Xg_s)
  if (K == 1) return(list(a0 = 0, A = matrix(0, q, 1), ok = TRUE, reason = "K = 1"))
  if (!slopes) {
    ## intercept-only gating: the MLE is the mean membership weight
    pk <- colMeans(tau)
    return(list(a0 = log(pk) - mean(log(pk)), A = matrix(0, q, K), ok = TRUE,
                reason = "intercept only"))
  }
  Xg <- if (q == 1) cbind(Xg_s, 0) else Xg_s      # glmnet needs >= 2 columns
  g <- tryCatch(glmnet::glmnet(Xg, tau, family = "multinomial", alpha = 0,
                               lambda = lambda_gate, standardize = FALSE,
                               intercept = TRUE, maxit = 1e5),
                error = function(e) NULL)
  if (is.null(g)) {
    if (is.null(prev)) stop("Gating fit failed at initialisation.")
    prev$ok <- FALSE
    prev$reason <- "gating fit failed; previous gating kept"
    return(prev)
  }
  cf <- coef(g, s = lambda_gate)
  a0 <- vapply(cf, function(m) m[1, 1], 0)
  A <- matrix(vapply(cf, function(m) as.numeric(m[-1, 1]), numeric(ncol(Xg))),
              ncol = K)
  if (q == 1) A <- A[1, , drop = FALSE]
  list(a0 = unname(a0), A = unname(A), ok = TRUE, reason = "ok")
}

# log P(Z = k | x) for scaled gating features (n x K).
gate_logprob <- function(gate, Xg_s) {
  L <- sweep(Xg_s %*% gate$A, 2, gate$a0, "+")
  L - rowLogSumExp(L)
}

# E-step evidence: log P(Z = k | x_i) + log f_k(t_i, d_i | x_i).
lcc_log_evidence <- function(Xg_s, time, status, eta_mat, baselines, gate) {
  L <- gate_logprob(gate, Xg_s)
  for (k in seq_len(ncol(L))) {
    L[, k] <- L[, k] + log_surv_density(time, status, eta_mat[, k], baselines[[k]])
  }
  L
}

# The penalties glmnet applies, on the log-likelihood scale: the gating fit
# uses unit weights (n lambda_gate / 2 sum_k ||a_k||^2); each Cox update
# normalises its weights tau_k, so its ridge term is n_k lambda / 2
# ||beta_k||^2 with n_k = sum_i tau_ik of the weights it was fitted with.
lcc_penalty <- function(gate, coxfit, tau_used, st, n) {
  nk <- colSums(tau_used)
  pen_beta <- sum(vapply(seq_along(coxfit), function(k) {
    nk[k] * st$lambda * (1 - st$alpha) / 2 * sum(coxfit[[k]]$beta_s^2)
  }, 0))
  n * st$lambda_gate / 2 * sum(gate$A^2) + pen_beta
}

# One EM run for the latent-class Cox model (structure of em_fit()).
lcc_em <- function(X_gmm, X_cox, time, status, st, seed, verbose = FALSE) {
  K <- st$K
  ctl <- st$control
  n <- nrow(X_gmm)

  gate_scaler <- make_x_scaler(X_gmm)
  Xg_s <- apply_x_scaler(X_gmm, gate_scaler)
  cox_scaler <- make_x_scaler(X_cox)
  Xs <- apply_x_scaler(X_cox, cox_scaler)

  global <- cox_glmnet_weighted(time, status, Xs, rep(1, n), st$lambda, st$alpha,
                                min_eff_events = 0, control = ctl)
  tau <- init_tau(Xg_s, K, st$init, seed)
  gate <- gate_fit(Xg_s, tau, st$lambda_gate, st$gate_slopes)
  coxfit <- lapply(seq_len(K), function(k) {
    cox_glmnet_weighted(time, status, Xs, tau[, k], st$lambda, st$alpha,
                        ctl$min_eff_events, ctl, fallback_fit = global)
  })
  if (st$baseline == "shared") coxfit <- cox_mstep(time, status, Xs, tau, coxfit, global, st)
  tau_used <- tau

  trace <- pen_trace <- numeric(st$max_iter)
  converged <- FALSE
  for (iter in seq_len(st$max_iter)) {
    E <- eta_matrix(Xs, coxfit, ctl$eta_cap)
    L <- sanitize_evidence(lcc_log_evidence(Xg_s, time, status, E,
                                            lapply(coxfit, `[[`, "baseline"), gate))
    trace[iter] <- sum(rowLogSumExp(L))
    pen_trace[iter] <- trace[iter] - lcc_penalty(gate, coxfit, tau_used, st, n)
    tau <- tau_from_evidence(L, 1)
    if (verbose) cat(sprintf("Iter %d: loglik = %.3f\n", iter, trace[iter]))
    if (iter > 1) {
      rel <- abs(trace[iter] - trace[iter - 1]) / (abs(trace[iter - 1]) + 1e-8)
      if (rel < st$tol) {
        converged <- TRUE
        break
      }
    }
    gate <- gate_fit(Xg_s, tau, st$lambda_gate, st$gate_slopes, prev = gate)
    coxfit <- cox_mstep(time, status, Xs, tau, coxfit, global, st)
    tau_used <- tau
  }
  trace <- trace[seq_len(iter)]
  pen_trace <- pen_trace[seq_len(iter)]

  ## final E-step at the returned parameters
  E <- eta_matrix(Xs, coxfit, ctl$eta_cap)
  L <- lcc_log_evidence(Xg_s, time, status, E, lapply(coxfit, `[[`, "baseline"), gate)
  if (any(!is.finite(L))) stop("Non-finite membership evidence at the final parameters.")
  tau <- tau_from_evidence(L, 1)

  cn <- paste0("C", seq_len(K))
  colnames(tau) <- cn
  B_s <- matrix(vapply(coxfit, `[[`, numeric(ncol(Xs)), "beta_s"), ncol(Xs),
                dimnames = list(colnames(X_cox), cn))
  dimnames(gate$A) <- list(colnames(X_gmm), cn)
  names(gate$a0) <- cn
  baselines <- lapply(coxfit, `[[`, "baseline")
  names(baselines) <- cn
  list(
    beta = B_s / cox_scaler$scale,
    beta_scaled = B_s,
    gate = gate,
    gate_slopes_orig = gate$A / gate_scaler$scale,
    baseline = baselines,
    tau = tau,
    cluster = max.col(tau, ties.method = "first"),
    eff_events = stats::setNames(colSums(tau * status), cn),
    loglik = trace,
    pen_loglik = pen_trace,
    final_score = sum(rowLogSumExp(L)),
    converged = converged,
    iterations = iter,
    cox_status = data.frame(cluster = cn,
                            ok = vapply(coxfit, function(cf) isTRUE(cf$ok), NA),
                            reason = vapply(coxfit, `[[`, "", "reason"),
                            stringsAsFactors = FALSE),
    gate_ok = isTRUE(gate$ok),
    global = list(beta = global$beta_s / cox_scaler$scale, baseline = global$baseline),
    gate_scaler = gate_scaler,
    cox_scaler = cox_scaler
  )
}

# Fit the latent-class Cox model (internal; same arguments and validation
# as gemcox() where they apply). Returns an object of class "lcc".
lcc_fit <- function(X_gmm, X_cox = X_gmm, time, status, K = 2, lambda = 0.05, alpha = 0,
                    lambda_gate = lambda, gate_slopes = TRUE,
                    baseline = c("shared", "cluster"), n_starts = 1,
                    max_iter = getOption("gemcox.max_iter", 1000L),
                    tol = getOption("gemcox.tol", 1e-8), init = c("kmeans", "random"), seed = NULL,
                    verbose = FALSE, control = gemcox_control()) {
  cl <- match.call()
  baseline <- match.arg(baseline)
  init <- match.arg(init)
  X_gmm <- as_feature_matrix(X_gmm, "X_gmm")
  X_cox <- as_feature_matrix(X_cox, "X_cox")
  n <- nrow(X_gmm)
  if (nrow(X_cox) != n) {
    stop(sprintf("X_gmm has %d rows but X_cox has %d.", n, nrow(X_cox)), call. = FALSE)
  }
  y <- validate_outcome(time, status, n)
  check_scalar(K, "K", lower = 1, upper = n, integer = TRUE)
  check_scalar(lambda, "lambda", lower = 0)
  check_scalar(lambda_gate, "lambda_gate", lower = 0)
  check_scalar(alpha, "alpha", lower = 0, upper = 1)
  check_scalar(n_starts, "n_starts", lower = 1, integer = TRUE)
  check_scalar(max_iter, "max_iter", lower = 1, integer = TRUE)
  check_scalar(tol, "tol", lower = 0)
  if (is.null(control$min_eff_events)) {
    control$min_eff_events <- if (baseline == "shared") 1 else 0.5
  }
  settings <- list(model = "latent-class Cox (logistic gating)", K = as.integer(K),
                   gamma = 1, lambda = lambda, alpha = alpha, lambda_gate = lambda_gate,
                   gate_alpha = 0, gate_slopes = gate_slopes, baseline = baseline,
                   temp = 1, n_starts = as.integer(n_starts),
                   max_iter = as.integer(max_iter), tol = tol, init = init, seed = seed,
                   control = control)
  seeds <- if (is.null(seed)) seq_len(n_starts) else seed + seq_len(n_starts) - 1L

  run <- collect_guards({
    fits <- vector("list", n_starts)
    errs <- rep(NA_character_, n_starts)
    for (s in seq_len(n_starts)) {
      fits[[s]] <- tryCatch(
        lcc_em(X_gmm, X_cox, y$time, y$status, settings, seeds[s],
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
  out <- c(fits[[best]], list(
    starts = data.frame(start = seq_len(n_starts), seed = seeds, final_score = score,
                        error = run$value$errs, stringsAsFactors = FALSE),
    best_start = best, guards = run$guards, settings = settings,
    data = list(X_gmm = X_gmm, X_cox = X_cox, time = y$time, status = y$status),
    call = cl))
  class(out) <- "lcc"
  warn_guard_summary(run$guards)
  warn_not_converged(out$converged, out$iterations, tol)
  out
}

# Membership probabilities for new subjects from the gating model only.
lcc_predict_tau <- function(object, X_gmm_new) {
  Xg <- match_columns(as_feature_matrix(X_gmm_new, "X_gmm_new"),
                      rownames(object$gate$A), "X_gmm_new")
  P <- exp(gate_logprob(object$gate, apply_x_scaler(Xg, object$gate_scaler)))
  colnames(P) <- colnames(object$beta)
  P
}
