## estep.R -- E-step: membership weights tau_ik.
##
##   log e_ik = log pi_k + s * log N(x_gmm_i | mu_k, Sigma_k) + gamma * log f_k(t_i, d_i | x_cox_i)
##   tau_ik   = softmax_k(log e_ik / temp)
##
## s is 1, or 1/q when normalize_gmm_by_dim = TRUE. With gamma = 0 the
## weights use the clustering features only.

# Evaluate log phi(X; mu, Sigma) for each row of X.
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

# n x K matrix of log pi_k + s * log N(x | mu_k, Sigma_k), on the scaled GMM features.
gmm_log_evidence <- function(X_gmm_s, pi_k, mu, Sigma, ridge, gmm_scale,
                             pi_min = 0) {
  K <- length(pi_k)
  L <- matrix(NA_real_, nrow(X_gmm_s), K)
  for (k in seq_len(K)) {
    L[, k] <- log(pmax(pi_k[k], pi_min)) +
      gmm_scale * dmvnorm_log(X_gmm_s, mu[, k], Sigma[[k]], ridge = ridge)
  }
  L
}

# Full E-step evidence: GMM part plus gamma times the survival log-density.
estep_log_evidence <- function(X_gmm_s, time, status, eta_mat, baselines,
                               pi_k, mu, Sigma, ridge, gmm_scale, gamma,
                               pi_min = 0) {
  L <- gmm_log_evidence(X_gmm_s, pi_k, mu, Sigma, ridge, gmm_scale, pi_min)
  for (k in seq_len(ncol(L))) {
    L[, k] <- L[, k] + gamma * log_surv_density(time, status, eta_mat[, k], baselines[[k]])
  }
  L
}

# Replace non-finite evidence by -700 (warning when it happens).
sanitize_evidence <- function(L) {
  bad <- !is.finite(L)
  check_guard(bad, "nonfinite_evidence")
  L[bad] <- -700
  L
}

# Tempered softmax over clusters, floored at 1e-12 and renormalised.
tau_from_evidence <- function(L, temp) {
  Lt <- L / temp
  tau <- exp(sweep(Lt, 1, rowLogSumExp(Lt), "-"))
  tau <- pmax(tau, 1e-12)
  tau / rowSums(tau)
}

# Objective recorded in the log-likelihood trace. With gamma = 1 and
# temp = 1 this is the observed-data log-likelihood.
evidence_score <- function(L, temp) {
  if (abs(temp - 1) < 1e-9) sum(rowLogSumExp(L)) else sum(rowLogSumExp(L / temp) * temp)
}
