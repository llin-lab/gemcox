## predict.R -- predict.gemcox(): tau from X_gmm only, linear predictors,
## and mixture survival curves. This is the single path for computing
## membership weights on new data.

#' Predict from a GeM-Cox fit
#'
#' Membership weights for new subjects are computed from the clustering
#' features `newX_gmm` only, because a new subject's outcome is unknown at
#' prediction time. They therefore differ from the training weights
#' `object$tau` whenever `gamma > 0`.
#'
#' @param object A `"gemcox"` fit.
#' @param newX_gmm Clustering features for new subjects (columns matched to
#'   the training `X_gmm` by name).
#' @param newX_cox Survival features for new subjects. May be omitted when
#'   the model was fitted with `X_cox = X_gmm`.
#' @param type `"tau"`: n x K membership weights. `"lp"`: the mixture linear
#'   predictor `sum_k tau_k (x - xbar)' beta_k`, centred at the training
#'   means of `X_cox` so that it pairs with the fitted baseline (see
#'   "Baseline and centring" in [gemcox()]). The uncentred version is
#'   `rowSums(predict(object, newX_gmm) * (newX_cox %*% coef(object)))`.
#'   `"survival"`:
#'   an n x length(times) matrix of mixture survival probabilities
#'   `sum_k tau_k S_k(t | x)`.
#' @param times Times at which to evaluate survival (for `type = "survival"`).
#' @param ... Unused.
#' @return A matrix (types `"tau"` and `"survival"`) or numeric vector
#'   (type `"lp"`).
#' @examples
#' d <- gemcox_simulate(n = 300, p = 4, mu_sep = 2, beta_sep = 2, seed = 2)
#' fit <- gemcox(d$X, time = d$time, status = d$status)
#' new <- gemcox_simulate(n = 5, p = 4, mu_sep = 2, beta_sep = 2, seed = 3)
#' predict(fit, new$X)
#' predict(fit, new$X, type = "lp")
#' predict(fit, new$X, type = "survival", times = c(50, 100))
#' @export
predict.gemcox <- function(object, newX_gmm, newX_cox,
                           type = c("tau", "lp", "survival"), times, ...) {
  type <- match.arg(type)
  if (missing(newX_gmm)) stop("newX_gmm is required.", call. = FALSE)
  Xg <- match_columns(as_feature_matrix(newX_gmm, "newX_gmm"),
                      rownames(object$mu), "newX_gmm")
  tau <- predict_tau(object, Xg)
  if (type == "tau") return(tau)

  if (missing(newX_cox)) {
    if (!isTRUE(object$x_cox_is_x_gmm)) {
      stop("newX_cox is required: the model was fitted with X_cox different from X_gmm.",
           call. = FALSE)
    }
    newX_cox <- newX_gmm
  }
  Xc <- match_columns(as_feature_matrix(newX_cox, "newX_cox"),
                      rownames(object$beta), "newX_cox")
  if (nrow(Xc) != nrow(Xg)) {
    stop("newX_gmm and newX_cox must have the same number of rows.", call. = FALSE)
  }
  E <- predict_eta(object, Xc)

  if (type == "lp") return(as.numeric(rowSums(tau * E)))

  if (missing(times)) stop("times is required for type = \"survival\".", call. = FALSE)
  tt <- as.numeric(times)
  S <- matrix(0, nrow(Xc), length(tt))
  for (k in seq_len(ncol(tau))) {
    Lambda_k <- get_cumhaz_at(object$baseline[[k]], tt)
    S <- S + tau[, k] * exp(-tcrossprod(exp(E[, k]), Lambda_k))
  }
  colnames(S) <- paste0("t=", signif(tt, 4))
  S
}

# Membership weights from the Gaussian mixture part only.
predict_tau <- function(object, Xg) {
  Xs <- apply_x_scaler(Xg, object$gmm_scaler)
  st <- object$settings
  L <- gmm_log_evidence(Xs, object$pi, object$mu, object$Sigma, st$control$gmm_ridge,
                        st$gmm_log_density_scale, pi_min = 1e-12)
  L <- sanitize_evidence(L)
  lse <- rowLogSumExp(L)
  tau <- exp(sweep(L, 1, lse, "-"))
  tau <- pmax(tau, 1e-12)
  tau <- tau / rowSums(tau)
  colnames(tau) <- names(object$pi)
  tau
}

# Per-cluster linear predictors (n x K), centred and capped as in fitting.
predict_eta <- function(object, Xc) {
  Xs <- apply_x_scaler(Xc, object$cox_scaler)
  E <- Xs %*% object$beta_scaled
  matrix(clamp_eta(as.numeric(E), object$settings$control$eta_cap, "eta_cap"),
         nrow(Xs), dimnames = list(NULL, colnames(object$beta)))
}
