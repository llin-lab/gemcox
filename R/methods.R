## methods.R -- print / summary / coef / plot methods for class "gemcox".

#' Methods for GeM-Cox fits
#'
#' `print()` shows a compact overview, `summary()` the per-cluster
#' coefficients, mixing proportions, effective events, the coefficient
#' contrast (K = 2), convergence and any numerical guards that bound.
#' `coef()` returns the p x K coefficient matrix on the original feature
#' scale. `plot()` draws the log-likelihood trace, the distribution of the
#' largest membership weight, and the coefficients by cluster.
#'
#' @param x,object A `"gemcox"` fit (or its summary, for printing).
#' @param digits Number of significant digits to print.
#' @param which Which plots to draw: any of `"trace"`, `"tau"`, `"coef"`.
#' @param ... Unused, or passed to plotting functions.
#' @return `print()` and `plot()` return `x` invisibly; `summary()` returns
#'   an object of class `"summary.gemcox"`; `coef()` a p x K matrix.
#' @examples
#' d <- gemcox_simulate(n = 300, p = 4, mu_sep = 1, beta_sep = 2, seed = 1)
#' fit <- gemcox(d$X, time = d$time, status = d$status)
#' summary(fit)
#' coef(fit)
#' plot(fit)
#' @name gemcox-methods
NULL

#' @rdname gemcox-methods
#' @export
coef.gemcox <- function(object, ...) object$beta

#' @rdname gemcox-methods
#' @export
print.gemcox <- function(x, digits = 3, ...) {
  st <- x$settings
  cat("GeM-Cox fit:", st$K, if (st$K == 1) "cluster" else "clusters",
      sprintf("(%s baseline, %s covariance)\n", st$baseline, st$covariance))
  cat(sprintf("  n = %d, events = %d, gamma = %g, lambda = %g, alpha = %g\n",
              length(x$data$time), sum(x$data$status), st$gamma, st$lambda, st$alpha))
  cat(sprintf("  EM: %s after %d iterations; final log-likelihood %.3f\n",
              if (x$converged) "converged" else "NOT converged", x$iterations,
              utils::tail(x$loglik, 1)))
  cat("  mixing proportions:", format(round(x$pi, digits)), "\n\n")
  cat("Coefficients by cluster:\n")
  print(signif(x$beta, digits))
  if (nrow(x$guards)) {
    cat("\nNumerical guards bound during fitting:",
        paste(x$guards$guard, collapse = ", "), "\n")
  }
  invisible(x)
}

#' @rdname gemcox-methods
#' @export
summary.gemcox <- function(object, ...) {
  K <- object$settings$K
  clusters <- data.frame(
    cluster = names(object$pi),
    pi = as.numeric(object$pi),
    weight_sum = colSums(object$tau),
    eff_events = as.numeric(object$eff_events),
    cox_update = object$cox_status$reason,
    row.names = NULL, stringsAsFactors = FALSE)
  contrast <- if (K == 2) object$beta[, 1] - object$beta[, 2] else NULL
  structure(list(
    call = object$call, settings = object$settings, clusters = clusters,
    coefficients = object$beta, contrast = contrast,
    sharpness = mean(apply(object$tau, 1, max)),
    converged = object$converged, iterations = object$iterations,
    loglik = utils::tail(object$loglik, 1), starts = object$starts,
    guards = object$guards,
    n = length(object$data$time), events = sum(object$data$status)),
    class = "summary.gemcox")
}

#' @rdname gemcox-methods
#' @export
print.summary.gemcox <- function(x, digits = 3, ...) {
  st <- x$settings
  cat("Call:\n")
  print(x$call)
  cat(sprintf("\nn = %d, events = %d; K = %d, %s baseline, %s covariance\n",
              x$n, x$events, st$K, st$baseline, st$covariance))
  cat(sprintf("gamma = %g, lambda = %g, alpha = %g, normalize_gmm_by_dim = %s, temp = %g\n",
              st$gamma, st$lambda, st$alpha, st$normalize_gmm_by_dim, st$temp))
  cat(sprintf("EM: %s after %d iterations (tol %g); final log-likelihood %.3f; %d start(s)\n",
              if (x$converged) "converged" else "NOT converged", x$iterations,
              st$tol, x$loglik, nrow(x$starts)))
  cat(sprintf("Mean largest membership weight (sharpness): %.3f\n\n", x$sharpness))
  cat("Clusters:\n")
  cl <- x$clusters
  cl[, c("pi", "weight_sum", "eff_events")] <-
    lapply(cl[, c("pi", "weight_sum", "eff_events")], signif, digits)
  print(cl, row.names = FALSE)
  cat("\nCoefficients (original feature scale):\n")
  cf <- x$coefficients
  if (!is.null(x$contrast)) cf <- cbind(cf, `C1 - C2` = x$contrast)
  print(signif(cf, digits))
  if (nrow(x$guards)) {
    cat("\nNumerical guards that bound during fitting:\n")
    print(x$guards, row.names = FALSE)
  }
  invisible(x)
}

#' @rdname gemcox-methods
#' @export
plot.gemcox <- function(x, which = c("trace", "tau", "coef"), ...) {
  which <- match.arg(which, several.ok = TRUE)
  op <- graphics::par(mfrow = c(1, length(which)))
  on.exit(graphics::par(op))
  if ("trace" %in% which) {
    graphics::plot(seq_along(x$loglik), x$loglik, type = "b", pch = 20,
                   xlab = "EM iteration", ylab = "log-likelihood", main = "EM trace")
  }
  if ("tau" %in% which) {
    graphics::hist(apply(x$tau, 1, max), breaks = 20, xlim = c(0, 1),
                   xlab = "largest membership weight", main = "Membership weights")
  }
  if ("coef" %in% which) {
    graphics::barplot(t(x$beta), beside = TRUE, las = 2,
                      legend.text = colnames(x$beta), main = "Coefficients", ...)
    graphics::abline(h = 0)
  }
  invisible(x)
}
