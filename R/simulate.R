## simulate.R -- the one data-generating mechanism.
##
## Ported verbatim from sim_data() in the simulation pipeline's 00_core.R,
## so that the package examples, gemcox_selftest() and the pipeline share a
## single DGP. The pipeline calls this function; nothing redefines it.

#' Simulate data from the two-cluster GeM-Cox design
#'
#' Two latent clusters with equal probability:
#' * `X | Z = k ~ N(mu_k, (1 - rho) I + rho J)` with `mu_k = -/+ (mu_sep / 2) u`;
#' * `beta_k = base -/+ (beta_sep / 2) v`;
#' * `T | X, Z` Weibull (shape 1.5) with Cox linear predictor `x' beta_k`;
#' * independent exponential censoring calibrated to `event_rate`, times
#'   floored at 0.1.
#'
#' `u` (profile direction) and `v` (coefficient direction) are orthonormal
#' and `base` is orthogonal to `v`. So `mu_sep` controls only how far apart
#' the feature profiles sit, `beta_sep` only how different the survival
#' mechanisms are, and the two clusters have equal coefficient norms.
#'
#' @param n Number of subjects.
#' @param p Number of features (at least 1; `u` and `v` need p >= 4 to be
#'   distinct as designed).
#' @param mu_sep Distance between cluster means along `u`.
#' @param beta_sep Distance between cluster coefficient vectors along `v`.
#' @param rho Exchangeable feature correlation.
#' @param K_true 1, 2 or 3 clusters. With 1, all subjects share `base` and mean
#'   0. With 3 (added after the frozen study; needs `p >= 8`), clusters sit at
#'   the vertices of equilateral triangles with side `mu_sep` in profile
#'   space and `beta_sep` in coefficient space, in two further orthogonal
#'   directions, so that K = 3 cannot be reduced to two collinear groups.
#' @param event_rate Target proportion of events (random censoring only).
#' @param censoring `"random"` (default, the frozen design): independent
#'   exponential censoring calibrated per dataset to `event_rate`.
#'   `"administrative"` (added after sim-freeze-v2 for experiment E6): every
#'   subject is followed up to the same time `followup`, and there is no
#'   other censoring. The event rate is then set by `followup`, not by
#'   `event_rate`.
#' @param followup Fixed follow-up time for `censoring = "administrative"`
#'   (`Inf` gives no censoring).
#' @param features `"gaussian"` (the design above; the default) or
#'   `"lognormal"` (added after sim-freeze-v2 for the misspecified-features
#'   experiment C5). With `"lognormal"`, each coordinate of the Gaussian
#'   noise `z` is replaced by the standardised log-normal
#'   `(exp(sdlog z) - m) / s`, with `m` and `s` the log-normal mean and
#'   standard deviation, so every feature has mean 0 and variance 1 about its
#'   cluster mean but is right-skewed. The cluster means, directions,
#'   coefficients and survival model are unchanged, and the same random
#'   numbers are drawn. Requires `rho = 0`.
#' @param sdlog Log-scale standard deviation for `features = "lognormal"`.
#' @param seed Random seed. The caller's random number stream is left
#'   unchanged.
#' @return A list with `X` (n x p, columns `f1..fp`), `time`, `status`, true
#'   cluster means `mu_list` (added after sim-freeze-v2), true
#'   labels `Z`, true linear predictor `eta`, `beta_list`, and `p`.
#' @examples
#' d <- gemcox_simulate(n = 200, p = 5, mu_sep = 0, beta_sep = 2, seed = 1)
#' mean(d$status)
#' @export
gemcox_simulate <- function(n, p, mu_sep, beta_sep, rho = 0, K_true = 2,
                            event_rate = 0.55, features = c("gaussian", "lognormal"),
                            sdlog = 0.8, censoring = c("random", "administrative"),
                            followup = NULL, seed) {
  features <- match.arg(features)
  censoring <- match.arg(censoring)
  if (censoring == "administrative") {
    if (is.null(followup) || length(followup) != 1 || is.na(followup) || followup <= 0) {
      stop("censoring = \"administrative\" needs a positive followup.", call. = FALSE)
    }
  }
  if (!K_true %in% 1:3) stop("K_true must be 1, 2 or 3.", call. = FALSE)
  if (K_true == 3 && p < 8) stop("K_true = 3 needs p >= 8.", call. = FALSE)
  if (features == "lognormal" && rho != 0) {
    stop("features = \"lognormal\" requires rho = 0.", call. = FALSE)
  }
  if (features == "lognormal") check_scalar(sdlog, "sdlog", lower = 0, lower_open = TRUE)
  with_local_seed(seed, {
    v <- rep(0, p); v[seq_len(min(3, p))] <- c(0.5, -0.8, 0.3)[seq_len(min(3, p))]
    v <- v / sqrt(sum(v^2))
    u <- rep(0, p); u[c(1, min(4, p))] <- c(0.8, 0.6)
    u <- u - as.numeric(crossprod(u, v)) * v; u <- u / sqrt(sum(u^2))
    base <- rep(0, p); base[seq_len(min(4, p))] <- c(0.5, 0.3, -0.4, 0.35)[seq_len(min(4, p))]
    base <- base - as.numeric(crossprod(base, v)) * v

    beta_list <- if (K_true == 1) list(base) else
      list(base - (beta_sep / 2) * v, base + (beta_sep / 2) * v)
    mu_list <- if (K_true == 1) list(rep(0, p)) else
      list(-(mu_sep / 2) * u, (mu_sep / 2) * u)
    if (K_true == 3) {
      ## Triangle geometry (added after sim-freeze-v1). Second profile and
      ## coefficient directions u2 (features 7-8) and v2 (features 5-6) are
      ## orthogonal to u, v and base, so {u, u2, v, v2} is orthonormal and
      ## base is orthogonal to span(v, v2): the three coefficient vectors
      ## have equal norms. Clusters sit at the vertices of equilateral
      ## triangles with side mu_sep (profiles) and beta_sep (mechanisms).
      ## Deterministic; draws no random numbers.
      unit <- function(x) x / sqrt(sum(x^2))
      v2 <- rep(0, p); v2[5:6] <- c(0.8, -0.6); v2 <- unit(v2)
      u2 <- rep(0, p); u2[7:8] <- c(0.6, 0.8); u2 <- unit(u2)
      th <- pi / 2 + 2 * pi * (0:2) / 3
      beta_list <- lapply(th, function(a) base + (beta_sep / sqrt(3)) * (cos(a) * v + sin(a) * v2))
      mu_list <- lapply(th, function(a) (mu_sep / sqrt(3)) * (cos(a) * u + sin(a) * u2))
    }

    Sig <- if (rho > 0) (1 - rho) * diag(p) + rho * matrix(1, p, p) else diag(p)
    Z <- sample(seq_len(K_true), n, replace = TRUE)
    noise <- MASS::mvrnorm(n, rep(0, p), Sig)
    if (features == "lognormal") {
      ## Added after sim-freeze-v2 (C5). Elementwise transform of the same
      ## draws; the Gaussian path above is unchanged.
      s2 <- sdlog^2
      noise <- (exp(sdlog * noise) - exp(s2 / 2)) / sqrt((exp(s2) - 1) * exp(s2))
    }
    X <- noise + do.call(rbind, mu_list[Z])
    colnames(X) <- paste0("f", seq_len(p))
    eta <- as.numeric(rowSums(X * t(sapply(Z, function(k) beta_list[[k]]))))

    Tt <- 200 * (-log(stats::runif(n)))^(1 / 1.5) * exp(-eta / 1.5)
    if (censoring == "random") {
      lo <- 1e-6; hi <- 10                     # calibrate censoring to event_rate
      for (i in 1:80) { mid <- (lo + hi) / 2
        if (mean(Tt <= pmin(stats::rexp(n, mid), 1e4)) < event_rate) hi <- mid else lo <- mid }
      Cc <- pmin(stats::rexp(n, mid), 1e4)
    } else {
      Cc <- rep(followup, n)                   # administrative: fixed follow-up
    }
    list(X = X, time = pmax(pmin(Tt, Cc), 0.1), status = as.integer(Tt <= Cc),
         Z = Z, eta = eta, beta_list = beta_list, p = p, mu_list = mu_list)
  })
}
