## inference.R -- the cross-validated partial likelihood (Verweij & van
## Houwelingen 1993) and gemcox_permutation_test(): K = 2 vs K = 1.
##
## Breslow jumps are NOT used as a density for held-out subjects: a jump is
## not a continuous-time density at an arbitrary held-out event time. The
## cross-validated partial likelihood avoids needing one.

# Cox log partial likelihood with Breslow ties for a given linear predictor.
cox_partial_loglik <- function(time, status, eta) {
  o <- order(time, decreasing = TRUE)
  tt <- time[o]
  d <- status[o]
  e <- eta[o]
  m <- max(e)
  cs <- cumsum(exp(e - m))
  last <- length(tt) + 1L - match(tt, rev(tt))  # risk set = all with time >= t
  sum(d * (e - (m + log(cs[last]))))
}

# Proper log-density of the fitted Gaussian mixture on the original X_gmm
# scale (log sum_k pi_k N(x | mu_k, Sigma_k) minus the scaling Jacobian).
gmm_mixture_logdens <- function(object, Xg) {
  Xs <- apply_x_scaler(Xg, object$gmm_scaler)
  L <- gmm_log_evidence(Xs, object$pi, object$mu, object$Sigma,
                        object$settings$control$gmm_ridge, gmm_scale = 1,
                        pi_min = 1e-12)
  rowLogSumExp(L) - sum(log(object$gmm_scaler$scale))
}

# Stratified fold labels (events and censored split separately).
stratified_folds <- function(status, nfolds, seed) {
  with_local_seed(seed, {
    fold <- integer(length(status))
    for (grp in c(0, 1)) {
      idx <- which(status == grp)
      fold[idx] <- sample(rep(seq_len(nfolds), length.out = length(idx)))
    }
    fold
  })
}

# Cross-validated criterion for each K, on fixed folds. The per-K partial
# likelihood and held-out GMM log-density parts are attached as the
# "components" attribute (the GMM part is computed for either criterion; it
# enters the total only for "joint").
cv_loglik_folds <- function(X_gmm, X_cox, time, status, folds, K, criterion, args) {
  out <- stats::setNames(numeric(length(K)), paste0("K", K))
  comp <- matrix(0, 2, length(K), dimnames = list(c("partial", "gmm"), names(out)))
  for (f in sort(unique(folds))) {
    tr <- folds != f
    for (j in seq_along(K)) {
      fit <- do.call(gemcox, c(list(X_gmm = X_gmm[tr, , drop = FALSE],
                                    X_cox = X_cox[tr, , drop = FALSE],
                                    time = time[tr], status = status[tr],
                                    K = K[j]), args))
      eta <- predict(fit, X_gmm, X_cox, type = "lp")
      val <- cox_partial_loglik(time, status, eta) -
        cox_partial_loglik(time[tr], status[tr], eta[tr])
      gmm <- sum(gmm_mixture_logdens(fit, X_gmm[!tr, , drop = FALSE]))
      comp["partial", j] <- comp["partial", j] + val
      comp["gmm", j] <- comp["gmm", j] + gmm
      if (criterion == "joint") val <- val + gmm
      out[j] <- out[j] + val
    }
  }
  attr(out, "components") <- comp
  out
}

# Pull X_gmm / X_cox / time / status and gemcox() settings from the inputs.
resolve_test_inputs <- function(object_or_data, dots) {
  fit_names <- setdiff(names(formals(gemcox)),
                       c("X_gmm", "X_cox", "time", "status", "K", "verbose"))
  if (inherits(object_or_data, "gemcox")) {
    d <- object_or_data$data
    st <- object_or_data$settings
    args <- st[intersect(names(st), fit_names)]
  } else {
    d <- if (is.null(object_or_data)) dots else object_or_data
    if (!is.list(d) || is.null(d$X_gmm) || is.null(d$time) || is.null(d$status)) {
      stop("Supply a gemcox fit, a list with X_gmm, time and status, or those as arguments.",
           call. = FALSE)
    }
    if (is.null(d$X_cox)) d$X_cox <- d$X_gmm
    args <- list()
  }
  extra <- dots[intersect(names(dots), fit_names)]
  args[names(extra)] <- extra
  X_gmm <- as_feature_matrix(d$X_gmm, "X_gmm")
  X_cox <- as_feature_matrix(d$X_cox, "X_cox")
  y <- validate_outcome(d$time, d$status, nrow(X_gmm))
  list(X_gmm = X_gmm, X_cox = X_cox, time = y$time, status = y$status, args = args)
}

#' Cross-validated log-likelihood of GeM-Cox fits
#'
#' Cross-validated partial likelihood of Verweij and van Houwelingen (1993):
#' for each fold k the model is fitted without fold k, and the fold
#' contributes `l(beta_-k) - l_-k(beta_-k)`, the full-data log partial
#' likelihood minus the training-data log partial likelihood, both at the
#' linear predictor of the fit without fold k. For a mixture the linear
#' predictor is `sum_k tau_k(x) x' beta_k` with `tau` from the clustering
#' features only (see [predict.gemcox()]). With `criterion = "joint"` the
#' held-out Gaussian-mixture log-density of fold k's clustering features is
#' added, which rewards structure in feature space as well.
#'
#' @param X_gmm,X_cox,time,status Data, as in [gemcox()].
#' @param K Vector of numbers of clusters to evaluate.
#' @param nfolds Number of folds, stratified by event status.
#' @param criterion `"partial"` or `"joint"`.
#' @param seed Seed for the fold assignment; `NULL` draws it from the
#'   current random number stream.
#' @param ... Further arguments to [gemcox()] (for example `gamma`, `lambda`).
#' @return Named numeric vector, one cross-validated log-likelihood per K
#'   (larger is better). Attribute `"components"` is a 2 x length(K) matrix
#'   with the partial-likelihood part and the held-out GMM log-density part
#'   (the "profile part"); the latter enters the total only for `"joint"`.
#' @references Verweij PJM, van Houwelingen HC (1993). Cross-validation in
#'   survival analysis. *Statistics in Medicine* 12:2305-2314.
#' @examples
#' d <- gemcox_simulate(n = 200, p = 3, mu_sep = 2, beta_sep = 2, seed = 4)
#' gemcox_cv_loglik(d$X, time = d$time, status = d$status, nfolds = 3, seed = 1)
#' @export
gemcox_cv_loglik <- function(X_gmm, X_cox = X_gmm, time, status, K = 1:2,
                             nfolds = 5, criterion = c("partial", "joint"),
                             seed = NULL, ...) {
  criterion <- match.arg(criterion)
  inp <- resolve_test_inputs(list(X_gmm = X_gmm, X_cox = X_cox, time = time,
                                  status = status), list(...))
  check_scalar(nfolds, "nfolds", lower = 2, upper = sum(inp$status), integer = TRUE)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)
  folds <- stratified_folds(inp$status, nfolds, seed)
  n_cap <- 0L
  out <- withCallingHandlers(
    cv_loglik_folds(inp$X_gmm, inp$X_cox, inp$time, inp$status, folds, K, criterion, inp$args),
    gemcox_not_converged = function(w) { n_cap <<- n_cap + 1L; invokeRestart("muffleWarning") })
  if (n_cap > 0) {
    warning(sprintf("%d of the %d fold fits stopped at max_iter without converging.",
                    n_cap, nfolds * length(K)), call. = FALSE)
  }
  out
}


# Parametric bootstrap generator under a fitted K = 1 model. Event times are
# drawn from S(t | x) = exp(-Lambda0(t) exp(eta)) using the fitted Breslow
# baseline (so they fall on the observed event-time grid, or beyond its end);
# censoring times are drawn from the Kaplan-Meier estimate of the censoring
# distribution. Follow-up is truncated at the largest observed time.
k1_bootstrap_sampler <- function(fit1, X_cox, time, status) {
  eta <- predict_eta(fit1, X_cox)[, 1]
  base <- fit1$baseline[[1]]
  km <- survival::survfit(Surv(time, 1 - status) ~ 1)
  cens_t <- km$time
  cens_s <- km$surv
  tmax <- max(time)
  n <- length(time)
  function() {
    ## event time: first grid time with Lambda0(t) >= E / exp(eta), E ~ Exp(1)
    H <- -log(stats::runif(n)) / exp(eta)
    j <- findInterval(H, base$cumhaz, left.open = TRUE) + 1L
    Tt <- rep(Inf, n)
    ok <- j <= length(base$time)
    Tt[ok] <- base$time[j[ok]]
    ## censoring time: first time with G(t) <= U, U ~ U(0, 1)
    U <- stats::runif(n)
    jc <- findInterval(-U, -cens_s, left.open = TRUE) + 1L
    Cc <- rep(Inf, n)
    okc <- jc <= length(cens_t)
    Cc[okc] <- cens_t[jc[okc]]
    list(time = pmin(Tt, Cc, tmax), status = as.integer(Tt <= Cc & Tt <= tmax))
  }
}

# Parametric bootstrap generator under a fitted model with K0 >= 2 clusters
# (added after sim-freeze-v1). Features are kept. Each subject's cluster is
# drawn from the fitted Gaussian mixture's membership probabilities given
# its clustering features (predict(type = "tau"), the model's P(Z | x));
# the event time from that cluster's Cox model with its fitted baseline;
# censoring as in the K = 1 sampler.
kK_bootstrap_sampler <- function(fitK, X_gmm, X_cox, time, status) {
  tau <- predict_tau(fitK, X_gmm)
  K <- ncol(tau)
  cum <- t(apply(tau, 1, cumsum))
  E <- predict_eta(fitK, X_cox)
  bases <- fitK$baseline
  km <- survival::survfit(Surv(time, 1 - status) ~ 1)
  cens_t <- km$time
  cens_s <- km$surv
  tmax <- max(time)
  n <- length(time)
  function() {
    Z <- 1L + as.integer(rowSums(stats::runif(n) > cum[, -K, drop = FALSE]))
    H <- -log(stats::runif(n)) / exp(E[cbind(seq_len(n), Z)])
    Tt <- rep(Inf, n)
    for (k in sort(unique(Z))) {
      idx <- which(Z == k)
      b <- bases[[k]]
      j <- findInterval(H[idx], b$cumhaz, left.open = TRUE) + 1L
      ok <- j <= length(b$time)
      Tt[idx[ok]] <- b$time[j[ok]]
    }
    U <- stats::runif(n)
    jc <- findInterval(-U, -cens_s, left.open = TRUE) + 1L
    Cc <- rep(Inf, n)
    okc <- jc <= length(cens_t)
    Cc[okc] <- cens_t[jc[okc]]
    list(time = pmin(Tt, Cc, tmax), status = as.integer(Tt <= Cc & Tt <= tmax))
  }
}

#' Test for heterogeneous survival mechanisms (K = 2 vs K = 1)
#'
#' Tests whether two latent subgroups with different Cox coefficient vectors
#' fit better than one. The null hypothesis is **one mechanism with a shared
#' beta**.
#'
#' **All calibration results below were obtained with the earlier EM
#' setting (`tol = 1e-4, max_iter = 100`), not with the current default
#' (`tol = 1e-8, max_iter = 1000`).**
#' - Calibration at the current default has not yet been established. A
#'   rerun was started but not completed.
#' - At the earlier setting both the observed and the null fits stopped
#'   early, so the numbers below are provisional for the current default.
#' - To use the setting they describe, set
#'   `options(gemcox.tol = 1e-4, gemcox.max_iter = 100)`.
#' - Each test with `B` replicates runs `2 B + 3` fits. At the current
#'   default that is roughly 1-10 s per fit at n = 400-800.
#'
#' **Warning: not validated when the features form distinct profile
#' clusters.** In experiment E5d (added after the frozen study), subgroups
#' differed only in their feature profiles and shared one Cox coefficient
#' vector.
#' * At mu_sep = 0.5 the default test rejected at 0.070 (MCSE 0.011; 35 of
#'   500 datasets), above the pre-registered tolerance of 0.0695.
#' * At mu_sep = 1, at mu_sep = 2, and at mu_sep = 1 with correlated
#'   features it rejected at 0.052, 0.046 and 0.060.
#'
#' By the rule recorded in advance, the test is not validated for such
#' data: a rejection may partly reflect profile structure rather than
#' different survival mechanisms. The excess is small, and chance cannot be
#' ruled out (probability about 0.09 of one such cell among three), but a
#' real inflation cannot be either.
#'
#' **Statistic.** The default, `statistic = "lrt"`, is the in-sample
#' likelihood ratio `T = 2 (l_2 - l_1)`. Here `l_K` is the maximised mixture
#' log-likelihood of the full-data fit with K clusters: the Gaussian mixture
#' part plus the survival part with the fitted Breslow baseline, at the
#' penalised estimates.
#' * **Requires `gamma = 1` and `temp = 1`.** Only then is the fitted
#'   objective the mixture log-likelihood; other values are refused.
#' * With `normalize_gmm_by_dim = TRUE` the objective down-weights the
#'   Gaussian part, so `T` is a likelihood ratio of that tempered objective.
#'   It is still a valid test statistic because it is calibrated by
#'   simulation, not by a reference distribution.
#' * The chi-square reference is never used: it is invalid for mixtures.
#'
#' `statistic = "cv"` uses instead the gain in cross-validated partial
#' likelihood, `T = CV(K = 2) - CV(K = 1)` (see [gemcox_cv_loglik()]).
#'
#' **Null distribution.** The default, `null = "bootstrap"`, is a parametric
#' bootstrap under the fitted K = 1 model. The features are kept; event
#' times are simulated from the fitted Cox model with its Breslow baseline,
#' and censoring times from the Kaplan-Meier estimate of the censoring
#' distribution, with follow-up truncated at the largest observed time.
#' This keeps the common Cox effect, so it targets "one mechanism" directly,
#' at the price of relying on the fitted model.
#'
#' `null = "permutation"` permutes `(time, status)` jointly across subjects.
#' That also destroys the common effect, so strictly it tests "no
#' association between features and outcome".
#'
#' The p-value is `(1 + #{T_null >= T_obs}) / (1 + n_valid)`, where
#' `n_valid` counts null replicates whose fits succeeded.
#'
#' **Calibration in simulation (earlier EM setting, see above).** Source:
#' `gemcox/inst/sim`, experiments E5a and E5c. Setting: n = 400 or 800,
#' p = 10, B = 19, K = 1 null data
#' with a nonzero shared coefficient vector. Rejection rates at 0.05:
#'
#' | test | null data, independent features (500 datasets) | exchangeable correlation 0.5 (200) | power at beta_sep = 2 (4 cells, 200 each) |
#' |---|---|---|---|
#' | LRT, bootstrap (default) | 0.052 (MCSE 0.010) | 0.030 (conservative) | 0.30-0.57 |
#' | LRT, permutation | 0.028 (conservative) | 0.030 | 0.21-0.57 |
#' | CV, bootstrap | 0.062 | 0.060 | 0.05-0.10 |
#' | CV, permutation | 0.112 (anti-conservative) | 0.055 | 0.05-0.10 |
#'
#' Further results for the default test:
#' * At beta_sep = 1 it had power 0.09-0.16.
#' * At CVIA078 scale (n = 117, about 52 events, p = 9; experiment E2b) it
#'   was calibrated (0.055 under the null) but had only 0.13-0.14 power at
#'   beta_sep = 2.
#'
#' With `normalize_gmm_by_dim = TRUE`:
#' * at n = 400-800 the same test was more powerful (0.71-0.82 at
#'   beta_sep = 2, null rejection 0.040);
#' * it was anti-conservative with few events: 0.135 under the null at
#'   n = 117 (about 52 events) and 0.080 at n = 200 (about 88 events);
#' * it reached the nominal level at n = 300 (0.050, about 132 events;
#'   experiment E5f).
#'
#' Power is modest when events are few, so a non-significant result is
#' inconclusive: it is not evidence that the mechanisms are homogeneous.
#'
#' @param object_or_data A `"gemcox"` fit (its data and settings are
#'   reused), a list with `X_gmm`, `time`, `status` (and optionally `X_cox`),
#'   or `NULL` with those supplied through `...`.
#' @param statistic `"lrt"` (default, in-sample likelihood ratio) or `"cv"`
#'   (cross-validated partial likelihood); see Details.
#' @param null `"bootstrap"` (default) or `"permutation"`; see Details.
#' @param B Number of null replicates (bootstrap samples or permutations).
#' @param nfolds,criterion Folds and criterion (`"partial"` or `"joint"`,
#'   see [gemcox_cv_loglik()]); used only by `statistic = "cv"`.
#' @param K0 Number of clusters under the null; the test compares `K0 + 1`
#'   against `K0` (default 1: two subgroups against one). `K0 >= 2` (added
#'   after the frozen study, for sequential selection of K) requires
#'   `statistic = "lrt"` and `null = "bootstrap"`. The bootstrap then draws
#'   each subject's cluster from the fitted K0 model's membership
#'   probabilities given the features, and an event time from that
#'   cluster's Cox model. Its calibration has not been established.
#'   Used sequentially (2 vs 1, then 3 vs 2, ...; experiment E5e), this was
#'   the only selector tested that detected differing mechanisms without
#'   clear profile separation. It selected K = 2 in 20-38% of
#'   mechanism-only datasets with two subgroups, and K = 1 in 93-98% of
#'   null datasets.
#' @param observed Optional result of an earlier call on the same data with
#'   the same `statistic`, `nfolds`, `criterion` and fit settings. Its
#'   observed statistic (and folds) are reused, so that both nulls can be
#'   compared against one observed statistic without refitting it.
#' @param seed Seed for the null replicates and folds; `NULL` uses the
#'   current random number stream.
#' @param verbose Print progress.
#' @param ... Data (`X_gmm`, `X_cox`, `time`, `status`) when
#'   `object_or_data` is `NULL`, and/or [gemcox()] settings overriding those
#'   of a supplied fit.
#' @return An object of class `"htest"` with the observed statistic,
#'   p-value, the null distribution (`null.distribution`), the numbers of
#'   null replicates run and valid, the per-K criterion values (`cv_loglik`:
#'   full-data mixture log-likelihoods for `"lrt"`, cross-validated values
#'   for `"cv"`), and the fold seed (`fold_seed`). For `"lrt"` it also
#'   counts the fits that stopped at `max_iter` without meeting the
#'   convergence tolerance (`n_not_converged`, out of `n_fits`: the null
#'   model fitted for the bootstrap and every K0 and K0 + 1 fit on the
#'   observed and null data). Added after the frozen study; `NA` for
#'   `"cv"`.
#' @references Verweij PJM, van Houwelingen HC (1993). Cross-validation in
#'   survival analysis. *Statistics in Medicine* 12:2305-2314.
#' @examples
#' \donttest{
#' d <- gemcox_simulate(n = 200, p = 3, mu_sep = 0, beta_sep = 2, seed = 5)
#' ht <- gemcox_heterogeneity_test(X_gmm = d$X, time = d$time, status = d$status,
#'                                 B = 9, seed = 1)
#' ht
#' ## the same observed statistic against the permutation null
#' gemcox_heterogeneity_test(X_gmm = d$X, time = d$time, status = d$status,
#'                           B = 9, seed = 1, null = "permutation", observed = ht)
#' }
#' @export
gemcox_heterogeneity_test <- function(object_or_data = NULL,
                                      statistic = c("lrt", "cv"),
                                      null = c("bootstrap", "permutation"),
                                      B = 200, nfolds = 5,
                                      criterion = c("partial", "joint"), K0 = 1,
                                      observed = NULL, seed = NULL, verbose = FALSE, ...) {
  heterogeneity_test_impl(object_or_data, n_perm = B, nfolds = nfolds,
                          criterion = match.arg(criterion), null_type = match.arg(null),
                          stat_type = match.arg(statistic), observed = observed,
                          seed = seed, verbose = verbose, dots = list(...), K0 = K0)
}

#' Deprecated: use gemcox_heterogeneity_test()
#'
#' `gemcox_permutation_test()` is the former name of
#' [gemcox_heterogeneity_test()]. It keeps its original signature and
#' defaults (`statistic = "cv"`, `n_perm`), so existing code returns the same
#' results, and it signals a deprecation warning. Note that the new function
#' defaults to the likelihood-ratio statistic.
#'
#' @inheritParams gemcox_heterogeneity_test
#' @param n_perm Number of null replicates.
#' @param statistic `"cv"` (the former default) or `"lrt"`.
#' @return As [gemcox_heterogeneity_test()].
#' @keywords internal
#' @export
gemcox_permutation_test <- function(object_or_data = NULL, n_perm = 200, nfolds = 5,
                                    criterion = c("partial", "joint"),
                                    null = c("bootstrap", "permutation"),
                                    statistic = c("cv", "lrt"), observed = NULL, seed = NULL,
                                    verbose = FALSE, ...) {
  .Deprecated("gemcox_heterogeneity_test", package = "gemcox",
              msg = paste("gemcox_permutation_test() is deprecated; use",
                          "gemcox_heterogeneity_test(), whose default statistic is \"lrt\"."))
  heterogeneity_test_impl(object_or_data, n_perm = n_perm, nfolds = nfolds,
                          criterion = match.arg(criterion), null_type = match.arg(null),
                          stat_type = match.arg(statistic), observed = observed,
                          seed = seed, verbose = verbose, dots = list(...))
}

heterogeneity_test_impl <- function(object_or_data, n_perm, nfolds, criterion, null_type,
                                    stat_type, observed, seed, verbose, dots, K0 = 1) {
  inp <- resolve_test_inputs(object_or_data, dots)
  check_scalar(n_perm, "the number of null replicates (B)", lower = 1, integer = TRUE)
  check_scalar(K0, "K0", lower = 1, integer = TRUE)
  if (K0 > 1 && (stat_type != "lrt" || null_type != "bootstrap")) {
    stop("K0 >= 2 needs statistic = \"lrt\" and null = \"bootstrap\".", call. = FALSE)
  }
  Ks <- c(K0, K0 + 1)
  if (stat_type == "cv") {
    check_scalar(nfolds, "nfolds", lower = 2, upper = sum(inp$status), integer = TRUE)
  } else {
    g <- if (is.null(inp$args$gamma)) 1 else inp$args$gamma
    tp <- if (is.null(inp$args$temp)) 1 else inp$args$temp
    if (g != 1 || tp != 1) {
      stop("statistic = \"lrt\" needs gamma = 1 and temp = 1: only then is the fitted ",
           "objective the mixture log-likelihood.", call. = FALSE)
    }
  }
  n <- length(inp$time)

  n_guard_fits <- 0L
  n_fits <- 0L; n_not_converged <- 0L
  count_fit <- function(fit) {
    n_fits <<- n_fits + 1L
    if (!isTRUE(fit$converged)) n_not_converged <<- n_not_converged + 1L
    fit
  }
  n_cap_warn <- 0L
  quiet_guards <- function(expr) {
    withCallingHandlers(expr, gemcox_guard_summary = function(w) {
      n_guard_fits <<- n_guard_fits + 1L
      invokeRestart("muffleWarning")
    }, gemcox_not_converged = function(w) {
      n_cap_warn <<- n_cap_warn + 1L
      invokeRestart("muffleWarning")
    })
  }

  sampler <- NULL
  if (null_type == "bootstrap") {
    fit0 <- count_fit(quiet_guards(do.call(gemcox, c(list(X_gmm = inp$X_gmm, X_cox = inp$X_cox,
                                                          time = inp$time, status = inp$status,
                                                          K = K0), inp$args))))
    sampler <- if (K0 == 1) {
      k1_bootstrap_sampler(fit0, inp$X_cox, inp$time, inp$status)
    } else {
      kK_bootstrap_sampler(fit0, inp$X_gmm, inp$X_cox, inp$time, inp$status)
    }
  }
  draw <- function() {
    list(fold_seed = sample.int(.Machine$integer.max, 1),
         null_data = lapply(seq_len(n_perm), function(b) {
           if (null_type == "permutation") {
             idx <- sample.int(n)
             list(time = inp$time[idx], status = inp$status[idx])
           } else {
             sampler()
           }
         }))
  }
  rnd <- if (is.null(seed)) draw() else with_local_seed(seed, draw())
  if (!is.null(observed)) {
    ok <- inherits(observed, "htest") && !is.null(observed$cv_loglik) &&
      !is.null(observed$fold_seed) && isTRUE(observed$nfolds == nfolds) &&
      identical(observed$criterion, criterion) &&
      identical(observed$statistic_type, stat_type) &&
      isTRUE((if (is.null(observed$null_K)) 1 else observed$null_K) == K0) &&
      identical(observed$data.name, sprintf("%d subjects, %d events", n, sum(inp$status)))
    if (!ok) {
      stop("`observed` must come from gemcox_permutation_test() on the same data, ",
           "with the same statistic, nfolds and criterion.", call. = FALSE)
    }
    rnd$fold_seed <- observed$fold_seed
  }

  stat_for <- function(time, status) {
    if (stat_type == "lrt") {
      ll <- vapply(Ks, function(K) {
        count_fit(quiet_guards(do.call(gemcox, c(list(X_gmm = inp$X_gmm, X_cox = inp$X_cox,
                                                      time = time, status = status, K = K),
                                                 inp$args))))$final_score
      }, 0)
      return(stats::setNames(ll, paste0("K", Ks)))
    }
    folds <- stratified_folds(status, nfolds, rnd$fold_seed)
    quiet_guards(cv_loglik_folds(inp$X_gmm, inp$X_cox, time, status, folds, 1:2,
                                 criterion, inp$args))
  }
  stat_of <- function(v) {
    d <- unname(v[2] - v[1])          # K0 + 1 minus K0
    if (stat_type == "lrt") 2 * d else d
  }
  min_events <- if (stat_type == "cv") nfolds else 2

  cv_obs <- if (is.null(observed)) stat_for(inp$time, inp$status) else observed$cv_loglik
  t_obs <- stat_of(cv_obs)
  null_stat <- rep(NA_real_, n_perm)
  for (b in seq_len(n_perm)) {
    yb <- rnd$null_data[[b]]
    cv_b <- if (sum(yb$status) >= min_events) {
      tryCatch(stat_for(yb$time, yb$status), error = function(e) NULL)
    }
    if (!is.null(cv_b)) null_stat[b] <- stat_of(cv_b)
    if (verbose && b %% 10 == 0) cat(sprintf("  null replicate %d / %d\n", b, n_perm))
  }
  n_valid <- sum(is.finite(null_stat))
  if (n_valid < n_perm) {
    warning(sprintf("%d of %d null replicates failed and were excluded.",
                    n_perm - n_valid, n_perm), call. = FALSE)
  }
  if (n_guard_fits > 0) {
    warning(sprintf("Numerical guards bound in %d of the fits (observed and null data).",
                    n_guard_fits), call. = FALSE)
  }
  if (n_cap_warn > 0) {
    warning(sprintf("%d of the fits (observed and null data) stopped at max_iter without converging.",
                    n_cap_warn), call. = FALSE)
  }
  stat_name <- if (stat_type == "lrt") sprintf("LRT 2(l%d - l%d)", K0 + 1, K0) else
    "CV log-lik gain (K=2 vs K=1)"
  structure(list(
    statistic = stats::setNames(t_obs, stat_name),
    parameter = c(`valid null replicates` = n_valid),
    p.value = (1 + sum(null_stat[is.finite(null_stat)] >= t_obs)) / (1 + n_valid),
    method = if (stat_type == "lrt") {
      sprintf("GeM-Cox test of K = %d vs K = %d: %s null, in-sample likelihood ratio",
              K0 + 1, K0,
              if (null_type == "permutation") "permutation" else "parametric bootstrap")
    } else {
      sprintf("GeM-Cox test of K = 2 vs K = 1: %s null, %s criterion, %d-fold CV",
              if (null_type == "permutation") "permutation" else "parametric bootstrap",
              if (criterion == "partial") "cross-validated partial likelihood" else "joint",
              nfolds)
    },
    data.name = sprintf("%d subjects, %d events", n, sum(inp$status)),
    null.distribution = null_stat, null = null_type, n_perm = n_perm, n_valid = n_valid,
    nfolds = nfolds, criterion = criterion, statistic_type = stat_type, null_K = K0,
    cv_loglik = cv_obs, fold_seed = rnd$fold_seed,
    n_fits = if (stat_type == "lrt") n_fits else NA_integer_,
    n_not_converged = if (stat_type == "lrt") n_not_converged else NA_integer_),
    class = "htest")
}
