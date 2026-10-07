###############################################################
## C_methods.R  --  estimators added after sim-freeze-v2, shared by
## 05_competitor.R and 07_convergence_rerun.R. Source after 00_core.R and
## C_common.R. Moved here unchanged from 05_competitor.R (commit 9b57ed8);
## m_gemcox_tight() gained gamma, normalize and temp arguments whose
## defaults reproduce its earlier behaviour exactly.
##
## ASCII only.
###############################################################

## ---- the competitor ------------------------------------------------------
## Latent-class Cox with multinomial logistic gating (gemcox:::lcc_fit,
## internal). GeM-Cox's estimator settings: lambda, alpha, shared baseline,
## k-means initialisation with seed 1, one start, max_iter 100, tol 1e-4.
## Gating ridge lambda_gate = lambda (primary) or lambda / 10 (sensitivity).
## New subjects are assigned by the gating, P(Z = k | x); the test-set
## predictor is the uncentred mixture sum_k P(Z = k | x) x' beta_k, as for
## GeM-Cox.
m_lcc <- function(d, te, cfg, lambda_gate = cfg$lambda, tol = 1e-4, max_iter = 100) {
  n_warn <- 0L
  f <- withCallingHandlers(
    gemcox:::lcc_fit(d$X, time = d$time, status = d$status, K = cfg$K_fit,
                     lambda = cfg$lambda, alpha = cfg$alpha, lambda_gate = lambda_gate,
                     seed = 1, tol = tol, max_iter = max_iter),
    warning = function(w) { n_warn <<- n_warn + 1L; invokeRestart("muffleWarning") })
  tn <- gemcox:::lcc_predict_tau(f, te$X)
  list(beta_list = lapply(seq_len(ncol(f$beta)), function(k) f$beta[, k]),
       eta = rowSums(tn * (te$X %*% f$beta)),
       tau = f$tau, labels = f$cluster, fit = f,
       diag = data.frame(converged = f$converged, iterations = f$iterations,
                         guard_bound = NROW(f$guards) > 0, n_warnings = n_warn,
                         gate_ok = f$gate_ok, lambda_gate = lambda_gate))
}
LCC_METHODS <- stats::setNames(list(
  function(d, te, cfg) m_lcc(d, te, cfg),
  function(d, te, cfg) m_lcc(d, te, cfg, lambda_gate = cfg$lambda / 10)), c(LCC, LCC_SENS))

## ---- GeM-Cox to a tight tolerance ------------------------------------------
## GeM-Cox's frozen method (m_gemcox in 00_core.R) passes no tolerance, so
## the tight fit is written out here with the same arguments plus tol and
## max_iter.
m_gemcox_tight <- function(d, te, cfg, gamma = 1, normalize = cfg$normalize, temp = 1) {
  n_warn <- 0L
  f <- withCallingHandlers(
    gemcox(d$X, time = d$time, status = d$status, K = cfg$K_fit, gamma = gamma,
           lambda = cfg$lambda, alpha = cfg$alpha, normalize_gmm_by_dim = normalize,
           temp = temp, seed = 1, tol = TIGHT_TOL, max_iter = TIGHT_MAXIT),
    warning = function(w) { n_warn <<- n_warn + 1L; invokeRestart("muffleWarning") })
  tn <- predict(f, te$X, type = "tau")
  list(beta_list = lapply(seq_len(ncol(f$beta)), function(k) f$beta[, k]),
       eta = rowSums(tn * (te$X %*% coef(f))),
       tau = f$tau, labels = f$cluster, fit = f,
       diag = data.frame(converged = f$converged, iterations = f$iterations,
                         guard_bound = nrow(f$guards) > 0, n_warnings = n_warn))
}
TIGHT_METHODS <- stats::setNames(list(
  function(d, te, cfg) m_lcc(d, te, cfg, tol = TIGHT_TOL, max_iter = TIGHT_MAXIT),
  function(d, te, cfg) m_gemcox_tight(d, te, cfg)), c(LCC_TIGHT, GEM_TIGHT))

## ---- added for the convergence rerun (R1) -----------------------------------
GEM0_TIGHT <- "GeM-Cox (gamma=0, tol 1e-8)"
GEM0_TIGHT_METHOD <- stats::setNames(list(
  function(d, te, cfg) m_gemcox_tight(d, te, cfg, gamma = 0)), GEM0_TIGHT)

## ---- matched-regularisation comparison (sensitivity analysis) ---------------
## GeM-Cox (gamma = 1): normalize_gmm_by_dim x temp; competitor: gating
## penalty lambda x {1/10, 1, 10}; all to tolerance 1e-8. The default
## configurations are GEM_TIGHT and LCC_TIGHT (fitted in the convergence
## arm, reused); the five and two others are defined here.
MATCHED_GEM <- subset(expand.grid(normalize = c(FALSE, TRUE), temp = c(1, 2, 5)),
                      !(normalize == FALSE & temp == 1))
MATCHED_GEM$method <- sprintf("GeM-Cox (gamma=1, normalize=%s, temp=%g, tol 1e-8)",
                              MATCHED_GEM$normalize, MATCHED_GEM$temp)
MATCHED_LCC <- data.frame(mult = c(0.1, 10),
                          method = c("Latent-class Cox (gating lambda/10, tol 1e-8)",
                                     "Latent-class Cox (gating 10 lambda, tol 1e-8)"))
MATCHED_METHODS <- c(
  stats::setNames(lapply(seq_len(nrow(MATCHED_GEM)), function(i) {
    nz <- MATCHED_GEM$normalize[i]; tp <- MATCHED_GEM$temp[i]
    function(d, te, cfg) m_gemcox_tight(d, te, cfg, normalize = nz, temp = tp)
  }), MATCHED_GEM$method),
  stats::setNames(lapply(MATCHED_LCC$mult, function(mu) {
    function(d, te, cfg) m_lcc(d, te, cfg, lambda_gate = cfg$lambda * mu,
                               tol = TIGHT_TOL, max_iter = TIGHT_MAXIT)
  }), MATCHED_LCC$method))
## configuration table for every matched-comparison arm, defaults included
MATCHED_CONFIGS <- rbind(
  data.frame(method = GEM_TIGHT, family = "GeM-Cox", normalize = FALSE, temp = 1,
             gate_mult = NA_real_, default = TRUE),
  data.frame(method = MATCHED_GEM$method, family = "GeM-Cox", normalize = MATCHED_GEM$normalize,
             temp = MATCHED_GEM$temp, gate_mult = NA_real_, default = FALSE),
  data.frame(method = LCC_TIGHT, family = "Latent-class Cox", normalize = NA, temp = NA,
             gate_mult = 1, default = TRUE),
  data.frame(method = MATCHED_LCC$method, family = "Latent-class Cox", normalize = NA,
             temp = NA, gate_mult = MATCHED_LCC$mult, default = FALSE))

## ---- row helpers -------------------------------------------------------------
hash_extra <- function(d, te, m, cfg)
  data.frame(hash_train = hash_train(d), hash_test = hash_test(te), stringsAsFactors = FALSE)

label_rows <- function(x, key, design) {
  x$experiment <- key
  x$K_true     <- if (is.null(design$K_true)) 2 else design$K_true
  x$event_rate <- if (is.null(design$event_rate)) CFG$event_rate else design$event_rate
  x$features   <- if (is.null(design$features)) "gaussian" else design$features
  x
}

## ---- E6: methods and per-fit records (added after sim-freeze-v2) -------------
## Every EM fit at tolerance 1e-8, cap 1000; two-stage and oracle use coxph.
E6_METHODS <- c(TIGHT_METHODS[c(GEM_TIGHT, LCC_TIGHT)],
                METHODS[c("Two-stage (GMM+Cox)", "Oracle (true labels)")])
## Profile recovery: |cos| between the estimated and true between-subgroup
## mean-difference directions in X (true: mu_2 - mu_1). GeM-Cox: the fitted
## Gaussian means (original scale). Latent-class Cox: the gating slopes,
## which are proportional to Sigma^-1 (mu_2 - mu_1) = mu_2 - mu_1 for the
## identity covariance of the DGP. Two-stage: membership-weighted means
## (mclust's means up to its prior shrinkage). Oracle: true-label means.
## NA when mu_sep = 0 (no true direction).
profile_direction <- function(m, d) {
  K <- length(m$beta_list)
  if (K != 2) return(NA_real_)
  if (!is.null(m$fit) && inherits(m$fit, "gemcox")) {
    sc <- m$fit$gmm_scaler
    mu <- m$fit$mu * sc$scale + sc$center
    return(mu[, 2] - mu[, 1])
  }
  if (!is.null(m$fit) && inherits(m$fit, "lcc")) {
    A <- m$fit$gate_slopes_orig
    return(A[, 2] - A[, 1])
  }
  w <- if (!is.null(m$tau)) m$tau else outer(d$Z, 1:2, "==") * 1
  mk <- sapply(1:2, function(k) colSums(d$X * w[, k]) / sum(w[, k]))
  mk[, 2] - mk[, 1]
}
e6_extra <- function(d, te, m, cfg) {
  ev_sub <- tapply(d$status, factor(d$Z, levels = 1:2), sum)
  truth <- d$mu_list[[2]] - d$mu_list[[1]]
  pr <- NA_real_
  if (sqrt(sum(truth^2)) > 0) {
    est <- profile_direction(m, d)
    if (all(is.finite(est)) && sqrt(sum(est^2)) > 0) {
      pr <- abs(sum(est * truth)) / (sqrt(sum(est^2)) * sqrt(sum(truth^2)))
    }
  }
  data.frame(hash_train = hash_train(d), hash_test = hash_test(te),
             events_sub1 = unname(ev_sub[1]), events_sub2 = unname(ev_sub[2]),
             events_sub_min = min(ev_sub), events_per_coef = min(ev_sub) / ncol(d$X),
             realized_rate = mean(d$status), profile_recovery = pr)
}
