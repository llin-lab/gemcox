###############################################################
## 02_experiments.R  --  the manuscript experiments
##
##  E1  Existence and scaling        Fig 1   (propositions P1, P2)
##  E1b Normalisation on E1's datasets (pre-specified after S1; separate)
##  E2  Comparison across scenarios  Fig 2
##  E3  Operating envelope           Fig 3   (n and p)
##  E4  Theory checks                Fig 4   (two scales, mechanism, softness)
##  E5  Inference                    Fig 5   (a: test calibration and power,
##                                            b: K selection by regime)
##  E5c In-sample LRT on E5a's datasets (added after the full run)
##  E2b Adaptive normalisation: (a) E2 cells, (b) LRT power, (c) CVIA078
##      scale. The final experiment; pre-specified with its decision rule.
##  S1  Normalisation                Suppl   (regime dependence)
##
## Run 01_pilot.R first. Run from gemcox/inst/sim. Select experiments with
##   GEMCOX_RUN="E1,E2" Rscript 02_experiments.R
## and rehearse quickly with GEMCOX_REPS=20 (E5 counts scale with it).
###############################################################

source("00_core.R")
REPS <- as.integer(Sys.getenv("GEMCOX_REPS", CFG$reps_full))
RUN  <- strsplit(Sys.getenv("GEMCOX_RUN", "E1,E1b,E2,E2b,E3,E4,E5,E5c,E5d,E5e,E5f,S1"), ",")[[1]]
E5_SCALE <- REPS / CFG$reps_full      # 1 in the full run

## ================= E1: existence and scaling ==========================
if ("E1" %in% RUN) {
  cells <- expand.grid(n = 400, p = 10, mu_sep = 0,
                       beta_sep = c(0.5, 1, 2, 3), KEEP.OUT.ATTRS = FALSE)
  run_experiment(1, "existence and scaling", cells, REPS, "E1_existence.rds")
}

## ================= E1b: normalisation on E1's datasets ================
## Pre-specified after S1 (S1: normalize = TRUE was better at mu_sep = 0 for
## both p = 10 and p = 80). Analyses E1's exact datasets (seed_exp = 1) and
## is reported separately from E1. Arms:
##   normalize = TRUE, and a regime-adaptive arm that fits with TRUE when the
##   PROFILE PART of the joint cross-validated criterion (held-out GMM
##   log-density, 5 folds, fits at the default settings) selects K = 1, and
##   with FALSE otherwise. The overall joint choice is recorded too.
## E1's methods are refitted on the same data as paired comparators.
m_adaptive <- function(d, te, cfg) {
  cvj <- suppressWarnings(gemcox_cv_loglik(
    d$X, time = d$time, status = d$status, K = 1:2, nfolds = 5, criterion = "joint",
    seed = d$seed, gamma = cfg$gamma, lambda = cfg$lambda, alpha = cfg$alpha,
    normalize_gmm_by_dim = cfg$normalize))
  comp <- attr(cvj, "components")
  profile_K <- unname(which.max(comp["gmm", ]))
  nz <- profile_K == 1
  m <- m_gemcox(d, te, cfg, gamma = 1, normalize = nz)
  m$diag <- cbind(m$diag, adaptive_normalize = nz, profile_K = profile_K,
                  joint_K = unname(which.max(cvj)),
                  profile_gain = unname(diff(comp["gmm", ])))
  m
}
E1B_METHODS <- c(METHODS, list(
  "GeM-Cox (gamma=1, normalize=TRUE)" = function(d, te, cfg) m_gemcox(d, te, cfg, 1, TRUE),
  "GeM-Cox (gamma=1, adaptive)"       = function(d, te, cfg) m_adaptive(d, te, cfg)))
if ("E1b" %in% RUN) {
  cells <- expand.grid(n = 400, p = 10, mu_sep = 0,
                       beta_sep = c(0.5, 1, 2, 3), KEEP.OUT.ATTRS = FALSE)
  run_experiment(11, "E1b: normalisation on E1's datasets", cells, REPS,
                 "E1b_normalisation.rds", methods = E1B_METHODS, seed_exp = 1)
}

## ================= E2: comparison across scenarios ====================
if ("E2" %in% RUN) {
  cells <- expand.grid(n = 800, p = 10,
                       mu_sep = c(0, 0.5, 1, 2), beta_sep = c(1, 2, 3),
                       KEEP.OUT.ATTRS = FALSE)
  run_experiment(2, "method comparison across separation regimes",
                 cells, REPS, "E2_comparison.rds")
}

## ================= E3: operating envelope =============================
if ("E3" %in% RUN) {
  cells_n <- expand.grid(n = c(200, 400, 800, 1600), p = 10,
                         mu_sep = 0.5, beta_sep = 2, KEEP.OUT.ATTRS = FALSE)
  run_experiment(3, "envelope: sample size", cells_n, REPS, "E3a_n.rds")
  cells_p <- expand.grid(n = 800, p = c(10, 20, 40, 80),
                         mu_sep = 0.5, beta_sep = 2, KEEP.OUT.ATTRS = FALSE)
  run_experiment(31, "envelope: dimension", cells_p, REPS, "E3b_p.rds")
}

## ================= E4: theory checks ==================================
if ("E4" %in% RUN) {
  ## (a) TWO SCALES: individual labels vs the contrast, as n grows.
  cells <- expand.grid(n = c(800, 1600, 3200, 6400), p = 10,
                       mu_sep = 0, beta_sep = 2, KEEP.OUT.ATTRS = FALSE)
  run_experiment(4, "two scales: labels vs contrast", cells, REPS,
                 "E4a_two_scales.rds")

  ## (b) MECHANISM: does the fitted log-odds track
  ##     (x' Delta) x (martingale residual)? Everything is centred at the
  ##     training means, which is where the fitted shared baseline lives.
  mech <- function(d, te, m, cfg) {
    if (is.null(m$fit) || is.null(m$tau)) return(data.frame(r_mech = NA_real_))
    f <- m$fit
    xc <- sweep(d$X, 2, f$cox_scaler$center)
    Dl <- f$beta[, 1] - f$beta[, 2]
    bb <- rowMeans(f$beta)
    M  <- d$status - cumhaz_at(f$baseline[[1]], d$time) * exp(as.numeric(xc %*% bb))
    pred <- as.numeric(xc %*% Dl) * M
    tt <- pmin(pmax(m$tau[, 1], 1e-6), 1 - 1e-6)
    data.frame(r_mech = suppressWarnings(
      cor(log(tt / (1 - tt)), pred, method = "spearman")))
  }
  cells <- expand.grid(n = 800, p = 10, mu_sep = 0, beta_sep = 2,
                       KEEP.OUT.ATTRS = FALSE)
  run_experiment(41, "mechanism: log-odds vs disagreement x residual",
                 cells, REPS, "E4b_mechanism.rds",
                 methods = METHODS[c("GeM-Cox (gamma=0)", "GeM-Cox (gamma=1)")],
                 extra = mech)

  ## (c) POSTERIOR SOFTNESS. Membership weights for the test set come from
  ## predict(), which does not apply temp.
  soft <- list(
    "base (gamma=1)"        = function(d, te, cfg) m_gemcox(d, te, cfg, 1,  FALSE, 1),
    "gamma=p"               = function(d, te, cfg) m_gemcox(d, te, cfg, 10, FALSE, 1),
    "normalize"             = function(d, te, cfg) m_gemcox(d, te, cfg, 1,  TRUE,  1),
    "gamma=p, temp=p"       = function(d, te, cfg) m_gemcox(d, te, cfg, 10, FALSE, 10),
    "normalize, sharpened"  = function(d, te, cfg) m_gemcox(d, te, cfg, 1,  TRUE,  0.05),
    "Oracle (true labels)"  = METHODS[["Oracle (true labels)"]])
  run_experiment(42, "posterior softness", cells, REPS, "E4c_softness.rds",
                 methods = soft)
}

## ================= E5: inference ======================================
## Per-dataset functions are defined at top level so they can be timed
## and rerun in isolation; the if-block below only runs them.

## (a) CALIBRATION AND POWER of the K = 2 vs K = 1 test, both nulls.
## B = 19 null replicates, 3 folds; the observed statistic is computed
## once per dataset and shared by both nulls. K = 1 null cells have a
## NONZERO shared coefficient vector (the base coefficient of the DGP).
B <- 19L; NFOLDS <- 3L
e5a_cells <- rbind(
  data.frame(role = "type I", n = 400, p = 10, mu_sep = 0, beta_sep = 0,
             rho = 0, K_true = 1, reps = 500),
  data.frame(role = "type I", n = 400, p = 10, mu_sep = 0, beta_sep = 0,
             rho = 0.5, K_true = 1, reps = 200),
  cbind(role = "power",
        expand.grid(n = c(400, 800), p = 10, mu_sep = c(0, 0.5), beta_sep = c(1, 2),
                    rho = 0, K_true = 2, reps = 200, KEEP.OUT.ATTRS = FALSE)))
e5a_cells$reps <- pmax(2L, round(e5a_cells$reps * E5_SCALE))

e5a_one <- function(r, cell_id, design, cfg) {
  pin_rng()
  s <- mk_seed(5, cell_id, r)
  d <- simulate(design, s)
  n_warn <- 0L
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(withCallingHandlers({
    args <- list(X_gmm = d$X, time = d$time, status = d$status, B = B,
                 statistic = "cv", nfolds = NFOLDS, seed = s, gamma = cfg$gamma,
                 lambda = cfg$lambda, alpha = cfg$alpha,
                 normalize_gmm_by_dim = cfg$normalize)
    bs <- do.call(gemcox_heterogeneity_test, c(args, list(null = "bootstrap")))
    pm <- do.call(gemcox_heterogeneity_test,
                  c(args, list(null = "permutation", observed = bs)))
    list(bs = bs, pm = pm)
  }, warning = function(w) { n_warn <<- n_warn + 1L; invokeRestart("muffleWarning") }),
  error = function(e) e)
  base <- data.frame(exp = 5, cell = cell_id, rep = r, seed = s, role = design$role,
                     n = design$n, p = design$p, mu_sep = design$mu_sep,
                     beta_sep = design$beta_sep, rho = design$rho,
                     K_true = design$K_true, events = sum(d$status),
                     secs = proc.time()[["elapsed"]] - t0, n_warnings = n_warn,
                     stringsAsFactors = FALSE)
  if (inherits(res, "error")) {
    return(cbind(base, failed = TRUE, error = conditionMessage(res), stat = NA_real_,
                 p_bootstrap = NA_real_, p_permutation = NA_real_,
                 valid_bootstrap = NA_integer_, valid_permutation = NA_integer_))
  }
  cbind(base, failed = FALSE, error = NA_character_,
        stat = unname(res$bs$statistic), p_bootstrap = res$bs$p.value,
        p_permutation = res$pm$p.value, valid_bootstrap = res$bs$n_valid,
        valid_permutation = res$pm$n_valid)
}

## (b) K SELECTION BY REGIME: choose K in {1, 2} by the cross-validated
## joint criterion (partial likelihood + held-out GMM log-density) and by
## the partial likelihood alone.
ksel_cells <- data.frame(
  regime = c("null", "profile only", "mechanism only", "both"),
  n = 400, p = 10, mu_sep = c(0, 3, 0, 3), beta_sep = c(0, 0, 3, 3),
  stringsAsFactors = FALSE)
e5b_one <- function(r, cell_id, design, cfg) {
  pin_rng()
  s <- mk_seed(6, cell_id, r)
  d <- simulate(design, s)
  args <- list(X_gmm = d$X, time = d$time, status = d$status, K = 1:2, nfolds = 5,
               seed = s, gamma = cfg$gamma, lambda = cfg$lambda, alpha = cfg$alpha,
               normalize_gmm_by_dim = cfg$normalize)
  res <- tryCatch(suppressWarnings(list(
    joint = do.call(gemcox_cv_loglik, c(args, list(criterion = "joint"))),
    partial = do.call(gemcox_cv_loglik, c(args, list(criterion = "partial"))))),
    error = function(e) e)
  base <- data.frame(exp = 6, cell = cell_id, rep = r, seed = s, regime = design$regime,
                     n = design$n, p = design$p, mu_sep = design$mu_sep,
                     beta_sep = design$beta_sep, events = sum(d$status),
                     stringsAsFactors = FALSE)
  if (inherits(res, "error")) {
    return(cbind(base, failed = TRUE, error = conditionMessage(res),
                 K_joint = NA_integer_, K_partial = NA_integer_,
                 gain_joint = NA_real_, gain_partial = NA_real_))
  }
  cbind(base, failed = FALSE, error = NA_character_,
        K_joint = unname(which.max(res$joint)), K_partial = unname(which.max(res$partial)),
        gain_joint = unname(diff(res$joint)), gain_partial = unname(diff(res$partial)))
}
## (c) IN-SAMPLE LRT (added after the full run): 2 (l2 - l1) on the mixture
## log-likelihood, calibrated by the same parametric bootstrap (features
## fixed, times from the fitted K = 1 Cox model, Kaplan-Meier censoring) and
## by the permutation null, sharing one observed statistic. B = 19. Same
## cells and the SAME DATASETS as E5a (seeds from experiment 5).
e5c_one <- function(r, cell_id, design, cfg) {
  pin_rng()
  s <- mk_seed(5, cell_id, r)
  d <- simulate(design, s)
  n_warn <- 0L
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(withCallingHandlers({
    args <- list(X_gmm = d$X, time = d$time, status = d$status, B = B,
                 statistic = "lrt", seed = s, gamma = cfg$gamma, lambda = cfg$lambda,
                 alpha = cfg$alpha, normalize_gmm_by_dim = cfg$normalize)
    bs <- do.call(gemcox_heterogeneity_test, c(args, list(null = "bootstrap")))
    pm <- do.call(gemcox_heterogeneity_test,
                  c(args, list(null = "permutation", observed = bs)))
    list(bs = bs, pm = pm)
  }, warning = function(w) { n_warn <<- n_warn + 1L; invokeRestart("muffleWarning") }),
  error = function(e) e)
  base <- data.frame(exp = 7, cell = cell_id, rep = r, seed = s, role = design$role,
                     n = design$n, p = design$p, mu_sep = design$mu_sep,
                     beta_sep = design$beta_sep, rho = design$rho,
                     K_true = design$K_true, events = sum(d$status),
                     secs = proc.time()[["elapsed"]] - t0, n_warnings = n_warn,
                     stringsAsFactors = FALSE)
  if (inherits(res, "error")) {
    return(cbind(base, failed = TRUE, error = conditionMessage(res), stat = NA_real_,
                 p_bootstrap = NA_real_, p_permutation = NA_real_,
                 valid_bootstrap = NA_integer_, valid_permutation = NA_integer_))
  }
  cbind(base, failed = FALSE, error = NA_character_,
        stat = unname(res$bs$statistic), p_bootstrap = res$bs$p.value,
        p_permutation = res$pm$p.value, valid_bootstrap = res$bs$n_valid,
        valid_permutation = res$pm$n_valid)
}
if ("E5c" %in% RUN) {
  cat(sprintf("\n[E5c] in-sample LRT: calibration and power (%d cells, B = %d)\n",
              nrow(e5a_cells), B))
  t0 <- Sys.time()
  E5c <- bind_rows(lapply(seq_len(nrow(e5a_cells)), function(j) {
    design <- as.list(e5a_cells[j, ])
    t1 <- Sys.time()
    x <- bind_rows(par_lapply(seq_len(design$reps), e5c_one, cell_id = j,
                              design = design, cfg = CFG))
    cat(sprintf("  cell %2d %-6s n=%d rho=%g mu=%g beta=%g: %d datasets, %d failed, reject(boot/perm) %.3f/%.3f (%.1f min)\n",
                j, design$role, design$n, design$rho, design$mu_sep, design$beta_sep,
                nrow(x), sum(x$failed), mean(x$p_bootstrap <= 0.05, na.rm = TRUE),
                mean(x$p_permutation <= 0.05, na.rm = TRUE),
                as.numeric(difftime(Sys.time(), t1, units = "mins"))))
    x
  }))
  save_results(E5c, "E5c_lrt.rds", t0)
}

if ("E5" %in% RUN) {
  cat(sprintf("\n[E5a] test calibration and power (%d cells, B = %d, %d folds)\n",
              nrow(e5a_cells), B, NFOLDS))
  t0 <- Sys.time()
  E5a <- bind_rows(lapply(seq_len(nrow(e5a_cells)), function(j) {
    design <- as.list(e5a_cells[j, ])
    t1 <- Sys.time()
    x <- bind_rows(par_lapply(seq_len(design$reps), e5a_one, cell_id = j,
                              design = design, cfg = CFG))
    cat(sprintf("  cell %2d %-6s n=%d rho=%g mu=%g beta=%g: %d datasets, %d failed, reject(boot/perm) %.3f/%.3f (%.1f min)\n",
                j, design$role, design$n, design$rho, design$mu_sep, design$beta_sep,
                nrow(x), sum(x$failed), mean(x$p_bootstrap <= 0.05, na.rm = TRUE),
                mean(x$p_permutation <= 0.05, na.rm = TRUE),
                as.numeric(difftime(Sys.time(), t1, units = "mins"))))
    x
  }))
  save_results(E5a, "E5a_calibration.rds", t0)

  cat("\n[E5b] K selection by regime\n")
  t0 <- Sys.time()
  E5b <- bind_rows(lapply(seq_len(nrow(ksel_cells)), function(j) {
    x <- bind_rows(par_lapply(seq_len(REPS), e5b_one, cell_id = j,
                              design = as.list(ksel_cells[j, ]), cfg = CFG))
    cat(sprintf("  %-15s P(K=2) joint %.2f, partial %.2f\n", ksel_cells$regime[j],
                mean(x$K_joint == 2, na.rm = TRUE), mean(x$K_partial == 2, na.rm = TRUE)))
    x
  }))
  save_results(E5b, "E5b_kselect.rds", t0)
}

## ================= E2b: adaptive normalisation (final experiment) =====
## Pre-specified before running, with its decision rule (see README):
## (a) adaptive rule vs normalize = FALSE on E2's 12 cells and datasets
##     (seed_exp = 2); normalize = TRUE is an auxiliary arm; gamma = 0 and
##     the oracle give the normalised gain. DECISION RULE for the package
##     default: adaptive replaces FALSE iff (i) in no cell is adaptive worse
##     than FALSE by more than 2 MCSE of the paired difference and (ii) the
##     paired difference averaged over the 12 cells exceeds 2 MCSE of that
##     average. Degenerate fits excluded (primary analysis).
## (b) LRT (bootstrap null, B = 19) with normalize chosen by the adaptive
##     rule, on E5a/E5c's cells and datasets (seeds from experiment 5),
##     including the two type I cells, so it pairs with E5c's default LRT.
##     The rule is applied once to the observed data and held fixed across
##     bootstrap replicates (features are fixed under the bootstrap).
## (c) CVIA078 scale: n = 117, event rate 0.44, p = 9, beta_sep = 2,
##     mu_sep in {0, 0.5}, plus a K = 1 null cell; LRT (bootstrap, B = 19)
##     under normalize = FALSE and under the adaptive rule.
adaptive_choice <- function(d, cfg) {
  cvj <- suppressWarnings(gemcox_cv_loglik(
    d$X, time = d$time, status = d$status, K = 1:2, nfolds = 5, criterion = "joint",
    seed = d$seed, gamma = cfg$gamma, lambda = cfg$lambda, alpha = cfg$alpha,
    normalize_gmm_by_dim = cfg$normalize))
  unname(which.max(attr(cvj, "components")["gmm", ])) == 1
}
lrt_boot <- function(d, s, cfg, normalize) {
  ht <- suppressWarnings(gemcox_heterogeneity_test(
    X_gmm = d$X, time = d$time, status = d$status, statistic = "lrt",
    null = "bootstrap", B = B, seed = s, gamma = cfg$gamma, lambda = cfg$lambda,
    alpha = cfg$alpha, normalize_gmm_by_dim = normalize))
  c(stat = unname(ht$statistic), p = ht$p.value, valid = ht$n_valid)
}
e2b_test_one <- function(r, cell_id, design, cfg, seed_exp, exp_id, arms) {
  pin_rng()
  s <- mk_seed(seed_exp, cell_id, r)
  d <- simulate(design, s)
  d$seed <- s
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch({
    nz <- adaptive_choice(d, cfg)
    out <- list(adaptive_normalize = nz)
    if ("adaptive" %in% arms) out$adaptive <- lrt_boot(d, s, cfg, nz)
    if ("false" %in% arms) out$false <- lrt_boot(d, s, cfg, FALSE)
    out
  }, error = function(e) e)
  base <- data.frame(exp = exp_id, cell = cell_id, rep = r, seed = s, role = design$role,
                     n = design$n, p = design$p, mu_sep = design$mu_sep,
                     beta_sep = design$beta_sep, rho = design$rho, K_true = design$K_true,
                     events = sum(d$status), secs = proc.time()[["elapsed"]] - t0,
                     stringsAsFactors = FALSE)
  if (inherits(res, "error")) {
    return(cbind(base, failed = TRUE, error = conditionMessage(res),
                 adaptive_normalize = NA, stat_adaptive = NA_real_, p_adaptive = NA_real_,
                 valid_adaptive = NA_real_, stat_false = NA_real_, p_false = NA_real_,
                 valid_false = NA_real_))
  }
  g <- function(arm, k) if (is.null(res[[arm]])) NA_real_ else unname(res[[arm]][k])
  cbind(base, failed = FALSE, error = NA_character_,
        adaptive_normalize = res$adaptive_normalize,
        stat_adaptive = g("adaptive", "stat"), p_adaptive = g("adaptive", "p"),
        valid_adaptive = g("adaptive", "valid"), stat_false = g("false", "stat"),
        p_false = g("false", "p"), valid_false = g("false", "valid"))
}
E2B_METHODS <- list(
  "GeM-Cox (gamma=1)"                 = METHODS[["GeM-Cox (gamma=1)"]],
  "GeM-Cox (gamma=1, adaptive)"       = function(d, te, cfg) m_adaptive(d, te, cfg),
  "GeM-Cox (gamma=1, normalize=TRUE)" = function(d, te, cfg) m_gemcox(d, te, cfg, 1, TRUE),
  "GeM-Cox (gamma=0)"                 = METHODS[["GeM-Cox (gamma=0)"]],
  "Oracle (true labels)"              = METHODS[["Oracle (true labels)"]])
e2bc_cells <- data.frame(role = c("type I", "power", "power"), n = 117, p = 9,
                         mu_sep = c(0, 0, 0.5), beta_sep = c(0, 2, 2), rho = 0,
                         K_true = c(1, 2, 2), event_rate = 0.44, reps = 200)
e2bc_cells$reps <- pmax(2L, round(e2bc_cells$reps * E5_SCALE))
if ("E2b" %in% RUN) {
  cells <- expand.grid(n = 800, p = 10, mu_sep = c(0, 0.5, 1, 2), beta_sep = c(1, 2, 3),
                       KEEP.OUT.ATTRS = FALSE)
  run_experiment(21, "E2b(a): adaptive vs FALSE on E2's datasets", cells, REPS,
                 "E2b_a_adaptive.rds", methods = E2B_METHODS, seed_exp = 2)

  run_tests <- function(cells, exp_id, seed_exp, arms, file, label) {
    cat(sprintf("\n[E%d] %s (%d cells, B = %d)\n", exp_id, label, nrow(cells), B))
    t0 <- Sys.time()
    X <- bind_rows(lapply(seq_len(nrow(cells)), function(j) {
      design <- as.list(cells[j, ])
      t1 <- Sys.time()
      x <- bind_rows(par_lapply(seq_len(design$reps), e2b_test_one, cell_id = j,
                                design = design, cfg = CFG, seed_exp = seed_exp,
                                exp_id = exp_id, arms = arms))
      cat(sprintf("  cell %2d %-6s n=%d mu=%g beta=%g: %d datasets, %d failed, chose TRUE %.2f, reject(adaptive/FALSE) %.3f/%.3f (%.1f min)\n",
                  j, design$role, design$n, design$mu_sep, design$beta_sep, nrow(x),
                  sum(x$failed), mean(x$adaptive_normalize, na.rm = TRUE),
                  mean(x$p_adaptive <= 0.05, na.rm = TRUE),
                  mean(x$p_false <= 0.05, na.rm = TRUE),
                  as.numeric(difftime(Sys.time(), t1, units = "mins"))))
      x
    }))
    save_results(X, file, t0)
  }
  run_tests(e5a_cells, 22, 5, "adaptive", "E2b_b_lrt_adaptive.rds",
            "E2b(b): LRT under the adaptive rule on E5a/E5c's datasets")
  run_tests(e2bc_cells, 23, 23, c("adaptive", "false"), "E2b_c_cvia078.rds",
            "E2b(c): LRT at CVIA078 scale, FALSE and adaptive")
}

## ================= ADDED AFTER sim-freeze-v1 ==========================
## E5d: type I error of the default test when profiles separate but
## mechanisms do not (profile-only null). Pre-specified in README.md.
e5d_cells <- data.frame(role = "type I (profile-only)", n = 400, p = 10,
                        mu_sep = c(0.5, 1, 2, 1), beta_sep = 0,
                        rho = c(0, 0, 0, 0.5), K_true = 2, reps = c(500, 500, 500, 200))
e5d_cells$reps <- pmax(2L, round(e5d_cells$reps * E5_SCALE))
e5d_one <- function(r, cell_id, design, cfg) {
  pin_rng()
  s <- mk_seed(24, cell_id, r)
  d <- simulate(design, s)
  t0 <- proc.time()[["elapsed"]]
  n_warn <- 0L
  res <- tryCatch(withCallingHandlers(
    gemcox_heterogeneity_test(X_gmm = d$X, time = d$time, status = d$status,
                              statistic = "lrt", null = "bootstrap", B = B, seed = s,
                              gamma = cfg$gamma, lambda = cfg$lambda, alpha = cfg$alpha,
                              normalize_gmm_by_dim = cfg$normalize),
    warning = function(w) { n_warn <<- n_warn + 1L; invokeRestart("muffleWarning") }),
    error = function(e) e)
  base <- data.frame(exp = 24, cell = cell_id, rep = r, seed = s, role = design$role,
                     n = design$n, p = design$p, mu_sep = design$mu_sep,
                     beta_sep = design$beta_sep, rho = design$rho, K_true = design$K_true,
                     events = sum(d$status), secs = proc.time()[["elapsed"]] - t0,
                     n_warnings = n_warn, stringsAsFactors = FALSE)
  if (inherits(res, "error")) {
    return(cbind(base, failed = TRUE, error = conditionMessage(res), stat = NA_real_,
                 p_bootstrap = NA_real_, valid_bootstrap = NA_integer_))
  }
  cbind(base, failed = FALSE, error = NA_character_, stat = unname(res$statistic),
        p_bootstrap = res$p.value, valid_bootstrap = res$n_valid)
}
if ("E5d" %in% RUN) {
  cat(sprintf("\n[E5d] profile-only null: LRT calibration (%d cells, B = %d)\n",
              nrow(e5d_cells), B))
  t0 <- Sys.time()
  E5d <- bind_rows(lapply(seq_len(nrow(e5d_cells)), function(j) {
    design <- as.list(e5d_cells[j, ])
    t1 <- Sys.time()
    x <- bind_rows(par_lapply(seq_len(design$reps), e5d_one, cell_id = j,
                              design = design, cfg = CFG))
    cat(sprintf("  cell %d mu=%g rho=%g: %d datasets, %d failed, reject %.3f (%.1f min)\n",
                j, design$mu_sep, design$rho, nrow(x), sum(x$failed),
                mean(x$p_bootstrap <= 0.05, na.rm = TRUE),
                as.numeric(difftime(Sys.time(), t1, units = "mins"))))
    x
  }))
  save_results(E5d, "E5d_profile_null.rds", t0)
}

## E5e: selecting K (K grid 1:4) by four selectors, K_true in {1, 2, 3}.
## Pre-specified in README.md. SEQ_REPS limits the sequential LRT to the
## first SEQ_REPS datasets per cell (set from the timing run; see README).
SEQ_REPS <- as.integer(Sys.getenv("GEMCOX_SEQ_REPS", "200"))
e5e_cells <- local({
  main <- expand.grid(K_true = c(2, 3), regime = c("profile only", "mechanism only", "both"),
                      separation = c("moderate", "strong"), n = c(400, 800),
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  mu <- ifelse(main$separation == "moderate", 1, 3)
  be <- ifelse(main$separation == "moderate", 2, 3)
  main$mu_sep <- ifelse(main$regime == "mechanism only", 0, mu)
  main$beta_sep <- ifelse(main$regime == "profile only", 0, be)
  main$set <- "main"; main$p <- 10; main$event_rate <- 0.55
  k1 <- data.frame(K_true = 1, regime = "null", separation = "none", n = c(400, 800),
                   mu_sep = 0, beta_sep = 0, set = "main", p = 10, event_rate = 0.55,
                   stringsAsFactors = FALSE)
  cv <- data.frame(K_true = c(1, 2), regime = c("null", "both"),
                   separation = c("none", "moderate"), n = 117, mu_sep = c(0, 1),
                   beta_sep = c(0, 2), set = "CVIA078 scale", p = 9, event_rate = 0.44,
                   stringsAsFactors = FALSE)
  x <- rbind(k1, main, cv)
  x$rho <- 0; x$reps <- pmax(2L, round(200 * E5_SCALE))
  x
})
E5E_KMAX <- 4L
e5e_bic_df <- function(K, q, p) (K - 1) + 2 * K * q + K * p
e5e_one <- function(r, cell_id, design, cfg, seq_reps = SEQ_REPS) {
  pin_rng()
  s <- mk_seed(25, cell_id, r)
  d <- simulate(design, s)
  d$seed <- s
  t0 <- proc.time()[["elapsed"]]
  fit_args <- list(gamma = cfg$gamma, lambda = cfg$lambda, alpha = cfg$alpha,
                   normalize_gmm_by_dim = cfg$normalize)
  Ks <- seq_len(E5E_KMAX)
  ## selectors 1 and 2: joint and partial cross-validated criteria (one call)
  cvj <- tryCatch(suppressWarnings(do.call(gemcox_cv_loglik, c(list(
    d$X, time = d$time, status = d$status, K = Ks, nfolds = 5, criterion = "joint",
    seed = s), fit_args))), error = function(e) e)
  cv_ok <- !inherits(cvj, "error")
  joint <- if (cv_ok) as.numeric(cvj) else rep(NA_real_, E5E_KMAX)
  partial <- if (cv_ok) as.numeric(attr(cvj, "components")["partial", ]) else rep(NA_real_, E5E_KMAX)
  ## selector 3: BIC of the full-data mixture fits
  ll <- vapply(Ks, function(K) {
    f <- tryCatch(suppressWarnings(do.call(gemcox, c(list(d$X, time = d$time,
                                                          status = d$status, K = K), fit_args))),
                  error = function(e) NULL)
    if (is.null(f)) NA_real_ else f$final_score
  }, 0)
  q <- ncol(d$X)
  bic <- -2 * ll + e5e_bic_df(Ks, q, q) * log(length(d$time))
  ## selector 4: sequential bootstrap LRT (K0 + 1 vs K0, stop at first p > 0.05)
  seq_p <- rep(NA_real_, E5E_KMAX - 1)
  K_seq <- NA_integer_; seq_failed <- NA
  if (r <= seq_reps) {
    K_seq <- E5E_KMAX; seq_failed <- FALSE
    for (k0 in seq_len(E5E_KMAX - 1)) {
      ht <- tryCatch(suppressWarnings(do.call(gemcox_heterogeneity_test, c(list(
        X_gmm = d$X, time = d$time, status = d$status, statistic = "lrt",
        null = "bootstrap", B = B, K0 = k0, seed = s + k0), fit_args))),
        error = function(e) e)
      if (inherits(ht, "error")) { K_seq <- k0; seq_failed <- TRUE; break }
      seq_p[k0] <- ht$p.value
      if (ht$p.value > 0.05) { K_seq <- k0; break }
    }
  }
  pick_max <- function(v) if (any(is.finite(v))) which.max(v) else NA_integer_
  pick_min <- function(v) if (any(is.finite(v))) which.min(v) else NA_integer_
  data.frame(exp = 25, cell = cell_id, rep = r, seed = s, set = design$set,
             K_true = design$K_true, regime = design$regime, separation = design$separation,
             n = design$n, p = design$p, mu_sep = design$mu_sep, beta_sep = design$beta_sep,
             events = sum(d$status), secs = proc.time()[["elapsed"]] - t0,
             cv_failed = !cv_ok,
             K_joint = pick_max(joint), K_partial = pick_max(partial), K_bic = pick_min(bic),
             K_seq = K_seq, seq_failed = seq_failed,
             p_2v1 = seq_p[1], p_3v2 = seq_p[2], p_4v3 = seq_p[3],
             joint_1 = joint[1], joint_2 = joint[2], joint_3 = joint[3], joint_4 = joint[4],
             partial_1 = partial[1], partial_2 = partial[2], partial_3 = partial[3],
             partial_4 = partial[4], loglik_1 = ll[1], loglik_2 = ll[2], loglik_3 = ll[3],
             loglik_4 = ll[4], stringsAsFactors = FALSE)
}
if ("E5e" %in% RUN) {
  cat(sprintf("\n[E5e] selecting K (grid 1:%d; %d cells; sequential LRT on the first %d datasets per cell)\n",
              E5E_KMAX, nrow(e5e_cells), SEQ_REPS))
  t0 <- Sys.time()
  E5e <- bind_rows(lapply(seq_len(nrow(e5e_cells)), function(j) {
    design <- as.list(e5e_cells[j, ])
    t1 <- Sys.time()
    x <- bind_rows(par_lapply(seq_len(design$reps), e5e_one, cell_id = j,
                              design = design, cfg = CFG))
    cat(sprintf("  cell %2d %-13s K=%d %-14s %-8s n=%4d: correct joint/partial/BIC/seq %.2f/%.2f/%.2f/%.2f (%.1f min)\n",
                j, design$set, design$K_true, design$regime, design$separation, design$n,
                mean(x$K_joint == design$K_true, na.rm = TRUE),
                mean(x$K_partial == design$K_true, na.rm = TRUE),
                mean(x$K_bic == design$K_true, na.rm = TRUE),
                mean(x$K_seq == design$K_true, na.rm = TRUE),
                as.numeric(difftime(Sys.time(), t1, units = "mins"))))
    x
  }))
  save_results(E5e, "E5e_kselect.rds", t0)
}

## E5f: where the adaptive test becomes calibrated (optional; added after
## sim-freeze-v1). Pre-specified in README.md. Same datasets for the
## adaptive and the default (normalize = FALSE) LRT, bootstrap null, B = 19.
## Null cells are K_true = 1 with a nonzero shared coefficient, as in the
## E2b(c) n = 117 row they are compared with. Documentation only.
e5f_cells <- data.frame(role = rep(c("type I", "power"), 2), n = rep(c(200, 300), each = 2),
                        p = 10, mu_sep = 0, beta_sep = rep(c(0, 2), 2), rho = 0,
                        K_true = rep(c(1, 2), 2), event_rate = 0.44, reps = 200)
e5f_cells$reps <- pmax(2L, round(e5f_cells$reps * E5_SCALE))
if ("E5f" %in% RUN) {
  cat(sprintf("\n[E5f] adaptive vs default LRT at n = 200, 300 (%d cells, B = %d)\n",
              nrow(e5f_cells), B))
  t0 <- Sys.time()
  E5f <- bind_rows(lapply(seq_len(nrow(e5f_cells)), function(j) {
    design <- as.list(e5f_cells[j, ])
    t1 <- Sys.time()
    x <- bind_rows(par_lapply(seq_len(design$reps), e2b_test_one, cell_id = j,
                              design = design, cfg = CFG, seed_exp = 26, exp_id = 26,
                              arms = c("adaptive", "false")))
    cat(sprintf("  cell %d %-6s n=%d beta=%g: %d datasets, %d failed, chose TRUE %.2f, reject(adaptive/FALSE) %.3f/%.3f (%.1f min)\n",
                j, design$role, design$n, design$beta_sep, nrow(x), sum(x$failed),
                mean(x$adaptive_normalize, na.rm = TRUE),
                mean(x$p_adaptive <= 0.05, na.rm = TRUE), mean(x$p_false <= 0.05, na.rm = TRUE),
                as.numeric(difftime(Sys.time(), t1, units = "mins"))))
    x
  }))
  save_results(E5f, "E5f_adaptive_scale.rds", t0)
}

## ================= S1: normalisation (supplement) =====================
if ("S1" %in% RUN) {
  meths <- list(
    "GeM-Cox, normalize = FALSE" = function(d, te, cfg) m_gemcox(d, te, cfg, 1, FALSE),
    "GeM-Cox, normalize = TRUE"  = function(d, te, cfg) m_gemcox(d, te, cfg, 1, TRUE),
    "Oracle (true labels)"       = METHODS[["Oracle (true labels)"]])
  cells <- expand.grid(n = 800, p = c(10, 80),
                       mu_sep = c(0, 2), beta_sep = 2, KEEP.OUT.ATTRS = FALSE)
  run_experiment(9, "normalisation is regime-dependent", cells, REPS,
                 "S1_normalisation.rds", methods = meths)
}

cat("\nAll requested experiments complete. Next: 03_summarise.R\n")
