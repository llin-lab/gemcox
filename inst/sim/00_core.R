###############################################################
## 00_core.R  --  shared foundation for every simulation
##
## Run every script from this directory (gemcox/inst/sim), with the gemcox
## package installed:  R CMD INSTALL gemcox;  cd gemcox/inst/sim;
##                      Rscript 01_pilot.R
##
## One data-generating mechanism, one set of estimators, one set of
## metrics. The DGP is gemcox::gemcox_simulate(); nothing here or in the
## experiment scripts redefines it. Membership weights for new subjects come
## from predict.gemcox() only.
##
## ASCII only.
###############################################################

suppressMessages({
  library(gemcox); library(survival); library(mclust)
  library(dplyr); library(tidyr); library(future.apply)
})

## ================= CONFIG =============================================
CFG <- list(
  ## estimator settings, fixed across every experiment
  lambda    = 0.05,      # small ridge; lambda = 0 gave all-zero components at p = 80
  alpha     = 0,         # ridge, not elastic net
  gamma     = 1,         # generative value; no data-driven criterion selects it defensibly
  normalize = FALSE,     # varied only in the supplementary experiment S1
  ## design defaults
  event_rate = 0.55,
  n_test     = 2000,
  K_fit      = 2,
  ## replication
  reps_pilot = 20L,
  reps_full  = 200L,
  seed_base  = 20260101L,
  ## parallel workers (PSOCK / multisession); 1 = serial
  workers    = as.integer(Sys.getenv("GEMCOX_WORKERS",
                                     max(1L, parallel::detectCores() - 1L)))
)

## EM settings of the frozen study. The package default became tol = 1e-8,
## max_iter = 1000 in gemcox 0.2.0. The pipeline keeps the settings it was
## frozen with, so stored results reproduce; scripts that run converged
## fits (C_methods.R, 07_, 09_) pass tol and max_iter explicitly. They are
## set here and again in every worker task (par_lapply).
FROZEN_EM <- list(gemcox.tol = 1e-4, gemcox.max_iter = 100L)
options(FROZEN_EM)

## GEMCOX_RESULTS redirects output (e.g. "rehearsal/results") so a
## rehearsal never overwrites the manuscript results; figures go alongside.
RESULTS <- Sys.getenv("GEMCOX_RESULTS", "results")
FIGURES <- file.path(dirname(RESULTS), "figures")
dir.create(RESULTS, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURES, showWarnings = FALSE, recursive = TRUE)

## RANDOM NUMBER GENERATOR. Every random draw in the pipeline (simulated
## data, initialisation, folds, null replicates, mclust subsampling) uses
## L'Ecuyer-CMRG, the generator future.seed = TRUE installs in workers.
## set.seed(s) gives different streams under different generators, so each
## per-dataset function pins it; a replicate then reproduces from
## (experiment, cell, replicate) in any R session. NOTE: calling
## gemcox_simulate(seed = s) in a session using the default Mersenne-Twister
## generator gives DIFFERENT data from the pipeline's.
pin_rng <- function() RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")

## deterministic, collision-free seeds (cells < 500, reps < 1000)
mk_seed <- function(exp_id, cell, rep)
  CFG$seed_base + exp_id * 1e6 + cell * 1e3 + rep

## the one DGP, with the pipeline's event rate
simulate <- function(design, seed, n = design$n) {
  gemcox_simulate(n = n, p = design$p, mu_sep = design$mu_sep,
                  beta_sep = design$beta_sep,
                  rho = if (is.null(design$rho)) 0 else design$rho,
                  K_true = if (is.null(design$K_true)) 2 else design$K_true,
                  event_rate = if (is.null(design$event_rate)) CFG$event_rate
                               else design$event_rate,
                  ## added after sim-freeze-v2 for C5; absent means "gaussian",
                  ## the frozen DGP
                  features = if (is.null(design$features)) "gaussian" else design$features,
                  ## added for E6; absent means the frozen random censoring
                  censoring = if (is.null(design$censoring)) "random" else design$censoring,
                  followup = design$followup,
                  seed = seed)
}

## ================= PARALLEL SETUP =====================================
## PSOCK/multisession workers, never forked: glmnet uses OpenMP and forked
## (mclapply) workers previously failed to return. OMP_NUM_THREADS = 1 is
## set before the workers start, so they inherit it, and again inside each
## task. Results do not depend on the worker count: every random draw is
## seeded from (experiment, cell, replicate), and future.seed = TRUE gives
## reproducible streams for anything else.
setup_parallel <- function(workers = CFG$workers) {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  PLAN_NAME <<- if (workers <= 1) "sequential" else "multisession"
  if (workers <= 1) future::plan(future::sequential)
  else future::plan(future::multisession, workers = workers)
  invisible(workers)
}
par_lapply <- function(X, FUN, ...) {
  future_lapply(X, function(x, ...) {
    Sys.setenv(OMP_NUM_THREADS = "1")
    options(FROZEN_EM)
    FUN(x, ...)
  }, ..., future.seed = TRUE,
  ## mclust must be ATTACHED in workers: Mclust() evaluates mclustBIC() by
  ## name in the caller's environment and fails when it is only loaded.
  future.packages = c("gemcox", "survival", "mclust"))
}

## ================= METRICS ============================================
## Primary estimand: the DIRECTION of the between-cluster coefficient
## contrast. For general K the target is the span of the centred
## coefficient vectors, compared by principal angles; at K = 2 this reduces
## to the absolute cosine between the two contrast vectors.
span_basis <- function(bl) {
  if (length(bl) < 2) return(NULL)
  B <- do.call(cbind, bl); B <- B - rowMeans(B)
  s <- svd(B); k <- s$d > 1e-8 * max(s$d, 1e-12)
  if (!any(k)) return(NULL)
  s$u[, k, drop = FALSE]
}
contrast_recovery <- function(bhat, btrue, p) {
  A <- span_basis(bhat); B <- span_basis(btrue)
  if (is.null(A) || is.null(B)) return(random_floor(p))  # no contrast => chance
  d <- min(ncol(A), ncol(B))
  mean(pmin(pmax(svd(t(A) %*% B)$d[seq_len(d)], 0), 1))
}
## expected |cos| between two random directions in R^p
random_floor <- function(p) sqrt(2 / (pi * p))
pct_of_oracle <- function(x, oracle, p) {
  fl <- random_floor(p)
  100 * (x - fl) / (oracle - fl)
}
ari <- function(a, b) mclust::adjustedRandIndex(a, b)
mcse <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))

## MSE of the linear predictor up to a constant: the overall level of a
## Cox linear predictor is not identified separately from the baseline, so
## the mean of the difference is removed before squaring.
mse_lp <- function(pred, truth) {
  e <- pred - truth
  mean((e - mean(e))^2)
}

## Cumulative hazard of a fitted step-function baseline, via stepfun
## (independent of the package internals).
cumhaz_at <- function(b, t) {
  if (!length(b$time)) return(rep(0, length(t)))
  stats::stepfun(b$time, c(0, b$cumhaz))(t)
}

## ================= ESTIMATORS =========================================
## Linear predictors on the test set are UNCENTRED mixtures
## sum_k w_k x' beta_k, compared with the true eta by mse_lp().
safe_cox <- function(time, status, X, idx = seq_along(time)) {
  if (length(idx) < 10 || sum(status[idx]) < 5) return(rep(0, ncol(X)))
  f <- suppressWarnings(tryCatch(
    coxph(Surv(time[idx], status[idx]) ~ X[idx, , drop = FALSE]),
    error = function(e) NULL))
  if (is.null(f)) return(rep(0, ncol(X)))
  b <- as.numeric(coef(f)); b[!is.finite(b)] <- 0
  if (max(abs(b)) > 20) b <- rep(0, ncol(X))   # divergent fit scored as null
  b
}

m_single <- function(d, te, cfg) {           # standard practice
  b <- safe_cox(d$time, d$status, d$X)
  list(beta_list = list(b, b), eta = as.numeric(te$X %*% b), tau = NULL)
}
## Two-stage comparator: diagonal, varying-volume Gaussian mixture (VVI)
## with mclust's conjugate prior, in EVERY cell, so it is one method along
## the whole p-axis and matches GeM-Cox's regularised diagonal covariance.
## Without the prior VVI returned no fit at p >= 40 (n = 800).
m_twostage <- function(d, te, cfg) {         # cluster antibodies, then model
  m <- tryCatch(mclust::Mclust(d$X, G = cfg$K_fit, modelNames = "VVI",
                               prior = mclust::priorControl(), verbose = FALSE),
                error = function(e) NULL)
  if (is.null(m)) stop("Mclust failed")
  bl <- lapply(seq_len(cfg$K_fit), function(k)
    safe_cox(d$time, d$status, d$X, which(m$classification == k)))
  z <- predict(m, te$X)$z
  list(beta_list = bl, eta = rowSums(z * (te$X %*% do.call(cbind, bl))),
       tau = m$z, labels = m$classification)
}
m_gemcox <- function(d, te, cfg, gamma = NULL, normalize = NULL, temp = 1) {
  g  <- if (is.null(gamma)) cfg$gamma else gamma
  nz <- if (is.null(normalize)) cfg$normalize else normalize
  n_warn <- 0L
  f <- withCallingHandlers(
    gemcox(d$X, time = d$time, status = d$status, K = cfg$K_fit, gamma = g,
           lambda = cfg$lambda, alpha = cfg$alpha, normalize_gmm_by_dim = nz,
           temp = temp, seed = 1),
    warning = function(w) { n_warn <<- n_warn + 1L; invokeRestart("muffleWarning") })
  tn <- predict(f, te$X, type = "tau")
  list(beta_list = lapply(seq_len(ncol(f$beta)), function(k) f$beta[, k]),
       eta = rowSums(tn * (te$X %*% coef(f))),
       tau = f$tau, labels = f$cluster, fit = f,
       diag = data.frame(converged = f$converged, iterations = f$iterations,
                         guard_bound = nrow(f$guards) > 0, n_warnings = n_warn))
}
m_oracle <- function(d, te, cfg) {           # true labels: the ceiling
  bl <- lapply(seq_len(max(d$Z)), function(k)
    safe_cox(d$time, d$status, d$X, which(d$Z == k)))
  eta <- rowSums(sapply(seq_along(bl), function(k)
    (te$Z == k) * (te$X %*% bl[[k]])))
  list(beta_list = bl, eta = eta, tau = NULL, labels = d$Z)
}

METHODS <- list(
  "Single Cox model"      = function(d, te, cfg) m_single(d, te, cfg),
  "Two-stage (GMM+Cox)"   = function(d, te, cfg) m_twostage(d, te, cfg),
  "GeM-Cox (gamma=0)"     = function(d, te, cfg) m_gemcox(d, te, cfg, gamma = 0),
  "GeM-Cox (gamma=1)"     = function(d, te, cfg) m_gemcox(d, te, cfg, gamma = 1),
  "Oracle (true labels)"  = function(d, te, cfg) m_oracle(d, te, cfg)
)
METHOD_LEVELS <- names(METHODS)

## ================= CLUSTER DEGENERACY =================================
## A fitted cluster is DEGENERATE if it has fewer than 10 subjects or fewer
## than 5 events (hard labels: the assigned cluster, or argmax tau for
## GeM-Cox). A single Cox model is one cluster, the whole sample. Recorded
## for every method; 03_summarise.R excludes degenerate fits in the primary
## analysis and scores them at the random-direction floor in a sensitivity
## analysis. Also recorded: tau-weighted sizes/events where the method has
## membership weights, and how many cluster fits returned all-zero
## coefficients (safe_cox scores degenerate or divergent fits as zero).
DEGEN_MIN_SIZE <- 10; DEGEN_MIN_EVENTS <- 5
cluster_stats <- function(m, d) {
  K <- length(m$beta_list)
  lab <- if (is.null(m$labels)) factor(rep(1L, length(d$time))) else
    factor(m$labels, levels = seq_len(K))
  size <- as.integer(table(lab))
  ev <- as.integer(tapply(d$status, lab, sum)); ev[is.na(ev)] <- 0L
  eff_size <- eff_ev <- NA_real_
  if (!is.null(m$tau)) {
    eff_size <- min(colSums(m$tau)); eff_ev <- min(colSums(m$tau * d$status))
  }
  data.frame(min_size = min(size), min_events = min(ev),
             degenerate = any(size < DEGEN_MIN_SIZE | ev < DEGEN_MIN_EVENTS),
             min_eff_size = eff_size, min_eff_events = eff_ev,
             n_zero_fits = sum(vapply(m$beta_list, function(b) all(b == 0), NA)))
}

## ================= HARNESS ============================================
## Raw per-replicate rows only. Summaries and Monte Carlo SEs are computed
## in 03_summarise.R from the saved rows. A method that errors is recorded
## as a failed row (metrics NA, error message kept), never dropped.
## seed_exp: the experiment id used for SEEDS (defaults to exp_id). E1b sets
## it to 1 so that it analyses exactly E1's datasets.
one_rep <- function(r, exp_id, cell_id, design, methods, cfg, extra, seed_exp = exp_id) {
  pin_rng()
  s  <- mk_seed(seed_exp, cell_id, r)
  d  <- simulate(design, s)
  d$seed <- s
  te <- simulate(design, s + 500000L, n = cfg$n_test)
  rows <- lapply(names(methods), function(nm) {
    t0 <- proc.time()[["elapsed"]]
    ## Seed the stream from the dataset seed before EVERY method: some
    ## methods draw random numbers (mclust initialises from a random subset
    ## when n > 2000), and future.seed streams are not tied to (exp, cell, rep).
    set.seed(s)
    m <- tryCatch(methods[[nm]](d, te, cfg), error = function(e) e)
    secs <- proc.time()[["elapsed"]] - t0
    base <- data.frame(
      exp = exp_id, cell = cell_id, rep = r, seed = s, method = nm,
      n = design$n, p = design$p, mu_sep = design$mu_sep,
      beta_sep = design$beta_sep,
      rho = if (is.null(design$rho)) 0 else design$rho,
      events = sum(d$status), secs = secs, stringsAsFactors = FALSE)
    if (inherits(m, "error")) {
      return(cbind(base, failed = TRUE, error = conditionMessage(m),
                   recovery = NA_real_, ari = NA_real_, mse_eta = NA_real_,
                   sharpness = NA_real_, contrast_norm = NA_real_,
                   min_size = NA_integer_, min_events = NA_integer_, degenerate = NA,
                   min_eff_size = NA_real_, min_eff_events = NA_real_,
                   n_zero_fits = NA_integer_))
    }
    row <- cbind(base, failed = FALSE, error = NA_character_,
      recovery = contrast_recovery(m$beta_list, d$beta_list, design$p),
      ari = if (is.null(m$labels)) NA_real_ else ari(m$labels, d$Z),
      mse_eta = mse_lp(m$eta, te$eta),
      sharpness = if (is.null(m$tau)) NA_real_ else mean(apply(m$tau, 1, max)),
      contrast_norm = sqrt(sum((m$beta_list[[1]] -
                                m$beta_list[[length(m$beta_list)]])^2)),
      cluster_stats(m, d))
    if (!is.null(m$diag)) row <- cbind(row, m$diag)
    if (!is.null(extra)) row <- cbind(row, extra(d, te, m, cfg))
    row
  })
  bind_rows(rows)
}

run_cell <- function(exp_id, cell_id, design, reps, methods = METHODS,
                     cfg = CFG, extra = NULL, seed_exp = exp_id) {
  bind_rows(par_lapply(seq_len(reps), one_rep, exp_id = exp_id, cell_id = cell_id,
                       design = design, methods = methods, cfg = cfg, extra = extra,
                       seed_exp = seed_exp))
}

## Provenance stored in every results file.
run_meta <- function(t0) {
  git <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
                  error = function(e) NA_character_, warning = function(w) NA_character_)
  ## generated outputs (results/, figures/) are excluded: only code counts
  dirty <- tryCatch(length(system2("git", c("status", "--porcelain", "--untracked-files=no", "--",
                                            "../..", "':(exclude)results'", "':(exclude)figures'"),
                                   stdout = TRUE, stderr = FALSE)) > 0,
                    error = function(e) NA, warning = function(w) NA)
  list(created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
       gemcox_version = as.character(utils::packageVersion("gemcox")),
       git_commit = if (length(git)) git[1] else NA_character_,
       git_dirty = dirty,   # TRUE: tracked code under gemcox/ (not results/) had uncommitted changes
       config = CFG,
       plan = PLAN_NAME,
       workers = future::nbrOfWorkers(),
       wall_minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")),
       sessionInfo = utils::capture.output(utils::sessionInfo()))
}
save_results <- function(rows, file, t0) {
  saveRDS(list(rows = rows, meta = run_meta(t0)), file.path(RESULTS, file))
}
read_results <- function(file) {
  p <- file.path(RESULTS, file)
  if (!file.exists(p)) { cat("  [skip]", file, "\n"); return(NULL) }
  readRDS(p)$rows
}

run_experiment <- function(exp_id, label, cells, reps, file, ...) {
  cat(sprintf("\n[E%d] %s  (%d cells x %d reps, %d workers)\n", exp_id, label,
              nrow(cells), reps, future::nbrOfWorkers()))
  t0 <- Sys.time()
  res <- bind_rows(lapply(seq_len(nrow(cells)), function(j) {
    design <- as.list(cells[j, , drop = FALSE])
    design <- design[!vapply(design, is.na, TRUE)]
    t1 <- Sys.time()
    x <- run_cell(exp_id, j, design, reps, ...)
    cat(sprintf("  cell %2d: %s  %d rows, %d failed  (%.1f min)\n", j,
                paste(sprintf("%s=%s", names(design), unlist(design)), collapse = " "),
                nrow(x), sum(x$failed),
                as.numeric(difftime(Sys.time(), t1, units = "mins"))))
    x
  }))
  save_results(res, file, t0)
  cat(sprintf("  saved %s  (%.1f min)\n", file,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  invisible(res)
}

setup_parallel()
cat("00_core.R loaded.\n")
cat(sprintf("  gemcox %s; gamma=%g lambda=%g alpha=%g normalize=%s event_rate=%.2f; %d workers\n",
            utils::packageVersion("gemcox"), CFG$gamma, CFG$lambda, CFG$alpha,
            CFG$normalize, CFG$event_rate, future::nbrOfWorkers()))
