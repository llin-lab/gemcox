###############################################################
## 09_tolerance_study.R  --  evidence for a proposed default EM tolerance
## (added after sim-freeze-v2; a diagnostic for a proposal, NOT applied to
## the package; the rule below was fixed in README.md before it was run)
##
## GeM-Cox (gamma = 1, all other settings at the package defaults) is fitted
## at tolerances 1e-4 (the current default), 1e-5, 1e-6, 1e-7 and 1e-8, each
## with at most 1000 iterations, on the first 20 datasets of each cell
## below (10 for E4a). The 1e-8 fit is the reference.
##   E1 (C1) 4 cells; E2 (C2) 12 cells; E3b (C3b) p = 40, 80; C4 2 cells;
##   C5 2 cells; E4a 4 cells (n 800-6400).
## Recorded per fit: seconds, iterations, convergence, contrast recovery,
## the log-likelihood gap to the reference, and the coefficient distance to
## the reference (max over clusters of |beta_k - beta_ref,k|, relative to
## the largest |beta_ref,k|, clusters matched).
##
## Proposal rule: the loosest tolerance at which at least 99% of fits,
## pooled, and at least 95% in every cell, have recovery within 0.01 of the
## reference; with max_iter the smallest of 200, 500, 1000 within which at
## least 99.5% of fits at that tolerance converge. If no tolerance looser
## than 1e-8 qualifies, 1e-8 is proposed. Timing: median and 90th
## percentile seconds per fit, relative to the current default.
##
##   Rscript 09_tolerance_study.R     (sourced by 07's timing mode, it runs nothing)
## ASCII only.
###############################################################

if (!exists("CFG")) { source("00_core.R"); source("C_common.R"); source("C_methods.R") }
TOLS <- c(1e-4, 1e-5, 1e-6, 1e-7, 1e-8)
T_MAXIT <- 1000L
T_DESIGN <- rbind(
  data.frame(key = "C1", cell = 1:4, reps = 20),
  data.frame(key = "C2", cell = 1:12, reps = 20),
  data.frame(key = "C3b", cell = 3:4, reps = 20),
  data.frame(key = "C4", cell = 2:3, reps = 20),
  data.frame(key = "C5", cell = 1:2, reps = 20),
  data.frame(key = "E4a", cell = 1:4, reps = 10))
T_DESIGN$reps <- pmin(T_DESIGN$reps, as.integer(Sys.getenv("GEMCOX_T_REPS", "1000")))  # rehearsal only

t_one <- function(r, row) {
  E <- if (row$key == "E4a") R_EXPS$E4a else C_EXPS[[row$key]]
  j <- which(E$cells$cell == row$cell)
  design <- c_design(E, j)
  pin_rng()
  s <- mk_seed(E$seed_exp, row$cell, r)
  d <- simulate(design, s)
  fits <- lapply(TOLS, function(tol) {
    t0 <- proc.time()[["elapsed"]]
    f <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, K = CFG$K_fit,
                                 gamma = 1, lambda = CFG$lambda, alpha = CFG$alpha,
                                 normalize_gmm_by_dim = CFG$normalize, seed = 1,
                                 tol = tol, max_iter = T_MAXIT))
    list(f = f, secs = proc.time()[["elapsed"]] - t0)
  })
  ref <- fits[[length(TOLS)]]$f
  bref <- ref$beta
  scale <- max(sqrt(colSums(bref^2)))
  dist <- function(b) min(vapply(list(1:2, 2:1), function(o)
    max(sqrt(colSums((b[, o, drop = FALSE] - bref)^2))), 0)) / scale
  data.frame(key = row$key, cell = row$cell, rep = r, seed = s, n = design$n, p = design$p,
             mu_sep = design$mu_sep, beta_sep = design$beta_sep, tol = TOLS,
             secs = vapply(fits, `[[`, 0, "secs"),
             iterations = vapply(fits, function(x) x$f$iterations, 0),
             converged = vapply(fits, function(x) isTRUE(x$f$converged), NA),
             recovery = vapply(fits, function(x) contrast_recovery(
               lapply(1:2, function(k) x$f$beta[, k]), d$beta_list, design$p), 0),
             loglik_gap = ref$final_score - vapply(fits, function(x) x$f$final_score, 0),
             coef_dist = vapply(fits, function(x) dist(x$f$beta), 0),
             hash_train = hash_train(d), stringsAsFactors = FALSE)
}

if (sys.nframe() == 0L) {
  jobs <- do.call(rbind, lapply(seq_len(nrow(T_DESIGN)), function(i)
    data.frame(i = i, r = seq_len(T_DESIGN$reps[i]))))
  cat(sprintf("\n[T1] tolerance study: %d datasets x %d tolerances (%d workers)\n",
              nrow(jobs), length(TOLS), future::nbrOfWorkers()))
  t0 <- Sys.time()
  X <- bind_rows(par_lapply(seq_len(nrow(jobs)), function(k)
    t_one(jobs$r[k], T_DESIGN[jobs$i[k], ])))
  X <- X %>% group_by(key, cell, rep) %>%
    mutate(ref_recovery = recovery[tol == min(TOLS)], rec_dev = abs(recovery - ref_recovery),
           ref_converged = converged[tol == min(TOLS)]) %>% ungroup()
  save_results(X, "T1_tolerance.rds", t0)

  q90 <- function(x) unname(stats::quantile(x, 0.9))
  base <- X %>% filter(tol == 1e-4)
  S <- X %>% group_by(tol) %>%
    summarise(fits = n(), within_0.01 = mean(rec_dev <= 0.01),
              worst_cell = min(tapply(rec_dev <= 0.01, paste(key, cell), mean)),
              mean_abs_rec_dev = mean(rec_dev), max_coef_dist = max(coef_dist),
              median_loglik_gap = median(loglik_gap),
              cap_hits = sum(!converged), iters_median = median(iterations),
              over_100 = mean(iterations > 100), over_200 = mean(iterations > 200),
              over_500 = mean(iterations > 500),
              secs_median = median(secs), secs_p90 = q90(secs), .groups = "drop") %>%
    mutate(time_vs_default = secs_median / secs_median[tol == 1e-4])
  cat("\n=== tolerance study (reference: tolerance 1e-8, max_iter 1000) ===\n")
  print(as.data.frame(S %>% mutate(across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)
  cat(sprintf("  reference fits that themselves hit the cap: %d of %d\n",
              sum(!base$ref_converged), nrow(base)))
  ok <- S$tol[S$tol > min(TOLS) & S$within_0.01 >= 0.99 & S$worst_cell >= 0.95]
  prop_tol <- if (length(ok)) max(ok) else min(TOLS)
  at <- X %>% filter(tol == prop_tol)
  prop_maxit <- c(200, 500, 1000)[which(c(mean(at$converged & at$iterations <= 200),
                                          mean(at$converged & at$iterations <= 500),
                                          mean(at$converged)) >= 0.995)[1]]
  if (is.na(prop_maxit)) prop_maxit <- 1000
  cat(sprintf("\n  PROPOSAL (by the pre-stated rule; not applied): tol = %g, max_iter = %d\n",
              prop_tol, prop_maxit))
  cat(sprintf("  at that setting: %.1f%% of fits within 0.01 of the reference (worst cell %.1f%%);\n",
              100 * S$within_0.01[S$tol == prop_tol], 100 * S$worst_cell[S$tol == prop_tol]))
  cat(sprintf("  median %.2f s per fit (%.1fx the current default), 90th percentile %.2f s\n",
              S$secs_median[S$tol == prop_tol], S$time_vs_default[S$tol == prop_tol],
              S$secs_p90[S$tol == prop_tol]))
  byc <- X %>% group_by(key, cell, n, p, tol) %>%
    summarise(within_0.01 = mean(rec_dev <= 0.01), secs_median = median(secs),
              iters_median = median(iterations), .groups = "drop")
  wcsv <- function(x, f) write.csv(x, file.path(RESULTS, f), row.names = FALSE)
  wcsv(S, "tab_T1_tolerance.csv"); wcsv(byc, "tab_T1_tolerance_cells.csv")
}
