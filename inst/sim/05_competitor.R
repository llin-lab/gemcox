###############################################################
## 05_competitor.R  --  latent-class Cox competitor (added after
## sim-freeze-v2; pre-registered in README.md before it was run)
##
##  C1   E1's 4 cells and datasets (mu_sep = 0): mechanism only
##  C2   E2's 12 cells and datasets: main comparison
##  C3a  E3a's 4 cells and datasets: sample size
##  C3b  E3b's 4 cells and datasets: dimension
##  C4   E2b(c)'s two power cells and datasets: CVIA078 scale
##  C5   NEW: E2's mu_sep in {0.5, 1} x beta_sep = 2 with log-normal features
##
## C1-C3 fit only the competitor and are paired with the frozen rows of the
## same datasets; every row records md5 hashes of its training and test
## data, checked against results/C0_dataset_hashes.rds (05_hash_manifest.R)
## before any paired difference is computed. C4 and C5 have no frozen
## recovery results, so every method is fitted. The sensitivity arm (gating
## penalty lambda / 10) runs on C1's cells and on C2's cells 6 and 7.
##
## GEMCOX_ARMS="matched" runs the matched-regularisation comparison (a
## sensitivity analysis, pre-registered in README.md): the non-default
## configurations in C_methods.R on C1, C2 cells 6-7 and C5, tolerance 1e-8,
## saved to <file>_matched.rds.
##
## GEMCOX_ARMS="convergence" runs the pre-registered convergence-sensitivity
## arm instead: the competitor and GeM-Cox (gamma = 1), both with EM
## tolerance 1e-8 and max_iter 1000, on every cell's datasets, saved to
## <file>_convergence.rds. The frozen rows are not touched.
##
##   Rscript 05_hash_manifest.R            # dataset identity, first
##   GEMCOX_RUN="C1,C2" Rscript 05_competitor.R
##   GEMCOX_ARMS=convergence Rscript 05_competitor.R
##   GEMCOX_TIME=1 Rscript 05_competitor.R # time one dataset per cell; saves nothing
## ASCII only.
###############################################################

source("00_core.R")
source("C_common.R")
REPS <- as.integer(Sys.getenv("GEMCOX_REPS", CFG$reps_full))
RUN  <- strsplit(Sys.getenv("GEMCOX_RUN", "C1,C2,C3a,C3b,C4,C5"), ",")[[1]]
TIME_ONLY <- Sys.getenv("GEMCOX_TIME", "0") == "1"
ARMS <- match.arg(Sys.getenv("GEMCOX_ARMS", "primary"), c("primary", "convergence", "matched"))
if (ARMS == "matched") RUN <- intersect(RUN, names(MATCHED_CELLS))

## Estimators (m_lcc, m_gemcox_tight, ...) and row helpers are in
## C_methods.R, shared with 07_convergence_rerun.R.
source("C_methods.R")

c_methods <- function(E, cell) {
  if (ARMS == "convergence") return(TIGHT_METHODS)
  if (ARMS == "matched") return(MATCHED_METHODS)
  m <- LCC_METHODS[LCC]
  if (cell %in% E$sens_cells) m <- c(m, LCC_METHODS[LCC_SENS])
  if (E$fit == "all") m <- c(METHODS, m)
  m
}
## ---- timing: one dataset per cell, serial ----------------------------------
## For C1-C3 the frozen GeM-Cox (gamma = 1) fit is also refitted on the timed
## dataset and its recovery compared with the stored value.
if (TIME_ONLY) {
  cat("\n[timing] one dataset (replicate 1) per cell, serial\n")
  TT <- list()
  for (key in RUN) {
    E <- C_EXPS[[key]]
    fz <- if (E$fit == "competitor") read_results(E$frozen) else NULL
    for (j in seq_len(nrow(E$cells))) {
      design <- c_design(E, j); cid <- E$cells$cell[j]
      meth <- c_methods(E, cid)
      if (!is.null(fz)) meth <- c(meth, METHODS["GeM-Cox (gamma=1)"])
      t1 <- proc.time()[["elapsed"]]
      x <- one_rep(1, E$id, cid, design, meth, CFG, hash_extra, seed_exp = E$seed_exp)
      wall <- proc.time()[["elapsed"]] - t1
      same <- NA
      if (!is.null(fz)) {
        ref <- fz$recovery[fz$cell == cid & fz$rep == 1 & fz$method == "GeM-Cox (gamma=1)"]
        same <- identical(x$recovery[x$method == "GeM-Cox (gamma=1)"], ref)
        wall <- wall - x$secs[x$method == "GeM-Cox (gamma=1)"]   # not part of the run
      }
      TT[[length(TT) + 1]] <- data.frame(
        experiment = key, cell = cid, n = design$n, p = design$p, mu_sep = design$mu_sep,
        beta_sep = design$beta_sep, methods = length(c_methods(E, cid)), secs = wall,
        lcc_secs = x$secs[x$method == LCC], lcc_iters = x$iterations[x$method == LCC],
        failed = sum(x$failed), gamma1_refit_identical = same)
      cat(sprintf("  %-4s cell %2d n=%4d p=%2d mu=%.1f beta=%.1f: %5.1f s (competitor %5.1f s, %3d iter)%s\n",
                  key, cid, design$n, design$p, design$mu_sep, design$beta_sep, wall,
                  x$secs[x$method == LCC], x$iterations[x$method == LCC],
                  if (is.na(same)) "" else sprintf("  gamma=1 refit identical: %s", same)))
    }
  }
  TT <- do.call(rbind, TT)
  serial_h <- sum(TT$secs) * REPS / 3600
  cat(sprintf("\n  serial: %.1f s per replicate set; x %d replicates = %.1f CPU hours\n",
              sum(TT$secs), REPS, serial_h))
  print(aggregate(secs ~ experiment, TT, function(s) round(sum(s) * REPS / 3600, 2)))
  saveRDS(TT, file.path(tempdir(), "c_timing.rds"))
  quit(save = "no")
}

## ---- the runs -------------------------------------------------------------
for (key in RUN) {
  E <- C_EXPS[[key]]
  cat(sprintf("\n[%s, %s arms] exp %d, seeds from exp %d, %d cells x %d reps, fit = %s (%d workers)\n",
              key, ARMS, E$id, E$seed_exp, nrow(E$cells), REPS, E$fit, future::nbrOfWorkers()))
  t0 <- Sys.time()
  run_j <- if (ARMS == "matched") which(E$cells$cell %in% MATCHED_CELLS[[key]]) else
    seq_len(nrow(E$cells))
  res <- bind_rows(lapply(run_j, function(j) {
    design <- c_design(E, j); cid <- E$cells$cell[j]
    t1 <- Sys.time()
    x <- run_cell(E$id, cid, design, REPS, methods = c_methods(E, cid), extra = hash_extra,
                  seed_exp = E$seed_exp)
    x <- label_rows(x, key, design)
    lc <- x[x$method == switch(ARMS, primary = LCC, convergence = LCC_TIGHT,
                                matched = MATCHED_LCC$method[1]), ]
    cat(sprintf("  cell %2d n=%4d p=%2d mu=%.1f beta=%.1f: %d rows, %d failed; competitor recovery %.3f, degenerate %.3f (%.1f min)\n",
                cid, design$n, design$p, design$mu_sep, design$beta_sep, nrow(x), sum(x$failed),
                mean(lc$recovery, na.rm = TRUE), mean(lc$degenerate, na.rm = TRUE),
                as.numeric(difftime(Sys.time(), t1, units = "mins"))))
    x
  }))
  out <- switch(ARMS, primary = E$file, convergence = conv_file(E), matched = matched_file(E))
  save_results(res, out, t0)
  cat(sprintf("  saved %s (%.1f min)\n", out, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat("\nCompetitor runs complete. Next: 06_competitor_summarise.R\n")
