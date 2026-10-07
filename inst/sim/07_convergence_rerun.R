###############################################################
## 07_convergence_rerun.R  --  the core frozen experiments at convergence
## (added after sim-freeze-v2; pre-registered in README.md before it was run)
##
## Purpose: establish which frozen claims hold for converged estimators,
## especially GeM-Cox gamma = 1 versus gamma = 0.
##
##  R1  E1, E2, E3a, E3b, E4a: GeM-Cox gamma = 0, GeM-Cox gamma = 1 and the
##      competitor, all to EM tolerance 1e-8 with at most 1000 iterations.
##      - E1-E3b: gamma = 0 is fitted here. gamma = 1 and the competitor at
##        this tolerance were already fitted on the same datasets with the
##        same code in the convergence arm (C*_convergence.rds) and are
##        reused; the timing run refits them on one dataset per cell and
##        checks that they reproduce the stored rows exactly.
##      - E4a: all three are fitted here.
##      Single Cox, two-stage and the oracle are the frozen rows, unchanged.
##  R2  E5c at convergence: the in-sample LRT with the bootstrap null
##      (B = 19), every fit to tolerance 1e-8, on E5c's datasets; type I
##      and power. Fits that stop at 1000 iterations are counted.
##
## Every row records the md5 of its data; 08_convergence_summarise.R checks
## them against C0_dataset_hashes.rds (E1-E3b) and R0_dataset_hashes.rds
## (E4a, E5c) before any comparison with the frozen rows.
##
##   GEMCOX_MANIFEST_SET=R Rscript 05_hash_manifest.R   # E4a/E5c identity, first
##   Rscript 07_convergence_rerun.R                    # GEMCOX_RUN="E1,E4a,E5c" to select
##   GEMCOX_TIME=1 Rscript 07_convergence_rerun.R      # time one dataset per cell of
##                                                     # everything; saves nothing
## ASCII only.
###############################################################

source("00_core.R")
source("C_common.R")
source("C_methods.R")
REPS <- as.integer(Sys.getenv("GEMCOX_REPS", CFG$reps_full))
RUN  <- strsplit(Sys.getenv("GEMCOX_RUN", "E1,E2,E3a,E3b,E4a,E5c"), ",")[[1]]
TIME_ONLY <- Sys.getenv("GEMCOX_TIME", "0") == "1"
E5_SCALE <- REPS / CFG$reps_full
B <- 19L                                   # as in 02_experiments.R

## ---- R1 designs ----------------------------------------------------------------
## key: the C experiment whose datasets these are (NA for E4a); reuse: the
## convergence-arm file with gamma = 1 and the competitor at 1e-8.
R1 <- list(
  E1  = list(id = 201, key = "C1",  file = "R1_E1_convergence.rds"),
  E2  = list(id = 202, key = "C2",  file = "R1_E2_convergence.rds"),
  E3a = list(id = 203, key = "C3a", file = "R1_E3a_convergence.rds"),
  E3b = list(id = 204, key = "C3b", file = "R1_E3b_convergence.rds"),
  E4a = list(id = 205, key = NA,    file = "R1_E4a_convergence.rds"))
r1_design_src <- function(nm) if (is.na(R1[[nm]]$key)) R_EXPS$E4a else C_EXPS[[R1[[nm]]$key]]
r1_methods <- function(nm) if (is.na(R1[[nm]]$key)) c(GEM0_TIGHT_METHOD, TIGHT_METHODS) else
  GEM0_TIGHT_METHOD

## ---- R2: E5c at convergence ------------------------------------------------------
e5c_tight_one <- function(r, cell_id, design, cfg) {
  pin_rng()
  s <- mk_seed(5, cell_id, r)
  d <- simulate(design, s)
  n_warn <- 0L
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(withCallingHandlers(
    gemcox_heterogeneity_test(X_gmm = d$X, time = d$time, status = d$status, B = B,
                              statistic = "lrt", null = "bootstrap", seed = s,
                              gamma = cfg$gamma, lambda = cfg$lambda, alpha = cfg$alpha,
                              normalize_gmm_by_dim = cfg$normalize,
                              tol = TIGHT_TOL, max_iter = TIGHT_MAXIT),
    warning = function(w) { n_warn <<- n_warn + 1L; invokeRestart("muffleWarning") }),
    error = function(e) e)
  base <- data.frame(exp = 207, cell = cell_id, rep = r, seed = s, role = design$role,
                     n = design$n, p = design$p, mu_sep = design$mu_sep,
                     beta_sep = design$beta_sep, rho = design$rho, K_true = design$K_true,
                     events = sum(d$status), secs = proc.time()[["elapsed"]] - t0,
                     n_warnings = n_warn, hash_train = hash_train(d), stringsAsFactors = FALSE)
  if (inherits(res, "error")) {
    return(cbind(base, failed = TRUE, error = conditionMessage(res), stat = NA_real_,
                 p_bootstrap = NA_real_, valid_bootstrap = NA_integer_, n_fits = NA_integer_,
                 n_not_converged = NA_integer_))
  }
  cbind(base, failed = FALSE, error = NA_character_, stat = unname(res$statistic),
        p_bootstrap = res$p.value, valid_bootstrap = res$n_valid, n_fits = res$n_fits,
        n_not_converged = res$n_not_converged)
}
e5c_cells <- R_EXPS$E5c$cells
e5c_cells$reps <- pmax(2L, round(e5c_cells$reps * E5_SCALE))

## ---- timing: one dataset per cell of everything, serial -------------------------
## Also checks that refitting gamma = 1 and the competitor at 1e-8 on
## replicate 1 reproduces the convergence-arm rows that R1 reuses.
if (TIME_ONLY) {
  source("09_tolerance_study.R")     # defines T_DESIGN and t_one(); runs nothing when sourced
  TT <- list()
  add <- function(part, key, cid, secs, reps, extra = "") {
    TT[[length(TT) + 1]] <<- data.frame(part = part, key = key, cell = cid, secs = secs,
                                        reps = reps, stringsAsFactors = FALSE)
    cat(sprintf("  %-9s %-4s cell %2d: %6.1f s x %3d reps%s\n", part, key, cid, secs, reps, extra))
  }
  cat("\n[timing] one dataset (replicate 1) per cell, serial\n")
  for (nm in names(R1)) {
    E <- r1_design_src(nm)
    stored <- if (!is.na(R1[[nm]]$key)) read_results(conv_file(E)) else NULL
    for (j in seq_len(nrow(E$cells))) {
      design <- c_design(E, j); cid <- E$cells$cell[j]
      x <- one_rep(1, R1[[nm]]$id, cid, design, r1_methods(nm), CFG, hash_extra,
                   seed_exp = E$seed_exp)
      chk <- ""
      if (!is.null(stored)) {
        y <- one_rep(1, E$id, cid, design, TIGHT_METHODS, CFG, hash_extra, seed_exp = E$seed_exp)
        s1 <- stored[stored$cell == cid & stored$rep == 1, ]
        same <- all(vapply(c(LCC_TIGHT, GEM_TIGHT), function(m)
          identical(y$recovery[y$method == m], s1$recovery[s1$method == m]) &&
            identical(y$iterations[y$method == m], s1$iterations[s1$method == m]), NA))
        chk <- sprintf("  reused tight rows reproduce: %s", same)
      }
      add("R1", nm, cid, sum(x$secs), REPS, chk)
    }
  }
  for (j in seq_len(nrow(e5c_cells))) {
    t1 <- proc.time()[["elapsed"]]
    x <- e5c_tight_one(1, e5c_cells$cell[j], as.list(e5c_cells[j, ]), CFG)
    add("R2 (E5c)", "E5c", e5c_cells$cell[j], proc.time()[["elapsed"]] - t1, e5c_cells$reps[j],
        sprintf("  (%d fits, %d at the cap)", x$n_fits, x$n_not_converged))
  }
  for (key in names(MATCHED_CELLS)) {
    E <- C_EXPS[[key]]
    for (j in which(E$cells$cell %in% MATCHED_CELLS[[key]])) {
      x <- one_rep(1, E$id, E$cells$cell[j], c_design(E, j), MATCHED_METHODS, CFG, hash_extra,
                   seed_exp = E$seed_exp)
      add("matched", key, E$cells$cell[j], sum(x$secs), REPS)
    }
  }
  for (i in seq_len(nrow(T_DESIGN))) {
    t1 <- proc.time()[["elapsed"]]
    invisible(t_one(1, T_DESIGN[i, ]))
    add("tolerance", T_DESIGN$key[i], T_DESIGN$cell[i], proc.time()[["elapsed"]] - t1,
        T_DESIGN$reps[i])
  }
  TT <- do.call(rbind, TT)
  TT$cpu_h <- TT$secs * TT$reps / 3600
  S <- aggregate(cpu_h ~ part, TT, sum)
  ## wall = serial CPU hours x factor. The convergence arm's actual wall time
  ## was 0.34 x its serial projection on 7 workers; 0.5 is the upper figure
  ## used before (3-4x slowdown under load).
  S$wall_h_central <- S$cpu_h * 0.34
  S$wall_h_upper   <- S$cpu_h * 0.50
  cat("\n  projection (hours):\n")
  print(transform(S, cpu_h = round(cpu_h, 1), wall_h_central = round(wall_h_central, 1),
                  wall_h_upper = round(wall_h_upper, 1)), row.names = FALSE)
  cat(sprintf("  total: %.1f CPU hours; wall %.1f h (central) to %.1f h (upper)\n",
              sum(S$cpu_h), sum(S$wall_h_central), sum(S$wall_h_upper)))
  saveRDS(TT, file.path(tempdir(), "r_timing.rds"))
  quit(save = "no")
}

## ---- R1 runs ---------------------------------------------------------------------
for (nm in intersect(RUN, names(R1))) {
  E <- r1_design_src(nm)
  meth <- r1_methods(nm)
  cat(sprintf("\n[R1 %s] exp %d, seeds from exp %d, %d cells x %d reps, methods: %s (%d workers)\n",
              nm, R1[[nm]]$id, E$seed_exp, nrow(E$cells), REPS, paste(names(meth), collapse = "; "),
              future::nbrOfWorkers()))
  t0 <- Sys.time()
  res <- bind_rows(lapply(seq_len(nrow(E$cells)), function(j) {
    design <- c_design(E, j); cid <- E$cells$cell[j]
    t1 <- Sys.time()
    x <- run_cell(R1[[nm]]$id, cid, design, REPS, methods = meth, extra = hash_extra,
                  seed_exp = E$seed_exp)
    x <- label_rows(x, nm, design)
    g0 <- x[x$method == GEM0_TIGHT, ]
    cat(sprintf("  cell %d n=%4d p=%2d mu=%.1f beta=%.1f: %d rows, %d failed, %d at the cap (%.1f min)\n",
                cid, design$n, design$p, design$mu_sep, design$beta_sep, nrow(x), sum(x$failed),
                sum(!x$failed & !x$converged, na.rm = TRUE),
                as.numeric(difftime(Sys.time(), t1, units = "mins"))))
    x
  }))
  save_results(res, R1[[nm]]$file, t0)
  cat(sprintf("  saved %s (%.1f min)\n", R1[[nm]]$file,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

## ---- R2 run ------------------------------------------------------------------------
if ("E5c" %in% RUN) {
  cat(sprintf("\n[R2 E5c] LRT, bootstrap null, B = %d, every fit to tolerance 1e-8 (%d cells, %d workers)\n",
              B, nrow(e5c_cells), future::nbrOfWorkers()))
  t0 <- Sys.time()
  res <- bind_rows(lapply(seq_len(nrow(e5c_cells)), function(j) {
    design <- as.list(e5c_cells[j, ])
    t1 <- Sys.time()
    x <- bind_rows(par_lapply(seq_len(design$reps), e5c_tight_one, cell_id = design$cell,
                              design = design, cfg = CFG))
    cat(sprintf("  cell %2d %-6s n=%d rho=%g mu=%g beta=%g: %d datasets, %d failed, %d of %d fits at the cap (%.1f min)\n",
                design$cell, design$role, design$n, design$rho, design$mu_sep, design$beta_sep,
                nrow(x), sum(x$failed), sum(x$n_not_converged, na.rm = TRUE),
                sum(x$n_fits, na.rm = TRUE), as.numeric(difftime(Sys.time(), t1, units = "mins"))))
    x
  }))
  save_results(res, R_EXPS$E5c$file, t0)
  cat(sprintf("  saved %s (%.1f min)\n", R_EXPS$E5c$file,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat("\nConvergence rerun complete. Next: 08_convergence_summarise.R\n")
