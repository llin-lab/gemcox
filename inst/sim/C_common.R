###############################################################
## C_common.R  --  shared by 05_hash_manifest.R, 05_competitor.R and
## 06_competitor_summarise.R (latent-class Cox competitor, added after
## sim-freeze-v2; pre-registered in README.md).
##
## Self-contained: 05_hash_manifest.R also sources this file in R sessions
## running the FROZEN versions of the package and of 00_core.R.
##
## ASCII only.
###############################################################

## ---- dataset hashes --------------------------------------------------------
## md5 of the serialised data. Serialisation format 2 writes no R version
## into the header, so the hash depends only on the data.
data_hash <- function(x) {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f))
  saveRDS(x, f, version = 2, compress = FALSE)
  unname(tools::md5sum(f))
}
HASH_FIELDS <- c("X", "time", "status", "Z", "eta")
hash_train <- function(d) data_hash(d[HASH_FIELDS])
hash_test  <- function(te) data_hash(te[HASH_FIELDS])

## ---- the experiments -------------------------------------------------------
## id: experiment id recorded in the rows. seed_exp: experiment id the SEEDS
## come from, i.e. the frozen experiment whose datasets are reused (C5 is new
## and uses its own). cell: the frozen cell index, which the seed depends on.
## frozen / commit: the frozen results file and the commit that produced it.
## fit: "competitor" (pair with the frozen rows) or "all" (no frozen recovery
## results exist, so every method is fitted). sens_cells: cells that also
## run the sensitivity arm (gating penalty lambda / 10).
C_EXPS <- list(
  C1 = list(id = 101, seed_exp = 1, file = "C1_mechanism_only.rds",
            frozen = "E1_existence.rds", commit = "de83301",
            cells = data.frame(cell = 1:4, n = 400, p = 10, mu_sep = 0,
                               beta_sep = c(0.5, 1, 2, 3)),
            fit = "competitor", sens_cells = 1:4),
  C2 = list(id = 102, seed_exp = 2, file = "C2_comparison.rds",
            frozen = "E2_comparison.rds", commit = "f965278",
            cells = cbind(cell = 1:12, expand.grid(n = 800, p = 10, mu_sep = c(0, 0.5, 1, 2),
                                                   beta_sep = c(1, 2, 3), KEEP.OUT.ATTRS = FALSE)),
            fit = "competitor", sens_cells = c(6, 7)),      # mu_sep 0.5 and 1, beta_sep 2
  C3a = list(id = 103, seed_exp = 3, file = "C3a_n.rds",
             frozen = "E3a_n.rds", commit = "f965278",
             cells = data.frame(cell = 1:4, n = c(200, 400, 800, 1600), p = 10, mu_sep = 0.5,
                                beta_sep = 2),
             fit = "competitor", sens_cells = integer(0)),
  C3b = list(id = 104, seed_exp = 31, file = "C3b_p.rds",
             frozen = "E3b_p.rds", commit = "f965278",
             cells = data.frame(cell = 1:4, n = 800, p = c(10, 20, 40, 80), mu_sep = 0.5,
                                beta_sep = 2),
             fit = "competitor", sens_cells = integer(0)),
  C4 = list(id = 105, seed_exp = 23, file = "C4_cvia078.rds",
            frozen = "E2b_c_cvia078.rds", commit = "30595fa",
            cells = data.frame(cell = 2:3, n = 117, p = 9, mu_sep = c(0, 0.5), beta_sep = 2,
                               K_true = 2, event_rate = 0.44),
            fit = "all", sens_cells = integer(0)),
  C5 = list(id = 106, seed_exp = 106, file = "C5_lognormal.rds",
            frozen = NA_character_, commit = NA_character_,
            cells = data.frame(cell = 1:2, n = 800, p = 10, mu_sep = c(0.5, 1), beta_sep = 2,
                               features = "lognormal", stringsAsFactors = FALSE),
            fit = "all", sens_cells = integer(0)))

## design list for one cell (the frozen fields; rho = 0 throughout)
c_design <- function(E, j) {
  row <- E$cells[j, , drop = FALSE]
  as.list(row[setdiff(names(row), "cell")])
}

LCC        <- "Latent-class Cox"
LCC_SENS   <- "Latent-class Cox (gating lambda/10)"
## convergence-sensitivity arm (pre-registered; see README): both models
## with a tight EM tolerance, same datasets, results in *_convergence.rds
TIGHT_TOL  <- 1e-8
TIGHT_MAXIT <- 1000L
LCC_TIGHT  <- "Latent-class Cox (tol 1e-8)"
GEM_TIGHT  <- "GeM-Cox (gamma=1, tol 1e-8)"
conv_file  <- function(E) sub("\\.rds$", "_convergence.rds", E$file)
## matched-regularisation comparison (sensitivity analysis; see README)
MATCHED_CELLS <- list(C1 = 1:4, C2 = c(6, 7), C5 = 1:2)
matched_file <- function(E) sub("\\.rds$", "_matched.rds", E$file)

## ---- convergence rerun of frozen experiments (R-series; see README) ---------
## E1, E2, E3a and E3b are C1, C2, C3a and C3b's datasets (verified in
## C0_dataset_hashes.rds). E4a and E5c are verified in R0_dataset_hashes.rds.
## E5c's cells are E5a's (02_experiments.R, e5a_cells), in the same order.
E5C_CELLS <- rbind(
  data.frame(cell = 1:2, role = "type I", n = 400, p = 10, mu_sep = 0, beta_sep = 0,
             rho = c(0, 0.5), K_true = 1, reps = c(500, 200), stringsAsFactors = FALSE),
  cbind(cell = 3:10, role = "power",
        expand.grid(n = c(400, 800), p = 10, mu_sep = c(0, 0.5), beta_sep = c(1, 2),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE),
        rho = 0, K_true = 2, reps = 200, stringsAsFactors = FALSE))
R_EXPS <- list(
  E4a = list(id = 205, seed_exp = 4, file = "R1_E4a_convergence.rds",
             frozen = "E4a_two_scales.rds", commit = "de83301",
             cells = data.frame(cell = 1:4, n = c(800, 1600, 3200, 6400), p = 10, mu_sep = 0,
                                beta_sep = 2),
             test_sets = TRUE),
  E5c = list(id = 207, seed_exp = 5, file = "R2_E5c_convergence.rds",
             frozen = "E5c_lrt.rds", commit = "f965278", cells = E5C_CELLS,
             test_sets = FALSE))
## Whether the frozen experiment generated test sets (C: those it pairs with)
frozen_tests <- function(E) if (!is.null(E$test_sets)) E$test_sets else identical(E$fit, "competitor")

## ---- E6: events vs sample size (added after sim-freeze-v2; see README) -------
## Administrative censoring at a fixed follow-up. The follow-up for each
## (p, mu_sep, event rate) is the event-time quantile at that rate, from one
## reference sample of 200,000 (gemcox_simulate(censoring = "administrative",
## followup = Inf, seed = 20261006), L'Ecuyer-CMRG, beta_sep = 2). It is
## fixed before the run and does not depend on n.
E6_FOLLOWUP <- data.frame(
  p = c(rep(10, 12), 20, 40, 20, 40, 20, 40, 20, 40),
  mu_sep = c(rep(c(0, 1), 6), 0, 0, 1, 1, 0, 0, 1, 1),
  event_rate = c(rep(c(0.1, 0.175, 0.2, 0.275, 0.35, 0.55), each = 2), rep(0.2, 4), rep(0.55, 4)),
  followup = c(29.199, 28.4579, 46.7501, 45.2791, 52.4465, 51.336, 71.7102, 70.154, 92.9026,
               91.8562, 168.71, 168.234, 52.132, 52.5913, 50.949, 51.444, 168.711, 168.886,
               168.318, 168.731))
E6_CELLS <- local({
  fac <- expand.grid(n = c(200, 400, 800, 1600), event_rate = c(0.1, 0.2, 0.35, 0.55),
                     mu_sep = c(0, 1), p = 10, KEEP.OUT.ATTRS = FALSE)
  fac$set <- "factorial"
  mat <- expand.grid(n = c(400, 800, 1600), event_rate = c(0.175, 0.275), mu_sep = c(0, 1),
                     p = 10, KEEP.OUT.ATTRS = FALSE)
  mat$set <- "matched events"      # pairs with (n / 2, 2 x rate): rates 0.35 and 0.55
  pax <- expand.grid(p = c(20, 40), n = c(400, 1600), event_rate = c(0.2, 0.55), mu_sep = c(0, 1),
                     KEEP.OUT.ATTRS = FALSE)
  pax$set <- "p axis"
  x <- rbind(fac[, c("n", "p", "mu_sep", "event_rate", "set")],
             mat[, c("n", "p", "mu_sep", "event_rate", "set")],
             pax[, c("n", "p", "mu_sep", "event_rate", "set")])
  x <- merge(x, E6_FOLLOWUP, by = c("p", "mu_sep", "event_rate"), sort = FALSE)
  x <- x[order(match(x$set, c("factorial", "matched events", "p axis")), x$p, x$mu_sep,
               x$event_rate, x$n), ]
  stopifnot(!anyNA(x$followup), nrow(x) == 60)
  data.frame(cell = seq_len(nrow(x)), n = x$n, p = x$p, mu_sep = x$mu_sep, beta_sep = 2,
             event_rate = x$event_rate, followup = x$followup, censoring = "administrative",
             set = x$set, stringsAsFactors = FALSE)
})
E6_EXP_ID <- 27
E6_FILE <- "E6_events.rds"
