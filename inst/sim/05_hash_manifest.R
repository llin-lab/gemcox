###############################################################
## 05_hash_manifest.R  --  dataset identity for the competitor experiments
## (added after sim-freeze-v2; pre-registered in README.md)
##
## C1-C4 reuse frozen datasets, and C1-C3 pair the competitor with frozen
## rows. This script establishes that the datasets are the same:
##  1. For each frozen results file it installs the package AT THE COMMIT
##     THAT PRODUCED THE FILE (recorded in its provenance) into a temporary
##     library and, in a separate R session, sources that commit's 00_core.R
##     and regenerates every training dataset (and, for C1-C3, every test
##     set) from the seeds stored in the frozen rows. Event counts are
##     compared with the stored ones.
##  2. With the current package and 00_core.R it generates the same
##     datasets exactly as 05_competitor.R will (seeds from mk_seed()).
##  3. It saves both sets of md5 hashes, and the comparison, to
##     results/C0_dataset_hashes.rds.
## 06_competitor_summarise.R refuses to compute a paired difference for any
## replicate whose hashes, seed or event count do not match.
##
##   Rscript 05_hash_manifest.R        (from gemcox/inst/sim)
##   GEMCOX_MANIFEST_SET=R Rscript 05_hash_manifest.R
## Set "C" (default): C1-C4, saved as C0_dataset_hashes.rds. Set "R" (added
## for the convergence rerun): E4a and E5c, saved as R0_dataset_hashes.rds.
## ASCII only.
###############################################################

args <- commandArgs(trailingOnly = TRUE)

## ================= CHILD: runs under a frozen commit ===================
if (length(args) && args[1] == "--child") {
  commit <- args[2]; frozen_sim <- args[3]; results_dir <- args[4]; out <- args[5]
  SET <- if (length(args) >= 6) args[6] else "C"
  lib <- Sys.getenv("R_LIBS")
  source(file.path(frozen_sim, "00_core.R"), chdir = TRUE)   # the FROZEN 00_core.R
  source("C_common.R")                                       # hashes and designs only
  EXPS <- if (SET == "C") C_EXPS else R_EXPS
  pkg <- normalizePath(find.package("gemcox"))
  if (!startsWith(pkg, normalizePath(lib))) stop("gemcox not loaded from the frozen library: ", pkg)
  ## The frozen runs drew their data in future workers, which use
  ## L'Ecuyer-CMRG; the commits before pin_rng() relied on that.
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  rows <- list()
  for (key in names(EXPS)) {
    E <- EXPS[[key]]
    if (is.na(E$commit) || E$commit != commit) next
    R <- readRDS(file.path(results_dir, E$frozen))$rows
    for (j in seq_len(nrow(E$cells))) {
      design <- c_design(E, j)
      cid <- E$cells$cell[j]
      F <- unique(R[R$cell == cid, intersect(c("rep", "seed", "n", "p", "mu_sep", "beta_sep", "events"),
                                              names(R))])
      stopifnot(!anyDuplicated(F$rep))
      ## the frozen rows must describe the same cell
      for (v in c("n", "p", "mu_sep", "beta_sep")) stopifnot(all(F[[v]] == design[[v]]))
      for (v in intersect(c("rho", "K_true"), intersect(names(R), names(design)))) {
        stopifnot(all(R[R$cell == cid, v] == design[[v]]))
      }
      for (i in seq_len(nrow(F))) {
        s <- F$seed[i]
        d <- simulate(design, s)
        ht <- NA_character_
        if (frozen_tests(E)) ht <- hash_test(simulate(design, s + 500000L, n = CFG$n_test))
        rows[[length(rows) + 1]] <- data.frame(
          experiment = key, frozen_file = E$frozen, cell = cid, rep = F$rep[i], seed_frozen = s,
          events_stored = F$events[i], events_frozen_code = sum(d$status),
          hash_train_frozen = hash_train(d), hash_test_frozen = ht, stringsAsFactors = FALSE)
      }
    }
  }
  saveRDS(list(rows = do.call(rbind, rows), commit = commit, package_path = pkg,
               gemcox_version = as.character(utils::packageVersion("gemcox")),
               sessionInfo = utils::capture.output(utils::sessionInfo())), out)
  quit(save = "no")
}

## ================= PARENT ===============================================
Sys.setenv(GEMCOX_WORKERS = "1")
source("00_core.R")
source("C_common.R")
SET <- Sys.getenv("GEMCOX_MANIFEST_SET", "C")
EXPS <- if (SET == "C") C_EXPS else R_EXPS
t0 <- Sys.time()
REPS <- as.integer(Sys.getenv("GEMCOX_REPS", CFG$reps_full))
TMP <- Sys.getenv("GEMCOX_MANIFEST_TMP", file.path(tempdir(), "frozen"))
dir.create(TMP, showWarnings = FALSE, recursive = TRUE)
top <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
results_dir <- normalizePath(RESULTS)

## 1. frozen hashes, one R session per producing commit
commits <- unique(na.omit(vapply(EXPS, `[[`, "", "commit")))
frozen <- list(); frozen_info <- list()
for (cm in commits) {
  full <- system2("git", c("-C", shQuote(top), "rev-parse", cm), stdout = TRUE)
  src <- file.path(TMP, cm); lib <- file.path(TMP, paste0(cm, "_lib"))
  unlink(c(src, lib), recursive = TRUE)
  dir.create(src, recursive = TRUE); dir.create(lib)
  tarf <- file.path(TMP, paste0(cm, ".tar"))
  stopifnot(system2("git", c("-C", shQuote(top), "archive", "-o", shQuote(tarf), full, "gemcox")) == 0)
  utils::untar(tarf, exdir = src)
  stopifnot(system2(file.path(R.home("bin"), "R"),
                    c("CMD", "INSTALL", "--no-docs", paste0("--library=", shQuote(lib)),
                      shQuote(file.path(src, "gemcox"))), stdout = FALSE, stderr = FALSE) == 0)
  out <- file.path(TMP, paste0(cm, "_hashes.rds"))
  cat(sprintf("  frozen commit %s: regenerating datasets ...\n", cm))
  st <- system2(file.path(R.home("bin"), "Rscript"),
                c("05_hash_manifest.R", "--child", cm, shQuote(file.path(src, "gemcox", "inst", "sim")),
                  shQuote(results_dir), shQuote(out), SET),
                env = c(paste0("R_LIBS=", lib), "GEMCOX_WORKERS=1",
                        paste0("GEMCOX_RESULTS=", file.path(TMP, paste0(cm, "_results")))),
                stdout = FALSE, stderr = file.path(TMP, paste0(cm, "_child.log")))
  if (st != 0) stop("child session failed for ", cm, "; see ", file.path(TMP, paste0(cm, "_child.log")))
  x <- readRDS(out)
  frozen[[cm]] <- cbind(x$rows, frozen_commit = full, stringsAsFactors = FALSE)
  frozen_info[[cm]] <- x[c("commit", "package_path", "gemcox_version", "sessionInfo")]
}
FZ <- do.call(rbind, frozen)

## 2. current hashes, generated exactly as one_rep() in 00_core.R does
cur <- list()
for (key in names(EXPS)) {
  E <- EXPS[[key]]
  if (is.na(E$frozen)) next
  for (j in seq_len(nrow(E$cells))) {
    design <- c_design(E, j); cid <- E$cells$cell[j]
    reps <- if (!is.null(E$cells$reps)) E$cells$reps[j] else REPS
    for (r in seq_len(reps)) {
      pin_rng()
      s <- mk_seed(E$seed_exp, cid, r)
      d <- simulate(design, s)
      ## test sets only where the rerun uses them (E5c's does not)
      hte <- if (isFALSE(E$test_sets)) NA_character_ else
        hash_test(simulate(design, s + 500000L, n = CFG$n_test))
      cur[[length(cur) + 1]] <- data.frame(
        experiment = key, cell = cid, rep = r, seed_current = s,
        events_current = sum(d$status), hash_train_current = hash_train(d),
        hash_test_current = hte, stringsAsFactors = FALSE)
    }
  }
}
CU <- do.call(rbind, cur)

## 3. compare
M <- merge(CU, FZ, by = c("experiment", "cell", "rep"), all = TRUE)
M$seed_matches     <- M$seed_current == M$seed_frozen
M$events_match     <- M$events_current == M$events_stored & M$events_frozen_code == M$events_stored
M$identical_train  <- M$hash_train_current == M$hash_train_frozen
M$identical_test   <- ifelse(is.na(M$hash_test_frozen), NA, M$hash_test_current == M$hash_test_frozen)
M$verified <- M$seed_matches & M$events_match & M$identical_train &
  (is.na(M$identical_test) | M$identical_test)
M$verified[is.na(M$verified)] <- FALSE

cat("\n=== dataset identity: current code vs the commit that produced each frozen file ===\n")
S <- do.call(rbind, lapply(split(M, M$experiment), function(x) data.frame(
  experiment = x$experiment[1], frozen_file = x$frozen_file[1],
  frozen_commit = substr(x$frozen_commit[1], 1, 7), datasets = nrow(x),
  seeds_match = sum(x$seed_matches, na.rm = TRUE), events_match = sum(x$events_match, na.rm = TRUE),
  train_identical = sum(x$identical_train, na.rm = TRUE),
  test_identical = if (all(is.na(x$identical_test))) NA else sum(x$identical_test, na.rm = TRUE),
  verified = sum(x$verified))))
print(S, row.names = FALSE)
cat(sprintf("\n  ALL %d DATASETS VERIFIED: %s\n", nrow(M), all(M$verified)))

meta <- run_meta(t0)
meta$frozen_sessions <- frozen_info
out_file <- sprintf("%s0_dataset_hashes.rds", SET)
saveRDS(list(rows = M, meta = meta), file.path(RESULTS, out_file))
cat("  saved", out_file, "\n")
