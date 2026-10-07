###############################################################
## run_cell.R  --  run ONE cell of a registered experiment, checkpointed
## (added after sim-freeze-v2 for cluster runs; tasks/HPC.md, step 3)
##
##   Rscript run_cell.R --exp <id> --cell <k>     run the cell (skips if complete)
##   Rscript run_cell.R --exp <id> --cell <k> --time-only
##                                                time replicate 1, print a suggested
##                                                wall-time limit; writes nothing
##   Rscript run_cell.R --exp <id> --list         list the cells and their status
##
## Output: results/<id>/cell_<k>.rds, written atomically (temporary file in
## the same directory, then rename). It holds list(rows, meta, exp, cell,
## reps), with the standard provenance. A cell whose file already exists and
## has the requested number of replicates is skipped, so resubmission is safe.
## aggregate_cells.R assembles the cell files into the per-experiment object
## (list(rows, meta)) that 06_ and 08_ summaries read.
##
## The per-cell computation is the existing harness's, called unchanged:
##   R2          e5c_tight_one() from 07_convergence_rerun.R (E5c at convergence)
##   matched_C1  00_core.R run_cell() with MATCHED_METHODS (C_methods.R), as in
##   matched_C2  05_competitor.R with GEMCOX_ARMS=matched
##   matched_C5
##   E6          00_core.R run_cell() with E6_METHODS and e6_extra (C_methods.R)
## Seeds are deterministic in (experiment, cell, replicate), so a cell gives
## identical rows whether it is run here or by the original script.
## Workers: GEMCOX_WORKERS (on SLURM, set it to SLURM_CPUS_PER_TASK).
## ASCII only.
###############################################################

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag) { i <- match(flag, args); if (is.na(i)) NA_character_ else args[i + 1] }
EXP <- arg_val("--exp")
CELL <- suppressWarnings(as.integer(arg_val("--cell")))
TIME_ONLY_CELL <- "--time-only" %in% args
LIST <- "--list" %in% args
if (is.na(EXP)) stop("usage: Rscript run_cell.R --exp <id> (--cell <k> [--time-only] | --list)")

Sys.setenv(GEMCOX_RUN = "none", GEMCOX_TIME = "0")   # 07 then defines its functions and runs nothing
if (TIME_ONLY_CELL || LIST) Sys.setenv(GEMCOX_WORKERS = "1")
suppressMessages(source("07_convergence_rerun.R"))

## ---- registry ------------------------------------------------------------------
matched_entry <- function(key) {
  E <- C_EXPS[[key]]
  cells <- E$cells[E$cells$cell %in% MATCHED_CELLS[[key]], , drop = FALSE]
  list(file = matched_file(E), cells = cells$cell, reps = function(cid) REPS,
       design = function(cid) c_design(E, which(E$cells$cell == cid)),
       run = function(cid, reps) {
         design <- c_design(E, which(E$cells$cell == cid))
         x <- run_cell(E$id, cid, design, reps, methods = MATCHED_METHODS, extra = hash_extra,
                       seed_exp = E$seed_exp)
         label_rows(x, key, design)
       },
       one = function(cid) {
         design <- c_design(E, which(E$cells$cell == cid))
         one_rep(1, E$id, cid, design, MATCHED_METHODS, CFG, hash_extra, seed_exp = E$seed_exp)
       })
}
REGISTRY <- list(
  R2 = list(file = R_EXPS$E5c$file, cells = e5c_cells$cell,
            reps = function(cid) e5c_cells$reps[e5c_cells$cell == cid],
            design = function(cid) as.list(e5c_cells[e5c_cells$cell == cid, ]),
            run = function(cid, reps) {
              design <- as.list(e5c_cells[e5c_cells$cell == cid, ])
              bind_rows(par_lapply(seq_len(reps), e5c_tight_one, cell_id = cid, design = design, cfg = CFG))
            },
            one = function(cid) e5c_tight_one(1, cid, as.list(e5c_cells[e5c_cells$cell == cid, ]), CFG)),
  E6 = list(file = E6_FILE, cells = E6_CELLS$cell, reps = function(cid) REPS,
            design = function(cid) as.list(E6_CELLS[E6_CELLS$cell == cid, ]),
            run = function(cid, reps) {
              design <- as.list(E6_CELLS[E6_CELLS$cell == cid, ])
              x <- run_cell(E6_EXP_ID, cid, design, reps, methods = E6_METHODS, extra = e6_extra,
                            seed_exp = E6_EXP_ID)
              x <- label_rows(x, "E6", design)
              x$set <- design$set; x$followup <- design$followup
              x
            },
            one = function(cid) {
              design <- as.list(E6_CELLS[E6_CELLS$cell == cid, ])
              one_rep(1, E6_EXP_ID, cid, design, E6_METHODS, CFG, e6_extra, seed_exp = E6_EXP_ID)
            }),
  matched_C1 = matched_entry("C1"),
  matched_C2 = matched_entry("C2"),
  matched_C5 = matched_entry("C5"))
if (!EXP %in% names(REGISTRY)) stop("unknown --exp ", EXP, "; registered: ", paste(names(REGISTRY), collapse = ", "))
X <- REGISTRY[[EXP]]
DIR <- file.path(RESULTS, EXP)
dir.create(DIR, recursive = TRUE, showWarnings = FALSE)
cell_file <- function(cid) file.path(DIR, sprintf("cell_%d.rds", cid))
cell_done <- function(cid) {
  f <- cell_file(cid)
  if (!file.exists(f)) return(FALSE)
  o <- tryCatch(readRDS(f), error = function(e) NULL)
  !is.null(o) && identical(as.integer(o$reps), as.integer(X$reps(cid)))
}

## ---- list -----------------------------------------------------------------------
if (LIST) {
  for (cid in X$cells) cat(sprintf("%s cell %d: %d replicates, %s\n", EXP, cid, X$reps(cid),
                                   if (cell_done(cid)) "complete" else "missing"))
  quit(save = "no")
}
if (is.na(CELL) || !CELL %in% X$cells) {
  stop("--cell must be one of ", paste(X$cells, collapse = ", "), " for ", EXP)
}

## ---- time one replicate ----------------------------------------------------------
if (TIME_ONLY_CELL) {
  t1 <- proc.time()[["elapsed"]]
  invisible(X$one(CELL))
  secs <- proc.time()[["elapsed"]] - t1
  cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", Sys.getenv("GEMCOX_CPUS", "16")))
  ## suggested limit: replicates x seconds / cpus, x 3 for slowdown under load
  ## and slower datasets, plus 10 minutes of start-up
  mins <- ceiling(X$reps(CELL) * secs / cpus * 3 / 60 + 10)
  cat(sprintf("%s\t%d\t%d\t%.1f\t%d\t%d\n", EXP, CELL, X$reps(CELL), secs, cpus, mins))
  quit(save = "no")
}

## ---- run ------------------------------------------------------------------------------
if (cell_done(CELL)) {
  cat(sprintf("%s cell %d already complete: %s (skipped)\n", EXP, CELL, cell_file(CELL)))
  quit(save = "no")
}
t0 <- Sys.time()
reps <- X$reps(CELL)
cat(sprintf("[%s cell %d] %d replicates, %d workers, started %s\n", EXP, CELL, reps,
            future::nbrOfWorkers(), format(t0)))
rows <- X$run(CELL, reps)
obj <- list(rows = rows, meta = run_meta(t0), exp = EXP, cell = CELL, reps = reps)
tmp <- tempfile(pattern = sprintf(".cell_%d_", CELL), tmpdir = DIR, fileext = ".rds.tmp")
saveRDS(obj, tmp)
stopifnot(identical(readRDS(tmp)$reps, reps))
if (!file.rename(tmp, cell_file(CELL))) stop("could not rename ", tmp)
cat(sprintf("[%s cell %d] saved %s: %d rows, %d failed (%.1f min)\n", EXP, CELL, cell_file(CELL),
            nrow(rows), sum(rows$failed), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
