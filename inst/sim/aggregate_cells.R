###############################################################
## aggregate_cells.R  --  assemble results/<id>/cell_<k>.rds (run_cell.R)
## into the per-experiment results file the summaries read
## (added after sim-freeze-v2; tasks/HPC.md, step 3)
##
##   Rscript aggregate_cells.R --exp <id>
##
## It refuses to write if any registered cell is missing or incomplete.
## Rows are bound in cell order, which is the order the original scripts
## produce. The meta keeps every cell's provenance, plus the set of commits
## and whether any cell ran with uncommitted changes.
## ASCII only.
###############################################################

args <- commandArgs(trailingOnly = TRUE)
EXP <- if ("--exp" %in% args) args[match("--exp", args) + 1] else stop("usage: --exp <id>")
Sys.setenv(GEMCOX_WORKERS = "1")
## the registry lives in run_cell.R; load it without running a cell
src <- parse("run_cell.R", keep.source = FALSE)
Sys.setenv(GEMCOX_RUN = "none", GEMCOX_TIME = "0")
suppressMessages(source("07_convergence_rerun.R"))
for (e in src) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% c("matched_entry", "REGISTRY")) eval(e, globalenv())
}
if (!EXP %in% names(REGISTRY)) stop("unknown --exp ", EXP)
X <- REGISTRY[[EXP]]
DIR <- file.path(RESULTS, EXP)
objs <- lapply(X$cells, function(cid) {
  f <- file.path(DIR, sprintf("cell_%d.rds", cid))
  if (!file.exists(f)) return(NULL)
  o <- readRDS(f)
  if (!identical(as.integer(o$reps), as.integer(X$reps(cid)))) return(NULL)
  o
})
missing <- X$cells[vapply(objs, is.null, NA)]
if (length(missing)) stop(EXP, ": missing or incomplete cells: ", paste(missing, collapse = ", "))
rows <- bind_rows(lapply(objs, `[[`, "rows"))
metas <- lapply(objs, `[[`, "meta")
commits <- unique(vapply(metas, function(m) as.character(m$git_commit), ""))
dirty <- any(vapply(metas, function(m) isTRUE(m$git_dirty), NA))
meta <- list(created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
             assembled_by = "aggregate_cells.R", exp = EXP, cells = X$cells,
             gemcox_version = metas[[1]]$gemcox_version,
             git_commit = if (length(commits) == 1) commits else commits, git_dirty = dirty,
             config = metas[[1]]$config,
             wall_minutes = sum(vapply(metas, function(m) m$wall_minutes, 0)),
             cell_meta = metas)
if (length(commits) > 1) warning("cells were run at different commits: ", paste(commits, collapse = ", "))
if (dirty) warning("at least one cell ran with uncommitted tracked changes")
saveRDS(list(rows = rows, meta = meta), file.path(RESULTS, X$file))
cat(sprintf("%s: %d cells, %d rows, %d failed -> %s (commits: %s; dirty: %s)\n", EXP, length(objs),
            nrow(rows), sum(rows$failed), file.path(RESULTS, X$file),
            paste(substr(commits, 1, 7), collapse = ","), dirty))
