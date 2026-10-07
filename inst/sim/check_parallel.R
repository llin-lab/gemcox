###############################################################
## check_parallel.R  --  serial and parallel runs must be identical
##
## Runs two 5-replicate cells (all five methods; n = 400 and n = 2400)
## serially and on PSOCK workers, and compares every result column except
## timings. Run from gemcox/inst/sim before any experiment.
###############################################################

source("00_core.R")
## Two cells: n = 400, and n = 2400, where mclust initialises from a random
## subset (n > 2000) and results depend on seeding the stream per method.
designs <- list(list(n = 400, p = 10, mu_sep = 0.5, beta_sep = 2),
                list(n = 2400, p = 10, mu_sep = 0.5, beta_sep = 2))
strip <- function(x) x[, setdiff(names(x), "secs"), drop = FALSE]
run_both <- function() bind_rows(lapply(seq_along(designs), function(j)
  run_cell(99, j, designs[[j]], 5)))

setup_parallel(1)
t0 <- Sys.time(); serial <- run_both(); ts <- Sys.time() - t0
workers_used <- max(2L, CFG$workers)
setup_parallel(workers_used)
omp <- unlist(par_lapply(1:workers_used, function(i) Sys.getenv("OMP_NUM_THREADS")))
t0 <- Sys.time(); parallel <- run_both(); tp <- Sys.time() - t0

same <- identical(strip(serial), strip(parallel))
cat(sprintf("\ntwo 5-replicate cells (n = 400, 2400): %d rows; serial %.1f s, parallel (%d workers) %.1f s\n",
            nrow(serial), as.numeric(ts, units = "secs"), workers_used,
            as.numeric(tp, units = "secs")))
cat("OMP_NUM_THREADS in workers:", paste(unique(omp), collapse = ", "), "\n")
cat("identical results (all columns except timings):", same, "\n")
cat("failed method fits (serial / parallel):", sum(serial$failed), "/", sum(parallel$failed), "\n")
if (any(serial$failed) || any(parallel$failed)) {
  print(unique(rbind(serial, parallel)[serial$failed | parallel$failed, c("method", "error")]))
  stop("Method fits failed in the check cell; identical failures would hide a problem.")
}
if (!same) {
  diffcols <- names(strip(serial))[!mapply(identical, strip(serial), strip(parallel))]
  cat("columns that differ:", paste(diffcols, collapse = ", "), "\n")
  stop("Serial and parallel runs differ.")
}
