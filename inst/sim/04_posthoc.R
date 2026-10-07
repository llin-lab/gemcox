###############################################################
## 04_posthoc.R  --  POST HOC diagnostic (NOT pre-specified)
##
## Added after the pilot, which found the two-stage estimator's fitted
## contrast aligned with the shared coefficient ("base", orthogonal to the
## true contrast by design), placing it below the isotropic random-direction
## floor. Question: does that alignment require mechanism heterogeneity?
##
## Design: mu_sep = 0, n = 400, p = 10, 200 replicates.
##   beta_sep = 0  the diagnostic (no heterogeneity: both clusters = base)
##   beta_sep = 2  a heterogeneous reference, same code and seeds scheme
## Methods: two-stage (GMM + Cox) and GeM-Cox (gamma = 0).
## Recorded: contrast norm ||beta1_hat - beta2_hat|| (the harness column
## contrast_norm) and |cos(beta1_hat - beta2_hat, base)|.
##
## Results go to results/PH1_alignment.rds and are summarised in
## 03_summarise.R under a separate POST HOC heading.
###############################################################

source("00_core.R")
REPS <- as.integer(Sys.getenv("GEMCOX_REPS", CFG$reps_full))

align <- function(d, te, m, cfg) {
  base <- Reduce(`+`, d$beta_list) / length(d$beta_list)   # shared coefficient
  ctr <- m$beta_list[[1]] - m$beta_list[[length(m$beta_list)]]
  nrm <- sqrt(sum(ctr^2))
  data.frame(cos_base = if (nrm > 0) abs(sum(ctr * base)) / (nrm * sqrt(sum(base^2)))
                        else NA_real_)
}

cells <- expand.grid(n = 400, p = 10, mu_sep = 0, beta_sep = c(0, 2),
                     KEEP.OUT.ATTRS = FALSE)
run_experiment(90, "POST HOC: two-stage alignment with the shared coefficient",
               cells, REPS, "PH1_alignment.rds",
               methods = METHODS[c("Two-stage (GMM+Cox)", "GeM-Cox (gamma=0)")],
               extra = align)
