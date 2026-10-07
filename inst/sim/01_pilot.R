###############################################################
## 01_pilot.R  --  run this BEFORE committing compute
##
## Checks that the unified DGP reproduces the qualitative findings the
## paper rests on. Four checks, each with its expectation recorded in
## advance (unchanged from the pre-package pilot):
##
##  P1  IMPOSSIBILITY. With mu_sep = 0 the covariates carry no information
##      about cluster membership, so any covariate-based estimator must
##      sit at the random-direction floor. Expect: single Cox, two-stage
##      and GeM-Cox(gamma=0) at the floor; GeM-Cox(gamma=1) above it.
##      NOTE: the single Cox model returns the same coefficients for both
##      "clusters", which the metric scores exactly at the floor, so that
##      part of P1 holds by construction and is not evidence.
##
##  P2  QUADRATIC SCALING. Expect gain / beta_sep^2 roughly constant.
##
##  P3  TWO SCALES. Expect ARI at chance while recovery exceeds the floor.
##
##  P4  ENGINE. Expect gemcox_selftest() to pass, including K = 1
##      reproducing coxph and the K = 2 E-step survival term.
##
## Run from gemcox/inst/sim:  Rscript 01_pilot.R
###############################################################

source("00_core.R")

cat("\n================ P4: engine self-test ================\n")
P4 <- gemcox_selftest(verbose = TRUE)
if (!isTRUE(P4)) stop("Engine self-test failed. Fix before running any experiment.")

REPS <- CFG$reps_pilot
cells <- expand.grid(n = 400, p = 10, mu_sep = 0,
                     beta_sep = c(0.5, 1, 2, 3),
                     KEEP.OUT.ATTRS = FALSE)
pilot <- run_experiment(0, "pilot: mechanism-only regime", cells, REPS, "pilot.rds")

cat(sprintf("\nfailed method fits: %d of %d\n", sum(pilot$failed), nrow(pilot)))
S <- pilot %>% filter(!failed) %>% group_by(beta_sep, method) %>%
  summarise(reps = n(),
            rec_mcse = mcse(recovery),        # dispersion before the mean
            rec = mean(recovery, na.rm = TRUE),
            ari = mean(ari, na.rm = TRUE), .groups = "drop") %>%
  group_by(beta_sep) %>%
  mutate(orc = rec[method == "Oracle (true labels)"],
         pct = pct_of_oracle(rec, orc, 10)) %>% ungroup()

FL <- random_floor(10)
cat(sprintf("\nrandom-direction floor at p = 10: %.3f\n\n", FL))
print(as.data.frame(S %>%
  transmute(beta_sep, method, reps, recovery = round(rec, 3),
            mcse = round(rec_mcse, 3), pct_of_oracle = round(pct, 1),
            ARI = round(ari, 3))), row.names = FALSE)

## ---- P1 ---------------------------------------------------------------
cat("\n================ P1: impossibility ================\n")
p1 <- S %>% filter(beta_sep >= 1) %>% group_by(method) %>%
  summarise(rec = mean(rec), .groups = "drop") %>%
  mutate(above_floor = round(rec - FL, 3))
print(as.data.frame(p1 %>% mutate(rec = round(rec, 3))), row.names = FALSE)
blind <- p1 %>% filter(method %in% c("Single Cox model", "Two-stage (GMM+Cox)",
                                     "GeM-Cox (gamma=0)"))
g1 <- p1 %>% filter(method == "GeM-Cox (gamma=1)")
cat(sprintf("\n  covariate-based methods: max excess over floor = %+.3f (expect ~0)\n",
            max(blind$above_floor)))
cat(sprintf("  GeM-Cox (gamma=1)     : excess over floor = %+.3f (expect > 0)\n",
            g1$above_floor))
P1 <- max(blind$above_floor) < 0.05 && g1$above_floor > 0.05
cat("  VERDICT: ", if (P1) "P1 reproduced.\n" else
      "P1 NOT reproduced under the unified DGP. Investigate before writing.\n", sep = "")

## ---- P2 ---------------------------------------------------------------
cat("\n================ P2: quadratic scaling ================\n")
W <- S %>% filter(method %in% c("GeM-Cox (gamma=0)", "GeM-Cox (gamma=1)")) %>%
  select(beta_sep, method, rec) %>%
  pivot_wider(names_from = method, values_from = rec) %>%
  mutate(gain = `GeM-Cox (gamma=1)` - `GeM-Cox (gamma=0)`,
         scaled = gain / beta_sep^2)
print(as.data.frame(W %>% mutate(across(where(is.numeric), ~round(.x, 4)))),
      row.names = FALSE)
cv <- sd(W$scaled[W$beta_sep >= 1]) / mean(W$scaled[W$beta_sep >= 1])
cat(sprintf("\n  coefficient of variation of gain/beta_sep^2 (beta_sep >= 1): %.2f\n", cv))
P2 <- is.finite(cv) && cv < 0.5
cat("  VERDICT: ", if (P2) "P2 reproduced (roughly constant).\n" else
    "scaling not clean at pilot size; recheck at full replication.\n", sep = "")

## ---- P3 ---------------------------------------------------------------
cat("\n================ P3: two scales ================\n")
p3 <- S %>% filter(method == "GeM-Cox (gamma=1)") %>%
  transmute(beta_sep, ARI = round(ari, 3), recovery = round(rec, 3),
            excess_over_floor = round(rec - FL, 3))
print(as.data.frame(p3), row.names = FALSE)
P3 <- max(p3$ARI, na.rm = TRUE) < 0.10 && max(p3$excess_over_floor) > 0.05
cat("  VERDICT: ", if (P3)
      "P3 reproduced: labels at chance while the contrast is recovered.\n" else
      "P3 pattern unclear at pilot size.\n", sep = "")

cat("\n================ summary ================\n")
print(data.frame(check = c("P1", "P2", "P3", "P4"),
                 reproduced = c(P1, P2, P3, isTRUE(P4))), row.names = FALSE)
cat("\n##### If all four verdicts are positive, proceed to 02_experiments.R.\n")
cat("##### If any is negative, the unified DGP has changed a headline\n")
cat("##### result and that must be resolved before drafting.\n")
