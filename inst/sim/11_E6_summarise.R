###############################################################
## 11_E6_summarise.R  --  E6: does recovery follow events or sample size?
## (added after sim-freeze-v2; the analysis pre-registered in README.md,
## written before E6 was run). Reads results/E6_events.rds only.
##
## Outcomes:
##   mechanism recovery  frozen contrast-recovery metric (|cos| of the
##                       coefficient contrast), every cell
##   profile recovery    |cos| of the between-subgroup mean-difference
##                       direction in X (C_methods.R profile_direction()),
##                       mu_sep = 1 cells only
## For each method x outcome x mu_sep: weighted least squares of
## logit(cell-mean recovery) on log(mean events), log(n) and log(p), with
## weights 1 / Var(logit mean) from the delta method. Primary: degenerate
## fits excluded (the frozen convention). Sensitivity: degenerate fits
## scored at the random-direction floor; and a replicate-level fit.
## Classification rule (README): which of log events and log n "carries"
## the effect. Matched-events pairs: 2n at rate r/2 minus n at rate r.
## ASCII only.
###############################################################

suppressMessages(library(ggplot2))
Sys.setenv(GEMCOX_WORKERS = "1")
source("00_core.R"); source("C_common.R"); source("C_methods.R")
OUTP <- function(f) file.path(RESULTS, f)
wcsv <- function(x, f) write.csv(x, OUTP(f), row.names = FALSE)
X <- read_results(E6_FILE)
if (is.null(X)) stop("E6 results not found")

## ---- verification: seeds, and one dataset per replicate across methods ---------
seed_ok <- all(X$seed == mk_seed(E6_EXP_ID, X$cell, X$rep))
same_data <- X %>% filter(!failed) %>% group_by(cell, rep) %>%
  summarise(k = n_distinct(hash_train), .groups = "drop")
cat(sprintf("seeds match mk_seed(27, cell, rep): %s; replicates where methods saw different data: %d\n",
            seed_ok, sum(same_data$k > 1)))
stopifnot(seed_ok, all(same_data$k == 1))

X <- X %>% left_join(E6_CELLS %>% select(cell, set), by = "cell", suffix = c("", ".design")) %>%
  mutate(method = factor(method, c(GEM_TIGHT, LCC_TIGHT, "Two-stage (GMM+Cox)", "Oracle (true labels)")),
         degenerate = !is.na(degenerate) & degenerate,
         capped = !failed & !is.na(converged) & !converged)

## ---- per-cell table: events, events per subgroup and per coefficient, recovery ---
cellsum <- function(R, mode) {
  R <- R %>% filter(!failed)
  if (mode == "exclude") R <- R %>% filter(!degenerate)
  if (mode == "floor") R <- R %>% mutate(recovery = ifelse(degenerate, random_floor(p), recovery),
                                         profile_recovery = ifelse(degenerate & mu_sep > 0,
                                                                   random_floor(p), profile_recovery))
  R %>% group_by(cell, set, n, p, mu_sep, event_rate, method) %>%
    summarise(fits = n(), rec_mcse = mcse(recovery), recovery = mean(recovery),
              prof_mcse = mcse(profile_recovery), profile = mean(profile_recovery),
              .groups = "drop")
}
D <- X %>% group_by(cell, set, n, p, mu_sep, event_rate, followup, method) %>%
  summarise(datasets = n(), failed = sum(failed), degenerate = sum(degenerate & !failed),
            capped = sum(capped), events = mean(events), events_sub_min = mean(events_sub_min),
            events_per_coef = mean(events_per_coef), realized_rate = mean(realized_rate),
            .groups = "drop")
S <- left_join(D, cellsum(X, "exclude"), by = c("cell", "set", "n", "p", "mu_sep", "event_rate", "method")) %>%
  left_join(cellsum(X, "floor") %>% select(cell, method, recovery_floor = recovery, profile_floor = profile),
            by = c("cell", "method"))
wcsv(S, "tab_E6_cells.csv")
cat("\n=== E6 cells: events, events per subgroup (smaller true subgroup), events per coefficient ===\n")
print(as.data.frame(S %>% filter(method == GEM_TIGHT) %>%
  transmute(cell, set, n, p, mu_sep, rate = event_rate, followup = round(followup, 1),
            events = round(events, 1), per_subgroup = round(events_sub_min, 1),
            per_coef = round(events_per_coef, 2))), row.names = FALSE)
cat("\n=== failures, degenerate fits and cap hits by method (all cells) ===\n")
print(as.data.frame(S %>% group_by(method) %>%
  summarise(fits = sum(datasets), failed = sum(failed), degenerate = sum(degenerate),
            capped = sum(capped), .groups = "drop")), row.names = FALSE)

## ---- regressions on log events, log n, log p ------------------------------------
## minimum usable fits for a cell to enter the fit: 20 (pre-registered);
## GEMCOX_E6_MINFITS lowers it for rehearsals only
MIN_FITS <- as.integer(Sys.getenv("GEMCOX_E6_MINFITS", "20"))
fit_law <- function(Sx, yvar, sevar) {
  d <- Sx %>% filter(is.finite(.data[[yvar]]), is.finite(.data[[sevar]]), fits >= MIN_FITS) %>%
    mutate(y = pmin(pmax(.data[[yvar]], 0.01), 0.99),
           se_logit = .data[[sevar]] / (y * (1 - y)), w = 1 / pmax(se_logit, 1e-4)^2,
           logit_y = qlogis(y), lE = log(events), lN = log(n), lP = log(p))
  if (nrow(d) < 6) return(NULL)
  m <- lm(logit_y ~ lE + lN + lP, data = d, weights = w)
  ci <- suppressMessages(confint(m))
  data.frame(term = c("log events", "log n", "log p"),
             estimate = unname(coef(m)[c("lE", "lN", "lP")]),
             lo = unname(ci[c("lE", "lN", "lP"), 1]), hi = unname(ci[c("lE", "lN", "lP"), 2]),
             cells = nrow(d))
}
classify <- function(L) {
  g <- function(t) L[L$term == t, ]
  e <- g("log events"); nn <- g("log n")
  sigE <- e$lo > 0; sigN <- nn$lo > 0
  ## rule (README): one term carries the effect if its CI is above 0 and
  ## the other's CI includes 0 or its estimate is under a third of it
  if (sigE && (!(nn$lo > 0 || nn$hi < 0) || abs(nn$estimate) < abs(e$estimate) / 3)) return("events")
  if (sigN && (!(e$lo > 0 || e$hi < 0) || abs(e$estimate) < abs(nn$estimate) / 3)) return("n")
  if (sigE && sigN) return("both")
  "neither"
}
LAW <- list(); CLS <- list()
for (mode in c("primary (degenerate excluded)", "sensitivity (degenerate at floor)")) {
  yv <- if (grepl("primary", mode)) c(mech = "recovery", prof = "profile") else c(mech = "recovery_floor", prof = "profile_floor")
  for (meth in levels(S$method)) for (mu in c(0, 1)) for (out in c("mechanism", "profile")) {
    if (out == "profile" && mu == 0) next
    Sx <- S %>% filter(method == meth, mu_sep == mu)
    L <- fit_law(Sx, if (out == "mechanism") yv["mech"] else yv["prof"],
                 if (out == "mechanism") "rec_mcse" else "prof_mcse")
    if (is.null(L)) next
    L <- cbind(analysis = mode, method = meth, mu_sep = mu, outcome = out, L)
    LAW[[length(LAW) + 1]] <- L
    CLS[[length(CLS) + 1]] <- data.frame(analysis = mode, method = meth, mu_sep = mu, outcome = out,
                                         carries = classify(L),
                                         expected = if (out == "mechanism") "events" else "n")
  }
}
if (!length(LAW)) stop("no method x outcome had enough cells with >= ", MIN_FITS, " usable fits")
LAW <- do.call(rbind, LAW); CLS <- do.call(rbind, CLS)
CLS$agrees <- CLS$carries == CLS$expected
wcsv(LAW, "tab_E6_law.csv"); wcsv(CLS, "tab_E6_carries.csv")
cat("\n=== logit(recovery) ~ log events + log n + log p (WLS over cells; 95% CI) ===\n")
print(as.data.frame(LAW %>% mutate(across(c(estimate, lo, hi), ~ round(.x, 3)))), row.names = FALSE)
cat("\n=== which carries the effect (expectation: mechanism -> events, profile -> n) ===\n")
print(CLS, row.names = FALSE)
dis <- CLS %>% filter(!agrees)
cat(sprintf("\nRESULTS THAT DISAGREE WITH THE EXPECTATION: %d of %d\n", nrow(dis), nrow(CLS)))
if (nrow(dis)) print(dis, row.names = FALSE)

## replicate-level sensitivity: each fit's own event count
RL <- bind_rows(lapply(levels(X$method), function(meth) bind_rows(lapply(c(0, 1), function(mu) {
  d <- X %>% filter(method == meth, mu_sep == mu, !failed, !degenerate, events > 0) %>%
    mutate(y = qlogis(pmin(pmax(recovery, 0.01), 0.99)), lE = log(events), lN = log(n), lP = log(p))
  if (nrow(d) < 50) return(NULL)
  m <- lm(y ~ lE + lN + lP, data = d); ci <- suppressMessages(confint(m))
  data.frame(method = meth, mu_sep = mu, outcome = "mechanism", term = c("log events", "log n", "log p"),
             estimate = unname(coef(m)[2:4]), lo = ci[2:4, 1], hi = ci[2:4, 2], fits = nrow(d))
}))))
wcsv(RL, "tab_E6_law_replicate.csv")

## ---- matched-events pairs ---------------------------------------------------------
pairs <- E6_CELLS %>% filter(p == 10) %>% select(cell, n, mu_sep, event_rate) %>%
  inner_join(E6_CELLS %>% filter(p == 10) %>% select(cell2 = cell, n2 = n, mu_sep, rate2 = event_rate),
             by = "mu_sep", relationship = "many-to-many") %>%
  filter(n2 == 2 * n, abs(rate2 - event_rate / 2) < 1e-9)
MP <- bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
  a <- S %>% filter(cell == pairs$cell[i]); b <- S %>% filter(cell == pairs$cell2[i])
  inner_join(a, b, by = "method", suffix = c("_n", "_2n")) %>%
    transmute(method, mu_sep = pairs$mu_sep[i], expected_events = pairs$n[i] * pairs$event_rate[i],
              n = pairs$n[i], rate = pairs$event_rate[i], n2 = pairs$n2[i], rate2 = pairs$rate2[i],
              events_n = events_n, events_2n = events_2n,
              mech_diff = recovery_2n - recovery_n, mech_se = sqrt(rec_mcse_n^2 + rec_mcse_2n^2),
              prof_diff = profile_2n - profile_n, prof_se = sqrt(prof_mcse_n^2 + prof_mcse_2n^2))
})) %>% mutate(mech_verdict = ifelse(abs(mech_diff) < 2 * mech_se, "no difference (as expected)",
                                     ifelse(mech_diff > 0, "2n higher (disagrees)", "2n lower (disagrees)")),
               prof_verdict = ifelse(mu_sep == 0, NA,
                                     ifelse(prof_diff > 2 * prof_se, "2n higher (as expected)",
                                            "not higher (disagrees)")))
wcsv(MP, "tab_E6_matched_pairs.csv")
cat("\n=== matched-events pairs: (2n, rate/2) minus (n, rate); expectation: mechanism equal, profile higher at 2n ===\n")
print(as.data.frame(MP %>% mutate(across(c(mech_diff, mech_se, prof_diff, prof_se), ~ round(.x, 3)))),
      row.names = FALSE)

## ---- figures ----------------------------------------------------------------------
g1 <- ggplot(S %>% filter(p == 10), aes(events, recovery, colour = factor(n))) +
  geom_point() + geom_line(aes(group = interaction(n, method))) + scale_x_log10() +
  facet_grid(mu_sep ~ method, labeller = label_both) +
  labs(title = "E6: mechanism recovery against expected events (p = 10)",
       subtitle = "Administrative censoring at a fixed follow-up; EM to tolerance 1e-8. Degenerate fits excluded.",
       x = "Mean events (log scale)", y = "Contrast recovery", colour = "n") +
  theme_minimal(base_size = 10)
ggsave(file.path(FIGURES, "figE6_mechanism_events.pdf"), g1, width = 13, height = 6)
g2 <- ggplot(S %>% filter(p == 10, mu_sep == 1), aes(n, profile, colour = factor(event_rate))) +
  geom_point() + geom_line(aes(group = interaction(event_rate, method))) + scale_x_log10() +
  facet_wrap(~ method, nrow = 1) +
  labs(title = "E6: profile recovery against n (mu_sep = 1, p = 10)",
       x = "n (log scale)", y = "Profile-direction recovery", colour = "event rate") +
  theme_minimal(base_size = 10)
ggsave(file.path(FIGURES, "figE6_profile_n.pdf"), g2, width = 13, height = 4)
cat("\nDone.\n")
