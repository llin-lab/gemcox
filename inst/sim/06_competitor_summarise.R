###############################################################
## 06_competitor_summarise.R  --  latent-class Cox competitor: tables,
## verdicts and the C2 figure (added after sim-freeze-v2)
##
## Implements the analysis pre-registered in README.md ("Added after
## sim-freeze-v2"); written and committed before the C runs. Reads only
## results/; nothing is refitted.
##
## Dataset identity: every C row carries md5 hashes of its training and test
## data. They are checked against results/C0_dataset_hashes.rds (the frozen
## commits' datasets) before any paired difference is computed; a replicate
## that fails the check is dropped from every comparison and counted.
##
## ASCII only.
###############################################################

suppressMessages(library(ggplot2))
Sys.setenv(GEMCOX_WORKERS = "1")   # no refitting here, so no workers
source("00_core.R")
source("C_common.R")

## The frozen summary helpers, evaluated from 03_summarise.R itself (parsed,
## not copied), so their definitions cannot drift.
import_03 <- function(nms) {
  for (e in parse("03_summarise.R", keep.source = FALSE)) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
        as.character(e[[2]]) %in% nms) eval(e, globalenv())
  }
  miss <- nms[!vapply(nms, exists, NA, envir = globalenv(), inherits = FALSE)]
  if (length(miss)) stop("not found in 03_summarise.R: ", paste(miss, collapse = ", "))
}
import_03(c("BLU", "ORG", "GRY", "SKY", "YEL", "INK", "PAL", "th", "wcsv", "binom_mcse",
            "TARGET", "COMPARATORS", "is_degenerate", "score_degenerate",
            "summarise_cells", "wide_recovery", "paired_diffs", "normalised_gain"))

BY <- c("n", "p", "mu_sep", "beta_sep")
GRN <- "#009E73"
PAL_C <- c(PAL, stats::setNames(GRN, LCC))
GAUSSIAN <- c("C1", "C2", "C3a", "C3b", "C4")
regime_of <- function(key, mu_sep) dplyr::case_when(
  key == "C1" ~ "mechanism only",
  key == "C2" & mu_sep == 0 ~ "mechanism only",
  key == "C2" ~ "profile separation",
  key %in% c("C3a", "C3b") ~ "sample size and dimension",
  key == "C4" ~ "application scale",
  key == "C5" ~ "non-Gaussian features")
REGIMES <- c("mechanism only", "profile separation", "sample size and dimension",
             "application scale", "non-Gaussian features")

## ---- verdict rules (README) ------------------------------------------------
verdict_cell <- function(diff, se) dplyr::case_when(
  !is.finite(diff) | !is.finite(se) ~ NA_character_,
  diff > 2 * se ~ "beats", diff < -2 * se ~ "loses", TRUE ~ "matches")
verdict_regime <- function(v) {
  v <- v[!is.na(v)]
  if (!length(v)) return(NA_character_)
  b <- any(v == "beats"); l <- any(v == "loses")
  if (b && !l) "beats" else if (l && !b) "loses" else if (!b && !l) "matches" else "mixed"
}

## ---- loading and dataset verification --------------------------------------
MAN <- readRDS(file.path(RESULTS, "C0_dataset_hashes.rds"))$rows
if (!all(MAN$verified)) stop("C0_dataset_hashes.rds contains unverified datasets")

## Returns the rows and the replicates (cell, rep) that fail verification.
## ref: rows whose hashes the file must match (C5's convergence file is
## checked against C5's primary rows; C5 has no frozen datasets).
load_c <- function(key, file, ref = NULL) {
  E <- C_EXPS[[key]]
  X <- read_results(file)
  if (is.null(X)) return(NULL)
  j <- X[, c("cell", "rep", "method", "failed", "seed", "events", "hash_train", "hash_test")]
  j$seed_ok <- j$seed == mk_seed(E$seed_exp, j$cell, j$rep)
  if (!is.na(E$frozen)) {
    m <- MAN[MAN$experiment == key, c("cell", "rep", "seed_frozen", "events_stored",
                                      "hash_train_frozen", "hash_test_frozen")]
    j <- merge(j, m, by = c("cell", "rep"), all.x = TRUE)
    j$good <- with(j, seed_ok & !is.na(seed_frozen) & seed == seed_frozen &
                     events == events_stored &
                     (failed | (hash_train == hash_train_frozen &
                                (is.na(hash_test_frozen) | hash_test == hash_test_frozen))))
  } else {
    if (!is.null(ref)) {
      r <- unique(ref[!ref$failed, c("cell", "rep", "hash_train", "hash_test")])
      names(r)[3:4] <- c("ref_train", "ref_test")
      j <- merge(j, r, by = c("cell", "rep"), all.x = TRUE)
      j$good <- with(j, seed_ok & (failed | (hash_train == ref_train & hash_test == ref_test)))
    } else {
      j <- j %>% group_by(cell, rep) %>%
        mutate(good = seed_ok & n_distinct(hash_train[!failed]) <= 1 &
                 n_distinct(hash_test[!failed]) <= 1) %>% ungroup()
    }
  }
  j$good[is.na(j$good)] <- FALSE
  bad <- unique(paste(j$cell, j$rep)[!j$good])
  list(rows = X, bad = bad, checked = length(unique(paste(j$cell, j$rep))))
}

ALL <- list(); CONV <- list(); VER <- list()
for (key in names(C_EXPS)) {
  E <- C_EXPS[[key]]
  a <- load_c(key, E$file)
  if (is.null(a)) next
  base <- if (E$fit == "competitor") read_results(E$frozen) else NULL
  R <- bind_rows(base, a$rows) %>% filter(!paste(cell, rep) %in% a$bad) %>%
    mutate(experiment = key, regime = regime_of(key, mu_sep))
  ALL[[key]] <- R
  VER[[length(VER) + 1]] <- data.frame(experiment = key, file = E$file, replicates = a$checked,
                                       failing_verification = length(a$bad))
  cv <- load_c(key, conv_file(E), ref = if (is.na(E$frozen)) a$rows else NULL)
  if (!is.null(cv)) {
    CONV[[key]] <- bind_rows(R, cv$rows) %>% filter(!paste(cell, rep) %in% c(a$bad, cv$bad)) %>%
      mutate(experiment = key, regime = regime_of(key, mu_sep))
    VER[[length(VER) + 1]] <- data.frame(experiment = key, file = conv_file(E),
                                         replicates = cv$checked,
                                         failing_verification = length(cv$bad))
  }
}
if (!length(ALL)) stop("no competitor results found")
VER <- do.call(rbind, VER)
wcsv(VER, "tab_C_verification.csv")
cat("\n=== dataset verification (hashes, seeds, event counts vs C0_dataset_hashes.rds) ===\n")
print(VER, row.names = FALSE)

## ---- paired differences ----------------------------------------------------
pd <- function(R, target, comparators, contrast, mode = "exclude") {
  if (!target %in% R$method || !any(comparators %in% R$method)) return(NULL)
  paired_diffs(R, BY, target = target, comparators = comparators, mode = mode) %>%
    filter(pairs > 0) %>% mutate(target = target, contrast = contrast, scoring = mode)
}
cell_of <- function(R) distinct(R, across(all_of(c(BY, "cell", "experiment", "regime"))))
P <- bind_rows(lapply(names(ALL), function(key) {
  R <- ALL[[key]]
  x <- bind_rows(
    pd(R, TARGET, LCC, "primary"),
    pd(R, TARGET, LCC, "primary", mode = "floor"),
    pd(R, TARGET, LCC_SENS, "gating lambda/10"),
    pd(R, LCC, c("Two-stage (GMM+Cox)", "GeM-Cox (gamma=0)", "Single Cox model"), "ranking"))
  if (!is.null(CONV[[key]])) x <- bind_rows(x,
    pd(CONV[[key]], GEM_TIGHT, LCC_TIGHT, "convergence (both tol 1e-8)"),
    pd(CONV[[key]], LCC_TIGHT, LCC, "competitor: tight minus primary"),
    pd(CONV[[key]], GEM_TIGHT, TARGET, "GeM-Cox: tight minus frozen"))
  left_join(x, cell_of(R), by = BY)
}))
P <- P %>% mutate(verdict = verdict_cell(diff, diff_mcse))
wcsv(P, "tab_C_paired.csv")

## per-cell verdicts, GeM-Cox minus competitor
VC <- P %>% filter(contrast %in% c("primary", "gating lambda/10", "convergence (both tol 1e-8)"),
                   scoring == "exclude") %>%
  select(experiment, regime, cell, all_of(BY), contrast, pairs, diff, diff_mcse, win, verdict)
wcsv(VC, "tab_C_verdicts_cell.csv")
fmt_cells <- function(X) as.data.frame(X %>% transmute(experiment, cell, n, p, mu_sep, beta_sep,
  pairs, diff = round(diff, 3), mcse = round(diff_mcse, 3), win = round(win, 2), verdict))
cat("\n=== PRIMARY: paired difference in recovery, GeM-Cox (gamma = 1) minus the competitor ===\n")
cat("    (degenerate fits excluded; verdict: beats / loses if beyond 2 MCSE, else matches)\n")
print(fmt_cells(VC %>% filter(contrast == "primary")), row.names = FALSE)

## per-regime verdicts
VR <- VC %>% group_by(contrast, regime) %>%
  summarise(cells = n(), beats = sum(verdict == "beats", na.rm = TRUE),
            loses = sum(verdict == "loses", na.rm = TRUE),
            matches = sum(verdict == "matches", na.rm = TRUE),
            verdict = verdict_regime(verdict), .groups = "drop") %>%
  mutate(regime = factor(regime, REGIMES)) %>% arrange(contrast, regime)
wcsv(VR, "tab_C_verdicts_regime.csv")
cat("\n=== PRIMARY: verdict per regime (GeM-Cox gamma = 1 vs the competitor) ===\n")
print(as.data.frame(VR %>% filter(contrast == "primary") %>% select(-contrast)), row.names = FALSE)

## sensitivity arm 1: gating penalty lambda / 10, over the cells where it ran
S1 <- VC %>% filter(contrast == "gating lambda/10")
if (nrow(S1)) {
  on <- paste(S1$experiment, S1$cell)
  cmp <- VC %>% filter(contrast %in% c("primary", "gating lambda/10"),
                       paste(experiment, cell) %in% on) %>%
    group_by(regime, contrast) %>% summarise(verdict = verdict_regime(verdict), cells = n(),
                                             .groups = "drop") %>%
    tidyr::pivot_wider(names_from = contrast, values_from = verdict) %>%
    mutate(changed = primary != `gating lambda/10`)
  wcsv(cmp, "tab_C_sensitivity_gating.csv")
  cat("\n=== SENSITIVITY 1: gating penalty lambda / 10 (cells where the arm ran) ===\n")
  print(fmt_cells(S1), row.names = FALSE)
  print(as.data.frame(cmp), row.names = FALSE)
  cat(sprintf("  Any per-regime verdict changed: %s\n", any(cmp$changed, na.rm = TRUE)))
}

## sensitivity arm 2: convergence
S2 <- VC %>% filter(contrast == "convergence (both tol 1e-8)")
if (nrow(S2)) {
  cmp2 <- VR %>% filter(contrast %in% c("primary", "convergence (both tol 1e-8)")) %>%
    select(regime, contrast, verdict) %>%
    tidyr::pivot_wider(names_from = contrast, values_from = verdict) %>%
    mutate(changed = primary != `convergence (both tol 1e-8)`)
  wcsv(cmp2, "tab_C_sensitivity_convergence.csv")
  cat("\n=== SENSITIVITY 2: both methods to tolerance 1e-8 (GeM-Cox tight minus competitor tight) ===\n")
  print(fmt_cells(S2), row.names = FALSE)
  print(as.data.frame(cmp2), row.names = FALSE)
  cat(sprintf("  Any per-regime verdict changed: %s\n", any(cmp2$changed, na.rm = TRUE)))
  D2 <- P %>% filter(contrast %in% c("competitor: tight minus primary", "GeM-Cox: tight minus frozen"))
  cat("\n  descriptive: the effect of the tolerance on each method (paired, degenerate excluded):\n")
  print(as.data.frame(D2 %>% transmute(contrast, experiment, cell, mu_sep, beta_sep, n, p, pairs,
        diff = round(diff, 3), mcse = round(diff_mcse, 3))), row.names = FALSE)
}

## floor-scored sensitivity of the primary comparison
FL <- P %>% filter(contrast == "primary") %>%
  select(experiment, cell, scoring, diff, diff_mcse, verdict) %>%
  tidyr::pivot_wider(names_from = scoring, values_from = c(diff, diff_mcse, verdict))
changed_fl <- FL %>% filter(verdict_exclude != verdict_floor)
cat(sprintf("\n  floor-scored degenerate fits: %d cell verdict(s) differ from the primary scoring\n",
            nrow(changed_fl)))
if (nrow(changed_fl)) print(as.data.frame(changed_fl), row.names = FALSE)

## ranking: the competitor against the other comparators
cat("\n=== RANKING: the competitor minus each comparator (paired, degenerate excluded) ===\n")
print(as.data.frame(P %>% filter(contrast == "ranking") %>%
  transmute(experiment, cell, mu_sep, beta_sep, n, p, comparator, pairs,
            diff = round(diff, 3), mcse = round(diff_mcse, 3))), row.names = FALSE)

## normalised gains: GeM-Cox and the competitor over gamma = 0, relative to the oracle
NG <- bind_rows(lapply(names(ALL), function(key) bind_rows(
  normalised_gain(ALL[[key]], BY, arm = TARGET) %>% mutate(arm = TARGET),
  normalised_gain(ALL[[key]], BY, arm = LCC) %>% mutate(arm = LCC)) %>%
  mutate(experiment = key)))
wcsv(NG, "tab_C_normalised.csv")
cat("\n=== normalised gain (arm - gamma0) / (oracle - gamma0) ===\n")
print(as.data.frame(NG %>% transmute(experiment, mu_sep, beta_sep, n, p, arm, pairs,
      norm_gain = round(norm_gain, 3), mcse = round(norm_mcse, 3))), row.names = FALSE)

## ---- secondary metrics -----------------------------------------------------
paired_metric <- function(R, target, comparator, metric) {
  score_degenerate(R, "exclude") %>% filter(!failed, method %in% c(target, comparator)) %>%
    select(all_of(c(BY, "rep", "method", metric))) %>%
    tidyr::pivot_wider(names_from = method, values_from = all_of(metric)) %>%
    mutate(d = .data[[target]] - .data[[comparator]]) %>%
    group_by(across(all_of(BY))) %>%
    summarise(pairs = sum(is.finite(d)), diff_mcse = mcse(d), diff = mean(d, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(target = target, comparator = comparator, metric = metric)
}
PM <- bind_rows(lapply(names(ALL), function(key) bind_rows(lapply(c("mse_eta", "ari"), function(mt)
  paired_metric(ALL[[key]], TARGET, LCC, mt))) %>% mutate(experiment = key)))
wcsv(PM, "tab_C_paired_secondary.csv")
DS <- bind_rows(lapply(names(ALL), function(key)
  summarise_cells(ALL[[key]], c(BY, "method")) %>% mutate(experiment = key)))
wcsv(DS, "tab_C_descriptive.csv")
cat("\n=== descriptive: recovery, ARI, MSE of the linear predictor, sharpness (degenerate excluded) ===\n")
print(as.data.frame(DS %>% filter(method %in% c(TARGET, LCC, "Oracle (true labels)")) %>%
  transmute(experiment, mu_sep, beta_sep, n, p, method, reps, recovery = round(recovery, 3),
            rec_mcse = round(rec_mcse, 3), ARI = round(ari, 3), mse_eta = round(mse_eta, 3),
            sharpness = round(sharpness, 3))), row.names = FALSE)
cat("\n  paired MSE and ARI, GeM-Cox (gamma = 1) minus the competitor:\n")
print(as.data.frame(PM %>% transmute(experiment, metric, mu_sep, beta_sep, n, p, pairs,
      diff = round(diff, 3), mcse = round(diff_mcse, 3))), row.names = FALSE)

## ---- degeneracy, failure and early stopping --------------------------------
DGC <- bind_rows(lapply(names(ALL), function(key) {
  R <- if (!is.null(CONV[[key]])) CONV[[key]] else ALL[[key]]
  R %>% group_by(experiment, cell, across(all_of(BY)), method) %>%
    summarise(fits = n(), failed = sum(failed),
              degenerate = sum(!failed & !is.na(degenerate) & degenerate),
              early_stop = sum(!failed & !is.na(iterations) & iterations <= 3),
              not_converged = sum(!failed & !is.na(converged) & !converged),
              gate_failed = sum(!failed & !is.na(gate_ok) & !gate_ok),
              .groups = "drop")
})) %>% mutate(ok = fits - failed,
               failure_rate = failed / fits, failure_mcse = binom_mcse(failed, fits),
               degenerate_rate = degenerate / ok, degenerate_mcse = binom_mcse(degenerate, ok),
               early_stop_rate = early_stop / ok)
wcsv(DGC, "tab_C_degeneracy.csv")
cat("\n=== failure, degeneracy and early stopping (<= 3 EM iterations), GeM-Cox and the competitor ===\n")
print(as.data.frame(DGC %>% filter(method %in% c(TARGET, LCC, GEM_TIGHT, LCC_TIGHT)) %>%
  transmute(experiment, cell, method, fits, failure = round(failure_rate, 3),
            degenerate = round(degenerate_rate, 3), early_stop = round(early_stop_rate, 3),
            not_converged, gate_failed)), row.names = FALSE)

## ---- expectations (README) ---------------------------------------------------
cat("\n=== pre-registered expectations ===\n")
X1 <- VC %>% filter(contrast == "primary", experiment == "C2", mu_sep > 0)
if (nrow(X1)) {
  fails <- X1 %>% filter(verdict == "loses")
  cat(sprintf("  1. C2, mu_sep > 0: GeM-Cox >= competitor in all %d cells. Cells failing (difference < -2 MCSE): %d. Expectation %s.\n",
              nrow(X1), nrow(fails), if (nrow(fails)) "FAILED" else "met"))
}
if (!is.null(ALL$C4)) {
  s <- summarise_cells(ALL$C4, c(BY, "method")) %>%
    filter(method %in% c(TARGET, LCC, "Oracle (true labels)")) %>%
    group_by(across(all_of(BY))) %>%
    mutate(oracle = recovery[method == "Oracle (true labels)"], floor = random_floor(p),
           near_floor = recovery < (floor + oracle) / 2) %>% ungroup() %>%
    filter(method != "Oracle (true labels)")
  wcsv(s, "tab_C_expectation4.csv")
  cat("  4. C4, both near the floor (mean recovery closer to sqrt(2/(9 pi)) = 0.266 than to the oracle's):\n")
  print(as.data.frame(s %>% transmute(mu_sep, method, recovery = round(recovery, 3),
        mcse = round(rec_mcse, 3), oracle = round(oracle, 3), near_floor)), row.names = FALSE)
  cat(sprintf("     Expectation %s.\n", if (all(s$near_floor)) "met" else "not met"))
}
flag <- VC %>% filter(contrast == "primary", experiment %in% GAUSSIAN, verdict == "loses")
flag2 <- VC %>% filter(contrast == "convergence (both tol 1e-8)", experiment %in% GAUSSIAN,
                       verdict == "loses")
cat(sprintf("  FLAG: competitor beats GeM-Cox in a Gaussian cell (primary): %s%s\n",
            if (nrow(flag)) "YES in " else "no",
            if (nrow(flag)) paste(paste(flag$experiment, "cell", flag$cell), collapse = ", ") else ""))
if (nrow(S2)) {
  cat(sprintf("        ... and under the convergence arm: %s%s\n",
              if (nrow(flag2)) "YES in " else "no",
              if (nrow(flag2)) paste(paste(flag2$experiment, "cell", flag2$cell), collapse = ", ") else ""))
}

## ---- consistency with the frozen E2 table ------------------------------------
if (!is.null(ALL$C2) && file.exists(file.path(RESULTS, "tab_E2_paired.csv"))) {
  old <- read.csv(file.path(RESULTS, "tab_E2_paired.csv"))
  new <- paired_diffs(ALL$C2, BY) %>% as.data.frame()
  k <- merge(old, new, by = c(BY, "comparator"))
  cat(sprintf("\n  consistency: frozen comparators' paired differences in C2 reproduce tab_E2_paired.csv: %s (%d rows)\n",
              isTRUE(all.equal(k$diff.x, k$diff.y, tolerance = 1e-12)) &&
                isTRUE(all.equal(k$pairs.x, k$pairs.y)), nrow(k)))
}

## ---- figure: C2 as a four-panel comparison matching fig2_paired.pdf ---------
if (!is.null(ALL$C2)) {
  CMP <- c(COMPARATORS, LCC)
  P2 <- paired_diffs(ALL$C2, BY, comparators = CMP)
  g <- ggplot(P2 %>% mutate(comparator = factor(comparator, CMP)),
              aes(factor(beta_sep), diff, colour = comparator, group = comparator)) +
    geom_hline(yintercept = 0, linetype = 2, colour = GRY) +
    geom_line(linewidth = 1) +
    geom_pointrange(aes(ymin = diff - diff_mcse, ymax = diff + diff_mcse), size = .35) +
    scale_colour_manual(values = PAL_C[CMP], labels = paste("minus", CMP)) +
    facet_wrap(~ mu_sep, nrow = 1, labeller = label_both) +
    labs(title = "Paired gain in contrast recovery, with the latent-class Cox competitor",
         subtitle = paste("C2 (E2's datasets): n = 800, p = 10. Bars: +/- 1 Monte Carlo SE",
                          "of the paired difference. Degenerate fits excluded."),
         x = "Difference between survival mechanisms (beta_sep)",
         y = "Paired difference in contrast recovery\n(GeM-Cox gamma = 1 minus comparator)") + th
  ggsave(file.path(FIGURES, "figC2_paired.pdf"), g, width = 12, height = 4.6)
  cat("\n  saved figC2_paired.pdf\n")
}
cat("\nDone.\n")
