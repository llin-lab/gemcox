###############################################################
## 03_summarise.R  --  summaries, Monte Carlo SEs, manuscript figures
##
## Reads only from results/. Nothing here refits anything.
##
## PRIMARY METRIC (decided after the pilot, before the full run):
## paired per-replicate differences in contrast recovery, GeM-Cox
## (gamma = 1) minus each comparator, computed on the same simulated
## dataset. MCSE = SD of the paired differences / sqrt(number of pairs).
## Normalised version: (gamma1 - gamma0) / (oracle - gamma0), a ratio of
## paired means with a delta-method MCSE.
##
## The random-direction floor sqrt(2 / (pi p)) is DESCRIPTIVE CONTEXT ONLY.
## The estimators' null behaviour is not isotropic: in the pilot the
## two-stage contrast aligned with the shared coefficient (orthogonal to
## the true contrast by design) and fell below the floor, while GeM-Cox
## gamma = 0 rose slightly above it at beta_sep = 3. Percentages "of oracle
## above the floor" are therefore reported only as description.
##
## Monte Carlo standard errors are reported for every quantity (ADEMP;
## Morris, White & Crowther, Stat Med 2019). Inside summarise() every
## dispersion is computed BEFORE its column is overwritten by the mean.
## Figure titles describe what is plotted; they do not state conclusions.
##
## DEGENERATE FITS: a fitted cluster with fewer than 10 subjects or fewer
## than 5 events (see 00_core.R). PRIMARY analysis excludes degenerate fits
## (a pair is dropped if either member is degenerate) and the degeneracy
## rate is reported for every method in every cell. SENSITIVITY analysis
## scores degenerate fits at the random-direction floor instead of the
## zero-coefficient score they receive from the harness.
##
## ASCII only.
###############################################################

suppressMessages(library(ggplot2))
Sys.setenv(GEMCOX_WORKERS = "1")   # no refitting here, so no workers
source("00_core.R")

BLU <- "#0072B2"; ORG <- "#D55E00"; GRY <- "#8C9196"
SKY <- "#56B4E9"; YEL <- "#E69F00"; INK <- "#1A1A1A"
PAL <- c("Single Cox model" = GRY, "Two-stage (GMM+Cox)" = SKY,
         "GeM-Cox (gamma=0)" = YEL, "GeM-Cox (gamma=1)" = ORG,
         "Oracle (true labels)" = BLU)
th <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "#555555", size = 9),
        legend.position = "top", legend.title = element_blank(),
        panel.grid.minor = element_blank())
wcsv <- function(x, f) write.csv(x, file.path(RESULTS, f), row.names = FALSE)
binom_mcse <- function(k, n) sqrt((k / n) * (1 - k / n) / n)
FLOOR_NOTE <- paste0("Dashed line: isotropic random-direction reference sqrt(2/(pi p)), ",
                     "descriptive only; the estimators' null behaviour is not isotropic.")
TARGET <- "GeM-Cox (gamma=1)"
COMPARATORS <- c("GeM-Cox (gamma=0)", "Two-stage (GMM+Cox)", "Single Cox model")

## ---- degenerate fits ------------------------------------------------------
is_degenerate <- function(R) !is.na(R$degenerate) & R$degenerate
## mode: "exclude" (primary), "floor" (sensitivity), "as_scored" (reference)
score_degenerate <- function(R, mode = c("exclude", "floor", "as_scored")) {
  mode <- match.arg(mode)
  if (!"degenerate" %in% names(R) || mode == "as_scored") return(R)
  dg <- is_degenerate(R)
  if (mode == "exclude") return(R[!dg, , drop = FALSE])
  R$recovery[dg] <- random_floor(R$p[dg])
  R
}

## ---- descriptive recovery (SD first, then the mean) ----------------------
## Failed and (primary analysis) degenerate fits are excluded here; both
## are counted in the diagnostics and degeneracy tables.
summarise_cells <- function(R, by = c("n", "p", "mu_sep", "beta_sep", "method")) {
  by <- intersect(by, names(R))
  score_degenerate(R, "exclude") %>% filter(!failed) %>% group_by(across(all_of(by))) %>%
    summarise(reps          = n(),
              rec_mcse      = mcse(recovery),
              ari_mcse      = mcse(ari),
              mse_mcse      = mcse(mse_eta),
              sharp_mcse    = mcse(sharpness),
              recovery      = mean(recovery, na.rm = TRUE),
              ari           = mean(ari, na.rm = TRUE),
              mse_eta       = mean(mse_eta, na.rm = TRUE),
              sharpness     = mean(sharpness, na.rm = TRUE),
              events        = mean(events),
              .groups = "drop")
}
## descriptive only: % of oracle recovery above the isotropic floor
add_pct <- function(S, group = c("n", "p", "mu_sep", "beta_sep")) {
  group <- intersect(group, names(S))
  S %>% group_by(across(all_of(group))) %>%
    mutate(oracle   = recovery[method == "Oracle (true labels)"],
           pct      = pct_of_oracle(recovery, oracle, p),
           pct_mcse = 100 * rec_mcse / (oracle - random_floor(p))) %>%
    ungroup()
}

## ---- PRIMARY: paired per-replicate differences ---------------------------
wide_recovery <- function(R, by, mode = "exclude") {
  score_degenerate(R, mode) %>% filter(!failed) %>%
    select(all_of(c(by, "rep", "method", "recovery"))) %>%
    pivot_wider(names_from = method, values_from = recovery)
}
## recovery(target) - recovery(comparator) on the same replicate; the MCSE
## is the SD of the paired differences over sqrt(pairs)
paired_diffs <- function(R, by, target = TARGET, comparators = COMPARATORS,
                         mode = "exclude") {
  W <- wide_recovery(R, by, mode)
  bind_rows(lapply(intersect(comparators, names(W)), function(cm) {
    W %>% mutate(d = .data[[target]] - .data[[cm]]) %>%
      group_by(across(all_of(by))) %>%
      summarise(comparator = cm,
                pairs = sum(is.finite(d)),
                diff_mcse = mcse(d),
                win = mean(d > 0, na.rm = TRUE),
                diff = mean(d, na.rm = TRUE),
                .groups = "drop")
  }))
}
## (gamma1 - gamma0) / (oracle - gamma0) as a ratio of paired means over
## replicates where all three fits succeeded. Delta-method MCSE:
##   var(R) ~ (var(a) - 2 R cov(a, b) + R^2 var(b)) / (pairs * mean(b)^2)
normalised_gain <- function(R, by, arm = "GeM-Cox (gamma=1)") {
  W <- wide_recovery(R, by)
  W %>% mutate(a = .data[[arm]] - .data[["GeM-Cox (gamma=0)"]],
               b = .data[["Oracle (true labels)"]] - .data[["GeM-Cox (gamma=0)"]]) %>%
    filter(is.finite(a), is.finite(b)) %>%
    group_by(across(all_of(by))) %>%
    summarise(pairs = n(), va = var(a), vb = var(b), cab = cov(a, b),
              ma = mean(a), mb = mean(b), .groups = "drop") %>%
    mutate(norm_gain = ma / mb,
           norm_mcse = sqrt(pmax(va - 2 * norm_gain * cab + norm_gain^2 * vb, 0) / pairs) /
             abs(mb)) %>%
    select(-va, -vb, -cab)
}
print_paired <- function(P, by, digits = 3) {
  print(as.data.frame(P %>% transmute(across(all_of(by)), comparator, pairs,
        diff = round(diff, digits), mcse = round(diff_mcse, digits), win = round(win, 2))),
        row.names = FALSE)
}
paired_plot <- function(P, x, facet = NULL, title, subtitle, xlab) {
  g <- ggplot(P %>% mutate(comparator = factor(comparator, COMPARATORS)),
              aes(factor(.data[[x]]), diff, colour = comparator, group = comparator)) +
    geom_hline(yintercept = 0, linetype = 2, colour = GRY) +
    geom_line(linewidth = 1) +
    geom_pointrange(aes(ymin = diff - diff_mcse, ymax = diff + diff_mcse), size = .35) +
    scale_colour_manual(values = PAL[COMPARATORS],
                        labels = paste("minus", COMPARATORS)) +
    labs(title = title, subtitle = subtitle, x = xlab,
         y = "Paired difference in contrast recovery\n(GeM-Cox gamma = 1 minus comparator)") + th
  if (!is.null(facet)) g <- g + facet_wrap(as.formula(paste("~", facet)), nrow = 1,
                                           labeller = label_both)
  g
}

## ================= diagnostics across experiments =====================
diag_files <- c("E1_existence.rds", "E1b_normalisation.rds", "E2_comparison.rds", "E2b_a_adaptive.rds", "E3a_n.rds", "E3b_p.rds",
                "E4a_two_scales.rds", "E4b_mechanism.rds", "E4c_softness.rds",
                "S1_normalisation.rds", "PH1_alignment.rds")
DG <- bind_rows(lapply(diag_files, function(f) {
  R <- read_results(f)
  if (is.null(R)) return(NULL)
  R %>% group_by(method) %>%
    summarise(file = f, rows = n(), failed = sum(failed),
              not_converged = sum(!is.na(converged) & !converged),
              guard_bound = sum(!is.na(guard_bound) & guard_bound),
              mean_secs = mean(secs), .groups = "drop")
}))
if (nrow(DG)) {
  cat("\n=== fit diagnostics (failed fits are excluded from all summaries) ===\n")
  print(as.data.frame(DG %>% mutate(mean_secs = round(mean_secs, 2))), row.names = FALSE)
  wcsv(DG, "tab_diagnostics.csv")
}

## ================= degeneracy rates, every method, every cell ==========
degen_files <- c(E1 = "E1_existence.rds", E1b = "E1b_normalisation.rds",
                 E2 = "E2_comparison.rds", E2b = "E2b_a_adaptive.rds",
                 E3a = "E3a_n.rds", E3b = "E3b_p.rds",
                 E4a = "E4a_two_scales.rds", E4b = "E4b_mechanism.rds",
                 E4c = "E4c_softness.rds", S1 = "S1_normalisation.rds",
                 PH1 = "PH1_alignment.rds")
DR <- bind_rows(lapply(names(degen_files), function(e) {
  R <- read_results(degen_files[[e]])
  if (is.null(R) || !"degenerate" %in% names(R)) return(NULL)
  R %>% filter(!failed) %>% group_by(n, p, mu_sep, beta_sep, method) %>%
    summarise(experiment = e, fits = n(), degenerate = sum(!is.na(degenerate) & degenerate),
              zero_coef_fits = sum(n_zero_fits > 0, na.rm = TRUE), .groups = "drop") %>%
    mutate(rate = degenerate / fits, rate_mcse = binom_mcse(degenerate, fits))
}))
if (nrow(DR)) {
  wcsv(DR, "tab_degeneracy.csv")
  cat("\n=== degeneracy rate (cluster < 10 subjects or < 5 events), every method, every cell ===\n")
  print(as.data.frame(DR %>% select(experiment, n, p, mu_sep, beta_sep, method, rate) %>%
        mutate(rate = round(rate, 3)) %>%
        pivot_wider(names_from = method, values_from = rate)), row.names = FALSE)
  cat("  (full table with counts, MCSEs and zero-coefficient fits: tab_degeneracy.csv)\n")
}

## paired gains over two-stage under the three scorings of degenerate fits
two_stage_gains <- function(R, by, label) {
  bind_rows(lapply(c(primary_exclude = "exclude", sensitivity_floor = "floor",
                     as_scored_zero = "as_scored"), function(md) {
    paired_diffs(R, by, comparators = "Two-stage (GMM+Cox)", mode = md) %>%
      mutate(scoring = names(which(c(primary_exclude = "exclude", sensitivity_floor = "floor",
                                     as_scored_zero = "as_scored") == md)))
  })) %>% mutate(experiment = label)
}

## ================= E1: existence ======================================
E1 <- read_results("E1_existence.rds")
if (!is.null(E1)) {
  by <- c("n", "p", "mu_sep", "beta_sep")
  P <- paired_diffs(E1, by); NG <- normalised_gain(E1, by)
  wcsv(P, "tab_E1_paired.csv"); wcsv(NG, "tab_E1_normalised.csv")
  cat("\n=== E1 existence (mu_sep = 0): PRIMARY, paired differences in recovery ===\n")
  print_paired(P, "beta_sep")
  cat("\n  normalised gain (gamma1 - gamma0) / (oracle - gamma0):\n")
  print(as.data.frame(NG %>% transmute(beta_sep, pairs, norm_gain = round(norm_gain, 3),
        mcse = round(norm_mcse, 3))), row.names = FALSE)

  ## P2 as recorded in advance (NOT reproduced in the pilot; E1 unchanged).
  ## Reported here for completeness with the full replication.
  G <- P %>% filter(comparator == "GeM-Cox (gamma=0)") %>%
    mutate(scaled = diff / beta_sep^2, scaled_mcse = diff_mcse / beta_sep^2)
  cvs <- with(G[G$beta_sep >= 1, ], sd(scaled) / mean(scaled))
  cat("\n  P2 (gain / beta_sep^2 roughly constant), recorded in advance:\n")
  print(as.data.frame(G %>% transmute(beta_sep, gain = round(diff, 4),
        scaled = round(scaled, 4), scaled_mcse = round(scaled_mcse, 4))), row.names = FALSE)
  cat(sprintf("  CV of gain/beta_sep^2 over beta_sep >= 1: %.2f (criterion < 0.5; %d pairs per cell)\n",
              cvs, min(G$pairs)))
  gb <- G[G$beta_sep >= 1, ]
  cat(sprintf(paste0(
    "  P2 depends on replication. At pilot size the criterion recorded in advance failed\n",
    "  (CV 0.60, 20 replicates). At %d replicates it is %s (CV %.2f %s 0.5). Even so, the\n",
    "  gain divided by beta_sep^2 %s: %s at beta_sep = %s (MCSE %s), and it is %.3f at\n",
    "  beta_sep = 0.5.\n"),
    min(G$pairs), if (cvs < 0.5) "met" else "not met", cvs, if (cvs < 0.5) "<" else ">=",
    if (all(diff(gb$scaled) < 0)) "still declines" else "is not monotone",
    paste(sprintf("%.3f", gb$scaled), collapse = ", "), paste(gb$beta_sep, collapse = ", "),
    paste(sprintf("%.3f", gb$scaled_mcse), collapse = ", "), G$scaled[G$beta_sep == 0.5]))

  S <- add_pct(summarise_cells(E1))
  wcsv(S, "tab_E1_recovery.csv")
  cat("\n  descriptive recovery (floor-relative % is context only):\n")
  print(as.data.frame(S %>% transmute(beta_sep, method, reps,
    recovery = round(recovery, 3), mcse = round(rec_mcse, 3),
    pct_of_oracle_above_floor = round(pct, 1), ARI = round(ari, 3),
    ari_mcse = round(ari_mcse, 3))), row.names = FALSE)

  ggsave(file.path(FIGURES, "fig1_paired.pdf"),
         paired_plot(P, "beta_sep",
                     title = "Paired gain in contrast recovery when feature profiles are identical",
                     subtitle = "mu_sep = 0, n = 400, p = 10. Bars: +/- 1 Monte Carlo SE of the paired difference.",
                     xlab = "Difference between survival mechanisms (beta_sep)"),
         width = 8, height = 4.8)
  FL <- random_floor(10)
  g <- ggplot(S %>% mutate(method = factor(method, METHOD_LEVELS)),
              aes(factor(beta_sep), recovery, colour = method, group = method)) +
    geom_hline(yintercept = FL, linetype = 2, colour = INK) +
    geom_line(linewidth = 1) +
    geom_pointrange(aes(ymin = recovery - rec_mcse, ymax = recovery + rec_mcse), size = .4) +
    scale_colour_manual(values = PAL) + ylim(0, 1) +
    labs(title = "Contrast recovery by method (descriptive)",
         subtitle = paste("mu_sep = 0, n = 400, p = 10.", FLOOR_NOTE),
         x = "Difference between survival mechanisms (beta_sep)",
         y = "Recovery of the true contrast") + th
  ggsave(file.path(FIGURES, "fig1_recovery_context.pdf"), g, width = 8, height = 4.8)
}

## ================= E1b: normalisation on E1's datasets (separate) ======
E1b <- read_results("E1b_normalisation.rds")
if (!is.null(E1b)) {
  by <- c("n", "p", "mu_sep", "beta_sep")
  ARMS <- c("GeM-Cox (gamma=1, normalize=TRUE)", "GeM-Cox (gamma=1, adaptive)")
  CMP_B <- c("GeM-Cox (gamma=1)", COMPARATORS)
  P <- bind_rows(lapply(ARMS, function(a) paired_diffs(E1b, by, target = a,
                                                       comparators = CMP_B) %>%
                          mutate(arm = a)))
  NG <- bind_rows(lapply(c("GeM-Cox (gamma=1)", ARMS), function(a)
    normalised_gain(E1b, by, arm = a) %>% mutate(arm = a)))
  wcsv(P, "tab_E1b_paired.csv"); wcsv(NG, "tab_E1b_normalised.csv")
  cat("\n=== E1b (reported separately from E1): normalize = TRUE and adaptive arms on E1's datasets ===\n")
  cat("  paired differences, arm minus comparator (primary: degenerate fits excluded):\n")
  print(as.data.frame(P %>% transmute(arm, beta_sep, comparator, pairs,
        diff = round(diff, 3), mcse = round(diff_mcse, 3), win = round(win, 2))),
        row.names = FALSE)
  cat("\n  normalised gain (arm - gamma0) / (oracle - gamma0):\n")
  print(as.data.frame(NG %>% transmute(arm, beta_sep, pairs, norm_gain = round(norm_gain, 3),
        mcse = round(norm_mcse, 3))), row.names = FALSE)
  AD <- E1b %>% filter(!failed, method == "GeM-Cox (gamma=1, adaptive)") %>%
    group_by(beta_sep) %>%
    summarise(datasets = n(), chose_TRUE = mean(adaptive_normalize),
              profile_K1 = mean(profile_K == 1), joint_K1 = mean(joint_K == 1),
              rules_agree = mean(profile_K == joint_K), .groups = "drop")
  wcsv(AD, "tab_E1b_adaptive_choice.csv")
  cat("\n  adaptive rule: share choosing normalize = TRUE (profile part selects K = 1),\n")
  cat("  and agreement with the overall joint criterion:\n")
  print(as.data.frame(AD %>% mutate(across(where(is.double), ~round(.x, 3)))), row.names = FALSE)
  if (!is.null(E1)) {
    key <- c("cell", "rep", "method")
    chk <- inner_join(E1 %>% select(all_of(key), recovery),
                      E1b %>% select(all_of(key), recovery), by = key)
    cat(sprintf("  consistency: E1 methods refitted in E1b reproduce E1 exactly: %s (%d fits)\n",
                identical(chk$recovery.x, chk$recovery.y), nrow(chk)))
  }
}

## ================= E2: comparison =====================================
E2 <- read_results("E2_comparison.rds")
if (!is.null(E2)) {
  by <- c("n", "p", "mu_sep", "beta_sep")
  P <- paired_diffs(E2, by); NG <- normalised_gain(E2, by)
  wcsv(P, "tab_E2_paired.csv"); wcsv(NG, "tab_E2_normalised.csv")
  cat("\n=== E2 comparison (n = 800, p = 10): PRIMARY, paired differences in recovery ===\n")
  print_paired(P, c("mu_sep", "beta_sep"))
  cat("\n  normalised gain (gamma1 - gamma0) / (oracle - gamma0):\n")
  print(as.data.frame(NG %>% transmute(mu_sep, beta_sep, pairs,
        norm_gain = round(norm_gain, 3), mcse = round(norm_mcse, 3))), row.names = FALSE)

  S <- add_pct(summarise_cells(E2))
  wcsv(S, "tab_E2_recovery.csv")
  cat("\n  descriptive recovery, ARI and MSE of the linear predictor:\n")
  print(as.data.frame(S %>% transmute(mu_sep, beta_sep, method,
    recovery = round(recovery, 3), mcse = round(rec_mcse, 3),
    pct_of_oracle_above_floor = round(pct, 1), ARI = round(ari, 3),
    mse_eta = round(mse_eta, 3), mse_mcse = round(mse_mcse, 3))), row.names = FALSE)

  ggsave(file.path(FIGURES, "fig2_paired.pdf"),
         paired_plot(P, "beta_sep", facet = "mu_sep",
                     title = "Paired gain in contrast recovery across profile and mechanism separation",
                     subtitle = "n = 800, p = 10. Bars: +/- 1 Monte Carlo SE of the paired difference.",
                     xlab = "Difference between survival mechanisms (beta_sep)"),
         width = 12, height = 4.6)
}

## ================= E3: envelope =======================================
E3a <- read_results("E3a_n.rds"); E3b <- read_results("E3b_p.rds")
if (!is.null(E3a) && !is.null(E3b)) {
  Pa <- paired_diffs(E3a, c("n", "p")) %>% left_join(
    E3a %>% group_by(n) %>% summarise(events = mean(events), .groups = "drop"), by = "n")
  Pb <- paired_diffs(E3b, c("n", "p"))
  NGa <- normalised_gain(E3a, c("n", "p")); NGb <- normalised_gain(E3b, c("n", "p"))
  wcsv(bind_rows(Pa, Pb), "tab_E3_paired.csv"); wcsv(bind_rows(NGa, NGb), "tab_E3_normalised.csv")
  wcsv(bind_rows(add_pct(summarise_cells(E3a)), add_pct(summarise_cells(E3b))),
       "tab_E3_recovery.csv")
  cat("\n=== E3 envelope: PRIMARY, paired differences in recovery ===\n")
  cat("  (a) sample size (p = 10):\n"); print_paired(Pa, "n")
  cat("  (b) dimension (n = 800):\n"); print_paired(Pb, "p")
  cat("\n  normalised gain (gamma1 - gamma0) / (oracle - gamma0):\n")
  print(as.data.frame(bind_rows(NGa, NGb) %>% transmute(n, p, pairs,
        norm_gain = round(norm_gain, 3), mcse = round(norm_mcse, 3))), row.names = FALSE)
  G0 <- bind_rows(
    Pa %>% filter(comparator == "GeM-Cox (gamma=0)") %>%
      transmute(x = events, diff, diff_mcse, facet = "Number of events (n varied, p = 10)"),
    Pb %>% filter(comparator == "GeM-Cox (gamma=0)") %>%
      transmute(x = p, diff, diff_mcse, facet = "Number of features (p varied, n = 800)"))
  g <- ggplot(G0, aes(x, diff)) +
    geom_hline(yintercept = 0, linetype = 2, colour = GRY) +
    geom_line(colour = ORG, linewidth = 1.1) +
    geom_pointrange(aes(ymin = diff - diff_mcse, ymax = diff + diff_mcse), colour = ORG, size = .4) +
    facet_wrap(~ facet, scales = "free_x") + scale_x_log10() +
    labs(title = "Paired gain from survival weighting by events and dimension",
         subtitle = "GeM-Cox gamma = 1 minus gamma = 0, same datasets. mu_sep = 0.5, beta_sep = 2. Bars: +/- 1 MCSE.",
         x = NULL, y = "Paired difference in recovery") + th
  ggsave(file.path(FIGURES, "fig3_envelope.pdf"), g, width = 9, height = 4.2)
}

## ================= paired gains over two-stage, three scorings =========
TS <- bind_rows(
  if (!is.null(E1)) two_stage_gains(E1, c("n", "p", "mu_sep", "beta_sep"), "E1"),
  if (!is.null(E2)) two_stage_gains(E2, c("n", "p", "mu_sep", "beta_sep"), "E2"),
  if (!is.null(E3a)) two_stage_gains(E3a, c("n", "p", "mu_sep", "beta_sep"), "E3a"),
  if (!is.null(E3b)) two_stage_gains(E3b, c("n", "p", "mu_sep", "beta_sep"), "E3b"))
if (nrow(TS)) {
  wcsv(TS, "tab_two_stage_scorings.csv")
  cat("\n=== paired gain of GeM-Cox (gamma=1) over two-stage: degenerate fits excluded (primary),\n")
  cat("    scored at the floor (sensitivity), or scored as zero coefficients (as run) ===\n")
  print(as.data.frame(TS %>%
    transmute(experiment, n, p, mu_sep, beta_sep, scoring,
              est = sprintf("%.3f (%.3f) n=%d", diff, diff_mcse, pairs)) %>%
    pivot_wider(names_from = scoring, values_from = est)), row.names = FALSE)
}

## ================= E4: theory =========================================
E4a <- read_results("E4a_two_scales.rds")
if (!is.null(E4a)) {
  S <- add_pct(summarise_cells(E4a))
  wcsv(S, "tab_E4a.csv")
  FL <- random_floor(10)
  cat("\n=== E4a two scales ===\n")
  print(as.data.frame(S %>% filter(grepl("gamma=1|Oracle", method)) %>%
        transmute(n, events = round(events), method, recovery = round(recovery, 3),
                  mcse = round(rec_mcse, 3), ARI = round(ari, 3),
                  ari_mcse = round(ari_mcse, 3))), row.names = FALSE)
  D <- S %>% filter(grepl("gamma=1", method))
  g <- ggplot() +
    geom_hline(yintercept = FL, linetype = 2, colour = GRY) +
    geom_line(data = D, aes(events, recovery, colour = "Mechanism contrast"), linewidth = 1.1) +
    geom_pointrange(data = D, aes(events, recovery, colour = "Mechanism contrast",
                    ymin = recovery - rec_mcse, ymax = recovery + rec_mcse), size = .4) +
    geom_line(data = D, aes(events, ari, colour = "Individual labels (ARI)"), linewidth = 1.1) +
    geom_pointrange(data = D, aes(events, ari, colour = "Individual labels (ARI)",
                    ymin = ari - ari_mcse, ymax = ari + ari_mcse), size = .4) +
    scale_colour_manual(values = c("Mechanism contrast" = ORG,
                                   "Individual labels (ARI)" = GRY)) +
    scale_x_log10() + ylim(-0.02, 1) +
    labs(title = "Label recovery and contrast recovery as the sample grows",
         subtitle = paste("GeM-Cox gamma = 1, mu_sep = 0, beta_sep = 2. Bars: +/- 1 Monte Carlo SE.", FLOOR_NOTE),
         x = "Number of events", y = "Recovery") + th
  ggsave(file.path(FIGURES, "fig4a_two_scales.pdf"), g, width = 8, height = 4.6)
}

E4b <- read_results("E4b_mechanism.rds")
if (!is.null(E4b) && "r_mech" %in% names(E4b)) {
  S <- E4b %>% filter(!failed) %>% group_by(method) %>%
    summarise(reps = sum(is.finite(r_mech)), mcse = mcse(r_mech),
              r = mean(r_mech, na.rm = TRUE), .groups = "drop")
  wcsv(S, "tab_E4b.csv")
  cat("\n=== E4b mechanism check: Spearman(fitted log-odds, (x'Delta) x martingale residual) ===\n")
  print(as.data.frame(S %>% mutate(r = round(r, 3), mcse = round(mcse, 3))), row.names = FALSE)
}

E4c <- read_results("E4c_softness.rds")
if (!is.null(E4c)) {
  ## Degeneracy rule applied as everywhere else; configurations whose fits
  ## are all degenerate are reported separately. No correlation across
  ## configurations: per-configuration recovery +/- MCSE only.
  CL <- E4c %>% filter(!failed, !grepl("Oracle", method)) %>% group_by(method) %>%
    summarise(fits = n(), degenerate = sum(!is.na(degenerate) & degenerate),
              min_size_median = stats::median(min_size),
              sharp_mcse = mcse(sharpness), sharpness_all = mean(sharpness),
              .groups = "drop")
  S <- summarise_cells(E4c, by = "method") %>% filter(!grepl("Oracle", method))
  wcsv(S, "tab_E4c.csv"); wcsv(CL, "tab_E4c_collapse.csv")
  cat("\n=== E4c posterior softness: recovery +/- MCSE by configuration (degenerate fits excluded) ===\n")
  print(as.data.frame(S %>% arrange(sharpness) %>%
    transmute(method, fits = reps, mean_max_tau = round(sharpness, 3),
              sharp_mcse = round(sharp_mcse, 3), recovery = round(recovery, 3),
              mcse = round(rec_mcse, 3))), row.names = FALSE)
  collapsed <- CL %>% filter(degenerate == fits)
  if (nrow(collapsed)) {
    cat("\n  Configurations that collapsed (every fit degenerate; reported separately):\n")
    print(as.data.frame(collapsed %>% transmute(method, fits, degenerate,
          median_smallest_cluster = min_size_median, mean_max_tau = round(sharpness_all, 3))),
          row.names = FALSE)
  }
  g <- ggplot(S %>% mutate(method = stats::reorder(method, sharpness)),
              aes(method, recovery)) +
    geom_pointrange(aes(ymin = recovery - rec_mcse, ymax = recovery + rec_mcse),
                    colour = ORG, size = .5) +
    geom_text(aes(label = sprintf("max tau %.2f", sharpness)), vjust = -1.2, size = 3,
              colour = INK) +
    coord_flip() +
    labs(title = "Contrast recovery by E-step configuration",
         subtitle = paste0("n = 800, p = 10, mu_sep = 0, beta_sep = 2. Bars: +/- 1 MCSE. ",
                           "Degenerate fits excluded",
                           if (nrow(collapsed)) paste0("; collapsed: ",
                             paste(collapsed$method, collapse = ", ")) else "", "."),
         x = NULL, y = "Recovery of the contrast") + th
  ggsave(file.path(FIGURES, "fig4c_softness.pdf"), g, width = 8, height = 4.5)
}

## ================= E5a: test calibration and power ====================
E5a <- read_results("E5a_calibration.rds")
if (!is.null(E5a)) {
  L <- E5a %>% filter(!failed) %>%
    pivot_longer(c(p_bootstrap, p_permutation), names_to = "null", values_to = "pval",
                 names_prefix = "p_")
  S <- L %>% group_by(role, n, rho, mu_sep, beta_sep, null) %>%
    summarise(datasets = sum(is.finite(pval)), rejections = sum(pval <= 0.05, na.rm = TRUE),
              mean_p = mean(pval, na.rm = TRUE), .groups = "drop") %>%
    mutate(rate = rejections / datasets, mcse = binom_mcse(rejections, datasets))
  F5 <- E5a %>% group_by(role, n, rho, mu_sep, beta_sep) %>%
    summarise(failed = sum(failed), events = mean(events), secs = mean(secs), .groups = "drop")
  wcsv(S, "tab_E5a.csv"); wcsv(F5, "tab_E5a_failures.csv")
  side <- S %>% select(role, n, rho, mu_sep, beta_sep, datasets, null, rate, mcse) %>%
    pivot_wider(names_from = null, values_from = c(rate, mcse))
  cat("\n=== E5a rejection rate at 0.05 (B = 19, so the exact null rate is 1/20 = 0.05) ===\n")
  cat("  Type I rows: K = 1 with a nonzero shared beta. Read them first.\n")
  print(as.data.frame(side %>% arrange(desc(role == "type I"), rho, n, mu_sep, beta_sep) %>%
        mutate(across(starts_with(c("rate", "mcse")), ~round(.x, 3)))), row.names = FALSE)
  print(as.data.frame(F5 %>% mutate(events = round(events), secs = round(secs, 1))),
        row.names = FALSE)
  H <- L %>% filter(role == "type I") %>%
    mutate(cell = sprintf("K = 1, rho = %g (%s null)", rho, null))
  g <- ggplot(H, aes(pval)) +
    geom_histogram(breaks = seq(0, 1, by = 0.05), fill = BLU, colour = "white") +
    geom_hline(data = H %>% count(cell) %>% mutate(expected = n / 20),
               aes(yintercept = expected), linetype = 2, colour = INK) +
    facet_wrap(~ cell, scales = "free_y") +
    labs(title = "Null p-value distributions of the K = 2 vs K = 1 test",
         subtitle = "Dashed line: count expected under uniformity (B = 19, p-values on a 1/20 grid).",
         x = "p-value", y = "Datasets") + th
  ggsave(file.path(FIGURES, "fig5a_null_pvalues.pdf"), g, width = 9, height = 6)
  P <- S %>% filter(role == "power")
  g <- ggplot(P, aes(factor(beta_sep), rate, colour = null, group = null)) +
    geom_hline(yintercept = 0.05, linetype = 2, colour = GRY) +
    geom_line(linewidth = 1) +
    geom_pointrange(aes(ymin = rate - mcse, ymax = rate + mcse), size = .4) +
    facet_grid(paste("mu_sep =", mu_sep) ~ paste("n =", n)) +
    scale_colour_manual(values = c(bootstrap = ORG, permutation = BLU)) + ylim(0, 1) +
    labs(title = "Power of the K = 2 vs K = 1 test, both null distributions",
         subtitle = "Rejection at 0.05; p = 10. Bars: +/- 1 Monte Carlo SE.",
         x = "beta_sep", y = "Rejection rate") + th
  ggsave(file.path(FIGURES, "fig5a_power.pdf"), g, width = 8, height = 6)
}

## ================= E5c alongside E5a: CV statistic vs in-sample LRT =====
E5c <- read_results("E5c_lrt.rds")
if (!is.null(E5c)) {
  ## label, not `stat`: the results have a column called `stat` (the observed
  ## statistic), which data masking would pick up instead of the argument
  tidy_tests <- function(R, label) {
    R %>% filter(!failed) %>%
      pivot_longer(c(p_bootstrap, p_permutation), names_to = "null", values_to = "pval",
                   names_prefix = "p_") %>%
      mutate(statistic = .env$label)
  }
  L <- bind_rows(if (!is.null(E5a)) tidy_tests(E5a, "CV"), tidy_tests(E5c, "LRT"))
  S <- L %>% group_by(role, n, rho, mu_sep, beta_sep, statistic, null) %>%
    summarise(datasets = sum(is.finite(pval)), rejections = sum(pval <= 0.05, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(rate = rejections / datasets, mcse = binom_mcse(rejections, datasets))
  wcsv(S, "tab_E5ac.csv")
  cat("\n=== E5a (CV statistic) and E5c (in-sample LRT) side by side: rejection at 0.05 (MCSE) ===\n")
  cat("  Same datasets for both statistics. Type I rows: K = 1 with a nonzero shared beta.\n")
  print(as.data.frame(S %>%
    transmute(role, n, rho, mu_sep, beta_sep, test = paste(statistic, null, sep = "-"),
              est = sprintf("%.3f (%.3f)", rate, mcse)) %>%
    pivot_wider(names_from = test, values_from = est) %>%
    arrange(desc(role == "type I"), rho, n, mu_sep, beta_sep)), row.names = FALSE)
  KS <- L %>% filter(role == "type I") %>% group_by(statistic, null, rho) %>%
    summarise(ks_p = suppressWarnings(ks.test(pval, "punif")$p.value),
              mean_p = mean(pval), .groups = "drop")
  cat("  null p-value uniformity (type I cells):\n")
  print(as.data.frame(KS %>% mutate(ks_p = signif(ks_p, 2), mean_p = round(mean_p, 3))),
        row.names = FALSE)
  F5c <- E5c %>% group_by(role, n, rho, mu_sep, beta_sep) %>%
    summarise(failed = sum(failed), secs = mean(secs), .groups = "drop")
  cat(sprintf("  E5c failed datasets: %d of %d\n", sum(F5c$failed), nrow(E5c)))
  H <- L %>% filter(role == "type I", statistic == "LRT") %>%
    mutate(cell = sprintf("K = 1, rho = %g (LRT, %s null)", rho, null))
  g <- ggplot(H, aes(pval)) +
    geom_histogram(breaks = seq(0, 1, by = 0.05), fill = BLU, colour = "white") +
    geom_hline(data = H %>% count(cell) %>% mutate(expected = n / 20),
               aes(yintercept = expected), linetype = 2, colour = INK) +
    facet_wrap(~ cell, scales = "free_y") +
    labs(title = "Null p-value distributions of the in-sample LRT",
         subtitle = "Dashed line: count expected under uniformity (B = 19).",
         x = "p-value", y = "Datasets") + th
  ggsave(file.path(FIGURES, "fig5c_lrt_null_pvalues.pdf"), g, width = 9, height = 6)
  P <- S %>% filter(role == "power")
  g <- ggplot(P, aes(factor(beta_sep), rate, colour = statistic, linetype = null,
                     group = interaction(statistic, null))) +
    geom_hline(yintercept = 0.05, linetype = 2, colour = GRY) +
    geom_line(linewidth = 0.9) +
    geom_pointrange(aes(ymin = rate - mcse, ymax = rate + mcse), size = .3,
                    position = position_dodge(width = 0.25)) +
    facet_grid(paste("mu_sep =", mu_sep) ~ paste("n =", n)) +
    scale_colour_manual(values = c(CV = BLU, LRT = ORG)) + ylim(0, 1) +
    labs(title = "Power of the K = 2 vs K = 1 tests: CV statistic and in-sample LRT",
         subtitle = "Rejection at 0.05; p = 10; same datasets. Bars: +/- 1 Monte Carlo SE.",
         x = "beta_sep", y = "Rejection rate") + th
  ggsave(file.path(FIGURES, "fig5c_power_cv_vs_lrt.pdf"), g, width = 8.5, height = 6)
}

## ================= E5b: K selection by regime =========================
E5b <- read_results("E5b_kselect.rds")
if (!is.null(E5b)) {
  S <- E5b %>% filter(!failed) %>% group_by(regime, mu_sep, beta_sep) %>%
    summarise(datasets = n(), k2_joint = sum(K_joint == 2), k2_partial = sum(K_partial == 2),
              .groups = "drop") %>%
    mutate(P_K2_joint = k2_joint / datasets, mcse_joint = binom_mcse(k2_joint, datasets),
           P_K2_partial = k2_partial / datasets, mcse_partial = binom_mcse(k2_partial, datasets))
  wcsv(S, "tab_E5b.csv")
  cat("\n=== E5b K selection: proportion choosing K = 2 ===\n")
  print(as.data.frame(S %>% select(regime, mu_sep, beta_sep, datasets, P_K2_joint,
        mcse_joint, P_K2_partial, mcse_partial) %>%
        mutate(across(where(is.double), ~round(.x, 3)))), row.names = FALSE)
  cat(sprintf("  failed datasets: %d\n", sum(E5b$failed)))
}


## ================= E2b: adaptive normalisation (final experiment) =====
E2ba <- read_results("E2b_a_adaptive.rds")
if (!is.null(E2ba)) {
  by <- c("n", "p", "mu_sep", "beta_sep")
  ARM <- "GeM-Cox (gamma=1, adaptive)"; AUX <- "GeM-Cox (gamma=1, normalize=TRUE)"
  P <- bind_rows(lapply(c(ARM, AUX), function(a) paired_diffs(
    E2ba, by, target = a, comparators = "GeM-Cox (gamma=1)") %>% mutate(arm = a)))
  NG <- bind_rows(lapply(c("GeM-Cox (gamma=1)", ARM, AUX), function(a)
    normalised_gain(E2ba, by, arm = a) %>% mutate(arm = a)))
  CH <- E2ba %>% filter(!failed, method == ARM) %>% group_by(mu_sep, beta_sep) %>%
    summarise(chose_TRUE = mean(adaptive_normalize), .groups = "drop")
  wcsv(P, "tab_E2b_a_paired.csv"); wcsv(NG, "tab_E2b_a_normalised.csv")
  cat("\n=== E2b(a): adaptive rule vs normalize = FALSE on E2's datasets (paired, degenerate excluded) ===\n")
  print(as.data.frame(P %>% select(arm, mu_sep, beta_sep, pairs, diff, diff_mcse) %>%
    mutate(est = sprintf("%.3f (%.3f)", diff, diff_mcse)) %>%
    select(-diff, -diff_mcse, -pairs) %>%
    pivot_wider(names_from = arm, values_from = est) %>% left_join(CH, by = c("mu_sep", "beta_sep")) %>%
    mutate(chose_TRUE = round(chose_TRUE, 3))), row.names = FALSE)
  cat("\n  normalised gain (arm - gamma0) / (oracle - gamma0):\n")
  print(as.data.frame(NG %>% transmute(arm, mu_sep, beta_sep,
        est = sprintf("%.3f (%.3f)", norm_gain, norm_mcse)) %>%
        pivot_wider(names_from = arm, values_from = est)), row.names = FALSE)
  ## pre-registered decision rule
  A <- P %>% filter(arm == ARM)
  worse <- A %>% filter(diff < -2 * diff_mcse)
  avg <- mean(A$diff); avg_mcse <- sqrt(sum(A$diff_mcse^2)) / nrow(A)
  adopt <- nrow(worse) == 0 && avg > 2 * avg_mcse
  cat(sprintf(paste0("\n  DECISION RULE (pre-registered): (i) cells where adaptive is worse than FALSE by > 2 MCSE: %d;\n",
                     "  (ii) mean paired difference over %d cells = %.4f (MCSE %.4f; threshold 2 MCSE = %.4f).\n",
                     "  DECISION: %s\n"),
              nrow(worse), nrow(A), avg, avg_mcse, 2 * avg_mcse,
              if (adopt) "the adaptive rule becomes the default." else
                "normalize_gmm_by_dim = FALSE stays the default."))
  if (nrow(worse)) print(as.data.frame(worse %>% select(mu_sep, beta_sep, diff, diff_mcse)),
                         row.names = FALSE)
  wcsv(data.frame(cells_worse = nrow(worse), mean_diff = avg, mean_diff_mcse = avg_mcse,
                  adopt_adaptive = adopt), "tab_E2b_decision.csv")
  key <- c("cell", "rep", "method")
  if (!is.null(E2)) {
    chk <- inner_join(E2 %>% select(all_of(key), recovery),
                      E2ba %>% filter(method == "GeM-Cox (gamma=1)") %>% select(all_of(key), recovery),
                      by = key)
    cat(sprintf("  consistency: GeM-Cox (gamma=1) refitted in E2b reproduces E2 exactly: %s (%d fits)\n",
                identical(chk$recovery.x, chk$recovery.y), nrow(chk)))
  }
}

test_rates <- function(R, cols) {
  R %>% filter(!failed) %>% group_by(role, n, p, rho, mu_sep, beta_sep) %>%
    summarise(across(all_of(cols), list(
      rate = ~ mean(.x <= 0.05, na.rm = TRUE),
      mcse = ~ binom_mcse(sum(.x <= 0.05, na.rm = TRUE), sum(is.finite(.x)))),
      .names = "{.col}_{.fn}"),
      datasets = n(), chose_TRUE = mean(adaptive_normalize, na.rm = TRUE), .groups = "drop")
}
E2bb <- read_results("E2b_b_lrt_adaptive.rds")
if (!is.null(E2bb)) {
  S <- test_rates(E2bb, "p_adaptive")
  if (!is.null(E5c)) {
    D <- E5c %>% filter(!failed) %>% group_by(role, n, rho, mu_sep, beta_sep) %>%
      summarise(E5c_default_rate = mean(p_bootstrap <= 0.05),
                E5c_default_mcse = binom_mcse(sum(p_bootstrap <= 0.05), n()), .groups = "drop")
    S <- S %>% left_join(D, by = c("role", "n", "rho", "mu_sep", "beta_sep"))
  }
  wcsv(S, "tab_E2b_b.csv")
  cat("\n=== E2b(b): LRT (bootstrap null) under the adaptive rule vs under FALSE (E5c), same datasets ===\n")
  print(as.data.frame(S %>% transmute(role, n, rho, mu_sep, beta_sep, datasets,
        chose_TRUE = round(chose_TRUE, 3),
        adaptive = sprintf("%.3f (%.3f)", p_adaptive_rate, p_adaptive_mcse),
        default_FALSE = if ("E5c_default_rate" %in% names(S))
          sprintf("%.3f (%.3f)", E5c_default_rate, E5c_default_mcse) else NA) %>%
        arrange(desc(role == "type I"), rho, n, mu_sep, beta_sep)), row.names = FALSE)
  KS <- E2bb %>% filter(!failed, role == "type I") %>% group_by(rho) %>%
    summarise(ks_p = suppressWarnings(ks.test(p_adaptive, "punif")$p.value), .groups = "drop")
  cat("  null p-value uniformity (adaptive):", paste(sprintf("rho=%g KS p=%.2g", KS$rho, KS$ks_p),
                                                   collapse = "; "), "\n")
  cat(sprintf("  failed datasets: %d of %d\n", sum(E2bb$failed), nrow(E2bb)))
}

E2bc <- read_results("E2b_c_cvia078.rds")
if (!is.null(E2bc)) {
  S <- test_rates(E2bc, c("p_false", "p_adaptive"))
  EV <- E2bc %>% group_by(role, mu_sep, beta_sep) %>%
    summarise(events = mean(events), failed = sum(failed), .groups = "drop")
  wcsv(S, "tab_E2b_c.csv")
  cat("\n=== E2b(c): CVIA078 scale (n = 117, event rate 0.44, p = 9): LRT, bootstrap null ===\n")
  print(as.data.frame(S %>% left_join(EV, by = c("role", "mu_sep", "beta_sep")) %>%
        transmute(role, mu_sep, beta_sep, datasets, events = round(events, 1), failed,
                  chose_TRUE = round(chose_TRUE, 3),
                  FALSE_ = sprintf("%.3f (%.3f)", p_false_rate, p_false_mcse),
                  adaptive = sprintf("%.3f (%.3f)", p_adaptive_rate, p_adaptive_mcse))),
        row.names = FALSE)
}

## ================= ADDED AFTER sim-freeze-v1: E5d ======================
E5d <- read_results("E5d_profile_null.rds")
if (!is.null(E5d)) {
  S <- E5d %>% filter(!failed) %>% group_by(role, n, rho, mu_sep, beta_sep) %>%
    summarise(datasets = sum(is.finite(p_bootstrap)),
              rejections = sum(p_bootstrap <= 0.05, na.rm = TRUE),
              ks_p = suppressWarnings(ks.test(p_bootstrap, "punif")$p.value),
              mean_p = mean(p_bootstrap, na.rm = TRUE), .groups = "drop") %>%
    mutate(rate = rejections / datasets, mcse = binom_mcse(rejections, datasets),
           threshold = 0.05 + 2 * sqrt(0.05 * 0.95 / datasets),
           verdict = ifelse(rate > threshold, "FAILS (exceeds 0.05 + 2 MCSE)", "calibrated"))
  REF <- NULL
  if (!is.null(E5c)) {
    REF <- E5c %>% filter(!failed, role == "type I") %>% group_by(role, n, rho, mu_sep, beta_sep) %>%
      summarise(datasets = n(), rejections = sum(p_bootstrap <= 0.05),
                ks_p = suppressWarnings(ks.test(p_bootstrap, "punif")$p.value),
                mean_p = mean(p_bootstrap), .groups = "drop") %>%
      mutate(role = "type I (E5c, mu_sep = 0)", rate = rejections / datasets,
             mcse = binom_mcse(rejections, datasets),
             threshold = 0.05 + 2 * sqrt(0.05 * 0.95 / datasets), verdict = "(reference)")
  }
  T5 <- bind_rows(REF, S)
  wcsv(T5, "tab_E5d.csv")
  cat("\n=== E5d (added after sim-freeze-v1): profile-only null, LRT with bootstrap null, B = 19 ===\n")
  print(as.data.frame(T5 %>% transmute(role, n, rho, mu_sep, datasets,
        rate = sprintf("%.3f (%.3f)", rate, mcse), threshold = round(threshold, 4),
        ks_p = signif(ks_p, 2), mean_p = round(mean_p, 3), verdict)), row.names = FALSE)
  cat(sprintf("  E5d VERDICT: %s\n", if (any(S$rate > S$threshold))
    "the test is NOT valid for data with profile structure (a cell exceeds 0.05 + 2 MCSE)." else
    "calibrated in every profile-only cell (no cell exceeds 0.05 + 2 MCSE)."))
  H <- E5d %>% filter(!failed) %>% mutate(cell = sprintf("mu_sep = %g, rho = %g", mu_sep, rho))
  g <- ggplot(H, aes(p_bootstrap)) +
    geom_histogram(breaks = seq(0, 1, by = 0.05), fill = BLU, colour = "white") +
    geom_hline(data = H %>% count(cell) %>% mutate(expected = n / 20),
               aes(yintercept = expected), linetype = 2, colour = INK) +
    facet_wrap(~ cell, scales = "free_y") +
    labs(title = "E5d: null p-values of the LRT when profiles separate but mechanisms do not",
         subtitle = "Added after sim-freeze-v1. Dashed line: count expected under uniformity (B = 19).",
         x = "p-value", y = "Datasets") + th
  ggsave(file.path(FIGURES, "fig5d_profile_null_pvalues.pdf"), g, width = 9, height = 6)
}

## ================= ADDED AFTER sim-freeze-v1: E5e ======================
E5e <- read_results("E5e_kselect.rds")
if (!is.null(E5e)) {
  SEL <- c(joint = "K_joint", partial = "K_partial", BIC = "K_bic", `sequential LRT` = "K_seq")
  cell_by <- c("set", "K_true", "regime", "separation", "n")
  L <- bind_rows(lapply(names(SEL), function(nm) {
    E5e %>% transmute(across(all_of(cell_by)), selector = nm, K_hat = .data[[SEL[[nm]]]])
  })) %>% filter(!is.na(K_hat))
  CM <- L %>% group_by(across(all_of(c(cell_by, "selector")))) %>%
    summarise(datasets = n(), k1 = sum(K_hat == 1), k2 = sum(K_hat == 2),
              k3 = sum(K_hat == 3), k4 = sum(K_hat == 4), correct = sum(K_hat == K_true),
              .groups = "drop") %>%
    pivot_longer(c(k1, k2, k3, k4, correct), names_to = "k", values_to = "count") %>%
    mutate(prob = count / datasets, mcse = binom_mcse(count, datasets))
  wcsv(CM, "tab_E5e_confusion.csv")
  cat("\n=== E5e (added after sim-freeze-v1): selecting K, grid 1:4. P(K_hat = k | cell) (MCSE) ===\n")
  for (sel in names(SEL)) {
    cat(sprintf("\n  -- selector: %s --\n", sel))
    W <- CM %>% filter(selector == sel) %>%
      mutate(est = sprintf("%.2f (%.2f)", prob, mcse)) %>%
      select(all_of(cell_by), datasets, k, est) %>%
      pivot_wider(names_from = k, values_from = est) %>%
      arrange(desc(set == "main"), K_true, regime, separation, n)
    print(as.data.frame(W %>% rename(`K=1` = k1, `K=2` = k2, `K=3` = k3, `K=4` = k4)),
          row.names = FALSE)
  }
  ## which question each selector answers: P(K_hat >= 2) by regime
  Q <- L %>% filter(set == "main") %>%
    mutate(group = ifelse(K_true == 1, "K_true = 1 (null)", paste0(regime, ", K_true >= 2"))) %>%
    group_by(selector, group) %>%
    summarise(datasets = n(), more_than_one = sum(K_hat >= 2), .groups = "drop") %>%
    mutate(est = sprintf("%.2f (%.3f)", more_than_one / datasets,
                         binom_mcse(more_than_one, datasets))) %>%
    select(selector, group, est) %>% pivot_wider(names_from = group, values_from = est)
  wcsv(Q, "tab_E5e_question.csv")
  cat("\n  P(K_hat >= 2) by regime (main cells): which question each selector answers\n")
  print(as.data.frame(Q), row.names = FALSE)
  S1R <- E5e %>% filter(!is.na(p_2v1)) %>% group_by(across(all_of(cell_by))) %>%
    summarise(datasets = n(), rej = sum(p_2v1 <= 0.05), .groups = "drop") %>%
    mutate(est = sprintf("%.3f (%.3f)", rej / datasets, binom_mcse(rej, datasets)))
  wcsv(S1R, "tab_E5e_seq_step1.csv")
  cat("\n  sequential LRT, first step (K = 2 vs 1): rejection rate at 0.05\n")
  print(as.data.frame(S1R %>% select(all_of(cell_by), datasets, est) %>%
        arrange(desc(set == "main"), K_true, regime, separation, n)), row.names = FALSE)
  cat(sprintf("  failures: CV criterion %d, sequential LRT step %d, BIC fits %d (of %d datasets)\n",
              sum(E5e$cv_failed), sum(E5e$seq_failed, na.rm = TRUE),
              sum(is.na(E5e$K_bic)), nrow(E5e)))
  F <- CM %>% filter(k != "correct") %>%
    mutate(K_hat = as.integer(sub("k", "", k)),
           cell = sprintf("%s K=%d %s %s n=%d", ifelse(set == "main", "", "CVIA"), K_true,
                          regime, separation, n))
  g <- ggplot(F, aes(factor(K_hat), cell, fill = prob)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", prob)), size = 2.2) +
    facet_wrap(~ selector, nrow = 1) +
    scale_fill_gradient(low = "white", high = BLU, limits = c(0, 1)) +
    labs(title = "E5e: P(selected K | true design)",
         subtitle = "Added after sim-freeze-v1. K grid 1:4; 200 datasets per cell.",
         x = "Selected K", y = NULL, fill = "P") + th +
    theme(axis.text.y = element_text(size = 6), legend.position = "right")
  ggsave(file.path(FIGURES, "fig5e_confusion.pdf"), g, width = 13, height = 8)
}

## ================= ADDED AFTER sim-freeze-v1: E5f (optional) ===========
E5f <- read_results("E5f_adaptive_scale.rds")
if (!is.null(E5f)) {
  S <- test_rates(E5f, c("p_false", "p_adaptive"))
  if (!is.null(E2bc)) S <- bind_rows(test_rates(E2bc, c("p_false", "p_adaptive")), S)
  KS <- bind_rows(E2bc, E5f) %>% filter(!failed, role == "type I") %>% group_by(n) %>%
    summarise(ks_false = suppressWarnings(ks.test(p_false, "punif")$p.value),
              ks_adaptive = suppressWarnings(ks.test(p_adaptive, "punif")$p.value),
              .groups = "drop")
  wcsv(S, "tab_E5f.csv")
  cat("\n=== E5f (added after sim-freeze-v1): adaptive vs default LRT (bootstrap, B = 19), event rate 0.44 ===\n")
  cat("  n = 117 rows are E2b(c) (p = 9), for context; n = 200, 300 are E5f (p = 10). Same datasets within a row.\n")
  print(as.data.frame(S %>% arrange(n, desc(role == "type I")) %>%
        transmute(role, n, p, beta_sep, datasets, chose_TRUE = round(chose_TRUE, 3),
                  default_FALSE = sprintf("%.3f (%.3f)", p_false_rate, p_false_mcse),
                  adaptive = sprintf("%.3f (%.3f)", p_adaptive_rate, p_adaptive_mcse))),
        row.names = FALSE)
  cat("  null p-value uniformity (KS p):\n")
  print(as.data.frame(KS %>% mutate(across(starts_with("ks"), ~signif(.x, 2)))), row.names = FALSE)
  cat(sprintf("  failed datasets: %d of %d\n", sum(E5f$failed), nrow(E5f)))
}

## ================= S1: normalisation ==================================
S1 <- read_results("S1_normalisation.rds")
if (!is.null(S1)) {
  by <- c("n", "p", "mu_sep", "beta_sep")
  P <- paired_diffs(S1, by, target = "GeM-Cox, normalize = FALSE",
                    comparators = "GeM-Cox, normalize = TRUE")
  wcsv(P, "tab_S1_paired.csv")
  S <- add_pct(summarise_cells(S1))
  wcsv(S, "tab_S1_recovery.csv")
  cat("\n=== S1 normalisation: paired difference, normalize = FALSE minus TRUE ===\n")
  print_paired(P, c("p", "mu_sep"))
  cat("\n  descriptive recovery:\n")
  print(as.data.frame(S %>% transmute(p, mu_sep, method, recovery = round(recovery, 3),
        mcse = round(rec_mcse, 3), ARI = round(ari, 3))), row.names = FALSE)
  cat("  Recorded expectation: normalize = FALSE better whenever mu_sep > 0;\n")
  cat("  normalize = TRUE better only at mu_sep = 0 with large p.\n")
}

## ================= POST HOC diagnostic (not pre-specified) ============
## Added after the pilot. Does the two-stage contrast align with the shared
## coefficient only when mechanisms differ? beta_sep = 0 is the diagnostic;
## beta_sep = 2 is a heterogeneous reference. See 04_posthoc.R.
PH <- read_results("PH1_alignment.rds")
if (!is.null(PH)) {
  S <- PH %>% filter(!failed) %>% group_by(beta_sep, method) %>%
    summarise(reps = n(),
              norm_mcse = mcse(contrast_norm), cos_mcse = mcse(cos_base),
              contrast_norm = mean(contrast_norm, na.rm = TRUE),
              cos_base = mean(cos_base, na.rm = TRUE), .groups = "drop")
  wcsv(S, "tab_PH1_alignment.csv")
  cat("\n=== POST HOC (not pre-specified): alignment of the fitted contrast with the shared coefficient ===\n")
  cat("  mu_sep = 0, n = 400, p = 10. cos_base = |cos(beta1_hat - beta2_hat, shared coefficient)|.\n")
  print(as.data.frame(S %>% mutate(across(c(contrast_norm, norm_mcse, cos_base, cos_mcse),
                                          ~round(.x, 3)))), row.names = FALSE)
  cat(sprintf("  isotropic reference for |cos| with a fixed direction in R^10: %.3f (descriptive)\n",
              random_floor(10)))
  cat(sprintf("  failed fits: %d of %d\n", sum(PH$failed), nrow(PH)))
}

cat("\nFigures written to", FIGURES, "and tables to", RESULTS, "\n")
