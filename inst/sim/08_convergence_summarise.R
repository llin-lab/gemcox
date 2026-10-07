###############################################################
## 08_convergence_summarise.R  --  which frozen claims hold for converged
## estimators (R1, R2), and the matched-regularisation comparison (a
## sensitivity analysis). Added after sim-freeze-v2; implements the analysis
## pre-registered in README.md and was written before the runs. Reads only
## results/; nothing is refitted.
##
## Every new row is checked against C0_dataset_hashes.rds (E1-E3b, C1, C2)
## or R0_dataset_hashes.rds (E4a, E5c), or against C5's primary rows (C5,
## which has no frozen datasets), before it is compared with anything.
##
## Each frozen claim is evaluated twice with the SAME code and criterion:
## on the frozen fits (tolerance 1e-4, max_iter 100) and on the converged
## fits (tolerance 1e-8, max_iter 1000).
##
## ASCII only.
###############################################################

suppressMessages(library(ggplot2))
Sys.setenv(GEMCOX_WORKERS = "1")
source("00_core.R")
source("C_common.R")
source("C_methods.R")

import_03 <- function(nms) {
  for (e in parse("03_summarise.R", keep.source = FALSE)) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
        as.character(e[[2]]) %in% nms) eval(e, globalenv())
  }
  miss <- nms[!vapply(nms, exists, NA, envir = globalenv(), inherits = FALSE)]
  if (length(miss)) stop("not found in 03_summarise.R: ", paste(miss, collapse = ", "))
}
import_03(c("BLU", "ORG", "GRY", "SKY", "YEL", "INK", "PAL", "th", "wcsv", "binom_mcse",
            "is_degenerate", "score_degenerate", "summarise_cells", "wide_recovery",
            "paired_diffs", "normalised_gain"))

BY <- c("n", "p", "mu_sep", "beta_sep")
G1F <- "GeM-Cox (gamma=1)"; G0F <- "GeM-Cox (gamma=0)"; TS <- "Two-stage (GMM+Cox)"
SC <- "Single Cox model"; OR <- "Oracle (true labels)"
GRN <- "#009E73"
verdict_cell <- function(diff, se) dplyr::case_when(
  !is.finite(diff) | !is.finite(se) ~ NA_character_,
  diff > 2 * se ~ "beats", diff < -2 * se ~ "loses", TRUE ~ "matches")
verdict_regime <- function(v) {
  v <- v[!is.na(v)]
  if (!length(v)) return(NA_character_)
  b <- any(v == "beats"); l <- any(v == "loses")
  if (b && !l) "beats" else if (l && !b) "loses" else if (!b && !l) "matches" else "mixed"
}
pd <- function(R, target, comparators, contrast, mode = "exclude") {
  if (!target %in% R$method || !any(comparators %in% R$method)) return(NULL)
  paired_diffs(R, BY, target = target, comparators = comparators, mode = mode) %>%
    filter(pairs > 0) %>% mutate(target = target, contrast = contrast)
}
paired_metric <- function(R, target, comparator, metric) {
  score_degenerate(R, "exclude") %>% filter(!failed, method %in% c(target, comparator)) %>%
    select(all_of(c(BY, "rep", "method", metric))) %>%
    tidyr::pivot_wider(names_from = method, values_from = all_of(metric)) %>%
    mutate(d = .data[[target]] - .data[[comparator]]) %>%
    group_by(across(all_of(BY))) %>%
    summarise(pairs = sum(is.finite(d)), diff_mcse = mcse(d), diff = mean(d, na.rm = TRUE),
              .groups = "drop")
}
FIT_DIAG <- function(R, label) R %>% group_by(across(all_of(BY)), method) %>%
  summarise(experiment = label, fits = n(), failed = sum(failed),
            degenerate = sum(!failed & !is.na(degenerate) & degenerate),
            cap_hits = sum(!failed & !is.na(converged) & !converged),
            iters_median = median(iterations[!failed]), .groups = "drop")

## ---- verification ------------------------------------------------------------
MAN <- bind_rows(readRDS(file.path(RESULTS, "C0_dataset_hashes.rds"))$rows,
                 readRDS(file.path(RESULTS, "R0_dataset_hashes.rds"))$rows)
if (!all(MAN$verified)) stop("a dataset manifest contains unverified datasets")
## replicates (cell, rep) failing verification against manifest key mkey, or
## against the hashes of reference rows ref
bad_reps <- function(X, seed_exp, mkey = NULL, ref = NULL) {
  j <- X
  if (!"hash_test" %in% names(j)) j$hash_test <- NA_character_
  j <- j[, c("cell", "rep", "failed", "seed", "events", "hash_train", "hash_test")]
  j$seed_ok <- j$seed == mk_seed(seed_exp, j$cell, j$rep)
  if (!is.null(mkey)) {
    m <- MAN[MAN$experiment == mkey, c("cell", "rep", "seed_frozen", "events_stored",
                                       "hash_train_frozen", "hash_test_frozen")]
    j <- merge(j, m, by = c("cell", "rep"), all.x = TRUE)
    j$good <- with(j, seed_ok & !is.na(seed_frozen) & seed == seed_frozen &
                     events == events_stored &
                     (failed | (hash_train == hash_train_frozen &
                                (is.na(hash_test_frozen) | hash_test == hash_test_frozen))))
  } else {
    r <- unique(ref[!ref$failed, c("cell", "rep", "hash_train", "hash_test")])
    names(r)[3:4] <- c("ref_train", "ref_test")
    j <- merge(j, r, by = c("cell", "rep"), all.x = TRUE)
    j$good <- with(j, seed_ok & (failed | (hash_train == ref_train & hash_test == ref_test)))
  }
  j$good[is.na(j$good)] <- FALSE
  unique(paste(j$cell, j$rep)[!j$good])
}
VER <- list()
note_ver <- function(part, file, X, bad) VER[[length(VER) + 1]] <<- data.frame(
  part = part, file = file, replicates = length(unique(paste(X$cell, X$rep))),
  failing_verification = length(bad))

## ================= R1: E1, E2, E3a, E3b, E4a at convergence =================
R1_SPEC <- list(
  E1  = list(frozen = "E1_existence.rds", key = "C1",  g0 = "R1_E1_convergence.rds"),
  E2  = list(frozen = "E2_comparison.rds", key = "C2", g0 = "R1_E2_convergence.rds"),
  E3a = list(frozen = "E3a_n.rds", key = "C3a", g0 = "R1_E3a_convergence.rds"),
  E3b = list(frozen = "E3b_p.rds", key = "C3b", g0 = "R1_E3b_convergence.rds"),
  E4a = list(frozen = "E4a_two_scales.rds", key = NA, g0 = "R1_E4a_convergence.rds"))
R1 <- list()
for (nm in names(R1_SPEC)) {
  S <- R1_SPEC[[nm]]
  F <- read_results(S$frozen); X <- read_results(S$g0)
  if (is.null(F) || is.null(X)) next
  if (is.na(S$key)) {
    mkey <- "E4a"; seed_exp <- R_EXPS$E4a$seed_exp; tight <- X
    bad <- bad_reps(X, seed_exp, mkey = mkey)
    note_ver("R1", S$g0, X, bad)
  } else {
    E <- C_EXPS[[S$key]]; mkey <- S$key; seed_exp <- E$seed_exp
    cv <- read_results(conv_file(E))
    b1 <- bad_reps(X, seed_exp, mkey = mkey); b2 <- bad_reps(cv, seed_exp, mkey = mkey)
    note_ver("R1", S$g0, X, b1); note_ver("R1 (reused)", conv_file(E), cv, b2)
    bad <- union(b1, b2)
    tight <- bind_rows(X, cv %>% filter(method %in% c(GEM_TIGHT, LCC_TIGHT)))
  }
  R1[[nm]] <- bind_rows(F, tight) %>% filter(!paste(cell, rep) %in% bad) %>%
    mutate(experiment = nm)
}

## ================= R2: E5c at convergence ===================================
F5 <- read_results("E5c_lrt.rds"); N5 <- read_results(R_EXPS$E5c$file)
if (!is.null(N5)) {
  b5 <- bad_reps(N5, 5, mkey = "E5c")
  note_ver("R2", R_EXPS$E5c$file, N5, b5)
  N5 <- N5 %>% filter(!paste(cell, rep) %in% b5)
}

## ================= matched regularisation (sensitivity analysis) =============
MR <- list()
for (key in names(MATCHED_CELLS)) {
  E <- C_EXPS[[key]]
  cv <- read_results(conv_file(E)); mt <- read_results(matched_file(E))
  if (is.null(cv) || is.null(mt)) next
  cv <- cv %>% filter(cell %in% MATCHED_CELLS[[key]])
  if (is.na(E$frozen)) {                     # C5: check against its primary rows
    ref <- read_results(E$file)
    b1 <- bad_reps(cv, E$seed_exp, ref = ref); b2 <- bad_reps(mt, E$seed_exp, ref = ref)
  } else {
    b1 <- bad_reps(cv, E$seed_exp, mkey = key); b2 <- bad_reps(mt, E$seed_exp, mkey = key)
  }
  note_ver("matched", matched_file(E), mt, b2)
  MR[[key]] <- bind_rows(cv, mt) %>% filter(!paste(cell, rep) %in% union(b1, b2)) %>%
    mutate(experiment = key)
}

VER <- do.call(rbind, VER)
wcsv(VER, "tab_R_verification.csv")
cat("\n=== dataset verification ===\n")
print(VER, row.names = FALSE)

## ================= R1 tables =================================================
if (length(R1)) {
  P <- bind_rows(lapply(names(R1), function(nm) {
    R <- R1[[nm]]
    bind_rows(
      pd(R, G1F, G0F, "gamma1 - gamma0: frozen"),
      pd(R, GEM_TIGHT, GEM0_TIGHT, "gamma1 - gamma0: converged"),
      pd(R, G1F, c(TS, SC), "gamma1 - comparator: frozen"),
      pd(R, GEM_TIGHT, c(TS, SC), "gamma1 - comparator: converged"),
      pd(R, GEM_TIGHT, LCC_TIGHT, "gamma1 - competitor: converged"),
      pd(R, LCC_TIGHT, c(GEM0_TIGHT, TS, SC), "competitor - comparator: converged"),
      pd(R, GEM_TIGHT, G1F, "gamma1: converged - frozen"),
      pd(R, GEM0_TIGHT, G0F, "gamma0: converged - frozen")) %>% mutate(experiment = nm)
  })) %>% mutate(verdict = verdict_cell(diff, diff_mcse))
  wcsv(P, "tab_R1_paired.csv")
  show <- function(ct) {
    x <- P %>% filter(contrast == ct)
    if (!nrow(x)) return(invisible())
    cat(sprintf("\n--- %s ---\n", ct))
    print(as.data.frame(x %>% transmute(experiment, n, p, mu_sep, beta_sep, comparator, pairs,
          diff = round(diff, 3), mcse = round(diff_mcse, 3), verdict)), row.names = FALSE)
  }
  cat("\n=== R1: paired differences in recovery (degenerate excluded) ===\n")
  for (ct in unique(P$contrast)) show(ct)

  ## normalised gains, frozen and converged (the frozen function, with the
  ## converged gamma = 0 in the gamma = 0 slot)
  NG <- bind_rows(lapply(names(R1), function(nm) {
    R <- R1[[nm]]
    Rc <- R %>% filter(method != G0F) %>% mutate(method = ifelse(method == GEM0_TIGHT, G0F, method))
    bind_rows(normalised_gain(R, BY, arm = G1F) %>% mutate(fits = "frozen", arm = "GeM-Cox gamma=1"),
              normalised_gain(Rc, BY, arm = GEM_TIGHT) %>% mutate(fits = "converged", arm = "GeM-Cox gamma=1"),
              normalised_gain(Rc, BY, arm = LCC_TIGHT) %>% mutate(fits = "converged", arm = "competitor")) %>%
      mutate(experiment = nm)
  }))
  wcsv(NG, "tab_R1_normalised.csv")
  cat("\n=== R1: normalised gain (arm - gamma0) / (oracle - gamma0) ===\n")
  print(as.data.frame(NG %>% transmute(experiment, n, p, mu_sep, beta_sep, fits, arm, pairs,
        norm_gain = round(norm_gain, 3), mcse = round(norm_mcse, 3)) %>%
        tidyr::pivot_wider(names_from = c(fits, arm), values_from = c(norm_gain, mcse, pairs))),
        row.names = FALSE)

  DS <- bind_rows(lapply(names(R1), function(nm) summarise_cells(R1[[nm]], c(BY, "method")) %>%
                           mutate(experiment = nm)))
  wcsv(DS, "tab_R1_descriptive.csv")
  FD <- bind_rows(lapply(names(R1), function(nm) FIT_DIAG(R1[[nm]] %>%
                           filter(method %in% c(GEM0_TIGHT, GEM_TIGHT, LCC_TIGHT)), nm)))
  wcsv(FD, "tab_R1_fits.csv")
  cat("\n=== R1: converged fits: failures, degenerate fits, cap hits (1000 iterations) ===\n")
  print(as.data.frame(FD %>% group_by(experiment, method) %>%
        summarise(fits = sum(fits), failed = sum(failed), degenerate = sum(degenerate),
                  cap_hits = sum(cap_hits), .groups = "drop")), row.names = FALSE)

  ## ---- the frozen claims, frozen fits vs converged fits ----------------------
  CL <- list()
  claim <- function(id, claim, source, criterion, frozen, converged, holds_f, holds_c)
    CL[[length(CL) + 1]] <<- data.frame(id = id, claim = claim, source = source,
                                        criterion = criterion, frozen = frozen,
                                        converged = converged, holds_frozen = holds_f,
                                        holds_converged = holds_c, stringsAsFactors = FALSE)
  cnt <- function(ct, exps, comp = NULL) {
    x <- P %>% filter(contrast == ct, experiment %in% exps)
    if (!is.null(comp)) x <- x %>% filter(comparator == comp)
    c(beats = sum(x$verdict == "beats", na.rm = TRUE), cells = nrow(x))
  }
  k <- function(v) sprintf("beats in %d of %d cells", v[["beats"]], v[["cells"]])
  a_f <- cnt("gamma1 - gamma0: frozen", c("E1", "E2", "E3a", "E3b"))
  a_c <- cnt("gamma1 - gamma0: converged", c("E1", "E2", "E3a", "E3b"))
  claim("A", "gamma = 1 beats gamma = 0", "E1, E2, E3a, E3b", "paired difference > 2 MCSE in every cell",
        k(a_f), k(a_c), a_f[["beats"]] == a_f[["cells"]], a_c[["beats"]] == a_c[["cells"]])
  for (cm in c(TS, SC)) {
    f <- cnt("gamma1 - comparator: frozen", c("E1", "E2"), cm)
    cc <- cnt("gamma1 - comparator: converged", c("E1", "E2"), cm)
    claim(if (cm == TS) "B" else "C", paste("gamma = 1 beats", cm), "E1, E2",
          "paired difference > 2 MCSE in every cell with pairs", k(f), k(cc),
          f[["beats"]] == f[["cells"]], cc[["beats"]] == cc[["cells"]])
  }
  if (!is.null(R1$E1)) {
    S1 <- summarise_cells(R1$E1, c(BY, "method")) %>% filter(beta_sep >= 1)
    fl <- random_floor(10)
    p1 <- function(g1, blind) {
      m <- S1 %>% group_by(method) %>% summarise(rec = mean(recovery), .groups = "drop")
      ex <- m$rec - fl; names(ex) <- m$method
      c(blind = max(ex[blind]), g1 = ex[[g1]])
    }
    pf <- p1(G1F, c(SC, TS, G0F)); pc <- p1(GEM_TIGHT, c(SC, TS, GEM0_TIGHT))
    fmt1 <- function(v) sprintf("covariate-only max excess %+.3f; gamma = 1 excess %+.3f", v[["blind"]], v[["g1"]])
    claim("D", "P1: covariate-only methods at the floor, gamma = 1 above it", "E1 (beta_sep >= 1)",
          "max covariate-only excess < 0.05 and gamma = 1 excess > 0.05", fmt1(pf), fmt1(pc),
          pf[["blind"]] < 0.05 && pf[["g1"]] > 0.05, pc[["blind"]] < 0.05 && pc[["g1"]] > 0.05)
    p2 <- function(ct) {
      G <- P %>% filter(experiment == "E1", contrast == ct) %>% mutate(scaled = diff / beta_sep^2)
      gb <- G[G$beta_sep >= 1, ]
      list(cv = sd(gb$scaled) / mean(gb$scaled), scaled = gb$scaled)
    }
    qf <- p2("gamma1 - gamma0: frozen"); qc <- p2("gamma1 - gamma0: converged")
    fmt2 <- function(q) sprintf("CV %.2f; gain/beta^2 %s", q$cv, paste(sprintf("%.3f", q$scaled), collapse = ", "))
    claim("E", "P2: gain / beta_sep^2 roughly constant", "E1 (beta_sep >= 1)", "CV < 0.5",
          fmt2(qf), fmt2(qc), is.finite(qf$cv) && qf$cv < 0.5, is.finite(qc$cv) && qc$cv < 0.5)
  }
  if (!is.null(R1$E4a)) {
    S4 <- summarise_cells(R1$E4a, c(BY, "method"))
    p3 <- function(m) {
      x <- S4 %>% filter(method == m) %>% arrange(n)
      list(ari = max(x$ari), exc = min(x$recovery - random_floor(10)),
           rise = x$recovery[x$n == max(x$n)] - x$recovery[x$n == min(x$n)],
           rise_se = sqrt(x$rec_mcse[x$n == max(x$n)]^2 + x$rec_mcse[x$n == min(x$n)]^2),
           rec = x$recovery)
    }
    f <- p3(G1F); cc <- p3(GEM_TIGHT)
    fmt3 <- function(q) sprintf("max ARI %.3f; min excess over floor %+.3f", q$ari, q$exc)
    claim("F", "P3: labels at chance while the contrast is recovered", "E4a (all n)",
          "every cell: ARI < 0.10 and recovery - floor > 0.05", fmt3(f), fmt3(cc),
          f$ari < 0.10 && f$exc > 0.05, cc$ari < 0.10 && cc$exc > 0.05)
    fmt4 <- function(q) sprintf("recovery %s; n 6400 - n 800 = %+.3f (SE %.3f)",
                                paste(sprintf("%.3f", q$rec), collapse = ", "), q$rise, q$rise_se)
    claim("G", "recovery rises with n (two scales)", "E4a", "n = 6400 minus n = 800 > 2 SE",
          fmt4(f), fmt4(cc), f$rise > 2 * f$rise_se, cc$rise > 2 * cc$rise_se)
  }
  gain_trend <- function(nm, v, ct) {
    x <- P %>% filter(experiment == nm, contrast == ct) %>% arrange(.data[[v]])
    d <- x$diff[nrow(x)] - x$diff[1]; se <- sqrt(x$diff_mcse[nrow(x)]^2 + x$diff_mcse[1]^2)
    slope <- if (all(x$diff > 0)) unname(coef(lm(log(x$diff) ~ log(x[[v]])))[2]) else NA_real_
    list(diff = x$diff, d = d, se = se, slope = slope)
  }
  if (!is.null(R1$E3a)) {
    f <- gain_trend("E3a", "n", "gamma1 - gamma0: frozen"); cc <- gain_trend("E3a", "n", "gamma1 - gamma0: converged")
    fmt5 <- function(q) sprintf("gain %s; n 1600 - n 200 = %+.3f (SE %.3f)",
                                paste(sprintf("%.3f", q$diff), collapse = ", "), q$d, q$se)
    claim("H", "gain over gamma = 0 grows with n", "E3a", "n = 1600 minus n = 200 > 2 SE",
          fmt5(f), fmt5(cc), f$d > 2 * f$se, cc$d > 2 * cc$se)
  }
  if (!is.null(R1$E3b)) {
    f <- gain_trend("E3b", "p", "gamma1 - gamma0: frozen"); cc <- gain_trend("E3b", "p", "gamma1 - gamma0: converged")
    fmt6 <- function(q) sprintf("gain %s; log-log slope %s; p 80 - p 10 = %+.3f (SE %.3f)",
                                paste(sprintf("%.3f", q$diff), collapse = ", "),
                                if (is.na(q$slope)) "undefined" else sprintf("%.2f", q$slope), q$d, q$se)
    claim("I", "gain over gamma = 0 decays with p (roughly p^-1.15)", "E3b",
          "p = 80 minus p = 10 < -2 SE (slope reported)", fmt6(f), fmt6(cc),
          f$d < -2 * f$se, cc$d < -2 * cc$se)
  }
  if (!is.null(R1$E2)) {
    mm <- function(g1) paired_metric(R1$E2 %>% filter(mu_sep == 0), g1, SC, "mse_eta")
    f <- mm(G1F); cc <- mm(GEM_TIGHT)
    fmt7 <- function(x) paste(sprintf("%+.3f (%.3f)", x$diff, x$diff_mcse), collapse = "; ")
    claim("J", "prediction: at mu_sep = 0 gamma = 1's MSE is no smaller than single Cox's",
          "E2 (mu_sep = 0)", "no cell with gamma = 1 minus single Cox < -2 MCSE",
          fmt7(f), fmt7(cc), !any(f$diff < -2 * f$diff_mcse), !any(cc$diff < -2 * cc$diff_mcse))
  }
}

## ================= R2 tables =================================================
if (!is.null(N5) && !is.null(F5)) {
  J5 <- merge(F5[, c("cell", "rep", "p_bootstrap", "failed")] %>%
                rename(p_frozen = p_bootstrap, failed_frozen = failed),
              N5[, c("cell", "rep", "role", "n", "rho", "mu_sep", "beta_sep", "p_bootstrap",
                     "failed", "n_fits", "n_not_converged", "valid_bootstrap")],
              by = c("cell", "rep")) %>% filter(!failed, !failed_frozen)
  T5 <- J5 %>% group_by(cell, role, n, rho, mu_sep, beta_sep) %>%
    summarise(datasets = n(), rej_frozen = mean(p_frozen <= 0.05),
              rej_converged = mean(p_bootstrap <= 0.05),
              diff_mcse = sd((p_bootstrap <= 0.05) - (p_frozen <= 0.05)) / sqrt(n()),
              diff = rej_converged - rej_frozen,
              ks_p_converged = suppressWarnings(stats::ks.test(p_bootstrap, "punif")$p.value),
              fits = sum(n_fits), cap_hits = sum(n_not_converged),
              min_valid = min(valid_bootstrap), .groups = "drop") %>%
    mutate(mcse_frozen = binom_mcse(rej_frozen * datasets, datasets),
           mcse_converged = binom_mcse(rej_converged * datasets, datasets),
           threshold = ifelse(role == "type I", 0.05 + 2 * sqrt(0.05 * 0.95 / datasets), NA),
           changed = abs(diff) > 2 * diff_mcse)
  wcsv(T5, "tab_R2_E5c.csv")
  cat("\n=== R2: E5c at convergence: LRT, bootstrap null, rejection at 0.05 (same datasets) ===\n")
  print(as.data.frame(T5 %>% transmute(cell, role, n, rho, mu_sep, beta_sep, datasets,
        frozen = sprintf("%.3f (%.3f)", rej_frozen, mcse_frozen),
        converged = sprintf("%.3f (%.3f)", rej_converged, mcse_converged),
        paired_change = sprintf("%+.3f (%.3f)", diff, diff_mcse),
        threshold = round(threshold, 4), cap_hits = sprintf("%d / %d", cap_hits, fits))),
        row.names = FALSE)
  if (exists("claim")) {
    ti <- T5 %>% filter(role == "type I")
    fmt8 <- function(r) paste(sprintf("%.3f", r), collapse = ", ")
    claim("K", "the LRT with the bootstrap null is calibrated (K = 1 null, nonzero beta)", "E5c type I cells",
          "rejection <= 0.05 + 2 sqrt(0.05 x 0.95 / N) in both cells",
          fmt8(ti$rej_frozen), fmt8(ti$rej_converged),
          all(ti$rej_frozen <= ti$threshold), all(ti$rej_converged <= ti$threshold))
    pw <- T5 %>% filter(role == "power", beta_sep == 2)
    claim("L", "power 0.30-0.57 at beta_sep = 2 (n = 400-800)", "E5c power cells, beta_sep = 2",
          "descriptive: range, and whether any cell changes by > 2 MCSE (paired)",
          sprintf("%.3f-%.3f", min(pw$rej_frozen), max(pw$rej_frozen)),
          sprintf("%.3f-%.3f; changed in %d of %d cells", min(pw$rej_converged), max(pw$rej_converged),
                  sum(pw$changed), nrow(pw)), NA, NA)
  }
}

if (exists("CL") && length(CL)) {
  CL <- do.call(rbind, CL)
  wcsv(CL, "tab_R_claims.csv")
  cat("\n=== WHICH FROZEN CLAIMS HOLD FOR CONVERGED ESTIMATORS ===\n")
  for (i in seq_len(nrow(CL))) {
    cat(sprintf("\n[%s] %s  (%s; criterion: %s)\n      frozen:    %s  -> %s\n      converged: %s  -> %s\n",
                CL$id[i], CL$claim[i], CL$source[i], CL$criterion[i], CL$frozen[i],
                if (is.na(CL$holds_frozen[i])) "descriptive" else if (CL$holds_frozen[i]) "holds" else "does not hold",
                CL$converged[i],
                if (is.na(CL$holds_converged[i])) "descriptive" else if (CL$holds_converged[i]) "HOLDS" else "DOES NOT HOLD"))
  }
}

## ================= matched regularisation (sensitivity analysis) =============
if (length(MR)) {
  REG <- c(C1 = "mechanism only", C2 = "profile separation", C5 = "non-Gaussian features")
  MS <- bind_rows(lapply(names(MR), function(key)
    summarise_cells(MR[[key]], c(BY, "cell", "method")) %>% mutate(experiment = key))) %>%
    left_join(MATCHED_CONFIGS, by = "method") %>% mutate(regime = REG[experiment])
  MF <- bind_rows(lapply(names(MR), function(key) FIT_DIAG(MR[[key]], key)))
  MS <- MS %>% left_join(MF %>% select(experiment, all_of(BY), method, fits, failed, degenerate,
                                       cap_hits), by = c("experiment", BY, "method"))
  wcsv(MS, "tab_M_configs.csv")
  cat("\n=== SENSITIVITY ANALYSIS: matched regularisation, every configuration (tolerance 1e-8) ===\n")
  print(as.data.frame(MS %>% arrange(experiment, cell, family, method) %>%
        transmute(experiment, cell, mu_sep, beta_sep, family,
                  config = ifelse(family == "GeM-Cox", sprintf("normalize=%s temp=%g", normalize, temp),
                                  sprintf("gating x%g", gate_mult)),
                  recovery = round(recovery, 3), mcse = round(rec_mcse, 3), ARI = round(ari, 3),
                  mse = round(mse_eta, 3), sharpness = round(sharpness, 3), degenerate, cap_hits)),
        row.names = FALSE)
  best <- MS %>% group_by(experiment, cell, family) %>% slice_max(recovery, n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(experiment, cell, family, method)
  MV <- bind_rows(lapply(names(MR), function(key) {
    R <- MR[[key]]
    bind_rows(lapply(unique(R$cell), function(cid) {
      Rc <- R %>% filter(cell == cid)
      bg <- best$method[best$experiment == key & best$cell == cid & best$family == "GeM-Cox"]
      bl <- best$method[best$experiment == key & best$cell == cid & best$family == "Latent-class Cox"]
      bind_rows(pd(Rc, GEM_TIGHT, LCC_TIGHT, "default vs default"),
                pd(Rc, bg, bl, "best vs best (oracle-selected upper bound)")) %>%
        mutate(experiment = key, cell = cid, best_gem = bg, best_lcc = bl)
    }))
  })) %>% mutate(verdict = verdict_cell(diff, diff_mcse), regime = REG[experiment])
  wcsv(MV, "tab_M_verdicts.csv")
  cat("\n  GeM-Cox minus the competitor (paired, degenerate excluded):\n")
  print(as.data.frame(MV %>% transmute(contrast, experiment, cell, mu_sep, beta_sep, pairs,
        diff = round(diff, 3), mcse = round(diff_mcse, 3), verdict,
        best_gem = ifelse(contrast == "default vs default", "", best_gem),
        best_lcc = ifelse(contrast == "default vs default", "", best_lcc))), row.names = FALSE)
  MRV <- MV %>% group_by(contrast, regime) %>%
    summarise(cells = n(), verdict = verdict_regime(verdict), .groups = "drop")
  wcsv(MRV, "tab_M_verdicts_regime.csv")
  cat("\n  per regime:\n"); print(as.data.frame(MRV), row.names = FALSE)
}

## ================= figure: E2 at convergence (four panels, as fig2) ===========
if (!is.null(R1$E2)) {
  CMP <- c(GEM0_TIGHT, TS, SC, LCC_TIGHT)
  cols <- stats::setNames(c(YEL, SKY, GRY, GRN), CMP)
  P2 <- paired_diffs(R1$E2, BY, target = GEM_TIGHT, comparators = CMP)
  g <- ggplot(P2 %>% mutate(comparator = factor(comparator, CMP)),
              aes(factor(beta_sep), diff, colour = comparator, group = comparator)) +
    geom_hline(yintercept = 0, linetype = 2, colour = GRY) +
    geom_line(linewidth = 1) +
    geom_pointrange(aes(ymin = diff - diff_mcse, ymax = diff + diff_mcse), size = .35) +
    scale_colour_manual(values = cols, labels = paste("minus", CMP)) +
    facet_wrap(~ mu_sep, nrow = 1, labeller = label_both) +
    labs(title = "Paired gain in contrast recovery at convergence (E2's datasets)",
         subtitle = paste("n = 800, p = 10. GeM-Cox, gamma = 0 and the competitor to tolerance 1e-8;",
                          "two-stage and single Cox as frozen. Bars: +/- 1 MCSE. Degenerate fits excluded."),
         x = "Difference between survival mechanisms (beta_sep)",
         y = "Paired difference in contrast recovery\n(GeM-Cox gamma = 1, converged, minus comparator)") + th
  ggsave(file.path(FIGURES, "figR1_E2_paired_converged.pdf"), g, width = 12, height = 4.6)
  cat("\n  saved figR1_E2_paired_converged.pdf\n")
}
cat("\nDone.\n")
