## utils.R -- input validation, feature scaling, internal helpers.

#' gemcox: Gaussian mixture of Cox models for mechanism subgroups
#'
#' Fits latent subgroups whose survival mechanisms (Cox coefficient
#' vectors) differ, with a shared baseline hazard. The main entry points
#' are [gemcox()], [predict.gemcox()] and [gemcox_heterogeneity_test()].
#'
#' @keywords internal
#' @importFrom stats coef predict kmeans sd
#' @importFrom survival Surv coxph
#' @importFrom glmnet glmnet
"_PACKAGE"

## ---- input validation ------------------------------------------------

## Coerce a feature matrix while keeping the ORIGINAL column names.
## Names are never passed through make.names(): as.data.frame()-style
## mangling ("IgG3 AMA-1 MFI | V20" -> "IgG3.AMA.1.MFI...V20") silently
## broke feature classification in earlier application code.
as_feature_matrix <- function(X, what) {
  if (is.data.frame(X)) {
    num <- vapply(X, is.numeric, logical(1))
    if (!all(num)) {
      stop(sprintf("%s: all columns must be numeric; non-numeric column(s): %s",
                   what, paste(names(X)[!num], collapse = ", ")), call. = FALSE)
    }
    nm <- names(X)
    X <- as.matrix(X)
    colnames(X) <- nm
  } else if (is.numeric(X) && is.null(dim(X))) {
    X <- matrix(X, ncol = 1)
  }
  if (!is.matrix(X) || !is.numeric(X)) {
    stop(sprintf("%s must be a numeric matrix or data frame.", what), call. = FALSE)
  }
  if (nrow(X) == 0 || ncol(X) == 0) {
    stop(sprintf("%s has no rows or no columns.", what), call. = FALSE)
  }
  if (anyNA(X)) {
    bad <- colnames(X)[colSums(is.na(X)) > 0]
    stop(sprintf("%s contains NA values%s. Impute or remove them before fitting.",
                 what, if (length(bad)) paste0(" (columns: ", paste(bad, collapse = ", "), ")") else ""),
         call. = FALSE)
  }
  if (any(!is.finite(X))) {
    stop(sprintf("%s contains non-finite values (Inf/NaN).", what), call. = FALSE)
  }
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
  if (anyDuplicated(colnames(X))) {
    stop(sprintf("%s has duplicated column names: %s", what,
                 paste(unique(colnames(X)[duplicated(colnames(X))]), collapse = ", ")),
         call. = FALSE)
  }
  storage.mode(X) <- "double"
  X
}

validate_outcome <- function(time, status, n) {
  if (!is.numeric(time) || length(time) != n) {
    stop(sprintf("time must be a numeric vector of length %d (the number of rows of X_gmm); got length %d.",
                 n, length(time)), call. = FALSE)
  }
  if (anyNA(time) || any(!is.finite(time))) {
    stop("time contains NA or non-finite values.", call. = FALSE)
  }
  if (any(time < 0)) {
    stop(sprintf("time must be non-negative; %d value(s) are negative.", sum(time < 0)),
         call. = FALSE)
  }
  if (is.logical(status)) status <- as.integer(status)
  if (!is.numeric(status) || length(status) != n) {
    stop(sprintf("status must be a vector of length %d coded 0 (censored) / 1 (event); got length %d.",
                 n, length(status)), call. = FALSE)
  }
  if (anyNA(status) || !all(status %in% c(0, 1))) {
    stop("status must be binary, coded 0 (censored) / 1 (event), with no NA.",
         call. = FALSE)
  }
  if (sum(status) == 0) stop("status contains no events.", call. = FALSE)
  list(time = as.numeric(time), status = as.integer(status))
}

check_scalar <- function(x, name, lower = -Inf, upper = Inf, integer = FALSE,
                         lower_open = FALSE) {
  ok <- is.numeric(x) && length(x) == 1 && is.finite(x) &&
    (if (lower_open) x > lower else x >= lower) && x <= upper &&
    (!integer || x == round(x))
  if (!ok) {
    rng <- sprintf("%s%s, %s]", if (lower_open) "(" else "[", format(lower), format(upper))
    stop(sprintf("%s must be a single %s in %s.", name,
                 if (integer) "integer" else "number", rng), call. = FALSE)
  }
  invisible(x)
}

## Reorder new data to the training columns, matching on original names.
match_columns <- function(Xnew, ref_names, what) {
  if (ncol(Xnew) != length(ref_names)) {
    stop(sprintf("%s has %d columns; the model was fitted with %d.",
                 what, ncol(Xnew), length(ref_names)), call. = FALSE)
  }
  if (!all(ref_names %in% colnames(Xnew))) {
    if (all(grepl("^x[0-9]+$", colnames(Xnew)))) return(Xnew)  # unnamed input
    missing <- setdiff(ref_names, colnames(Xnew))
    stop(sprintf("%s is missing column(s) used in fitting: %s", what,
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
  Xnew[, ref_names, drop = FALSE]
}

## ---- scaling ---------------------------------------------------------

# Store column centers/scales for leakage-safe standardization.
make_x_scaler <- function(X) {
  X <- as.matrix(X)
  center <- colMeans(X, na.rm = TRUE)
  scalev <- apply(X, 2, stats::sd, na.rm = TRUE)
  scalev[!is.finite(scalev) | scalev <= 0] <- 1
  list(center = center, scale = scalev)
}

# Apply a previously learned scaler to train, validation, or test data.
apply_x_scaler <- function(X, scaler) {
  X <- as.matrix(X)
  sweep(sweep(X, 2, scaler$center, "-"), 2, scaler$scale, "/")
}

## ---- numerical helpers -----------------------------------------------

# Stable log(sum(exp(.))) by row, used in responsibility calculations.
rowLogSumExp <- function(A) {
  m <- apply(A, 1, max)
  m[!is.finite(m)] <- 0
  s <- rowSums(exp(A - m))
  s <- pmax(s, 1e-300)
  m + log(s)
}

## ---- numerical guards ------------------------------------------------
## Every clipping bound in the engine reports when it binds. Low-level
## functions signal a classed warning ("gemcox_guard_warning") each time;
## gemcox() collects them during fitting and re-issues one summary warning,
## so a fit reports clipping once rather than once per EM iteration.

guard_warning <- function(guard, n_bound, n_total) {
  msg <- sprintf(paste0("gemcox: numerical guard '%s' bound for %d of %d values; ",
                        "estimates may be attenuated."), guard, n_bound, n_total)
  cond <- structure(
    class = c("gemcox_guard_warning", "warning", "condition"),
    list(message = msg, call = NULL, guard = guard,
         n_bound = n_bound, n_total = n_total))
  warning(cond)
}

check_guard <- function(bound, guard) {
  nb <- sum(bound)
  if (nb > 0) guard_warning(guard, nb, length(bound))
  invisible(nb)
}

# Clamp a linear predictor to [-cap, cap], warning when the bound binds.
clamp_eta <- function(eta, cap, guard) {
  check_guard(abs(eta) > cap, guard)
  pmin(pmax(eta, -cap), cap)
}

# Compute a capped Cox linear predictor from raw features and scaled betas.
eta_from_scaled <- function(X, beta_s, scaler, eta_cap = 12) {
  Xs <- apply_x_scaler(X, scaler)
  eta <- as.numeric(Xs %*% beta_s)
  clamp_eta(eta, eta_cap, "eta_cap")
}

# Evaluate expr, muffling guard warnings and returning their tally.
collect_guards <- function(expr) {
  tally <- list()
  value <- withCallingHandlers(expr, gemcox_guard_warning = function(w) {
    cur <- tally[[w$guard]]
    if (is.null(cur)) cur <- c(events = 0, bound = 0, total = 0)
    tally[[w$guard]] <<- cur + c(1, w$n_bound, w$n_total)
    invokeRestart("muffleWarning")
  })
  guards <- if (length(tally)) {
    data.frame(guard = names(tally),
               events = vapply(tally, `[[`, 0, 1),
               values_bound = vapply(tally, `[[`, 0, 2),
               values_checked = vapply(tally, `[[`, 0, 3),
               row.names = NULL, stringsAsFactors = FALSE)
  } else {
    data.frame(guard = character(), events = numeric(), values_bound = numeric(),
               values_checked = numeric(), stringsAsFactors = FALSE)
  }
  list(value = value, guards = guards)
}

warn_guard_summary <- function(guards, where = "fitting") {
  if (!nrow(guards)) return(invisible())
  txt <- paste(sprintf("%s (%g of %g values)", guards$guard, guards$values_bound,
                       guards$values_checked), collapse = "; ")
  msg <- sprintf(paste0("gemcox: numerical guards bound during %s: %s. ",
                        "Estimates may be attenuated; see the 'guards' element of the fit ",
                        "and gemcox_control()."), where, txt)
  warning(structure(class = c("gemcox_guard_summary", "warning", "condition"),
                    list(message = msg, call = NULL, guards = guards)))
}

# Warn when EM stopped at max_iter without meeting the tolerance. Classed,
# so that functions running many fits (the heterogeneity test, the CV
# criterion) can count these warnings and report them once.
warn_not_converged <- function(converged, iterations, tol) {
  if (isTRUE(converged)) return(invisible())
  msg <- sprintf(paste0("gemcox: EM stopped at max_iter = %d without meeting tol = %g; ",
                        "estimates may not be at the optimum. Increase max_iter."),
                 as.integer(iterations), tol)
  warning(structure(class = c("gemcox_not_converged", "warning", "condition"),
                    list(message = msg, call = NULL)))
}

## ---- reproducible randomness without touching the caller's stream ------

# Evaluate expr after set.seed(seed), then restore the global RNG state.
with_local_seed <- function(seed, expr) {
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = genv)
    } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      rm(".Random.seed", envir = genv)
    }
  })
  set.seed(seed)
  expr
}
