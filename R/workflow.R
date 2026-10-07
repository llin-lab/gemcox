## workflow.R -- a step-by-step workflow in the style of the GeMLR package
## (read_data -> fit_model -> runCV -> finalModel -> plot_beta_heatmap).
## These are thin wrappers around gemcox(), gemcox_cv_loglik()'s
## cross-validation and predict(); the estimator is unchanged.

#' Read a dataset and prepare the GeM-Cox inputs
#'
#' Reads a file, or takes a data frame, and splits it into the biomarker
#' matrix `X`, the survival outcome (`time`, `status`) and optional indicator
#' covariates `Indi` (for example vaccine group). It mirrors `read_data()` in
#' the GeMLR package. Indicator covariates enter the Cox models only, not the
#' Gaussian mixture.
#'
#' @param dat_path Path to a data file (`.csv`, `.tsv`, `.txt`, `.rds`,
#'   `.rda`/`.RData` holding one data frame, or `.xlsx` if readxl is
#'   installed), or a data frame or matrix.
#' @param time_col Column (name or index) with follow-up times.
#' @param status_col Column (name or index) with the event indicator
#'   (1 = event, 0 = censored).
#' @param Indi_col Optional column(s) (names or indices) of covariates used
#'   only in the Cox models, for example a 0/1 vaccine indicator. `NULL`
#'   (default) means none.
#' @param encoding Optional file encoding for text files.
#' @return A list with `dim` (number of biomarkers), `numdata` (number of
#'   subjects), `rawdat` (the data as read), `X` (biomarkers, data frame),
#'   `Xs` (biomarkers standardised to mean 0, SD 1; 0/1 columns kept),
#'   `time`, `status`, and `Indi` (matrix or `NULL`).
#' @examples
#' f <- system.file("extdata", "gemcox_example.csv", package = "gemcox")
#' dat <- read_data(f, time_col = "time", status_col = "status", Indi_col = "vaccine")
#' dat$dim
#' head(dat$X)
#' @export
read_data <- function(dat_path, time_col, status_col, Indi_col = NULL, encoding = NULL) {
  if (is.data.frame(dat_path) || is.matrix(dat_path)) {
    rawdat <- as.data.frame(dat_path, stringsAsFactors = FALSE)
  } else {
    if (!is.character(dat_path) || length(dat_path) != 1 || !file.exists(dat_path)) {
      stop("dat_path must be an existing file or a data frame.", call. = FALSE)
    }
    ext <- tolower(tools::file_ext(dat_path))
    rd <- function(f, ...) {
      if (is.null(encoding)) f(dat_path, ...) else f(dat_path, ..., fileEncoding = encoding)
    }
    rawdat <- switch(ext,
      csv = rd(utils::read.csv, stringsAsFactors = FALSE, check.names = FALSE),
      tsv = , tab = rd(utils::read.delim, stringsAsFactors = FALSE, check.names = FALSE),
      txt = , dat = rd(utils::read.table, header = TRUE, stringsAsFactors = FALSE,
                       check.names = FALSE),
      rds = readRDS(dat_path),
      rda = , rdata = {
        e <- new.env()
        nm <- load(dat_path, envir = e)
        dfs <- nm[vapply(nm, function(x) is.data.frame(e[[x]]), NA)]
        if (length(dfs) != 1) stop("The .RData file must hold exactly one data frame.", call. = FALSE)
        e[[dfs]]
      },
      xlsx = , xls = {
        if (!requireNamespace("readxl", quietly = TRUE)) {
          stop("Reading Excel files needs the readxl package.", call. = FALSE)
        }
        as.data.frame(readxl::read_excel(dat_path))
      },
      stop("Unsupported file type: .", ext, call. = FALSE))
    rawdat <- as.data.frame(rawdat, stringsAsFactors = FALSE)
  }
  col_index <- function(col, what) {
    if (is.character(col)) {
      bad <- setdiff(col, names(rawdat))
      if (length(bad)) stop(what, ": column(s) not found: ", paste(bad, collapse = ", "), call. = FALSE)
      match(col, names(rawdat))
    } else {
      col <- as.integer(col)
      if (any(col < 1 | col > ncol(rawdat))) stop(what, ": column index out of range.", call. = FALSE)
      col
    }
  }
  if (missing(time_col) || missing(status_col)) {
    stop("Give time_col and status_col (column names or indices).", call. = FALSE)
  }
  it <- col_index(time_col, "time_col"); ist <- col_index(status_col, "status_col")
  ii <- if (is.null(Indi_col)) integer(0) else col_index(Indi_col, "Indi_col")
  if (length(intersect(ii, c(it, ist)))) stop("Indi_col overlaps time_col or status_col.", call. = FALSE)
  time <- rawdat[[it]]; status <- rawdat[[ist]]
  y <- validate_outcome(time, status, nrow(rawdat))
  xi <- setdiff(seq_len(ncol(rawdat)), c(it, ist, ii))
  if (!length(xi)) stop("No biomarker columns left after removing time, status and Indi.", call. = FALSE)
  X <- rawdat[, xi, drop = FALSE]
  num <- vapply(X, is.numeric, NA)
  if (!all(num)) {
    stop("Biomarker columns must be numeric; non-numeric: ", paste(names(X)[!num], collapse = ", "),
         ". Recode them (for example with model.matrix()) or pass them as Indi_col.", call. = FALSE)
  }
  na <- names(X)[colSums(is.na(X)) > 0]
  if (length(na)) stop("Missing values in: ", paste(na, collapse = ", "), ". Impute or remove first.",
                       call. = FALSE)
  is01 <- function(v) all(v %in% c(0, 1))
  Xs <- X
  for (j in names(X)) if (!is01(X[[j]])) Xs[[j]] <- as.numeric(scale(X[[j]]))
  Indi <- NULL
  if (length(ii)) {
    Indi <- as_feature_matrix(rawdat[, ii, drop = FALSE], "Indi")
  }
  list(dim = ncol(X), numdata = nrow(rawdat), rawdat = rawdat, X = X, Xs = Xs,
       time = y$time, status = y$status, Indi = Indi)
}

# Resolve variable selections given as names or indices of X.
resolve_vars <- function(v, X, what) {
  if (is.null(v)) return(seq_len(ncol(X)))
  if (is.character(v)) {
    bad <- setdiff(v, colnames(X))
    if (length(bad)) stop(what, ": not found in X: ", paste(bad, collapse = ", "), call. = FALSE)
    return(match(v, colnames(X)))
  }
  v <- as.integer(v)
  if (any(v < 1 | v > ncol(X))) stop(what, ": indices out of range [1, ", ncol(X), "].", call. = FALSE)
  v
}
workflow_inputs <- function(X, Indi, vargmm, varcox) {
  X <- as_feature_matrix(X, "X")
  ig <- resolve_vars(vargmm, X, "vargmm"); ic <- resolve_vars(varcox, X, "varcox")
  X_gmm <- X[, ig, drop = FALSE]
  X_cox <- X[, ic, drop = FALSE]
  if (!is.null(Indi)) {
    Indi <- as_feature_matrix(Indi, "Indi")
    if (is.null(colnames(Indi))) colnames(Indi) <- paste0("Indi", seq_len(ncol(Indi)))
    X_cox <- cbind(X_cox, Indi)
  }
  list(X_gmm = X_gmm, X_cox = X_cox, vargmm = colnames(X)[ig], varcox = colnames(X)[ic])
}

#' Fit a GeM-Cox model (quick fit for a given K)
#'
#' Fits GeM-Cox with `K` subgroups, the counterpart of `fit_model()` in
#' GeMLR. Biomarkers in `vargmm` define the subgroups through a Gaussian
#' mixture. Biomarkers in `varcox`, plus any indicator covariates `Indi`,
#' enter each subgroup's Cox model. The result is a `"gemcox"` object
#' ([gemcox()]), so `summary()`, `predict()` and `plot()` work on it. It
#' carries the extra elements listed under Value.
#'
#' @param X Biomarker matrix or data frame (n x p), raw scale (for example
#'   `read_data()$X`). gemcox standardises internally.
#' @param time,status Follow-up time and event indicator (1 = event).
#' @param Indi Optional covariates for the Cox models only (for example
#'   `read_data()$Indi`).
#' @param K Number of subgroups (default 2).
#' @param vargmm Biomarkers (names or indices of `X`) for the Gaussian
#'   mixture. `NULL` = all.
#' @param varcox Biomarkers (names or indices of `X`) for the Cox models.
#'   `NULL` = all.
#' @param lambda Ridge penalty on the Cox coefficients (default 0.05).
#' @param nseeds Number of EM starts; the best is kept (default 5).
#' @param seed Seed for the starting partitions (default 1).
#' @param verbose 0 = quiet (default), 1 = report the fit, 2 = also print the
#'   EM trace.
#' @param ... Further arguments to [gemcox()], for example `gamma`,
#'   `baseline` or `max_iter`.
#' @return A `"gemcox"` fit with, in addition:
#' \describe{
#'   \item{beta_sd}{Coefficients per standard deviation of each Cox covariate
#'     (rows) by subgroup (columns), comparable across covariates; used by
#'     [plot_beta_heatmap()].}
#'   \item{vargmm, varcox, Indi_names}{The variables used.}
#'   \item{metrics}{A list: `loglik` (final log-likelihood), `cindex`
#'     (in-sample concordance of the mixture risk score, with membership
#'     from the biomarkers only; optimistic because it is computed on the
#'     fitting data), `events` (effective events per subgroup), and
#'     `converged`.}
#' }
#' @seealso [runCV()] and [finalModel()] to choose K; [gemcox()] for all
#'   settings.
#' @examples
#' dat <- read_data(gemcox_example, time_col = "time", status_col = "status",
#'                  Indi_col = "vaccine")
#' fit <- fit_model(dat$X, dat$time, dat$status, Indi = dat$Indi, K = 2, nseeds = 2)
#' fit$metrics
#' fit$beta_sd
#' @export
fit_model <- function(X, time, status, Indi = NULL, K = 2, vargmm = NULL, varcox = NULL,
                      lambda = 0.05, nseeds = 5, seed = 1, verbose = 0, ...) {
  inp <- workflow_inputs(X, Indi, vargmm, varcox)
  fit <- gemcox(inp$X_gmm, inp$X_cox, time = time, status = status, K = K, lambda = lambda,
                n_starts = nseeds, seed = seed, verbose = verbose >= 2, ...)
  fit$call <- match.call()
  fit$beta_sd <- fit$beta_scaled
  fit$vargmm <- inp$vargmm; fit$varcox <- inp$varcox
  fit$Indi_names <- if (is.null(Indi)) character(0) else setdiff(colnames(inp$X_cox), inp$varcox)
  lp <- stats::predict(fit, inp$X_gmm, inp$X_cox, type = "lp")
  cidx <- tryCatch(survival::concordance(survival::Surv(fit$data$time, fit$data$status) ~ lp,
                                         reverse = TRUE)$concordance, error = function(e) NA_real_)
  fit$metrics <- list(loglik = utils::tail(fit$loglik, 1), cindex = unname(cidx),
                      events = fit$eff_events, converged = fit$converged)
  if (verbose >= 1) {
    message(sprintf("GeM-Cox, K = %d: %s after %d iterations; in-sample C-index %.3f.", K,
                    if (fit$converged) "converged" else "NOT converged", fit$iterations, cidx))
  }
  fit
}

# Held-out (Verweij & van Houwelingen) log-likelihood per fold and K: the
# per-fold version of cv_loglik_folds(); column sums equal its totals.
cv_loglik_by_fold <- function(X_gmm, X_cox, time, status, folds, K, criterion, args, verbose) {
  fl <- sort(unique(folds))
  out <- matrix(0, length(fl), length(K))
  for (i in seq_along(fl)) {
    tr <- folds != fl[i]
    for (j in seq_along(K)) {
      fit <- do.call(gemcox, c(list(X_gmm = X_gmm[tr, , drop = FALSE], X_cox = X_cox[tr, , drop = FALSE],
                                    time = time[tr], status = status[tr], K = K[j]), args))
      eta <- stats::predict(fit, X_gmm, X_cox, type = "lp")
      val <- cox_partial_loglik(time, status, eta) - cox_partial_loglik(time[tr], status[tr], eta[tr])
      if (criterion == "joint") val <- val + sum(gmm_mixture_logdens(fit, X_gmm[!tr, , drop = FALSE]))
      out[i, j] <- val
    }
    if (verbose >= 1) message(sprintf("  fold %d of %d done", i, length(fl)))
  }
  out
}

#' Cross-validation over the number of subgroups
#'
#' The counterpart of `runCV()` in GeMLR. Each fold is held out in turn;
#' GeM-Cox is fitted on the rest for every K in `ncmp`, and the held-out
#' fold is scored. The score is the cross-validated partial log-likelihood
#' of Verweij and van Houwelingen (1993); with `criterion = "joint"` the
#' held-out Gaussian-mixture log-density is added. Larger is better.
#'
#' GeMLR chooses K by cross-validated AUC. For GeM-Cox the C-index is not
#' used to choose K, by a settled project decision: it measures how well
#' risk is ranked, not whether subgroups differ in how biomarkers act.
#'
#' @inheritParams fit_model
#' @param ncmp Numbers of subgroups to compare (default `1:3`; include 1 to
#'   compare against a single Cox model).
#' @param k Number of folds (default 5). Folds are stratified by event
#'   status.
#' @param criterion `"partial"` (default; the survival part only) or
#'   `"joint"` (also counts how well the feature profiles are modelled;
#'   tends to count feature-profile clusters).
#' @param nseeds Number of EM starts per fit (default 1, for speed).
#' @param verbose 1 (default) reports progress per fold; 0 is quiet.
#' @return A list: `cvLLfinal`, a k x length(ncmp) matrix of held-out
#'   log-likelihoods (rows = folds, columns = `"cluster=K"`); `mean`, the
#'   column means; `best_K`, the K with the largest mean; and `criterion`,
#'   `ncmp` and `folds`.
#' @references Verweij PJM, van Houwelingen HC (1993). Cross-validation in
#'   survival analysis. *Statistics in Medicine* 12:2305-2314.
#' @examples
#' \donttest{
#' cv <- runCV(gemcox_example[, paste0("marker", 1:6)], gemcox_example$time,
#'             gemcox_example$status, ncmp = 1:2, k = 3, verbose = 0)
#' cv$cvLLfinal
#' cv$best_K
#' }
#' @export
runCV <- function(X, time, status, Indi = NULL, ncmp = 1:3, k = 5, vargmm = NULL, varcox = NULL,
                  criterion = c("partial", "joint"), lambda = 0.05, nseeds = 1, seed = 1,
                  verbose = 1, ...) {
  criterion <- match.arg(criterion)
  inp <- workflow_inputs(X, Indi, vargmm, varcox)
  y <- validate_outcome(time, status, nrow(inp$X_gmm))
  check_scalar(k, "k", lower = 2, upper = sum(y$status), integer = TRUE)
  folds <- stratified_folds(y$status, k, seed)
  args <- c(list(lambda = lambda, n_starts = nseeds, seed = seed), list(...))
  n_cap <- 0L
  M <- withCallingHandlers(
    cv_loglik_by_fold(inp$X_gmm, inp$X_cox, y$time, y$status, folds, ncmp, criterion, args, verbose),
    gemcox_not_converged = function(w) { n_cap <<- n_cap + 1L; invokeRestart("muffleWarning") })
  if (n_cap > 0) warning(sprintf("%d of the %d fold fits stopped at max_iter without converging.",
                                 n_cap, k * length(ncmp)), call. = FALSE)
  dimnames(M) <- list(paste(seq_len(k), "fold"), paste0("cluster=", ncmp))
  m <- colMeans(M)
  list(cvLLfinal = M, mean = m, best_K = ncmp[which.max(m)], criterion = criterion,
       ncmp = ncmp, folds = folds)
}

#' Fit the final model at the cross-validated number of subgroups
#'
#' The counterpart of `finalModel()` in GeMLR. It takes the K with the
#' largest mean held-out log-likelihood from [runCV()], and fits it on all
#' the data with [fit_model()].
#'
#' @inheritParams fit_model
#' @param cv The result of [runCV()].
#' @param nseeds Number of EM starts for the final fit (default 10).
#' @return The [fit_model()] result for the selected K, with `K_selected`
#'   and the cross-validation result `cv` attached.
#' @examples
#' \donttest{
#' X <- gemcox_example[, paste0("marker", 1:6)]
#' cv <- runCV(X, gemcox_example$time, gemcox_example$status, ncmp = 1:2, k = 3, verbose = 0)
#' fit <- finalModel(cv, X, gemcox_example$time, gemcox_example$status, nseeds = 2)
#' fit$K_selected
#' }
#' @export
finalModel <- function(cv, X, time, status, Indi = NULL, vargmm = NULL, varcox = NULL,
                       lambda = 0.05, nseeds = 10, seed = 1, verbose = 0, ...) {
  if (!is.list(cv) || is.null(cv$best_K)) stop("cv must be the result of runCV().", call. = FALSE)
  if (verbose >= 1) {
    message(sprintf("Selected K = %d (largest mean held-out log-likelihood, %s criterion).",
                    cv$best_K, cv$criterion))
  }
  fit <- fit_model(X, time, status, Indi = Indi, K = cv$best_K, vargmm = vargmm, varcox = varcox,
                   lambda = lambda, nseeds = nseeds, seed = seed, verbose = verbose, ...)
  fit$K_selected <- cv$best_K
  fit$cv <- cv
  fit$call <- match.call()
  fit
}

#' Heatmap of subgroup-specific Cox coefficients
#'
#' The counterpart of `plot_beta_heatmap()` in GeMLR. Rows are covariates,
#' columns are subgroups, and colour gives the sign and size of each
#' coefficient (red = higher hazard, blue = lower), with values printed in
#' the cells. Uses base graphics only.
#'
#' @param beta A [fit_model()] result (its per-SD coefficients `beta_sd` are
#'   plotted), any `"gemcox"` fit, or a coefficient matrix (covariates x
#'   subgroups).
#' @param output_file `NULL` (default) draws on the current device; a `.png`
#'   or `.pdf` file name saves the plot there.
#' @param width,height Size in inches when saving; `res` is the PNG
#'   resolution.
#' @param res Resolution (dpi) for PNG output.
#' @param main Title.
#' @return The plotted matrix, invisibly.
#' @examples
#' dat <- read_data(gemcox_example, time_col = "time", status_col = "status")
#' fit <- fit_model(dat$X, dat$time, dat$status, K = 2, nseeds = 2)
#' plot_beta_heatmap(fit)
#' @export
plot_beta_heatmap <- function(beta, output_file = NULL, width = 5.5, height = 5, res = 300,
                              main = "Cox coefficients by subgroup (per SD)") {
  if (inherits(beta, "gemcox")) {
    beta <- if (!is.null(beta$beta_sd)) beta$beta_sd else beta$beta_scaled
  }
  beta <- as.matrix(beta)
  if (!is.null(output_file)) {
    if (grepl("\\.pdf$", output_file, ignore.case = TRUE)) {
      grDevices::pdf(output_file, width = width, height = height)
    } else {
      grDevices::png(output_file, width = width, height = height, units = "in", res = res)
    }
    on.exit(grDevices::dev.off())
  }
  m <- max(abs(beta), 1e-8)
  pal <- grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(101)
  nr <- nrow(beta); nc <- ncol(beta)
  op <- graphics::par(mar = c(3, max(6, max(nchar(rownames(beta)), 0) * 0.55), 3, 1))
  on.exit(graphics::par(op), add = TRUE)
  graphics::image(seq_len(nc), seq_len(nr), t(beta[nr:1, , drop = FALSE]), col = pal,
                  zlim = c(-m, m), axes = FALSE, xlab = "", ylab = "", main = main)
  graphics::axis(1, at = seq_len(nc), labels = colnames(beta), tick = FALSE)
  graphics::axis(2, at = seq_len(nr), labels = rev(rownames(beta)), las = 1, tick = FALSE)
  graphics::text(rep(seq_len(nc), each = nr), rep(nr:1, nc), sprintf("%.2f", as.vector(beta)),
                 cex = 0.8)
  graphics::box()
  invisible(beta)
}
