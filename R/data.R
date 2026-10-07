#' Example data: simulated biomarkers and time to infection
#'
#' A simulated dataset for trying out the package. It is not trial data. It
#' has 300 subjects in two latent subgroups. Their biomarker profiles
#' overlap (mean separation 1 SD), but the biomarkers relate to the hazard
#' differently in the two subgroups: the Cox coefficients differ mainly for
#' `marker1` to `marker3`. `vaccine` is a randomised 0/1 indicator with no
#' effect on the outcome. Everyone is followed to day 150 (administrative
#' censoring). The same data are in
#' `system.file("extdata", "gemcox_example.csv", package = "gemcox")`.
#'
#' @format A data frame with 300 rows and 9 columns:
#' \describe{
#'   \item{vaccine}{0/1 indicator (randomised; no effect).}
#'   \item{marker1, marker2, marker3, marker4, marker5, marker6}{Biomarker
#'     levels (standardised scale).}
#'   \item{time}{Follow-up time (days; at most 150).}
#'   \item{status}{1 = event (infection), 0 = censored.}
#' }
#' @source `gemcox_simulate(n = 300, p = 6, mu_sep = 1, beta_sep = 2,
#'   censoring = "administrative", followup = 150, seed = 2026)`; see
#'   `data-raw/make_gemcox_example.R`.
"gemcox_example"
