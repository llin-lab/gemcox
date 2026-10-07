# gemcox 0.3.0

For sharing within the research team; not a public release.

## A step-by-step workflow in the style of GeMLR

New functions, which are thin wrappers around the existing estimator:

* `read_data()` reads a file (`.csv`, `.tsv`, `.txt`, `.rds`, `.RData`,
  `.xlsx`) or a data frame into biomarkers `X`, standardised `Xs`,
  `time`, `status`, and indicator covariates `Indi` (used in the Cox models
  only).
* `fit_model()` fits GeM-Cox for a given K. It selects biomarkers for
  clustering (`vargmm`) and for the Cox models (`varcox`), and returns per-SD
  coefficients `beta_sd` and in-sample `metrics`.
* `runCV()` gives the held-out log-likelihood for each fold and each K
  (`cvLLfinal`, folds x K). Its column sums equal `gemcox_cv_loglik()`.
* `finalModel()` refits at the K with the largest mean held-out
  log-likelihood.
* `plot_beta_heatmap()` draws a heatmap of the subgroup coefficients
  (base graphics; it can save to PNG or PDF).
* `gemcox_example` is a simulated example dataset (300 subjects, 6
  biomarkers, a vaccine indicator, follow-up to day 150), also available
  as `inst/extdata/gemcox_example.csv`.

Differences from GeMLR: K is chosen by the held-out log-likelihood, not
by AUC or C-index; no packages are installed or attached by the
functions. The workflow functions share their names with GeMLR's; when
both packages are loaded, use `gemcox::`.

The README follows GeMLR's step-by-step format, with an "Understanding the
output" section.

# gemcox 0.2.0

For sharing within the research team; not a public release.

## Behaviour change: EM runs to convergence by default

* `gemcox()` (and the internal `lcc_fit()`) now default to `tol = 1e-8`
  and `max_iter = 1000`, instead of `tol = 1e-4` and `max_iter = 100`.
  - In the pre-registered tolerance study (`inst/sim/09_tolerance_study.R`;
    480 simulated datasets), only 32% of fits at the old setting gave a
    contrast within 0.01 of the converged fit. 1e-7 reached 90%, with the
    worst scenario at 65%; only 1e-8 met the rule.
  - The cost is about 7 times more fitting time (median 3.3 s per fit at
    n = 800, p = 10).
  - The defaults are read from `getOption("gemcox.tol")` and
    `getOption("gemcox.max_iter")`. `options(gemcox.tol = 1e-4,
    gemcox.max_iter = 100)` restores the old behaviour, and the simulation
    pipeline sets it, so the frozen results still reproduce exactly
    (checked on refitted replicates).
  - `gemcox_heterogeneity_test()` and `gemcox_cv_loglik()` inherit the new
    defaults.
* **The test's calibration** (type I error, power) was established at the
  old setting and is provisional at the new default; the help page, README
  and vignette say so.

## New warnings

* A fit that stops at `max_iter` without meeting `tol` now warns (class
  `gemcox_not_converged`).
* `gemcox_heterogeneity_test()` and `gemcox_cv_loglik()` collect these
  warnings and report one summary. The test also returns `n_fits` and
  `n_not_converged`.

## Simulation

* `gemcox_simulate(censoring = "administrative", followup = )`: every
  subject is followed to the same time, with no other censoring
  (experiment E6). The default, random censoring calibrated to
  `event_rate`, is unchanged; the frozen datasets' hashes still match.
  The result also includes the true cluster means, `mu_list`.

## Documentation

* `DESCRIPTION` describes the default test (likelihood ratio with a
  parametric bootstrap null), not the deprecated permutation test.
* The README has installation notes for the team, a "Status of the
  evidence" section (what is established at convergence and what is
  provisional), run-time notes and the new default.

# gemcox 0.1.0

First packaged version. The engine is ported from `GeMCox_combined_PATCHED.R`
(sha256 5e42cb0f...). `tests/testthat/test-equivalence.R` checks that
`gemcox()` reproduces that file's fitters to machine precision when given the
same arguments.

## Bug fixes relative to `GeMCox_combined.R`

* `get_cumhaz_at()`: the original ended with
  `ifelse(idx > 0, baseline$cumhaz[idx], 0)`. When a query time preceded
  the first event time, `idx` contained 0, `cumhaz[0]` dropped elements,
  and `ifelse()` recycled the shorter vector. Every subject after the first
  pre-grid time received another subject's cumulative hazard, which
  attenuated shared-baseline fits roughly 2x. Only valid positions are now
  indexed.
* The same defect was still present, inline, in `log_surv_density()`: the
  E-step survival term for both fitters. It fired whenever any subject was
  observed before the first event time (158 of 200 simulated datasets at
  n = 400), and then corrupted the term for almost every subject. It
  affected every `gamma > 0` fit made with `GeMCox_combined.R`. The
  survival density now reads cumulative hazards only through
  `get_cumhaz_at()`.
* Feature names are kept verbatim. The legacy fitters passed column names
  through `make.names()`, which mangles names such as
  `"IgG3 AMA-1 MFI | V20"`.

## Algorithm changes

* **The shared-baseline Poisson M-step is fitted without an intercept
  (`intercept = FALSE`).** `GeMCox_combined.R` fitted a per-cluster intercept
  and discarded it. A per-cluster intercept lets clusters differ in baseline
  risk, which contradicts the shared-baseline model, in which clusters
  differ only through `beta`. Estimates differ from `GeMCox_combined.R` for
  this reason as well as the bug fixes above.
* The returned membership weights `tau` are recomputed at the returned
  parameters (one final E-step), and multi-start selects the start with the
  highest final E-step score. (Both from PATCHED.)

## Defaults

* One fitter, `gemcox()`, replaces `gemcox_full()` (`baseline = "cluster"`),
  `gemcox_full_shared_v4()` (`baseline = "shared"`, default) and
  `gemcox_full_multistart_shared_v4()` (`n_starts > 1`).
* `normalize_gmm_by_dim = FALSE` everywhere. Legacy `gemcox_full()` and the
  `cv_select_K_*` functions defaulted to `TRUE`.
* `init = "kmeans"` on the scaled clustering features (or `"random"`).
  PATCHED's default `"supervised"` initialisation, which uses the outcome, is
  not offered: it would make the `gamma = 0` comparator outcome-informed.
* `lambda = 0.05`, `alpha = 0` (ridge). Legacy defaults were `alpha = 0.5`
  and `lambda = NULL`, which selected `lambda.1se` by `cv.glmnet` and gave
  all-zero components at p = 80.
* `gamma = 1`, `max_iter = 100` (legacy 60 or 50), `tol = 1e-4`.
* `normalize_gmm_by_dim` stays `FALSE` after a pre-specified comparison
  (simulation experiment E2b). The comparison was with an adaptive rule:
  use `TRUE` when the held-out Gaussian-mixture part of the joint
  cross-validated criterion prefers one feature cluster.
  * Averaged over 12 scenarios, the adaptive rule improved contrast
    recovery by 0.085 (MCSE 0.004).
  * It was worse in one scenario (mu_sep = 1, beta_sep = 1: -0.030, MCSE
    0.011), which the pre-registered rule did not allow.
  * `TRUE` alone was markedly better when feature profiles barely differed
    (mu_sep <= 0.5), and worse once they clearly differed (mu_sep = 2).

## Numerical guards now warn

* The Breslow caps (`max_jump = 10`, `max_cumhaz = 500`), the risk-set
  denominator floor, the linear-predictor caps (`eta_cap = 12`,
  `eta_clamp = 30`) and the replacement of non-finite E-step evidence all
  used to clip silently. They now warn when they bind. `gemcox()` collects
  these into one warning per fit and records them in `fit$guards`. The
  thresholds are unchanged and set through `gemcox_control()`.

## Testing K = 2 against K = 1

* `gemcox_heterogeneity_test()` is the test of heterogeneous survival
  mechanisms. Its default is the in-sample likelihood ratio (`statistic =
  "lrt"`) with the parametric bootstrap null. This choice follows the
  package simulations, experiments E5a and E5c in `inst/sim`:
  * Only this combination was both calibrated and powered: rejection 0.052
    under the null (MCSE 0.010), and power 0.30-0.57 at beta_sep = 2.
  * The cross-validated statistic was calibrated with the bootstrap (0.062)
    but had power 0.05-0.10 in the same cells.
  * The joint permutation null was anti-conservative for the
    cross-validated statistic (0.112) and conservative for the LRT (0.028).
* The LRT requires `gamma = 1` and `temp = 1`, the only settings in which
  the fitted objective is the mixture log-likelihood.
* `gemcox_permutation_test()` is deprecated. It keeps its original
  signature and defaults (`statistic = "cv"`, `n_perm`), so existing code
  returns the same results, and it signals a deprecation warning.

## Added after the frozen simulation study (sim-freeze-v1)

* `gemcox_simulate(K_true = 3)`: three clusters at the vertices of
  equilateral triangles in profile and coefficient space, in two further
  orthogonal directions. It needs p >= 8. Draws for `K_true = 1, 2` are
  unchanged; a test checks them against the frozen version.
* `gemcox_heterogeneity_test(K0 = )` tests `K0 + 1` against `K0` clusters,
  for sequential selection of K. The default `K0 = 1` is unchanged.
  `K0 >= 2` needs the LRT with the bootstrap null. The bootstrap then draws
  cluster membership from the fitted model's probabilities given the
  features. Its calibration for `K0 >= 2` has not been established.

## Findings added after the frozen study (sim-freeze-v2)

* **E5d: not validated for profile-only structure.** When subgroups
  differed only in feature profile, the default test rejected at 0.070
  (MCSE 0.011) at mu_sep = 0.5, above the pre-registered tolerance of
  0.0695. It was calibrated at mu_sep = 1 and 2 and with correlated
  features. By the rule recorded in advance, the test is not validated
  for data with profile structure; the help page and README say so.
* **E5e: selecting K on the grid 1:4.**
  - The joint cross-validated criterion and BIC count feature-profile
    clusters, and recover K = 3 only with strong separation.
  - The sequential bootstrap LRT was the only selector sensitive to
    mechanism differences without profile separation.
  - The partial cross-validated criterion selected K >= 2 in 13% of null
    datasets.
* **E5f: the adaptive-normalisation test** was anti-conservative up to
  about 88 events and at the nominal level by about 132 events.

## Internal competitor (added after sim-freeze-v2)

* `lcc_fit()` (internal, not exported): latent-class Cox with multinomial
  logistic gating, for comparison with GeM-Cox in the simulation study.
  Apart from the gating, it reuses GeM-Cox's machinery: the M-step,
  baseline, E-step density, ridge penalty, convergence rule,
  initialisation and guards. The gating is fitted by glmnet
  (multinomial, alpha = 0, lambda_gate = lambda by default, unpenalised
  intercept). No existing function or default changed.
* `gemcox_simulate(features = "lognormal", sdlog = 0.8)`: standardised,
  right-skewed features about the same cluster means, for the
  misspecified-features experiment C5. The default, `"gaussian"`, is the
  frozen DGP. Checksum tests, and a manifest of all 5,200 reused datasets
  (`inst/sim/05_hash_manifest.R`), confirm that the Gaussian datasets are
  unchanged.
* `gemcox_heterogeneity_test()` returns, for the LRT, `n_fits` and
  `n_not_converged`: how many of its fits stopped at `max_iter` without
  meeting the tolerance. This is additive; the statistic and p-value are
  unchanged.

## New functions

* `predict.gemcox()`: the single path for scoring new subjects. Membership
  weights come from the clustering features only. Like the legacy
  `predict_tau()`, it does not apply `temp`.
* `gemcox_cv_loglik()` and `gemcox_permutation_test()`: the cross-validated
  partial likelihood of Verweij and van Houwelingen (1993), with a mixture
  linear predictor and `tau` from the clustering features. It is optionally
  joined with the held-out Gaussian-mixture log-density. Breslow jumps are
  not used as a density for held-out subjects.
* `gemcox_permutation_test(null = )` (now `gemcox_heterogeneity_test()`)
  offers two null distributions for K = 2 vs K = 1:
  * `"bootstrap"` (the default) is a parametric bootstrap under the fitted
    K = 1 model (Breslow baseline, Kaplan-Meier censoring). It keeps a
    common Cox effect.
  * `"permutation"` permutes `(time, status)` jointly. Strictly this tests
    "no association between features and outcome".

  **Why the bootstrap is the default.** The question is whether survival
  mechanisms differ between latent subgroups, so the null is "one mechanism
  with a shared beta". The joint permutation also destroys the shared
  effect: its null data have beta = 0, a different hypothesis. It is
  calibrated for "one mechanism" only if the gain statistic happens to have
  the same distribution with and without a common effect, which is not
  guaranteed. The bootstrap keeps the fitted common effect, at the cost of
  relying on the fitted K = 1 model. The permutation null stays available;
  experiment E5a compares the two.
* `observed =` reuses the observed statistic and folds of an earlier call,
  so both nulls can be compared against one observed statistic.
* `gemcox_permutation_test(statistic = "lrt")`: the in-sample likelihood
  ratio 2 (l2 - l1) on the mixture log-likelihood, calibrated by either
  null. The chi-square reference is not used, since it is invalid for
  mixtures. Requires `gamma = 1` and `temp = 1`, the only setting in which
  the fitted objective is the log-likelihood. The default statistic stays
  `"cv"`.
* `gemcox_cv_loglik()` returns the per-K partial-likelihood and held-out
  GMM log-density parts as attribute `"components"`. Totals are unchanged.
* The shared baseline is documented as the hazard at the training means of
  `X_cox`. Linear predictors are compared up to a constant.
* `gemcox_simulate()`: the simulation pipeline's data-generating mechanism,
  moved here so there is one definition.
* `gemcox_selftest()`: fast engine checks run by the simulation pilot.

## Not ported

* The cross-validated C-index selectors (`cv_select_K_gamma*`,
  `cv_select_v5`, `gemcox_pipeline_v*`): predictive CV selects `gamma` away
  and is not appropriate for choosing K.
* `select_lambda()` (`lambda = NULL`), `heldout_loglik()` (disabled in
  PATCHED), `predict_risk()`, the preprocessing helpers
  (`preprocess_*_v4`, `screen_features_univariate`) and the CVIA078-specific
  helpers.
