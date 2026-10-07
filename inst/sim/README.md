# GeM-Cox simulation pipeline

Everything the manuscript reports comes from `results/`. Figures are made
from saved results, never from a refit.

## Layout

```
00_core.R          config, estimators, metrics, parallel harness (DGP is gemcox_simulate())
check_parallel.R   serial vs parallel identity check
01_pilot.R         20-replicate check of four expectations recorded in advance
02_experiments.R   E1-E5 and the supplement S1
03_summarise.R     summaries with Monte Carlo SEs, and figures
04_posthoc.R       POST HOC diagnostic added after the pilot (not pre-specified)
C_common.R         competitor experiments C1-C5: designs and dataset hashes (after sim-freeze-v2)
05_hash_manifest.R verifies that C1-C4 reuse the frozen datasets (hashes vs the frozen commits)
05_competitor.R    competitor experiments C1-C5
06_competitor_summarise.R  competitor tables, verdicts and figure
C_methods.R        estimators added after sim-freeze-v2 (shared by 05 and 07)
07_convergence_rerun.R     E1-E4a and E5c at convergence (R1, R2)
08_convergence_summarise.R frozen claims at convergence; matched regularisation
09_tolerance_study.R       evidence for a proposed default tolerance (not applied)
run_cell.R         one cell of a registered experiment, checkpointed (cluster runs)
aggregate_cells.R  cell files -> the per-experiment results file the summaries read
slurm/             SLURM array scripts: time_cells.sh, submit.sh, run_task.sh, status.sh
11_E6_summarise.R  E6 analysis: events vs n (pre-registered; added after sim-freeze-v2)
results/           raw per-replicate rows, one .rds per experiment
figures/           manuscript figures
```

## Run order

```sh
R CMD INSTALL gemcox          # from the repository root
cd gemcox/inst/sim
Rscript check_parallel.R      # serial and parallel must agree
Rscript 01_pilot.R
Rscript 02_experiments.R      # GEMCOX_RUN="E1,E2b" to select; GEMCOX_REPS=20 to rehearse
Rscript 04_posthoc.R          # post hoc diagnostic
Rscript 03_summarise.R
# added after sim-freeze-v2
Rscript 05_hash_manifest.R
Rscript 05_competitor.R
GEMCOX_ARMS=convergence Rscript 05_competitor.R
Rscript 06_competitor_summarise.R
GEMCOX_MANIFEST_SET=R Rscript 05_hash_manifest.R
Rscript 07_convergence_rerun.R
GEMCOX_ARMS=matched Rscript 05_competitor.R
Rscript 09_tolerance_study.R
Rscript 08_convergence_summarise.R
```

`GEMCOX_WORKERS` sets the number of PSOCK workers (default: cores - 1).

## Design decisions

**One data-generating mechanism.** `gemcox::gemcox_simulate()`. `u` (profile
direction) and `v` (coefficient direction) are orthonormal, and the base
coefficient is orthogonal to `v`. So `mu_sep` controls only how far apart the
feature profiles sit, `beta_sep` controls only how different the survival
mechanisms are, and the two clusters have equal coefficient norms.

**Fixed estimator settings.** gamma = 1, ridge lambda = 0.05, alpha = 0,
`normalize_gmm_by_dim = FALSE`, event rate 0.55, k-means initialisation
with seed 1.

**Comparators.**
- *Single Cox model:* `coxph` on all features.
- *Two-stage:*
  - The first stage is `mclust::Mclust(G = 2, modelNames = "VVI", prior =
    priorControl())`, a diagonal, varying-volume mixture with mclust's
    default conjugate prior. It is followed by `coxph` within each hard
    cluster.
  - The same model is used in every cell of every experiment. It matches
    GeM-Cox's regularised diagonal covariance, and it keeps the comparator
    one method along the whole p-axis.
  - Without the prior, VVI returned no fit at p >= 40 (n = 800). With it,
    VVI fitted 10/10 datasets at p = 40 and 80 for n = 400, 800 and 1600.
  - EEI in every cell was the pre-agreed fallback if VVI with the prior
    still failed at p >= 40; it was not needed.
  - A rare remaining failure (a mixture proportion estimated as zero; 1 of
    240 test datasets, at n = 400, p = 20) is recorded as a failed fit.
- *Oracle:* `coxph` within the true clusters.

**The estimand is the contrast, not membership.** Recovery is the mean
cosine of the principal angles between the span of the fitted centred
coefficients and the truth. At K = 2 that is the absolute cosine between
contrast vectors.

**Primary metric (fixed after the pilot, before the full run).**
- Paired per-replicate differences in recovery: GeM-Cox (gamma = 1) minus
  each comparator, computed on the same simulated dataset.
- MCSE = SD of the paired differences / sqrt(number of pairs).
- Normalised version: (gamma1 - gamma0) / (oracle - gamma0), a ratio of
  paired means with a delta-method MCSE.

**The random-direction floor sqrt(2 / (pi p)) is descriptive context only.**
The estimators' null behaviour is not isotropic. In the pilot, the
two-stage contrast aligned with the shared coefficient, which the design
makes orthogonal to the true contrast, and fell below the floor. GeM-Cox
with gamma = 0 rose slightly above it at beta_sep = 3. Percentages "of
oracle above the floor" are reported only as description.

**P2 depends on replication.** At pilot size the criterion recorded in
advance failed (CV 0.60). At 200 replicates it is met (CV 0.19 < 0.5).
Even so, the gain divided by beta_sep^2 still declines: 0.042, 0.036,
0.029 at beta_sep = 1, 2, 3 (MCSE 0.009, 0.003, 0.002), and it is 0.111 at
beta_sep = 0.5. E1 was not modified.

**E2b (the final experiment; pre-specified, with its decision rule,
before it was run).**

*(a) Adaptive rule vs normalize = FALSE on E2's 12 cells.*
- Uses E2's datasets, 200 replicates.
- normalize = TRUE is an auxiliary arm. gamma = 0 and the oracle give the
  normalised gain.
- **Decision rule for the package default:** the adaptive rule replaces
  FALSE if and only if both hold:
  1. in no cell is it worse than FALSE by more than 2 MCSE of the paired
     difference;
  2. the paired difference averaged over the 12 cells exceeds 2 MCSE of
     that average.

  Degenerate fits are excluded.

*(b) LRT power under the adaptive rule.*
- The LRT uses the bootstrap null, B = 19, on E5a/E5c's cells and
  datasets, so it pairs with E5c's LRT under FALSE.
- *Added:* the two type I cells, because the rule selects TRUE on K = 1
  data, a setting E5c did not calibrate.
- The rule is applied once to the observed data and held fixed across
  bootstrap replicates, which keep the features fixed.
- Under TRUE the statistic is a likelihood ratio of the tempered objective;
  the bootstrap calibrates it.

*(c) CVIA078 scale.*
- n = 117, event rate 0.44, p = 9, beta_sep = 2, mu_sep in {0, 0.5}.
- LRT (bootstrap, B = 19) under FALSE and under the adaptive rule.
- *Added:* a K = 1 null cell at the same scale (200 replicates).

## Added after sim-freeze-v1

These experiments were added after the study was frozen at
`sim-freeze-v1`. Nothing frozen is changed: not the DGP for K_true <= 2,
the estimator settings, the metrics, E1-E5c, or their results files. Each
experiment below was written here and committed before it was run.

E5b's K grid was {1, 2}, so it could not select K > 2.

**E5d: type I error when profiles separate but mechanisms do not.**
- *Question:* does the default test (LRT, bootstrap null, B = 19,
  normalize = FALSE) reject "one survival mechanism" merely because
  feature profiles differ?
- *Design:* profile-only null, K_true = 2 with beta_sep = 0 (both clusters
  have the base coefficient) and mu_sep in {0.5, 1, 2}; n = 400, p = 10,
  500 datasets per cell. Plus rho = 0.5 at mu_sep = 1 (200 datasets).
  Seeds from experiment id 24.
- *Report:* rejection at 0.05 with MCSE, p-value histogram and KS test,
  alongside the E5a/E5c type I rows.
- *Expectation:* rejection within 2 MCSE of 0.05 in every cell.
- *Failure criterion:* a cell fails if its rejection rate exceeds
  0.05 + 2 sqrt(0.05 x 0.95 / valid datasets), with the MCSE evaluated at
  the nominal rate: 0.0695 for 500 datasets, 0.081 for 200. If any cell
  fails, the test is not valid for data with profile structure. This will
  be stated prominently, with no fix attempted in this step.

**E5e: selecting K, and how accurately.**
- *Question:* how accurately can K be selected, by which criterion, and
  where does selection fail?
- *Design:* 200 datasets per cell, seeds from experiment id 25, K grid
  1:4.
  - K_true = 1: n in {400, 800}, p = 10.
  - K_true in {2, 3}: regime (profile only, mechanism only, both) x
    separation (moderate: mu_sep = 1, beta_sep = 2; strong: mu_sep = 3,
    beta_sep = 3, applied to whichever components the regime separates)
    x n in {400, 800}; p = 10.
  - K_true = 3 uses the triangle geometry of `gemcox_simulate(K_true = 3)`
    (added after sim-freeze-v1), so K = 3 is irreducible.
  - CVIA078 scale: n = 117, event rate 0.44, p = 9, K_true in {1, 2},
    "both" regime at moderate separation.
- *Selectors:*
  1. **Joint cross-validated criterion**: held-out GMM log-density plus
     cross-validated partial likelihood, 5 folds, as in E5b; K = argmax.
  2. **Partial criterion**: cross-validated partial likelihood only, from
     the same folds and fits; K = argmax.
  3. **BIC of the full-data mixture fit**: -2 l_K + df_K log n with
     df_K = (K - 1) + 2 K q + K p. That counts mixing proportions,
     diagonal Gaussian means and variances, and cluster coefficients. The
     shared Breslow baseline is identical across K and is excluded.
     K = argmin.
  4. **Sequential bootstrap LRT**: test 2 vs 1 with
     `gemcox_heterogeneity_test(K0 = 1)`; if p <= 0.05, test 3 vs 2
     (`K0 = 2`); and so on. Select the first K0 whose test is not
     rejected, or 4 if all three are. B = 19 per test (bootstrap seed:
     dataset seed + K0). If a step fails, select the current K0 and flag
     the failure.
- *Report:* confusion matrices P(K_hat = k | cell) with MCSEs for every
  selector and cell; P(K_hat >= 2) by regime, as the question each
  selector answers; the sequential LRT's first-step rejection rates;
  failure counts.
- *Expectations (from the brief, recorded before running):*
  - Mechanism-only rows: near-zero accuracy for every selector. This
    would confirm the two-scale result, not be a failure.
  - Profile-only rows: the joint criterion and BIC select the number of
    profile clusters; the partial criterion and the sequential LRT should
    not report distinct mechanisms.
  - Moderate separation with K_true = 3: unknown; reported as found.
- *Cost:* one dataset per cell, timed under full load, gave about 4.8
  hours with all selectors on 200 datasets per cell (3.4 hours with the
  sequential LRT on 100). That is below the 8-hour threshold, so all
  selectors run on 200 datasets per cell (`GEMCOX_SEQ_REPS = 200`).

**E5f (optional): where the adaptive test becomes calibrated.**
- *Question:* the adaptive rule roughly doubled the LRT's power at
  220-440 events but was anti-conservative at about 52 events (0.135
  under the null, E2b(c)). Where does it become calibrated?
- *Design:* n in {200, 300}, event rate 0.44, p = 10, mu_sep = 0. Null
  cells are K_true = 1 with a nonzero shared coefficient, as in E2b(c);
  power cells are K_true = 2, beta_sep = 2. 200 datasets per cell, seeds
  from experiment id 26.
- *Tests:* the adaptive and default (normalize = FALSE) LRTs on the same
  datasets, bootstrap null, B = 19.
- *Report:* type I error and power for both, with MCSEs and KS tests,
  plus the E2b(c) n = 117 rows for context.
- *Use of results:* documentation only; the package default is not
  changed on the basis of this experiment. No expectation or decision
  rule is attached.

**E4c reporting.** The degeneracy rule applies as in every other
experiment. A configuration whose fits are all degenerate (the "sharpened"
arm collapses to one cluster) is reported separately. Recovery is reported
per configuration with its MCSE; no correlation across configurations is
computed.

**Degenerate fits.**
- *Definition:* a fitted cluster with fewer than 10 subjects or fewer than
  5 events. Hard labels are used: the assigned cluster, or argmax tau for
  GeM-Cox. A single Cox model is one cluster (the whole sample).
- *Rates:* reported for every method in every cell.
- *Primary analysis:* excludes degenerate fits; a paired difference is
  dropped if either member is degenerate.
- *Sensitivity analysis:* scores them at the random-direction floor, not
  at the zero-coefficient score `safe_cox()` gives them.
- *Two-stage gains:* paired gains over two-stage are reported under both
  scorings, and as originally scored.
- *Why:* at mu_sep = 0 the two-stage mixture step often produced a tiny
  cluster whose zeroed fit made its "contrast" equal the other cluster's
  coefficients.

**E1b (pre-specified after S1, reported separately from E1).**
- Uses E1's exact datasets (same seeds).
- *Arms:*
  - GeM-Cox gamma = 1 with `normalize_gmm_by_dim = TRUE`.
  - A regime-adaptive arm: TRUE when the profile part of the joint
    cross-validated criterion (held-out GMM log-density, 5 folds, fits at
    the default settings) selects K = 1, and FALSE otherwise. The overall
    joint-criterion choice is recorded as well.
- E1's methods are refitted on the same data as paired comparators.

**E5c (added after the full run).**
- The in-sample likelihood ratio 2 (l2 - l1) on the mixture
  log-likelihood (gamma = 1).
- Calibrated by the same parametric bootstrap as E5a (features fixed,
  times from the fitted K = 1 Cox model, Kaplan-Meier censoring) and by
  the permutation null, sharing one observed statistic; B = 19.
- Uses E5a's exact cells and datasets, and is reported alongside E5a.

**Post hoc diagnostic (04_posthoc.R, not pre-specified).**
- Question: does the two-stage contrast's alignment with the shared
  coefficient require mechanism heterogeneity?
- Design: mu_sep = 0, n = 400, p = 10, 200 replicates, two-stage and
  GeM-Cox (gamma = 0).
- The diagnostic cell is beta_sep = 0; beta_sep = 2 is a heterogeneous
  reference.
- Records the contrast norm and |cos(contrast, shared coefficient)|.

**Secondary metrics.**
- ARI.
- MSE of the linear predictor on an independent test set of 2000. The
  predictor is the uncentred mixture `sum_k tau_k x'beta_k`, with `tau` from
  `predict()`, compared with the true eta after removing the mean of the
  difference, since the overall constant is not identified.
- Posterior sharpness.

**Inference (E5a).**
- The K = 2 vs K = 1 test uses the cross-validated partial likelihood, with
  B = 19 null replicates and 3 folds.
- The observed statistic is computed once per dataset and compared against
  both the parametric-bootstrap null (the package default) and the joint
  permutation null.
- Type I error is evaluated on K = 1 data with a nonzero shared beta.

**Raw rows are saved; summaries are computed separately**, with every
dispersion taken before its column is replaced by a mean.

**Failures are recorded, not dropped.** A method that errors gives a row
with `failed = TRUE` and the error message.

**Seeds** are deterministic in (experiment, cell, replicate).
- *Generator:* every random draw uses the L'Ecuyer-CMRG generator, the
  one `future.seed = TRUE` installs in workers. Each per-dataset function
  pins it, so any replicate reproduces exactly from any R session.
- *Per-method seeding:* the stream is re-seeded from the dataset seed
  before every method, because mclust initialises from a random subset
  when n > 2000.
- *Caveat:* `gemcox_simulate(seed = s)` called in a session using R's
  default Mersenne-Twister generator gives different data from the
  pipeline's; call `RNGkind("L'Ecuyer-CMRG")` first.

**Parallelism.**
- Workers are PSOCK (multisession) processes with `OMP_NUM_THREADS = 1`.
- `future.seed = TRUE` makes results identical to a serial run, which
  `check_parallel.R` verifies.
- Every results file stores the package version, git commit, configuration
  and `sessionInfo()`.

## Added after sim-freeze-v2: latent-class Cox competitor (C1-C5)

Written here and committed before any C experiment was run. Nothing frozen
changes: not the DGP (the new feature option is off by default, and the
datasets were verified identical; see below), the estimator settings, the
metrics, E1-E5f, or their results files and rows.

**The competitor: latent-class Cox with logistic gating** (a mixture of
experts; `gemcox:::lcc_fit()`, internal).
- Membership P(Z = k | x) = softmax_k(a_k0 + a_k' x) on the clustering
  features; experts h_k(t | x) = lambda0(t) exp(beta_k' x); E-step
  tau_ik proportional to P(Z = k | x_i) f_k(t_i, d_i | x_i).
- It differs from GeM-Cox only in the membership model: logistic in x
  (discriminative) instead of a Gaussian mixture for x (generative).
  Everything else is GeM-Cox's own code: the tau-weighted Cox M-step with
  the shared Breslow baseline, the same ridge on beta (lambda = 0.05,
  alpha = 0), the E-step survival density, k-means initialisation with
  seed 1, one start, the convergence rule (relative change in the
  log-likelihood below 1e-4, at most 100 iterations) and the guards.
- Gating: glmnet multinomial, alpha = 0, unpenalised intercept, gating
  penalty equal in number to the beta penalty (0.05). glmnet normalises
  the weights, so the gating ridge is n lambda / 2 sum_k |a_k|^2 while
  cluster k's Cox ridge is n_k lambda / 2 |beta_k|^2: the same number, not
  the same effective strength.
- New subjects are assigned by the gating. The test-set predictor is the
  uncentred mixture sum_k P(Z = k | x) x' beta_k, as for GeM-Cox.
- Verified in `tests/testthat/test-lcc.R`.

**Both membership models are correctly specified in C1-C4, and both are
misspecified in C5.**
- Under the frozen Gaussian DGP (rho = 0 in every C cell), X | Z = k is
  N(-/+ (mu_sep / 2) u, I) with equal proportions, so
  log P(Z = 2 | x) / P(Z = 1 | x) = mu_sep u' x, exactly linear in x. The
  competitor's gating is correctly specified, and so is GeM-Cox's
  diagonal Gaussian mixture. At mu_sep = 0 the truth is constant
  membership, which both models contain. **C1-C4 therefore compare
  generative against discriminative membership models when both are
  correct.**
- In C5 the features are standardised log-normal about the same cluster
  means. The Gaussian mixture is misspecified. The true log-odds is not
  linear in x: it involves log(x - shift) terms and is infinite near the
  lower edge of each cluster's support. So the logistic gate is
  misspecified too. **C5 compares the two when both are misspecified.**

**Data: the frozen datasets, verified by hash.**
- C1-C4 reuse frozen datasets: seeds come from the frozen experiment id
  and cell index.
- `05_hash_manifest.R` installs the package at the commit that produced
  each frozen results file (E1: de83301; E2, E3a, E3b: f965278; E2b(c):
  30595fa). In a separate R session it sources that commit's `00_core.R`
  and regenerates every training set from the seeds stored in the frozen
  rows. For C1-C3 it regenerates every test set as well. It compares event
  counts with the stored ones, then generates the same datasets with the
  current code, exactly as `05_competitor.R` does.
- Result, before any C run: all 5,200 datasets identical (seed, event
  count, md5 of X, time, status, Z and eta; test sets too for C1-C3).
  Saved as `results/C0_dataset_hashes.rds`.
- Every C row records the md5 of its training and test data. The summary
  computes no paired difference for a replicate whose hashes, seed or
  event count do not match the manifest.
- On the 24 timed datasets of C1-C3, refitting GeM-Cox (gamma = 1)
  reproduced the stored recovery exactly.

**Experiments.** 200 datasets per cell, rho = 0 throughout.

| id | cells | datasets | methods fitted |
|---|---|---|---|
| C1 | E1's 4: n 400, p 10, mu_sep 0, beta_sep {0.5, 1, 2, 3} | E1's (exp 1) | competitor and the lambda / 10 arm; paired with E1's rows |
| C2 | E2's 12: n 800, p 10, mu_sep {0, 0.5, 1, 2} x beta_sep {1, 2, 3} | E2's (exp 2) | competitor; lambda / 10 arm in cells 6 and 7 (mu_sep 0.5 and 1, beta_sep 2); paired with E2's rows |
| C3 | E3a's 4 (n 200-1600) and E3b's 4 (p 10-80); mu_sep 0.5, beta_sep 2 | E3a's (exp 3), E3b's (exp 31) | competitor; paired with the frozen rows |
| C4 | E2b(c)'s two power cells: n 117, p 9, event rate 0.44, beta_sep 2, mu_sep {0, 0.5} | E2b(c)'s (exp 23, cells 2 and 3) | all five frozen methods and the competitor (E2b(c) stored test results only) |
| C5 | new: n 800, p 10, mu_sep {0.5, 1}, beta_sep 2, log-normal features | new (exp 106) | all five frozen methods and the competitor |

E2b(c)'s K = 1 null cell is not run: without a contrast, recovery is
undefined.

**C5 feature model.** `gemcox_simulate(features = "lognormal", sdlog =
0.8)`, a new argument whose default (`"gaussian"`) is the frozen DGP.
- Each coordinate of the Gaussian noise z is replaced by
  (exp(0.8 z) - e^0.32) / sqrt((e^0.64 - 1) e^0.64): mean 0, variance 1,
  skewness 3.69.
- The cluster means, directions, coefficients, survival and censoring
  models are unchanged, and the same random numbers are drawn.

**Metrics** (the frozen definitions).
- Contrast recovery, ARI, and MSE of the linear predictor on the 2000
  test subjects (uncentred mixture, mean difference removed).
- Posterior sharpness and degeneracy (a cluster with fewer than 10
  subjects or fewer than 5 events, hard labels = argmax tau).
- Failure status, iterations, convergence, guard warnings, and whether
  the gating fit succeeded.

**Primary comparison.**
- Per cell, the paired per-replicate difference in recovery: GeM-Cox
  (gamma = 1) minus the competitor. MCSE = SD of the differences /
  sqrt(pairs).
- Degenerate fits are excluded (a pair is dropped if either member is
  degenerate), as in every frozen experiment. The floor-scored
  sensitivity is reported alongside.
- Also reported:
  - competitor minus two-stage and competitor minus GeM-Cox (gamma = 0),
    paired the same way;
  - the competitor's normalised gain (competitor - gamma0) /
    (oracle - gamma0), next to GeM-Cox's;
  - ARI, MSE and sharpness, with MCSEs;
  - degeneracy and failure rates for every method.

**Verdict rule.**
- *Per cell:* GeM-Cox **beats** the competitor if the paired difference
  exceeds 2 MCSE, **loses** if it is below -2 MCSE, and **matches**
  otherwise.
- *Per regime:* **beats** if it beats in at least one cell and loses in
  none; **loses** if it loses in at least one cell and beats in none;
  **matches** if every cell matches; **mixed** otherwise.
- *Regimes:*
  - mechanism only: C1, and C2's mu_sep = 0 cells;
  - profile separation: C2's mu_sep > 0 cells;
  - sample size and dimension: C3;
  - application scale: C4;
  - non-Gaussian features: C5.

  All but C5 have Gaussian features.
- *Multiplicity:* the 2-MCSE threshold is applied per cell, without
  correction. Across 28 cells, one or two spurious verdicts can occur by
  chance where the true difference is zero. Verdicts are reported next to
  the differences.

**Expectations (from the brief).**
1. *Profile separation, Gaussian features (C2, mu_sep > 0):* GeM-Cox
   (gamma = 1) recovery is at least the competitor's. The expectation
   fails in a cell if the paired difference is below -2 MCSE, and fails
   overall if it fails in any of the 9 cells.
2. *Mechanism only (C1):* genuinely unknown; reported as found.
3. *Non-Gaussian features (C5):* genuinely unknown. The competitor may do
   better.
4. *CVIA078 scale (C4):* both near the floor. Operationalised
   descriptively: both methods' mean recovery is closer to the
   random-direction reference sqrt(2 / (9 pi)) = 0.266 than to the
   oracle's. Degeneracy and failure rates are reported for both; they may
   differ.

- *Flag:* if the competitor beats GeM-Cox (difference below -2 MCSE) in
  any cell with Gaussian features (C1-C4), this is stated at the top of
  the report section, because it changes the paper's claims.

**Sensitivity arm 1: gating penalty lambda / 10.**
- Gating penalty 0.005, on C1's 4 cells and C2's cells 6 and 7, same
  datasets.
- Verdicts are recomputed with this arm in those cells.
- Reported: whether any per-regime verdict changes. The regime verdict is
  computed over the cells where the arm ran, for both arms.

**Sensitivity arm 2: convergence (added after timing, before any C run).**
- *Why.* Timing found that under the shared rule (relative change below
  1e-4) the competitor often stops after one EM step.
  - From the k-means start, its gating is fitted to hard, linearly
    separable labels, so the first step barely changes the likelihood.
  - In 20 datasets per cell (iterations and time recorded, recovery not),
    it stopped within 3 iterations in 24% of datasets (0-70% by cell).
    GeM-Cox did so in none of the frozen E2 fits.
  - At tolerance 1e-8 the competitor needs a median of 39-326 iterations
    and GeM-Cox (gamma = 1) 28-258, against GeM-Cox's frozen 4-18. The
    frozen rule stops both early, and the competitor much more often.
- *Disclosure.* The diagnostic that found this also computed recovery on
  four datasets (replicate 1 of C1 cell 4, C2 cells 3 and 5, and C3a
  cell 4). In each, the competitor's recovery rose substantially with the
  tight tolerance.
- *Design.* The competitor and GeM-Cox (gamma = 1), both with tolerance
  1e-8 and at most 1000 iterations, all other settings unchanged, on
  every C cell's datasets. None of the 84 fits of each method timed at
  this tolerance needed more than 1000 iterations. These are new arms,
  saved in `*_convergence.rds`; the frozen rows are not touched.
- *Report.*
  - Paired GeM-Cox (tight) minus competitor (tight) per cell, with
    verdicts, and whether any per-regime verdict differs from the
    primary.
  - Descriptive: competitor (tight) minus competitor (primary), and
    GeM-Cox (tight) minus GeM-Cox (frozen), per cell. If GeM-Cox's
    recovery moves materially with the tight tolerance, the report says
    so as a limitation of the frozen study.
  - No decision rule is attached. **The primary comparison remains the
    one specified in the brief.**

**Cost.** One dataset per cell, timed serially.
- Primary arms: 11 s per set of 28 cells, 0.6 CPU hours for 200 datasets
  per cell. Projected at about 20 minutes on 7 workers, given the 3-4x
  slowdown under parallel load seen before.
- Convergence arm (mean of 3 datasets per cell): 174 s per set, 9.7 CPU
  hours. Projected at 4-5.5 hours.
- A 3-replicate rehearsal of both arms and of `06_competitor_summarise.R`,
  in a scratch directory and then discarded, checked that the code runs
  and that the hash verification passes.

**Interpretation of the convergence arm (committed while it was running,
before any of its results files were written).**
- *At the time of this commit* the arm's progress log had printed the
  tight competitor's mean recovery in C1 cells 1 and 2 (0.310, 0.428). No
  GeM-Cox (tight) result and no paired comparison had been produced.
- *Operationalised* with the verdict rule above, applied to GeM-Cox
  (tight) minus the competitor (tight):
  - the competitor *still beats* GeM-Cox at convergence if GeM-Cox loses
    (difference below -2 MCSE) in at least one Gaussian cell (C1-C4);
  - the advantage *disappears* if it loses in none.

  The number of such cells and their regimes are reported with the
  branch, because a single cell can cross 2 MCSE by chance (see
  Multiplicity).
- **If the competitor still beats GeM-Cox at convergence in Gaussian
  cells**, GeM-Cox's generative membership model is inferior in these
  settings. The paper is reframed around outcome-informed mixture-of-Cox
  models and recommends the better-performing membership model.
- **If the advantage disappears at convergence**, the primary-arm result
  reflected early-stopping regularisation. Proceed to the
  matched-regularisation comparison below.
- **In either case**, C5 is reported as a limitation of Gaussian-mixture
  membership under skewed features.
- *Provenance note:* the convergence results files will record this
  commit, not 9b57ed8, because `run_meta()` reads HEAD when a file is
  saved. The two commits differ only in this README, so the code that
  runs is the same.

**Matched-regularisation comparison (pre-registered; not yet run).**
- *Question:* with each method's regularisation varied over a small
  grid, and both run to convergence, which membership model recovers the
  contrast better?
- *Configurations:*
  - GeM-Cox (gamma = 1): `normalize_gmm_by_dim` in {FALSE, TRUE} x `temp`
    in {1, 2, 5}: six configurations. The default is FALSE, 1.
  - Competitor: gating penalty in {lambda, lambda / 10, lambda x 10} =
    {0.05, 0.005, 0.5}: three configurations. The default is lambda.
  - All other settings as in the primary arms. EM tolerance 1e-8 and at
    most 1000 iterations, as in the convergence arm, so that early
    stopping does not act as an uncontrolled regulariser. (Tolerance
    chosen here; the request did not specify one.)
- *Cells and data:* C1 (4 cells), C2 at mu_sep in {0.5, 1} x
  beta_sep = 2 (cells 6 and 7), and C5 (2 cells); the same
  hash-verified datasets, 200 per cell.
- *Primary:* default versus default, i.e. the convergence arm's
  comparison in these cells, with the same verdict rule per cell and per
  regime.
- *Oracle upper bound: best versus best by true recovery.*
  - In each cell, each method's configuration with the highest mean
    recovery (degenerate fits excluded) is selected.
  - The paired difference between the two selected configurations is
    reported, with its MCSE and the same verdict rule.
  - This uses the truth, so it is an upper bound, not an achievable
    procedure. Selection on the same replicates biases it upward, more
    for GeM-Cox (six configurations) than for the competitor (three).
- *Reported for every configuration in every cell:*
  - recovery, ARI, MSE of the linear predictor and posterior sharpness,
    with MCSEs;
  - degeneracy, failure, convergence and iteration counts.
- *Cost (not yet timed):* about 5 GeM-Cox and 2 competitor fits per
  dataset beyond the convergence arm's two, 8 cells x 200 datasets.
  Before running, one dataset per cell will be timed and the projection
  reported.

## Added after sim-freeze-v2: the core experiments at convergence (R1, R2)

Written here and committed before any R run, after the convergence arm
showed that the frozen EM rule stops GeM-Cox early. At tolerance 1e-8,
GeM-Cox's recovery was up to 0.35 higher on the same datasets. Nothing
frozen changes; the frozen rows stay as they are.

**Purpose: establish which frozen claims hold for converged estimators,
especially GeM-Cox gamma = 1 versus gamma = 0.**

**R1: E1, E2, E3a, E3b and E4a at convergence.**
- *Methods:* GeM-Cox gamma = 0, GeM-Cox gamma = 1 and the competitor, all
  to EM tolerance 1e-8 with at most 1000 iterations, other settings
  frozen. Single Cox, two-stage and the oracle are the frozen rows,
  unchanged.
- *Data:* each frozen experiment's datasets and test sets, 200 per cell.
  - E1-E3b are C1-C3b's datasets, verified in `C0_dataset_hashes.rds`.
  - E4a's 800 datasets and test sets are verified in
    `R0_dataset_hashes.rds` against the commit that produced
    `E4a_two_scales.rds` (de83301): all identical.
- *Reuse:*
  - For E1-E3b, gamma = 1 and the competitor at tolerance 1e-8 were
    already fitted on these datasets, with this code, in the convergence
    arm. Those rows are reused; gamma = 0 is fitted here.
  - The timing run refitted both on replicate 1 of all 24 cells, and each
    reproduced the stored recovery and iteration count exactly.
  - For E4a all three methods are fitted here.
- *Cap hits:* every fit records whether it stopped at 1000 iterations.
  Counts are reported by method and cell.
- *Comparisons (per cell; same metrics, pairing and degenerate-fit rule
  as the frozen study; verdict rule as for C):*
  - gamma = 1 minus gamma = 0;
  - gamma = 1 minus two-stage, and minus single Cox;
  - the normalised gain (gamma1 - gamma0) / (oracle - gamma0), using the
    converged gamma = 0;
  - the competitor against each of these;
  - each estimator, converged minus frozen.

  Each comparison is shown next to its frozen value.

**Claims evaluated.** Each is computed by the same code on the frozen fits
and on the converged fits (`08_convergence_summarise.R`, table
`tab_R_claims.csv`). A claim *holds for converged estimators* if its
criterion is met with converged fits. REPORT.md will state, for every
claim, whether it holds at convergence, with both values.

| id | claim | source | criterion |
|---|---|---|---|
| A | gamma = 1 beats gamma = 0 | E1, E2, E3a, E3b (24 cells) | paired difference > 2 MCSE in every cell |
| B | gamma = 1 beats two-stage | E1, E2 | > 2 MCSE in every cell with pairs |
| C | gamma = 1 beats single Cox | E1, E2 | > 2 MCSE in every cell |
| D | P1: covariate-only methods at the floor, gamma = 1 above it | E1, beta_sep >= 1 | max covariate-only excess over sqrt(2 / (pi p)) < 0.05 and gamma = 1's > 0.05 (the pilot's rule; converged gamma = 0 counts as covariate-only) |
| E | P2: gain / beta_sep^2 roughly constant | E1, beta_sep >= 1 | CV < 0.5 |
| F | P3: labels at chance while the contrast is recovered | E4a, every n | ARI < 0.10 and recovery above the floor by > 0.05 in every cell |
| G | recovery rises with n (two scales) | E4a | n = 6400 minus n = 800 > 2 SE |
| H | gain over gamma = 0 grows with n | E3a | n = 1600 minus n = 200 > 2 SE |
| I | gain over gamma = 0 decays with p (roughly p^-1.15) | E3b | p = 80 minus p = 10 < -2 SE; log-log slope reported |
| J | at mu_sep = 0 gamma = 1's prediction error is no smaller than single Cox's | E2, mu_sep = 0 | no cell with MSE difference < -2 MCSE |
| K | the default test is calibrated (K = 1 null, nonzero beta) | E5c type I cells | rejection <= 0.05 + 2 sqrt(0.05 x 0.95 / N): 0.0695 (N = 500), 0.081 (N = 200) |
| L | power 0.30-0.57 at beta_sep = 2 | E5c power cells | descriptive: range, and cells whose rejection changes by > 2 MCSE (paired) |

The frozen E1 P2 value is recomputed with the full-run convention
(degenerate fits excluded), which is how REPORT.md reports it.

**R2: E5c at convergence.**
- *Test:* the in-sample LRT with the parametric bootstrap null, B = 19,
  seed as in E5c. Every fit inside the test, including the K = 1 model
  the bootstrap draws from, runs to tolerance 1e-8 with at most 1000
  iterations; the settings pass through `gemcox_heterogeneity_test(...)`.
  The permutation null is not rerun.
- *Data:* E5c's datasets (E5a's cells and seeds): 500 and 200 K = 1
  datasets for type I, and 200 per power cell. All 2,300 are verified in
  `R0_dataset_hashes.rds` against the commit that produced `E5c_lrt.rds`
  (f965278): all identical.
- *Report:*
  - rejection at 0.05 with MCSE, next to the frozen E5c value on the same
    datasets;
  - the paired change, with MCSE from the paired SD;
  - a KS test of the type I p-values against uniform (descriptive; the
    p-values are on a 1/20 grid);
  - fits that stopped at the cap.
- *Cap counts:* `gemcox_heterogeneity_test()` now returns `n_fits` and
  `n_not_converged` for the LRT. This is an additive change: the
  statistic and p-value are unchanged, and it is tested.
- The progress log prints no rejection rates.

**Matched-regularisation comparison: now a sensitivity analysis.**
- Pre-registered in 412f85f, and relabelled here as a sensitivity
  analysis.
- It runs regardless of which interpretation branch the convergence arm
  gave. The design is unchanged: configurations, cells, tolerance 1e-8,
  default-vs-default primary, best-vs-best by true recovery as an oracle
  upper bound, and posterior sharpness for every configuration.
- The default configurations (GeM-Cox FALSE, 1; competitor lambda) are
  the convergence arm's rows, reused. The seven others are fitted by
  `GEMCOX_ARMS=matched Rscript 05_competitor.R` and saved to
  `*_matched.rds`.
- Datasets are verified against C0 (C1, C2) and against C5's primary rows
  (C5).
- The per-cell and per-regime verdict rules are as for C.

**Proposed default tolerance (a diagnostic; not applied to the package).**
- *Design:* `09_tolerance_study.R`. GeM-Cox (gamma = 1) is fitted at
  tolerances 1e-4, 1e-5, 1e-6, 1e-7 and 1e-8, each with at most 1000
  iterations.
  - Datasets: the first 20 of each of E1's 4 cells, E2's 12, E3b's p = 40
    and 80, C4's 2 and C5's 2, and the first 10 of E4a's 4 cells.
  - The 1e-8 fit is the reference.
- *Recorded:* seconds, iterations, convergence, recovery, the
  log-likelihood gap and the coefficient distance to the reference.
- *Rule, fixed here:* propose the loosest tolerance at which
  - at least 99% of fits, pooled, have recovery within 0.01 of the
    reference, and
  - at least 95% do in every cell.

  With it, propose max_iter as the smallest of 200, 500 and 1000 within
  which at least 99.5% of fits at that tolerance converge. If no looser
  tolerance qualifies, 1e-8 is proposed. Timing is reported relative to
  the current default.
- A default change would affect every fit, the test and the
  cross-validation functions, so it is proposed, not applied.

**Order of runs.** R0 manifest (clean tree), then R1, R2, matched, and the
tolerance study.

**Cost.**
- One dataset per cell of every part was timed serially. Wall time is
  projected at 0.34 x the serial CPU hours (the convergence arm's actual
  ratio on 7 workers), up to 0.5x.

  | part | CPU hours | wall hours |
  |---|---|---|
  | R1 | 5.8 | 2.0-2.9 |
  | R2 (E5c) | 19.9 | 6.8-9.9 |
  | matched | 2.0 | 0.7-1.0 |
  | tolerance study | 0.3 | 0.1-0.2 |
  | total | 28.0 | 9.5-14 |

- None of the ten timed E5c tests hit the cap (41 fits each).
- A 2-replicate rehearsal of every script, in a scratch directory and
  then discarded, checked that the code runs and that every replicate
  passes verification.

**Cap-hit rule for summaries at convergence (registered 2026-10-01,
before any R1 summary).**
- *Primary:* includes every fit and reports the rate of fits stopped at
  the 1000-iteration cap, per method per cell.
- *Sensitivity:* excludes capped fits.
- Both are reported. This applies to R1, R2, the matched-regularisation
  arm and the talk figures (`talk/README.md`).

## Added after sim-freeze-v2: E6, does recovery follow events or sample size?

Written here and committed before E6 was run. Nothing frozen changes. The
DGP gains an administrative-censoring option, off by default; the frozen
datasets' hashes are unchanged.

**Question.** Does mechanism recovery follow the number of events, and
profile recovery the number of subjects?

**Expectation (recorded now):** mechanism recovery follows events;
profile recovery follows n. Any result that disagrees is reported.

**Design.** Seeds come from experiment id 27; 200 datasets per cell;
beta_sep = 2, K_true = 2, rho = 0; test sets of 2000.
- *Censoring:* administrative, at a fixed follow-up and with no other
  censoring (`gemcox_simulate(censoring = "administrative", followup = )`).
  - The follow-up for each (p, mu_sep, event rate) is the event-time
    quantile at that rate, from one reference sample of 200,000 (seed
    20261006). It is fixed before the run and does not depend on n
    (`E6_FOLLOWUP` in `C_common.R`).
  - It runs from 28.5 (10% events) to 168.9 (55% events).
- *Cells* (60; `E6_CELLS` in `C_common.R`):
  - **factorial:** n in {200, 400, 800, 1600} x event rate in {0.1, 0.2,
    0.35, 0.55} x mu_sep in {0, 1}, with p = 10 (32 cells);
  - **matched events:** event rate in {0.175, 0.275} at n in {400, 800,
    1600} x mu_sep in {0, 1}, with p = 10 (12 cells). Each pairs with
    (n / 2, twice the rate), so it has the same expected events and twice
    the subjects. With the factorial's own pairs (0.1 vs 0.2) there are
    matched pairs at every rate, from 40 to 440 expected events;
  - **p axis:** p in {20, 40} x n in {400, 1600} x event rate in {0.2,
    0.55} x mu_sep in {0, 1} (16 cells), so that log p is estimable.
- *Methods:*
  - GeM-Cox gamma = 1 and latent-class Cox, both run to tolerance 1e-8
    with a cap of 1000;
  - two-stage (mclust VVI with the prior, then coxph);
  - the oracle (true labels).

  The other settings are the frozen ones.

**Recorded for every fit:**
- events; events in each true subgroup and in the smaller one; events per
  coefficient (smaller subgroup's events / p); realised event rate;
- the frozen metrics (contrast recovery, ARI, MSE, sharpness,
  degeneracy, failure, iterations, cap hits);
- **profile recovery:** |cos| between the estimated and true
  between-subgroup mean-difference direction in X, for mu_sep = 1 only.
  The estimated direction is taken from:
  - GeM-Cox: the fitted Gaussian means;
  - latent-class Cox: the gating slopes, which point along mu_2 - mu_1
    under this DGP's identity covariance;
  - two-stage: the membership-weighted means;
  - the oracle: the true-label means.

**Analysis** (`11_E6_summarise.R`, written before the run).
- *Primary:* for each method x outcome (mechanism recovery, profile
  recovery) x mu_sep, a weighted least-squares fit over cells of
  logit(cell-mean recovery) on log(mean events) + log(n) + log(p). The
  weights are 1 / Var(logit mean), by the delta method. Only cells with at
  least 20 usable fits enter, and degenerate fits are excluded, as in
  every frozen summary.
- *Which term carries the effect:* the one whose 95% CI lies above 0,
  provided the other term's CI includes 0 or its estimate is under a
  third of the first's. If both terms are positive, it is "both"; if
  neither is, "neither". This is classified for every fit and compared
  with the expectation.
- *Matched-events pairs:* the difference in mean recovery, (2n, rate / 2)
  minus (n, rate), with SE sqrt(MCSE_a^2 + MCSE_b^2).
  - For mechanism recovery, the expectation is |difference| < 2 SE.
  - For profile recovery, a difference > 2 SE.
- *Sensitivity:*
  - degenerate fits scored at the random-direction floor, instead of
    being excluded. Low-event cells have many degenerate fits, so
    excluding them selects the easier datasets. Both analyses are
    reported, and any classification that changes is listed.
  - a replicate-level fit, with each fit's own event count.
- *Reported:*
  - the per-cell table of events, events per subgroup, events per
    coefficient, and failure, degeneracy and cap-hit counts;
  - the regression coefficients with CIs, the classifications, the
    matched pairs, and every disagreement with the expectation;
  - two figures: mechanism recovery against events, and profile recovery
    against n.

**After E6:** rewrite the scaling law and the information index in terms
of events. The user will supply their current definitions; the repository
does not define them.

**Cost.** One dataset per cell, timed serially on the laptop (on battery,
in Low Power Mode): 445 s per replicate set, 24.7 CPU hours for 200
replicates. On 7 laptop workers that projects to about 8-12 hours, and R1
took twice its projection. The run is checkpointed per cell
(`run_cell.R --exp E6`, then `aggregate_cells.R --exp E6`), so an
interruption loses at most one cell. It can run on the laptop or the
cluster (`slurm/submit.sh E6`).

## Reporting

Use ADEMP (Morris, White & Crowther, *Stat Med* 2019;38:2074-2102).

## Added after sim-freeze-v2: separate G/X contrast-recovery pilot (2026-10-01)

The user requested a new preliminary comparison with Fei et al.'s proportional-
baseline latent-class PH model. The preregistration is
`../../../benchmark/fei_pilot/SPECIFICATION.md`; the pre-execution audit is
`../../../benchmark/fei_pilot/AUDIT.md`. Both are committed before calibration,
smoke tests, timing or simulation. Core: n={120,300}, d={0,1.5,3}, b={0,.5,1},
200 replicates/cell; selected n=800 and seven one-factor sensitivity cells follow.
Primary estimand is label-invariant squared coefficient-contrast error; absent
contrasts have undefined cosine, never an assigned random-direction score.
This NEW DGP uses separate G and X and .475 target events (user-authorized
exception to frozen .55). Frozen code, metrics, settings and outputs stay intact.
No coefficient-heterogeneity p-values; null cells quantify spurious contrasts.
