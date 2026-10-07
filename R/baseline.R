## baseline.R -- Breslow baselines, get_cumhaz_at(), and the survival
## log-density log f_k(t, d | x) built on them.
##
## A baseline is a step function list(time, jump, cumhaz): event times,
## hazard jumps at those times, and their cumulative sum.

empty_baseline <- function() {
  list(time = numeric(0), jump = numeric(0), cumhaz = numeric(0))
}

# Apply the jump / cumulative-hazard caps, warning when either binds.
cap_baseline <- function(uniqT, jump, control) {
  check_guard(jump > control$max_jump, "max_jump")
  jump <- pmin(jump, control$max_jump)
  cumhaz <- cumsum(jump)
  check_guard(cumhaz > control$max_cumhaz, "max_cumhaz")
  list(time = uniqT, jump = jump, cumhaz = pmin(cumhaz, control$max_cumhaz))
}

# Floor the risk-set denominator. It only changes a jump when there are
# weighted events at that time (num > 0); with num = 0 the jump is 0 either
# way, so only the former is reported.
floor_denominator <- function(denom, num, control) {
  check_guard(denom < control$denom_floor & num > 0, "denom_floor")
  pmax(denom, control$denom_floor)
}

# Numerators and risk-set denominators of the Breslow estimator at the
# unique event times: num_m = sum of event weights at t_m, denom_m = sum of
# risk weights over subjects with time >= t_m. Vectorised with suffix sums
# over the sorted times; this replaced an O(n x events) loop (the result is
# unchanged to rounding, checked by test-equivalence.R).
breslow_terms <- function(time, status, event_w, risk_w) {
  ev <- status == 1
  uniqT <- sort(unique(time[ev]))
  num <- as.numeric(rowsum(event_w[ev], time[ev], reorder = TRUE))
  o <- order(time)
  suffix <- rev(cumsum(rev(risk_w[o])))
  denom <- suffix[match(uniqT, time[o])]
  list(time = uniqT, num = num, denom = denom)
}

# Weighted Breslow estimator for one cluster (baseline = "cluster").
breslow_weighted <- function(time, status, eta, w, control = gemcox_control()) {
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  status <- as.integer(status)
  eta <- clamp_eta(as.numeric(eta), control$eta_clamp, "eta_clamp")

  if (sum(w * status) < 1e-10) return(empty_baseline())

  bt <- breslow_terms(time, status, event_w = w, risk_w = w * exp(eta))
  cap_baseline(bt$time, bt$num / floor_denominator(bt$denom, bt$num, control), control)
}

#' Shared weighted Breslow baseline
#'
#' Breslow estimator of one baseline hazard shared by all clusters: at each
#' event time the jump is (number of events) / sum over the risk set of
#' sum_k tau_ik exp(eta_ik).
#'
#' @param time,status Follow-up time and event indicator (0/1).
#' @param eta_list List of K linear-predictor vectors, one per cluster.
#' @param tau n x K matrix of membership weights.
#' @param control Numerical guards, from [gemcox_control()].
#' @return A step function `list(time, jump, cumhaz)`.
#' @keywords internal
#' @noRd
shared_breslow <- function(time, status, eta_list, tau, control = gemcox_control()) {
  status <- as.integer(status)
  if (!any(status == 1)) return(empty_baseline())

  risk_w <- numeric(length(time))
  for (k in seq_len(ncol(tau))) {
    risk_w <- risk_w +
      tau[, k] * exp(clamp_eta(as.numeric(eta_list[[k]]), control$eta_clamp, "eta_clamp"))
  }
  bt <- breslow_terms(time, status, event_w = status, risk_w = risk_w)
  cap_baseline(bt$time, bt$num / floor_denominator(bt$denom, bt$num, control), control)
}

# Read cumulative hazard values from a step-function baseline.
#
# Known defect fixed here: the original ended with
#   ifelse(idx > 0, baseline$cumhaz[idx], 0)
# When idx contains zeros (times before the first event time), x[0] silently
# drops elements and ifelse() recycles the shorter vector, so every subject
# after the first pre-grid time received another subject's cumulative
# hazard. That attenuated every shared-baseline fit roughly 2x. Only valid
# positions are indexed. This is the ONLY place a cumulative hazard is read.
get_cumhaz_at <- function(baseline, times) {
  tt <- as.numeric(times)
  if (length(baseline$time) == 0) return(rep(0, length(tt)))
  idx <- findInterval(tt, baseline$time)
  out <- numeric(length(tt))
  ok <- idx > 0
  out[ok] <- baseline$cumhaz[idx[ok]]
  out
}

# Baseline hazard jump for each time: the jump at the last grid time <= t.
# A time before the first grid time uses the first jump (only possible for
# held-out subjects; training event times always lie on the grid).
get_jump_at <- function(baseline, times) {
  idx <- findInterval(as.numeric(times), baseline$time)
  baseline$jump[pmax(idx, 1L)]
}

# log f_k(t_i, d_i | x_i) = d_i (log dLambda0(t_i) + eta_i) - Lambda0(t_i) exp(eta_i)
log_surv_density <- function(time, status, eta, baseline) {
  status <- as.integer(status)
  n <- length(time)
  if (length(baseline$time) == 0) return(rep(0, n))

  Lambda <- get_cumhaz_at(baseline, time)
  haz0 <- rep(1, n)
  ev <- which(status == 1)
  if (length(ev) > 0) haz0[ev] <- pmax(get_jump_at(baseline, time[ev]), 1e-12)

  status * (log(haz0) + eta) - Lambda * exp(eta)
}
