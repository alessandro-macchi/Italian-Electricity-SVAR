## 09_conditional_forecast.R ---------------------------------------------------
## Conditional forecasts (thesis sec. 5.4) with every parameter FIXED at the
## 2006-2025 estimates: no VAR is re-estimated and no rotation redrawn on the
## 2026 data. Three exercises share one machinery:
##   A   ex post: one-step-ahead structural shocks over the 2026 months with
##       all four variables observed, and the path without `cf$attack_shocks`
##       from `cf$attack_month` on;
##   D1  TTF held at its origin level to `cf$horizon_end`, read three ways --
##       all shocks (minimum norm, Waggoner & Zha, 1999: the reference path),
##       gas-specific shocks only, industrial demand shocks only;
##   C1  the median-target gas-specific shocks of `cf$c1_window` replayed on
##       top of the reference path.
## Every path is the zero-residual forecast plus hd_contributions() of an
## h x K matrix of shocks, so the historical decomposition's convolution does
## the work here too. Bands are pointwise percentiles over the admissible
## rotations in ident$draws$b0 -- identification uncertainty only, no
## bootstrap. Every rotation enters, except that a single-shock D1 reading
## drops out where it is not numerically defined (drop_undefined()).
##
## Depends on load_all_series() (01), transform_data() (02), split_by_role()
## and monthly_dummies() (04), reduced_form_irf() and structural_irf() (05),
## extract_var_coefs() and simulate_var() (06), hd_contributions() (08).
## -----------------------------------------------------------------------

#' Exogenous block for `dates`, which run past the data: each exogenous series
#' carried beyond its last observation by its last monthly step (RES has no
#' 2026 values), then the monthly dummies -- the estimation's `exo` columns,
#' in the same order.
extend_exo <- function(data, dates, variables_cfg) {
  x <- split_by_role(data, variables_cfg)$exogenous
  ext <- apply(x, 2, function(v) {
    last <- max(which(!is.na(v)))
    c(v[seq_len(last)], v[last] + (v[last] - v[last - 1]) * seq_len(length(dates) - last))
  })
  cbind(ext, monthly_dummies(dates))
}

#' Zero-residual forecast of the `h` months after row `origin` of `y`:
#' simulate_var() from the p realised months up to `origin`, with x's rows for
#' the same months (x and y share row indices). h = 1 is the one-step forecast
#' behind every forecast error below.
forecast_path <- function(coefs, y, x, origin, h) {
  p <- dim(coefs$A)[3]
  rows <- (origin - p + 1):(origin + h)
  sim <- simulate_var(coefs$A, coefs$const, init = y[rows[seq_len(p)], , drop = FALSE],
                      resid_mat = matrix(0, h, ncol(y)), B = coefs$B,
                      exo = x[rows, , drop = FALSE])
  sim[p + seq_len(h), , drop = FALSE]
}

#' One-step-ahead forecast errors u_t = y_t - E_{t-1} y_t for each row in
#' `rows`: fixed coefficients, realised lags. One row per month.
one_step_errors <- function(coefs, y, x, rows) {
  t(vapply(rows, function(t) y[t, ] - forecast_path(coefs, y, x, t - 1, 1)[1, ],
           numeric(ncol(y))))
}

#' Ragged edge: variables missing from `y_t` (NA) get their expectation given
#' the month's OBSERVED one-step errors, E[u_m | u_o] = Sigma_mo Sigma_oo^-1 u_o,
#' added to the one-step forecast `yhat_t`. A complete row passes through.
nowcast_row <- function(y_t, yhat_t, sigma_u) {
  m <- is.na(y_t)
  u_o <- (y_t - yhat_t)[!m]
  y_t[m] <- yhat_t[m] + sigma_u[m, !m, drop = FALSE] %*%
    solve(sigma_u[!m, !m, drop = FALSE], u_o)
  y_t
}

## With perfectly correlated errors (u_2 = 2 u_1) the nowcast must hand back
## the missing value exactly.
stopifnot(isTRUE(all.equal(
  nowcast_row(c(1.5, NA), c(1, 3), matrix(c(1, 2, 2, 4), 2)), c(1.5, 4)
)))

#' L_j[t, s] = Theta_{t-s}[v, j] for s <= t, zero above the diagonal: how shock
#' j in forecast month s moves variable v in month t.
irf_toeplitz <- function(theta, v, j, h) {
  l <- stats::toeplitz(theta[v, j, seq_len(h)])
  l[upper.tri(l)] <- 0
  l
}

#' Minimum-norm shocks putting variable v on the gap r with all K shocks free
#' (Waggoner & Zha, 1999): w = R'(RR')^-1 r, R = [L_1 ... L_K]. Returned as an
#' h x K matrix, month by shock. Rotating B0 rotates w but leaves its norm, and
#' so every implied path, unchanged.
min_norm_shocks <- function(theta, v, r) {
  h <- length(r)
  R <- do.call(cbind, lapply(seq_len(dim(theta)[2]), irf_toeplitz,
                             theta = theta, v = v, h = h))
  matrix(crossprod(R, solve(tcrossprod(R), r)), h)
}

#' Shocks of type j alone putting variable v on the gap r, every other shock
#' zero: L_j w_j = r solved forward, month by month. The diagonal
#' Theta_0[v, j] is non-negative under the sign normalisation.
one_shock_shocks <- function(theta, v, j, r) {
  w <- matrix(0, length(r), dim(theta)[2])
  w[, j] <- forwardsolve(irf_toeplitz(theta, v, j, length(r)), r)
  w
}

#' Deviation of every variable from its forecast caused by the h x K shock
#' matrix w, row 1 = first month: the historical decomposition's convolution,
#' summed over shocks.
shock_path <- function(theta, w) {
  apply(hd_contributions(theta[, , seq_len(nrow(w)), drop = FALSE], w), c(1, 2), sum)
}

## Both solvers must put variable 1 exactly on the gap through shock_path(),
## which pins the Toeplitz layout to hd_contributions()' convolution.
stopifnot({
  set.seed(1)
  th <- array(stats::rnorm(16), c(2, 2, 4))
  th[1, , 1] <- abs(th[1, , 1])
  r <- stats::rnorm(4)
  max(abs(shock_path(th, min_norm_shocks(th, 1, r))[, 1] - r)) < 1e-10 &&
    max(abs(shock_path(th, one_shock_shocks(th, 1, 2, r))[, 1] - r)) < 1e-10
})

#' Everything that depends on the rotation, for one impact matrix b0: A's
#' shocks and the path without the attack shocks, the benchmark month's
#' shocks, and the four forward paths (month x variable x reading, logs) with
#' each reading's largest implied shock. Run on every admissible rotation, so
#' the median-target results are simply element mt_index of the list. `s`
#' holds the rotation-free inputs built in conditional_forecast().
cf_rotation <- function(b0, phi, s) {
  theta <- structural_irf(phi, b0)
  eps <- function(u) t(solve(b0, t(u)))

  e_A <- eps(s$u_A)
  off <- e_A
  off[!s$after_attack, ] <- 0
  off[, -s$attack_cols] <- 0

  ## C1 replays THIS rotation's gas-specific shocks of the window.
  w <- list(reference = min_norm_shocks(theta, s$gas_var, s$r),
            gas = one_shock_shocks(theta, s$gas_var, s$gas_shock, s$r),
            demand = one_shock_shocks(theta, s$gas_var, s$demand_shock, s$r),
            c1 = matrix(0, length(s$r), ncol(b0)))
  w$c1[seq_len(nrow(s$u_c1)), s$gas_shock] <- eps(s$u_c1)[, s$gas_shock]

  ref <- s$yhat + shock_path(theta, w$reference)
  paths <- simplify2array(list(
    reference = ref,
    gas = s$yhat + shock_path(theta, w$gas),
    demand = s$yhat + shock_path(theta, w$demand),
    c1 = ref + shock_path(theta, w$c1)
  ))

  list(eps_A = e_A, y_cf = s$y_A - shock_path(theta, off),
       eps_bench = eps(s$u_bench)[1, ], paths = paths,
       max_shock = vapply(w, function(z) max(abs(z)), numeric(1)))
}

#' D1's single-shock readings divide by Theta_0[TTF, j] month after month, and
#' the sign restriction bounds it only at zero: in the admissible rotations
#' where it is nearly zero the implied shocks explode and TTF leaves the
#' target. A reading is kept for a rotation only where it is numerically
#' defined -- TTF on the target to 1e-8 and every path finite in levels;
#' otherwise its path and shock size become NA and drop out of the bands.
drop_undefined <- function(z, v, target, readings = c("gas", "demand")) {
  for (r in readings) {
    ok <- isTRUE(max(abs(z$paths[, v, r] - target)) < 1e-8) &&
      all(is.finite(exp(z$paths[, , r])))
    if (!ok) {
      z$paths[, , r] <- NA
      z$max_shock[[r]] <- NA
    }
  }
  z
}

#' Pointwise percentile band at `conf_level` over the admissible rotations,
#' from a list holding one array (or vector) per rotation. NAs -- readings
#' drop_undefined() removed -- are skipped.
rotation_band <- function(draws, conf_level) {
  a <- simplify2array(draws)
  keep <- seq_len(length(dim(a)) - 1)
  alpha <- 1 - conf_level
  list(lower = apply(a, keep, stats::quantile, probs = alpha / 2, names = FALSE, na.rm = TRUE),
       upper = apply(a, keep, stats::quantile, probs = 1 - alpha / 2, names = FALSE,
                     na.rm = TRUE))
}

#' Comparison-table numbers for one rotation's forward paths, one row per
#' reading: PUN means over Q4 2026 and 2027 and its peak (EUR/MWh), IPI and
#' consumption at the horizon end in % deviation from the reference path
#' (100 x log difference, as in the IRF and HD figures), and the largest
#' implied shock (s.d.).
scenario_stats <- function(paths, max_shock, fdates) {
  pun <- exp(paths[, "pun", ])
  h <- dim(paths)[1]
  dev <- function(v) 100 * (paths[h, v, ] - paths[h, v, "reference"])
  q4 <- quarters(fdates) == "Q4" & format(fdates, "%Y") == "2026"
  cbind(pun_q4_2026 = colMeans(pun[q4, , drop = FALSE]),
        pun_2027 = colMeans(pun[format(fdates, "%Y") == "2027", , drop = FALSE]),
        pun_peak = apply(pun, 2, max),
        ipi_dev = dev("ipi"), cons_dev = dev("energy_consumption"),
        max_shock = max_shock)
}

#' "MT [lower, upper]", or MT alone where every rotation gives the same value.
fmt_band <- function(mt, lo, hi, digits = 2) {
  f <- paste0("%.", digits, "f")
  ifelse(hi - lo < 1e-8, sprintf(f, mt), sprintf(paste0(f, " [", f, ", ", f, "]"), mt, lo, hi))
}

#' Section 5.4 end to end, on the fitted model, its admissible rotations and
#' its historical decomposition (whose median-target shocks are A's
#' benchmark). Stops if the one-step errors do not reproduce the VAR's last
#' residual, if TTF misses the target under any D1 reading that is kept, if
#' the median-target model loses a single-shock reading, or if the reference
#' path depends on the rotation.
conditional_forecast <- function(var_out, ident, hd, config) {
  cf <- config$cf
  fit <- var_out$fit
  coefs <- extract_var_coefs(fit, var_out$p, var_out$exo)
  shock_names <- c(ident$shock_names, "Unlabelled")

  ## Every series to cf$data_end, same transforms; the exogenous block runs on
  ## to cf$horizon_end so the forecast has regressors.
  cfg <- config
  cfg$sample$end <- cf$data_end
  data <- transform_data(load_all_series(cfg), config$variables)
  y <- split_by_role(data, config$variables)$endogenous
  origin <- nrow(y)
  fdates <- seq(data$date[origin], as.Date(cf$horizon_end), by = "month")[-1]
  x <- extend_exo(data, c(data$date, fdates), config$variables)
  stopifnot(identical(colnames(y), ident$var_names),
            identical(colnames(x), colnames(var_out$exo)))

  ## Check 1: the last estimation month's one-step error is its VAR residual.
  est_end <- match(as.Date(config$sample$end), data$date)
  err_resid <- max(abs(one_step_errors(coefs, y, x, est_end) -
                         utils::tail(stats::residuals(fit), 1)))
  stopifnot(err_resid < 1e-8)

  ## A uses the months after the sample with all four variables observed; at
  ## most one later month -- the last -- is ragged and gets the nowcast.
  post <- (est_end + 1):origin
  a_rows <- post[stats::complete.cases(y[post, , drop = FALSE])]
  stopifnot(length(a_rows) > 0, a_rows == post[seq_along(a_rows)],
            length(post) - length(a_rows) <= 1)
  nowcast <- colnames(y)[is.na(y[origin, ])]
  y[origin, ] <- nowcast_row(y[origin, ], forecast_path(coefs, y, x, origin - 1, 1)[1, ],
                             summary(fit)$covres)

  a_dates <- data$date[a_rows]
  attack <- as.Date(cf$attack_month)
  attack_cols <- match(cf$attack_shocks, shock_names)
  stopifnot(attack %in% a_dates, !anyNA(attack_cols))

  u_hist <- stats::residuals(fit)
  c1_rows <- which(hd$dates >= as.Date(cf$c1_window[1]) & hd$dates <= as.Date(cf$c1_window[2]))
  h <- length(fdates)
  stopifnot(length(c1_rows) <= h)

  gas_var <- match("gas_price", colnames(y))
  target <- y[origin, gas_var]
  yhat <- forecast_path(coefs, y, x, origin, h)
  s <- list(u_A = one_step_errors(coefs, y, x, a_rows), y_A = y[a_rows, , drop = FALSE],
            after_attack = a_dates >= attack, attack_cols = attack_cols,
            u_bench = u_hist[hd$dates == as.Date("2022-03-01"), , drop = FALSE],
            u_c1 = u_hist[c1_rows, , drop = FALSE],
            yhat = yhat, r = target - yhat[, gas_var], gas_var = gas_var,
            gas_shock = match("Gas-Specific_Shock", shock_names),
            demand_shock = match("Industrial_Demand_Shock", shock_names))

  phi <- reduced_form_irf(fit, max(h, length(a_rows)) - 1)
  rot <- lapply(ident$draws$b0, cf_rotation, phi = phi, s = s)
  rot <- lapply(rot, drop_undefined, v = gas_var, target = target)
  mt <- rot[[ident$mt$mt_index]]
  band <- function(f) rotation_band(lapply(rot, f), config$bootstrap$conf_level)
  n_defined <- rowSums(vapply(rot, function(z) !is.na(z$max_shock[c("gas", "demand")]),
                              logical(2)))

  ## Check 2: TTF on the target under every D1 reading -- the reference in
  ## every rotation, a single-shock reading wherever it is kept -- and the
  ## median-target model keeps both single-shock readings.
  ## Check 3: the all-shocks path is the same under every admissible rotation.
  d1 <- c("reference", "gas", "demand")
  err_target <- max(vapply(rot, function(z)
    max(abs(z$paths[, gas_var, d1] - target), na.rm = TRUE), numeric(1)))
  err_rotation <- max(vapply(rot, function(z)
    max(abs(z$paths[, , "reference"] - mt$paths[, , "reference"])), numeric(1)))
  stopifnot(length(rot) >= 2, !anyNA(mt$max_shock), err_target < 1e-8, err_rotation < 1e-8)

  ## A: shocks of the attack month against March 2022 and against the sample
  ## distribution of |eps| (median-target), then PUN without the attack shocks.
  eps_band <- band(function(z) z$eps_A)
  bench_band <- band(function(z) z$eps_bench)
  a <- match(attack, a_dates)
  shock_table <- tibble::tibble(
    shock = gsub("_", " ", shock_names),
    attack = fmt_band(mt$eps_A[a, ], eps_band$lower[a, ], eps_band$upper[a, ]),
    bench = fmt_band(mt$eps_bench, bench_band$lower, bench_band$upper),
    pct = vapply(seq_along(shock_names), function(j)
      100 * stats::ecdf(abs(hd$eps[, j]))(abs(mt$eps_A[a, j])), numeric(1))
  )
  names(shock_table) <- c("Shock", format(attack, "%b %Y"), "Mar 2022",
                          "Percentile of |eps| in sample (MT)")

  pun <- match("pun", colnames(y))
  actual <- exp(s$y_A[, pun])
  pun_cf <- c(list(actual = actual, mt = exp(mt$y_cf[, pun])),
              band(function(z) exp(z$y_cf[, pun])))
  gap_band <- band(function(z) actual - exp(z$y_cf[, pun]))
  after <- s$after_attack
  gap_table <- tibble::tibble(
    Month = format(a_dates[after], "%b %Y"),
    `PUN actual` = sprintf("%.1f", actual[after]),
    `PUN without attack shocks` = fmt_band(pun_cf$mt[after], pun_cf$lower[after],
                                           pun_cf$upper[after], 1),
    `Attack contribution` = fmt_band((actual - pun_cf$mt)[after], gap_band$lower[after],
                                     gap_band$upper[after], 1)
  )

  ## Forward block: the comparison table, bands over the per-rotation numbers.
  stats <- lapply(rot, function(z) scenario_stats(z$paths, z$max_shock, fdates))
  stats_band <- rotation_band(stats, config$bootstrap$conf_level)
  st <- stats[[ident$mt$mt_index]]
  col <- function(k, d) fmt_band(st[, k], stats_band$lower[, k], stats_band$upper[, k], d)
  comparison <- tibble::tibble(
    Reading = c("Reference (all shocks)", "D1 gas-specific only",
                "D1 industrial demand only", "C1 escalation"),
    `PUN mean Q4 2026` = col("pun_q4_2026", 1),
    `PUN mean 2027` = col("pun_2027", 1),
    `PUN peak` = col("pun_peak", 1),
    `Peak month (MT)` = format(fdates[apply(mt$paths[, pun, ], 2, which.max)], "%b %Y"),
    `IPI vs reference (%)` = col("ipi_dev", 2),
    `Consumption vs reference (%)` = col("cons_dev", 2),
    `Max |shock| (s.d.)` = col("max_shock", 2)
  )

  list(dates = data$date, y = y, nowcast = nowcast, a_dates = a_dates, attack = attack,
       shock_names = shock_names, eps = c(list(mt = mt$eps_A), eps_band),
       shock_table = shock_table, pun_cf = pun_cf, gap_table = gap_table,
       fdates = fdates, paths = c(list(mt = mt$paths), band(function(z) z$paths)),
       comparison = comparison, n_defined = n_defined, n_rotations = length(rot),
       checks = c(one_step_vs_residual = err_resid, ttf_on_target = err_target,
                  reference_across_rotations = err_rotation))
}
