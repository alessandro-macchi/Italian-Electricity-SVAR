## 03_sanity_checks.R -----------------------------------------------------------
## Pre-estimation diagnostics: ADF and KPSS unit-root tests, one row per
## variable with a p-value, reported in two tables (level, and after
## `02_transform.R`'s configured transform) -- and a Johansen cointegration
## test on the endogenous variables' level-form series, controlling for the
## crisis and seasonal dummies. This file only computes numbers, it never
## branches the pipeline on their values -- no
## pass/fail flags, just statistics against p = 0.05 (or the 5% critical
## value, for Johansen).
## -----------------------------------------------------------------------

#' A variable's "level" series for unit-root/cointegration testing: logged
#' when its configured transform is log-based, raw otherwise -- i.e. the
#' pre-differencing space, so the "level" table checks it before
#' `02_transform.R`'s transform is applied.
level_form <- function(x, transform) {
  if (transform %in% c("log", "logdiff")) return(log(x))
  if (transform == "log1p") return(log1p(x))
  x
}

#' ADF test (H0: unit root): statistic and p-value (urca::punitroot,
#' drift case) for one series. `lags` is ur.df()'s MAXIMUM lag, not a fixed
#' one -- selectlags = "AIC" searches 1:lags, so leaving urca's default
#' lags = 1 collapses the search to a single point and runs every test with
#' one augmentation lag. On monthly data that leaves the seasonal
#' autocorrelation in the test regression and inflates the DF t-ratio, so
#' the maximum is set to one full annual cycle.
adf_stats <- function(x) {
  x <- stats::na.omit(x)
  fit <- tryCatch(urca::ur.df(x, type = "drift", lags = 12, selectlags = "AIC"), error = function(e) NULL)
  if (is.null(fit)) return(c(adf_stat = NA_real_, adf_p = NA_real_))
  stat <- unname(fit@teststat[1])
  c(adf_stat = stat, adf_p = urca::punitroot(stat, N = length(x), trend = "c", statistic = "t"))
}

#' KPSS test (H0: stationary): statistic and p-value for one series,
#' interpolated (stats::approx) from urca::ur.kpss's own critical values --
#' the same table lookup tseries::kpss.test uses internally.
kpss_stats <- function(x) {
  x <- stats::na.omit(x)
  fit <- tryCatch(urca::ur.kpss(x, type = "mu", lags = "short"), error = function(e) NULL)
  if (is.null(fit)) return(c(kpss_stat = NA_real_, kpss_p = NA_real_))
  stat <- unname(fit@teststat[1])
  p <- stats::approx(fit@cval[1, ], c(0.10, 0.05, 0.025, 0.01), xout = stat, rule = 2)$y
  c(kpss_stat = stat, kpss_p = p)
}

#' Unit-root tests, one row per variable: ADF + KPSS statistic and p-value
#' side by side. `space = "level"` tests `level_form()`; `space =
#' "transform"` tests the series as `02_transform.R::apply_transform()`
#' actually feeds it to the model.
unit_root_table <- function(raw_data, variables_cfg, space = c("level", "transform")) {
  space <- match.arg(space)
  purrr::map_dfr(variables_cfg, function(v) {
    series <- if (space == "level") level_form(raw_data[[v$id]], v$transform)
              else apply_transform(raw_data[[v$id]], v$transform)
    adf <- adf_stats(series)
    kpss <- kpss_stats(series)
    tibble::tibble(
      variable  = v$id,
      adf_stat  = adf[["adf_stat"]],   adf_p  = adf[["adf_p"]],
      kpss_stat = kpss[["kpss_stat"]], kpss_p = kpss[["kpss_p"]]
    )
  })
}

#' Johansen trace test (H0: cointegration rank <= r) on the endogenous
#' variables' level-form series, using the VAR's own lag order p. Controls
#' for the same structural breaks the VAR conditions on (energy_crisis,
#' covid_crisis) plus monthly seasonals, via ca.jo()'s `dumvar` -- otherwise
#' those level shifts can register as spurious cointegration.
johansen_test <- function(raw_data, variables_cfg, p) {
  endo_cfg <- Filter(function(v) v$role == "endogenous", variables_cfg)

  levels_mat <- sapply(endo_cfg, function(v) level_form(raw_data[[v$id]], v$transform))
  colnames(levels_mat) <- purrr::map_chr(endo_cfg, "id")

  dumvar <- cbind(
    as.matrix(raw_data[, c("energy_crisis")]),
    monthly_dummies(raw_data$date)
  )

  keep <- stats::complete.cases(levels_mat)
  levels_mat <- levels_mat[keep, , drop = FALSE]
  dumvar <- dumvar[keep, , drop = FALSE]

  jo <- urca::ca.jo(levels_mat, type = "trace", ecdet = "const", K = p,
                     spec = "transitory")

  tibble::tibble(
    rank      = trimws(gsub("\\|", "", rownames(jo@cval))),
    statistic = as.numeric(jo@teststat),
    cv_5pct   = jo@cval[, "5pct"]
  )
}

## Deterministic misspecification diagnostics -----------------------------------
## ADF and KPSS disagreeing is not evidence of a unit root: it is what a
## misspecified deterministic component looks like. Each endogenous variable
## therefore has its own `det` terms in 00_config.R, which are partialled out
## before testing, and the test settings below follow from that choice.

#' Covid-period indicator read straight from data/raw/covid_crisis_dummy.csv
#' and aligned to `dates`. Loaded through read_one_series() rather than
#' config$variables on purpose: here it is a diagnostic regressor, and putting
#' it in config would also put it back into the VAR exogenous block.
covid_dummy <- function(dates, raw_dir) {
  cfg <- list(id = "covid", file = "covid_crisis_dummy.csv",
              date_col = "Date", date_format = "%m-%Y",
              value_col = "covid_crisis_dummy", delim = ",", decimal = ".")
  series <- read_one_series(cfg, raw_dir)
  series$covid[match(dates, series$date)]
}

#' Deterministic design matrix. "const" never appears here -- it is the lm
#' intercept. Monthly dummies are centred (each minus 1/12) so they stay
#' orthogonal to that intercept and the trend coefficient remains pure drift
#' instead of absorbing seasonal means.
det_matrix <- function(dates, terms, covid = NULL) {
  X <- NULL
  if ("trend"  %in% terms) X <- cbind(X, trend = seq_along(dates))
  if ("season" %in% terms) X <- cbind(X, monthly_dummies(dates) - 1 / 12)
  if ("covid"  %in% terms) X <- cbind(X, covid = covid)
  X
}

#' The series with its configured deterministic component partialled out.
det_residuals <- function(y, dates, terms, covid = NULL) {
  X <- det_matrix(dates, terms, covid)
  if (is.null(X)) return(y - mean(y))
  unname(stats::residuals(stats::lm(y ~ X)))
}

#' ADF statistic plus the critical value it must be read against. `type` is
#' "trend" whenever a linear trend was partialled out beforehand: pre-removing
#' a trend shifts the Dickey-Fuller null distribution onto tau_tau (-3.42 at
#' 5%) rather than tau_mu (-2.87). Pre-removing seasonal dummies does not
#' shift it (Dickey, Bell and Miller, 1986), which is why they can be taken
#' out for free.
adf_at <- function(x, type, lags = 12) {
  fit <- urca::ur.df(x, type = type, lags = lags, selectlags = "AIC")
  c(stat = unname(fit@teststat[1]), cv5 = unname(fit@cval[1, "5pct"]),
    lags = fit@lags)
}

#' KPSS statistic on an already deterministic-adjusted series, hence type =
#' "mu". `lags` selects the Newey-West bandwidth: "short" is 4(n/100)^0.25,
#' "long" is 12(n/100)^0.25. Reported at both because a statistic that falls
#' as the bandwidth widens is leftover short-run autocorrelation inflating the
#' long-run variance, not non-stationarity.
kpss_at <- function(x, lags) {
  unname(urca::ur.kpss(x, type = "mu", lags = lags)@teststat[1])
}

#' Zivot-Andrews with an endogenously dated level shift. Two urca details
#' matter for reading the output: model = "intercept" still carries a linear
#' trend in the test regression, and `bpoint` is the last observation BEFORE
#' the shift (its dummy is 0 over 1:bpoint), so the first affected month is
#' bpoint + 1. urca also searches every candidate break in 1:(n-1) with no
#' trimming, so `frac` is returned to show whether the break landed in the
#' unreliable sample tails.
za_at <- function(y, dates, lag = 12) {
  fit <- urca::ur.za(y, model = "intercept", lag = lag)
  b <- fit@bpoint + 1
  list(stat = unname(fit@teststat), cval = unname(fit@cval),
       break_date = dates[b], frac = b / length(y))
}

#' One row per endogenous variable: ADF and KPSS on the deterministic-adjusted
#' level and on its first difference (the I(2) check), plus Zivot-Andrews for
#' the variables config marks za = TRUE. Differencing the adjusted series
#' rather than the raw one keeps both halves of a row on the same footing.
stationarity_table <- function(raw_data, config) {
  endo <- Filter(function(v) v$role == "endogenous", config$variables)
  covid <- covid_dummy(raw_data$date, config$paths$raw_dir)

  purrr::map_dfr(endo, function(v) {
    y <- level_form(raw_data[[v$id]], v$transform)
    ok <- !is.na(y)
    y <- y[ok]
    dates <- raw_data$date[ok]

    u <- det_residuals(y, dates, v$det, covid[ok])
    du <- diff(u)
    adf_lev <- adf_at(u, if ("trend" %in% v$det) "trend" else "drift")
    adf_dif <- adf_at(du, "drift")
    za <- if (isTRUE(v$za)) za_at(y, dates) else NULL

    tibble::tibble(
      variable         = v$id,
      n                = length(y),
      deterministics   = paste(v$det, collapse = " + "),
      adf_level        = adf_lev[["stat"]],
      adf_level_cv5    = adf_lev[["cv5"]],
      kpss_level_short = kpss_at(u, "short"),
      kpss_level_long  = kpss_at(u, "long"),
      adf_diff         = adf_dif[["stat"]],
      kpss_diff_short  = kpss_at(du, "short"),
      kpss_diff_long   = kpss_at(du, "long"),
      za_stat          = if (is.null(za)) NA_real_ else za$stat,
      za_break         = if (is.null(za)) as.Date(NA) else za$break_date
    )
  })
}

#' Row-by-row note for stationarity_table(), so the rows are comparable: what
#' was partialled out of each row, and the critical values every column is
#' read against.
stationarity_note <- function(tbl) {
  rows <- paste0(tbl$variable, " = ", tbl$deterministics, collapse = "; ")
  paste0(
    "Unit-root tests on deterministic-adjusted series. Terms removed before ",
    "testing, by row: ", rows, ". Here const is the regression intercept, ",
    "trend a linear trend, season 11 centred monthly dummies (Feb-Dec, each ",
    "minus 1/12), covid the covid_crisis_dummy series (2020-02 to 2022-05). ",
    "ADF H0: unit root -- levels use type = trend where a trend was removed ",
    "and type = drift otherwise, 5% critical value in adf_level_cv5; first ",
    "differences use type = drift, 5% critical value -2.87. Both search up ",
    "to 12 augmentation lags by AIC. KPSS H0: stationary, type = mu, 5% ",
    "critical value 0.463 (10%: 0.347). Zivot-Andrews H0: unit root with no ",
    "break, model = intercept, 5% critical value -4.80; za_break is the ",
    "first month of the estimated shift."
  )
}

#' Zivot-Andrews detail for the variables config marks za = TRUE: statistic,
#' all three critical values, the estimated break month, and its position in
#' the sample (urca does not trim the candidate breaks, so a fraction near 0
#' or 1 is a warning sign rather than a finding).
za_table <- function(raw_data, config) {
  vars <- Filter(function(v) v$role == "endogenous" && isTRUE(v$za), config$variables)
  purrr::map_dfr(vars, function(v) {
    y <- level_form(raw_data[[v$id]], v$transform)
    ok <- !is.na(y)
    za <- za_at(y[ok], raw_data$date[ok])
    tibble::tibble(
      variable = v$id, statistic = za$stat,
      cv_1pct = za$cval[1], cv_5pct = za$cval[2], cv_10pct = za$cval[3],
      break_date = za$break_date, break_fraction = za$frac
    )
  })
}
