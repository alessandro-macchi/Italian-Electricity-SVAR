## 03_sanity_checks.R -----------------------------------------------------------
## Pre-estimation diagnostics, one compact table per test: ADF and KPSS
## unit-root tests (level and first difference, every variable) and a
## Johansen cointegration test on the endogenous variables' level-form
## series. This file only computes numbers, it never branches the pipeline
## on their values -- no pass/fail flags, just statistics against their own
## native critical values.
## -----------------------------------------------------------------------

#' A variable's "level" series for unit-root/cointegration testing: logged
#' when its configured transform is log-based, raw otherwise -- i.e. the
#' same space `02_transform.R` would difference, so the test verifies that
#' choice.
level_form <- function(x, transform) {
  if (transform %in% c("log", "logdiff")) return(log(x))
  if (transform == "log1p") return(log1p(x))
  x
}

#' ADF test (H0: unit root) on one series, one row: statistic against its
#' own 10%/5%/1% critical values (urca::ur.df, drift case).
adf_row <- function(x, variable, spec) {
  x <- stats::na.omit(x)
  fit <- tryCatch(urca::ur.df(x, type = "drift", selectlags = "AIC"), error = function(e) NULL)
  if (is.null(fit)) {
    return(tibble::tibble(variable = variable, spec = spec, statistic = NA_real_,
                           cv_10pct = NA_real_, cv_5pct = NA_real_, cv_1pct = NA_real_))
  }
  cv <- fit@cval[1, ]
  tibble::tibble(
    variable  = variable, spec = spec,
    statistic = fit@teststat[1],
    cv_10pct  = unname(cv["10pct"]), cv_5pct = unname(cv["5pct"]), cv_1pct = unname(cv["1pct"])
  )
}

#' KPSS test (H0: stationary) on one series, one row: statistic against its
#' own 10%/5%/2.5%/1% critical values (urca::ur.kpss, level case).
kpss_row <- function(x, variable, spec) {
  x <- stats::na.omit(x)
  fit <- tryCatch(urca::ur.kpss(x, type = "mu", lags = "short"), error = function(e) NULL)
  if (is.null(fit)) {
    return(tibble::tibble(variable = variable, spec = spec, statistic = NA_real_,
                           cv_10pct = NA_real_, cv_5pct = NA_real_, cv_2.5pct = NA_real_, cv_1pct = NA_real_))
  }
  cv <- fit@cval[1, ]
  tibble::tibble(
    variable  = variable, spec = spec,
    statistic = fit@teststat[1],
    cv_10pct  = unname(cv["10pct"]), cv_5pct = unname(cv["5pct"]),
    cv_2.5pct = unname(cv["2.5pct"]), cv_1pct = unname(cv["1pct"])
  )
}

#' ADF, level and first difference, for every configured variable.
run_adf_tests <- function(raw_data, variables_cfg) {
  purrr::map_dfr(variables_cfg, function(v) {
    lvl <- level_form(raw_data[[v$id]], v$transform)
    dplyr::bind_rows(adf_row(lvl, v$id, "level"), adf_row(diff(lvl), v$id, "1st diff"))
  })
}

#' KPSS, level and first difference, for every configured variable.
run_kpss_tests <- function(raw_data, variables_cfg) {
  purrr::map_dfr(variables_cfg, function(v) {
    lvl <- level_form(raw_data[[v$id]], v$transform)
    dplyr::bind_rows(kpss_row(lvl, v$id, "level"), kpss_row(diff(lvl), v$id, "1st diff"))
  })
}

#' Johansen trace test (H0: cointegration rank <= r) on the endogenous
#' variables' level-form series, using the VAR's own lag order p.
johansen_test <- function(raw_data, variables_cfg, p) {
  endo_cfg <- Filter(function(v) v$role == "endogenous", variables_cfg)

  levels_mat <- sapply(endo_cfg, function(v) level_form(raw_data[[v$id]], v$transform))
  colnames(levels_mat) <- purrr::map_chr(endo_cfg, "id")
  levels_mat <- stats::na.omit(levels_mat)

  jo <- urca::ca.jo(levels_mat, type = "trace", ecdet = "const", K = p, spec = "transitory")

  tibble::tibble(
    rank      = trimws(gsub("\\|", "", rownames(jo@cval))),
    statistic = as.numeric(jo@teststat),
    cv_10pct  = jo@cval[, "10pct"],
    cv_5pct   = jo@cval[, "5pct"],
    cv_1pct   = jo@cval[, "1pct"]
  )
}
