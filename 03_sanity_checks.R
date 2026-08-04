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
#' drift case) for one series.
adf_stats <- function(x) {
  x <- stats::na.omit(x)
  fit <- tryCatch(urca::ur.df(x, type = "drift", selectlags = "AIC"), error = function(e) NULL)
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
