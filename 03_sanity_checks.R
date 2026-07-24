## 03_sanity_checks.R -----------------------------------------------------------
## Minimal sanity checks: one compact table (NA count, ADF screening
## p-value, stationarity flag at 5%) and one summary print line. This is a
## quick screening signal only -- the formal unit-root testing for the
## thesis is done separately.
## -----------------------------------------------------------------------

#' Approximate an ADF p-value by interpolating the Dickey-Fuller critical
#' value table returned by urca::ur.df(). Screening use only.
approx_adf_pvalue <- function(stat, cval) {
  probs <- c(0.01, 0.05, 0.10)
  cv <- as.numeric(cval[c("1pct", "5pct", "10pct")])
  if (stat <= min(cv)) return(0.01)
  if (stat >= max(cv)) return(1.00)
  stats::approx(x = cv, y = probs, xout = stat, rule = 2)$y
}

#' Run the ADF screening test on one series and summarize it in one row.
sanity_check_series <- function(x, name) {
  n <- length(x)
  n_na <- sum(is.na(x))
  x_clean <- stats::na.omit(x)

  adf_p <- NA_real_
  if (length(x_clean) >= 10) {
    test <- tryCatch(urca::ur.df(x_clean, type = "drift", selectlags = "AIC"),
                      error = function(e) NULL)
    if (!is.null(test)) adf_p <- approx_adf_pvalue(test@teststat[1], test@cval[1, ])
  }

  tibble::tibble(
    variable        = name,
    n               = n,
    n_na            = n_na,
    adf_p           = adf_p,
    stationary_5pct = !is.na(adf_p) & adf_p < 0.05
  )
}

#' Run sanity checks on every configured variable and print one summary line.
run_sanity_checks <- function(data, variables_cfg) {
  ids <- purrr::map_chr(variables_cfg, "id")
  results <- purrr::map_dfr(ids, ~ sanity_check_series(data[[.x]], .x))

  cat(sprintf(
    "Sanity check: %d variables, %d NA cell(s) total, %d/%d stationary at 5%% (ADF screening).\n",
    nrow(results), sum(results$n_na), sum(results$stationary_5pct), nrow(results)
  ))

  results
}
