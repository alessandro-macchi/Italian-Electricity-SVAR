## 04_var_model.R --------------------------------------------------------------
## Reduced-form VAR estimation. The endogenous/exogenous split happens
## automatically from config$variables[[i]]$role, and the lag order is
## either fixed or auto-selected per config$lag_selection.
## -----------------------------------------------------------------------

#' Split the modeling dataset's columns into endogenous/exogenous matrices
#' based on each variable's `role` in config.
split_by_role <- function(data, variables_cfg) {
  ids <- purrr::map_chr(variables_cfg, "id")
  roles <- purrr::map_chr(variables_cfg, "role")

  endo_ids <- ids[roles == "endogenous"]
  exo_ids <- ids[roles == "exogenous"]

  list(
    endogenous = as.matrix(data[, endo_ids, drop = FALSE]),
    exogenous  = if (length(exo_ids) > 0) as.matrix(data[, exo_ids, drop = FALSE]) else NULL
  )
}

#' 11 monthly seasonal dummy columns (Feb..Dec; Jan is the omitted
#' reference month, avoiding collinearity with the VAR's intercept), one
#' row per date -- passed through the exogenous block instead of
#' vars::VAR()'s built-in `season` argument, so they flow through the same
#' `exo` matrix as every other exogenous variable.
monthly_dummies <- function(dates) {
  m <- as.integer(format(dates, "%m"))
  d <- sapply(2:12, function(mm) as.numeric(m == mm))
  colnames(d) <- month.abb[2:12]
  d
}
stopifnot({
  test_dates <- seq(as.Date("2026-01-01"), by = "month", length.out = 12)
  d <- monthly_dummies(test_dates)
  ncol(d) == 11 && all(rowSums(d) == c(0, rep(1, 11)))
})

#' Fixed or automatically selected (vars::VARselect) lag order.
select_lag <- function(endo, lag_cfg) {
  if (lag_cfg$method == "fixed") return(lag_cfg$p_fixed)

  sel <- vars::VARselect(endo, lag.max = lag_cfg$lag_max, type = "const")
  crit_row <- paste0(lag_cfg$criterion, "(n)")
  as.integer(sel$selection[crit_row])
}

#' Estimate the reduced-form VAR, keeping only rows with no NA among the
#' endogenous variables (and the matching exogenous rows). Monthly
#' seasonal dummies are always appended to the exogenous block.
estimate_var <- function(data, config) {
  split <- split_by_role(data, config$variables)

  complete_idx <- stats::complete.cases(split$endogenous)
  endo <- split$endogenous[complete_idx, , drop = FALSE]
  exo <- if (!is.null(split$exogenous)) split$exogenous[complete_idx, , drop = FALSE] else NULL

  dummies <- monthly_dummies(data$date[complete_idx])
  exo <- if (is.null(exo)) dummies else cbind(exo, dummies)

  p <- select_lag(endo, config$lag_selection)
  fit <- vars::VAR(endo, p = p, type = config$var$type, exogen = exo)

  ## `dates` matches `endo` row for row -- the endogenous matrix itself
  ## carries no dates, and the historical decomposition needs them for its
  ## time axis.
  list(fit = fit, p = p, endo = endo, exo = exo, dates = data$date[complete_idx])
}

#' Residual autocorrelation (Portmanteau/LM test) and normality
#' (Jarque-Bera) for the fitted VAR -- one table, native statistic/df/
#' p-value for each, no pass/fail flag.
residual_diagnostics <- function(fit, p) {
  pt <- vars::serial.test(fit, lags.pt = p + 12, type = "PT.asymptotic")$serial
  jb <- vars::normality.test(fit, multivariate.only = TRUE)$jb.mul$JB

  tibble::tibble(
    test      = c("Portmanteau (LM)", "Jarque-Bera"),
    statistic = c(as.numeric(pt$statistic), as.numeric(jb$statistic)),
    df        = c(as.numeric(pt$parameter), as.numeric(jb$parameter)),
    p_value   = c(as.numeric(pt$p.value), as.numeric(jb$p.value))
  )
}

#' Companion-matrix root moduli (all of them -- stability requires every
#' one strictly below 1), sorted descending and wrapped into a compact
#' grid instead of one long column.
stability_table <- function(fit) {
  r <- sort(vars::roots(fit), decreasing = TRUE)
  ncol <- min(6, length(r))
  pad <- (ncol - length(r) %% ncol) %% ncol
  m <- matrix(c(r, rep(NA_real_, pad)), ncol = ncol, byrow = TRUE)
  colnames(m) <- as.character(seq_len(ncol))
  m
}
