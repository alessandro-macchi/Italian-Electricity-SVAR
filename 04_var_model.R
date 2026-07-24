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

#' Fixed or automatically selected (vars::VARselect) lag order.
select_lag <- function(endo, lag_cfg) {
  if (lag_cfg$method == "fixed") return(lag_cfg$p_fixed)

  sel <- vars::VARselect(endo, lag.max = lag_cfg$lag_max, type = "const")
  crit_row <- paste0(lag_cfg$criterion, "(n)")
  as.integer(sel$selection[crit_row])
}

#' Estimate the reduced-form VAR, keeping only rows with no NA among the
#' endogenous variables (and the matching exogenous rows).
estimate_var <- function(data, config) {
  split <- split_by_role(data, config$variables)

  complete_idx <- stats::complete.cases(split$endogenous)
  endo <- split$endogenous[complete_idx, , drop = FALSE]
  exo <- if (!is.null(split$exogenous)) split$exogenous[complete_idx, , drop = FALSE] else NULL

  p <- select_lag(endo, config$lag_selection)
  fit <- vars::VAR(endo, p = p, type = config$var$type, exogen = exo)

  list(fit = fit, p = p, endo = endo, exo = exo)
}
