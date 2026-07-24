## 02_transform.R --------------------------------------------------------------
## Reversible per-variable transformations, driven entirely by the
## `transform` field in config$variables. Switching a variable between
## level/log/diff/logdiff is a one-line change in 00_config.R; this file
## never needs to be touched.
## -----------------------------------------------------------------------

#' Apply one of the four supported transformations to a numeric vector.
apply_transform <- function(x, type) {
  switch(
    type,
    level   = x,
    log     = log(x),
    diff    = c(NA_real_, diff(x)),
    logdiff = c(NA_real_, diff(log(x))),
    stop("Unknown transform '", type, "' -- expected level/log/diff/logdiff.")
  )
}

#' Apply each variable's configured transformation to the loaded dataset.
transform_data <- function(data, variables_cfg) {
  for (v in variables_cfg) {
    data[[v$id]] <- apply_transform(data[[v$id]], v$transform)
  }
  data
}
