# Project-wide paths and current experiment settings.

project_paths <- list(
  data_dir = "data",
  plots_dir = "plots",
  tables_dir = file.path("output", "tables"),
  models_dir = file.path("output", "models")
)

ensure_project_dirs <- function(paths = project_paths) {
  dir.create(paths$plots_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(paths$tables_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(paths$models_dir, showWarnings = FALSE, recursive = TRUE)

  invisible(paths)
}

experiment_config <- list(
  sample_start = as.Date("2005-01-01"),
  sample_end = as.Date("2025-12-01"),
  estimation_start = as.Date("2005-01-01"),
  estimation_end = as.Date("2021-08-01"),
  test_start = as.Date("2021-09-01"),
  test_end = as.Date("2025-12-01"),
  adf_lag = 12,
  lag_max = 12,
  selected_lag = 12,
  include_month_dummies = TRUE,
  forecast_target_variable = "pun",
  forecast_target_name = "pun"
)

require_package <- function(package, purpose = NULL) {
  if (!requireNamespace(package, quietly = TRUE)) {
    detail <- if (is.null(purpose)) "" else paste0(" for ", purpose)
    stop(
      "Package '", package, "' is needed", detail,
      ". Install it with install.packages('", package, "').",
      call. = FALSE
    )
  }

  invisible(TRUE)
}
