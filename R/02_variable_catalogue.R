# One place to describe raw variables and the current experiment.
#
# To try a different specification, usually edit only these columns:
# - model_variant: level, log, d_level, or d_log
# - role: endogenous, exogenous, or unused
# - svar_order: Cholesky order, only for endogenous variables

valid_model_variants <- c("level", "log", "d_level", "d_log")
valid_variable_roles <- c("endogenous", "exogenous", "unused")

variable_catalogue <- data.frame(
  variable_id = c(
    "pun",
    "natural_gas",
    "res_capacity",
    "heating_degree_days",
    "cooling_degree_days",
    "electricity_consumption"
  ),
  file = c(
    "pun.csv",
    "natural_gas.csv",
    "italy_res_capacity_monthly.csv",
    "hdd.csv",
    "cdd.csv",
    "energy_consumption.csv"
  ),
  separator = c(";", ";", ",", ",", ",", ";"),
  value_column = c(
    "PUN (€/MWh)",
    "Natural gas Europe (€/MWh)",
    "RES Installed Capacity (MW)",
    "heating_degree_days",
    "cooling_degree_days",
    "Volumi MWh"
  ),
  number_format = c("euro", "euro", "plain", "plain", "plain", "euro"),
  display_name = c(
    "PUN",
    "Natural gas",
    "RES capacity",
    "Heating degree days",
    "Cooling degree days",
    "Electricity consumption"
  ),
  model_variant = c(
    "d_log",
    "log",
    "d_log",
    "level",
    "level",
    "d_log"
  ),
  role = c(
    "endogenous",
    "endogenous",
    "endogenous",
    "exogenous",
    "exogenous",
    "endogenous"
  ),
  svar_order = c(4, 2, 1, NA, NA, 3),
  stringsAsFactors = FALSE
)

variant_column_name <- function(variable_id, model_variant) {
  if (model_variant == "level") {
    return(variable_id)
  }

  if (model_variant == "log") {
    return(paste0("log_", variable_id))
  }

  if (model_variant == "d_level") {
    return(paste0("d_", variable_id))
  }

  if (model_variant == "d_log") {
    return(paste0("d_log_", variable_id))
  }

  stop("Unknown model_variant: ", model_variant)
}

get_active_model_column <- function(catalogue_rows) {
  if (nrow(catalogue_rows) == 0) {
    return(character(0))
  }

  mapply(
    variant_column_name,
    catalogue_rows$variable_id,
    catalogue_rows$model_variant,
    USE.NAMES = FALSE
  )
}

get_active_model_columns <- function(catalogue = variable_catalogue,
                                     include_unused = FALSE) {
  rows <- catalogue

  if (!include_unused) {
    rows <- rows[rows$role != "unused", ]
  }

  unname(get_active_model_column(rows))
}

get_catalogue_model_columns <- function(catalogue = variable_catalogue) {
  get_active_model_columns(catalogue)
}

get_difference_columns <- function(catalogue = variable_catalogue) {
  rows <- catalogue[
    catalogue$role != "unused" &
      catalogue$model_variant %in% c("d_level", "d_log"),
  ]

  if (nrow(rows) == 0) {
    return(character(0))
  }

  unname(get_active_model_column(rows))
}

get_svar_variables <- function(catalogue = variable_catalogue) {
  rows <- catalogue[catalogue$role == "endogenous", ]

  if (nrow(rows) == 0) {
    stop("At least one endogenous variable is required.")
  }

  if (any(is.na(rows$svar_order))) {
    stop("Every endogenous variable needs a svar_order in the catalogue.")
  }

  if (any(duplicated(rows$svar_order))) {
    stop("Endogenous variables cannot share the same svar_order.")
  }

  rows <- rows[order(rows$svar_order), ]
  unname(get_active_model_column(rows))
}

get_forecast_target_info <- function(catalogue = variable_catalogue,
                                     target_variable = experiment_config$forecast_target_variable) {
  rows <- catalogue[catalogue$variable_id == target_variable, ]

  if (nrow(rows) != 1) {
    stop("The forecast target must match exactly one catalogue row.")
  }

  data.frame(
    variable_id = rows$variable_id,
    display_name = rows$display_name,
    model_variant = rows$model_variant,
    model_column = get_active_model_column(rows),
    level_column = rows$variable_id,
    log_column = variant_column_name(rows$variable_id, "log"),
    stringsAsFactors = FALSE
  )
}

forecast_required_columns <- function(catalogue = variable_catalogue,
                                      target_variable = experiment_config$forecast_target_variable) {
  target <- get_forecast_target_info(catalogue, target_variable)

  required <- c(target$level_column)

  if (target$model_variant == "d_log") {
    required <- c(required, target$log_column)
  }

  unique(required)
}

validate_variable_catalogue <- function(catalogue = variable_catalogue) {
  unknown_variants <- setdiff(catalogue$model_variant, valid_model_variants)
  unknown_roles <- setdiff(catalogue$role, valid_variable_roles)

  if (length(unknown_variants) > 0) {
    stop(
      "Unknown model_variant values: ",
      paste(unknown_variants, collapse = ", ")
    )
  }

  if (length(unknown_roles) > 0) {
    stop("Unknown role values: ", paste(unknown_roles, collapse = ", "))
  }

  invisible(TRUE)
}

describe_variable_catalogue <- function(catalogue = variable_catalogue) {
  validate_variable_catalogue(catalogue)

  data.frame(
    variable_id = catalogue$variable_id,
    file = catalogue$file,
    display_name = catalogue$display_name,
    model_variant = catalogue$model_variant,
    active_model_column = get_active_model_column(catalogue),
    role = catalogue$role,
    svar_order = catalogue$svar_order,
    stringsAsFactors = FALSE
  )
}
