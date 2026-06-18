# VAR and SVAR helpers.

get_monthly_dummy_variables <- function(data) {
  grep("^month_", names(data), value = TRUE)
}

get_exogenous_variables <- function(data,
                                    catalogue = variable_catalogue,
                                    include_month_dummies = TRUE) {
  rows <- catalogue[catalogue$role == "exogenous", ]
  variables <- if (nrow(rows) == 0) {
    character(0)
  } else {
    get_active_model_column(rows)
  }

  if (isTRUE(include_month_dummies)) {
    variables <- c(variables, get_monthly_dummy_variables(data))
  }

  unname(variables)
}

build_svar_input <- function(data,
                             catalogue = variable_catalogue,
                             config = experiment_config) {
  validate_variable_catalogue(catalogue)

  svar_variables <- get_svar_variables(catalogue)
  exogenous_variables <- get_exogenous_variables(
    data,
    catalogue = catalogue,
    include_month_dummies = config$include_month_dummies
  )
  target_required_columns <- forecast_required_columns(
    catalogue,
    target_variable = config$forecast_target_variable
  )

  required_columns <- unique(c(
    "Date",
    target_required_columns,
    svar_variables,
    exogenous_variables
  ))

  check_required_columns(data, required_columns, context = "SVAR dataset")

  data[
    complete.cases(data[, required_columns]),
    required_columns
  ]
}

make_model_variable_table <- function(svar_variables, exogenous_variables) {
  endogenous_table <- data.frame(
    block = rep("endogenous", length(svar_variables)),
    variable = svar_variables,
    stringsAsFactors = FALSE
  )

  exogenous_table <- data.frame(
    block = rep("exogenous", length(exogenous_variables)),
    variable = exogenous_variables,
    stringsAsFactors = FALSE
  )

  rbind(endogenous_table, exogenous_table)
}

make_model_samples <- function(svar_input,
                               svar_variables,
                               exogenous_variables,
                               config = experiment_config) {
  estimation_sample <- svar_input[
    svar_input$Date >= config$estimation_start &
      svar_input$Date <= config$estimation_end,
  ]

  test_sample <- svar_input[
    svar_input$Date >= config$test_start &
      svar_input$Date <= config$test_end,
  ]

  check_non_empty_sample(estimation_sample, "estimation")
  check_non_empty_sample(test_sample, "test")

  endog_estimation <- as.matrix(estimation_sample[, svar_variables])
  exog_estimation <- if (length(exogenous_variables) == 0) {
    NULL
  } else {
    as.matrix(estimation_sample[, exogenous_variables])
  }

  list(
    estimation_sample = estimation_sample,
    test_sample = test_sample,
    endog_estimation = endog_estimation,
    exog_estimation = exog_estimation
  )
}

select_var_lags <- function(endog, exog, lag_max, type = "const") {
  require_package("vars", "lag selection")

  lag_selection <- vars::VARselect(
    y = endog,
    lag.max = lag_max,
    type = type,
    exogen = exog
  )

  criteria_table <- as.data.frame(t(lag_selection$criteria))
  criteria_table$lag <- as.integer(row.names(criteria_table))
  criteria_table <- criteria_table[
    c("lag", row.names(lag_selection$criteria))
  ]

  list(
    raw = lag_selection,
    criteria_table = criteria_table
  )
}

fit_reduced_form_var <- function(endog, exog, selected_lag, type = "const") {
  require_package("vars", "the VAR model")

  vars::VAR(
    y = endog,
    p = selected_lag,
    type = type,
    exogen = exog
  )
}

make_cholesky_amat <- function(svar_variables) {
  k_endog <- length(svar_variables)

  amat <- matrix(NA, nrow = k_endog, ncol = k_endog)
  rownames(amat) <- svar_variables
  colnames(amat) <- svar_variables
  amat[upper.tri(amat)] <- 0

  amat
}

check_cholesky_identification <- function(amat) {
  k_endog <- nrow(amat)
  needed_restrictions <- k_endog * (k_endog - 1) / 2
  zero_restrictions <- sum(amat[upper.tri(amat)] == 0)
  free_lower_off_diagonal <- sum(is.na(amat[lower.tri(amat)]))

  identification_check <- data.frame(
    k_endog = k_endog,
    needed_zero_restrictions = needed_restrictions,
    imposed_zero_restrictions = zero_restrictions,
    free_lower_off_diagonal_parameters = free_lower_off_diagonal,
    exactly_identified = zero_restrictions == needed_restrictions,
    stringsAsFactors = FALSE
  )

  if (!identification_check$exactly_identified) {
    stop("The Cholesky A-matrix is not exactly identified.")
  }

  identification_check
}

fit_svar_cholesky <- function(var_fit, amat) {
  require_package("vars", "the SVAR model")

  vars::SVAR(
    x = var_fit,
    estmethod = "scoring",
    Amat = amat,
    Bmat = NULL
  )
}
