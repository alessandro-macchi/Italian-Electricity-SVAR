# Forecast evaluation helpers.

build_recursive_regressors <- function(data,
                                       forecast_index,
                                       svar_variables,
                                       exogenous_variables,
                                       selected_lag) {
  if (forecast_index <= selected_lag) {
    stop("Forecast date does not have enough lagged observations.")
  }

  new_regressors <- data.frame(matrix(nrow = 1, ncol = 0))

  for (lag in seq_len(selected_lag)) {
    lagged_values <- data[
      forecast_index - lag,
      svar_variables,
      drop = FALSE
    ]
    names(lagged_values) <- paste0(names(lagged_values), ".l", lag)
    new_regressors <- cbind(new_regressors, lagged_values)
  }

  new_regressors$const <- 1

  if (length(exogenous_variables) == 0) {
    return(new_regressors)
  }

  next_exog <- data[forecast_index, exogenous_variables, drop = FALSE]
  cbind(new_regressors, next_exog)
}

forecast_one_step <- function(data,
                              forecast_index,
                              svar_variables,
                              exogenous_variables,
                              selected_lag,
                              catalogue = variable_catalogue,
                              config = experiment_config) {
  require_package("vars", "recursive VAR forecasting")

  forecast_date <- data$Date[forecast_index]
  target <- get_forecast_target_info(
    catalogue,
    target_variable = config$forecast_target_variable
  )

  if (!target$model_column %in% svar_variables) {
    stop("The forecast target must be included as an endogenous variable.")
  }

  recursive_train <- data[
    data$Date >= config$estimation_start &
      data$Date < forecast_date,
  ]

  if (nrow(recursive_train) <= selected_lag) {
    stop("Recursive training sample is too short for ", forecast_date, ".")
  }

  recursive_endog <- as.matrix(recursive_train[, svar_variables])
  recursive_exog <- if (length(exogenous_variables) == 0) {
    NULL
  } else {
    as.matrix(recursive_train[, exogenous_variables])
  }

  recursive_var <- vars::VAR(
    y = recursive_endog,
    p = selected_lag,
    type = "const",
    exogen = recursive_exog
  )

  new_regressors <- build_recursive_regressors(
    data = data,
    forecast_index = forecast_index,
    svar_variables = svar_variables,
    exogenous_variables = exogenous_variables,
    selected_lag = selected_lag
  )

  endog_forecast <- sapply(
    svar_variables,
    function(variable) {
      coefficients <- coef(recursive_var$varresult[[variable]])
      coefficients[is.na(coefficients)] <- 0

      missing_regressors <- setdiff(names(coefficients), names(new_regressors))

      if (length(missing_regressors) > 0) {
        stop(
          "Missing recursive forecast regressors: ",
          paste(missing_regressors, collapse = ", ")
        )
      }

      sum(
        as.numeric(new_regressors[1, names(coefficients)]) *
          coefficients
      )
    }
  )

  previous_index <- forecast_index - 1
  target_forecast <- endog_forecast[[target$model_column]]
  forecast_level <- convert_target_forecast_to_level(
    data = data,
    previous_index = previous_index,
    target_forecast = target_forecast,
    target = target
  )

  data.frame(
    Date = forecast_date,
    forecast_pun = forecast_level,
    actual_pun = data[[target$level_column]][forecast_index],
    rw_pun = data[[target$level_column]][previous_index],
    stringsAsFactors = FALSE
  )
}

convert_target_forecast_to_level <- function(data,
                                             previous_index,
                                             target_forecast,
                                             target) {
  if (target$model_variant == "level") {
    return(target_forecast)
  }

  if (target$model_variant == "log") {
    return(exp(target_forecast))
  }

  if (target$model_variant == "d_level") {
    return(data[[target$level_column]][previous_index] + target_forecast)
  }

  if (target$model_variant == "d_log") {
    return(exp(data[[target$log_column]][previous_index] + target_forecast))
  }

  stop("Cannot convert unknown target model_variant: ", target$model_variant)
}

run_recursive_forecasts <- function(data,
                                    svar_variables,
                                    exogenous_variables,
                                    selected_lag,
                                    catalogue = variable_catalogue,
                                    config = experiment_config) {
  forecast_indices <- which(
    data$Date >= config$test_start &
      data$Date <= config$test_end
  )

  if (length(forecast_indices) == 0) {
    stop("No recursive forecast dates found.")
  }

  do.call(
    rbind,
    lapply(forecast_indices, function(i) {
      forecast_one_step(
        data = data,
        forecast_index = i,
        svar_variables = svar_variables,
        exogenous_variables = exogenous_variables,
        selected_lag = selected_lag,
        catalogue = catalogue,
        config = config
      )
    })
  )
}

forecast_metrics <- function(actual, model_forecast, rw_forecast) {
  model_error <- actual - model_forecast
  rw_error <- actual - rw_forecast

  rmse_model <- sqrt(mean(model_error^2, na.rm = TRUE))
  rmse_rw <- sqrt(mean(rw_error^2, na.rm = TRUE))

  data.frame(
    rmse_svar = rmse_model,
    rmse_rw = rmse_rw,
    mae_svar = mean(abs(model_error), na.rm = TRUE),
    mae_rw = mean(abs(rw_error), na.rm = TRUE),
    theils_u = rmse_model / rmse_rw,
    stringsAsFactors = FALSE
  )
}

clark_west_test <- function(actual, model_forecast, rw_forecast) {
  model_error <- actual - model_forecast
  rw_error <- actual - rw_forecast
  adjustment <- (model_forecast - rw_forecast)^2
  cw_series <- rw_error^2 - (model_error^2 - adjustment)
  cw_series <- cw_series[is.finite(cw_series)]

  if (length(cw_series) < 2 || sd(cw_series) == 0) {
    return(
      data.frame(
        clark_west_statistic = NA_real_,
        clark_west_p_value_one_sided = NA_real_
      )
    )
  }

  statistic <- mean(cw_series) / (sd(cw_series) / sqrt(length(cw_series)))

  data.frame(
    clark_west_statistic = statistic,
    clark_west_p_value_one_sided = 1 - pnorm(statistic),
    stringsAsFactors = FALSE
  )
}

forecast_accuracy_table <- function(forecasts, target_name = "pun") {
  cbind(
    variable = target_name,
    forecast_metrics(
      forecasts$actual_pun,
      forecasts$forecast_pun,
      forecasts$rw_pun
    ),
    clark_west_test(
      forecasts$actual_pun,
      forecasts$forecast_pun,
      forecasts$rw_pun
    )
  )
}
