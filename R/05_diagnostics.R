# Lightweight statistical checks used by the notebook.

adf_conclusion <- function(p_value) {
  ifelse(
    p_value < 0.05,
    "Reject unit root: stationary",
    "Do not reject unit root: non-stationary"
  )
}

run_adf_tests <- function(data, variables, lags) {
  require_package("tseries", "the ADF test")

  if (length(variables) == 0) {
    return(
      data.frame(
        variable = character(0),
        lags = integer(0),
        adf_statistic = numeric(0),
        p_value = numeric(0),
        conclusion = character(0)
      )
    )
  }

  do.call(
    rbind,
    lapply(variables, function(v) {
      if (!v %in% names(data)) {
        stop("ADF variable not found: ", v)
      }

      test <- tseries::adf.test(
        na.omit(data[[v]]),
        k = lags
      )

      data.frame(
        variable = v,
        lags = lags,
        adf_statistic = as.numeric(test$statistic),
        p_value = test$p.value,
        conclusion = adf_conclusion(test$p.value),
        stringsAsFactors = FALSE
      )
    })
  )
}

check_required_columns <- function(data, required_columns, context = "dataset") {
  missing_columns <- setdiff(required_columns, names(data))

  if (length(missing_columns) > 0) {
    stop(
      "The ", context, " is missing these columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  invisible(TRUE)
}

check_non_empty_sample <- function(data, sample_name) {
  if (nrow(data) == 0) {
    stop("The ", sample_name, " sample is empty.")
  }

  invisible(TRUE)
}

summarise_sample_sizes <- function(estimation_sample, test_sample) {
  data.frame(
    sample = c("estimation", "test"),
    start = c(min(estimation_sample$Date), min(test_sample$Date)),
    end = c(max(estimation_sample$Date), max(test_sample$Date)),
    observations = c(nrow(estimation_sample), nrow(test_sample)),
    stringsAsFactors = FALSE
  )
}

summarise_var_diagnostics <- function(var_fit, selected_lag) {
  require_package("vars", "VAR diagnostics")

  diagnostic_tests <- list(
    Portmanteau = vars::serial.test(
      var_fit,
      lags.pt = selected_lag + 12,
      type = "PT.asymptotic"
    )$serial,
    "Jarque-Bera" = vars::normality.test(
      var_fit,
      multivariate.only = TRUE
    )$jb.mul$JB,
    ARCH = vars::arch.test(
      var_fit,
      lags.multi = 5,
      multivariate.only = TRUE
    )$arch.mul
  )

  diagnostic_p_values <- sapply(
    diagnostic_tests,
    function(test) as.numeric(test$p.value)
  )

  data.frame(
    test = names(diagnostic_tests),
    statistic = round(
      sapply(diagnostic_tests, function(test) as.numeric(test$statistic)),
      3
    ),
    df = sapply(diagnostic_tests, function(test) {
      parameter <- as.numeric(test$parameter)
      if (length(parameter) == 0) NA_real_ else parameter
    }),
    p_value = signif(diagnostic_p_values, 3),
    decision_5pct = ifelse(
      diagnostic_p_values < 0.05,
      "reject H0",
      "do not reject H0"
    ),
    stringsAsFactors = FALSE
  )
}
