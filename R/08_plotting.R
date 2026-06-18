# Plotting and table helpers used by the report.

make_plot_name <- function(variable, plot_type, plots_dir = project_paths$plots_dir) {
  variable_name <- tolower(variable)
  variable_name <- gsub("[^a-z0-9]+", "_", variable_name)
  variable_name <- gsub("^_|_$", "", variable_name)

  file.path(plots_dir, paste(variable_name, plot_type, sep = "_"))
}

save_plot <- function(file_name, plot_code) {
  png(file_name, width = 10, height = 5, units = "in", res = 300)
  on.exit(dev.off())
  force(plot_code)
}

compact_table <- function(x, row_names = FALSE) {
  knitr::kable(x, row.names = row_names)
}

plot_level <- function(x, variable) {
  plot(
    x$Date,
    x$value,
    type = "l",
    main = paste(variable, "- level"),
    xlab = "Date",
    ylab = variable,
    col = "steelblue",
    lwd = 2
  )
  grid()
}

plot_acf <- function(x, variable) {
  acf(
    x$value,
    lag.max = 36,
    main = paste(variable, "- ACF in levels")
  )
}

plot_pacf <- function(x, variable) {
  pacf(
    x$value,
    lag.max = 36,
    main = paste(variable, "- PACF in levels")
  )
}

plot_level_diagnostics <- function(plot_data, variables = unique(plot_data$variable)) {
  level_plot_functions <- list(
    level_time_series = plot_level,
    acf_levels = plot_acf,
    pacf_levels = plot_pacf
  )

  for (v in variables) {
    x <- plot_data[plot_data$variable == v, ]

    for (plot_type in names(level_plot_functions)) {
      level_plot_functions[[plot_type]](x, v)
      save_plot(
        make_plot_name(v, paste0(plot_type, ".png")),
        level_plot_functions[[plot_type]](x, v)
      )
    }
  }
}

plot_forecast_comparison <- function(forecasts) {
  y_range <- range(
    forecasts$actual_pun,
    forecasts$forecast_pun,
    forecasts$rw_pun,
    na.rm = TRUE
  )

  plot(
    forecasts$Date,
    forecasts$actual_pun,
    type = "l",
    ylim = y_range,
    xlab = "Date",
    ylab = "PUN",
    main = "PUN: VAR Forecast vs Random Walk",
    col = "black",
    lwd = 2
  )
  lines(forecasts$Date, forecasts$forecast_pun, col = "steelblue", lwd = 2)
  lines(forecasts$Date, forecasts$rw_pun, col = "firebrick", lwd = 2, lty = 2)
  grid()
  legend(
    "topleft",
    legend = c("Actual", "VAR/SVAR forecast", "Random walk"),
    col = c("black", "steelblue", "firebrick"),
    lwd = 2,
    lty = c(1, 1, 2),
    bty = "n"
  )
}
