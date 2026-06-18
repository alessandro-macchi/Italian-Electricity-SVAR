# Prepare long and wide monthly datasets while keeping modelling choices explicit.

filter_common_sample <- function(long_data, start_date, end_date) {
  out <- long_data[
    long_data$Date >= start_date &
      long_data$Date <= end_date,
  ]

  out$year <- as.numeric(format(out$Date, "%Y"))
  out$month <- as.numeric(format(out$Date, "%m"))
  out$month_name <- factor(
    months(out$Date),
    levels = months(as.Date(paste0("2025-", sprintf("%02d", 1:12), "-01")))
  )

  out
}

summarise_sample_counts <- function(long_data) {
  sample_counts <- as.data.frame(
    table(long_data$variable),
    responseName = "observations"
  )
  names(sample_counts)[1] <- "variable"
  row.names(sample_counts) <- NULL

  sample_counts
}

wide_from_long <- function(long_data) {
  wide <- reshape(
    long_data[, c("Date", "variable", "value")],
    idvar = "Date",
    timevar = "variable",
    direction = "wide"
  )

  names(wide) <- sub("^value\\.", "", names(wide))
  wide <- wide[order(wide$Date), ]
  row.names(wide) <- NULL

  wide
}

prepare_candidate_data <- function(long_data,
                                   catalogue = variable_catalogue,
                                   config = experiment_config) {
  plot_data <- filter_common_sample(
    long_data,
    start_date = config$sample_start,
    end_date = config$sample_end
  )

  base_long <- make_base_long(plot_data)
  svar_data <- wide_from_long(base_long)
  svar_data <- add_catalogue_variants(svar_data, catalogue)

  if (isTRUE(config$include_month_dummies)) {
    svar_data <- add_month_dummies(svar_data)
  }

  list(
    plot_data = plot_data,
    base_long = base_long,
    svar_data = svar_data
  )
}
