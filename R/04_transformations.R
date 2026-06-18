# Transformation helpers driven by model_variant in the variable catalogue.

make_base_long <- function(long_data) {
  base_long <- long_data
  base_long$variable <- base_long$variable_id
  base_long
}

add_catalogue_variants <- function(data, catalogue = variable_catalogue) {
  validate_variable_catalogue(catalogue)

  for (i in seq_len(nrow(catalogue))) {
    row <- catalogue[i, ]
    level_column <- row$variable_id
    log_column <- variant_column_name(row$variable_id, "log")
    d_level_column <- variant_column_name(row$variable_id, "d_level")
    d_log_column <- variant_column_name(row$variable_id, "d_log")

    if (!level_column %in% names(data)) {
      stop("Missing raw level column for variable_id: ", row$variable_id)
    }

    active_uses_log <- row$role != "unused" &&
      row$model_variant %in% c("log", "d_log")

    if (any(data[[level_column]] <= 0, na.rm = TRUE)) {
      if (active_uses_log) {
        stop(
          "The active model_variant requires strictly positive values for ",
          row$display_name,
          "."
        )
      }

      data[[log_column]] <- NA_real_
    } else {
      data[[log_column]] <- log(data[[level_column]])
    }

    data[[d_level_column]] <- c(NA, diff(data[[level_column]]))
    data[[d_log_column]] <- c(NA, diff(data[[log_column]]))
  }

  data
}

add_month_dummies <- function(data) {
  data$month <- factor(
    format(data$Date, "%m"),
    levels = sprintf("%02d", 1:12),
    labels = month.abb
  )

  month_dummies <- model.matrix(~ month, data = data)[, -1, drop = FALSE]
  colnames(month_dummies) <- paste0("month_", month.abb[-1])

  data <- cbind(data, month_dummies)
  data$month <- NULL

  data
}
