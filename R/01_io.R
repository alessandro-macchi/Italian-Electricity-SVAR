# Input helpers for monthly CSV series.

parse_month <- function(x) {
  as.Date(paste0("01-", trimws(x)), format = "%d-%m-%Y")
}

parse_euro_number <- function(x) {
  x <- gsub("\\.", "", as.character(x))
  x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
}

parse_plain_number <- function(x) {
  suppressWarnings(as.numeric(x))
}

parse_catalogue_number <- function(x, number_format) {
  if (number_format == "euro") {
    return(parse_euro_number(x))
  }

  if (number_format == "plain") {
    return(parse_plain_number(x))
  }

  stop("Unknown number format in variable catalogue: ", number_format)
}

read_catalogue_series <- function(catalogue_row, data_dir = project_paths$data_dir) {
  file_path <- file.path(data_dir, catalogue_row[["file"]])

  raw <- read.csv(
    file_path,
    sep = catalogue_row[["separator"]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  value_column <- catalogue_row[["value_column"]]

  if (!"Date" %in% names(raw)) {
    stop("The file ", file_path, " does not contain a Date column.")
  }

  if (!value_column %in% names(raw)) {
    stop("The file ", file_path, " does not contain column ", value_column, ".")
  }

  data.frame(
    Date = parse_month(raw$Date),
    value = parse_catalogue_number(
      raw[[value_column]],
      catalogue_row[["number_format"]]
    ),
    variable = catalogue_row[["display_name"]],
    variable_id = catalogue_row[["variable_id"]],
    stringsAsFactors = FALSE
  )
}

load_catalogue_data <- function(catalogue, data_dir = project_paths$data_dir) {
  do.call(
    rbind,
    lapply(seq_len(nrow(catalogue)), function(i) {
      read_catalogue_series(catalogue[i, ], data_dir = data_dir)
    })
  )
}

summarise_data_coverage <- function(long_data) {
  coverage <- do.call(
    rbind,
    lapply(split(long_data, long_data$variable), function(x) {
      data.frame(
        variable = x$variable[1],
        start = min(x$Date, na.rm = TRUE),
        end = max(x$Date, na.rm = TRUE),
        observations = nrow(x),
        missing_values = sum(is.na(x$value)),
        stringsAsFactors = FALSE
      )
    })
  )

  row.names(coverage) <- NULL
  coverage
}
