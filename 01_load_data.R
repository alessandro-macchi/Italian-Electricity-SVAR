## 01_load_data.R -------------------------------------------------------------
## Reads every raw CSV listed in config$variables using ITS OWN delimiter,
## date format and decimal mark, then aligns all series onto one common
## monthly grid so missing months become explicit NA rows instead of
## silently shifting data out of alignment.
## -----------------------------------------------------------------------

#' Locale-independent parser for free-text "Mon-YY" / "Mon YYYY" dates
#' (e.g. "feb-05", "May 2005", "Jul 2026*"). Used when date_format ==
#' "monthname". Anything that doesn't match (section headers, footnotes,
#' blank lines mixed into a raw export) parses to NA and is dropped by the
#' caller -- this is what lets read_one_series() tolerate a messy raw file
#' without hardcoding its junk rows anywhere.
parse_month_name_date <- function(date_strings) {
  month_lookup <- c(jan = 1, feb = 2, mar = 3, apr = 4, may = 5, jun = 6,
                     jul = 7, aug = 8, sep = 9, oct = 10, nov = 11, dec = 12)

  cleaned <- trimws(gsub("\\*+$", "", date_strings))
  matches <- regmatches(cleaned, regexec("^([A-Za-z]{3})[A-Za-z]*[-\\s]+(\\d{2,4})$", cleaned, perl = TRUE))

  out <- as.Date(rep(NA_character_, length(date_strings)))
  for (i in seq_along(matches)) {
    mi <- matches[[i]]
    if (length(mi) != 3) next
    mon <- month_lookup[tolower(mi[2])]
    if (is.na(mon)) next
    yr <- as.integer(mi[3])
    if (yr < 100) yr <- if (yr < 70) 2000L + yr else 1900L + yr
    out[i] <- as.Date(sprintf("%04d-%02d-01", yr, mon))
  }
  out
}

#' Parse a vector of date strings using a per-variable strptime format, or
#' the "monthname" sentinel for free-text month names. Numeric formats
#' without a day component (e.g. "%Y-%m", "%m-%Y") are anchored to the 1st
#' of the month before parsing.
parse_var_date <- function(date_strings, date_format) {
  if (identical(date_format, "monthname")) {
    return(parse_month_name_date(date_strings))
  }

  has_day <- grepl("%d", date_format, fixed = TRUE)
  if (!has_day) {
    date_strings <- paste0(date_strings, "-01")
    date_format <- paste0(date_format, "-%d")
  }
  as.Date(date_strings, format = date_format)
}

#' Read one variable's raw CSV and return a two-column tibble(date, <id>),
#' with dates floored to the first of the month.
read_one_series <- function(var_cfg, raw_dir) {
  path <- file.path(raw_dir, var_cfg$file)
  loc <- readr::locale(decimal_mark = var_cfg$decimal)
  raw <- readr::read_delim(path, delim = var_cfg$delim, locale = loc,
                            show_col_types = FALSE, progress = FALSE)

  month <- lubridate::floor_date(
    parse_var_date(raw[[var_cfg$date_col]], var_cfg$date_format), "month"
  )
  value <- suppressWarnings(as.numeric(raw[[var_cfg$value_col]]))

  series <- tibble::tibble(date = month, value = value)
  series <- dplyr::filter(series, !is.na(date))

  if (anyDuplicated(series$date)) {
    stop("Duplicate months in ", var_cfg$file,
         " after parsing -- check `date_format` in 00_config.R.")
  }

  names(series)[2] <- var_cfg$id
  series
}

#' Build the common monthly grid spanning every series (or the configured
#' start/end, when given).
build_monthly_grid <- function(all_series, sample_cfg) {
  min_date <- min(purrr::map_dbl(all_series, ~ as.numeric(min(.x$date, na.rm = TRUE))))
  max_date <- max(purrr::map_dbl(all_series, ~ as.numeric(max(.x$date, na.rm = TRUE))))

  start <- if (!is.null(sample_cfg$start)) as.Date(sample_cfg$start) else as.Date(min_date, origin = "1970-01-01")
  end   <- if (!is.null(sample_cfg$end))   as.Date(sample_cfg$end)   else as.Date(max_date, origin = "1970-01-01")

  tibble::tibble(date = seq(start, end, by = "month"))
}

#' Load and align every configured variable onto one common monthly grid.
#' Returns a tibble(date, <id_1>, <id_2>, ...).
load_all_series <- function(config) {
  all_series <- purrr::map(config$variables, read_one_series, raw_dir = config$paths$raw_dir)
  grid <- build_monthly_grid(all_series, config$sample)
  data <- purrr::reduce(all_series, dplyr::left_join, by = "date", .init = grid)
  dplyr::arrange(data, date)
}
