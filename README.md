# Thesis SVAR Project

This project studies whether a Structural Vector Autoregression (SVAR) can help
forecast Italian electricity prices, measured by PUN.

The project is still experimental: the final variable set, transformations, lag
length, and SVAR restrictions are not fixed yet. The code is therefore organised
so that the notebook remains readable and the reusable mechanics live in small R
scripts.

## Project Structure

```text
thesis_svar/
├── SVAR_model.Rmd
├── README.md
├── data/
├── plots/
├── R/
│   ├── 00_config.R
│   ├── 01_io.R
│   ├── 02_variable_catalogue.R
│   ├── 03_prepare_data.R
│   ├── 04_transformations.R
│   ├── 05_diagnostics.R
│   ├── 06_var_svar.R
│   ├── 07_forecasting.R
│   └── 08_plotting.R
└── output/
    ├── tables/
    └── models/
```

## Main Workflow

- `SVAR_model.Rmd` is the report layer. It explains the analysis and calls
  reusable functions.
- `R/02_variable_catalogue.R` is the first file to edit when adding variables
  or trying a different model specification.
- `R/00_config.R` stores the current experiment dates, lag choices, and forecast
  target.
- `data/` stores raw CSV files.
- `plots/` stores generated diagnostic and forecast plots.
- `output/` is reserved for generated tables and saved model objects.

## Adding A New Variable

Start in `R/02_variable_catalogue.R` and add one row to `variable_catalogue`.

The key fields are:

- `file`: raw CSV file inside `data/`;
- `separator`: CSV separator, usually `","` or `";"`;
- `value_column`: raw data column to read;
- `number_format`: `"plain"` or `"euro"`;
- `display_name`: readable name used in plots;
- `model_variant`: active transformation for the model;
- `role`: `"endogenous"` or `"exogenous"` for the current experiment;
- `svar_order`: Cholesky order for endogenous variables.

After that, rerun `SVAR_model.Rmd`. The notebook should not need a new
hard-coded CSV reader for the variable.

## Trying A Different Specification

The fast experimentation fields are `model_variant`, `role`, and `svar_order`.

Allowed `model_variant` values are:

- `"level"`: use the raw level series;
- `"log"`: use the logged level series;
- `"d_level"`: use the first difference of the level series;
- `"d_log"`: use the first difference of the logged series.

Allowed `role` values are:

- `"endogenous"`: include in the VAR/SVAR system;
- `"exogenous"`: include as an external regressor;
- `"unused"`: keep loaded and plotted, but exclude from the model.

For example, to experiment with renewable capacity:

```r
# Logged exogenous renewable capacity
model_variant = "log"
role = "exogenous"
svar_order = NA

# First-differenced endogenous renewable capacity
model_variant = "d_log"
role = "endogenous"
svar_order = 2.5

# Drop renewable capacity from the VAR/SVAR model
role = "unused"
svar_order = NA
```

The transformation columns are generated automatically from `variable_id`, so
you do not need to manually create names such as `log_res_capacity` or
`d_log_res_capacity`.

`svar_order` only needs to sort the endogenous variables. Decimal values are
allowed, which is useful when inserting a variable between two existing
variables without renumbering the whole Cholesky order.

## Current Analysis Flow

1. Load raw monthly series from the variable catalogue.
2. Keep the common sample from January 2005 to December 2025.
3. Plot variables in levels and inspect ACF/PACF.
4. Apply preliminary transformations.
5. Run ADF stationarity tests.
6. Difference selected logged variables.
7. Define the current SVAR dataset and Cholesky ordering.
8. Select lags and estimate the reduced-form VAR.
9. Estimate the SVAR with recursive Cholesky restrictions.
10. Run recursive one-step PUN forecast evaluation.

## R Script Map

- `00_config.R`: paths, experiment dates, lag settings, forecast target.
- `01_io.R`: date parsing, number parsing, catalogue-driven CSV loading.
- `02_variable_catalogue.R`: one-row-per-variable metadata and model role.
- `03_prepare_data.R`: common sample filtering, long-to-wide preparation.
- `04_transformations.R`: log transforms, month dummies, first differences.
- `05_diagnostics.R`: ADF tests, required-column checks, VAR diagnostics.
- `06_var_svar.R`: model input construction, VAR estimation, Cholesky SVAR.
- `07_forecasting.R`: recursive forecasts, forecast metrics, Clark-West test.
- `08_plotting.R`: report tables and saved plots.

## Required R Packages

The notebook uses base R plus:

- `knitr`;
- `tseries`;
- `vars`.

Install missing packages with:

```r
install.packages(c("knitr", "tseries", "vars"))
```
