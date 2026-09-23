# SVAR of the Italian Electricity Market

## Structure

```
R/00_config.R              Central configuration — the only file you should
                          normally need to edit.
R/01_load_data.R           Reads each raw CSV with its own delimiter/date
                          format/decimal mark, aligns everything onto one
                          common monthly grid (missing months -> explicit NA).
R/02_transform.R           Reversible level/log/diff/logdiff transforms.
R/03_sanity_checks.R       Unit-root (ADF/KPSS) and Johansen cointegration
                          tests.
R/04_var_model.R           Reduced-form VAR: endogenous/exogenous split and
                          lag selection, both driven by config.
R/05_sign_restrictions.R   Identification: Haar-uniform rotations (QR),
                          B0 = P %*% Q, generic sign checking against
                          config$shocks, admissible-set search, Fry-Pagan
                          (2011) median-target selection.
R/06_bootstrap.R           Inoue & Kilian (2013) residual bootstrap: resample
                          residuals, rebuild the sample recursively,
                          re-estimate the VAR, and re-run the FULL
                          identification search inside every replication.
R/07_fevd.R                Forecast error variance decomposition under
                          partial identification: identified shares sum to
                          less than one, remainder left to unidentified
                          shocks.
R/08_historical_decomposition.R  Historical decomposition over the
                          estimation sample using the single Fry-Pagan
                          rotation; no VAR re-estimated, no rotation redrawn.
R/09_conditional_forecast.R  Conditional forecasts (thesis sec. 5.4) with
                          every parameter fixed at the estimation-sample
                          values: ex-post counterfactual, TTF-held-flat
                          scenario, and a replayed gas-specific-shock
                          scenario.
main.Rmd                  Orchestrator: sources the modules above and
                          produces the report, section by section.
data/raw/                 Raw CSVs (one per configured variable).
output/                   Rendered report + saved results (.rds).
```

## Running it

```r
rmarkdown::render("main.Rmd")
```

If R complains it can't find pandoc (common outside RStudio on Windows),
point it at RStudio's bundled copy first:

```r
Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")
rmarkdown::render("main.Rmd")
```

Opening the project in RStudio and clicking "Knit" does this automatically.

## Current data

`data/raw/` holds five series wired into `config$variables`: `pun.csv` (PUN,
the target), `ttf.csv` (TTF gas price, monthly close), `energy_consumption.csv`
(electricity consumption), `ita_ipi.csv` (Italian industrial production
index), and `ita_res_capacity.csv` (RES installed capacity, exogenous). They
arrive with different delimiters, decimal marks, and date-column formats
(`"%m-%Y"` for most, `"%Y-%m"` for `ita_ipi.csv`) — exactly the kind of
mismatch point 2 below is designed to absorb without touching the loader.
The loader also supports a locale-independent `"monthname"` sentinel for
free-text "Mon-YY"/"Mon YYYY" date columns (see `parse_month_name_date()` in
`01_load_data.R`), though none of the currently configured series need it.

A few other raw CSVs (`cdd.csv`, `hdd.csv`, `gas_storage_monthly.csv`,
`covid_hospitalizations.csv`, `covid_crisis_dummy.csv`,
`energy_crisis_dummy.csv`) sit in `data/raw/` but are not yet referenced by
`config$variables`.

## How to modify each of the 8 required points

**1. Scalability — adding a variable.** Append one block to
`config$variables` in `00_config.R` (id, file, date_col, date_format,
value_col, delim, decimal, label, transform, role). Nothing else changes:
the loader, transformer, sanity checks, and VAR split all iterate over this
list generically. (Restricting the new variable's shock response, if
desired, is a separate, deliberate edit to `config$shocks`.)

**2. Date standardization.** Each variable's `date_format` is a strptime
string (e.g. `"%Y-%m-%d"`, `"%d/%m/%Y"`, `"%m-%Y"`), or the sentinel
`"monthname"` for free-text "Mon-YY"/"Mon YYYY" columns like `ttf.csv`.
`load_all_series()` parses every series with its own format, floors dates
to the 1st of the month, and left-joins them all onto one common monthly
grid (`build_monthly_grid()`) — any month a source is missing becomes an
explicit `NA` row instead of shifting the others out of alignment.

**3. Diagnostics.** One compact table per check, in `main.Rmd`'s
"Diagnostics" section: `unit_root_table()` in `03_sanity_checks.R` (ADF +
KPSS, one row per variable with a p-value, in two tables -- level and
after `02_transform.R`'s configured transform), `johansen_test()` (same
file; cointegration rank test on the endogenous variables' level-form
series, against its 5% critical value), and
`residual_diagnostics()`/`stability_table()` in `04_var_model.R`
(Portmanteau/LM and Jarque-Bera tests on the fitted VAR's residuals, and
companion-matrix root moduli). No pass/fail flags anywhere -- read
statistics against p = 0.05 (or the 5% critical value) directly.

**4. Reversible transformations.** Change `transform` in a variable's
config block to `"level"`, `"log"`, `"diff"`, or `"logdiff"` and re-run —
`02_transform.R`'s `apply_transform()` dispatches on that string alone.

**5. Endogenous/exogenous split.** Change `role` to `"endogenous"` or
`"exogenous"` in a variable's config block. `split_by_role()` in
`04_var_model.R` reads it automatically; `estimate_var()` passes the
exogenous block to `vars::VAR(..., exogen = ...)`. `monthly_dummies()`
always appends 11 monthly seasonal dummies (Feb-Dec; Jan is the omitted
reference month) to that block, in place of `vars::VAR()`'s built-in
`season` argument.

**6. Lag selection.** `config$lag_selection$method` is `"fixed"` (uses
`p_fixed`) or `"auto"` (uses `vars::VARselect()` with
`config$lag_selection$criterion`, one of `AIC`/`HQ`/`SC`/`FPE`). Switching
is a one-word change.

**7. Sign restrictions.** `config$shocks` is a list, one entry per shock:
`name` plus a `restrictions` list mapping a variable id to `sign` (`"+"` or
`"-"`) and `horizons` (impulse-response horizons, 0 = impact). Add a shock
by appending another entry — `05_sign_restrictions.R` loops over
`config$shocks` generically and never needs to change. Shock *i* (in list
order) is checked against column *i* of the random rotation `B0`; because
`Q` is Haar-uniform over the whole orthogonal group, this is not a loss of
generality — the ensemble of random draws already explores every column
assignment. The current table is SATURATED — three shocks (industrial
demand, non-industrial demand, gas-specific), each restricting all four
endogenous variables, impact-only — and PAIRWISE DISJOINT, so no admissible
rotation can satisfy two shocks' restrictions in the same column and the
shock labels cannot switch across draws (see the table and discussion in
`00_config.R`). That saturation costs acceptance rate: roughly 0.14% of
draws are admissible, against `config$identification$n_draws = 300000`. If
you tighten restrictions further and get "No admissible rotations found",
that error message means the restrictions are infeasible for this dataset;
loosen the horizons/signs or raise `n_draws`.

**8. Readability.** Each module file has one job and small functions
(cholesky factor, one rotation draw, one restriction check, one bootstrap
replication, etc.) — see the file list above.

## Bootstrap and diagnostics

- `config$identification$n_draws` / `config$bootstrap$draws_per_boot`
  control how many candidate rotations are tried in the main search and in
  each bootstrap replication, respectively.
- The **acceptance rate** (admissible rotations / draws tried) is reported
  for the main identification and, as a distribution, across bootstrap
  replications — a low rate means the restrictions are informative
  (non-refutability diagnostic).
- A bootstrap replication that finds zero admissible rotations in
  `draws_per_boot` draws is dropped and counted in `n_failed`; the
  percentile bands use only the successful replications (`n_ok`).

## FEVD, historical decomposition, and conditional forecast

- `07_fevd.R` — forecast error variance decomposition under **partial**
  identification: the denominator is the full reduced-form forecast error
  variance, so identified shares sum to less than one and the remainder is
  reported explicitly as left to unidentified shocks. Horizons come from
  `config$fevd$horizons`.
- `08_historical_decomposition.R` — how much of each variable's deviation
  from its deterministic baseline each structural shock accounts for,
  period by period, using the single Fry-Pagan median-target rotation;
  nothing is re-estimated or redrawn.
- `09_conditional_forecast.R` (thesis sec. 5.4) — conditional forecasts with
  every parameter fixed at the estimation-sample values: an ex-post
  counterfactual switching off `cf$attack_shocks` from `cf$attack_month`
  on, a TTF-held-flat scenario (D1) read three ways, and a scenario (C1)
  replaying the median-target gas-specific shocks of `cf$c1_window`.
  `cf$attack_month` is currently a placeholder pending confirmation.
