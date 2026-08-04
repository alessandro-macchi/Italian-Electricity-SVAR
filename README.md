# SVAR of the Italian Electricity Market — Sign Restrictions

Frequentist SVAR identified by sign restrictions, following Kilian &
Lütkepohl (2017), *Structural Vector Autoregressive Analysis*, Ch. 13, with
frequentist bootstrap inference (Inoue & Kilian, 2013).

## Structure

```
00_config.R              Central configuration — the only file you should
                          normally need to edit.
01_load_data.R            Reads each raw CSV with its own delimiter/date
                          format/decimal mark, aligns everything onto one
                          common monthly grid (missing months -> explicit NA).
02_transform.R            Reversible level/log/diff/logdiff transforms.
03_sanity_checks.R        Unit-root (ADF/KPSS) and Johansen cointegration
                          tests.
04_var_model.R             Reduced-form VAR: endogenous/exogenous split and
                          lag selection, both driven by config.
05_sign_restrictions.R    Identification: Haar-uniform rotations (QR),
                          B0 = P %*% Q, generic sign checking against
                          config$shocks, admissible-set search, Fry-Pagan
                          (2011) median-target selection.
06_bootstrap.R             Inoue & Kilian (2013) residual bootstrap: resample
                          residuals, rebuild the sample recursively,
                          re-estimate the VAR, and re-run the FULL
                          identification search inside every replication.
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

`data/raw/` holds four real series: `pun.csv` (PUN, the target), `ita_res_capacity.csv`
(RES installed capacity), `ttf.csv` (TTF gas price, settlement), and
`energy_consumption.csv` (electricity consumption). They arrived in three
different date layouts and two different decimal conventions — exactly the
kind of mismatch point 2 below is designed to absorb without touching the
loader.

`ttf.csv` in particular is a messy export (looks like a pasted chat
response, complete with a disclaimer and a "copy button" instruction at the
bottom) that mixes two date styles in the same column (`feb-05` vs.
`May 2005`) and has ~20 trailing non-data rows (a "Series Summary" and
narrative notes). `01_load_data.R` handles this via the `date_format =
"monthname"` sentinel (a locale-independent month-name parser — see
`parse_month_name_date()`) and by dropping any row whose date can't be
parsed. **Before using this in the actual thesis, re-source `ttf.csv` from
an actual data provider** (e.g. ICE/Refinitiv/a Eurostat-adjacent series) —
the current file's provenance is not a market-data export and its numbers
should not be trusted for real analysis.

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
exogenous block to `vars::VAR(..., exogen = ...)` when non-empty.

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
assignment. The current restrictions (gas-price shock: impact-only on
`gas_price`; demand shock: impact-only on `energy_consumption`) were chosen
empirically: tighter windows (e.g. 0:3 on both shocks)
turned out to admit **zero** rotations out of 20000 for this dataset — a
real informativeness finding worth discussing, not a bug. If you tighten
restrictions and get "No admissible rotations found", that error message is
telling you exactly this; loosen the horizons or check `config$identification$n_draws`.

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
