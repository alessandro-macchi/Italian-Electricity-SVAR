## 00_config.R ---------------------------------------------------------------
## Central configuration for the SVAR pipeline.
##
## This is the ONLY file that should need editing for day-to-day changes:
##   - add/remove a variable                -> edit `variables`
##   - change a variable's transformation    -> edit `variables[[i]]$transform`
##   - change a variable's role              -> edit `variables[[i]]$role`
##   - change lag selection                  -> edit `lag_selection`
##   - add/edit a structural shock           -> edit `shocks`
##   - change bootstrap settings             -> edit `bootstrap`
## -----------------------------------------------------------------------

config <- list(

  paths = list(
    raw_dir    = "data/raw",
    output_dir = "output"
  ),

  ## Common monthly grid. NULL bounds default to the min/max date observed
  ## across all raw series (their union), so months missing at the start or
  ## end of any individual series surface as explicit NA rather than being
  ## silently dropped or misaligned. Fixed here to Feb 2005 - Dec 2025 (the
  ## thesis's analysis window): ttf.csv only starts Feb 2005, and Jan 2026
  ## onward is dropped deliberately.
  sample = list(
    start = "2005-02-01",
    end   = "2025-12-01"
  ),

  ## One block per variable. To add a variable, append one list() here with
  ## these fields -- no other file needs to change:
  ##   id          internal name; becomes the column name everywhere downstream
  ##   file        CSV file name inside paths$raw_dir
  ##   date_col    name of the date column in the raw CSV
  ##   date_format strptime-style format string for THIS file's date column
  ##               (e.g. "%Y-%m-%d", "%d/%m/%Y", "%Y-%m", "%m-%Y"), OR the
  ##               sentinel "monthname" for locale-independent "Mon-YY" /
  ##               "Mon YYYY" style text dates (see 01_load_data.R)
  ##   value_col   name of the value column in the raw CSV
  ##   delim       column delimiter used in the raw CSV ("," or ";")
  ##   decimal     decimal mark used in the value column ("." or ",")
  ##   label       human-readable label used in tables/plots
  ##   transform   "level" | "log" | "log1p" | "diff" | "logdiff"
  ##   role        "endogenous" | "exogenous"
  ##   header      optional; FALSE if the raw CSV has no header row (column
  ##               names default to X1, X2, ... and date_col/value_col must
  ##               match). Omit for the normal case (TRUE).
  variables = list(
    list(
      id = "pun", file = "pun.csv",
      date_col = "Date", date_format = "%m-%Y",
      value_col = "PUN (€/MWh)", delim = ";", decimal = ",",
      label = "PUN (EUR/MWh)", transform = "logdiff", role = "endogenous"
    ),
    list(
      id = "gas_price", file = "ttf.csv",
      date_col = "Date", date_format = "%d/%m/%Y",
      value_col = "Price", delim = ";", decimal = ".",
      label = "TTF Gas Price, Settlement (EUR/MWh)", transform = "logdiff", role = "endogenous"
    ),
    list(
      id = "energy_consumption", file = "energy_consumption.csv",
      date_col = "Date", date_format = "%m-%Y",
      value_col = "Volumi MWh", delim = ";", decimal = ",",
      label = "Electricity Consumption (MWh)", transform = "logdiff", role = "endogenous"
    ),
    list(
      id = "res_capacity", file = "ita_res_capacity.csv",
      date_col = "Date", date_format = "%m-%Y",
      value_col = "RES Installed Capacity (MW)", delim = ",", decimal = ".",
      label = "RES Installed Capacity (MW)", transform = "logdiff", role = "endogenous"
    ),
    list(
      id = "ipi", file = "eu_ipi.csv",
      date_col = "date", date_format = "%d/%m/%Y",
      value_col = "ipi", delim = ";", decimal = ".",
      label = "EU Industrial Production Index", transform = "level", role = "exogenous"
    ),
    list(
      id = "covid_hospitalizations", file = "covid_hospitalizations.csv",
      date_col = "Date", date_format = "%m-%Y",
      value_col = "covid_hospitalizations", delim = ",", decimal = ".",
      label = "COVID-19 New Hospitalizations (IT, monthly)", transform = "level", role = "exogenous"
    ),
    list(
      id = "energy_crisis", file = "energy_crisis_dummy.csv",
      date_col = "Date", date_format = "%m-%Y",
      value_col = "energy_crisis_dummy", delim = ",", decimal = ".",
      label = "Energy Crisis Dummy (2021-2022)", transform = "level", role = "exogenous"
    )
),

  ## Lag order. method = "fixed" uses p_fixed. method = "auto" uses
  ## vars::VARselect() and picks `criterion` (one of "AIC","HQ","SC","FPE").
  ## AIC with lag_max = 13, per the thesis's model-validation requirements.
  lag_selection = list(
    method    = "fixed",
    p_fixed   = 12,
    criterion = "AIC",
    lag_max   = 13
  ),

  ## Deterministic term passed to vars::VAR().
  var = list(
    type = "const"
  ),

  ## Structural shocks for sign-restriction identification. One entry per
  ## shock; `restrictions` maps a variable id (must match `variables[[i]]$id`)
  ## to a sign ("+"/"-") and the impulse-response horizons (0 = impact) over
  ## which that sign must hold. The identification code loops over this list
  ## generically -- adding a shock or a restriction requires no code changes.
  ## Restriction horizons were picked by checking the admissible-set
  ## acceptance rate for this dataset (see README, point 7): tighter
  ## horizon windows (e.g. 0:3 on both shocks) turned out to be
  ## infeasible for this particular reduced-form VAR (0 rotations
  ## admissible in 20000 draws) -- a real diagnostic finding, not a bug.
  shocks = list(
    list(
      name = "gas_price_shock",
      restrictions = list(
        gas_price = list(sign = "+", horizons = 0)
      )
    ),
    list(
      name = "electricity_demand_shock",
      restrictions = list(
        energy_consumption = list(sign = "+", horizons = 0)
      )
    )
  ),

  ## Sign-restriction search.
  identification = list(
    n_draws = 10000,
    horizon = 24,
    seed    = 6
  ),

  ## Frequentist bootstrap (Inoue & Kilian, 2013). The FULL identification
  ## search (draws_per_boot rotations) is re-run inside every replication.
  bootstrap = list(
    n_boot         = 100,
    draws_per_boot = 500,
    conf_level     = 0.90,
    seed           = 6
  )
)
