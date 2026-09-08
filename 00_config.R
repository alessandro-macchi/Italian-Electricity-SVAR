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
    output_dir = "output",
    ## Per-variable exploratory plots (time series, ACF, PACF) are written to
    ## plots_dir/<variable id>/ as PNGs; they are not shown in the report.
    plots_dir  = "plots"
  ),

  ## Common monthly grid. NULL bounds default to the min/max date observed
  ## across all raw series (their union), so months missing at the start or
  ## end of any individual series surface as explicit NA rather than being
  ## silently dropped or misaligned. Fixed here to Feb 2005 - Dec 2025: the
  ## start is not actually binding until gas_storage_monthly.csv (endogenous,
  ## only starts 2011-01), and estimate_var() drops any row with an NA
  ## endogenous value anyway -- Jan 2026 onward (every series runs past Dec
  ## 2025) is dropped deliberately.
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
  ##   offset      optional; constant added to the raw value column at load
  ##               time, before any transform. Use it when a source file stores
  ##               a series shifted off its natural scale. Currently unused.
  ##   det         endogenous only; deterministic terms partialled out before
  ##               the unit-root tests in 03_sanity_checks.R -- any of "const",
  ##               "trend", "season" (centred monthly dummies)
  ##   za          endogenous only; TRUE to also run a Zivot-Andrews test
  ##   header      optional; FALSE if the raw CSV has no header row (column
  ##               names default to X1, X2, ... and date_col/value_col must
  ##               match). Omit for the normal case (TRUE).
  variables = list(
    list(
      id = "pun", file = "pun.csv",
      date_col = "Date", date_format = "%m-%Y",
      value_col = "PUN (€/MWh)", delim = ";", decimal = ",",
      label = "PUN (EUR/MWh)", transform = "log", role = "endogenous",
      det = "const", za = TRUE
    ),
    list(
      id = "gas_price", file = "ttf.csv",
      date_col = "Date", date_format = "%d/%m/%Y",
      value_col = "Price", delim = ";", decimal = ".",
      label = "TTF Gas Price, Settlement (EUR/MWh)", transform = "log", role = "endogenous",
      det = "const", za = TRUE
    ),
    list(
      id = "energy_consumption", file = "energy_consumption.csv",
      date_col = "Date", date_format = "%m-%Y",
      value_col = "Volumi MWh", delim = ";", decimal = ",",
      label = "Electricity Consumption (MWh)", transform = "log", role = "endogenous",
      det = c("const", "season"), za = FALSE
    ),
    list(
      id = "ipi", file = "ita_ipi.csv",
      date_col = "time", date_format = "%Y-%m",
      value_col = "value", delim = ",", decimal = ".",
      label = "Italian Industrial Production Index", transform = "log", role = "endogenous",
      det = c("const", "season"), za = FALSE
    ),
    list(
      id = "res_capacity", file = "ita_res_capacity.csv",
      date_col = "Date", date_format = "%m-%Y",
      value_col = "RES Installed Capacity (MW)", delim = ",", decimal = ".",
      label = "RES Installed Capacity (MW)", transform = "level", role = "exogenous"
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
  ## shock; `restrictions` maps a variable id (must match
  ## `variables[[i]]$id`) to a sign ("+"/"-") and the impulse-response
  ## horizons (0 = impact) over which that sign must hold. The
  ## identification code loops over this list generically -- adding a shock
  ## or a restriction requires no code changes.
  ##
  ## The table below is SATURATED: every one of the 4 variables carries a
  ## sign under every one of the 3 shocks, 12 restrictions in all. That is
  ## what makes the labels stable (see below), and it is also what costs the
  ## acceptance rate -- 0.14%, against 2.5% for the earlier half-empty table.
  ## Restriction horizons are impact-only. Tighter windows (e.g. 0:3) turned
  ## out to be infeasible for this reduced-form VAR (0 rotations admissible
  ## in 20000 draws) -- a real informativeness finding, not a bug.
  ##
  ##                        Industrial   Gas-specific   Non-industrial
  ##   pun                       +             +              +
  ##   gas_price                 +             +              +
  ##   energy_consumption        +             -              +
  ##   ipi                       +             -              -
  ##
  ## `gas_price > 0` under all three shocks is a SIGN NORMALISATION, not an
  ## economic assumption: each column of B0 is defined up to sign, and fixing
  ## it so gas rises makes the three shocks comparable (Guntner et al., 2024,
  ## normalise every shock to raise the real gas price the same way).
  ##
  ## `pun > 0` under all three IS an economic restriction, and it comes from
  ## the merit order: gas-fired plants set the marginal price in the large
  ## majority of Italian hours, so any shock that raises TTF raises PUN, and
  ## a contraction of domestic demand-side pressure raises it directly. It
  ## binds hard -- roughly half the rotations admissible without it have PUN
  ## falling on impact while gas rises, which no account of the Italian merit
  ## order permits. Note what it does and does not settle: the SIGN of the
  ## impact pass-through is imposed, its MAGNITUDE, persistence and share of
  ## forecast error variance are estimated.
  ##
  ## The three restriction sets are PAIRWISE DISJOINT, so no admissible
  ## rotation can satisfy two of them in the same column and the shock labels
  ## cannot switch across draws (the labelling problem of Fry & Pagan, 2011,
  ## sec. 4). `warn_if_labels_overlap()` in 05_sign_restrictions.R checks this
  ## on the draws actually used. Two margins do the separating:
  ##   - `energy_consumption` splits the gas-specific shock (quantity down)
  ##     from both demand shocks (quantity up). Price and quantity moving
  ##     together is a demand shift, in opposite directions a supply shift --
  ##     the 2x2 logic of section 3.1 applied to the electricity market, the
  ##     one market here whose price AND quantity are both observed.
  ##   - `ipi` splits the two demand shocks. Industrial demand raises euro
  ##     area production by construction; non-industrial demand (residential
  ##     and weather-driven heating and cooling) raises energy prices while
  ##     the higher energy bill weighs on industrial output.
  ##
  ## Weather is deliberately NOT controlled for in the exogenous block: the
  ## non-industrial demand shock is meant to carry that content, and putting
  ## HDD/CDD among the exogenous regressors would partial it out of the
  ## residuals and leave the shock with nothing to identify.
  ##
  ## Sign normalisation: a POSITIVE realisation of `Gas-Specific_Shock` is an
  ## ADVERSE disturbance (prices up, activity and quantity down). Read every
  ## IRF and historical-decomposition bar accordingly.
  shocks = list(
    list(
      name = "Industrial_Demand_Shock",
      restrictions = list(
        pun = list(sign = "+", horizons = 0),
        gas_price = list(sign = "+", horizons = 0),
        ipi = list(sign = "+", horizons = 0),
        energy_consumption = list(sign = "+", horizons = 0)
      )
    ),
    list(
      name = "Non-Industrial_Demand_Shock",
      restrictions = list(
        pun = list(sign = "+", horizons = 0),
        gas_price = list(sign = "+", horizons = 0),
        ipi = list(sign = "-", horizons = 0),
        energy_consumption = list(sign = "+", horizons = 0)
      )
    ),
    list(
      name = "Gas-Specific_Shock",
      restrictions = list(
        pun = list(sign = "+", horizons = 0),
        gas_price = list(sign = "+", horizons = 0),
        ipi = list(sign = "-", horizons = 0),
        energy_consumption = list(sign = "-", horizons = 0)
      )
    )
  ),
  ## Sign-restriction search.
  ## n_draws is set from the measured acceptance rate, not by habit: three
  ## shocks and twelve restrictions over a 4x4 rotation leave ~0.16% of draws
  ## admissible, so 10000 draws yielded an admissible SET of 18 -- too few to
  ## call a set. 300000 puts it near 500 (~45 s).
  identification = list(
    n_draws = 300000,
    horizon = 36,
    seed    = 6
  ),

  ## Forecast error variance decomposition. `horizons` are the h reported
  ## in the FEVD table, in months after impact (h = 0 is the impact period).
  ## They are also the only horizons for which per-draw shares are kept in
  ## the bootstrap, so lengthening this list costs pooled memory.
  fevd = list(
    horizons = c(1, 6, 12, 24, 36)
  ),

  ## Frequentist bootstrap (Inoue & Kilian, 2013). The FULL identification
  ## search (draws_per_boot rotations) is re-run inside every replication.
  ## draws_per_boot must clear the same ~0.16% acceptance rate: at 1000 the
  ## expected admissible count per replication was 1.6, 8.6% of replications
  ## found none at all, and half of those that succeeded rested on <= 2
  ## rotations. 5000 puts the expectation near 8. Cost is linear -- the full
  ## run is roughly 40 minutes.
  bootstrap = list(
    n_boot         = 2000,
    draws_per_boot = 5000,
    conf_level     = 0.68,
    seed           = 6
  )
)
