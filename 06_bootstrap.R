## 06_bootstrap.R --------------------------------------------------------------
## Frequentist bootstrap inference (Inoue & Kilian, 2013): resample the
## reduced-form VAR residuals, rebuild the sample recursively, re-estimate
## the VAR, and re-run the FULL sign-restriction search (not just a
## coefficient resample) inside every replication. The exogenous block is
## held FIXED across replications (a fixed-regressor bootstrap): its
## estimated contribution is added back when the sample is rebuilt and the
## same block is passed to every re-estimation, so each replication
## reproduces the specification that was actually estimated rather than a
## seasonless, trendless variant of it. Bands are the pointwise
## percentiles of ALL admissible draws pooled across replications, so they
## span estimation uncertainty (reduced-form coefficients and Sigma_u) AND
## identification uncertainty (the whole set of rotations Q consistent with
## the restrictions).
## -----------------------------------------------------------------------

#' Extract per-equation coefficient matrices A_1..A_p, the intercept, and --
#' when `exo` is supplied -- the k x ncol(exo) matrix B of exogenous
#' coefficients, from a fitted vars::VAR object (deterministic term "const"
#' or "none" only). Exogenous coefficients are named after the exogen
#' matrix's columns, so they are read back by name.
extract_var_coefs <- function(fit, p, exo = NULL) {
  var_names <- names(fit$varresult)
  k <- length(var_names)
  coef_list <- lapply(fit$varresult, stats::coef)

  const <- vapply(coef_list, function(cf) {
    if ("const" %in% names(cf)) unname(cf[["const"]]) else 0
  }, numeric(1))

  A <- array(0, dim = c(k, k, p))
  for (l in seq_len(p)) {
    for (i in seq_len(k)) {
      cf <- coef_list[[i]]
      A[i, , l] <- vapply(var_names, function(vn) unname(cf[[paste0(vn, ".l", l)]]), numeric(1))
    }
  }

  B <- NULL
  if (!is.null(exo)) {
    B <- matrix(0, nrow = k, ncol = ncol(exo),
                dimnames = list(var_names, colnames(exo)))
    for (i in seq_len(k)) {
      cf <- coef_list[[i]]
      B[i, ] <- vapply(colnames(exo), function(xn) unname(cf[[xn]]), numeric(1))
    }
  }

  list(A = A, const = const, B = B, var_names = var_names)
}

#' Recursively simulate a VAR(p) path of the same length as `init` +
#' `nrow(resid_mat)`, using fixed presample values `init` and a given
#' (resampled) residual matrix. `B`/`exo` add the exogenous block's
#' contribution B %*% x_t, with `exo` indexed on the SAME rows as the
#' simulated path (its row t matches sim row t, presample included), so the
#' rebuilt sample carries the seasonality, the RES trend and the crisis
#' dummy that the estimated model attributes to them.
simulate_var <- function(A, const, init, resid_mat, B = NULL, exo = NULL) {
  p <- dim(A)[3]
  k <- dim(A)[1]
  n_new <- nrow(resid_mat)
  Tt <- n_new + p

  sim <- matrix(NA_real_, nrow = Tt, ncol = k)
  sim[seq_len(p), ] <- init

  for (t in (p + 1):Tt) {
    yt <- const
    for (l in seq_len(p)) yt <- yt + A[, , l] %*% sim[t - l, ]
    if (!is.null(B)) yt <- yt + B %*% exo[t, ]
    sim[t, ] <- yt + resid_mat[t - p, ]
  }

  colnames(sim) <- colnames(init)
  sim
}

## The exogenous block must actually reach the simulated path: with a zero
## VAR part, zero intercept and zero residuals, each observation is exactly
## that period's exogenous contribution.
stopifnot({
  s <- simulate_var(
    A = array(0, dim = c(2, 2, 1)), const = c(0, 0),
    init = matrix(0, nrow = 1, ncol = 2), resid_mat = matrix(0, nrow = 2, ncol = 2),
    B = matrix(1, nrow = 2, ncol = 1), exo = matrix(c(0, 5, 7), ncol = 1)
  )
  all(s[2, ] == 5) && all(s[3, ] == 7)
})

#' One bootstrap replication: resample residuals, rebuild the sample,
#' re-estimate the VAR, and re-run the sign-restriction search. Returns the
#' replication's ENTIRE admissible set as a [draw, var, shock, horizon]
#' array -- collapsing it to a single central model here (e.g. the
#' median-target draw) would discard identification uncertainty and leave
#' the bands too narrow. The FEVD is computed here too, from this
#' replication's own IRF array and reduced form, so one bootstrap run yields
#' both IRF and FEVD bands. Returns NULL if no admissible rotation is found.
run_one_bootstrap_replication <- function(fit, endo, exo, p, config) {
  coefs <- extract_var_coefs(fit, p, exo)
  resid_mat <- stats::residuals(fit)
  boot_resid <- resid_mat[sample(seq_len(nrow(resid_mat)), replace = TRUE), , drop = FALSE]

  sim <- simulate_var(coefs$A, coefs$const, init = endo[seq_len(p), , drop = FALSE],
                      resid_mat = boot_resid, B = coefs$B, exo = exo)
  boot_fit <- vars::VAR(sim, p = p, type = config$var$type, exogen = exo)

  ident <- tryCatch(
    identify_sign_restrictions(boot_fit, config,
                                n_draws = config$bootstrap$draws_per_boot, seed = NULL),
    error = function(e) NULL
  )
  if (is.null(ident)) return(NULL)

  ## ponytail: only config$fevd$horizons are kept per draw. Keeping all 37
  ## would add ~1.5x the pooled IRF array's memory for horizons nothing
  ## reports -- the plotted curve is the median-target point estimate, which
  ## keeps every horizon. Widen the subset here if bands are ever needed
  ## across the whole horizon.
  fevd <- fevd_shares(
    ident$irf_array,
    fevd_denominator(reduced_form_irf(boot_fit, config$identification$horizon),
                     summary(boot_fit)$covres)
  )[, , , config$fevd$horizons + 1L, drop = FALSE]

  list(irf = ident$irf_array, fevd = fevd,
       acceptance_rate = ident$draws$acceptance_rate)
}

#' Pool the per-replication [draw, var, shock, horizon] arrays into one
#' array by concatenating along the draw dimension. Replications contribute
#' different numbers of admissible draws, so each is weighted by the size of
#' its own admissible set.
## ponytail: holds every pooled draw in memory (~90 MB at n_boot = 1000 and a
## 2.5% acceptance rate). Accumulate per-cell quantiles incrementally if the
## draw count ever grows past what fits.
stack_replications <- function(irf_arrays) {
  n <- sum(vapply(irf_arrays, function(a) dim(a)[1], integer(1)))
  arr <- array(NA_real_, dim = c(n, dim(irf_arrays[[1]])[-1]))
  j <- 0L
  for (a in irf_arrays) {
    arr[j + seq_len(dim(a)[1]), , , ] <- a
    j <- j + dim(a)[1]
  }
  arr
}

#' Pointwise percentile bootstrap bands at `conf_level` for every
#' (var, shock, horizon) triple, taken over the pooled draw dimension.
bootstrap_bands <- function(draw_array, conf_level) {
  alpha <- 1 - conf_level
  list(
    lower = apply(draw_array, c(2, 3, 4), stats::quantile, probs = alpha / 2, na.rm = TRUE),
    upper = apply(draw_array, c(2, 3, 4), stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  )
}

## Pooling must concatenate along the draw dimension (not overwrite or
## collapse it), and the bands must aggregate over that dimension alone.
stopifnot({
  a <- array(as.numeric(1:24), dim = c(2, 2, 3, 2))
  b <- array(as.numeric(25:48), dim = c(2, 2, 3, 2))
  s <- stack_replications(list(a, b))
  dim(s)[1] == 4L && identical(s[1, , , ], a[1, , , ]) && identical(s[3, , , ], b[1, , , ]) &&
    identical(dim(bootstrap_bands(s, 0.68)$lower), c(2L, 3L, 2L))
})

#' Run the full Inoue-Kilian (2013) bootstrap: `n_boot` replications, each
#' re-running the entire identification search, then pointwise percentile
#' bands over every admissible draw pooled across replications.
run_bootstrap <- function(fit, endo, exo, p, config) {
  set.seed(config$bootstrap$seed)

  replications <- purrr::map(seq_len(config$bootstrap$n_boot), function(b) {
    run_one_bootstrap_replication(fit, endo, exo, p, config)
  })
  ok <- purrr::compact(replications)

  draw_array <- stack_replications(purrr::map(ok, "irf"))
  bands <- bootstrap_bands(draw_array, config$bootstrap$conf_level)

  ## Same shape ([draw, var, shock, horizon]), so pooling and the pointwise
  ## percentiles reuse the IRF machinery unchanged.
  fevd_array <- stack_replications(purrr::map(ok, "fevd"))
  fevd_bands <- bootstrap_bands(fevd_array, config$bootstrap$conf_level)

  list(
    n_boot = config$bootstrap$n_boot,
    n_ok = length(ok),
    n_failed = config$bootstrap$n_boot - length(ok),
    n_pooled_draws = dim(draw_array)[1],
    bands = bands,
    fevd_bands = fevd_bands,
    fevd_horizons = config$fevd$horizons,
    boot_acceptance_rates = purrr::map_dbl(ok, "acceptance_rate")
  )
}
