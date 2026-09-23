## 08_historical_decomposition.R ------------------------------------------------
## Historical decomposition over the estimation sample: how much of each
## variable's deviation from its deterministic baseline each structural
## shock accounts for, period by period. Everything comes from objects the
## pipeline has already built -- the fitted reduced form, its residuals, and
## the single Fry-Pagan median-target rotation reported everywhere else. No
## VAR is re-estimated, no rotation redrawn, no second bootstrap run.
##
## Identification is PARTIAL, exactly as in 07_fevd.R. B0 is a full
## invertible K x K matrix, so all K structural shocks are recoverable for
## the chosen rotation, but only the first length(config$shocks) columns
## carry an economic label; the rest are reported as one lumped "other
## shocks" series and never rescaled away.
##
## Depends on reduced_form_irf() and structural_irf() from
## 05_sign_restrictions.R.
## -----------------------------------------------------------------------

#' Contribution of every structural shock to every variable, period by
#' period: contrib[t, i, j] = sum_{h=0}^{t-1} Theta[i, j, h+1] * eps[t-h, j],
#' as a [period, variable, shock] array.
#'
#' That sum is a causal convolution of shock j's path with the (i, j) series
#' of structural MA weights, so it is one stats::filter() call per
#' (variable, shock) pair -- K^2 filters, not a T x K x K triple loop. The
#' n-1 leading zeros truncate the sum at t = 1 (nothing before the
#' estimation sample contributes) and the last n filtered values are the
#' sample periods themselves.
hd_contributions <- function(theta, eps) {
  n <- nrow(eps)
  k <- ncol(eps)

  contrib <- array(NA_real_, dim = c(n, k, k))
  for (j in seq_len(k)) {
    x <- c(rep(0, n - 1), eps[, j])
    for (i in seq_len(k)) {
      contrib[, i, j] <- stats::filter(x, theta[i, j, ], method = "convolution",
                                       sides = 1)[n:(2 * n - 1)]
    }
  }
  contrib
}

#' Historical decomposition of the median-target model, over the T_eff =
#' nrow(endo) - p periods the VAR residuals cover.
#'
#' Structural shocks are the reduced-form residuals run through the same B0
#' the IRFs and the FEVD use (eps_t = B0^-1 u_t); draw_admissible_rotations()
#' appends `b0` and `irf` in lockstep and stack_irf_array() preserves that
#' order, so mt_index indexes the b0 list too. `y_dev` sums ALL K shock
#' contributions and `baseline` is the residual y - y_dev, so the
#' decomposition adds up to the data exactly by construction. The baseline is
#' therefore everything the structural shocks do not move: the intercept, the
#' monthly dummies, the exogenous block, and the p presample values.
historical_decomposition <- function(var_out, ident) {
  fit <- var_out$fit
  b0 <- ident$draws$b0[[ident$mt$mt_index]]

  u <- stats::residuals(fit)
  eps <- t(solve(b0, t(u)))
  theta <- structural_irf(reduced_form_irf(fit, nrow(u) - 1L), b0)

  contrib <- hd_contributions(theta, eps)
  y <- var_out$endo[-seq_len(var_out$p), , drop = FALSE]
  y_dev <- apply(contrib, c(1, 2), sum)

  list(contrib = contrib, theta = theta, eps = eps,
       y = y, y_dev = y_dev, baseline = y - y_dev,
       dates = var_out$dates[-seq_len(var_out$p)],
       var_names = ident$var_names, shock_names = ident$shock_names)
}

## Mechanics on a synthetic K = 2, p = 1 VAR with a known coefficient matrix
## and constant: the decomposition must add up to the data, the impact
## period must hand back the reduced-form residual (Theta_0 eps_1 = B0 B0^-1
## u_1), and a lone shock of size 2 at t = 1 must trace out that shock's IRF
## scaled by 2.
stopifnot({
  set.seed(1)
  a <- matrix(c(0.5, 0.1, -0.2, 0.3), 2, 2)
  e <- matrix(stats::rnorm(200), 100, 2)
  y <- matrix(0, 101, 2, dimnames = list(NULL, c("y1", "y2")))
  for (t in 2:101) y[t, ] <- c(1, -1) + a %*% y[t - 1, ] + e[t - 1, ]

  vo <- list(fit = vars::VAR(y, p = 1, type = "const"), p = 1, endo = y,
             dates = seq(as.Date("2000-01-01"), by = "month", length.out = 101))
  id <- list(draws = list(b0 = list(matrix(c(1, 0.4, -0.2, 0.8), 2, 2))),
             mt = list(mt_index = 1), var_names = colnames(y), shock_names = c("s1", "s2"))

  hd <- historical_decomposition(vo, id)
  one <- hd_contributions(hd$theta, rbind(c(2, 0), matrix(0, 99, 2)))

  max(abs(hd$y_dev + hd$baseline - hd$y)) < 1e-10 &&
    max(abs(hd$y_dev[1, ] - stats::residuals(vo$fit)[1, ])) < 1e-10 &&
    max(abs(one[, , 1] - 2 * t(hd$theta[, 1, ]))) < 1e-10
})
