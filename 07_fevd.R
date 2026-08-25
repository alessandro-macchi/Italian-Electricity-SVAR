## 07_fevd.R -------------------------------------------------------------------
## Forecast error variance decomposition under PARTIAL identification. Only
## the shocks in config$shocks are identified, so the denominator is the
## reduced-form forecast error variance built from Sigma_u -- not the sum
## over the identified columns of B0. The identified shares therefore sum to
## less than one and the remainder is reported explicitly as the share left
## to the unidentified shocks; rescaling it away would silently claim a full
## identification the sign restrictions do not deliver. vars::fevd() is not
## usable here: it imposes a recursive ordering and ignores the rotation Q.
## -----------------------------------------------------------------------

#' Cumulative reduced-form forecast error variance, sum_{h<=H} diag(Phi_h
#' Sigma_u Phi_h'), as a [var, horizon] matrix. It depends on the reduced
#' form alone, so it is computed once per fitted VAR and recycled across
#' every rotation draw of that fit instead of being rebuilt per draw.
fevd_denominator <- function(phi, sigma_u) {
  per_h <- vapply(seq_len(dim(phi)[3]),
                  function(h) diag(phi[, , h] %*% sigma_u %*% t(phi[, , h])),
                  numeric(dim(phi)[1]))
  t(apply(per_h, 1, cumsum))
}

#' FEVD shares for every draw at once, straight from the already-stacked
#' [draw, var, shock, horizon] structural IRF array: the numerator is just
#' the cumulative sum of squared structural IRFs along the horizon
#' dimension, so no rotation is redrawn and no Phi is recomputed. Returns
#' the same array shape with one extra shock column holding the
#' unidentified remainder 1 - sum_j share_j.
fevd_shares <- function(irf_array, denom) {
  d <- dim(irf_array)
  num <- aperm(apply(irf_array^2, c(1, 2, 3), cumsum), c(2, 3, 4, 1))
  shares <- sweep(num, c(2, 4), denom, "/")

  out <- array(NA_real_, dim = c(d[1], d[2], d[3] + 1L, d[4]))
  out[, , seq_len(d[3]), ] <- shares
  out[, , d[3] + 1L, ] <- 1 - apply(shares, c(1, 2, 4), sum)
  out
}

## With a FULL K-column B0 satisfying B0 B0' = Sigma_u nothing is left
## unidentified: the K shares must sum to exactly 1 at every variable and
## every horizon, all lie in [0, 1], and the remainder column must vanish.
stopifnot({
  a <- matrix(c(0.5, 0.1, -0.2, 0.3), 2, 2)
  phi <- array(c(diag(2), a, a %*% a), dim = c(2, 2, 3))
  b0 <- matrix(c(1, 0.4, -0.2, 0.8), 2, 2)

  irf <- array(NA_real_, dim = c(1, 2, 2, 3))
  for (h in 1:3) irf[1, , , h] <- phi[, , h] %*% b0
  sh <- fevd_shares(irf, fevd_denominator(phi, b0 %*% t(b0)))

  max(abs(apply(sh[1, , 1:2, ], c(1, 3), sum) - 1)) < 1e-12 &&
    all(sh[, , 1:2, ] >= -1e-12 & sh[, , 1:2, ] <= 1 + 1e-12) &&
    max(abs(sh[1, , 3, ])) < 1e-12
})

#' FEVD of the Fry-Pagan median-target model, as a [var, shock + 1,
#' horizon] array over h = 0..horizon. The MT draw rather than a pointwise
#' median of shares, so the table and the reported IRF figure describe the
#' same single admissible structural model.
fevd_median_target <- function(fit, ident, horizon) {
  denom <- fevd_denominator(reduced_form_irf(fit, horizon), summary(fit)$covres)
  fevd_shares(ident$irf_array[ident$mt$mt_index, , , , drop = FALSE], denom)[1, , , ]
}
