## 05_sign_restrictions.R -------------------------------------------------------
## Identification by sign restrictions (Kilian & Lutkepohl, 2017, Ch. 13):
## Haar-uniform random rotations Q drawn via QR decomposition, structural
## impact matrix B0 = P %*% Q, generic sign checking driven by
## config$shocks, admissible-set search, and Fry & Pagan (2011)
## median-target (MT) selection (avoids reporting the incorrect pointwise
## median IRF).
## -----------------------------------------------------------------------

#' Lower-triangular Cholesky factor P of the VAR's residual covariance,
#' with P %*% t(P) = Sigma_u.
cholesky_P <- function(fit) {
  sigma_u <- summary(fit)$covres
  t(chol(sigma_u))
}

#' Reduced-form MA coefficient matrices Phi_0..Phi_H as a K x K x (H+1) array.
reduced_form_irf <- function(fit, horizon) {
  vars::Phi(fit, nstep = horizon)
}

#' One Haar-uniformly distributed K x K orthogonal matrix, via QR
#' decomposition of a random Gaussian matrix with the standard sign fix
#' (Rubio-Ramirez, Waggoner & Zha, 2010).
random_orthogonal_matrix <- function(k) {
  x <- matrix(stats::rnorm(k * k), nrow = k)
  qr_x <- qr(x)
  q <- qr.Q(qr_x)
  r <- qr.R(qr_x)
  q %*% diag(sign(diag(r)), k, k)
}

#' Structural IRFs Phi_h %*% B0 for h = 0..H, as a K x K x (H+1) array.
structural_irf <- function(phi, b0) {
  irf <- array(NA_real_, dim = dim(phi))
  for (h in seq_len(dim(phi)[3])) irf[, , h] <- phi[, , h] %*% b0
  irf
}

#' Check one shock's restrictions (variable -> sign over horizons) against
#' one column of the structural IRF array.
#'
#' A restriction is keyed EITHER by a variable id -- a sign on that single
#' response -- OR by a free name carrying `weights`, a named vector of
#' variable ids to weights defining a linear contrast of responses. The
#' relative-price restrictions use the latter: c(pun = 1, gas_price = -1) is
#' the electricity-price response net of the response of its marginal fuel
#' cost, which is what separates a gas shock from an electricity-specific
#' supply shock (see the `shocks` comment in 00_config.R).
shock_restrictions_hold <- function(irf_struct, col, restrictions, var_index) {
  for (key in names(restrictions)) {
    r <- restrictions[[key]]
    w <- if (is.null(r$weights)) stats::setNames(1, key) else r$weights
    rows <- var_index[names(w)]

    ## [length(rows), 1, length(horizons)]. The column dimension is a
    ## singleton, so R recycles `w` down the first dimension exactly once per
    ## horizon; summing over that dimension leaves the contrast at each
    ## horizon. The single-variable case is w = 1 and reduces to the plain
    ## response, so both kinds of restriction go through one code path.
    block <- irf_struct[rows, col, r$horizons + 1, drop = FALSE] * w
    values <- apply(block, 3, sum)

    ok <- if (r$sign == "+") all(values >= 0) else all(values <= 0)
    if (!ok) return(FALSE)
  }
  TRUE
}

#' Check every configured shock against the structural IRF array. Shock i
#' (in config$shocks order) is checked against column i of B0 -- valid
#' because Q is Haar-uniform over the full orthogonal group, so the ensemble
#' of random draws already explores every column assignment.
all_shocks_hold <- function(irf_struct, shocks_cfg, var_index) {
  for (s in seq_along(shocks_cfg)) {
    if (!shock_restrictions_hold(irf_struct, col = s,
                                  restrictions = shocks_cfg[[s]]$restrictions,
                                  var_index = var_index)) {
      return(FALSE)
    }
  }
  TRUE
}

#' Warn if two labelled shocks' restriction sets overlap on the admissible
#' draws. Pairwise disjointness is what keeps the shock labels stable from
#' draw to draw: where two sets intersect, the same economic disturbance lands
#' in different columns in different draws and the pooled IRF distribution
#' silently mixes them (Fry & Pagan, 2011, sec. 4). The config's restrictions
#' are built to be disjoint; this catches an edit that breaks that without
#' anyone noticing. Checked on the draws actually used rather than
#' analytically, which is exact for the reported set.
warn_if_labels_overlap <- function(draws, shocks_cfg, var_index) {
  n <- length(shocks_cfg)
  for (i in seq_len(n)) for (j in setdiff(seq_len(n), i)) {
    bad <- sum(vapply(draws$irf, shock_restrictions_hold, logical(1),
                      col = i, restrictions = shocks_cfg[[j]]$restrictions,
                      var_index = var_index))
    if (bad > 0) {
      warning(sprintf(
        paste("%d of %d admissible draws have column %d ('%s') also satisfying",
              "'%s' -- the restriction sets overlap, so the shock labels are",
              "not stable across draws."),
        bad, draws$n_accepted, i, shocks_cfg[[i]]$name, shocks_cfg[[j]]$name),
        call. = FALSE)
    }
  }
  invisible(NULL)
}

#' Draw candidate rotations and keep the admissible ones (those satisfying
#' every shock in config$shocks). Reads config$shocks generically -- adding
#' a shock or a restriction never requires changing this file.
draw_admissible_rotations <- function(fit, shocks_cfg, horizon, n_draws, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  var_names <- colnames(fit$y)
  k <- length(var_names)
  var_index <- stats::setNames(seq_len(k), var_names)

  p_chol <- cholesky_P(fit)
  phi <- reduced_form_irf(fit, horizon)

  accepted_b0 <- vector("list", n_draws)
  accepted_irf <- vector("list", n_draws)
  n_accepted <- 0L

  for (i in seq_len(n_draws)) {
    q <- random_orthogonal_matrix(k)
    b0 <- p_chol %*% q
    irf_struct <- structural_irf(phi, b0)

    if (all_shocks_hold(irf_struct, shocks_cfg, var_index)) {
      n_accepted <- n_accepted + 1L
      accepted_b0[[n_accepted]] <- b0
      accepted_irf[[n_accepted]] <- irf_struct
    }
  }

  list(
    b0 = accepted_b0[seq_len(n_accepted)],
    irf = accepted_irf[seq_len(n_accepted)],
    n_draws = n_draws,
    n_accepted = n_accepted,
    acceptance_rate = n_accepted / n_draws,
    var_names = var_names
  )
}

#' Stack the admissible draws' structural IRFs (only the identified shock
#' columns) into one [draw, var, shock, horizon] array.
stack_irf_array <- function(irf_list, n_shocks) {
  n_acc <- length(irf_list)
  dims <- dim(irf_list[[1]])
  k <- dims[1]
  horizon <- dims[3]

  arr <- array(NA_real_, dim = c(n_acc, k, n_shocks, horizon))
  for (i in seq_len(n_acc)) {
    arr[i, , , ] <- irf_list[[i]][, seq_len(n_shocks), , drop = FALSE]
  }
  arr
}

#' Fry & Pagan (2011) median-target selection: pick the single admissible
#' draw whose IRF is closest (in standardized squared distance) to the
#' pointwise median IRF, instead of reporting the pointwise median itself
#' (which need not correspond to any single admissible structural model).
fry_pagan_median_target <- function(irf_array) {
  med <- apply(irf_array, c(2, 3, 4), stats::median)
  spread <- apply(irf_array, c(2, 3, 4), stats::mad)
  spread[spread == 0] <- 1

  n_acc <- dim(irf_array)[1]
  dist <- vapply(seq_len(n_acc), function(i) {
    sum(((irf_array[i, , , ] - med) / spread)^2)
  }, numeric(1))

  best <- which.min(dist)
  list(median_pointwise = med, mt_index = best, mt_irf = irf_array[best, , , ])
}

#' Full identification pipeline: draw admissible rotations, then pick the
#' Fry-Pagan median-target model. `n_draws`/`seed` default to
#' config$identification but can be overridden (the bootstrap uses this to
#' run a smaller search per replication with an unseeded, independently
#' advancing RNG stream).
identify_sign_restrictions <- function(fit, config,
                                        n_draws = config$identification$n_draws,
                                        seed = config$identification$seed) {
  draws <- draw_admissible_rotations(
    fit, config$shocks, config$identification$horizon, n_draws, seed
  )
  if (draws$n_accepted == 0) {
    stop("No admissible rotations found -- restrictions may be too tight, ",
         "or n_draws too low.")
  }

  ## Only on the reported run: the bootstrap calls this with seed = NULL a
  ## couple of thousand times, and the restriction sets it checks are the
  ## same ones every replication.
  if (!is.null(seed)) {
    warn_if_labels_overlap(draws, config$shocks,
                            stats::setNames(seq_along(draws$var_names),
                                            draws$var_names))
  }

  n_shocks <- length(config$shocks)
  irf_array <- stack_irf_array(draws$irf, n_shocks)
  mt <- fry_pagan_median_target(irf_array)

  list(
    draws = draws,
    irf_array = irf_array,
    mt = mt,
    shock_names = purrr::map_chr(config$shocks, "name"),
    var_names = draws$var_names
  )
}
