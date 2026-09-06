#' Coarse-to-fine dynamic (space-time) spatial GLMMs (CF-DGLMMs)
#'
#' Prediction and regression via a separable space-time cascade. Given the
#' scales selected by \code{\link{cf_dglm_hv}}, the model is refitted on the
#' full sample and predictions (with standard deviations) are produced at sample
#' and, optionally, prediction sites. The link-scale linear predictor is
#' \eqn{g(\mu_{i,t}) = x_{i,t}'\beta + \sum_k f_k(s_i,t) + } offset, where each
#' scale-\eqn{k} field \eqn{f_k} couples a per-knot AR(1) Kalman smoother in time
#' with kernel kriging in space.
#'
#' The full-sample fit is a SINGLE coarse-to-fine cascade sweep, mirroring the
#' relationship between \code{cf_glm} and \code{cf_glm_hv}: it reuses the same
#' single greedy sweep that \code{\link{cf_dglm_hv}} performs for scale
#' selection, plus prediction. Within the sweep, for each band (coarse to fine)
#' the GLM working response/weights are refreshed (IRLS folded into the sweep, as
#' \code{cf_glm}'s per-band \code{glm()} does), the scale is fit and accumulated,
#' and the constant and time-varying coefficients are backfit. (The earlier
#' outer-IRLS implementation is archived as \code{cf_dglm_iter} under misc/.)
#'
#' @param y Vector of response variables (N x 1).
#' @param x Matrix of covariates (N x K).
#' @param coords Matrix of 2-dimensional point coordinates (N x 2). The
#'   space-time panel may be unbalanced (observed locations may differ across
#'   time points).
#' @param time Vector of time indices (N x 1); must use the same time points as
#'   in \code{\link{cf_dglm_hv}}.
#' @param offset Optional. Offset variable (N x 1), consistent with \code{glm}.
#' @param x0 Optional. Matrix of covariates at prediction sites (N0 x K).
#' @param coords0 Optional. Coordinates at prediction sites (N0 x 2).
#' @param time0 Optional. Time indices at prediction sites (N0 x 1). May include
#'   time points with no observations: interior time points absent from the
#'   training data are interpolated, and time points beyond the last observed one
#'   are forecast, via the per-knot AR(1) predict step (the Kalman gain is zero
#'   where a time column carries no data). Such time points are added to the
#'   working time grid, so predicting at interior gaps slightly re-spaces the
#'   AR(1) grid; forecasting beyond the last observed point leaves the
#'   training-time fit unchanged.
#' @param offset0 Optional. Offset at prediction sites (N0 x 1).
#' @param mod_hv Output object of \code{\link{cf_dglm_hv}}.
#' @param robust_se Logical; if \code{TRUE} (default), the constant-coefficient
#'   standard errors (and the coefficient-uncertainty term of the predictive SE)
#'   use a spatial-block cluster-robust sandwich that accounts for the cascade
#'   field being a correlated random component. The naive model-based covariance
#'   treats the field as a known offset and severely understates the SEs; the
#'   robust version restores near-nominal coverage. Set \code{FALSE} for the
#'   naive \code{vcov(glm)} SEs.
#' @param sill_cap Logical; if \code{TRUE} (default), the field component of the
#'   predictive variance is capped at the marginal variance of the fitted total
#'   field (the "sill"), \code{var(sum_r z_r)} on the link scale. This prevents
#'   the gPoE variance from diverging in deep extrapolation, mirroring a
#'   stationary GP that reverts to the prior marginal variance far from data. The
#'   cap is applied to the total field variance (never per scale) and leaves the
#'   predictive mean, RMSE, and the coefficient-uncertainty term unchanged. It is
#'   disabled automatically for \code{binomial} responses. Set \code{FALSE} to
#'   leave the field variance uncapped.
#'
#' @return A list (class \code{"cf_dglm"}) mirroring \code{\link{cf_glm}}:
#'   \code{beta}, \code{sd_summary}, \code{e_summary}, \code{pred}, \code{pred0},
#'   \code{pred_q}, \code{pred0_q}, \code{bands}, \code{Z}, \code{Z_sd},
#'   \code{Z0}, \code{Z0_sd}, \code{other}, \code{call}.
#'
#' @seealso \code{\link{cf_dglm_hv}}, \code{\link{cf_glm}}
#' @author Daisuke Murakami
#'
#' @examples
#' ############### Example 1: Space-time disease mapping (Poisson + offset)
#' set.seed(1234)
#' require(CARBayesdata); require(sf)
#' data(pollutionhealthdata); data(GGHB.IZ)
#'
#' ### Space-time panel: 271 areas observed over 2007-2011
#' cent   <- st_coordinates(st_centroid(GGHB.IZ))       # area centroids
#' id     <- match(pollutionhealthdata$IZ, GGHB.IZ$IZ)
#' coords <- cent[id, ]                                  # coordinates per row
#' time   <- pollutionhealthdata$year                   # time index per row
#' y      <- pollutionhealthdata$observed               # disease counts
#' x      <- pollutionhealthdata[, c("pm10", "jsa", "price")]
#' offset <- log(pollutionhealthdata$expected)          # log expected counts
#'
#' ### Holdout validation optimizing the number of spatial scales
#' mod_hv <- cf_dglm_hv(y = y, x = x, coords = coords, time = time,
#'                      offset = offset, family = poisson())
#'
#' ### Space-time modeling and prediction
#' mod    <- cf_dglm(y = y, x = x, coords = coords, time = time,
#'                   offset = offset, mod_hv = mod_hv)
#' mod                                                   # coefficients, rho, Q, tau
#'
#' ### Mapping the fitted relative risk in 2011
#' GGHB.IZ$rr2011 <- mod$pred$pred[time == 2011] /
#'                   pollutionhealthdata$expected[time == 2011]
#' plot(GGHB.IZ[, "rr2011"], lwd = 0.2, axes = TRUE, key.pos = 4, nbreaks = 50)
#'
#' ### Multiscale spatial pattern extraction (averaged over 2007-2011)
#' mod_s1 <- sp_scalewise(mod, bw_range = c(4000, Inf)) # large scale
#' mod_s2 <- sp_scalewise(mod, bw_range = c(0, 4000))   # small scale
#'
#'
#' ############### Example 2: Gaussian panel with sites differing every time point
#' set.seed(1)
#' ns  <- 150; nt <- 10                             # 150 fresh sites at each of 10 times
#' K   <- 25; cen <- cbind(runif(K), runif(K))      # fixed latent centres
#' Wc  <- exp(-as.matrix(dist(cen)) / 0.3); Wc <- Wc / sqrt(rowSums(Wc^2))
#' a   <- matrix(0, nt, K)                          # AR(1) process on the centres
#' a[1, ] <- as.numeric(Wc %*% rnorm(K))
#' for (t in 2:nt) a[t, ] <- 0.7 * a[t - 1, ] + as.numeric(Wc %*% rnorm(K))
#' coords <- cbind(runif(ns * nt), runif(ns * nt))  # new random locations every time
#' time   <- rep(1:nt, each = ns)
#' D   <- sqrt(outer(coords[, 1], cen[, 1], "-")^2 + outer(coords[, 2], cen[, 2], "-")^2)
#' Wp  <- exp(-D / 0.3); Wp <- Wp / rowSums(Wp)
#' field <- rowSums(Wp * a[time, ])                 # continuous field at each (site, time)
#' x1  <- rnorm(ns * nt)
#' y   <- 1 + 2 * x1 + field + rnorm(ns * nt, sd = 0.5)
#'
#' ### The space-time panel is unbalanced (no location is observed twice)
#' mod_hv <- cf_dglm_hv(y = y, x = cbind(x1 = x1), coords = coords, time = time)
#' mod    <- cf_dglm(y = y, x = cbind(x1 = x1), coords = coords, time = time,
#'                   mod_hv = mod_hv)
#' mod
#'
#' ### Coarse vs fine spatial process at time point 5
#' thr     <- stats::median(mod$bands)
#' s_large <- sp_scalewise(mod, bw_range = c(thr, Inf), time_range = c(5, 5))
#' s_small <- sp_scalewise(mod, bw_range = c(0, thr),   time_range = c(5, 5))
#'
#' @importFrom fields rdist
#' @importFrom stats glm gaussian predict vcov qnorm sd as.formula glm.fit lm.wfit
#' @export
cf_dglm <- function(y, x = NULL, coords, time, offset = NULL,
                    x0 = NULL, coords0 = NULL, time0 = NULL, offset0 = NULL,
                    mod_hv, robust_se = TRUE, sill_cap = TRUE) {

  family <- mod_hv$other$family
  bands  <- mod_hv$other$bands
  kernel <- mod_hv$other$kernel
  rho    <- mod_hv$other$rho; Q <- mod_hv$other$Q
  x_sel  <- mod_hv$other$x_sel; xname <- mod_hv$other$xname
  lev    <- mod_hv$other$time_levels
  sk     <- ifelse(is.null(mod_hv$other$seed), 4321, mod_hv$other$seed)
  tau    <- mod_hv$other$tau; if (is.null(tau) || !is.finite(tau) || tau <= 0) tau <- 1
  tv_cols <- mod_hv$other$tv_cols           # design-column indices (in X) with time-varying coefficient
  if (is.null(tv_cols)) tv_cols <- integer(0)
  q_tvc  <- mod_hv$other$q_tvc              # drift variance for the time-varying coefficients
  has_tv <- length(tv_cols) > 0
  n      <- length(y)
  coords <- as.matrix(coords)
  if (is.null(offset)) offset <- rep(0, n)
  has0   <- !is.null(coords0)

  if (has0) {
    if (!is.null(offset) && is.null(offset0) && any(offset != 0))
      stop("offset0 must be provided when offset is specified")
    if (!is.null(x) && is.null(x0)) stop("x0 must be provided when x is specified")
  }

  ## working time grid: extend to include any requested prediction times so that
  ## prediction is possible at time points with no observations. Columns absent
  ## from the training data carry no observation, so each per-knot AR(1) Kalman
  ## uses gain 0 there (predict only): interior gaps are interpolated and times
  ## beyond the last observed point are forecast (variance grows with horizon).
  lev_work <- lev
  if (has0 && !is.null(time0)) lev_work <- sort(unique(c(lev, time0)))

  ## ---- design matrices (intercept + selected covariates)
  if (is.null(x)) X <- matrix(1, n, 1) else { x <- as.matrix(x); X <- cbind(1, x[, x_sel, drop = FALSE]) }
  nx <- ncol(X)
  tv_cols <- tv_cols[tv_cols >= 1 & tv_cols <= nx & tv_cols != 1]  # never the intercept
  has_tv  <- length(tv_cols) > 0
  const_cols <- setdiff(seq_len(nx), tv_cols)                      # intercept + non-tv covariates

  ## ---- panels (training)
  pn  <- .dglm_panel(coords, time, time_levels = lev_work)
  nL  <- pn$nL; T <- pn$T
  Ctr <- pn$C

  ## ---- prediction-site panel (locations x same time grid)
  if (has0) {
    n0 <- nrow(coords0); coords0 <- as.matrix(coords0)
    if (is.null(offset0)) offset0 <- rep(0, n0)
    if (is.null(time0)) stop("time0 must be provided for prediction sites")
    X0  <- if (is.null(x)) matrix(1, n0, 1) else cbind(1, as.matrix(x0)[, x_sel, drop = FALSE])
    pn0 <- .dglm_panel(coords0, time0, time_levels = lev_work)
    Cpr <- pn0$C
  } else { Cpr <- NULL }

  f_obs <- function(field) field[cbind(pn$lk, pn$tk)]
  Xc <- X[, const_cols, drop = FALSE]                     # columns whose coef is constant
  Xtv <- X[, tv_cols, drop = FALSE]                       # columns whose coef is time-varying
  ## Knots and neighbourhoods depend only on the coordinates, so build them ONCE
  ## (train + grid) and reuse; the single sweep then only re-runs the residual-
  ## dependent C++ kernel per band.
  setups <- lapply(seq_along(bands), function(k)
    .dglm_scale_setup(Ctr, bands[k], kernel, sk, Cpr = Cpr))
  tvpart_of <- function(tvbeta) if (has_tv) rowSums(Xtv * tvbeta[pn$tk, , drop = FALSE]) else rep(0, n)

  ## ---- one coarse-to-fine cascade sweep (cf_glm-style, single pass). For each
  ## band, coarse to fine: (i) refresh the GLM working response/weights at the
  ## current linear predictor (IRLS folded into the sweep, as cf_glm's per-band
  ## glm() does); (ii) fit that scale to the working residual and accumulate its
  ## field; (iii) backfit the constant coefficients; (iv) re-estimate the
  ## time-varying coefficients. No outer iteration. `predict` toggles the
  ## prediction-site recombination (TRUE for the single final pass).
  casc_sweep <- function(beta, tvbeta, q_cur, predict = FALSE) {
    tvpart  <- tvpart_of(tvbeta)
    f       <- matrix(0, nL, T)                            # cumulative field (link scale)
    Ftr_sum <- matrix(0, nL, T)
    Fpr_sum <- if (has0 && predict) matrix(0, nrow(Cpr), T) else NULL
    sc_list <- vector("list", length(bands)); z <- w <- NULL
    for (k in seq_along(bands)) {
      fobs <- f_obs(f)
      eta  <- .dglm_clip_l(drop(X %*% beta) + tvpart + fobs + offset, family)
      zw   <- .dglm_work(family, eta, y, offset); z <- zw$z; w <- zw$w
      resid <- z - drop(X %*% beta) - tvpart - fobs  # working residual (z is already offset-free)
      Rp <- matrix(NA_real_, nL, T); Rp[cbind(pn$lk, pn$tk)] <- resid
      Wp <- matrix(NA_real_, nL, T); Wp[cbind(pn$lk, pn$tk)] <- w
      sc <- .dglm_scale_apply(setups[[k]], Rp, Wp, rho, Q, predict = predict)
      f <- f + sc$Ftr; Ftr_sum <- Ftr_sum + sc$Ftr
      if (has0 && predict) Fpr_sum <- Fpr_sum + sc$Fpr
      sc_list[[k]] <- sc
      ## constant-coefficient backfit on the peeled residual
      robs <- (Rp - sc$Ftr)[cbind(pn$lk, pn$tk)]
      ba   <- stats::lm.wfit(Xc, robs, w)$coefficients; ba[!is.finite(ba)] <- 0
      beta[const_cols] <- beta[const_cols] + ba
      ## time-varying-coefficient update (field + constant part removed)
      if (has_tv) {
        r_tv <- z - drop(X %*% beta) - f_obs(f)   # z is already offset-free
        dr   <- .dglm_dynreg(r_tv, Xtv, w, pn$tk, T, q = if (anyNA(q_cur)) NULL else q_cur)
        tvbeta <- dr$beta; if (anyNA(q_cur)) q_cur <- dr$q
        tvpart <- tvpart_of(tvbeta)
      }
    }
    list(beta = beta, tvbeta = tvbeta, q_cur = q_cur,
         Ftr = Ftr_sum, Fpr = Fpr_sum, scales = sc_list, z = z, w = w)
  }

  ## ---- initialize and run the single sweep
  beta <- stats::glm.fit(X, y, offset = offset, family = family)$coefficients
  tvbeta <- if (has_tv) matrix(0, T, length(tv_cols)) else NULL
  if (has_tv) beta[tv_cols] <- 0                          # constant part holds const_cols only
  q_cur <- if (has_tv && !is.null(q_tvc) && all(is.finite(q_tvc)) && all(q_tvc > 0)) q_tvc else NA_real_

  Z <- Z_sd <- matrix(0, n, max(length(bands), 1L))
  Z0 <- Z0_sd <- if (has0) matrix(0, n0, max(length(bands), 1L)) else NULL
  f_tr <- rep(0, n); f0_obs <- if (has0) rep(0, n0) else NULL
  tvpart <- tvpart_of(tvbeta); z <- w <- NULL
  if (length(bands) == 0) {
    eta <- .dglm_clip_l(drop(X %*% beta) + tvpart + offset, family)
    zw  <- .dglm_work(family, eta, y, offset); z <- zw$z; w <- zw$w
    beta[const_cols] <- stats::lm.wfit(Xc, z - tvpart, w)$coefficients
  } else {
    sw <- casc_sweep(beta, tvbeta, q_cur, predict = TRUE)
    beta <- sw$beta; tvbeta <- sw$tvbeta; q_cur <- sw$q_cur; z <- sw$z; w <- sw$w
    for (k in seq_along(bands)) {
      Z[, k]    <- sw$scales[[k]]$Ftr[cbind(pn$lk, pn$tk)]
      Z_sd[, k] <- sqrt(sw$scales[[k]]$Vtr[cbind(pn$lk, pn$tk)])
      if (has0) {
        Z0[, k]    <- sw$scales[[k]]$Fpr[cbind(pn0$lk, pn0$tk)]
        Z0_sd[, k] <- sqrt(sw$scales[[k]]$Vpr[cbind(pn0$lk, pn0$tk)])
      }
    }
    ## center each scale to zero mean (folding the bias into the intercept via
    ## the final GLM), as cf_dglm / cf_glm do. Prediction is unchanged.
    zmean <- colMeans(Z); Z <- sweep(Z, 2, zmean)
    if (has0) Z0 <- sweep(Z0, 2, zmean)   # center by TRAINING means for consistency
    f_tr <- rowSums(Z); if (has0) f0_obs <- rowSums(Z0)
  }

  ## ---- final time-varying coefficients (smoothed value + per-time covariance)
  tvV <- NULL
  if (has_tv) {
    dr <- .dglm_dynreg(z - drop(X %*% beta) - f_tr, Xtv, w, pn$tk, T,
                       q = if (anyNA(q_cur)) NULL else q_cur)
    tvbeta <- dr$beta; tvV <- dr$V; tvpart <- tvpart_of(tvbeta); q_tvc <- dr$q
  }

  ## ---- final GLM with the cascade field (and the time-varying part) as offset
  const_cov <- const_cols[const_cols != 1]                 # constant covariate columns (no intercept)
  ncv <- length(const_cov)
  cov_df <- if (ncv > 0) as.data.frame(X[, const_cov, drop = FALSE]) else NULL
  off_tr <- .dglm_clip_l(f_tr, family) + tvpart + offset
  dat <- data.frame(y = y, .off = off_tr)
  if (!is.null(cov_df)) { names(cov_df) <- xname[const_cov]; dat <- cbind(dat, cov_df) }
  form <- if (ncv > 0) as.formula(paste0("y ~ offset(.off) + ", paste(xname[const_cov], collapse = "+"))) else as.formula("y ~ offset(.off)")
  gmod <- stats::glm(form, data = dat, family = family)
  beta_int <- matrix(gmod$coefficients)
  ## design ordered as gmod (intercept, constant covariates) for variance propagation
  Xg <- if (ncv > 0) cbind(1, X[, const_cov, drop = FALSE]) else matrix(1, n, 1)
  Vbeta <- vcov(gmod)
  ## spatial-block cluster-robust covariance (default): the model-based vcov treats
  ## the cascade field as a known offset and badly understates Var(beta) because the
  ## residual is a correlated random field; .dglm_clusterSE puts the field back into
  ## the working residual and clusters over spatial blocks. Used for both the
  ## coefficient SEs and the coefficient-uncertainty term of the predictive SE.
  V_rob <- Vbeta; G_block <- NA_integer_
  if (robust_se && length(bands) > 0) {
    cse <- tryCatch(.dglm_clusterSE(y, Xg, beta_int, f_tr, tvpart, offset, family,
                                    coords, mod_hv$other$bands),
                    error = function(e) NULL)
    if (!is.null(cse)) { V_rob <- cse$V; G_block <- cse$G }
  }
  Vbeta <- V_rob
  beta_int_se <- sqrt(diag(Vbeta))
  beta_summ <- data.frame(coef = beta_int, coef_se = beta_int_se,
                          lower_95CI = beta_int - 1.96 * beta_int_se,
                          upper_95CI = beta_int + 1.96 * beta_int_se)
  row.names(beta_summ) <- c("Intercept", xname[const_cov])

  ## time-varying-coefficient variance contribution x_tv' V_t x_tv at each obs
  tvvar <- function(Xt, tk_) {
    if (!has_tv) return(rep(0, nrow(Xt)))
    v <- vapply(seq_len(nrow(Xt)), function(i) drop(Xt[i, ] %*% tvV[[tk_[i]]] %*% Xt[i, ]), numeric(1))
    pmax(v, 0)
  }

  ## sill cap: the field predictive variance cannot exceed the marginal variance
  ## of the fitted total field (link scale). Without it the gPoE variance
  ## V = 1/sum(phi/P) diverges where neighbours vanish (deep extrapolation). The
  ## cap is on the TOTAL field variance (never per scale -- scales are positively
  ## correlated, so sum_k var(Z_k) << var(sum_k Z_k) and per-scale caps undershoot).
  ## Disabled for binomial (bounded response, weak logit field) and when no scale.
  sill <- if (isTRUE(sill_cap) && family$family != "binomial" && length(bands) > 0) {
    sv <- stats::var(rowSums(Z)); if (is.finite(sv) && sv > 0) sv else Inf
  } else Inf
  pred     <- predict(gmod, type = "response")
  pred_lin <- predict(gmod, type = "link")
  pred_lin_sd <- sqrt(pmax(rowSums((Xg %*% Vbeta) * Xg) + tvvar(Xtv, pn$tk) +
                           pmin(tau * rowSums(Z_sd^2), sill), 0))
  pred_sd  <- abs(family$mu.eta(pred_lin)) * pred_lin_sd
  qs <- c(0.005, 0.025, 0.05, seq(0.1, 0.9, 0.1), 0.95, 0.975, 0.995)
  pred_q <- data.frame(family$linkinv(pred_lin + outer(pred_lin_sd, qnorm(qs), "*")))
  names(pred_q) <- paste0("q", qs)
  pred_ms <- data.frame(pred = pred, pred_sd = pred_sd)

  pred0_ms <- pred0_q <- NULL
  if (has0) {
    X0tv    <- X0[, tv_cols, drop = FALSE]
    tvpart0 <- if (has_tv) rowSums(X0tv * tvbeta[pn0$tk, , drop = FALSE]) else rep(0, n0)
    off0 <- .dglm_clip_l(f0_obs, family) + tvpart0 + offset0
    cov0 <- if (ncv > 0) as.data.frame(X0[, const_cov, drop = FALSE]) else NULL
    dat0 <- data.frame(.off = off0); if (!is.null(cov0)) { names(cov0) <- xname[const_cov]; dat0 <- cbind(dat0, cov0) }
    pred0     <- predict(gmod, newdata = dat0, type = "response")
    pred0_lin <- predict(gmod, newdata = dat0, type = "link")
    Xg0 <- if (ncv > 0) cbind(1, X0[, const_cov, drop = FALSE]) else matrix(1, n0, 1)
    pred0_lin_sd <- sqrt(pmax(rowSums((Xg0 %*% Vbeta) * Xg0) + tvvar(X0tv, pn0$tk) +
                              pmin(tau * rowSums(Z0_sd^2), sill), 0))
    pred0_sd  <- abs(family$mu.eta(pred0_lin)) * pred0_lin_sd
    pred0_q <- data.frame(family$linkinv(pred0_lin + outer(pred0_lin_sd, qnorm(qs), "*")))
    names(pred0_q) <- paste0("q", qs)
    pred0_ms <- data.frame(pred = pred0, pred_sd = pred0_sd)
  }

  ## ---- spatial-process objects (per scale)
  Zdf <- Zsd_df <- Z0df <- Z0sd_df <- NULL
  if (length(bands) > 0) {
    Zdf <- as.data.frame(Z); Zsd_df <- as.data.frame(Z_sd)
    names(Zdf) <- names(Zsd_df) <- paste0("scale", seq_along(bands))
    if (has0) {
      Z0df <- as.data.frame(Z0); Z0sd_df <- as.data.frame(Z0_sd)
      names(Z0df) <- names(Z0sd_df) <- paste0("scale", seq_along(bands))
    }
  }

  ## ---- time-varying coefficients (smoothed value and per-time SD), if any
  beta_tv <- beta_tv_sd <- NULL
  if (has_tv) {
    d_tv <- length(tv_cols)
    beta_tv    <- as.data.frame(matrix(tvbeta, nrow = T, ncol = d_tv))
    beta_tv_sd <- as.data.frame(matrix(unlist(lapply(tvV, function(V) sqrt(pmax(diag(V), 0)))),
                                       nrow = T, ncol = d_tv, byrow = TRUE))
    names(beta_tv) <- names(beta_tv_sd) <- xname[tv_cols]
    beta_tv$time <- beta_tv_sd$time <- lev_work
  }

  ## ---- sd summary
  if (length(bands) > 0) {
    elements <- c("xb", paste0("spatial_scale", seq_along(bands)))
    standard_deviation <- c(sd(Xg %*% beta_int), apply(Z, 2, sd))
  } else { elements <- "xb"; standard_deviation <- sd(Xg %*% beta_int) }
  if (has_tv) {
    elements <- c(elements, paste0("tv_", xname[tv_cols]))
    standard_deviation <- c(standard_deviation, apply(tvbeta, 2, sd))
  }
  sd_summary <- data.frame(elements, standard_deviation); row.names(sd_summary) <- NULL

  ## ---- validation error statistics (holdout from mod_hv)
  idt <- mod_hv$id_train
  yt <- y[-idt]; yp <- pred[-idt]
  yp_c <- switch(family$family,
                 binomial = pmin(pmax(yp, 1e-6), 1 - 1e-6),
                 poisson  = pmax(yp, 1e-8),
                 yp)
  gmod_null <- stats::glm(yt ~ 1, family = family)
  gmod_fix  <- stats::glm(yt ~ 0 + offset(family$linkfun(yp_c)), family = family)
  r2  <- 1 - gmod_fix$deviance / gmod_null$null.deviance
  rmse <- sqrt(mean((yt - yp)^2)); mae <- abs(mean(yt - yp))
  e_summary <- data.frame(stat = c("validation_Pseudo-R2", "validation_RMSE", "validation_MAE"),
                          value = c(r2, rmse, mae))

  other <- list(n = n, n0 = if (has0) n0 else NA, nx = nx, y = y,
                coords = coords, coords0 = coords0, rho = rho, Q = Q,
                kernel = kernel, beta_int_vmat = Vbeta, loss_hv = mod_hv$loss_hv,
                tau = tau, tv_cols = tv_cols, q_tvc = q_tvc,
                time = time, time0 = if (has0) time0 else NULL,
                time_levels = lev_work, time_levels_train = lev,
                robust_se = robust_se, se_blocks = G_block)
  result <- list(beta = beta_summ, beta_tv = beta_tv, beta_tv_sd = beta_tv_sd,
                 sd_summary = sd_summary, e_summary = e_summary,
                 pred = pred_ms, pred0 = pred0_ms, pred_q = pred_q, pred0_q = pred0_q,
                 bands = bands, Z = Zdf, Z_sd = Zsd_df, Z0 = Z0df, Z0_sd = Z0sd_df,
                 other = other, call = match.call())
  class(result) <- "cf_dglm"
  result
}

#' @noRd
#' @export
print.cf_dglm <- function(x, ...) {
  cat("Call:\n"); print(x$call)
  cat("\n---- Coefficients -------------------------------------\n")
  print(x$beta)
  cat("\n---- Standard deviations (model elements) -------------\n")
  print(x$sd_summary)
  cat("\n---- Error statistics ---------------------------------\n")
  print(x$e_summary)
  invisible(x)
}
