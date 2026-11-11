## Diagnostic for approximate W-orthogonality between the fitted CF-STM field
## and the fixed-effect design. Run from the sim_repro root.

suppressMessages({
  library(FNN)
  library(fields)
  library(dbscan)
  library(withr)
})
source("engine/internal_utils_dglm.R")
source("engine/cf_dglm_hv.R")
source("engine/cf_dglm.R")
.dglm_load_cpp()

NSEED <- as.integer(Sys.getenv("NSEED", "10"))
RHO   <- as.numeric(Sys.getenv("RHO", "0.7"))
OUT   <- Sys.getenv("OUT", "results/orthogonality_rho07.csv")

gen <- function(n_pop = 400, T = 40, rho = 0.7, range = 0.10,
                fam = "gaussian", seed = 1) {
  set.seed(seed)
  s <- cbind(runif(n_pop), runif(n_pop))
  D <- fields::rdist(s)
  W <- exp(-D / range)
  W <- W / sqrt(rowSums(W^2))
  b <- matrix(0, T, n_pop)
  b[1, ] <- as.numeric(W %*% rnorm(n_pop)) / sqrt(1 - rho^2)
  for (t in 2:T) b[t, ] <- rho * b[t - 1, ] + as.numeric(W %*% rnorm(n_pop))
  x1 <- x2 <- matrix(0, T, n_pop)
  for (t in 1:T) {
    x1[t, ] <- as.numeric(W %*% rnorm(n_pop))
    x2[t, ] <- as.numeric(W %*% rnorm(n_pop))
  }
  pt <- rep(1:T, each = n_pop)
  loc <- rep(1:n_pop, T)
  coords <- s[loc, ]
  al <- b[cbind(pt, loc)]
  X1 <- x1[cbind(pt, loc)]
  X2 <- x2[cbind(pt, loc)]
  if (fam == "gaussian")
    y <- 2 * X1 - 1.5 * X2 + al + rnorm(n_pop * T, sd = 0.5)
  else if (fam == "poisson")
    y <- rpois(n_pop * T, exp(0.5 + 0.3 * X1 - 0.2 * X2 + 0.7 * al / sd(al)))
  else
    y <- rbinom(n_pop * T, 1, plogis(0.5 * X1 - 0.4 * X2 + 1.2 * al / sd(al)))
  list(y = y, coords = coords, pt = pt, x = cbind(x1 = X1, x2 = X2), latent = al)
}

fam_obj <- function(f)
  if (f == "gaussian") gaussian() else if (f == "poisson") poisson() else binomial()

wcov <- function(a, b, w) {
  a <- a - weighted.mean(a, w)
  b <- b - weighted.mean(b, w)
  sum(w * a * b) / sqrt(sum(w * a^2) * sum(w * b^2))
}

## Selected-model local linearization of the fitted-field map. With the chosen
## bands, knots, rho and Q held fixed, apply the same coarse-to-fine cascade to
## an infinitesimal working-residual perturbation v. Thus H(v) approximates the
## change in fhat caused by changing X beta by -v. IRLS-weight derivatives and
## discrete scale-selection changes are deliberately excluded.
apply_selected_H <- function(v, w, coords, time, md, seed) {
  if (length(md$bands) == 0) return(rep(0, length(v)))
  lev <- md$other$time_levels
  pn <- .dglm_panel(coords, time, time_levels = lev)
  setups <- lapply(seq_along(md$bands), function(k)
    .dglm_scale_setup(pn$C, md$bands[k], md$other$kernel, seed, Cpr = NULL))
  f <- matrix(0, pn$nL, pn$T)
  obs <- cbind(pn$lk, pn$tk)
  for (k in seq_along(md$bands)) {
    resid <- v - f[obs]
    Rp <- matrix(NA_real_, pn$nL, pn$T); Rp[obs] <- resid
    Wp <- matrix(NA_real_, pn$nL, pn$T); Wp[obs] <- w
    sc <- .dglm_scale_apply(setups[[k]], Rp, Wp, md$other$rho,
                            md$other$Q, predict = FALSE)
    f <- f + sc$Ftr
  }
  out <- f[obs]
  out - weighted.mean(out, w)
}

diagnose_one <- function(seed, fam) {
  d <- gen(fam = fam, seed = seed, rho = RHO)
  ff <- fam_obj(fam)
  invisible(utils::capture.output(
    mh <- cf_dglm_hv(y = d$y, x = d$x, coords = d$coords, time = d$pt,
                     family = ff, seed = seed)
  ))
  md <- cf_dglm(y = d$y, x = d$x, coords = d$coords, time = d$pt,
                mod_hv = mh)
  fhat <- if (is.null(md$Z)) rep(0, length(d$y)) else rowSums(md$Z)
  X <- cbind(Intercept = 1, d$x)
  bhat <- md$beta[c("Intercept", "x1", "x2"), "coef"]
  eta <- drop(X %*% bhat) + fhat
  zw <- .dglm_work(ff, eta, d$y, rep(0, length(d$y)))
  w <- zw$w
  econd <- zw$z - drop(X %*% bhat) - fhat

  ## Weighted correlations are the scale-free form of x_j' W f after centering.
  c1 <- wcov(d$x[, "x1"], fhat, w)
  c2 <- wcov(d$x[, "x2"], fhat, w)
  ct1 <- wcov(d$x[, "x1"], d$latent, w)
  ct2 <- wcov(d$x[, "x2"], d$latent, w)

  ## Fraction of fitted-field weighted variation explained by the covariates.
  fw <- stats::lm.wfit(X, fhat, w)
  f0 <- fhat - weighted.mean(fhat, w)
  r2 <- 1 - sum(w * fw$residuals^2) / sum(w * f0^2)

  ## Contribution of the restored field to the BLOCK slope scores in Eq. (13).
  ## Using the global score would be uninformative: the final offset-GLM makes
  ## the conditional score sum exactly zero by construction.
  Xc <- sweep(d$x, 2, apply(d$x, 2, weighted.mean, w = w), "-")
  Xcs <- sweep(Xc, 2, sqrt(colSums(w * Xc^2)), "/")
  rng <- as.numeric(stats::quantile(md$bands, 0.5))
  span <- apply(d$coords, 2, function(z) diff(range(z)))
  Gxy <- pmax(2L, pmin(8L, as.integer(floor(span / rng))))
  qx <- unique(stats::quantile(d$coords[, 1], seq(0, 1, length.out = Gxy[1] + 1)))
  qy <- unique(stats::quantile(d$coords[, 2], seq(0, 1, length.out = Gxy[2] + 1)))
  blk <- interaction(cut(d$coords[, 1], qx, include.lowest = TRUE),
                     cut(d$coords[, 2], qy, include.lowest = TRUE), drop = TRUE)
  S_field <- rowsum(Xcs * (w * fhat), blk)
  S_cond  <- rowsum(Xcs * (w * econd), blk)
  score_ratio <- sqrt(sum(S_field^2)) /
    max(sqrt(sum((S_field + S_cond)^2)), .Machine$double.eps)
  score_cosine <- sum(S_field * S_cond) /
    sqrt(sum(S_field^2) * sum(S_cond^2))

  ## Local sensitivity: X' W H X after W-normalizing and centering X columns.
  ## Its largest singular value is the strongest covariate direction that leaks
  ## through the selected latent-field operator and back into the beta score.
  HX <- sapply(seq_len(ncol(Xcs)), function(j)
    apply_selected_H(Xcs[, j], w, d$coords, d$pt, md, seed))
  local_mat <- crossprod(Xcs, w * HX)
  local_leak <- max(svd(local_mat, nu = 0, nv = 0)$d)
  local_response <- max(sqrt(colSums(w * HX^2)))

  data.frame(seed = seed, family = fam, nscale = length(md$bands),
             corr_x1_fhat = c1, corr_x2_fhat = c2,
             max_abs_corr_fhat = max(abs(c1), abs(c2)),
             field_on_x_R2 = r2, score_field_fraction = score_ratio,
             score_field_cond_cosine = score_cosine,
             local_score_leakage = local_leak,
             local_field_response = local_response,
             corr_x1_truefield = ct1, corr_x2_truefield = ct2)
}

tasks <- expand.grid(seed = seq_len(NSEED),
                     family = c("gaussian", "poisson", "binary"),
                     stringsAsFactors = FALSE)
rows <- lapply(seq_len(nrow(tasks)), function(i) {
  cat(sprintf("seed=%d family=%s\n", tasks$seed[i], tasks$family[i]))
  tryCatch(diagnose_one(tasks$seed[i], tasks$family[i]),
           error = function(e) {
             warning(conditionMessage(e))
             NULL
           })
})
res <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
write.csv(res, OUT, row.names = FALSE)

cat("\nSummary by family (median [90th percentile])\n")
for (fm in unique(res$family)) {
  z <- res[res$family == fm, ]
  cat(sprintf("%-8s |corr|max %.3f [%.3f], R2 %.3f [%.3f], score fraction %.3f [%.3f], local leakage %.3f [%.3f]\n",
              fm,
              median(z$max_abs_corr_fhat), quantile(z$max_abs_corr_fhat, .9),
              median(z$field_on_x_R2), quantile(z$field_on_x_R2, .9),
              median(z$score_field_fraction), quantile(z$score_field_fraction, .9),
              median(z$local_score_leakage), quantile(z$local_score_leakage, .9)))
}
cat(sprintf("\nwrote %d rows -> %s\n", nrow(res), OUT))
