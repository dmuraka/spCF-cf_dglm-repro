## PAPER TABLE: Table 2 (predictive accuracy, irregular moving sites)  (see ../README.md)
###############################################################################
## MOVING-LOCATION benchmark: observation sites are RESAMPLED every time point
## (no persistent monitoring network). Same parameters/sizes as the fixed-site
## benchmark: n_pop=400 fresh sites per time, T=40, rho via RHO, range=0.10.
## Families: gaussian, poisson, binary. Methods: GLM, GAM, cf_dglm, KFAS-A
## (K=120 / K=240, now with TIME-VARYING loadings), sdmTMB. 100 reps each.
## Metrics vs observed test y: RMSE, CRPS, time(s). Output CSV.
###############################################################################
suppressMessages({
  library(FNN); library(fields); library(dbscan); library(withr)
  library(KFAS); library(sdmTMB); library(mgcv); library(scoringRules); library(parallel)
})
source("engine/internal_utils_dglm.R"); source("engine/cf_dglm_hv.R"); source("engine/cf_dglm.R")
.dglm_load_cpp()  # precompile fused C++ kernel in parent before mclapply forks (avoids per-worker build race)

NSEED <- as.integer(Sys.getenv("NSEED", "100"))
NCORE <- as.integer(Sys.getenv("NCORE", "20"))
OUT   <- Sys.getenv("OUT", "results/table2_moving_rho07.csv")
RHO   <- as.numeric(Sys.getenv("RHO","0.7"))

## ---- moving-location DGP -----------------------------------------------------
## Latent AR(1) field lives on FIXED anchors (temporal coherence); evaluated at
## fresh per-time sites by kernel interpolation. Covariates regenerated on each
## time's sites (independent across time, as in the fixed-site DGP).
gen <- function(n_pop=400, T=40, rho=0.7, range=0.10, fam="gaussian", seed=1){
  set.seed(seed)
  na <- n_pop
  anc <- cbind(runif(na), runif(na))
  Wa <- exp(-fields::rdist(anc)/range); Wa <- Wa/sqrt(rowSums(Wa^2))
  b  <- matrix(0,T,na); b[1,] <- as.numeric(Wa%*%rnorm(na))/sqrt(1-rho^2)
  for(t in 2:T) b[t,] <- rho*b[t-1,] + as.numeric(Wa%*%rnorm(na))
  N <- n_pop*T
  coords <- matrix(0,N,2); al <- numeric(N); X1 <- numeric(N); X2 <- numeric(N)
  pt <- rep(1:T, each=n_pop)
  for(t in 1:T){
    st <- cbind(runif(n_pop), runif(n_pop))
    Wi <- exp(-fields::rdist(st, anc)/range); Wi <- Wi/rowSums(Wi)         # field interpolation
    Wx <- exp(-fields::rdist(st)/range);      Wx <- Wx/sqrt(rowSums(Wx^2)) # covariate field
    idx <- ((t-1)*n_pop+1):(t*n_pop)
    coords[idx,] <- st
    al[idx] <- as.numeric(Wi %*% b[t,])
    X1[idx] <- as.numeric(Wx %*% rnorm(n_pop)); X2[idx] <- as.numeric(Wx %*% rnorm(n_pop))
  }
  al <- al * (sd(as.vector(b)) / sd(al))            # match field marginal sd to anchor field
  loc <- seq_len(N)                                  # every cell is a unique location
  if(fam=="gaussian")      y <- 2*X1 - 1.5*X2 + al + rnorm(N, sd=0.5)
  else if(fam=="poisson")  y <- rpois(N, exp(0.5 + 0.3*X1 - 0.2*X2 + 0.7*al/sd(al)))
  else                     y <- rbinom(N,1, plogis(0.0 + 0.5*X1 - 0.4*X2 + 1.2*al/sd(al)))
  list(y=y, coords=coords, pt=pt, loc=loc, x=cbind(x1=X1,x2=X2), n_pop=n_pop, T=T)
}
split <- function(d, seed=1){
  set.seed(seed+10000); N <- length(d$y)
  ix <- logical(N); ix[sample.int(N, round(N*0.30))] <- TRUE     # random 30% of cells held out
  list(ix=ix, d=d,
       y_tr=d$y[!ix], x_tr=d$x[!ix,,drop=FALSE], c_tr=d$coords[!ix,], p_tr=d$pt[!ix], loc_tr=d$loc[!ix],
       y_te=d$y[ix],  x_te=d$x[ix,,drop=FALSE],  c_te=d$coords[ix,],  p_te=d$pt[ix],  loc_te=d$loc[ix])
}
fam_obj <- function(f) if(f=="gaussian") gaussian() else if(f=="poisson") poisson() else binomial()
inv_link <- function(f) if(f=="gaussian") identity else if(f=="poisson") exp else plogis

score <- function(s, fam, mu, tr_res_sd){
  rmse <- sqrt(mean((s$y_te - mu)^2))
  crps <- if(fam=="gaussian") mean(scoringRules::crps_norm(s$y_te, mu, pmax(tr_res_sd,1e-6)))
          else if(fam=="poisson") mean(scoringRules::crps_pois(s$y_te, pmax(mu,1e-8)))
          else { p <- pmin(pmax(mu,1e-6),1-1e-6); mean(scoringRules::crps_binom(s$y_te, size=1, prob=p)) }
  c(rmse=rmse, crps=crps)
}

## ---------------- methods (return list(mu, mu_tr)) ---------------------------
m_glm <- function(s, fam){
  g <- glm(y_tr ~ x1 + x2, data=data.frame(y_tr=s$y_tr, s$x_tr), family=fam_obj(fam))
  list(mu=predict(g, data.frame(s$x_te), type="response"), mu_tr=predict(g, type="response"))
}
m_gam <- function(s, fam){
  dtr <- data.frame(y_tr=s$y_tr, x1=s$x_tr[,1], x2=s$x_tr[,2], cx=s$c_tr[,1], cy=s$c_tr[,2], tt=s$p_tr)
  dte <- data.frame(x1=s$x_te[,1], x2=s$x_te[,2], cx=s$c_te[,1], cy=s$c_te[,2], tt=s$p_te)
  g <- mgcv::gam(y_tr ~ x1 + x2 + te(cx,cy,tt, d=c(2,1), k=c(40,8), bs=c("tp","cr")),
                 data=dtr, family=fam_obj(fam), method="REML")
  list(mu=as.numeric(predict(g, dte, type="response")), mu_tr=as.numeric(predict(g, type="response")))
}
m_cfdglm <- function(s, fam, seed){
  mh <- cf_dglm_hv(y=s$y_tr, x=s$x_tr, coords=s$c_tr, time=s$p_tr, family=fam_obj(fam), seed=seed)
  md <- cf_dglm(y=s$y_tr, x=s$x_tr, coords=s$c_tr, time=s$p_tr,
                x0=s$x_te, coords0=s$c_te, time0=s$p_te, mod_hv=mh)
  list(mu=md$pred0$pred, mu_tr=md$pred$pred)
}
## KFAS construction A with TIME-VARYING loadings: balanced count n_pop per time,
## Z[,,t] = kernel(sites_t, knots), test cells set NA. State = K knot AR(1) coefs.
m_kfas <- function(s, fam, bands, nk, seed){
  d <- s$d; ix <- s$ix; n_pop <- d$n_pop; T <- d$T
  if(fam=="gaussian"){ resp <- d$y;                       inv <- identity; h0 <- 0.25 }
  else if(fam=="poisson"){ resp <- log(d$y + 0.5);        inv <- exp;      h0 <- 1 }
  else { resp <- log((d$y + 0.5)/(1.5 - d$y));            inv <- plogis;   h0 <- 1 }
  Xall <- cbind(1, d$x)
  beta <- lm.fit(Xall[!ix,,drop=FALSE], resp[!ix])$coefficients
  r <- resp - drop(Xall %*% beta)
  kn <- lapply(seq_along(bands), function(j) with_seed(seed, kmeans(d$coords[!ix,], nk[j], iter.max=20)$centers))
  K <- sum(nk)
  Zarr <- array(0, c(n_pop, K, T)); Y <- matrix(NA, T, n_pop)
  for(t in 1:T){
    idx <- ((t-1)*n_pop+1):(t*n_pop); st <- d$coords[idx,]
    Zarr[,,t] <- do.call(cbind, lapply(seq_along(bands), function(j) exp(-rdist(st, kn[[j]])/bands[j])))
    yt <- r[idx]; yt[ix[idx]] <- NA; Y[t,] <- yt
  }
  mod <- SSModel(Y ~ -1 + SSMcustom(Z=Zarr, T=diag(K), R=diag(K), Q=diag(K),
                 a1=rep(0,K), P1=diag(K), state_names=paste0("a",1:K)), H=diag(n_pop))
  upd <- function(p,model){ rho<-tanh(p[1]); q<-exp(p[2]); s2<-exp(p[3])
    model$T[,,1]<-rho*diag(K); model$Q[,,1]<-q*diag(K); model$P1[]<-(q/(1-rho^2))*diag(K); model$H[,,1]<-s2*diag(n_pop); model }
  fit <- fitSSM(mod, inits=c(atanh(0.7),log(1),log(h0)), updatefn=upd, method="BFGS", control=list(maxit=20))
  ks <- KFS(fit$model, smoothing="state"); ah <- ks$alphahat          # T x K
  mu <- numeric(length(d$y))
  for(t in 1:T){ idx <- ((t-1)*n_pop+1):(t*n_pop)
    mu[idx] <- inv(drop(Xall[idx,] %*% beta) + drop(Zarr[,,t] %*% ah[t,])) }
  list(mu = mu[ix], mu_tr = mu[!ix])
}
m_sdmtmb <- function(s, fam, cutoff=0.05){
  tr <- data.frame(y=s$y_tr, x1=s$x_tr[,1], x2=s$x_tr[,2], cx=s$c_tr[,1], cy=s$c_tr[,2], pt=s$p_tr)
  te <- data.frame(x1=s$x_te[,1], x2=s$x_te[,2], cx=s$c_te[,1], cy=s$c_te[,2], pt=s$p_te)
  mesh <- make_mesh(tr, xy_cols=c("cx","cy"), cutoff=cutoff)
  fit <- sdmTMB(y ~ x1 + x2, data=tr, mesh=mesh, time="pt", spatiotemporal="ar1",
                spatial="on", family=fam_obj(fam), silent=TRUE)
  inv <- inv_link(fam)
  list(mu=inv(predict(fit, newdata=te)$est), mu_tr=inv(predict(fit)$est))
}

run_task <- function(seed, fam){
  d <- gen(fam=fam, seed=seed, rho=RHO); s <- split(d, seed=seed)
  trsd <- function(mu_tr) sqrt(mean((s$y_tr - mu_tr)^2))
  methods <- list(
    GLM      = function() m_glm(s, fam),
    GAM      = function() m_gam(s, fam),
    cf_dglm  = function() m_cfdglm(s, fam, seed),
    KFAS_K120= function() m_kfas(s, fam, bands=0.12,              nk=120,          seed=seed),
    KFAS_K240= function() m_kfas(s, fam, bands=c(0.20,0.10,0.05), nk=c(40,80,120), seed=seed),
    sdmTMB   = function() m_sdmtmb(s, fam)
  )
  rows <- list()
  for(nm in names(methods)){
    t0 <- proc.time()["elapsed"]
    res <- tryCatch(suppressMessages(suppressWarnings({ utils::capture.output(rr <- methods[[nm]]()); rr })),
                    error=function(e) e)
    el <- as.numeric(proc.time()["elapsed"] - t0)
    if(inherits(res,"error") || is.null(res$mu) || any(!is.finite(res$mu))){
      rows[[nm]] <- data.frame(seed=seed, family=fam, method=nm, rmse=NA, crps=NA, time=el)
    } else {
      sc <- score(s, fam, res$mu, trsd(res$mu_tr))
      rows[[nm]] <- data.frame(seed=seed, family=fam, method=nm, rmse=sc["rmse"], crps=sc["crps"], time=el)
    }
  }
  do.call(rbind, rows)
}

tasks <- expand.grid(seed=1:NSEED, fam=c("gaussian","poisson","binary"), stringsAsFactors=FALSE)
cat(sprintf("MV methods: %d tasks on %d cores (rho=%.2f) -> %s\n", nrow(tasks), NCORE, RHO, OUT))
res_list <- mclapply(seq_len(nrow(tasks)), function(i)
                       tryCatch(run_task(tasks$seed[i], tasks$fam[i]), error=function(e) NULL),
                     mc.cores=NCORE, mc.preschedule=FALSE)
res <- do.call(rbind, res_list[!sapply(res_list, is.null)]); rownames(res) <- NULL
write.csv(res, OUT, row.names=FALSE)
cat(sprintf("Wrote %d rows to %s\n", nrow(res), OUT))
cat("\n=== mean over seeds (NA-dropped) ===\n")
agg <- aggregate(cbind(rmse,crps,time) ~ family + method, data=res, FUN=function(z) mean(z, na.rm=TRUE))
print(agg[order(agg$family, agg$rmse),], digits=4, row.names=FALSE)
