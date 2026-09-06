## PAPER TABLE: Table A2 (Appendix 2) -- predictive accuracy under the Gaussian-process DGP  (see ../README.md)
###############################################################################
## Benchmark: GLM, GAM, cf_dglm, KFAS-A (single K=120 / multi K=240), sdmTMB.
## Data: n_pop=400, T=40, Gaussian, Poisson and Binomial. 100 reps each.
## DGP: identical to Table 1 EXCEPT that the latent field and the covariates are
## drawn from a Matern Gaussian process (smoothness 1, marginal variance 1) that
## evolves as a temporal AR(1) -- i.e. exactly the covariance the SPDE assumes.
## Metrics per seed: RMSE, CRPS (vs observed test y), computation time (s).
## Output: CSV (seed, family, method, rmse, crps, time).
###############################################################################
suppressMessages({
  library(FNN); library(fields); library(dbscan); library(withr)
  library(KFAS); library(sdmTMB); library(mgcv); library(scoringRules); library(parallel)
})
source("engine/internal_utils_dglm.R"); source("engine/cf_dglm_hv.R"); source("engine/cf_dglm.R")
.dglm_load_cpp()  # precompile fused C++ kernel in parent before mclapply forks (avoids per-worker build race)

NSEED <- as.integer(Sys.getenv("NSEED", "100"))
NCORE <- as.integer(Sys.getenv("NCORE", "20"))
OUT   <- Sys.getenv("OUT", "results/tableA2_gp_rho07.csv")
RHO   <- as.numeric(Sys.getenv("RHO","0.7"))

## ---- data: spatiotemporal AR(1) latent field + covariates ----
gen <- function(n_pop=400, T=40, rho=0.7, range=0.10, fam="gaussian", seed=1){
  set.seed(seed)
  s <- cbind(runif(n_pop), runif(n_pop)); D <- fields::rdist(s)
  Sig <- fields::Matern(D, range=range, smoothness=1)      # Matern GP (nu = 1), marginal var 1
  Rc  <- chol(Sig + diag(1e-8, n_pop)); gp <- function() as.numeric(crossprod(Rc, rnorm(n_pop)))
  b <- matrix(0,T,n_pop); b[1,] <- gp()/sqrt(1-rho^2)      # AR(1) in time, GP innovations
  for(t in 2:T) b[t,] <- rho*b[t-1,] + gp()
  x1 <- matrix(0,T,n_pop); x2 <- matrix(0,T,n_pop)
  for(t in 1:T){ x1[t,]<-gp(); x2[t,]<-gp() }              # covariates from the same GP
  pt <- rep(1:T,each=n_pop); loc <- rep(1:n_pop,T); coords <- s[loc,]
  al <- b[cbind(pt,loc)]; X1 <- x1[cbind(pt,loc)]; X2 <- x2[cbind(pt,loc)]
  if(fam=="gaussian"){
    y <- 2*X1 - 1.5*X2 + al + rnorm(n_pop*T, sd=0.5)
  } else if(fam=="poisson"){                 # poisson: moderate counts
    y <- rpois(n_pop*T, exp(0.5 + 0.3*X1 - 0.2*X2 + 0.7*al/sd(al)))
  } else {                                   # binomial (Bernoulli)
    y <- rbinom(n_pop*T, 1, plogis(0.0 + 0.5*X1 - 0.4*X2 + 1.2*al/sd(al)))
  }
  list(y=y, coords=coords, pt=pt, loc=loc, x=cbind(x1=X1,x2=X2), n_pop=n_pop, T=T)
}
split <- function(d, seed=1){
  set.seed(seed+10000); te <- sort(sample.int(d$n_pop, round(d$n_pop*0.30))); ix <- d$loc %in% te
  list(ix=ix,
       y_tr=d$y[!ix], x_tr=d$x[!ix,,drop=FALSE], c_tr=d$coords[!ix,], p_tr=d$pt[!ix], loc_tr=d$loc[!ix],
       y_te=d$y[ix],  x_te=d$x[ix,,drop=FALSE],  c_te=d$coords[ix,],  p_te=d$pt[ix],  loc_te=d$loc[ix])
}
fam_obj <- function(f) if(f=="gaussian") gaussian() else if(f=="poisson") poisson() else binomial()

## ---- metrics: mu = predicted response mean at test; tr_res_sd = in-sample residual SD ----
score <- function(s, mu, tr_res_sd, fam){
  rmse <- sqrt(mean((s$y_te - mu)^2))
  crps <- if(fam=="gaussian") mean(scoringRules::crps_norm(s$y_te, mu, pmax(tr_res_sd,1e-6)))
          else if(fam=="poisson") mean(scoringRules::crps_pois(s$y_te, pmax(mu,1e-8)))
          else { p <- pmin(pmax(mu,1e-6),1-1e-6); mean(scoringRules::crps_binom(s$y_te, size=1, prob=p)) }
  c(rmse=rmse, crps=crps)
}

## ---------------- method implementations (return list(mu, mu_tr, ok)) -------
m_glm <- function(s, fam){
  g <- glm(y_tr ~ x1 + x2, data=data.frame(y_tr=s$y_tr, s$x_tr), family=fam_obj(fam))
  list(mu=predict(g, data.frame(s$x_te), type="response"),
       mu_tr=predict(g, type="response"))
}
m_gam <- function(s, fam){
  dtr <- data.frame(y_tr=s$y_tr, x1=s$x_tr[,1], x2=s$x_tr[,2], cx=s$c_tr[,1], cy=s$c_tr[,2], tt=s$p_tr)
  dte <- data.frame(x1=s$x_te[,1], x2=s$x_te[,2], cx=s$c_te[,1], cy=s$c_te[,2], tt=s$p_te)
  g <- mgcv::gam(y_tr ~ x1 + x2 + te(cx,cy,tt, d=c(2,1), k=c(40,8), bs=c("tp","cr")),
                 data=dtr, family=fam_obj(fam), method="REML")
  list(mu=as.numeric(predict(g, dte, type="response")),
       mu_tr=as.numeric(predict(g, type="response")))
}
m_cfdglm <- function(s, fam, seed){
  mh <- cf_dglm_hv(y=s$y_tr, x=s$x_tr, coords=s$c_tr, time=s$p_tr, family=fam_obj(fam), seed=seed)
  md <- cf_dglm(y=s$y_tr, x=s$x_tr, coords=s$c_tr, time=s$p_tr,
                x0=s$x_te, coords0=s$c_te, time0=s$p_te, mod_hv=mh)
  list(mu=md$pred0$pred, mu_tr=md$pred$pred)
}
## KFAS construction A: state = K knot coefficients (vector AR1), Z = kernel loadings.
## gaussian: model residual y - Xbeta. poisson: model z=log(y+0.5) (Gaussian approx), exp back.
m_kfas <- function(s, fam, bands, nk, seed){
  if(fam=="gaussian"){ resp <- s$y_tr; link_inv <- identity; h0 <- 0.25 }
  else if(fam=="poisson"){ resp <- log(s$y_tr + 0.5); link_inv <- exp; h0 <- 1 }
  else { resp <- log((s$y_tr + 0.5)/(1.5 - s$y_tr)); link_inv <- plogis; h0 <- 1 }
  beta <- lm.fit(cbind(1, s$x_tr), resp)$coefficients
  r <- resp - drop(cbind(1, s$x_tr) %*% beta)
  ul <- sort(unique(s$loc_tr)); n <- length(ul); T <- max(s$p_tr); lk <- match(s$loc_tr, ul)
  Y <- matrix(NA, T, n); Y[cbind(s$p_tr, lk)] <- r
  Cu <- matrix(0,n,2); fi <- which(!duplicated(lk)); Cu[lk[fi],] <- s$c_tr[fi,]
  ute <- sort(unique(s$loc_te)); nt <- length(ute); lkt <- match(s$loc_te, ute)
  Ct <- matrix(0,nt,2); ft <- which(!duplicated(lkt)); Ct[lkt[ft],] <- s$c_te[ft,]
  Zl<-list(); Z0l<-list()
  for(j in seq_along(bands)){ kn <- with_seed(seed, kmeans(Cu, nk[j], iter.max=20)$centers)
    Zl[[j]] <- exp(-rdist(Cu,kn)/bands[j]); Z0l[[j]] <- exp(-rdist(Ct,kn)/bands[j]) }
  Z <- do.call(cbind,Zl); Z0 <- do.call(cbind,Z0l); K <- ncol(Z)
  mod <- SSModel(Y ~ -1 + SSMcustom(Z=Z, T=diag(K), R=diag(K), Q=diag(K),
                 a1=rep(0,K), P1=diag(K), state_names=paste0("a",1:K)), H=diag(n))
  upd <- function(p,model){ rho<-tanh(p[1]); q<-exp(p[2]); s2<-exp(p[3])
    model$T[,,1]<-rho*diag(K); model$Q[,,1]<-q*diag(K); model$P1[]<-(q/(1-rho^2))*diag(K); model$H[,,1]<-s2*diag(n); model }
  fit <- fitSSM(mod, inits=c(atanh(0.7),log(1),log(h0)), updatefn=upd, method="BFGS",
                control=list(maxit=20))
  ks <- KFS(fit$model, smoothing="state"); ah <- ks$alphahat
  predf <- function(Zb, lkz){ fld <- ah %*% t(Zb); out <- numeric(0); out }
  mk <- function(Zb, lkz, p_idx, x_idx){
    fld <- ah %*% t(Zb); v <- numeric(length(p_idx))
    for(t in 1:T){ it <- which(p_idx==t); if(length(it)) v[it] <- drop(cbind(1,x_idx[it,,drop=FALSE])%*%beta) + fld[t,lkz[it]] }
    link_inv(v)
  }
  list(mu = mk(Z0, lkt, s$p_te, s$x_te), mu_tr = mk(Z, lk, s$p_tr, s$x_tr))
}
m_sdmtmb <- function(s, fam, cutoff=0.05){
  tr <- data.frame(y=s$y_tr, x1=s$x_tr[,1], x2=s$x_tr[,2], cx=s$c_tr[,1], cy=s$c_tr[,2], pt=s$p_tr)
  te <- data.frame(x1=s$x_te[,1], x2=s$x_te[,2], cx=s$c_te[,1], cy=s$c_te[,2], pt=s$p_te)
  mesh <- make_mesh(tr, xy_cols=c("cx","cy"), cutoff=cutoff)
  fit <- sdmTMB(y ~ x1 + x2, data=tr, mesh=mesh, time="pt", spatiotemporal="ar1",
                spatial="on", family=fam_obj(fam), silent=TRUE)
  pr  <- predict(fit, newdata=te); prt <- predict(fit)
  inv <- if(fam=="gaussian") identity else if(fam=="poisson") exp else plogis
  list(mu=inv(pr$est), mu_tr=inv(prt$est))
}

## ---------------- run one (seed, family) task ----------------
run_task <- function(seed, fam){
  d <- gen(fam=fam, seed=seed, rho=RHO); s <- split(d, seed=seed)
  trsd <- function(mu_tr) sqrt(mean((s$y_tr - mu_tr)^2))
  methods <- list(
    GLM      = function() m_glm(s, fam),
    GAM      = function() m_gam(s, fam),
    cf_dglm  = function() m_cfdglm(s, fam, seed),
    KFAS_K120= function() m_kfas(s, fam, bands=0.12,               nk=120,            seed=seed),
    KFAS_K240= function() m_kfas(s, fam, bands=c(0.20,0.10,0.05),  nk=c(40,80,120),   seed=seed),
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
      sc <- score(s, res$mu, trsd(res$mu_tr), fam)
      rows[[nm]] <- data.frame(seed=seed, family=fam, method=nm, rmse=sc["rmse"], crps=sc["crps"], time=el)
    }
  }
  do.call(rbind, rows)
}

## ---------------- driver (parallel over seed x family) ----------------
tasks <- expand.grid(seed=1:NSEED, fam=c("gaussian","poisson","binary"), stringsAsFactors=FALSE)
cat(sprintf("Running %d tasks on %d cores -> %s\n", nrow(tasks), NCORE, OUT))
res_list <- mclapply(seq_len(nrow(tasks)), function(i)
                       tryCatch(run_task(tasks$seed[i], tasks$fam[i]), error=function(e) NULL),
                     mc.cores=NCORE, mc.preschedule=FALSE)
res <- do.call(rbind, res_list[!sapply(res_list, is.null)])
rownames(res) <- NULL
write.csv(res, OUT, row.names=FALSE)
cat(sprintf("Wrote %d rows to %s\n", nrow(res), OUT))

## summary
cat("\n=== mean over seeds (NA-dropped) ===\n")
agg <- aggregate(cbind(rmse,crps,time) ~ family + method, data=res, FUN=function(z) mean(z, na.rm=TRUE))
print(agg[order(agg$family, agg$rmse),], digits=4, row.names=FALSE)
