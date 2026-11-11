## Model-based predictive distributions: CRPS, interval coverage, PIT for y, and
## coverage for the latent linear predictor eta.  Settings of Tables 1 (panel),
## 2 (moving) and A2 (gp): n_pop = 400, T = 40, Gaussian / Poisson / binary,
## rho via RHO, 70/30 split by site (panel, gp) or by cell (moving).
## Methods: GLM, GAM, CF-STM (default sigma^2; Gaussian also with the optional
## generalized-df correction), SPDE (sdmTMB, link-scale draws via nsim), and, for the
## Gaussian response only, GSSM / GSSM-MS (KFAS, same construction as sim_repro; predictive
## variance Z V_t Z' + H from the smoothed state covariance).
## Per task (seed x family): metrics row per method -> CSV; observation-level
## predictions -> preds/<design>_rho<..>_<family>_<seed>.rds when SAVE_PREDS=1
## Env: DESIGN=panel|moving|gp  RHO  NSEED  SEEDS  FAMS  NCORE  NSIM  ENGINE_MOD  DGLM_CPP  OUTDIR  PREDDIR  SAVE_PREDS  METHODS (comma list; default GLM,GAM,CF-STM,CF-STM-gdf,SPDE)
suppressMessages({
  library(FNN); library(fields); library(dbscan); library(withr)
  library(sdmTMB); library(mgcv); library(scoringRules); library(parallel); library(KFAS)
})
EM <- Sys.getenv("ENGINE_MOD"); source(file.path(EM,"internal_utils_dglm.R")); source(file.path(EM,"cf_dglm_hv.R")); source(file.path(EM,"cf_dglm.R")); .dglm_load_cpp()
DESIGN <- Sys.getenv("DESIGN","panel"); RHO <- as.numeric(Sys.getenv("RHO","0.7"))
NSEED <- as.integer(Sys.getenv("NSEED","100")); NCORE <- as.integer(Sys.getenv("NCORE","30")); NSIM <- as.integer(Sys.getenv("NSIM","500"))
SEEDS <- if(nzchar(Sys.getenv("SEEDS"))) as.integer(strsplit(Sys.getenv("SEEDS"),",")[[1]]) else 1:NSEED
FAMS  <- if(nzchar(Sys.getenv("FAMS"))) strsplit(Sys.getenv("FAMS"),",")[[1]] else c("gaussian","poisson","binary")
OUTDIR <- Sys.getenv("OUTDIR","results"); PREDDIR <- Sys.getenv("PREDDIR", file.path(dirname(OUTDIR),"preds"))
SAVE_PREDS <- identical(Sys.getenv("SAVE_PREDS", "0"), "1")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
if(SAVE_PREDS) dir.create(PREDDIR, recursive=TRUE, showWarnings=FALSE)
METHODS <- if(nzchar(Sys.getenv("METHODS"))) strsplit(Sys.getenv("METHODS"),",")[[1]] else NULL
tag <- sprintf("%s_rho%02d", DESIGN, round(RHO*10))

## ---------------- DGPs (verbatim from sim_repro; eta = true linear predictor added) ----------------
resp_of <- function(fam, X1, X2, al, N){
  if(fam=="gaussian"){ eta <- 2*X1 - 1.5*X2 + al; y <- eta + rnorm(N, sd=0.5) }
  else if(fam=="poisson"){ eta <- 0.5 + 0.3*X1 - 0.2*X2 + 0.7*al/sd(al); y <- rpois(N, exp(eta)) }
  else { eta <- 0.0 + 0.5*X1 - 0.4*X2 + 1.2*al/sd(al); y <- rbinom(N, 1, plogis(eta)) }
  list(y=y, eta=eta)
}
gen_panel <- function(n_pop=400, T=40, rho=0.7, range=0.10, fam="gaussian", seed=1){
  set.seed(seed)
  s <- cbind(runif(n_pop), runif(n_pop)); D <- fields::rdist(s); W <- exp(-D/range); W <- W/sqrt(rowSums(W^2))
  b <- matrix(0,T,n_pop); b[1,] <- as.numeric(W%*%rnorm(n_pop))/sqrt(1-rho^2)
  for(t in 2:T) b[t,] <- rho*b[t-1,] + as.numeric(W%*%rnorm(n_pop))
  x1 <- matrix(0,T,n_pop); x2 <- matrix(0,T,n_pop)
  for(t in 1:T){ x1[t,]<-as.numeric(W%*%rnorm(n_pop)); x2[t,]<-as.numeric(W%*%rnorm(n_pop)) }
  pt <- rep(1:T,each=n_pop); loc <- rep(1:n_pop,T); coords <- s[loc,]
  al <- b[cbind(pt,loc)]; X1 <- x1[cbind(pt,loc)]; X2 <- x2[cbind(pt,loc)]
  r <- resp_of(fam, X1, X2, al, n_pop*T)
  list(y=r$y, eta=r$eta, coords=coords, pt=pt, loc=loc, x=cbind(x1=X1,x2=X2), n_pop=n_pop, T=T)
}
gen_gp <- function(n_pop=400, T=40, rho=0.7, range=0.10, fam="gaussian", seed=1){
  set.seed(seed)
  s <- cbind(runif(n_pop), runif(n_pop)); D <- fields::rdist(s)
  Sig <- fields::Matern(D, range=range, smoothness=1)
  Rc  <- chol(Sig + diag(1e-8, n_pop)); gp <- function() as.numeric(crossprod(Rc, rnorm(n_pop)))
  b <- matrix(0,T,n_pop); b[1,] <- gp()/sqrt(1-rho^2)
  for(t in 2:T) b[t,] <- rho*b[t-1,] + gp()
  x1 <- matrix(0,T,n_pop); x2 <- matrix(0,T,n_pop)
  for(t in 1:T){ x1[t,]<-gp(); x2[t,]<-gp() }
  pt <- rep(1:T,each=n_pop); loc <- rep(1:n_pop,T); coords <- s[loc,]
  al <- b[cbind(pt,loc)]; X1 <- x1[cbind(pt,loc)]; X2 <- x2[cbind(pt,loc)]
  r <- resp_of(fam, X1, X2, al, n_pop*T)
  list(y=r$y, eta=r$eta, coords=coords, pt=pt, loc=loc, x=cbind(x1=X1,x2=X2), n_pop=n_pop, T=T)
}
gen_moving <- function(n_pop=400, T=40, rho=0.7, range=0.10, fam="gaussian", seed=1){
  set.seed(seed)
  na <- n_pop; anc <- cbind(runif(na), runif(na))
  Wa <- exp(-fields::rdist(anc)/range); Wa <- Wa/sqrt(rowSums(Wa^2))
  b  <- matrix(0,T,na); b[1,] <- as.numeric(Wa%*%rnorm(na))/sqrt(1-rho^2)
  for(t in 2:T) b[t,] <- rho*b[t-1,] + as.numeric(Wa%*%rnorm(na))
  N <- n_pop*T; coords <- matrix(0,N,2); al <- numeric(N); X1 <- numeric(N); X2 <- numeric(N)
  pt <- rep(1:T, each=n_pop)
  for(t in 1:T){
    st <- cbind(runif(n_pop), runif(n_pop))
    Wi <- exp(-fields::rdist(st, anc)/range); Wi <- Wi/rowSums(Wi)
    Wx <- exp(-fields::rdist(st)/range);      Wx <- Wx/sqrt(rowSums(Wx^2))
    idx <- ((t-1)*n_pop+1):(t*n_pop); coords[idx,] <- st
    al[idx] <- as.numeric(Wi %*% b[t,])
    X1[idx] <- as.numeric(Wx %*% rnorm(n_pop)); X2[idx] <- as.numeric(Wx %*% rnorm(n_pop))
  }
  al <- al * (sd(as.vector(b)) / sd(al)); loc <- seq_len(N)
  r <- resp_of(fam, X1, X2, al, N)
  list(y=r$y, eta=r$eta, coords=coords, pt=pt, loc=loc, x=cbind(x1=X1,x2=X2), n_pop=n_pop, T=T)
}
split_site <- function(d, seed=1){
  set.seed(seed+10000); te <- sort(sample.int(d$n_pop, round(d$n_pop*0.30))); ix <- d$loc %in% te
  list(ix=ix, y_tr=d$y[!ix], x_tr=d$x[!ix,,drop=FALSE], c_tr=d$coords[!ix,], p_tr=d$pt[!ix],
       y_te=d$y[ix], eta_te=d$eta[ix], x_te=d$x[ix,,drop=FALSE], c_te=d$coords[ix,], p_te=d$pt[ix],
       loc_tr=d$loc[!ix], loc_te=d$loc[ix], d=d)
}
split_cell <- function(d, seed=1){
  set.seed(seed+10000); N <- length(d$y); ix <- logical(N); ix[sample.int(N, round(N*0.30))] <- TRUE
  list(ix=ix, y_tr=d$y[!ix], x_tr=d$x[!ix,,drop=FALSE], c_tr=d$coords[!ix,], p_tr=d$pt[!ix],
       y_te=d$y[ix], eta_te=d$eta[ix], x_te=d$x[ix,,drop=FALSE], c_te=d$coords[ix,], p_te=d$pt[ix],
       loc_tr=d$loc[!ix], loc_te=d$loc[ix], d=d)
}
gen <- switch(DESIGN, panel=gen_panel, gp=gen_gp, moving=gen_moving)
split <- if(DESIGN=="moving") split_cell else split_site
fam_obj <- function(f) if(f=="gaussian") gaussian() else if(f=="poisson") poisson() else binomial()

## ---------------- methods: return list(lin, lin_sd, mu, sigma, mu_tr) (link scale lin/lin_sd) ----------------
m_glm <- function(s, fam){
  g <- glm(y_tr ~ x1 + x2, data=data.frame(y_tr=s$y_tr, s$x_tr), family=fam_obj(fam))
  p <- predict(g, data.frame(s$x_te), se.fit=TRUE)
  list(lin=as.numeric(p$fit), lin_sd=as.numeric(p$se.fit), mu=as.numeric(predict(g, data.frame(s$x_te), type="response")),
       sigma=if(fam=="gaussian") sqrt(summary(g)$dispersion) else NA, mu_tr=as.numeric(predict(g, type="response")))
}
m_gam <- function(s, fam){
  dtr <- data.frame(y_tr=s$y_tr, x1=s$x_tr[,1], x2=s$x_tr[,2], cx=s$c_tr[,1], cy=s$c_tr[,2], tt=s$p_tr)
  dte <- data.frame(x1=s$x_te[,1], x2=s$x_te[,2], cx=s$c_te[,1], cy=s$c_te[,2], tt=s$p_te)
  g <- mgcv::gam(y_tr ~ x1 + x2 + te(cx,cy,tt, d=c(2,1), k=c(40,8), bs=c("tp","cr")), data=dtr, family=fam_obj(fam), method="REML")
  p <- predict(g, dte, se.fit=TRUE)
  list(lin=as.numeric(p$fit), lin_sd=as.numeric(p$se.fit), mu=as.numeric(predict(g, dte, type="response")),
       sigma=if(fam=="gaussian") sqrt(g$sig2) else NA, mu_tr=as.numeric(predict(g, type="response")))
}
m_cfstm <- function(s, fam, seed, mh, sigma_est="insample"){
  fo <- fam_obj(fam)
  md <- cf_dglm(y=s$y_tr, x=s$x_tr, coords=s$c_tr, time=s$p_tr, x0=s$x_te, coords0=s$c_te, time0=s$p_te, mod_hv=mh, sigma_est=sigma_est)
  mu <- md$pred0$pred; lin <- fo$linkfun(mu); lin_sd <- md$pred0$pred_sd / pmax(abs(fo$mu.eta(lin)), 1e-12)
  sig <- if(fam!="gaussian") NA else if(sigma_est=="gdf") md$other$sigma_gdf else md$other$sigma_insample
  list(lin=lin, lin_sd=lin_sd, mu=mu, sigma=sig, mu_tr=md$pred$pred, R=length(md$bands), tau=md$other$tau)
}
## GSSM / GSSM-MS: KFAS construction A of sim_repro (state = K knot AR(1) coefficients, Z = kernel
## loadings, OLS beta plugged in), Gaussian only.  Predictive mean and variance of the field at a
## test site from the smoothed state (alphahat, V): lin = x'beta + Z0 a_t, lin_sd^2 = Z0 V_t Z0';
## sigma^2 = H (estimated).  The uncertainty of the plug-in beta is ignored, as in sim_repro.
## Fixed sites (panel, gp): time-invariant Z (sim_repro table1/tableA2); moving: time-varying Z (table2).
kfas_fit <- function(Y, Zmod, K, nobs, h0=0.25){
  mod <- SSModel(Y ~ -1 + SSMcustom(Z=Zmod, T=diag(K), R=diag(K), Q=diag(K),
                 a1=rep(0,K), P1=diag(K), state_names=paste0("a",1:K)), H=diag(nobs))
  upd <- function(p,model){ rho<-tanh(p[1]); q<-exp(p[2]); s2<-exp(p[3])
    model$T[,,1]<-rho*diag(K); model$Q[,,1]<-q*diag(K); model$P1[]<-(q/(1-rho^2))*diag(K); model$H[,,1]<-s2*diag(nobs); model }
  fit <- fitSSM(mod, inits=c(atanh(0.7),log(1),log(h0)), updatefn=upd, method="BFGS", control=list(maxit=20))
  ks <- KFS(fit$model, smoothing="state")
  list(ah=ks$alphahat, V=ks$V, s2=fit$model$H[1,1,1])
}
m_kfas <- function(s, fam, bands, nk, seed){
  if(fam!="gaussian") return(NULL)      # transformed-response approximation: no predictive distribution for y
  if(DESIGN=="moving") return(m_kfas_moving(s, bands, nk, seed))
  Xtr <- cbind(1, s$x_tr); beta <- lm.fit(Xtr, s$y_tr)$coefficients; r <- s$y_tr - drop(Xtr %*% beta)
  ul <- sort(unique(s$loc_tr)); n <- length(ul); T <- max(s$p_tr); lk <- match(s$loc_tr, ul)
  Y <- matrix(NA, T, n); Y[cbind(s$p_tr, lk)] <- r
  Cu <- matrix(0,n,2); fi <- which(!duplicated(lk)); Cu[lk[fi],] <- s$c_tr[fi,]
  ute <- sort(unique(s$loc_te)); nt <- length(ute); lkt <- match(s$loc_te, ute)
  Ct <- matrix(0,nt,2); ft <- which(!duplicated(lkt)); Ct[lkt[ft],] <- s$c_te[ft,]
  Zl<-list(); Z0l<-list()
  for(j in seq_along(bands)){ kn <- with_seed(seed, kmeans(Cu, nk[j], iter.max=20)$centers)
    Zl[[j]] <- exp(-rdist(Cu,kn)/bands[j]); Z0l[[j]] <- exp(-rdist(Ct,kn)/bands[j]) }
  Z <- do.call(cbind,Zl); Z0 <- do.call(cbind,Z0l); K <- ncol(Z)
  f <- kfas_fit(Y, Z, K, n)
  fld0 <- f$ah %*% t(Z0); fld <- f$ah %*% t(Z)                     # T x nt, T x n
  fsd0 <- sapply(1:T, function(t) sqrt(pmax(rowSums((Z0 %*% f$V[,,t]) * Z0), 0)))   # nt x T
  lin <- drop(cbind(1, s$x_te) %*% beta) + fld0[cbind(s$p_te, lkt)]
  list(lin=lin, lin_sd=fsd0[cbind(lkt, s$p_te)], mu=lin, sigma=sqrt(f$s2),
       mu_tr=drop(Xtr %*% beta) + fld[cbind(s$p_tr, lk)])
}
m_kfas_moving <- function(s, bands, nk, seed){
  d <- s$d; ix <- s$ix; n_pop <- d$n_pop; T <- d$T
  Xall <- cbind(1, d$x); beta <- lm.fit(Xall[!ix,,drop=FALSE], d$y[!ix])$coefficients; r <- d$y - drop(Xall %*% beta)
  kn <- lapply(seq_along(bands), function(j) with_seed(seed, kmeans(d$coords[!ix,], nk[j], iter.max=20)$centers))
  K <- sum(nk); Zarr <- array(0, c(n_pop, K, T)); Y <- matrix(NA, T, n_pop)
  for(t in 1:T){
    idx <- ((t-1)*n_pop+1):(t*n_pop); st <- d$coords[idx,]
    Zarr[,,t] <- do.call(cbind, lapply(seq_along(bands), function(j) exp(-rdist(st, kn[[j]])/bands[j])))
    yt <- r[idx]; yt[ix[idx]] <- NA; Y[t,] <- yt
  }
  f <- kfas_fit(Y, Zarr, K, n_pop)
  lin <- numeric(length(d$y)); lsd <- numeric(length(d$y))
  for(t in 1:T){ idx <- ((t-1)*n_pop+1):(t*n_pop); Zt <- Zarr[,,t]
    lin[idx] <- drop(Xall[idx,] %*% beta) + drop(Zt %*% f$ah[t,]); lsd[idx] <- sqrt(pmax(rowSums((Zt %*% f$V[,,t]) * Zt), 0)) }
  list(lin=lin[ix], lin_sd=lsd[ix], mu=lin[ix], sigma=sqrt(f$s2), mu_tr=lin[!ix])
}
m_sdmtmb <- function(s, fam, cutoff=0.05){
  tr <- data.frame(y=s$y_tr, x1=s$x_tr[,1], x2=s$x_tr[,2], cx=s$c_tr[,1], cy=s$c_tr[,2], pt=s$p_tr)
  te <- data.frame(x1=s$x_te[,1], x2=s$x_te[,2], cx=s$c_te[,1], cy=s$c_te[,2], pt=s$p_te)
  mesh <- make_mesh(tr, xy_cols=c("cx","cy"), cutoff=cutoff)
  fit <- sdmTMB(y ~ x1 + x2, data=tr, mesh=mesh, time="pt", spatiotemporal="ar1", spatial="on", family=fam_obj(fam), silent=TRUE)
  pr <- predict(fit, newdata=te); sims <- predict(fit, newdata=te, nsim=NSIM)
  inv <- fam_obj(fam)$linkinv
  sig <- NA; if(fam=="gaussian"){ rp <- tidy(fit, "ran_pars"); sig <- rp$estimate[rp$term=="phi"] }
  list(lin=as.numeric(pr$est), lin_sd=apply(sims,1,sd), mu=inv(as.numeric(pr$est)), sigma=sig, mu_tr=inv(as.numeric(predict(fit)$est)))
}

## ---------------- scoring ----------------
pit_hist <- function(u) as.numeric(table(cut(u, seq(0,1,0.1), include.lowest=TRUE)))
score_all <- function(s, fam, r, seed){
  y <- s$y_te; n <- length(y)
  rmse <- sqrt(mean((y - r$mu)^2)); trsd <- sqrt(mean((s$y_tr - r$mu_tr)^2))
  crps_paper <- if(fam=="gaussian") mean(crps_norm(y, r$mu, pmax(trsd,1e-6))) else if(fam=="poisson") mean(crps_pois(y, pmax(r$mu,1e-8))) else mean(crps_binom(y, 1, pmin(pmax(r$mu,1e-6),1-1e-6)))
  ## latent linear predictor
  z <- qnorm(0.975); eta_cov <- mean(abs(s$eta_te - r$lin) <= z*r$lin_sd); eta_w <- mean(2*z*r$lin_sd); eta_rmse <- sqrt(mean((s$eta_te - r$lin)^2))
  if(fam=="gaussian"){
    sd_y <- sqrt(r$lin_sd^2 + r$sigma^2)
    crps <- mean(crps_norm(y, r$lin, sd_y)); u <- pnorm(y, r$lin, sd_y)
    lo95 <- r$lin - z*sd_y; hi95 <- r$lin + z*sd_y; z90 <- qnorm(0.95); lo90 <- r$lin - z90*sd_y; hi90 <- r$lin + z90*sd_y
    cov95 <- mean(y>=lo95 & y<=hi95); w95 <- mean(hi95-lo95); cov90 <- mean(y>=lo90 & y<=hi90); w90 <- mean(hi90-lo90); mu_pd <- r$lin
  } else {
    set.seed(seed + 777L)                                            # common random numbers across methods
    E <- matrix(rnorm(n*NSIM), n, NSIM) * r$lin_sd + r$lin
    if(fam=="poisson"){
      Y <- matrix(rpois(n*NSIM, exp(E)), n, NSIM)
      crps <- mean(crps_sample(y, Y))
      q <- t(apply(Y, 1, quantile, probs=c(0.025,0.975,0.05,0.95), type=1))
      cov95 <- mean(y>=q[,1] & y<=q[,2]); w95 <- mean(q[,2]-q[,1]); cov90 <- mean(y>=q[,3] & y<=q[,4]); w90 <- mean(q[,4]-q[,3])
      Fy <- rowMeans(Y <= y); Fy1 <- rowMeans(Y <= y-1); u <- Fy1 + runif(n)*(Fy - Fy1); mu_pd <- rowMeans(exp(E))
    } else {
      p <- rowMeans(plogis(E)); p <- pmin(pmax(p,1e-6),1-1e-6)
      crps <- mean(crps_binom(y, 1, p)); cov95 <- NA; w95 <- NA; cov90 <- NA; w90 <- NA
      Fy <- ifelse(y==1, 1, 1-p); Fy1 <- ifelse(y==1, 1-p, 0); u <- Fy1 + runif(n)*(Fy - Fy1); mu_pd <- p
    }
  }
  c(rmse=rmse, rmse_pd=sqrt(mean((y-mu_pd)^2)), crps_paper=crps_paper, crps=crps, cov90=cov90, width90=w90, cov95=cov95, width95=w95,
    eta_rmse=eta_rmse, eta_cov95=eta_cov, eta_width95=eta_w, sigma_hat=if(is.null(r$sigma)) NA else r$sigma, mean_lin_sd=mean(r$lin_sd),
    setNames(pit_hist(u), paste0("pit", 1:10)))
}

run_task <- function(seed, fam){
  d <- gen(fam=fam, seed=seed, rho=RHO); s <- split(d, seed=seed)
  mh <- NULL
  meths <- list(
    GLM   = function() m_glm(s, fam),
    GAM   = function() m_gam(s, fam),
    `CF-STM` = function(){ mh <<- cf_dglm_hv(y=s$y_tr, x=s$x_tr, coords=s$c_tr, time=s$p_tr, family=fam_obj(fam), seed=seed); m_cfstm(s, fam, seed, mh) },
    `CF-STM-gdf` = function() if(fam=="gaussian") m_cfstm(s, fam, seed, mh, "gdf") else NULL,
    SPDE  = function() m_sdmtmb(s, fam),
    GSSM     = function() m_kfas(s, fam, bands=0.12,              nk=120,          seed=seed),
    `GSSM-MS`= function() m_kfas(s, fam, bands=c(0.20,0.10,0.05), nk=c(40,80,120), seed=seed))
  rows <- list(); keep <- list(y=s$y_te, eta=s$eta_te, pt=s$p_te)
  meths <- meths[names(meths) %in% (if(is.null(METHODS)) c("GLM","GAM","CF-STM","CF-STM-gdf","SPDE") else METHODS)]
  for(nm in names(meths)){
    t0 <- proc.time()[["elapsed"]]
    r <- tryCatch(suppressMessages(suppressWarnings({ utils::capture.output(rr <- meths[[nm]]()); rr })), error=function(e) e)
    el <- proc.time()[["elapsed"]] - t0
    if(is.null(r)) next
    if(inherits(r,"error") || any(!is.finite(r$mu))){ message(nm, " seed ", seed, " ", fam, " failed: ", if(inherits(r,"error")) conditionMessage(r) else "non-finite"); next }
    sc <- score_all(s, fam, r, seed)
    rows[[nm]] <- data.frame(design=DESIGN, rho=RHO, seed=seed, family=fam, method=nm, time=el, t(sc),
                             R=if(is.null(r$R)) NA else r$R, tau=if(is.null(r$tau)) NA else r$tau)
    keep[[nm]] <- list(lin=r$lin, lin_sd=r$lin_sd, mu=r$mu, sigma=r$sigma)
  }
  if(SAVE_PREDS)
    saveRDS(keep, file.path(PREDDIR, sprintf("%s_%s_%03d.rds", tag, fam, seed)), compress="xz")
  do.call(rbind, rows)
}
tasks <- expand.grid(seed=SEEDS, fam=FAMS, stringsAsFactors=FALSE)
cat(sprintf("%s: %d tasks on %d cores\n", tag, nrow(tasks), NCORE))
res <- mclapply(seq_len(nrow(tasks)), function(i) tryCatch(run_task(tasks$seed[i], tasks$fam[i]), error=function(e){message("task ",i," failed: ",conditionMessage(e)); NULL}), mc.cores=NCORE, mc.preschedule=FALSE)
cat("failed tasks:", sum(sapply(res,is.null)), "\n")
res <- do.call(rbind, res[!sapply(res,is.null)]); rownames(res) <- NULL
write.csv(res, file.path(OUTDIR, sprintf("coverage_%s.csv", tag)), row.names=FALSE)
cat("\n== mean over seeds ==\n")
print(aggregate(cbind(rmse,crps_paper,crps,cov90,cov95,width95,eta_cov95,eta_width95,time)~family+method, res, mean), digits=3, row.names=FALSE)
