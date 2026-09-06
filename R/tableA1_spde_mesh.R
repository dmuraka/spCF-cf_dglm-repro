## PAPER TABLE: Table A1 (Appendix 1) -- sensitivity of the SPDE benchmark to the
## mesh resolution.  (see ../README.md)
## sdmTMB mesh-cutoff sweep on the Table-2 (moving-site) DGP.
## Compares RMSE / CRPS / time / mesh-nodes across cutoff values; CF-STM shown
## as a (mesh-free) reference point. Gaussian family. Results saved incrementally.
suppressMessages({
  library(FNN); library(fields); library(withr); library(sdmTMB); library(scoringRules); library(parallel)
})
source("engine/internal_utils_dglm.R"); source("engine/cf_dglm_hv.R"); source("engine/cf_dglm.R")
.dglm_load_cpp()

NSEED  <- as.integer(Sys.getenv("NSEED","5"))
NCORE  <- as.integer(Sys.getenv("NCORE","5"))
RHO    <- as.numeric(Sys.getenv("RHO","0.7"))
FAM    <- "gaussian"
## Cutoffs reported in Table A1.  A finer cutoff of 0.01 is NOT run by default:
## a single fit was projected to exceed ten hours (see Appendix 1).  Override with
## e.g.  CUTS="0.2,0.1,0.05,0.025,0.01"
CUTS   <- as.numeric(strsplit(Sys.getenv("CUTS","0.20,0.10,0.05,0.025"),",")[[1]])
OUT    <- Sys.getenv("OUT","results/tableA1_spde_mesh.csv")

## ---- DGP / split / score (from table2_predictive_moving.R) ------------------
gen <- function(n_pop=400, T=40, rho=0.7, range=0.10, fam="gaussian", seed=1){
  set.seed(seed); na <- n_pop; anc <- cbind(runif(na), runif(na))
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
  al <- al * (sd(as.vector(b)) / sd(al))
  y <- 2*X1 - 1.5*X2 + al + rnorm(N, sd=0.5)
  list(y=y, coords=coords, pt=pt, x=cbind(x1=X1,x2=X2), n_pop=n_pop, T=T)
}
split <- function(d, seed=1){
  set.seed(seed+10000); N <- length(d$y); ix <- logical(N); ix[sample.int(N, round(N*0.30))] <- TRUE
  list(ix=ix, d=d,
       y_tr=d$y[!ix], x_tr=d$x[!ix,,drop=FALSE], c_tr=d$coords[!ix,], p_tr=d$pt[!ix],
       y_te=d$y[ix],  x_te=d$x[ix,,drop=FALSE],  c_te=d$coords[ix,],  p_te=d$pt[ix])
}
score <- function(s, mu, tr_res_sd){
  c(rmse=sqrt(mean((s$y_te - mu)^2)),
    crps=mean(scoringRules::crps_norm(s$y_te, mu, pmax(tr_res_sd,1e-6))))
}
trsd <- function(s, mu_tr) sqrt(mean((s$y_tr - mu_tr)^2))

m_sdmtmb <- function(s, cutoff){
  tr <- data.frame(y=s$y_tr, x1=s$x_tr[,1], x2=s$x_tr[,2], cx=s$c_tr[,1], cy=s$c_tr[,2], pt=s$p_tr)
  te <- data.frame(x1=s$x_te[,1], x2=s$x_te[,2], cx=s$c_te[,1], cy=s$c_te[,2], pt=s$p_te)
  mesh <- make_mesh(tr, xy_cols=c("cx","cy"), cutoff=cutoff)
  fit <- sdmTMB(y ~ x1 + x2, data=tr, mesh=mesh, time="pt", spatiotemporal="ar1",
                spatial="on", family=gaussian(), silent=TRUE)
  list(mu=predict(fit, newdata=te)$est, mu_tr=predict(fit)$est, nodes=mesh$mesh$n)
}
m_cfdglm <- function(s, seed){
  mh <- cf_dglm_hv(y=s$y_tr, x=s$x_tr, coords=s$c_tr, time=s$p_tr, family=gaussian(), seed=seed)
  md <- cf_dglm(y=s$y_tr, x=s$x_tr, coords=s$c_tr, time=s$p_tr,
                x0=s$x_te, coords0=s$c_te, time0=s$p_te, mod_hv=mh)
  list(mu=md$pred0$pred, mu_tr=md$pred$pred)
}

seeds <- 1:NSEED
all_rows <- list()
write_out <- function(){ res <- do.call(rbind, all_rows); write.csv(res, OUT, row.names=FALSE); res }

## ---- CF-STM reference (once per seed) ---------------------------------------
cat(sprintf("CF-STM reference: %d seeds\n", NSEED))
cf <- mclapply(seeds, function(sd){
  d <- gen(fam=FAM, seed=sd, rho=RHO); s <- split(d, sd)
  t0 <- proc.time()["elapsed"]
  r  <- tryCatch(m_cfdglm(s, sd), error=function(e) NULL)
  el <- as.numeric(proc.time()["elapsed"]-t0)
  if(is.null(r)) return(data.frame(method="CF-STM",cutoff=NA,nodes=NA,seed=sd,rmse=NA,crps=NA,time=el))
  sc <- score(s, r$mu, trsd(s, r$mu_tr))
  data.frame(method="CF-STM",cutoff=NA,nodes=NA,seed=sd,rmse=sc["rmse"],crps=sc["crps"],time=el)
}, mc.cores=NCORE, mc.preschedule=FALSE)
all_rows <- c(all_rows, cf); write_out()

## ---- sdmTMB per cutoff (coarse -> fine), incremental save -------------------
for(cut in CUTS){
  cat(sprintf("sdmTMB cutoff=%.3f (%d seeds)...\n", cut, NSEED)); flush.console()
  rows <- mclapply(seeds, function(sd){
    d <- gen(fam=FAM, seed=sd, rho=RHO); s <- split(d, sd)
    t0 <- proc.time()["elapsed"]
    r  <- tryCatch(suppressMessages(suppressWarnings(m_sdmtmb(s, cut))), error=function(e) NULL)
    el <- as.numeric(proc.time()["elapsed"]-t0)
    if(is.null(r) || any(!is.finite(r$mu)))
      return(data.frame(method="sdmTMB",cutoff=cut,nodes=if(is.null(r))NA else r$nodes,seed=sd,rmse=NA,crps=NA,time=el))
    sc <- score(s, r$mu, trsd(s, r$mu_tr))
    data.frame(method="sdmTMB",cutoff=cut,nodes=r$nodes,seed=sd,rmse=sc["rmse"],crps=sc["crps"],time=el)
  }, mc.cores=NCORE, mc.preschedule=FALSE)
  all_rows <- c(all_rows, rows); write_out()
  a <- do.call(rbind, rows)
  cat(sprintf("  nodes=%d  RMSE=%.3f  CRPS=%.3f  time=%.1fs\n",
      round(mean(a$nodes,na.rm=TRUE)), mean(a$rmse,na.rm=TRUE), mean(a$crps,na.rm=TRUE), mean(a$time,na.rm=TRUE)))
}

res <- write_out()
cat("\n=== mean over seeds ===\n")
agg <- aggregate(cbind(nodes,rmse,crps,time)~method+cutoff, data=res, FUN=function(z)mean(z,na.rm=TRUE), na.action=na.pass)
print(agg[order(agg$cutoff),], digits=4, row.names=FALSE)
cat(sprintf("\nWrote %s\n", OUT))
