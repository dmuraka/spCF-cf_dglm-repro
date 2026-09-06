## PAPER TABLE: Table 5 (computation-time comparison, all methods)  (see ../README.md)
###############################################################################
## Computation-time benchmark, ALL methods, on PANEL data (fixed sites).
## Grid: T in {40,80,200} x N in {400,800} = 6 cells, 5 reps each.
## Methods: GLM, cf_dglm, KFAS_K120, GAM, sdmTMB, KFAS_K240 (cheap->expensive).
## FIT-only (no holdout / no prediction), gaussian, rho=0.7, range=0.10.
## SEQUENTIAL: one fit at a time so timings are not polluted by contention.
## Methods ordered cheap->expensive so the 5 lighter methods finish first and
## KFAS_K240 (the heaviest) fills in last. Rows appended as completed.
## Output CSV: method, N, T, n_obs, rep, seed, time_fit.
###############################################################################
suppressMessages({ library(FNN); library(fields); library(dbscan); library(withr)
  library(KFAS); library(sdmTMB); library(mgcv) })
source("engine/internal_utils_dglm.R"); source("engine/cf_dglm_hv.R"); source("engine/cf_dglm.R"); .dglm_load_cpp()

OUT  <- Sys.getenv("OUT", "results/table5_timing.csv")
NREP <- as.integer(Sys.getenv("NREP", "5"))
RHO  <- 0.7; RANGE <- 0.10

gen_panel <- function(N, T, rho=0.7, range=0.10, seed=1){
  set.seed(seed)
  s <- cbind(runif(N), runif(N)); D <- fields::rdist(s); W <- exp(-D/range); W <- W/sqrt(rowSums(W^2))
  b <- matrix(0,T,N); b[1,] <- as.numeric(W%*%rnorm(N))/sqrt(1-rho^2)
  for(t in 2:T) b[t,] <- rho*b[t-1,] + as.numeric(W%*%rnorm(N))
  x1 <- matrix(0,T,N); x2 <- matrix(0,T,N)
  for(t in 1:T){ x1[t,]<-as.numeric(W%*%rnorm(N)); x2[t,]<-as.numeric(W%*%rnorm(N)) }
  pt <- rep(1:T,each=N); loc <- rep(1:N,T); s2 <- s[loc,]
  al <- b[cbind(pt,loc)]; X1 <- x1[cbind(pt,loc)]; X2 <- x2[cbind(pt,loc)]
  y  <- 2*X1 - 1.5*X2 + al + rnorm(N*T, sd=0.5)
  list(y=y, coords=s2, time=pt, x=cbind(x1=X1,x2=X2), s=s, N=N, T=T)
}

## fit-only timers (estimate the model; no prediction) -----------------------
fit_glm <- function(d) glm(y~x1+x2, data=data.frame(y=d$y,d$x), family=gaussian())
fit_cfdglm <- function(d){
  mh <- cf_dglm_hv(y=d$y, x=d$x, coords=d$coords, time=d$time, family=gaussian(), seed=1)
  cf_dglm(y=d$y, x=d$x, coords=d$coords, time=d$time, mod_hv=mh)
}
fit_gam <- function(d){ dd<-data.frame(y=d$y,x1=d$x[,1],x2=d$x[,2],cx=d$coords[,1],cy=d$coords[,2],tt=d$time)
  mgcv::gam(y~x1+x2+te(cx,cy,tt,d=c(2,1),k=c(40,8),bs=c("tp","cr")), data=dd, family=gaussian(), method="REML") }
fit_sdmtmb <- function(d){ tr<-data.frame(y=d$y,x1=d$x[,1],x2=d$x[,2],cx=d$coords[,1],cy=d$coords[,2],pt=d$time)
  mesh<-make_mesh(tr, xy_cols=c("cx","cy"), cutoff=0.05)
  sdmTMB(y~x1+x2, data=tr, mesh=mesh, time="pt", spatiotemporal="ar1", spatial="on", family=gaussian(), silent=TRUE) }
fit_kfas <- function(d, bands, nk, seed=1){
  N<-d$N; T<-d$T; resp<-d$y
  beta<-lm.fit(cbind(1,d$x),resp)$coefficients; r<-resp-drop(cbind(1,d$x)%*%beta)
  Y<-matrix(r, T, N, byrow=TRUE)                      # panel: row=time, col=site
  kn<-lapply(seq_along(bands), function(j) with_seed(seed, kmeans(d$s, nk[j], iter.max=20)$centers))
  Z<-do.call(cbind, lapply(seq_along(bands), function(j) exp(-rdist(d$s,kn[[j]])/bands[j]))); K<-ncol(Z)
  mod<-SSModel(Y ~ -1 + SSMcustom(Z=Z, T=diag(K), R=diag(K), Q=diag(K),
               a1=rep(0,K), P1=diag(K), state_names=paste0("a",1:K)), H=diag(N))
  upd<-function(p,model){ rho<-tanh(p[1]); q<-exp(p[2]); s2<-exp(p[3])
    model$T[,,1]<-rho*diag(K); model$Q[,,1]<-q*diag(K); model$P1[]<-(q/(1-rho^2))*diag(K); model$H[,,1]<-s2*diag(N); model }
  fit<-fitSSM(mod, inits=c(atanh(0.7),log(1),log(0.25)), updatefn=upd, method="BFGS", control=list(maxit=20))
  KFS(fit$model, smoothing="state")
}

## cheap -> expensive so light methods finish first; KFAS_K240 last
methods <- list(
  GLM       = function(d) fit_glm(d),
  cf_dglm   = function(d) fit_cfdglm(d),
  KFAS_K120 = function(d) fit_kfas(d, 0.12, 120),
  GAM       = function(d) fit_gam(d),
  sdmTMB    = function(d) fit_sdmtmb(d),
  KFAS_K240 = function(d) fit_kfas(d, c(0.20,0.10,0.05), c(40,80,120))
)
cells <- expand.grid(N=c(400,800), T=c(40,80,200))
cells <- cells[order(cells$N*cells$T), ]               # smallest first

hdr <- !file.exists(OUT)
cat(sprintf("timing bench (panel): %d methods x %d cells x %d reps -> %s\n",
            length(methods), nrow(cells), NREP, OUT))
for(nm in names(methods)){
  for(ci in seq_len(nrow(cells))){
    N<-cells$N[ci]; T<-cells$T[ci]
    for(rep in seq_len(NREP)){
      d <- gen_panel(N, T, rho=RHO, range=RANGE, seed=rep)
      t0 <- proc.time()["elapsed"]
      ok <- tryCatch({ invisible(capture.output(suppressWarnings(suppressMessages(methods[[nm]](d))))); TRUE },
                     error=function(e) conditionMessage(e))
      el <- as.numeric(proc.time()["elapsed"]-t0)
      row <- data.frame(method=nm, N=N, T=T, n_obs=N*T, rep=rep, seed=rep,
                        time_fit=if(isTRUE(ok)) round(el,3) else NA,
                        note=if(isTRUE(ok)) "" else substr(ok,1,40))
      write.table(row, OUT, sep=",", append=!hdr, col.names=hdr, row.names=FALSE); hdr<-FALSE
      cat(sprintf("  %-10s N=%4d T=%3d (n=%6d) rep %d/%d : %s\n", nm, N, T, N*T, rep, NREP,
                  if(isTRUE(ok)) sprintf("%.1fs",el) else paste("ERR",substr(ok,1,30))))
      rm(d); gc(verbose=FALSE)
    }
  }
}
cat("timing bench done.\n")
