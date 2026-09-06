## PAPER TABLE: Table 6 (CF-STM computation-time scaling)  (see ../README.md)
###############################################################################
## cf_dglm computation-time scaling study (latest single-pass cf_dglm / cf_dglm_hv
## ONLY -- no other methods). PANEL data (fixed sites): N sites x T times.
## Grid: N in {400,800,2000,5000} x T in {40,80,200,500} = 16 cases, 5 reps each.
## Gaussian DGP, rho=0.7, range=0.10, 2 covariates (no time-varying coefficients).
##
## Timing is split THREE ways per (N,T,rep):
##   t_hv   : cf_dglm_hv (scale selection)
##   t_fit  : cf_dglm WITHOUT prediction (full-sample model training)
##   t_pred : cf_dglm WITH a fresh prediction grid, MINUS t_fit
##            (the incremental cost of predicting at new locations)
## Prediction set = N fresh fixed locations x T times (same size/structure as the
## training panel). Runs SEQUENTIALLY (one fit at a time) for clean timings; rows
## are appended as completed.
## Output CSV: N, T, n_obs, n_pred, rep, seed, n_scales, t_hv, t_fit, t_pred, t_total.
###############################################################################
suppressMessages({ library(FNN); library(fields); library(dbscan); library(withr) })
source("engine/internal_utils_dglm.R"); source("engine/cf_dglm_hv.R"); source("engine/cf_dglm.R")
.dglm_load_cpp()                                   # compile fused C++ kernel once

OUT  <- Sys.getenv("OUT", "results/table6_scaling.csv")
NREP <- as.integer(Sys.getenv("NREP", "5"))
RHO  <- as.numeric(Sys.getenv("RHO", "0.7"))
RANGE<- 0.10

## fixed-location panel DGP; returns training panel + a fresh prediction panel
gen_panel <- function(N, T, rho=0.7, range=0.10, seed=1){
  set.seed(seed)
  mkfield <- function(S){
    D <- fields::rdist(S); W <- exp(-D/range); W <- W/sqrt(rowSums(W^2))
    b <- matrix(0,T,nrow(S)); b[1,] <- as.numeric(W%*%rnorm(nrow(S)))/sqrt(1-rho^2)
    for(t in 2:T) b[t,] <- rho*b[t-1,] + as.numeric(W%*%rnorm(nrow(S)))
    x1 <- matrix(0,T,nrow(S)); x2 <- matrix(0,T,nrow(S))
    for(t in 1:T){ x1[t,]<-as.numeric(W%*%rnorm(nrow(S))); x2[t,]<-as.numeric(W%*%rnorm(nrow(S))) }
    pt <- rep(1:T,each=nrow(S)); loc <- rep(1:nrow(S),T)
    list(coords=S[loc,], time=pt, x=cbind(x1=x1[cbind(pt,loc)], x2=x2[cbind(pt,loc)]),
         al=b[cbind(pt,loc)])
  }
  Str <- cbind(runif(N), runif(N)); Spr <- cbind(runif(N), runif(N))   # train + new sites
  tr <- mkfield(Str); pr <- mkfield(Spr)
  tr$y <- 2*tr$x[,1] - 1.5*tr$x[,2] + tr$al + rnorm(N*T, sd=0.5)
  list(tr=tr, pr=pr)
}

grid <- expand.grid(N=c(400,800,2000,5000), T=c(40,80,200,500))
grid <- grid[order(grid$N*grid$T), ]                # smallest first -> early feedback

hdr <- !file.exists(OUT)
cat(sprintf("scaling cf_dglm split (panel): %d cases x %d reps -> %s\n", nrow(grid), NREP, OUT))

for(g in seq_len(nrow(grid))){
  N <- grid$N[g]; T <- grid$T[g]
  for(rep in seq_len(NREP)){
    seed <- rep
    d <- gen_panel(N, T, rho=RHO, range=RANGE, seed=seed)
    tr <- d$tr; pr <- d$pr
    ## (1) scale selection
    t0 <- proc.time()["elapsed"]
    invisible(capture.output(suppressWarnings(suppressMessages(
      mh <- cf_dglm_hv(y=tr$y, x=tr$x, coords=tr$coords, time=tr$time, family=gaussian(), seed=seed)))))
    t_hv <- as.numeric(proc.time()["elapsed"] - t0)
    ## (2) model training only (no prediction)
    t0 <- proc.time()["elapsed"]
    invisible(capture.output(suppressWarnings(suppressMessages(
      md0 <- cf_dglm(y=tr$y, x=tr$x, coords=tr$coords, time=tr$time, mod_hv=mh)))))
    t_fit <- as.numeric(proc.time()["elapsed"] - t0)
    ## (3) training + prediction at new locations; prediction cost = difference
    t0 <- proc.time()["elapsed"]
    invisible(capture.output(suppressWarnings(suppressMessages(
      md1 <- cf_dglm(y=tr$y, x=tr$x, coords=tr$coords, time=tr$time,
                     x0=pr$x, coords0=pr$coords, time0=pr$time, mod_hv=mh)))))
    t_fitpred <- as.numeric(proc.time()["elapsed"] - t0)
    t_pred <- max(t_fitpred - t_fit, 0)
    row <- data.frame(N=N, T=T, n_obs=N*T, n_pred=N*T, rep=rep, seed=seed,
                      n_scales=length(mh$other$bands),
                      t_hv=round(t_hv,3), t_fit=round(t_fit,3), t_pred=round(t_pred,3),
                      t_total=round(t_hv + t_fitpred,3))
    write.table(row, OUT, sep=",", append=!hdr, col.names=hdr, row.names=FALSE); hdr <- FALSE
    cat(sprintf("  N=%4d T=%3d (n=%7d) rep %d/%d : hv=%6.2f fit=%6.2f pred=%6.2f total=%6.2f s\n",
                N, T, N*T, rep, NREP, t_hv, t_fit, t_pred, t_hv + t_fitpred))
    rm(d, tr, pr, mh, md0, md1); gc(verbose=FALSE)
  }
}
cat("scaling split done.\n")
