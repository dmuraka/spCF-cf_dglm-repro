## PAPER TABLE: Table A3 (Appendix 2) -- coefficient estimation accuracy under the
## Gaussian-process DGP  (see ../README.md)
###############################################################################
## Regression-coefficient (slope) estimation accuracy, per replicate.
## DGP: Matern GP (smoothness 1, marginal variance 1) innovations evolving as a
## temporal AR(1); the covariates come from the same GP.  All other settings are
## those of Section 4.
## Methods: GLM, GAM, cf_dglm, sdmTMB, KFAS (construction-A pooled beta).
## True slopes: gaussian (2,-1.5), poisson (0.3,-0.2), binary (0.5,-0.4) [link].
## beta_rmse = sqrt(mean((bhat_j - btrue_j)^2)) over j in {x1,x2}.
## Output CSV: seed, family, method, beta_rmse, b1_hat, b2_hat.
## (KFAS beta = pooled regression ignoring the field; on log/logit-transformed
##  response for non-Gaussian -> approximate, flagged.)
###############################################################################
suppressMessages({
  library(FNN); library(fields); library(dbscan); library(withr)
  library(sdmTMB); library(mgcv); library(parallel)
})
source("engine/internal_utils_dglm.R"); source("engine/cf_dglm_hv.R"); source("engine/cf_dglm.R")
.dglm_load_cpp()  # precompile fused C++ kernel in parent before mclapply forks (avoids per-worker build race)

NSEED <- as.integer(Sys.getenv("NSEED","100"))
NCORE <- as.integer(Sys.getenv("NCORE","20"))
OUT   <- Sys.getenv("OUT","results/tableA3_gp_coef_rho07.csv")
RHO   <- as.numeric(Sys.getenv("RHO","0.7"))

## same DGP as bench_methods (gaussian/poisson) + binary
gen <- function(n_pop=400,T=40,rho=0.7,range=0.10,fam="gaussian",seed=1){
  set.seed(seed)
  s <- cbind(runif(n_pop),runif(n_pop)); D<-fields::rdist(s)
  Sig <- fields::Matern(D, range=range, smoothness=1)      # Matern GP (nu=1), var 1
  Rc  <- chol(Sig + diag(1e-8, n_pop)); gp <- function() as.numeric(crossprod(Rc, rnorm(n_pop)))
  b <- matrix(0,T,n_pop); b[1,]<-gp()/sqrt(1-rho^2)
  for(t in 2:T) b[t,]<-rho*b[t-1,]+gp()
  x1<-matrix(0,T,n_pop); x2<-matrix(0,T,n_pop)
  for(t in 1:T){ x1[t,]<-gp(); x2[t,]<-gp() }
  pt<-rep(1:T,each=n_pop); loc<-rep(1:n_pop,T); coords<-s[loc,]
  al<-b[cbind(pt,loc)]; X1<-x1[cbind(pt,loc)]; X2<-x2[cbind(pt,loc)]
  if(fam=="gaussian")      y <- 2*X1 - 1.5*X2 + al + rnorm(n_pop*T,sd=0.5)
  else if(fam=="poisson")  y <- rpois(n_pop*T, exp(0.5 + 0.3*X1 - 0.2*X2 + 0.7*al/sd(al)))
  else                     y <- rbinom(n_pop*T,1, plogis(0.0 + 0.5*X1 - 0.4*X2 + 1.2*al/sd(al)))
  list(y=y, coords=coords, pt=pt, loc=loc, x=cbind(x1=X1,x2=X2))
}
fam_obj <- function(f) if(f=="gaussian") gaussian() else if(f=="poisson") poisson() else binomial()
TRUE_B  <- list(gaussian=c(2,-1.5), poisson=c(0.3,-0.2), binary=c(0.5,-0.4))

## beta extractors (return c(b1,b2) for x1,x2)
b_glm <- function(d,fam){ g<-glm(y~x1+x2,data=data.frame(y=d$y,d$x),family=fam_obj(fam)); coef(g)[c("x1","x2")] }
b_gam <- function(d,fam){
  dd<-data.frame(y=d$y,x1=d$x[,1],x2=d$x[,2],cx=d$coords[,1],cy=d$coords[,2],tt=d$pt)
  g<-mgcv::gam(y~x1+x2+te(cx,cy,tt,d=c(2,1),k=c(40,8),bs=c("tp","cr")),data=dd,family=fam_obj(fam),method="REML")
  coef(g)[c("x1","x2")]
}
b_cfdglm <- function(d,fam,seed){
  mh<-cf_dglm_hv(y=d$y,x=d$x,coords=d$coords,time=d$pt,family=fam_obj(fam),seed=seed)
  md<-cf_dglm(y=d$y,x=d$x,coords=d$coords,time=d$pt,mod_hv=mh)
  bc<-md$beta$coef; names(bc)<-rownames(md$beta); bc[c("x1","x2")]
}
b_sdmtmb <- function(d,fam){
  tr<-data.frame(y=d$y,x1=d$x[,1],x2=d$x[,2],cx=d$coords[,1],cy=d$coords[,2],pt=d$pt)
  mesh<-make_mesh(tr,xy_cols=c("cx","cy"),cutoff=0.05)
  fit<-sdmTMB(y~x1+x2,data=tr,mesh=mesh,time="pt",spatiotemporal="ar1",spatial="on",family=fam_obj(fam),silent=TRUE)
  td<-sdmTMB::tidy(fit,effects="fixed"); setNames(td$estimate[match(c("x1","x2"),td$term)],c("x1","x2"))
}
b_kfas <- function(d,fam){   # pooled regression (pre-field), transformed resp for non-gaussian
  resp <- if(fam=="gaussian") d$y else if(fam=="poisson") log(d$y+0.5) else log((d$y+0.5)/(1.5-d$y))
  bb<-lm.fit(cbind(1,d$x),resp)$coefficients; bb[2:3]
}

run_task <- function(seed, fam){
  d <- gen(fam=fam, seed=seed, rho=RHO); tb <- TRUE_B[[fam]]
  ex <- list(GLM=function() b_glm(d,fam), GAM=function() b_gam(d,fam),
             cf_dglm=function() b_cfdglm(d,fam,seed), sdmTMB=function() b_sdmtmb(d,fam),
             KFAS=function() b_kfas(d,fam))
  rows<-list()
  for(nm in names(ex)){
    bh <- tryCatch(suppressMessages(suppressWarnings({utils::capture.output(v<-ex[[nm]]()); as.numeric(v)})),
                   error=function(e) c(NA,NA))
    rmse <- if(any(is.na(bh))) NA else sqrt(mean((bh-tb)^2))
    rows[[nm]] <- data.frame(seed=seed,family=fam,method=nm,beta_rmse=rmse,b1_hat=bh[1],b2_hat=bh[2])
  }
  do.call(rbind,rows)
}

tasks <- expand.grid(seed=1:NSEED, fam=c("gaussian","poisson","binary"), stringsAsFactors=FALSE)
cat(sprintf("coef bench: %d tasks on %d cores -> %s\n", nrow(tasks), NCORE, OUT))
rl <- mclapply(seq_len(nrow(tasks)), function(i) tryCatch(run_task(tasks$seed[i],tasks$fam[i]),error=function(e)NULL),
               mc.cores=NCORE, mc.preschedule=FALSE)
res <- do.call(rbind, rl[!sapply(rl,is.null)]); rownames(res)<-NULL
write.csv(res, OUT, row.names=FALSE)
cat(sprintf("wrote %d rows -> %s\n", nrow(res), OUT))
cat("\n=== beta_rmse mean over seeds (true slopes per family) ===\n")
agg<-aggregate(beta_rmse~family+method,res,function(z)mean(z,na.rm=TRUE))
print(agg[order(agg$family,agg$beta_rmse),],digits=4,row.names=FALSE)

## Table A3 summary: bias and empirical SD of the x1 slope (paper reports beta1)
cat("\n=== Table A3: bias / SD of the x1 slope (rho=", RHO, ") ===\n", sep="")
TB1 <- sapply(TRUE_B, function(z) z[1])
for(fm in c("gaussian","poisson","binary")){
  for(nm in c("GLM","GAM","sdmTMB","cf_dglm")){
    sub <- res[res$family==fm & res$method==nm, ]
    bh  <- sub$b1_hat[is.finite(sub$b1_hat)]
    if(length(bh)) cat(sprintf("%-9s %-8s bias=%7.3f  SD=%6.3f  n=%d\n",
                               fm, nm, mean(bh)-TB1[[fm]], sd(bh), length(bh)))
  }
}
