## PAPER TABLE: Tables 3 and 4 (coefficient bias, SD, SE, coverage)  (see ../README.md)
###############################################################################
## Standard-error (SE) estimation accuracy for the regression slopes.
## Methods: GLM, GAM, cf_dglm (CF-STM, robust SE), sdmTMB. KFAS excluded
## (its beta uses a plain pooled GLM, as in Table 3 which drops GSSM).
## For each replicate we record bhat and its reported SE for x1,x2.
## Aggregated per (family): Bias, empirical SD (= true sampling variability of
## bhat), mean reported SE, and 95% CI coverage. A method's SE is "accurate"
## when mean(SE) ~ empirical SD and coverage ~ 0.95.
## True slopes: gaussian (2,-1.5), poisson (0.3,-0.2), binary (0.5,-0.4) [link].
## Output CSV: seed, family, method, b1_hat, se1, b2_hat, se2.
###############################################################################
suppressMessages({
  library(FNN); library(fields); library(dbscan); library(withr)
  library(sdmTMB); library(mgcv); library(parallel)
})
source("engine/internal_utils_dglm.R"); source("engine/cf_dglm_hv.R"); source("engine/cf_dglm.R")
.dglm_load_cpp()

NSEED <- as.integer(Sys.getenv("NSEED","100"))
NCORE <- as.integer(Sys.getenv("NCORE","14"))
OUT   <- Sys.getenv("OUT","results/table3_4_coef_se_rho07.csv")
RHO   <- as.numeric(Sys.getenv("RHO","0.7"))

## same DGP as bench_coef
gen <- function(n_pop=400,T=40,rho=0.7,range=0.10,fam="gaussian",seed=1){
  set.seed(seed)
  s <- cbind(runif(n_pop),runif(n_pop)); D<-fields::rdist(s); W<-exp(-D/range); W<-W/sqrt(rowSums(W^2))
  b <- matrix(0,T,n_pop); b[1,]<-as.numeric(W%*%rnorm(n_pop))/sqrt(1-rho^2)
  for(t in 2:T) b[t,]<-rho*b[t-1,]+as.numeric(W%*%rnorm(n_pop))
  x1<-matrix(0,T,n_pop); x2<-matrix(0,T,n_pop)
  for(t in 1:T){ x1[t,]<-as.numeric(W%*%rnorm(n_pop)); x2[t,]<-as.numeric(W%*%rnorm(n_pop)) }
  pt<-rep(1:T,each=n_pop); loc<-rep(1:n_pop,T); coords<-s[loc,]
  al<-b[cbind(pt,loc)]; X1<-x1[cbind(pt,loc)]; X2<-x2[cbind(pt,loc)]
  if(fam=="gaussian")      y <- 2*X1 - 1.5*X2 + al + rnorm(n_pop*T,sd=0.5)
  else if(fam=="poisson")  y <- rpois(n_pop*T, exp(0.5 + 0.3*X1 - 0.2*X2 + 0.7*al/sd(al)))
  else                     y <- rbinom(n_pop*T,1, plogis(0.0 + 0.5*X1 - 0.4*X2 + 1.2*al/sd(al)))
  list(y=y, coords=coords, pt=pt, loc=loc, x=cbind(x1=X1,x2=X2))
}
fam_obj <- function(f) if(f=="gaussian") gaussian() else if(f=="poisson") poisson() else binomial()
TRUE_B  <- list(gaussian=c(2,-1.5), poisson=c(0.3,-0.2), binary=c(0.5,-0.4))

## each extractor returns c(b1,se1,b2,se2)
es_glm <- function(d,fam){ g<-glm(y~x1+x2,data=data.frame(y=d$y,d$x),family=fam_obj(fam))
  s<-summary(g)$coefficients[c("x1","x2"),1:2]; c(s[1,1],s[1,2],s[2,1],s[2,2]) }
es_gam <- function(d,fam){
  dd<-data.frame(y=d$y,x1=d$x[,1],x2=d$x[,2],cx=d$coords[,1],cy=d$coords[,2],tt=d$pt)
  g<-mgcv::gam(y~x1+x2+te(cx,cy,tt,d=c(2,1),k=c(40,8),bs=c("tp","cr")),data=dd,family=fam_obj(fam),method="REML")
  s<-summary(g); c(s$p.coeff["x1"],s$se["x1"],s$p.coeff["x2"],s$se["x2"]) }
es_cfdglm <- function(d,fam,seed){
  mh<-cf_dglm_hv(y=d$y,x=d$x,coords=d$coords,time=d$pt,family=fam_obj(fam),seed=seed)
  md<-cf_dglm(y=d$y,x=d$x,coords=d$coords,time=d$pt,mod_hv=mh)
  b<-md$beta[c("x1","x2"),c("coef","coef_se")]; c(b[1,1],b[1,2],b[2,1],b[2,2]) }
es_sdmtmb <- function(d,fam){
  tr<-data.frame(y=d$y,x1=d$x[,1],x2=d$x[,2],cx=d$coords[,1],cy=d$coords[,2],pt=d$pt)
  mesh<-make_mesh(tr,xy_cols=c("cx","cy"),cutoff=0.05)
  fit<-sdmTMB(y~x1+x2,data=tr,mesh=mesh,time="pt",spatiotemporal="ar1",spatial="on",family=fam_obj(fam),silent=TRUE)
  td<-sdmTMB::tidy(fit,effects="fixed"); i<-match(c("x1","x2"),td$term)
  c(td$estimate[i[1]],td$std.error[i[1]],td$estimate[i[2]],td$std.error[i[2]]) }

run_task <- function(seed, fam){
  d <- gen(fam=fam, seed=seed, rho=RHO)
  ex <- list(GLM=function() es_glm(d,fam), GAM=function() es_gam(d,fam),
             cf_dglm=function() es_cfdglm(d,fam,seed), sdmTMB=function() es_sdmtmb(d,fam))
  rows<-list()
  for(nm in names(ex)){
    v <- tryCatch(suppressMessages(suppressWarnings({utils::capture.output(r<-ex[[nm]]()); as.numeric(r)})),
                  error=function(e) rep(NA,4))
    rows[[nm]] <- data.frame(seed=seed,family=fam,method=nm,b1_hat=v[1],se1=v[2],b2_hat=v[3],se2=v[4])
  }
  do.call(rbind,rows)
}

tasks <- expand.grid(seed=1:NSEED, fam=c("gaussian","poisson","binary"), stringsAsFactors=FALSE)
cat(sprintf("coef-SE bench: %d tasks on %d cores, rho=%.1f -> %s\n", nrow(tasks), NCORE, RHO, OUT))
rl <- mclapply(seq_len(nrow(tasks)), function(i) tryCatch(run_task(tasks$seed[i],tasks$fam[i]),error=function(e)NULL),
               mc.cores=NCORE, mc.preschedule=FALSE)
res <- do.call(rbind, rl[!sapply(rl,is.null)]); rownames(res)<-NULL
write.csv(res, OUT, row.names=FALSE)
cat(sprintf("wrote %d rows -> %s\n", nrow(res), OUT))

## ---- aggregate: bias, empirical SD, mean reported SE, 95% coverage (x1 slope) ----
## Paper Tables 3-4 report the primary slope x1. The two slopes attenuate in
## OPPOSITE directions under spatial confounding (x1>0 downward, x2<0 upward), so
## pooling would cancel the bias; hence the x1 coefficient is reported here.
ord <- c("GLM","GAM","sdmTMB","cf_dglm")
cat("\n=== SE estimation accuracy (rho=", RHO, ", ", NSEED, " seeds; x1 slope) ===\n", sep="")
cat(sprintf("%-9s %-8s %8s %8s %8s %8s\n","family","method","Bias","EmpSD","MeanSE","Cov95"))
for(fm in c("gaussian","poisson","binary")){
  b1<-TRUE_B[[fm]][1]
  for(nm in ord){
    sub<-res[res$family==fm & res$method==nm,]
    bh<-sub$b1_hat; se<-sub$se1
    ok<-is.finite(bh)&is.finite(se)
    bias <- mean(bh[ok]-b1)
    empsd<- sd(bh[ok])
    meanse<-mean(se[ok]); cov<-mean(abs(bh[ok]-b1)<=1.96*se[ok])
    cat(sprintf("%-9s %-8s %8.4f %8.4f %8.4f %8.3f\n",fm,nm,bias,empsd,meanse,cov))
  }
}
