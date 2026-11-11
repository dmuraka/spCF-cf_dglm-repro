## Appendix sensitivity (C): scale-specific ML (rho_r, v_r^2) alternative vs shared. Gaussian, 100 seeds.
suppressMessages({library(FNN);library(fields);library(dbscan);library(withr);library(Matrix);library(parallel)})
EM <- Sys.getenv("ENGINE_MOD"); source(file.path(EM,"internal_utils_dglm.R")); source(file.path(EM,"cf_dglm_hv.R")); source(file.path(EM,"cf_dglm.R")); .dglm_load_cpp()
source(Sys.getenv("GENS"))
NSEED <- as.integer(Sys.getenv("NSEED","100")); NCORE <- as.integer(Sys.getenv("NCORE","12"))
one <- function(design, rho_dgp, seed){
  d <- if(design=="panel") gen_panel(seed=seed,rho=rho_dgp) else gen_moving(seed=seed,rho=rho_dgp); ix<-d$ix
  a<-list(y=d$y[!ix],x=d$x[!ix,],coords=d$coords[!ix,],time=d$pt[!ix],family=gaussian(),seed=seed)
  rows<-list()
  for(kind in c("scalewise_ml")){
    invisible(capture.output(mh<-do.call(cf_dglm_hv,c(a,list(ar1="scalewise",disp_est=FALSE)))))
    md<-cf_dglm(y=a$y,x=a$x,coords=a$coords,time=a$time,x0=d$x[ix,],coords0=d$coords[ix,],time0=d$pt[ix],mod_hv=mh)
    mu<-md$pred0$pred
    rows[[kind]]<-data.frame(design=design,rho_dgp=rho_dgp,seed=seed,kind=kind,R=length(md$bands),
      rmse_y=sqrt(mean((d$y[ix]-mu)^2)),rmse_eta=sqrt(mean((d$eta[ix]-mu)^2)))}
  do.call(rbind,rows)}
tasks<-expand.grid(design=c("panel","moving"),rho_dgp=c(0.2,0.7),seed=1:NSEED,stringsAsFactors=FALSE)
res<-do.call(rbind,mclapply(seq_len(nrow(tasks)),function(i) tryCatch(one(tasks$design[i],tasks$rho_dgp[i],tasks$seed[i]),error=function(e) NULL),mc.cores=NCORE,mc.preschedule=FALSE))
OUT<-Sys.getenv("OUTDIR","results"); dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
write.csv(res,file.path(OUT,"tableA5_scalewise_ml.csv"),row.names=FALSE)
print(aggregate(cbind(rmse_y,rmse_eta,R)~design+rho_dgp+kind,res,mean),digits=4,row.names=FALSE)
