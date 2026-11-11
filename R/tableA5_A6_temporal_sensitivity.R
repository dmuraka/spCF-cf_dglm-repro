## Appendix sensitivity: (A) shared-Q multiplier sweep; (B) scale-wise rho_r ML with 95% CI. Gaussian, 100 seeds.
suppressMessages({library(FNN);library(fields);library(dbscan);library(withr);library(Matrix);library(parallel)})
EM <- Sys.getenv("ENGINE_MOD"); source(file.path(EM,"internal_utils_dglm.R")); source(file.path(EM,"cf_dglm_hv.R")); source(file.path(EM,"cf_dglm.R")); .dglm_load_cpp()
source(Sys.getenv("GENS"))
NSEED <- as.integer(Sys.getenv("NSEED","100")); NCORE <- as.integer(Sys.getenv("NCORE","12")); FACS <- c(0.1,0.3,1,3,10,1e6)
one <- function(design, rho_dgp, seed){
  d <- if(design=="panel") gen_panel(seed=seed,rho=rho_dgp) else gen_moving(seed=seed,rho=rho_dgp); ix<-d$ix
  a<-list(y=d$y[!ix],x=d$x[!ix,],coords=d$coords[!ix,],time=d$pt[!ix],family=gaussian(),seed=seed)
  ## base HV (shared) -> rho1, Q1 ; also scale-wise rho_r diagnostics (rho_only mode, c=1)
  invisible(capture.output(mh0<-do.call(cf_dglm_hv,c(a,list(ar1="rho_only",disp_est=FALSE)))))
  tb<-mh0$other$ar1_table; tb$design<-design; tb$rho_dgp<-rho_dgp; tb$seed<-seed
  rows<-list()
  for(f in FACS){
    invisible(capture.output(mh<-do.call(cf_dglm_hv,c(a,list(rho=mh0$other$rho,Q=mh0$other$Q*f)))))
    md<-cf_dglm(y=a$y,x=a$x,coords=a$coords,time=a$time,x0=d$x[ix,],coords0=d$coords[ix,],time0=d$pt[ix],mod_hv=mh)
    mu<-md$pred0$pred
    rows[[length(rows)+1]]<-data.frame(design=design,rho_dgp=rho_dgp,seed=seed,factor=f,R=length(md$bands),
      rmse_y=sqrt(mean((d$y[ix]-mu)^2)),rmse_eta=sqrt(mean((d$eta[ix]-mu)^2)),rho1=mh0$other$rho,Q1=mh0$other$Q)}
  list(sweep=do.call(rbind,rows),tbl=tb)}
tasks<-expand.grid(design=c("panel","moving"),rho_dgp=c(0.2,0.7),seed=1:NSEED,stringsAsFactors=FALSE)
res<-mclapply(seq_len(nrow(tasks)),function(i) tryCatch(one(tasks$design[i],tasks$rho_dgp[i],tasks$seed[i]),error=function(e) NULL),mc.cores=NCORE,mc.preschedule=FALSE)
res<-res[!sapply(res,is.null)]
sweep<-do.call(rbind,lapply(res,`[[`,"sweep")); tbl<-do.call(rbind,lapply(res,`[[`,"tbl"))
OUT<-Sys.getenv("OUTDIR","results"); dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
write.csv(sweep,file.path(OUT,"tableA5_Q_sweep.csv"),row.names=FALSE)
write.csv(tbl,file.path(OUT,"tableA6_rho_scalewise.csv"),row.names=FALSE)
cat("\n== (A) shared-Q multiplier: mean RMSE(y) / RMSE(eta) / R ==\n")
print(aggregate(cbind(rmse_y,rmse_eta,R)~design+rho_dgp+factor,sweep,mean),digits=4,row.names=FALSE)
cat("\n== (B) scale-wise rho_r (ML, c=1): mean estimate, mean se_z, share identifiable (se_z<0.5) ==\n")
tbl$ident<-is.finite(tbl$se_z)&tbl$se_z<0.5
ag<-aggregate(cbind(band,rho_ml,se_z,ident)~design+rho_dgp+scale,tbl[tbl$ident,],mean)
cnt<-aggregate(ident~design+rho_dgp+scale,tbl,mean); names(cnt)[4]<-"frac_ident"; ag<-merge(ag[,c("design","rho_dgp","scale","band","rho_ml","se_z")],cnt)
print(ag[order(ag$design,ag$rho_dgp,ag$scale),],digits=3,row.names=FALSE)
