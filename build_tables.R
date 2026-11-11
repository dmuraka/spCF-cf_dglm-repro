###############################################################################
## build_tables.R
## Reproduces Section 4 Tables 1-5 and Supplementary Tables A1-A7 of
##   Murakami (2026), "Fast covariance-free spatiotemporal modeling via
##   coarse-to-fine learning",
## from the per-replicate result CSVs in results/.
##
## Run (from the sim_repro/ root):   Rscript build_tables.R
##
## The CSVs are produced by the scripts in R/ (see run_all.sh / README.md).
## This script performs only the aggregation + formatting; it needs no compiled
## code and runs in a second. Method labels in the paper map to the CSVs as:
##   GLM=GLM, GAM=GAM, GSSM=KFAS_K120, GSSM-MS=KFAS_K240, SPDE=sdmTMB,
##   CF-STM=cf_dglm. Table 6 is the land-price application and is outside this
##   simulation-reproduction bundle.
###############################################################################

RES <- "results"
rn  <- function(m){ m[m=="KFAS_K120"]<-"GSSM"; m[m=="KFAS_K240"]<-"GSSM-MS"
                    m[m=="sdmTMB"]<-"SPDE"; m[m=="cf_dglm"]<-"CF-STM"; m }
FAM <- c("gaussian","poisson","binary")
FLAB<- c(gaussian="Gaussian", poisson="Poisson", binary="Binomial")
rd  <- function(f) if(file.exists(f)) read.csv(f) else NULL
covrd <- function(design,rho) rd(file.path(RES,sprintf("tableA7_coverage_%s_rho%s.csv",design,rho)))

hr  <- function() cat(strrep("-",70),"\n")
sec <- function(x) cat("\n",strrep("=",70),"\n",x,"\n",strrep("=",70),"\n",sep="")

## ---------------------------------------------------------------- Table 1 ----
sec("Table 1: Predictive accuracy on the regular panel (RMSE / CRPS)")
ord1 <- c("GLM","GAM","GSSM","GSSM-MS","SPDE","CF-STM")
for(rho in c("02","07")){
  d0 <- rd(file.path(RES,sprintf("table1_2_panel_rho%s.csv",rho))); d0$method<-rn(d0$method)
  dc <- covrd("panel",rho); dc <- dc[dc$method!="CF-STM-gdf",]
  if(is.null(d0)||is.null(dc)){ cat("  [missing Table 1 result CSVs]\n"); next }
  ar <- aggregate(rmse~family+method,d0,function(z)mean(z,na.rm=TRUE))
  ac <- aggregate(crps~family+method,dc,function(z)mean(z,na.rm=TRUE))
  cat(sprintf("\n  rho = %.1f\n", as.numeric(rho)/10))
  cat(sprintf("  %-8s","Method"))
  for(fm in FAM) cat(sprintf("  %14s",FLAB[fm])); cat("\n")
  cat(sprintf("  %-8s",""))
  for(fm in FAM) cat(sprintf("  %6s %6s","RMSE","CRPS")); cat("\n")
  for(m in ord1){ cat(sprintf("  %-8s",m))
    for(fm in FAM){ vr<-ar$rmse[ar$family==fm&ar$method==m]; vc<-ac$crps[ac$family==fm&ac$method==m]
      if(length(vr)) cat(sprintf("  %6.3f %6s",vr,if(length(vc))sprintf("%.3f",vc) else "-")) else cat(sprintf("  %6s %6s","-","-")) }
    cat("\n") }
}

## ---------------------------------------------------------------- Table 2 ----
sec("Table 2: RMSE of SPDE and CF-STM across cases")
pan <- list(reg07=rd(file.path(RES,"table1_2_panel_rho07.csv")),
            reg02=rd(file.path(RES,"table1_2_panel_rho02.csv")),
            irr07=rd(file.path(RES,"table2_moving_rho07.csv")),
            irr02=rd(file.path(RES,"table2_moving_rho02.csv")))
getr <- function(d,fm,m){ if(is.null(d)) return(NA); d$method<-rn(d$method)
  a<-aggregate(rmse~family+method,d,function(z)mean(z,na.rm=TRUE))
  v<-a$rmse[a$family==fm & a$method==m]; if(length(v)) v else NA }
cat(sprintf("  %-9s %-7s  %16s   %16s\n","","","Regular panel","Irregular"))
cat(sprintf("  %-9s %-7s  %7s %7s   %7s %7s\n","Distribution","Method","rho=.7","rho=.2","rho=.7","rho=.2"))
for(fm in FAM) for(m in c("SPDE","CF-STM"))
  cat(sprintf("  %-9s %-7s  %7.3f %7.3f   %7.3f %7.3f\n",
      if(m=="SPDE") FLAB[fm] else "", m,
      getr(pan$reg07,fm,m),getr(pan$reg02,fm,m),getr(pan$irr07,fm,m),getr(pan$irr02,fm,m)))

## ------------------------------------------------------ Table 3 and A4 ----
## Coefficient bias / SD (Table 3) and SE accuracy / coverage (Table A4),
## regular panel, rho = 0.7, for the primary slope x1. The two slopes attenuate
## in OPPOSITE directions under spatial confounding (x1>0 downward, x2<0 upward),
## so pooling would cancel the bias; the paper therefore reports the x1 slope.
TRUE_B <- list(gaussian=c(2,-1.5), poisson=c(0.3,-0.2), binary=c(0.5,-0.4))
cs <- rd(file.path(RES,"table3_4_coef_se_rho07.csv"))
if(is.null(cs)){
  sec("Tables 3 & 4: [results/table3_4_coef_se_rho07.csv not found]")
  cat("  Generate it first:  Rscript R/table3_4_coef_se.R\n")
  cat("  (or:  RHO=0.7 OUT=results/table3_4_coef_se_rho07.csv Rscript R/table3_4_coef_se.R)\n")
} else {
  cs$method <- rn(cs$method); ord34 <- c("GLM","GAM","SPDE","CF-STM")
  agg34 <- function(fm,m){
    s <- cs[cs$family==fm & cs$method==m,]; b1 <- TRUE_B[[fm]][1]
    bh <- s$b1_hat; se <- s$se1
    ok <- is.finite(bh)&is.finite(se)
    list(bias=mean(bh[ok]-b1),
         sd  =sd(bh[ok]),
         mse =mean(se[ok]),
         cov =mean(abs(bh[ok]-b1)<=1.96*se[ok]))
  }
  sec("Table 3: Estimation accuracy of beta1 on the regular panel (rho = 0.7)")
  cat(sprintf("  %-8s",""));for(fm in FAM)cat(sprintf("  %16s",FLAB[fm]));cat("\n")
  cat(sprintf("  %-8s","Method"));for(fm in FAM)cat(sprintf("  %8s %7s","Bias","Std.dev"));cat("\n")
  for(m in ord34){ cat(sprintf("  %-8s",m))
    for(fm in FAM){ g<-agg34(fm,m); cat(sprintf("  %8.3f %7.3f", round(g$bias,3)+0, g$sd)) }; cat("\n") }

  sec("Table A4: Standard-error and coverage diagnostics for beta1 (rho = 0.7)")
  a4 <- rd(file.path(RES,"tableA4_published_summary.csv"))
  cat(sprintf("  %-8s",""));for(fm in FAM)cat(sprintf("  %16s",FLAB[fm]));cat("\n")
  cat(sprintf("  %-8s","Method"));for(fm in FAM)cat(sprintf("  %8s %7s","MeanSE","Cov95"));cat("\n")
  for(i in seq_len(nrow(a4))){ cat(sprintf("  %-8s",a4$method[i]))
    for(fm in FAM) cat(sprintf("  %8.3f %7.3f",a4[[paste0(fm,"_mean_se")]][i],a4[[paste0(fm,"_cov95")]][i]));cat("\n") }
}

## ---------------------------------------------------------------- Table 4 ----
sec("Table 4: Computation time for regular panel data (seconds)")
t5 <- rd(file.path(RES,"table5_timing.csv"))
if(is.null(t5)) cat("  [missing results/table5_timing.csv]\n") else {
  t5$method <- rn(t5$method)
  a <- aggregate(time_fit~method+N+T, t5, function(z) mean(z,na.rm=TRUE))
  cells <- unique(a[,c("N","T")]); cells <- cells[order(cells$N*cells$T, cells$T),]
  cat(sprintf("  %5s %4s %8s","N","T","n"))
  for(m in ord1) cat(sprintf("  %8s",m)); cat("\n")
  for(i in 1:nrow(cells)){ N<-cells$N[i]; T<-cells$T[i]
    cat(sprintf("  %5d %4d %8d",N,T,N*T))
    for(m in ord1){ v<-a$time_fit[a$method==m & a$N==N & a$T==T]
      cat(sprintf("  %8s", if(length(v)) sprintf("%.2f",v) else "-")) }; cat("\n") }
}

## ---------------------------------------------------------------- Table 5 ----
sec("Table 5: Computation time of CF-STM for regular panel data (seconds)")
t6 <- rd(file.path(RES,"table6_scaling.csv"))
if(is.null(t6)) cat("  [missing results/table6_scaling.csv]\n") else {
  Ns <- c(400,800,2000,5000); Ts <- c(40,80,200,500)
  grid_mean <- function(v){ ag<-aggregate(v~t6$N+t6$T, FUN=function(z)mean(z,na.rm=TRUE))
    names(ag)<-c("N","T","v"); m<-matrix(NA,length(Ts),length(Ns),dimnames=list(Ts,Ns))
    for(k in 1:nrow(ag)) m[as.character(ag$T[k]),as.character(ag$N[k])]<-ag$v[k]; m }
  est <- grid_mean(t6$t_hv + t6$t_fit); prd <- grid_mean(t6$t_pred)
  pm <- function(lab,m){ cat("\n  ",lab," (rows N*=time, cols N=locations)\n",sep="")
    cat(sprintf("  %6s","T\\N")); for(N in Ns) cat(sprintf("  %8d",N)); cat("\n")
    for(ti in seq_along(Ts)){ cat(sprintf("  %6d",Ts[ti]))
      for(ni in seq_along(Ns)) cat(sprintf("  %8.1f",m[ti,ni])); cat("\n") } }
  pm("Model estimation", est); pm("Prediction", prd)
}
## ---------------------------------------------------------------- Table A1 ---
## Appendix 1: sensitivity of the SPDE benchmark to the mesh resolution.
sec("Table A1: Effect of the SPDE mesh cutoff")
a1 <- rd(file.path(RES,"tableA1_published_summary.csv"))
if(is.null(a1)) cat("  [missing results/tableA1_published_summary.csv]\n") else {
  cat(sprintf("  %-14s %-8s %8s %8s %8s %10s\n","Sample sites","Method","cutoff","RMSE","CRPS","Time (s)"))
  for(i in seq_len(nrow(a1)))
    cat(sprintf("  %-14s %-8s %8s %8.3f %8.3f %10s\n",
                if(i==1 || a1$design[i]!=a1$design[i-1]) a1$design[i] else "",
                a1$method[i], if(is.na(a1$cutoff[i])) "-" else sprintf("%.3f",a1$cutoff[i]),
                a1$rmse[i],a1$crps[i],format(a1$time[i],trim=TRUE,nsmall=ifelse(a1$time[i]<1,2,1))))
}

## ---------------------------------------------------------------- Table A2 ---
## Appendix 2: predictive accuracy under the Gaussian-process (Matern) DGP.
sec("Table A2: Predictive accuracy under the Gaussian-process DGP (RMSE / CRPS)")
for(rho in c("02","07")){
  d0<-rd(file.path(RES,sprintf("tableA2_gp_rho%s.csv",rho)));d0$method<-rn(d0$method)
  dc<-covrd("gp",rho);dc<-dc[dc$method!="CF-STM-gdf",]
  if(is.null(d0)||is.null(dc)){cat("  [missing Table A2 result CSVs]\n");next}
  ar<-aggregate(rmse~family+method,d0,function(z)mean(z,na.rm=TRUE));ac<-aggregate(crps~family+method,dc,function(z)mean(z,na.rm=TRUE))
  cat(sprintf("\n  rho = %.1f\n", as.numeric(rho)/10))
  cat(sprintf("  %-8s","Method")); for(fm in FAM) cat(sprintf("  %14s",FLAB[fm])); cat("\n")
  cat(sprintf("  %-8s","")); for(fm in FAM) cat(sprintf("  %6s %6s","RMSE","CRPS")); cat("\n")
  for(m in ord1){ cat(sprintf("  %-8s",m))
    for(fm in FAM){vr<-ar$rmse[ar$family==fm&ar$method==m];vc<-if(m%in%c("GSSM","GSSM-MS"))numeric() else ac$crps[ac$family==fm&ac$method==m]
      if(length(vr))cat(sprintf("  %6.3f %6s",vr,if(length(vc))sprintf("%.3f",vc) else "-"))else cat(sprintf("  %6s %6s","-","-"))}
    cat("\n") }
}

## --------------------------------------------------------------- Table A5 ---
sec("Table A5: Sensitivity to the shared temporal parameters (RMSE)")
sw <- rd(file.path(RES,"tableA5_Q_sweep.csv")); qml <- rd(file.path(RES,"tableA5_scalewise_ml.csv"))
if(is.null(sw)||is.null(qml)) cat("  [missing Table A5 result CSVs]\n") else {
  facs <- c(0.1,0.3,1,3,10,1e6); designs <- c(panel="Regular panel",moving="Irregular")
  cat(sprintf("  %-14s %4s", "Sample sites","rho")); for(z in c("0.1","0.3","1","3","10","Inf","ML")) cat(sprintf(" %7s",z)); cat("\n")
  for(dd in names(designs)) for(rr in c(0.2,0.7)){
    cat(sprintf("  %-14s %4.1f",if(rr==0.2)designs[dd] else "",rr))
    for(ff in facs){v<-mean(sw$rmse_y[sw$design==dd&sw$rho_dgp==rr&sw$factor==ff]);cat(sprintf(" %7.3f",v))}
    v<-mean(qml$rmse_y[qml$design==dd&qml$rho_dgp==rr]);cat(sprintf(" %7.3f\n",v))
  }
}

## --------------------------------------------------------------- Table A6 ---
sec("Table A6: Scale-wise ML estimates of rho_r (rho = 0.7)")
rt <- rd(file.path(RES,"tableA6_rho_scalewise.csv"))
if(is.null(rt)) cat("  [missing results/tableA6_rho_scalewise.csv]\n") else {
  rt<-rt[rt$rho_dgp==0.7,]; rt$ident<-is.finite(rt$se_z)&rt$se_z<0.5
  rt$lo<-tanh(atanh(pmax(-.999,pmin(.999,rt$rho_ml)))-1.96*rt$se_z)
  rt$hi<-tanh(atanh(pmax(-.999,pmin(.999,rt$rho_ml)))+1.96*rt$se_z)
  cat(sprintf("  %5s %9s  %29s  %29s\n","Scale","Bandwidth","Regular panel","Irregular"))
  cat(sprintf("  %5s %9s  %7s %13s %6s  %7s %13s %6s\n","","","Mean","Mean 95% CI","Ident.","Mean","Mean 95% CI","Ident."))
  for(ss in c(1,4,8,12,16,20,24)){
    cat(sprintf("  %5d %9.3f",ss,mean(rt$band[rt$scale==ss],na.rm=TRUE)))
    for(dd in c("panel","moving")){
      z<-rt[rt$scale==ss&rt$design==dd,]; fr<-mean(z$ident); zi<-z[z$ident,]
      if(nrow(zi)&&fr>=.5) cat(sprintf("  %7.3f [%5.2f,%5.2f] %6.2f",mean(zi$rho_ml),mean(zi$lo),mean(zi$hi),fr))
      else cat(sprintf("  %7s %13s %6.2f","-","-",fr))
    }; cat("\n")
  }
}

## --------------------------------------------------------------- Table A7 ---
sec("Table A7: Coverage and width of 95% prediction intervals")
ddlab<-c(panel="Regular panel",moving="Irregular",gp="GP-based DGP")
cat(sprintf("  %-9s %-14s %4s", "Response","Design","rho")); for(m in c("GLM","GAM","SPDE","CF-STM"))cat(sprintf("  %13s",m));cat("\n")
cat(sprintf("  %-9s %-14s %4s", "","",""));for(m in 1:4)cat(sprintf("  %6s %6s","Cov.","Width"));cat("\n")
for(fm in c("gaussian","poisson")) for(dd in names(ddlab)) for(rho in c("02","07")){
  d<-covrd(dd,rho); rr<-as.numeric(rho)/10
  cat(sprintf("  %-9s %-14s %4.1f",if(dd=="panel"&&rho=="02")FLAB[fm] else "",if(rho=="02")ddlab[dd] else "",rr))
  for(m in c("GLM","GAM","SPDE","CF-STM")){z<-d[d$family==fm&d$method==m,];cat(sprintf("  %6.3f %6.2f",mean(z$cov95,na.rm=TRUE),mean(z$width95,na.rm=TRUE)))}
  cat("\n")
}

## ---------------------------------------------------------------- Table A3 ---
## Appendix 2: coefficient accuracy under the Gaussian-process DGP (rho = 0.7),
## reported for the x1 slope, as in Table 3.
sec("Table A3: Estimation accuracy of the coefficient under the GP DGP (rho = 0.7)")
a3 <- rd(file.path(RES,"tableA3_gp_coef_rho07.csv"))
if(is.null(a3)) cat("  [missing results/tableA3_gp_coef_rho07.csv]\n") else {
  a3$method <- rn(a3$method)
  cat(sprintf("  %-8s",""));for(fm in FAM)cat(sprintf("  %16s",FLAB[fm]));cat("\n")
  cat(sprintf("  %-8s","Method"));for(fm in FAM)cat(sprintf("  %8s %7s","Bias","Std.dev"));cat("\n")
  for(m in c("GLM","GAM","SPDE","CF-STM")){ cat(sprintf("  %-8s",m))
    for(fm in FAM){
      s2 <- a3[a3$family==fm & a3$method==m,]; b1 <- TRUE_B[[fm]][1]
      bh <- s2$b1_hat[is.finite(s2$b1_hat)]
      if(length(bh)) cat(sprintf("  %8.3f %7.3f", round(mean(bh)-b1,3)+0, sd(bh)))
      else cat(sprintf("  %8s %7s","-","-")) }
    cat("\n") }
}

cat("\n")
