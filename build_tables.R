###############################################################################
## build_tables.R
## Reproduces Section 4 (Monte Carlo experiments) Tables 1-6 of
##   Murakami et al. (2026), "Coarse-to-fine dynamic spatio-temporal models",
## from the per-replicate result CSVs in results/.
##
## Run (from the sim_repro/ root):   Rscript build_tables.R
##
## The CSVs are produced by the scripts in R/ (see run_all.sh / README.md).
## This script performs only the aggregation + formatting; it needs no compiled
## code and runs in a second. Method labels in the paper map to the CSVs as:
##   GLM=GLM, GAM=GAM, GSSM=KFAS_K120, GSSM-MS=KFAS_K240, SPDE=sdmTMB,
##   CF-STM=cf_dglm.
###############################################################################

RES <- "results"
rn  <- function(m){ m[m=="KFAS_K120"]<-"GSSM"; m[m=="KFAS_K240"]<-"GSSM-MS"
                    m[m=="sdmTMB"]<-"SPDE"; m[m=="cf_dglm"]<-"CF-STM"; m }
FAM <- c("gaussian","poisson","binary")
FLAB<- c(gaussian="Gaussian", poisson="Poisson", binary="Binomial")
rd  <- function(f) if(file.exists(f)) read.csv(f) else NULL

hr  <- function() cat(strrep("-",70),"\n")
sec <- function(x) cat("\n",strrep("=",70),"\n",x,"\n",strrep("=",70),"\n",sep="")

## ---------------------------------------------------------------- Table 1 ----
sec("Table 1: Predictive accuracy on the regular panel (RMSE / CRPS)")
ord1 <- c("GLM","GAM","GSSM","GSSM-MS","SPDE","CF-STM")
for(rho in c("02","07")){
  d <- rd(file.path(RES,sprintf("table1_2_panel_rho%s.csv",rho)))
  if(is.null(d)){ cat("  [missing results/table1_2_panel_rho",rho,".csv]\n",sep=""); next }
  d$method <- rn(d$method)
  a <- aggregate(cbind(rmse,crps)~family+method, d, function(z) mean(z,na.rm=TRUE))
  cat(sprintf("\n  rho = %.1f\n", as.numeric(rho)/10))
  cat(sprintf("  %-8s","Method"))
  for(fm in FAM) cat(sprintf("  %14s",FLAB[fm])); cat("\n")
  cat(sprintf("  %-8s",""))
  for(fm in FAM) cat(sprintf("  %6s %6s","RMSE","CRPS")); cat("\n")
  for(m in ord1){ cat(sprintf("  %-8s",m))
    for(fm in FAM){ r<-a[a$family==fm & a$method==m,]
      if(nrow(r)) cat(sprintf("  %6.3f %6.3f",r$rmse,r$crps)) else cat(sprintf("  %6s %6s","-","-")) }
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

## ---------------------------------------------------------- Tables 3 & 4 ----
## Coefficient bias / SD (Table 3) and SE accuracy / coverage (Table 4),
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

  sec("Table 4: Accuracy of the SEs of beta1 on the regular panel (rho = 0.7)")
  cat(sprintf("  %-8s",""));for(fm in FAM)cat(sprintf("  %16s",FLAB[fm]));cat("\n")
  cat(sprintf("  %-8s","Method"));for(fm in FAM)cat(sprintf("  %8s %7s","MeanSE","Cov95"));cat("\n")
  for(m in ord34){ cat(sprintf("  %-8s",m))
    for(fm in FAM){ g<-agg34(fm,m); cat(sprintf("  %8.3f %7.3f",g$mse,g$cov)) }; cat("\n") }
}

## ---------------------------------------------------------------- Table 5 ----
sec("Table 5: Computation time for regular panel data (seconds)")
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

## ---------------------------------------------------------------- Table 6 ----
sec("Table 6: Computation time of CF-STM for regular panel data (seconds)")
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
sec("Table A1: Effect of the SPDE mesh cutoff (moving-site Gaussian setting)")
a1 <- rd(file.path(RES,"tableA1_spde_mesh.csv"))
if(is.null(a1)) cat("  [missing results/tableA1_spde_mesh.csv]\n") else {
  a1$method <- rn(a1$method)
  ag <- aggregate(cbind(rmse,crps,time) ~ method + cutoff, a1, function(z) mean(z,na.rm=TRUE))
  ref <- aggregate(cbind(rmse,crps,time) ~ method, a1[is.na(a1$cutoff),], function(z) mean(z,na.rm=TRUE))
  cat(sprintf("  %-8s %8s %8s %8s %10s\n","Method","cutoff","RMSE","CRPS","Time (s)"))
  for(i in seq_len(nrow(ref)))
    cat(sprintf("  %-8s %8s %8.3f %8.3f %10.1f\n", ref$method[i], "-", ref$rmse[i], ref$crps[i], ref$time[i]))
  sp <- ag[ag$method!="CF-STM",]; sp <- sp[order(-sp$cutoff),]
  for(i in seq_len(nrow(sp)))
    cat(sprintf("  %-8s %8.3f %8.3f %8.3f %10.1f\n",
                if(i==1) sp$method[i] else "", sp$cutoff[i], sp$rmse[i], sp$crps[i], sp$time[i]))
}

## ---------------------------------------------------------------- Table A2 ---
## Appendix 2: predictive accuracy under the Gaussian-process (Matern) DGP.
sec("Table A2: Predictive accuracy under the Gaussian-process DGP (RMSE / CRPS)")
for(rho in c("02","07")){
  d <- rd(file.path(RES,sprintf("tableA2_gp_rho%s.csv",rho)))
  if(is.null(d)){ cat("  [missing results/tableA2_gp_rho",rho,".csv]\n",sep=""); next }
  d$method <- rn(d$method)
  a <- aggregate(cbind(rmse,crps)~family+method, d, function(z) mean(z,na.rm=TRUE))
  cat(sprintf("\n  rho = %.1f\n", as.numeric(rho)/10))
  cat(sprintf("  %-8s","Method")); for(fm in FAM) cat(sprintf("  %14s",FLAB[fm])); cat("\n")
  cat(sprintf("  %-8s","")); for(fm in FAM) cat(sprintf("  %6s %6s","RMSE","CRPS")); cat("\n")
  for(m in ord1){ cat(sprintf("  %-8s",m))
    for(fm in FAM){ r<-a[a$family==fm & a$method==m,]
      if(nrow(r)) cat(sprintf("  %6.3f %6.3f",r$rmse,r$crps)) else cat(sprintf("  %6s %6s","-","-")) }
    cat("\n") }
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
