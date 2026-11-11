#!/usr/bin/env bash
###############################################################################
## run_all.sh  --  regenerate every result CSV for Section 4 Tables 1-5 and
## Supplementary Tables A1-A7,
## then rebuild the formatted tables.
##
## Run from the sim_repro/ root:      bash run_all.sh
##
## CF-STM is compiled/loaded from the STANDALONE sources in engine/ + src/
## (the scripts do NOT load an installed spCF package).
##
## Requirements (R packages):
##   CF-STM engine : Rcpp, FNN, fields, dbscan, nloptr, withr, Matrix
##   competitors   : mgcv (GAM), KFAS (GSSM/GSSM-MS), sdmTMB (SPDE)
##   scoring/util  : scoringRules, parallel
##
## Cost warning: the competitor methods are heavy. On ~12 cores the predictive
## and coefficient studies (100 seeds x 3 families) take on the order of an hour
## EACH; the Table 4 timing grid includes single fits that run ~40 min apiece
## (GSSM-MS at N=800, T=200). Lower NSEED / NREP for a quick smoke run.
###############################################################################
set -euo pipefail
cd "$(dirname "$0")"                      # always run from the bundle root
export NSEED="${NSEED:-100}"              # replicates for the accuracy studies
export NCORE="${NCORE:-12}"              # parallel workers
NREP="${NREP:-5}"                         # reps per cell for the timing studies

echo "### Table 1 & 2: predictive accuracy, regular panel (rho = 0.2, 0.7) ###"
RHO=0.2 OUT=results/table1_2_panel_rho02.csv Rscript R/table1_predictive_panel.R
RHO=0.7 OUT=results/table1_2_panel_rho07.csv Rscript R/table1_predictive_panel.R

echo "### Table 2: predictive accuracy, irregular (moving) sites (rho = 0.2, 0.7) ###"
RHO=0.2 OUT=results/table2_moving_rho02.csv Rscript R/table2_predictive_moving.R
RHO=0.7 OUT=results/table2_moving_rho07.csv Rscript R/table2_predictive_moving.R

echo "### Table 3 and Table A4: coefficient bias / SD / SE / coverage (rho = 0.7) ###"
RHO=0.7 OUT=results/table3_4_coef_se_rho07.csv Rscript R/table3_4_coef_se.R

echo "### Table 4: computation-time comparison (append-mode; fresh file) ###"
rm -f results/table5_timing.csv
NREP="$NREP" OUT=results/table5_timing.csv Rscript R/table5_timing.R

echo "### Table 5: CF-STM scaling (append-mode; fresh file) ###"
rm -f results/table6_scaling.csv
NREP="$NREP" RHO=0.7 OUT=results/table6_scaling.csv Rscript R/table6_scaling.R

echo "### Table A1: SPDE mesh-cutoff sensitivity (5 reps; very slow) ###"
rm -f results/tableA1_panel_replicates.csv results/tableA1_moving_replicates.csv results/tableA1_spde_mesh.csv
DESIGN=panel NSEED=5 RHO=0.7 OUT=results/tableA1_panel_replicates.csv Rscript R/tableA1_spde_mesh.R
DESIGN=moving NSEED=5 RHO=0.7 OUT=results/tableA1_moving_replicates.csv Rscript R/tableA1_spde_mesh.R
Rscript -e 'd<-rbind(read.csv("results/tableA1_panel_replicates.csv"),read.csv("results/tableA1_moving_replicates.csv"));write.csv(d,"results/tableA1_spde_mesh.csv",row.names=FALSE)'

echo "### Table A2 (Appendix 2): predictive accuracy under the Gaussian-process DGP ###"
RHO=0.2 OUT=results/tableA2_gp_rho02.csv Rscript R/tableA2_gp_predictive.R
RHO=0.7 OUT=results/tableA2_gp_rho07.csv Rscript R/tableA2_gp_predictive.R

echo "### Table A3 (Appendix 2): coefficient accuracy under the Gaussian-process DGP ###"
RHO=0.7 OUT=results/tableA3_gp_coef_rho07.csv Rscript R/tableA3_gp_coef.R

echo "### Tables A5-A6: shared temporal-parameter sensitivity ###"
DGLM_CPP=sensitivity/src/dglm_chunk.cpp ENGINE_MOD=sensitivity/engine GENS=sensitivity/R/gens.R OUTDIR=results NSEED="$NSEED" NCORE="$NCORE" Rscript R/tableA5_A6_temporal_sensitivity.R
DGLM_CPP=sensitivity/src/dglm_chunk.cpp ENGINE_MOD=sensitivity/engine GENS=sensitivity/R/gens.R OUTDIR=results NSEED="$NSEED" NCORE="$NCORE" Rscript R/tableA5_scalewise_ml.R

echo "### Table A7: predictive-distribution calibration (very slow) ###"
for design in panel moving gp; do
  for rho in 0.2 0.7; do
    DESIGN="$design" RHO="$rho" NSEED="$NSEED" NCORE="$NCORE" NSIM=500 SAVE_PREDS=0 \
      ENGINE_MOD=sensitivity/engine DGLM_CPP=sensitivity/src/dglm_chunk.cpp OUTDIR=results \
      Rscript R/tableA7_predictive_calibration.R
  done
done
rm -rf results_gssm
mkdir -p results_gssm
for design in panel gp; do
  for rho in 0.2 0.7; do
    DESIGN="$design" RHO="$rho" NSEED="$NSEED" NCORE="$NCORE" NSIM=500 SAVE_PREDS=0 FAMS=gaussian METHODS=GSSM,GSSM-MS \
      ENGINE_MOD=sensitivity/engine DGLM_CPP=sensitivity/src/dglm_chunk.cpp OUTDIR=results_gssm \
      Rscript R/tableA7_predictive_calibration.R
  done
done
for f in results_gssm/coverage_*.csv; do
  base="${f#results_gssm/}"
  Rscript -e 'a<-read.csv(commandArgs(TRUE)[1]);b<-read.csv(commandArgs(TRUE)[2]);a<-a[!a$method%in%c("GSSM","GSSM-MS"),];write.csv(rbind(a,b),commandArgs(TRUE)[1],row.names=FALSE)' "results/$base" "$f"
done
for f in results/coverage_*.csv; do mv "$f" "results/tableA7_${f#results/}"; done

echo "### Building formatted Tables 1-5 and A1-A7 ###"
Rscript build_tables.R
