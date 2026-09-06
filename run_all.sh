#!/usr/bin/env bash
###############################################################################
## run_all.sh  --  regenerate every result CSV for Section 4 Tables 1-6,
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
## EACH; the Table 5 timing grid includes single fits that run ~40 min apiece
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

echo "### Table 3 & 4: coefficient bias / SD / SE / coverage (rho = 0.7) ###"
RHO=0.7 OUT=results/table3_4_coef_se_rho07.csv Rscript R/table3_4_coef_se.R

echo "### Table 5: computation-time comparison (append-mode; fresh file) ###"
rm -f results/table5_timing.csv
NREP="$NREP" OUT=results/table5_timing.csv Rscript R/table5_timing.R

echo "### Table 6: CF-STM scaling (append-mode; fresh file) ###"
rm -f results/table6_scaling.csv
NREP="$NREP" RHO=0.7 OUT=results/table6_scaling.csv Rscript R/table6_scaling.R

echo "### Table A1 (Appendix 1): SPDE mesh-cutoff sensitivity (5 reps; slow) ###"
rm -f results/tableA1_spde_mesh.csv
NSEED=5 RHO=0.7 OUT=results/tableA1_spde_mesh.csv Rscript R/tableA1_spde_mesh.R

echo "### Table A2 (Appendix 2): predictive accuracy under the Gaussian-process DGP ###"
RHO=0.2 OUT=results/tableA2_gp_rho02.csv Rscript R/tableA2_gp_predictive.R
RHO=0.7 OUT=results/tableA2_gp_rho07.csv Rscript R/tableA2_gp_predictive.R

echo "### Table A3 (Appendix 2): coefficient accuracy under the Gaussian-process DGP ###"
RHO=0.7 OUT=results/tableA3_gp_coef_rho07.csv Rscript R/tableA3_gp_coef.R

echo "### Building formatted Tables 1-6 and A1-A3 ###"
Rscript build_tables.R
