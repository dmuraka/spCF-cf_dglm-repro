# Reproducibility: Section 4 (Monte Carlo experiments), Tables 1–6

This bundle reproduces every table in Section 4 of

> Murakami et al. (2026), *Coarse-to-fine dynamic spatio-temporal models*
> (CF-STM / CF-DGLMM).

CF-STM is provided here as **standalone R + C++ source** in `engine/` and `src/`.
The scripts load it with `source(...)`; **no installed `spCF` package is
required or assumed** (CF-STM is not yet part of `spCF` at submission time).

---

## Quick start

Build the tables from the shipped per-replicate result CSVs (takes ~1 second,
compiles nothing):

```r
# from this directory
Rscript build_tables.R    # prints Tables 1-6 and A1-A3
```

Regenerate the result CSVs from scratch (heavy — see cost warning below), then
rebuild:

```sh
bash run_all.sh              # NSEED=100 NCORE=12 NREP=5 by default
# quick smoke run:
NSEED=10 NREP=1 bash run_all.sh
```

Both must be run from the `sim_repro/` root (the CF-STM C++ kernel is located at
`./src/dglm_chunk.cpp` relative to the working directory).

---

## Layout

```
sim_repro/
├── README.md
├── LICENSE                        # GPL-2-or-later (same terms as spCF)
├── CITATION.cff                   # machine-readable citation metadata
├── .zenodo.json                   # Zenodo deposit metadata
├── sessionInfo.txt                # package versions / environment record
├── record_env.R                   # one command that regenerates sessionInfo.txt
├── build_tables.R                 # CSVs in results/ -> formatted Tables 1-6, A1-A3
├── run_all.sh                     # regenerate all CSVs, then build_tables.R
├── engine/                        # STANDALONE CF-STM (no spCF package)
│   ├── cf_dglm.R                  #   fit + prediction (CF-STM)
│   ├── cf_dglm_hv.R               #   holdout scale/hyper-parameter selection
│   ├── internal_utils_dglm.R      #   internals + .dglm_load_cpp()
│   └── sp_scalewise.R             #   (scale-wise field extraction; not used here)
├── src/
│   └── dglm_chunk.cpp             # fused C++ scale operator (compiled on first use)
├── R/
│   ├── table1_predictive_panel.R  # Table 1 + Table 2 (regular panel)
│   ├── table2_predictive_moving.R # Table 2 (irregular / moving sites)
│   ├── table3_4_coef_se.R         # Table 3 + Table 4
│   ├── table5_timing.R            # Table 5
│   ├── table6_scaling.R           # Table 6
│   ├── tableA1_spde_mesh.R        # Table A1 (Appendix 1)
│   ├── tableA2_gp_predictive.R    # Table A2 (Appendix 2)
│   └── tableA3_gp_coef.R          # Table A3 (Appendix 2)
└── results/                       # per-replicate CSVs (inputs to build_tables.R)
    ├── table1_2_panel_rho02.csv   ├── table1_2_panel_rho07.csv
    ├── table2_moving_rho02.csv    ├── table2_moving_rho07.csv
    ├── table3_4_coef_se_rho07.csv
    ├── table5_timing.csv          ├── table6_scaling.csv
    ├── tableA1_spde_mesh.csv      (+ .log)
    ├── tableA2_gp_rho02.csv       ├── tableA2_gp_rho07.csv
    └── tableA3_gp_coef_rho07.csv
```

---

## Table → script → data map

| Paper table | What it reports | Script (`R/`) | Result CSV (`results/`) |
|---|---|---|---|
| **1** | Predictive RMSE/CRPS, regular panel, ρ∈{0.2,0.7} | `table1_predictive_panel.R` | `table1_2_panel_rho{02,07}.csv` |
| **2** | RMSE of SPDE vs CF-STM, regular **and** irregular | `table1_predictive_panel.R` (regular) + `table2_predictive_moving.R` (irregular) | `table1_2_panel_rho{02,07}.csv`, `table2_moving_rho{02,07}.csv` |
| **3** | Bias & std.dev of β̂, regular panel, ρ=0.7 | `table3_4_coef_se.R` | `table3_4_coef_se_rho07.csv` |
| **4** | Mean SE & 95% coverage of β̂, ρ=0.7 | `table3_4_coef_se.R` | `table3_4_coef_se_rho07.csv` |
| **5** | Computation time, all methods, regular panel | `table5_timing.R` | `table5_timing.csv` |
| **6** | CF-STM computation time (estimation + prediction), larger panels | `table6_scaling.R` | `table6_scaling.csv` |
| **A1** | SPDE mesh-cutoff sensitivity (RMSE/CRPS/time), moving-site Gaussian | `tableA1_spde_mesh.R` | `tableA1_spde_mesh.csv` |
| **A2** | Predictive RMSE/CRPS under the **Gaussian-process DGP**, ρ∈{0.2,0.7} | `tableA2_gp_predictive.R` | `tableA2_gp_rho{02,07}.csv` |
| **A3** | Bias & std.dev of β₁ under the **Gaussian-process DGP**, ρ=0.7 | `tableA3_gp_coef.R` | `tableA3_gp_coef_rho07.csv` |

### Method-label mapping

Paper labels vs the `method` column in the CSVs:

| Paper | CSV `method` | Implementation |
|---|---|---|
| GLM | `GLM` | `stats::glm` |
| GAM | `GAM` | `mgcv::gam` with `te(cx,cy,tt)` |
| GSSM | `KFAS_K120` | `KFAS` state-space, single scale (K=120 knots) |
| GSSM-MS | `KFAS_K240` | `KFAS` state-space, multiscale (K=40+80+120) |
| SPDE | `sdmTMB` | `sdmTMB`, `spatiotemporal="ar1"` |
| CF-STM | `cf_dglm` | this bundle's `engine/` |

---

## Data-generating process

A spatio-temporal AR(1) latent field on `n_pop = 400` sites over `T = 40` time
points, `range = 0.10`, plus two spatially structured covariates. True slopes
per family (on the link scale):

- Gaussian: `(2, −1.5)`,  additive N(0, 0.5²) noise
- Poisson:  `(0.3, −0.2)`, intercept 0.5, field loading 0.7
- Binomial: `(0.5, −0.4)`, field loading 1.2

**Appendix 2 (Tables A2, A3) replaces this field.** The innovations and the two
covariates are instead drawn from a zero-mean Matérn Gaussian process
(smoothness 1, marginal variance 1) evolving as a temporal AR(1) — i.e. exactly
the covariance the SPDE benchmark assumes. Everything else (sizes, slopes,
holdout, replicates, methods) is unchanged, so Tables A2/A3 are directly
comparable with Tables 1/3.

`RHO ∈ {0.2, 0.7}` sets the temporal AR(1) coefficient. The **regular panel**
(Tables 1, 3–6) keeps the same sites across time; the **irregular** case
(Table 2) resamples fresh sites at every time point (the latent field lives on
fixed anchors and is kernel-interpolated to each time's sites). Predictive
accuracy holds out a random 30% of cells. Every replicate is seeded by its
index, so runs are reproducible.

---

## How to cite

Please cite the paper:

> Murakami, D. (2026). Fast covariance-free spatiotemporal modeling via
> coarse-to-fine learning. *TODO: journal, volume, pages.*
> doi:TODO

and, where the exact computational artifact matters, the archived snapshot of
this bundle:

> Murakami, D. (2026). *Reproduction code for "Fast covariance-free
> spatiotemporal modeling via coarse-to-fine learning"* (v1.0.0) [Software].
> Zenodo. doi:10.5281/zenodo.TODO

The CF-STM implementation itself is distributed in the R package **spCF**
(<https://github.com/dmuraka/spCF>); the copy in `engine/` and `src/` is the
frozen snapshot that produced the results reported in the paper.

Repository: <https://github.com/dmuraka/spCF-cf_dglm-repro>

`CITATION.cff` carries the same information in machine-readable form. Before
archiving, fill in the remaining `TODO` placeholders (ORCID, Zenodo DOI, and the
journal reference) in `CITATION.cff`, `.zenodo.json`, and this section.

---

## Notes

- **License.** GPL-2 or later, matching the spCF package from which the CF-STM
  sources in `engine/` and `src/` are taken (see `LICENSE`).
- **Environment.** `sessionInfo.txt` records the R version, the version of every
  required package, the C++ toolchain, and the CPU / core count (the timing
  tables are wall-clock measurements). Regenerate it with one command, from the
  bundle root:

  ```sh
  Rscript record_env.R
  ```

  Run it on a machine where every required package is installed. If one is
  missing, the script still writes the file but puts a prominent warning at its
  head, so an incomplete record cannot be mistaken for a complete one. The
  shipped copy carries such a warning: it was captured where `KFAS` was absent,
  so re-record it before archiving.
- **Requirements.** R with `Rcpp`, `FNN`, `fields`, `dbscan`, `nloptr`,
  `withr`, `Matrix` (CF-STM engine); `mgcv`, `KFAS`, `sdmTMB`,
  `scoringRules`, `parallel` (competitors, scoring, parallelism). A C++
  toolchain is needed to compile `src/dglm_chunk.cpp` on first use.
- **Timing tables (5 & 6) are machine-dependent.** The shipped values are
  5-replicate means from a single workstation and will differ on other
  hardware; only the orders-of-magnitude gaps are meant to be reproducible.
  The most expensive competitor cell in Table 5 (GSSM-MS at N=800, T=200,
  ≈44 min per fit) was measured in a separate single run.
- All shipped `results/` CSVs reproduce the printed tables: Tables 1–6 and
  A2, A3 match the paper to the printed precision, and Table A1 matches on
  RMSE/CRPS. The Table A1 *timing* column differs from the paper by about 1%
  (49.1 vs 49.7 s at cutoff 0.10; 8789.6 vs 8790.5 s at 0.025), because the
  paper's timings were transcribed from an earlier run on the same machine;
  wall-clock times are machine- and run-dependent in any case.
- `R/tableA1_spde_mesh.R` does not run the cutoff 0.01 case by default: a single
  fit was projected to exceed ten hours (as stated in Appendix 1). Add it with
  `CUTS="0.2,0.1,0.05,0.025,0.01"`.
- `build_tables.R` performs only aggregation/formatting and prints a note for
  any missing CSV, so you can rebuild individual tables without rerunning the
  whole study.
