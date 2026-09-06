#!/usr/bin/env Rscript
###############################################################################
## record_env.R -- record the run environment of this reproduction bundle.
##
## Usage, from the bundle root:
##     Rscript record_env.R
##
## Overwrites sessionInfo.txt with the R version, the versions of every package
## the bundle needs, the C++ toolchain, and the CPU / core count (the timing
## tables 5, 6 and A1 are wall-clock measurements, so the hardware matters).
##
## Run this on a machine where ALL packages below are installed; the script
## writes a prominent warning into the file if any of them is missing, so an
## incomplete record can never be mistaken for a complete one.
###############################################################################

OUT  <- "sessionInfo.txt"
PKGS <- c("Rcpp","FNN","fields","dbscan","nloptr","withr","Matrix",
          "mgcv","KFAS","sdmTMB","scoringRules","parallel")

have    <- vapply(PKGS, function(p) requireNamespace(p, quietly = TRUE), logical(1))
missing <- PKGS[!have]

## Attach what is available so that sessionInfo() reports it explicitly.
invisible(lapply(PKGS[have], function(p)
  suppressMessages(suppressWarnings(library(p, character.only = TRUE)))))

## ---- helpers ---------------------------------------------------------------
sh <- function(cmd, args = character()) {
  tryCatch(suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = FALSE))[1],
           error = function(e) NA_character_)
}
cpu_model <- function() {
  os <- Sys.info()[["sysname"]]
  v <- if (identical(os, "Darwin")) sh("sysctl", c("-n", "machdep.cpu.brand_string"))
       else if (identical(os, "Linux")) {
         l <- tryCatch(grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[1],
                       error = function(e) NA_character_)
         if (!is.na(l)) sub("^model name\\s*:\\s*", "", l) else NA_character_
       } else NA_character_
  if (is.na(v) || !nzchar(v)) "(unknown)" else v
}
cxx_info <- function() {
  cxx <- tryCatch(sh(file.path(R.home("bin"), "R"), c("CMD", "config", "CXX")),
                  error = function(e) NA_character_)
  if (is.na(cxx) || !nzchar(cxx)) return(c("(unknown)", "(unknown)"))
  ver <- sh(strsplit(cxx, " ")[[1]][1], "--version")
  c(cxx, if (is.na(ver) || !nzchar(ver)) "(version unavailable)" else ver)
}

con <- file(OUT, "w"); w <- function(...) writeLines(paste0(...), con)

w("Environment record for the CF-STM reproduction bundle")
w("Recorded: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  "   (regenerate with: Rscript record_env.R)")
w(strrep("=", 78)); w("")

if (length(missing)) {
  w("WARNING -- INCOMPLETE ENVIRONMENT")
  w(strrep("-", 78))
  w("The following package(s) are NOT installed on this machine, so the bundle")
  w("cannot be run end to end here and this file is not a complete record of the")
  w("environment that produced results/:")
  for (p in missing) w("    ", p)
  w("")
  w("Re-run  Rscript record_env.R  on a machine where every package is installed.")
} else {
  w("Every package required to reproduce all tables is installed here, so this is")
  w("a complete record of a working run environment.")
}
w("")

w("Required packages")
w(strrep("-", 78))
for (p in PKGS)
  w(sprintf("  %-14s %s", p,
            if (have[[p]]) as.character(utils::packageVersion(p)) else "-- NOT INSTALLED --"))
w("")

w("Hardware and toolchain")
w(strrep("-", 78))
w(sprintf("  %-14s %s", "platform",  R.version$platform))
w(sprintf("  %-14s %s %s", "OS", Sys.info()[["sysname"]], Sys.info()[["release"]]))
w(sprintf("  %-14s %s", "CPU",       cpu_model()))
w(sprintf("  %-14s %s", "cores",     tryCatch(as.character(parallel::detectCores()),
                                              error = function(e) "(unknown)")))
ci <- cxx_info()
w(sprintf("  %-14s %s", "C++ compiler", ci[1]))
w(sprintf("  %-14s %s", "",            ci[2]))
w("")
w("  src/dglm_chunk.cpp is compiled on first use by Rcpp::sourceCpp() via")
w("  .dglm_load_cpp(), so a working C++ toolchain is required.")
w("")
w("  Tables 5, 6 and A1 report wall-clock times and are therefore hardware-")
w("  dependent; only the orders-of-magnitude gaps are meant to reproduce.")
w("")

w("sessionInfo()")
w(strrep("-", 78))
close(con)
cat(capture.output(sessionInfo()), file = OUT, sep = "\n", append = TRUE)

## ---- console summary -------------------------------------------------------
if (length(missing)) {
  message(sprintf("Wrote %s -- WARNING: %d package(s) missing: %s",
                  OUT, length(missing), paste(missing, collapse = ", ")))
} else {
  message(sprintf("Wrote %s -- complete environment (%d packages).", OUT, length(PKGS)))
}
