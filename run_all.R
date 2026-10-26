## Run every reproducibility script and save its output to output/<script>.txt.
## The two slow scripts are skipped unless requested:
##   Rscript run_all.R            # all except 09 and 10
##   Rscript run_all.R slow       # everything (several hours, mostly the redundancy check in 10)
args  <- commandArgs(trailingOnly = TRUE)
slow  <- length(args) > 0 && args[1] == "slow"
files <- sort(list.files(pattern = "^[0-9]{2}-.*[.]R$"))
if (!slow) files <- setdiff(files, c("09-tripartite-bell.R", "10-benchmarks.R"))
dir.create("output", showWarnings = FALSE)
for (f in files) {
  out <- file.path("output", sub("[.]R$", ".txt", f))
  cat("Running", f, "->", out, "\n")
  t <- system.time(res <- system2(file.path(R.home("bin"), "Rscript"), f, stdout = out, stderr = out))
  cat(sprintf("  finished in %.1f s (exit status %d)\n", t["elapsed"], res))
}
