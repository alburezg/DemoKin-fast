# Run every scenario against one implementation and save outputs + timings.
#
# Usage: Rscript run_scenarios.R <package_dir> <out_dir> [scenario_regex]
#
# The package is loaded from source with pkgload::load_all() so the exact
# source tree is what gets timed, not whatever happens to be installed.

args <- commandArgs(trailingOnly = TRUE)
pkg_dir <- args[1]
out_dir <- args[2]
pattern <- if (length(args) >= 3) args[3] else NULL

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(pkgload::load_all(pkg_dir, export_all = FALSE, quiet = TRUE))

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(script_file), "scenarios.R"))

inp <- build_inputs()
scen <- build_scenarios(inp)
if (!is.null(pattern)) scen <- scen[grepl(pattern, names(scen))]

timings <- data.frame(scenario = names(scen), seconds = NA_real_, status = NA_character_, stringsAsFactors = FALSE)
for (i in seq_along(scen)) {
  nm <- names(scen)[i]
  cat(sprintf("[%s] %-28s ", format(Sys.time(), "%H:%M:%S"), nm))
  gc()
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(scen[[nm]](), error = function(e) structure(list(message = conditionMessage(e)), class = "scenario_error"))
  el <- proc.time()[["elapsed"]] - t0
  timings$seconds[i] <- el
  timings$status[i] <- if (inherits(res, "scenario_error")) paste("ERROR:", res$message) else "ok"
  saveRDS(res, file.path(out_dir, paste0(nm, ".rds")))
  cat(sprintf("%8.2f s  %s\n", el, timings$status[i]))
  rm(res)
}
write.csv(timings, file.path(out_dir, "_timings.csv"), row.names = FALSE)
cat("done\n")
