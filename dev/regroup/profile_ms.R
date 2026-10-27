# Profile the multistage functions (development source tree).
#   Rscript dev/regroup/profile_ms.R <which>   which = parity | svk | ms
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))   # run from the package root
which <- commandArgs(TRUE)[1]
q <- function(e) suppressWarnings(suppressMessages(e))
f <- switch(which,
  parity = function() q(kin_multi_stage_time_variant_2sex(
    U_list_females = Female_parity_mortality_list_UK[1:3], U_list_males = Male_parity_mortality_list_UK[1:3],
    F_list_females = Female_parity_fert_list_UK[1:3], F_list_males = Female_parity_fert_list_UK[1:3],
    T_list_females = Parity_transfers_by_age_list_UK[1:3], T_list_males = Parity_transfers_by_age_list_UK[1:3],
    birth_female = 1 - 0.51, parity = TRUE, output_kin = c("d", "oa", "ys", "os"), summary_kin = TRUE,
    sex_Focal = "Female", initial_stage_Focal = 1, output_years = 1965:1967)),
  parity12 = function() q(kin_multi_stage_time_variant_2sex(
    U_list_females = Female_parity_mortality_list_UK[1:13], U_list_males = Male_parity_mortality_list_UK[1:13],
    F_list_females = Female_parity_fert_list_UK[1:13], F_list_males = Female_parity_fert_list_UK[1:13],
    T_list_females = Parity_transfers_by_age_list_UK[1:13], T_list_males = Parity_transfers_by_age_list_UK[1:13],
    birth_female = 1 - 0.51, parity = TRUE, output_kin = c("d", "oa", "ys", "os"), summary_kin = TRUE,
    sex_Focal = "Female", initial_stage_Focal = 1, output_years = c(1965, 1970, 1977))),
  svk = function() {
    Ff <- (1 / 0.49) * svk_fxs
    q(kin_multi_stage_time_variant_2sex(
      list(svk_pxs), list(svk_pxs), list(Ff), list(Ff), list(svk_Uxs), list(svk_Uxs),
      birth_female = 0.49, output_kin = NULL, parity = TRUE, summary_kin = FALSE,
      sex_Focal = "Female", initial_stage_Focal = 1, output_years = seq(1990, 1990)))
  },
  ms = function() q(kin_multi_stage(U = svk_Uxs, f = svk_fxs, D = svk_pxs, H = svk_Hxs,
                                    birth_female = 0.49, parity = TRUE)))
prof_file <- tempfile(fileext = ".prof")
invisible(tryCatch(f(), error = function(e) message("first call: ", conditionMessage(e))))
t0 <- proc.time()[[3]]
Rprof(prof_file, interval = 0.01)
x <- f()
Rprof(NULL)
cat(which, ": ", proc.time()[[3]] - t0, "s\n")
s <- summaryRprof(prof_file)
print(head(s$by.total[, c("total.time", "total.pct", "self.time")], 30))
print(head(s$by.self[, c("self.time", "self.pct")], 15))
