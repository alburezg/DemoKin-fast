# Compare DemoKin (original, ./DemoKin) with DemoKin_fast (./DemoKin_fast)
#
# Runs the same calls through both versions and reports, for each call, whether
# the outputs are identical, where they differ (if they do), and how long each
# version took.
#
# Both packages are called "DemoKin", so they cannot be loaded in the same R
# session. Each version therefore runs in its own R process (callr), loaded
# from source with pkgload::load_all(), and saves its outputs to disk. This
# session then reads both sets of outputs back and compares them.
#
# What "identical" means here: identical() on the whole object after removing
# data.table's ".internal.selfref" attribute (a memory pointer that is never
# equal across sessions). That checks every value bit for bit, plus column
# names and order, row order, column types, NA/NaN positions, and class
# (data.frame / data.table / tibble). Any difference is printed with waldo.
#
# Usage (from the DemoKin project folder):
#   Rscript compare_DemoKin_original_vs_fast.R          # all 44 calls, ~12 min
#   Rscript compare_DemoKin_original_vs_fast.R quick    # skip the slowest calls, ~3 min
#
# Output: comparison_original_vs_fast/summary.csv and report.txt
#
# Requires: callr, pkgload, waldo, plus the packages both versions import.

root    <- normalizePath(".", winslash = "/")
old_dir <- file.path(root, "DemoKin")
new_dir <- file.path(root, "DemoKin_fast")
out_dir <- file.path(root, "comparison_original_vs_fast")
quick   <- "quick" %in% commandArgs(trailingOnly = TRUE)
stopifnot(dir.exists(old_dir), dir.exists(new_dir))

# ---------------------------------------------------------------------------
# The calls to compare. Each element: the exported function being exercised,
# a short description, and a zero-argument function that makes the call.
# Inputs only use data bundled with the package. `slow` marks calls that take
# more than ~15 s with the original package (skipped in quick mode).
# ---------------------------------------------------------------------------
scenarios <- function() {
  q <- function(expr) suppressMessages(suppressWarnings(expr))

  # inputs -------------------------------------------------------------------
  p_2015 <- swe_px[, "2015"]
  f_2015 <- swe_asfr[, "2015"]
  pi_2015 <- swe_pop[, "2015"] * swe_asfr[, "2015"]; pi_2015 <- pi_2015 / sum(pi_2015)
  pi_mat <- t(t(swe_pop * swe_asfr) / colSums(swe_pop * swe_asfr))

  fra_pf <- fra_surv_sex[, "pf"]; fra_pm <- fra_surv_sex[, "pm"]
  fra_ff <- fra_asfr_sex[, "ff"]; fra_fm <- fra_asfr_sex[, "fm"]
  fra_pif <- fra_ff * (0.3 + seq_len(101) / 101); fra_pif <- fra_pif / sum(fra_pif)
  fra_pim <- fra_fm * (0.5 + seq_len(101) / 101); fra_pim <- fra_pim / sum(fra_pim)

  # two-sex time-varying rates, as built in vignette 1_4
  years <- ncol(swe_px); ages <- nrow(swe_px)
  swe_pm <- swe_px^1.5
  swe_fm <- rbind(matrix(0, 5, years), swe_asfr[-((ages - 4):ages), ]) * 1.05
  dimnames(swe_fm) <- dimnames(swe_asfr)

  # causes of death (vignette 1_4, plus a 3-cause and a time-varying variant)
  Hf2 <- matrix(c(.5, 1), nrow = 2, ncol = 101); Hm2 <- matrix(c(.33, 1), nrow = 2, ncol = 101)
  Hf2[, 1:50] <- Hm2[, 1:50] <- 1
  Hf3 <- matrix(c(.5, 1, .25), nrow = 3, ncol = 101); Hm3 <- matrix(c(.33, 1, .4), nrow = 3, ncol = 101)
  Hf3[, 1:40] <- Hm3[, 1:40] <- 1
  y_cod <- as.character(1990:2018)
  Hf2_list <- stats::setNames(lapply(seq_along(y_cod), function(i) Hf2 * (1 + i / 100)), y_cod)
  Hm2_list <- stats::setNames(lapply(seq_along(y_cod), function(i) Hm2 * (1 + i / 200)), y_cod)
  y_cod3 <- as.character(2000:2008)
  Hf3_list <- stats::setNames(lapply(seq_along(y_cod3), function(i) Hf3 * (1 + i / 50)), y_cod3)
  Hm3_list <- stats::setNames(lapply(seq_along(y_cod3), function(i) Hm3 * (1 + i / 60)), y_cod3)

  y00 <- as.character(2000:2018)
  s <- function(fun, desc, call, slow = FALSE) list(fun = fun, desc = desc, call = call, slow = slow)

  list(
    # one sex, time invariant
    ti_1sex_default      = s("kin", "Sweden 2015, defaults",
                             function() q(kin(p = p_2015, f = f_2015))),
    ti_1sex_kinsel       = s("kin", "grouped kin (s, a), selected Focal ages, birth_female",
                             function() q(kin(p = p_2015, f = f_2015, output_kin = c("m", "d", "s", "a", "gm"),
                                              output_age_focal = c(0, 30, 65), birth_female = 0.49))),
    ti_1sex_pi           = s("kin", "user-supplied pi",
                             function() q(kin(p = p_2015, f = f_2015, pi = pi_2015))),
    ti_1sex_nosum        = s("kin", "summary_kin = FALSE",
                             function() q(kin(p = p_2015, f = f_2015, summary_kin = FALSE))),
    ti_1sex_matrix_input = s("kin", "matrix input with time_invariant = TRUE",
                             function() q(kin(p = swe_px[, c("2010", "2011")], f = swe_asfr[, c("2010", "2011")]))),
    ti_1sex_direct_df    = s("kin_time_invariant", "data frame, selected kin",
                             function() q(kin_time_invariant(p = p_2015, f = f_2015, output_kin = c("m", "os")))),
    ti_1sex_direct_list  = s("kin_time_invariant", "list_output = TRUE",
                             function() q(kin_time_invariant(p = p_2015, f = f_2015, list_output = TRUE))),

    # one sex, time varying
    tv_1sex_cohort       = s("kin", "1900-2018, cohorts 1910 and 1960",
                             function() q(kin(p = swe_px, f = swe_asfr, time_invariant = FALSE,
                                              output_cohort = c(1910, 1960))), slow = TRUE),
    tv_1sex_period       = s("kin", "1900-2018, periods 1950, 2000, 2018",
                             function() q(kin(p = swe_px, f = swe_asfr, time_invariant = FALSE,
                                              output_period = c(1950, 2000, 2018)))),
    tv_1sex_period_n     = s("kin", "population n given, grouped sisters",
                             function() q(kin(p = swe_px, f = swe_asfr, n = swe_pop, time_invariant = FALSE,
                                              output_period = 2000, output_kin = c("m", "gm", "s")))),
    tv_1sex_all_small    = s("kin", "2000-2018, all periods",
                             function() q(kin(p = swe_px[, y00], f = swe_asfr[, y00], time_invariant = FALSE))),
    tv_1sex_yearlast     = s("kin", "output_year_last = TRUE",
                             function() q(kin(p = swe_px[, y00], f = swe_asfr[, y00], time_invariant = FALSE,
                                              output_period = c(2017, 2018), output_year_last = TRUE))),
    tv_1sex_bf_vec       = s("kin", "birth_female shorter than years (padded)",
                             function() q(kin(p = swe_px[, as.character(1950:2018)], f = swe_asfr[, as.character(1950:2018)],
                                              time_invariant = FALSE, birth_female = c(0.49, 0.485), output_cohort = 1950))),
    tv_1sex_pi_matrix    = s("kin", "user-supplied pi matrix",
                             function() q(kin(p = swe_px, f = swe_asfr, pi = pi_mat, time_invariant = FALSE,
                                              output_cohort = 1980))),
    tv_1sex_age_focal    = s("kin", "selected Focal ages, no summary",
                             function() q(kin(p = swe_px, f = swe_asfr, time_invariant = FALSE, output_cohort = 1960,
                                              output_age_focal = c(10, 40), summary_kin = FALSE))),
    tv_1sex_direct_list  = s("kin_time_variant", "list_output = TRUE",
                             function() q(kin_time_variant(p = swe_px[, y00], f = swe_asfr[, y00],
                                                           output_cohort = 2000, list_output = TRUE))),
    tv_1sex_direct_df    = s("kin_time_variant", "data frame, one period",
                             function() q(kin_time_variant(p = swe_px[, as.character(2005:2015)], f = swe_asfr[, as.character(2005:2015)],
                                                           output_period = 2010, output_kin = c("d", "m")))),

    # two sex, time invariant
    ti_2sex_default      = s("kin2sex", "France 2012, birth_female = .5",
                             function() q(kin2sex(fra_pf, fra_pm, fra_ff, fra_fm, birth_female = .5))),
    ti_2sex_male         = s("kin2sex", "male Focal, grouped kin, selected ages",
                             function() q(kin2sex(fra_pf, fra_pm, fra_ff, fra_fm, sex_focal = "m",
                                                  output_kin = c("d", "s", "c"), output_age_focal = c(20, 50)))),
    ti_2sex_pi           = s("kin2sex", "user-supplied pif and pim",
                             function() q(kin2sex(fra_pf, fra_pm, fra_ff, fra_fm, pif = fra_pif, pim = fra_pim))),
    ti_2sex_nosum        = s("kin2sex", "summary_kin = FALSE",
                             function() q(kin2sex(fra_pf, fra_pm, fra_ff, fra_fm, summary_kin = FALSE))),
    ti_2sex_direct_list  = s("kin_time_invariant_2sex", "male Focal, list_output = TRUE",
                             function() q(kin_time_invariant_2sex(fra_pf, fra_pm, fra_ff, fra_fm,
                                                                  sex_focal = "m", list_output = TRUE))),

    # two sex, time varying
    tv_2sex_cohort       = s("kin2sex", "1900-2018, cohort 1900 (vignette 1_4)",
                             function() q(kin2sex(pf = swe_px, pm = swe_pm, ff = swe_asfr, fm = swe_fm,
                                                  time_invariant = FALSE, birth_female = .5, output_cohort = 1900)), slow = TRUE),
    tv_2sex_period       = s("kin2sex", "male Focal, nf/nm given, periods 1960 and 2000",
                             function() q(kin2sex(pf = swe_px, pm = swe_pm, ff = swe_asfr, fm = swe_fm,
                                                  nf = swe_pop, nm = swe_pop * 0.98, sex_focal = "m",
                                                  time_invariant = FALSE, output_period = c(1960, 2000))), slow = TRUE),
    tv_2sex_all_small    = s("kin2sex", "2000-2018, all periods",
                             function() q(kin2sex(pf = swe_px[, y00], pm = swe_pm[, y00], ff = swe_asfr[, y00], fm = swe_fm[, y00],
                                                  time_invariant = FALSE)), slow = TRUE),
    tv_2sex_bf_vec_last  = s("kin2sex", "birth_female vector, output_year_last, grouped kin",
                             function() q(kin2sex(pf = swe_px[, y00], pm = swe_pm[, y00], ff = swe_asfr[, y00], fm = swe_fm[, y00],
                                                  time_invariant = FALSE, birth_female = c(.5, .49, .48), output_period = 2018,
                                                  output_year_last = TRUE, output_kin = c("m", "a", "n"))), slow = TRUE),
    tv_2sex_direct_list  = s("kin_time_variant_2sex", "list_output = TRUE",
                             function() q(kin_time_variant_2sex(pf = swe_px[, as.character(2010:2018)], pm = swe_pm[, as.character(2010:2018)],
                                                                ff = swe_asfr[, as.character(2010:2018)], fm = swe_fm[, as.character(2010:2018)],
                                                                output_cohort = 2010, list_output = TRUE))),

    # causes of death
    ti_cod_default       = s("kin2sex", "2 causes (vignette 1_4)",
                             function() q(kin2sex(fra_pf, fra_pm, fra_ff, fra_fm, Hf = Hf2, Hm = Hm2))),
    ti_cod_3causes_male  = s("kin2sex", "3 causes, male Focal, grouped kin",
                             function() q(kin2sex(fra_pf, fra_pm, fra_ff, fra_fm, Hf = Hf3, Hm = Hm3, sex_focal = "m",
                                                  output_kin = c("m", "d", "s"), output_age_focal = c(30, 70)))),
    ti_cod_direct_list   = s("kin_time_invariant_2sex_cod", "list_output = TRUE",
                             function() q(kin_time_invariant_2sex_cod(fra_pf, fra_pm, fra_ff, fra_fm,
                                                                      Hf = Hf2, Hm = Hm2, list_output = TRUE))),
    tv_cod_cohort        = s("kin2sex", "time-varying causes, 1990-2018, cohort 1990",
                             function() q(kin2sex(pf = swe_px[, y_cod], pm = swe_pm[, y_cod], ff = swe_asfr[, y_cod], fm = swe_fm[, y_cod],
                                                  Hf = Hf2_list, Hm = Hm2_list, time_invariant = FALSE, output_cohort = 1990)), slow = TRUE),
    tv_cod_period        = s("kin2sex", "time-varying 3 causes, periods 2003 and 2005",
                             function() q(kin2sex(pf = swe_px[, y_cod3], pm = swe_pm[, y_cod3], ff = swe_asfr[, y_cod3], fm = swe_fm[, y_cod3],
                                                  Hf = Hf3_list, Hm = Hm3_list, time_invariant = FALSE, output_period = c(2003, 2005),
                                                  output_kin = c("m", "s", "gm"))), slow = TRUE),

    # multistage, one sex (Caswell 2020)
    ms_parity_list       = s("kin_multi_stage", "Slovakia parity, list_output (package test)",
                             function() q(kin_multi_stage(U = svk_Uxs, f = svk_fxs, D = svk_pxs, H = svk_Hxs,
                                                          birth_female = 1, list_output = TRUE, parity = TRUE)), slow = TRUE),
    ms_df_kinsel         = s("kin_multi_stage", "data frame, selected kin, parity = FALSE",
                             function() q(kin_multi_stage(U = svk_Uxs, f = svk_fxs, D = svk_pxs, H = svk_Hxs,
                                                          birth_female = 1, parity = FALSE, output_kin = c("m", "d", "os")))),
    ms_df_all            = s("kin_multi_stage", "data frame, all kin, parity = TRUE",
                             function() q(kin_multi_stage(U = svk_Uxs, f = svk_fxs, D = svk_pxs, H = svk_Hxs,
                                                          birth_female = 0.49, parity = TRUE)), slow = TRUE),

    # multistage, two sex, time varying (Butterick)
    joe_svk_1yr          = s("kin_multi_stage_time_variant_2sex", "Slovakia 1990, one year (package test)",
                             function() q(kin_multi_stage_time_variant_2sex(
                               list(svk_pxs), list(svk_pxs), list(svk_fxs / 0.49), list(svk_fxs / 0.49), list(svk_Uxs), list(svk_Uxs),
                               birth_female = 0.49, parity = TRUE, summary_kin = FALSE, sex_Focal = "Female",
                               initial_stage_Focal = 1, output_years = 1990)), slow = TRUE),
    joe_edu              = s("kin_multi_stage_time_variant_2sex", "Singapore education 2020-2090 (vignette 2_2)",
                             function() q(kin_multi_stage_time_variant_2sex(
                               U_list_females = U_mat_fem_edu, U_list_males = U_mat_male_edu,
                               F_list_females = F_mat_fem_edu, F_list_males = F_mat_male_edu,
                               T_list_females = T_mat_fem_edu, T_list_males = T_mat_fem_edu,
                               birth_female = 1 / (1.06 + 1), parity = FALSE, summary_kin = TRUE,
                               sex_Focal = "Female", initial_stage_Focal = 1, output_years = 1:15)), slow = TRUE),
    joe_edu_male         = s("kin_multi_stage_time_variant_2sex", "male Focal, stage 2, selected kin",
                             function() q(kin_multi_stage_time_variant_2sex(
                               U_list_females = U_mat_fem_edu[1:4], U_list_males = U_mat_male_edu[1:4],
                               F_list_females = F_mat_fem_edu[1:4], F_list_males = F_mat_male_edu[1:4],
                               T_list_females = T_mat_fem_edu[1:4], T_list_males = T_mat_male_edu[1:4],
                               birth_female = c(.49, .5), parity = FALSE, summary_kin = FALSE, sex_Focal = "Male",
                               initial_stage_Focal = 2, output_kin = c("Focal", "d", "m", "coa"), output_years = 1:4))),
    joe_parity_uk_3yrs   = s("kin_multi_stage_time_variant_2sex", "UK parity 1965-1967 (vignette 2_2, shortened)",
                             function() q(kin_multi_stage_time_variant_2sex(
                               U_list_females = Female_parity_mortality_list_UK[1:3], U_list_males = Male_parity_mortality_list_UK[1:3],
                               F_list_females = Female_parity_fert_list_UK[1:3], F_list_males = Female_parity_fert_list_UK[1:3],
                               T_list_females = Parity_transfers_by_age_list_UK[1:3], T_list_males = Parity_transfers_by_age_list_UK[1:3],
                               birth_female = 1 - 0.51, parity = TRUE, output_kin = c("d", "oa", "ys", "os"), summary_kin = TRUE,
                               sex_Focal = "Female", initial_stage_Focal = 1, output_years = 1965:1967)), slow = TRUE),

    # helpers
    rename_kin_f         = s("rename_kin", "female labels on kin() summary",
                             function() q(rename_kin(kin(p = p_2015, f = f_2015)$kin_summary, sex = "f"))),
    rename_kin_2sex      = s("rename_kin", "two-sex labels on kin2sex() summary",
                             function() q(rename_kin(kin2sex(fra_pf, fra_pm, fra_ff, fra_fm)$kin_summary, sex = "2sex"))),
    apc_cohort          = s("output_period_cohort_combination", "cohorts 1950 and 1980",
                             function() q(output_period_cohort_combination(output_cohort = c(1950, 1980),
                                                                           age = 0:100, years_data = 1900:2018))),
    apc_period           = s("output_period_cohort_combination", "periods, 5-year data",
                             function() q(output_period_cohort_combination(output_period = c(1950, 1980),
                                                                           age = 0:100, years_data = seq(1900, 2015, 5))))
  )
}

# Runs every call with one version of the package, in a separate R process.
run_version <- function(pkg_dir, save_dir, scenarios, keep) {
  suppressPackageStartupMessages(pkgload::load_all(pkg_dir, export_all = FALSE, quiet = TRUE))
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  sc <- scenarios()[keep]
  secs <- setNames(numeric(length(sc)), names(sc))
  for (nm in names(sc)) {
    cat(sprintf("  %-22s", nm))
    gc()
    t0 <- proc.time()[["elapsed"]]
    res <- tryCatch(sc[[nm]]$call(), error = function(e) structure(conditionMessage(e), class = "call_error"))
    secs[nm] <- proc.time()[["elapsed"]] - t0
    cat(sprintf("%7.2f s%s\n", secs[nm], if (inherits(res, "call_error")) paste("  ERROR:", res) else ""))
    saveRDS(res, file.path(save_dir, paste0(nm, ".rds")), compress = FALSE)
  }
  secs
}

# data.table keeps a pointer to itself in ".internal.selfref"; it is never
# equal between sessions, so it is removed before comparing.
strip_selfref <- function(x) {
  if (data.table::is.data.table(x)) attr(x, ".internal.selfref") <- NULL
  if (is.list(x) && !is.data.frame(x)) x[] <- lapply(x, strip_selfref)
  x
}

# Largest absolute difference across all numeric columns / matrices.
max_abs_diff <- function(a, b) {
  la <- if (is.list(a)) rapply(a, function(z) z, classes = "numeric", how = "unlist") else as.numeric(a)
  lb <- if (is.list(b)) rapply(b, function(z) z, classes = "numeric", how = "unlist") else as.numeric(b)
  if (length(la) != length(lb)) return(NA_real_)
  d <- abs(la - lb)
  if (all(is.na(d))) 0 else max(d, na.rm = TRUE)
}

# ---------------------------------------------------------------------------
sc_info <- local({
  # describe the scenarios without running them (needs only the bundled data)
  suppressPackageStartupMessages(pkgload::load_all(new_dir, export_all = FALSE, quiet = TRUE, attach = TRUE))
  sc <- scenarios()
  pkgload::unload("DemoKin")
  data.frame(scenario = names(sc), fun = vapply(sc, `[[`, "", "fun"),
             desc = vapply(sc, `[[`, "", "desc"), slow = vapply(sc, `[[`, NA, "slow"))
})
keep <- if (quick) sc_info$scenario[!sc_info$slow] else sc_info$scenario

unlink(out_dir, recursive = TRUE)
dir.create(out_dir)
cat(sprintf("Comparing %d calls (%s mode)\n\nOriginal DemoKin (%s):\n", length(keep), if (quick) "quick" else "full", old_dir))
t_old <- callr::r(run_version, args = list(old_dir, file.path(out_dir, "original"), scenarios, keep), show = TRUE)
cat(sprintf("\nDemoKin_fast (%s):\n", new_dir))
t_new <- callr::r(run_version, args = list(new_dir, file.path(out_dir, "fast"), scenarios, keep), show = TRUE)

report <- character()
say <- function(...) { txt <- paste0(...); cat(txt, "\n"); report <<- c(report, txt) }
rows <- list()
say("\n================ Differences ================")
for (nm in keep) {
  a <- strip_selfref(readRDS(file.path(out_dir, "original", paste0(nm, ".rds"))))
  b <- strip_selfref(readRDS(file.path(out_dir, "fast", paste0(nm, ".rds"))))
  same <- identical(a, b)
  if (!same) {
    say("\n--- ", nm, " (", sc_info$desc[sc_info$scenario == nm], ") ---")
    say(paste(format(waldo::compare(a, b, x_arg = "original", y_arg = "fast", max_diffs = 10)), collapse = "\n"))
  }
  rows[[nm]] <- data.frame(scenario = nm,
                           fun = sc_info$fun[sc_info$scenario == nm],
                           description = sc_info$desc[sc_info$scenario == nm],
                           identical = same,
                           max_abs_diff = if (same) 0 else max_abs_diff(a, b),
                           sec_original = round(t_old[[nm]], 2),
                           sec_fast = round(t_new[[nm]], 2),
                           speedup = round(t_old[[nm]] / max(t_new[[nm]], 0.01), 1))
}
summary_tab <- do.call(rbind, rows)
rownames(summary_tab) <- NULL
if (all(summary_tab$identical)) say("None: every output is identical.")

# Self-check: the comparison must be able to see the smallest possible change.
# Take one output, move a single non-zero value by one unit in the last place,
# and confirm it is reported as different.
probe <- strip_selfref(readRDS(file.path(out_dir, "original", paste0(keep[1], ".rds"))))
probe_df <- if (is.data.frame(probe)) probe else probe[[1]]
num_col <- names(probe_df)[vapply(probe_df, is.double, NA)][1]
i <- which(is.finite(probe_df[[num_col]]) & probe_df[[num_col]] != 0)[1]
nudged <- probe_df
nudged[[num_col]][i] <- nudged[[num_col]][i] * (1 + .Machine$double.eps)
say(sprintf("\nSelf-check: changing %s[%d] in '%s' by 1 ulp (%.17g -> %.17g) is detected: %s",
            num_col, i, keep[1], probe_df[[num_col]][i], nudged[[num_col]][i], !identical(probe_df, nudged)))

say("\n================ Summary ================")
say(paste(utils::capture.output(print(summary_tab[, c("scenario", "fun", "identical", "max_abs_diff",
                                                      "sec_original", "sec_fast", "speedup")], row.names = FALSE)),
          collapse = "\n"))
say(sprintf("\n%d of %d calls identical. Total time: original %.1f s, fast %.1f s.",
            sum(summary_tab$identical), nrow(summary_tab), sum(summary_tab$sec_original), sum(summary_tab$sec_fast)))
say("Run on ", format(Sys.time(), "%Y-%m-%d %H:%M"), ", ", R.version.string)

utils::write.csv(summary_tab, file.path(out_dir, "summary.csv"), row.names = FALSE)
writeLines(report, file.path(out_dir, "report.txt"))
unlink(file.path(out_dir, c("original", "fast")), recursive = TRUE)   # outputs are large; keep only the summary
cat("\nSaved:", file.path(out_dir, "summary.csv"), "and report.txt\n")
