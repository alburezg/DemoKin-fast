# Fixtures for the grouped-kin regression tests, produced with the ORIGINAL
# DemoKin 1.0.3 (installed as `DemoKin`). Run from the package root:
#   Rscript dev/regroup/make_fixtures.R
library(DemoKin)
yr <- as.character(2000:2018)
pm <- swe_px[, yr]^1.5
fm <- rbind(matrix(0, 5, length(yr)), swe_asfr[-(97:101), yr]) * 1.05
dimnames(fm) <- dimnames(swe_asfr[, yr])

tv <- suppressWarnings(suppressMessages(
  DemoKin::kin2sex(pf = swe_px[, yr], pm = pm, ff = swe_asfr[, yr], fm = fm,
                   nf = swe_pop[, yr], nm = swe_pop[, yr] * 0.98, time_invariant = FALSE,
                   sex_focal = "m", output_kin = c("c", "d", "gd", "gm", "m", "n", "a", "s"),
                   output_period = c(2010, 2018), output_age_focal = c(0, 10, 40, 70))))
saveRDS(tv, "tests/testthat/fixtures/tv_2sex_grouped.rds")

ti <- suppressWarnings(suppressMessages(
  DemoKin::kin2sex(fra_surv_sex[, "pf"], fra_surv_sex[, "pm"], fra_asfr_sex[, "ff"],
                   fra_asfr_sex[, "fm"], output_kin = c("s", "a", "n", "m", "ys"),
                   output_age_focal = c(5, 35, 65))))
saveRDS(ti, "tests/testthat/fixtures/ti_2sex_grouped.rds")
