# Outputs must be identical (values, attributes, classes) to DemoKin 1.0.3.
# Fixtures were produced by running the same calls with the original package
# (see dev/equivalence/ for the full 44-scenario harness and how to regenerate).

fixture <- function(name) readRDS(test_path("fixtures", paste0(name, ".rds")))

strip_selfref <- function(x) {
  if (is.data.frame(x)) { attr(x, ".internal.selfref") <- NULL; return(x) }
  if (is.list(x)) return(lapply(x, strip_selfref))
  x
}

expect_same_as_v1 <- function(new, name) {
  ref <- fixture(name)
  expect_identical(strip_selfref(new), strip_selfref(ref))
}

test_that("one-sex time invariant kin() with grouped kin and selected ages matches v1", {
  res <- suppressMessages(kin(p = swe_px[, "2015"], f = swe_asfr[, "2015"],
                              output_kin = c("m", "d", "s", "a", "gm"),
                              output_age_focal = c(0, 30, 65), birth_female = 0.49))
  expect_same_as_v1(res, "ti_1sex_kinsel")
})

test_that("one-sex time varying kin_time_variant() period output matches v1", {
  yr <- as.character(2005:2015)
  res <- suppressMessages(kin_time_variant(p = swe_px[, yr], f = swe_asfr[, yr],
                                           output_period = 2010, output_kin = c("d", "m")))
  expect_same_as_v1(res, "tv_1sex_direct_df")
})

test_that("one-sex time varying kin() with output_year_last matches v1", {
  yr <- as.character(2000:2018)
  res <- suppressMessages(kin(p = swe_px[, yr], f = swe_asfr[, yr], time_invariant = FALSE,
                              output_period = c(2017, 2018), output_year_last = TRUE))
  expect_same_as_v1(res, "tv_1sex_yearlast")
})

test_that("two-sex time invariant kin2sex() with male Focal matches v1", {
  res <- suppressMessages(kin2sex(fra_surv_sex[, "pf"], fra_surv_sex[, "pm"],
                                  fra_asfr_sex[, "ff"], fra_asfr_sex[, "fm"], sex_focal = "m",
                                  output_kin = c("d", "s", "c"), output_age_focal = c(20, 50)))
  expect_same_as_v1(res, "ti_2sex_male")
})

test_that("two-sex time invariant kin2sex() with 3 causes of death matches v1", {
  Hf3 <- matrix(c(.5, 1, .25), nrow = 3, ncol = 101)
  Hm3 <- matrix(c(.33, 1, .4), nrow = 3, ncol = 101)
  Hf3[, 1:40] <- Hm3[, 1:40] <- 1
  res <- suppressMessages(kin2sex(fra_surv_sex[, "pf"], fra_surv_sex[, "pm"],
                                  fra_asfr_sex[, "ff"], fra_asfr_sex[, "fm"], Hf = Hf3, Hm = Hm3,
                                  sex_focal = "m", output_kin = c("m", "d", "s"), output_age_focal = c(30, 70)))
  expect_same_as_v1(res, "ti_cod_3causes_male")
})

test_that("two-sex time varying kin2sex() with causes of death matches v1", {
  yrs <- as.character(2000:2008)
  years <- ncol(swe_px); ages <- nrow(swe_px)
  swe_pm <- swe_px^1.5
  swe_fm <- rbind(matrix(0, 5, years), swe_asfr[-((ages - 4):ages), ]) * 1.05
  dimnames(swe_fm) <- dimnames(swe_asfr)
  Hf3 <- matrix(c(.5, 1, .25), nrow = 3, ncol = 101)
  Hm3 <- matrix(c(.33, 1, .4), nrow = 3, ncol = 101)
  Hf3[, 1:40] <- Hm3[, 1:40] <- 1
  Hf_list <- stats::setNames(lapply(seq_along(yrs), function(i) Hf3 * (1 + i / 50)), yrs)
  Hm_list <- stats::setNames(lapply(seq_along(yrs), function(i) Hm3 * (1 + i / 60)), yrs)
  res <- suppressMessages(kin2sex(pf = swe_px[, yrs], pm = swe_pm[, yrs], ff = swe_asfr[, yrs], fm = swe_fm[, yrs],
                                  Hf = Hf_list, Hm = Hm_list, time_invariant = FALSE,
                                  output_period = c(2003, 2005), output_kin = c("m", "s", "gm")))
  expect_same_as_v1(res, "tv_cod_period")
})

test_that("two-sex multistage time varying model (education, male Focal) matches v1", {
  # The education data are Matrix objects: nrow() of them is NULL until the
  # Matrix namespace is loaded, and neither DemoKin 1.0.3 nor this package
  # imports it, so (as in the tutorials, which attach Matrix) load it first.
  # Under devtools::test() it happens to be loaded already; under R CMD check
  # it is not.
  skip_if_not_installed("Matrix")
  loadNamespace("Matrix")
  res <- suppressMessages(kin_multi_stage_time_variant_2sex(
    U_list_females = U_mat_fem_edu[1:4], U_list_males = U_mat_male_edu[1:4],
    F_list_females = F_mat_fem_edu[1:4], F_list_males = F_mat_male_edu[1:4],
    T_list_females = T_mat_fem_edu[1:4], T_list_males = T_mat_male_edu[1:4],
    birth_female = c(.49, .5), parity = FALSE, summary_kin = FALSE,
    sex_Focal = "Male", initial_stage_Focal = 2, output_kin = c("Focal", "d", "m", "coa"), output_years = 1:4))
  expect_same_as_v1(res, "joe_edu_male")
})

test_that("commutation matrix equals matrixcalc::commutation.matrix", {
  skip_if_not_installed("matrixcalc")
  expect_identical(DemoKinFast:::commutation_matrix(6, 110), matrixcalc::commutation.matrix(6, 110))
  expect_identical(DemoKinFast:::commutation_matrix(3, 5), matrixcalc::commutation.matrix(3, 5))
})

# Grouped kin types (s, c, a, n): kin2sex() sums the older/younger pairs. 2.1.0
# replaced the dplyr regrouping with an exact C++ group sum; these pin it to v1.
test_that("two-sex time varying kin2sex() with all grouped kin types matches v1", {
  yr <- as.character(2000:2018)
  pm <- swe_px[, yr]^1.5
  fm <- rbind(matrix(0, 5, length(yr)), swe_asfr[-(97:101), yr]) * 1.05
  dimnames(fm) <- dimnames(swe_asfr[, yr])
  res <- suppressWarnings(suppressMessages(
    kin2sex(pf = swe_px[, yr], pm = pm, ff = swe_asfr[, yr], fm = fm,
            nf = swe_pop[, yr], nm = swe_pop[, yr] * 0.98, time_invariant = FALSE,
            sex_focal = "m", output_kin = c("c", "d", "gd", "gm", "m", "n", "a", "s"),
            output_period = c(2010, 2018), output_age_focal = c(0, 10, 40, 70))))
  expect_same_as_v1(res, "tv_2sex_grouped")
})

test_that("two-sex time invariant kin2sex() with grouped and ungrouped kin matches v1", {
  res <- suppressWarnings(suppressMessages(
    kin2sex(fra_surv_sex[, "pf"], fra_surv_sex[, "pm"], fra_asfr_sex[, "ff"],
            fra_asfr_sex[, "fm"], output_kin = c("s", "a", "n", "m", "ys"),
            output_age_focal = c(5, 35, 65))))
  expect_same_as_v1(res, "ti_2sex_grouped")
})
