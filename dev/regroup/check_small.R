# Grouped-kin kin2sex() and kin() calls: development DemoKinFast vs the original DemoKin.
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))   # run from the package root
strip <- function(x) { if (is.data.frame(x)) attr(x, ".internal.selfref") <- NULL
  if (is.list(x) && !is.data.frame(x)) x <- lapply(x, strip); x }
pf <- fra_surv_sex[, "pf"]; pm <- fra_surv_sex[, "pm"]; ff <- fra_asfr_sex[, "ff"]; fm <- fra_asfr_sex[, "fm"]
spf <- swe_px; sff <- swe_asfr; spm <- swe_px^1.5
sfm <- rbind(matrix(0, 5, ncol(swe_px)), swe_asfr[-(97:101), ]) * 1.05; dimnames(sfm) <- dimnames(swe_asfr)
Hf <- matrix(c(.5, 1), 2, 101); Hm <- matrix(c(.33, 1), 2, 101); Hf[, 1:50] <- Hm[, 1:50] <- 1
yrs <- as.character(1990:2018)
Hfl <- setNames(lapply(seq_along(yrs), \(i) Hf * (1 + i / 100)), yrs)
Hml <- setNames(lapply(seq_along(yrs), \(i) Hm * (1 + i / 200)), yrs)
calls <- list(
  ti_grouped    = quote(kin2sex(pf, pm, ff, fm, output_kin = c("s", "c", "a", "n", "m"))),
  ti_allgroup   = quote(kin2sex(pf, pm, ff, fm, sex_focal = "m", output_kin = c("c","d","gd","gm","m","n","a","s"))),
  ti_partial    = quote(kin2sex(pf, pm, ff, fm, output_kin = c("s", "os", "m"))),
  ti_cod        = quote(kin2sex(pf, pm, ff, fm, Hf = Hf, Hm = Hm, output_kin = c("s", "c", "m"))),
  tv_period     = quote(kin2sex(spf[, yrs], spm[, yrs], sff[, yrs], sfm[, yrs], time_invariant = FALSE,
                                output_kin = c("c","d","gd","gm","m","n","a","s"), output_period = c(2000, 2015))),
  tv_cohort     = quote(kin2sex(spf[, yrs], spm[, yrs], sff[, yrs], sfm[, yrs], time_invariant = FALSE,
                                sex_focal = "m", output_kin = c("s", "a", "m"), output_cohort = c(1990, 2000))),
  tv_agefocal   = quote(kin2sex(spf[, yrs], spm[, yrs], sff[, yrs], sfm[, yrs], time_invariant = FALSE,
                                output_kin = c("n", "c"), output_period = 2010, output_age_focal = c(0, 30, 60))),
  tv_nosummary  = quote(kin2sex(spf[, yrs], spm[, yrs], sff[, yrs], sfm[, yrs], time_invariant = FALSE,
                                output_kin = c("s", "c"), summary_kin = FALSE)),
  tv_lastyear   = quote(kin2sex(spf[, yrs], spm[, yrs], sff[, yrs], sfm[, yrs], time_invariant = FALSE,
                                output_kin = c("s", "gm"), output_period = 2018, output_year_last = TRUE)),
  tv_cod        = quote(kin2sex(spf[, yrs], spm[, yrs], sff[, yrs], sfm[, yrs], time_invariant = FALSE,
                                Hf = Hfl, Hm = Hml, output_kin = c("s", "c", "m"), output_period = c(2000, 2010)))
)
p15 <- swe_px[, "2015"]; f15 <- swe_asfr[, "2015"]
calls <- c(calls, list(
  k1_ti_grouped   = quote(kin(p15, f15, output_kin = c("s", "c", "a", "n", "m", "gm"))),
  k1_ti_agefocal  = quote(kin(p15, f15, output_kin = c("s", "c"), output_age_focal = c(0, 25, 80))),
  k1_ti_partial   = quote(kin(p15, f15, output_kin = c("s", "ys", "a"), summary_kin = FALSE)),
  k1_tv_period    = quote(kin(swe_px[, yrs], swe_asfr[, yrs], time_invariant = FALSE,
                              output_kin = c("s", "c", "a", "n", "d"), output_period = c(2000, 2015))),
  k1_tv_cohort    = quote(kin(swe_px[, yrs], swe_asfr[, yrs], time_invariant = FALSE,
                              output_kin = c("c", "m"), output_cohort = c(1990, 2000))),
  k1_tv_nosummary = quote(kin(swe_px[, yrs], swe_asfr[, yrs], time_invariant = FALSE,
                              output_kin = c("s", "n"), summary_kin = FALSE)),
  k1_tv_lastyear  = quote(kin(swe_px[, yrs], swe_asfr[, yrs], time_invariant = FALSE,
                              output_kin = c("a", "gm"), output_period = 2018, output_year_last = TRUE))
))
res <- lapply(names(calls), function(nm) {
  a <- suppressWarnings(suppressMessages(eval(calls[[nm]], list(kin2sex = DemoKin::kin2sex, kin = DemoKin::kin))))
  t1 <- system.time(b <- suppressWarnings(suppressMessages(eval(calls[[nm]], list(kin2sex = kin2sex, kin = kin)))))[[3]]
  ok <- identical(strip(a), strip(b))
  why <- if (ok) "" else paste(head(all.equal(strip(a), strip(b)), 2), collapse = " | ")
  data.frame(call = nm, identical = ok, fast_s = round(t1, 2), note = substr(why, 1, 110))
})
print(do.call(rbind, res), right = FALSE)
